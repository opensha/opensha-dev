package scratch.ned.nshm23;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.annotation.Nullable;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSource;
import scratch.ned.nshm23.CEUS_FSS_creator.FaultModelEnum;

public class AK_FSS_creator {
	
	final static boolean D = true;
	
	 public enum DeformationModelEnum {
		 ELLIOT(0.3333), 
		 GEO(0.6667), 
		 BOTH(1.0) ;
		 private double weight;
		 DeformationModelEnum(double weight) {
			 this.weight=weight;
		 }
		 double getWeight() {
			 return weight;
		 }
	 }
	
	public static ArrayList<FaultSystemSolution> getFaultSystemSolutionList(String nshmModelDirPath, DeformationModelEnum defModel) {
		
		ArrayList<FaultSystemSolution> fssList = new ArrayList<FaultSystemSolution>();
		
		HashMap<Integer,Double> elliotSlipRateMap = new HashMap<Integer,Double>();
		HashMap<Integer,Double> geoSlipRateMap = new HashMap<Integer,Double>();
		ArrayList<GeoJSONFaultSection> fltSectList = getIsolatedFaultSectionsList(nshmModelDirPath, elliotSlipRateMap, geoSlipRateMap);
		
		HashMap<Integer,String> fltNameFromIdMap = new HashMap<Integer,String>();
		for(GeoJSONFaultSection fltSect:fltSectList)
			fltNameFromIdMap.put(fltSect.getSectionId(), fltSect.getName());
		
		HashMap<Integer,Double> moRateFromIdMap = new HashMap<Integer,Double>();
		HashMap<Integer,Double> rateWtFromIdMap = new HashMap<Integer,Double>(); // to adjust ERF rates for the Elliot and Geo branches

		for(GeoJSONFaultSection fltSect:fltSectList) {
			
			double slipRate = Double.NaN;
			double wtAveSlipRate = DeformationModelEnum.ELLIOT.getWeight()*elliotSlipRateMap.get(fltSect.getSectionId()) + 
					DeformationModelEnum.GEO.getWeight()*geoSlipRateMap.get(fltSect.getSectionId());
			switch (defModel) {
			case ELLIOT:
				slipRate = elliotSlipRateMap.get(fltSect.getSectionId());
				break;  
			case GEO:
				slipRate = geoSlipRateMap.get(fltSect.getSectionId());
				break;  
			case BOTH:
				slipRate = wtAveSlipRate;
				break;
			default:
				throw new RuntimeException("Problem with deformation model");
			}
			
			fltSect.setAveSlipRate(slipRate);
			moRateFromIdMap.put(fltSect.getSectionId(), fltSect.calcMomentRate(false));
			rateWtFromIdMap.put(fltSect.getSectionId(), slipRate/wtAveSlipRate); // mult ERF rates by this value to get branch-specific rates
		}

		
		NshmErf erf = AK_FaultZones_creator.getNshmERF(nshmModelDirPath);
		
		HashMap<Integer,Double> erfMoRateMap = new HashMap<Integer,Double>();
		HashMap<Integer,Double> erfTotRateMap = new HashMap<Integer,Double>();
		HashMap<Integer,Integer> numSrcForFaultMap = new HashMap<Integer,Integer>();
		for(int s=0;s<erf.getNumSources();s++) {
			NshmSource src = (NshmSource)erf.getSource(s);
if(src.getName().equals("Unnamed fault system source")) // temp fix for Peters ID duplicates
	continue;
			int id = src.getNSHM_ID();
			if(elliotSlipRateMap.containsKey(id)) {
				double totRate =ERF_Calculator.getTotalRateAboveMagForSource(src, 1.0, 0.0) * rateWtFromIdMap.get(id);
				double moRate = ERF_Calculator.getTotalMomentRateForSource(src, 1.0) * rateWtFromIdMap.get(id);
				if(erfMoRateMap.containsKey(id)) {
					erfMoRateMap.replace(id, erfMoRateMap.get(id) + moRate);
					erfTotRateMap.replace(id, erfTotRateMap.get(id) + totRate);
					numSrcForFaultMap.put(id, 1+numSrcForFaultMap.get(id));
//					if(D) System.out.println(src.getName()+", "+id+" has more than one source in ERF");
				}
				else {
					erfMoRateMap.put(id, moRate);
					erfTotRateMap.put(id, totRate);
					numSrcForFaultMap.put(id, 1);
				}
				
//				double totEventRate=0;
				SummedMagFreqDist mfd = new SummedMagFreqDist(6.05,25,0.1);
				double mMinERF=10, mMaxERF=0;
				ArrayList<String> strList = new ArrayList<String>();
				for(int r=0;r<src.getNumRuptures();r++) {
					double mag = src.getRupture(r).getMag();
					double eventRate = src.getRupture(r).getMeanAnnualRate(1.0);
					mfd.add(mfd.getClosestXIndex(mag), eventRate*rateWtFromIdMap.get(id));
//					totEventRate += eventRate;
					if(mMinERF>mag) mMinERF=mag;
					if(mMaxERF<mag) mMaxERF=mag;
//					if(src.getNSHM_ID() == 4283) {
//						String str = "\t"+mag+"\t"+src.getRupture(r).getMeanAnnualRate(1.0);
//						if(!strList.contains(str)) {
//							System.out.println("\t"+mag+"\t"+src.getRupture(r).getMeanAnnualRate(1.0));
//							strList.add(str);
//						}
//					}
				}
				double mMin = mfd.getMinMagWithNonZeroRate();
				double mMax = mfd.getMaxMagWithNonZeroRate();
//				double b = mfd.compute_bValue(mMin, mMax);
				double bAlt = mfd.compute_bValueAlt(mMin, mMax);
				int numMag = (int) Math.round((mMax-mMin)/0.1 + 1); 
				if(numMag==1) bAlt=Double.NaN;
				double moRateFract = moRate/moRateFromIdMap.get(id);
				double mfdMoRateTest = mfd.getTotalMomentRate()/moRate;
				double mfdRateTest = totRate/mfd.getTotalIncrRate();
//				double moRateFractRounded = ((double)Math.round(moRateFract*10000d))/10000d;
//				if(src.getNSHM_ID() == 4283) {
//					System.out.println(s+"\t"+src.getNSHM_ID()+"\t"+fltNameFromIdMap.get(src.getNSHM_ID())+"\n"+mfd);
//					done=true;
//				}
				System.out.println(s+"\t"+src.getNSHM_ID()+"\t"+mMinERF+"\t"+mMaxERF+"\t"+(float)moRate+"\t"+(float)moRateFract+"\t"+
						(float)mfdMoRateTest+"\t"+(float)mfdRateTest+"\t"+(float)bAlt+"\t"+numMag+"\t"+fltNameFromIdMap.get(src.getNSHM_ID()));
			}
		}
		
		for(GeoJSONFaultSection fltSect:fltSectList) {

			FaultSystemSolution fss = CEUS_FSS_creator.getFaultSystemSolution(rateWtFromIdMap.get(fltSect.getSectionId()),fltSect, erf);
			fssList.add(fss);
			if(D) {
				double fssMoRate = fss.getTotalFaultSolutionMomentRate();
				double fssTotRate = fss.getTotalRateForAllFaultSystemRups();

				double moRate = fltSect.calcMomentRate(false);
				double erfMoRate = erfMoRateMap.get(fltSect.getSectionId());
				int id = fltSect.getSectionId();
				if(D) System.out.println(id+"\t"+numSrcForFaultMap.get(id)+" Srces\t"+(float)moRate+"\t"+(float)erfMoRate+
						"\t"+(float)(moRate/erfMoRate)+"\t"+fltSect.getName()+"\t"+(float)fltSect.getTraceLength());
			}

		}
		
		
		
		// GET BIG FSS
		fssList.add(getBigFSS(nshmModelDirPath,defModel,erf));
		
		// Set tectonic region type for each
		for(int s=0;s<fssList.size();s++) {
			FaultSystemSolution fss = fssList.get(s);
		    TectonicRegionType[] trForRupArray = new TectonicRegionType[fss.getRupSet().getNumRuptures()];
		    for(int t=0;t<trForRupArray.length;t++)
			    trForRupArray[t] = TectonicRegionType.ACTIVE_SHALLOW;
		    RupSetTectonicRegimes tectonicRegimes = new RupSetTectonicRegimes(fss.getRupSet(),trForRupArray);
		    fss.getRupSet().addModule(tectonicRegimes);
		}

		
		return fssList;
	}
	
	/**
	 * This assumes the sections.geojson files are identical between deformation models, except for the slip rates.
	 * @param nshmModelDirPath
	 * @param defModel
	 * @return
	 */
	private static ArrayList<GeoJSONFaultSection> getFaultSystemSectionsList(String nshmModelDirPath, DeformationModelEnum defModel) {

		ArrayList<GeoJSONFaultSection> listGeol = new ArrayList<GeoJSONFaultSection>();
		ArrayList<GeoJSONFaultSection> listGeod = new ArrayList<GeoJSONFaultSection>();

		FeatureCollection fcGeol =null, fcGeod =null;
		String filePathStringGeod = nshmModelDirPath+"active-crust/fault/system/geodetic/sections.geojson";
		String filePathStringGeol = nshmModelDirPath+"active-crust/fault/system/geologic/sections.geojson";
		try {
			fcGeod = FeatureCollection.read(new File(filePathStringGeod));
			fcGeol = FeatureCollection.read(new File(filePathStringGeol));
		} catch (IOException e) {
			System.out.println("Problem with input file: "+filePathStringGeod+" or "+filePathStringGeol);
			e.printStackTrace();
		}

		for(Feature feature : fcGeol) {
			listGeol.add(GeoJSONFaultSection.fromNSHMP_HazFeature(feature));
		}
		
		for(Feature feature : fcGeod) {
			listGeod.add(GeoJSONFaultSection.fromNSHMP_HazFeature(feature));
		}
		
		if(defModel == DeformationModelEnum.GEO)
			return listGeol;
		else if (defModel == DeformationModelEnum.ELLIOT)
			return listGeod;
		else { // wt average slip rates
			for(int i=0;i<listGeol.size(); i++) {
				double aveSlipRate = listGeol.get(i).getOrigAveSlipRate()*DeformationModelEnum.GEO.getWeight() +
						listGeod.get(i).getOrigAveSlipRate()*DeformationModelEnum.ELLIOT.getWeight();
				listGeol.get(i).setAveSlipRate(aveSlipRate);
			}
			return listGeol;
		}
	}

	
	private static ArrayList<GeoJSONFaultSection> getIsolatedFaultSectionsList(String nshmModelDirPath, HashMap<Integer,Double> elliotSlipRateMap,
			HashMap<Integer,Double> geoSlipRateMap) {

		ArrayList<GeoJSONFaultSection> list = new ArrayList<GeoJSONFaultSection>();
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Atsaksovluk.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bagley (east).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bagley (middle).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bagley (west).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bancas Point.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bendeleben.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bering Glacier.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Billy Creek.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Boundary.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Broad Pass Thrust.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bulchitna Lake.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Bunco Lake.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Camden Bay (east).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Camden Bay (west).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Cape Cleare.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Castle Mountain.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Cathedral Rapids.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Chedotlothna.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Chirikof.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Chugach - St.Elias.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Cordova.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Council.geojson", elliotSlipRateMap, geoSlipRateMap));
//		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/DOT T Johnson.geojson", elliotSlipRateMap, geoSlipRateMap)); // remove in version 3.0.0
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Dall Mountain.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Decoeli Mountain.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Denali (Holitna).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Esker Creek.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Etches.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Foreland Thrust (Khitrov).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Foreland Thrust (Pamplona).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Granite Point.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Gulf of Alaska.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Hanning Bay.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Icy Point.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Iditarod - Nixon Fork (center).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Iditarod - Nixon Fork (east).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Iditarod - Nixon Fork (west).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kahiltna River.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kaktovik (east).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kaktovik (west).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kaltag.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kantishna Hills (north).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kantishna Hills (south).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kanuti.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kayak.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Kigluaik.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Lewis River.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Malaspina.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Marsh Creek.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/McCallum Creek.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/McLeod Creek.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Middle Fork.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Minto Flats (north).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Minto Flats (south).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Montague Strait.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Narrow Cape.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Natazhat.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/North Cook Inlet - SRS.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Northern Alaska Range.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Paimute (east).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Paimute (west).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Pass Creek.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Patton Bay.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Peters Dome.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Pivot.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Purcell Mountains.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Purkeypile.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Ragged Mountain.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Rude River.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Stevens (east).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Stevens (west).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Suckling Hills.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Susitna Glacier.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Ten Fathom.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Tintina (Medicine Lake - Preacher).geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Transition.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Ugak.geojson", elliotSlipRateMap, geoSlipRateMap));
		list.add(getIsolatedFaultSection(nshmModelDirPath+"active-crust/fault/Yakutat Foothills.geojson", elliotSlipRateMap, geoSlipRateMap));
		return list;
	}
	
	static GeoJSONFaultSection getIsolatedFaultSection(String filePathString, HashMap<Integer,Double> elliotSlipRateMap,
			HashMap<Integer,Double> geoSlipRateMap) {
		
		Feature feature=null;
		try {
			feature = Feature.read(new File(filePathString));
		} catch (IOException e) {
			System.out.println("Problem with input file: "+filePathString);
			e.printStackTrace();
		}
		
		feature.properties.containsKey("rate-map");
		FeatureProperties rateMapProps = (FeatureProperties) feature.properties.get("rate-map");
		FeatureProperties rateProp = (FeatureProperties)rateMapProps.get("ELLIOT");
		double elliotSlipRate = rateProp.getDouble("rate", Double.NaN);
		rateProp = (FeatureProperties)rateMapProps.get("GEO");
		double geoSlipRate = rateProp.getDouble("rate", Double.NaN);
		
		feature.properties.containsKey("rate-map");

		String rateType = feature.properties.getString("rate-type", "");
		Double length = feature.properties.getDouble("length", Double.NaN);
		Double mag = feature.properties.getDouble("magnitude", Double.NaN);
				
		GeoJSONFaultSection fltSect = GeoJSONFaultSection.fromNSHMP_HazFeature(feature);
		
		double lengthDiff = Math.abs(length - fltSect.getTraceLength());
		if(lengthDiff>0.1) // differences greater than 0.1 km
			throw new RuntimeException("PROBLEM: Trace length difference)");
		
		if(D) {
			System.out.println("\n"+fltSect.getName()+", "+fltSect.getSectionId());
			System.out.println("\telliotSlipRate="+elliotSlipRate);
			System.out.println("\tgeoSlipRate="+geoSlipRate);
			System.out.println("\trateType="+rateType);
			System.out.println("\trake="+fltSect.getAveRake());
			System.out.println("\tdip="+fltSect.getAveDip());
			System.out.println("\tlength="+length);
			System.out.println("\tlengthFromTraceDiff="+(float)lengthDiff);
			System.out.println("\tupper-depth="+fltSect.getOrigAveUpperDepth());
			System.out.println("\tlower-depth="+fltSect.getAveLowerDepth());
			System.out.println("\tmagnitude="+mag);
//			double ddw = 1e3*(fltSect.getAveLowerDepth()-fltSect.getOrigAveUpperDepth()); ///Math.sin(fltSect.getAveDip()*Math.PI/180);
			double ddw = 15000;
			length *= 1e3;
			NSHM23_ScalingRelationships sh09_Mod = NSHM23_ScalingRelationships.WIDTH_LIMITED; // this is rake independent
			double magCalc = sh09_Mod.getMag(ddw*length, length, ddw, ddw, Double.NaN);
			if(Math.abs(mag-magCalc)>0.05)
				System.out.println("\tWC_MagCalcDiff="+(float)(mag-magCalc));

			
//		    WC1994_MagLengthRelationship wcMagLength = new WC1994_MagLengthRelationship();
//			System.out.println("\tWC_MagDiffAll="+(float)(mag-wcMagLength.getMedianMag(length)));
//			System.out.println("\tWC_MagDiffSS="+(float)(mag-wcMagLength.getMedianMag(length,0.0)));
//			System.out.println("\tWC_MagDiffRV="+(float)(mag-wcMagLength.getMedianMag(length, 90)));
//			System.out.println("\tWC_MagDiffN="+(float)(mag-wcMagLength.getMedianMag(length, -90)));

		}

		double slipRateDipCorr = 1.0;
		if(rateType.equals("VERTICAL_SLIP")) {
			slipRateDipCorr = 1.0/Math.sin(fltSect.getAveDip()*Math.PI/180d);
			System.out.println("\tslipRateDipCorr="+slipRateDipCorr);
		}
		
		elliotSlipRateMap.put(fltSect.getSectionId(), elliotSlipRate*slipRateDipCorr);
		geoSlipRateMap.put(fltSect.getSectionId(), geoSlipRate*slipRateDipCorr);

		return fltSect;
	}

	public static void forPeter() {
		
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-main_Jan10_2024/";

	    Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW);
	    NshmErf erf = new NshmErf(Path.of(nshmModelDirPath), trts, IncludeBackgroundOption.EXCLUDE);
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
	    
	    // chose a source that shouldn't have a rupture at M 7.7
	    int testID = 4251; // Lewis River with "magnitude": "6.27" in geojson file

	    // print out the following for each source: index    NSHM_ID    tmMin    tmMax    totRate    tsrc.getName()
	    System.out.println("\nindex\tNSHM_ID\tmMin\tmMax\ttotRate\tsrc.getName()");
		for(int s=0;s<erf.getNumSources();s++) {
			NshmSource src = (NshmSource)erf.getSource(s);
			if(src.getNSHM_ID() == testID) {
				double mMin = 10;
				double mMax = 0;
				double totRate=0;
				for(int r=0;r<src.getNumRuptures();r++) {
					double mag = src.getRupture(r).getMag();
					if(mMin>mag) mMin=mag;
					if(mMax<mag) mMax=mag;
					totRate += src.getRupture(r).getMeanAnnualRate(1.0)	;
				}
				System.out.println(s+"\t"+src.getNSHM_ID()+"\t"+mMin+"\t"+mMax+"\t"+(float)totRate+"\t"+src.getName());
			}
		}
	}
	
	/**
	 * This gets the actual AK FSS used in the 2023 NSHM
	 * 
	 * This assumes csv files for both deformation models are identical except for rupture rates.
	 * @param nshmModelDirPath
	 * @param defModel
	 */
	private static FaultSystemSolution getBigFSS(String nshmModelDirPath, DeformationModelEnum defModel, NshmErf erf) {

		ArrayList<GeoJSONFaultSection> sectionList = getFaultSystemSectionsList(nshmModelDirPath, defModel);
		System.out.println(sectionList.size());
		for(int s=0;s<sectionList.size();s++) {
			GeoJSONFaultSection fs = sectionList.get(s);
			if (D) System.out.println(s+"\t"+fs.getSectionId()+"\t"+fs.getParentSectionId()+"\t"+fs.getOrigAveSlipRate()+"\t"+fs.getAveRake()+
					"\t"+fs.getAseismicSlipFactor()+"\t"+fs.getAveDip()+"\t"+fs.getOrigAveUpperDepth()+"\t"+fs.getAveLowerDepth()+
					"\t"+fs.getOrigAveSlipRate()+"\t"+fs.getName());
		}

		// read csv rate file
		File fileGeod = new File(nshmModelDirPath+"active-crust/fault/system/geodetic/ruptures.csv");
		File fileGeol = new File(nshmModelDirPath+"active-crust/fault/system/geologic/ruptures.csv");
		CSVFile<String> csvFileGeod=null, csvFileGeol=null;
		try {
			csvFileGeod = CSVFile.readFile(fileGeod, true);
			csvFileGeol = CSVFile.readFile(fileGeol, true);
			if(csvFileGeod.getNumRows() != csvFileGeol.getNumRows())
				throw new RuntimeException("csvFileGeod and csvFileGeol files have different number of rows");
		} catch (IOException e) {
			e.printStackTrace();
		}

		int numRups=0;
		for(int i=1; i<csvFileGeod.getNumRows(); i++) {
			String firstCell = csvFileGeod.get(i, 0);
			if(firstCell.charAt(0) != '#')
				numRups+=1;
//			else
//				System.out.println(firstCell);
		}

		double[] mags = new double[numRups];
		double[] rakes = new double[numRups];
		double[] rates = new double[numRups];
		double[] rupAreas = new double[numRups];
		double[] rupLengths = new double[numRups];
		List<List<Integer>> sectionForRups = new ArrayList<List<Integer>>();

		int r=0; // rupIndex
		//			int testNum=0;
//System.out.println("Start Test");
		for(int i=3; i<csvFileGeod.getNumRows(); i++) {
			if(csvFileGeod.get(i, 0).charAt(0) == '#') // skip these lines {
				continue;
			mags[r] = csvFileGeod.getDouble(i, 0);
			rakes[r] = csvFileGeod.getDouble(i, 5);
			if(defModel == DeformationModelEnum.GEO)
				rates[r] = csvFileGeol.getDouble(i, 1);
			else if (defModel == DeformationModelEnum.ELLIOT)
				rates[r] = csvFileGeod.getDouble(i, 1);
			else
				rates[r] = csvFileGeod.getDouble(i, 1)*DeformationModelEnum.ELLIOT.getWeight() + csvFileGeol.getDouble(i, 1)*DeformationModelEnum.GEO.getWeight();
			List<Integer> indexList = parseSectionsInRupture(csvFileGeod.get(i, 6));
			sectionForRups.add(indexList);
			//				if(csvFile.get(i, 6).contains("-") && testNum < 100) {
			//					System.out.println(r+"\t"+indexList);
			//					testNum+=1;
			//				}
			double rupLength=0;
			double rupArea=0;
			for(int sectIndex:indexList) {
				rupLength += sectionList.get(sectIndex).getTraceLength();
				rupArea += sectionList.get(sectIndex).getArea(false);
			}
			rupLengths[r]=rupLength*1e3; // km to m
			rupAreas[r]=rupArea; 
			r+=1;
		}
		FaultSystemRupSet rupSet = new FaultSystemRupSet(
				sectionList,
				sectionForRups,
				mags,
				rakes,
				rupAreas,
				rupLengths); 

		FaultSystemSolution fss= new FaultSystemSolution(rupSet, rates);
		fss.setInfoString("Big AK Fault System Solution");
		
		if(D && erf != null) {
	    	SummedMagFreqDist fssTotalMFD = new SummedMagFreqDist(6.05,80,0.05);
	    	double fssTotalMomentRate = 0;
	    	for(int rup=0;rup<fss.getRupSet().getNumRuptures();rup++) {
	    		double rate = fss.getRateForRup(rup);
	    		double mag = fss.getRupSet().getMagForRup(rup);
	    		fssTotalMomentRate += MagUtils.magToMoment(mag)*rate;
	    		int iMag = fssTotalMFD.getClosestXIndex(mag);
	    		fssTotalMFD.add(iMag, rate);
//	    		fssTotalMFD.add(mag, rate);
	    	}
	    	SummedMagFreqDist erfTotalMFD = new SummedMagFreqDist(6.05,80,0.05);
	    	for(int s=0;s<erf.getNumSources();s++) {
	    		NshmSource src = (NshmSource)erf.getSource(s);
	    		if(src.getName().equals("Unnamed fault system source")) {
	    			for(r=0;r<src.getNumRuptures();r++) {
	    				double mag = src.getRupture(r).getMag();
	    	    		int iMag = erfTotalMFD.getClosestXIndex(mag);
	    				erfTotalMFD.add(iMag, src.getRupture(r).getMeanAnnualRate(1.0));
	    			}
	    		}
	    	}
	    	
	    	for(int i=0;i<fssTotalMFD.size();i++) {
	    		double val1 = fssTotalMFD.getY(i);
	    		double val2 = erfTotalMFD.getY(i);
	    		if(val1 == 0.0)
	    			if(erfTotalMFD.getY(i) != val1)
	    				throw new RuntimeException("fssTotalMFD.getY(i) is zero but erfTotalMFD.getY(i) = "+val1);
	    		else if(Math.abs((val1-val2)/val1) > 0.00001)
    				throw new RuntimeException("fssTotalMFD != erfTotalMFD at M="+erfTotalMFD.getX(i)+"; fssTotalMFD="+val1+"; erfTotalMFD="+val2);
	    	}
	    	System.out.println("TEST PASSED: fssTotalMFD equals erfTotalMFD");
	    	
//	    	// the following verifies that Fairweather North has zero participation rates for nshm-alaska-3.0.1
//	    	// test Fairweather North 14 area; [-138.76926, 59.66737], [-138.72111, 59.64056]
//	    	Region f14_region = new Region(new Location(59.66737,-138.76926,0d), new Location(59.64056,-138.72111,0d));
//	    	
//	    	for(FaultSection fs:fss.getRupSet().getFaultSectionDataList()) {
//	    		for(Location loc: fs.getFaultTrace()) {
//	    			if(f14_region.contains(loc))
//	    				System.out.println("INSIDE: "+fs.getName());
//	    		}
//	    	}
//	    	double rateFSS_InRegion=0;
//	    	for(r=0;r<fss.getRupSet().getNumRuptures();r++) {
//	    		double fractInside=fss.getRupSet().getSurfaceForRupture(r, 1.0).getFractionOfSurfaceInRegion(f14_region);
//	    		if(fractInside>0)
//	    			rateFSS_InRegion =+ fss.getRateForRup(r);
//	    	}
//	    	double rateERF_InRegion=0;
//	    	boolean gotOne=false;
//	    	boolean gotOne2=false;
//	    	int numERF_FSS_sources=0;
//	    	for(int s=0;s<erf.getNumSources();s++) {
//	    		NshmSource src = (NshmSource)erf.getSource(s);
//	    		if(src.getName().equals("Unnamed fault system source")) {
//	    			gotOne=true;
//	    			numERF_FSS_sources+=1;
//	    			for(r=0;r<src.getNumRuptures();r++) {
//	    	    		double fractInside=src.getRupture(r).getRuptureSurface().getFractionOfSurfaceInRegion(f14_region);
//	    	    		if(fractInside>0) {
//	    	    			gotOne2=true;
//	    	    			rateERF_InRegion += src.getRupture(r).getMeanAnnualRate(erf.getTimeSpan().getDuration());
//	    	    		}
//	    			}
//	    		}
//	    	}
//	    	System.out.println("rateFSS_InRegion="+rateFSS_InRegion);
//	    	System.out.println("rateERF_InRegion="+rateERF_InRegion+"\t gotOne="+gotOne+"\t gotOne2="+gotOne2);
//	    	System.out.println("numERF_FSS_sources="+numERF_FSS_sources);
//	    	System.out.println("fss.getRupSet().getNumRuptures()="+fss.getRupSet().getNumRuptures());

		}
		return fss;
	}
	
	
	private static List<Integer> parseSectionsInRupture(String cellString) {
        String[] dashSplit = cellString.split("-");
        Preconditions.checkState(dashSplit.length >= 1);
        List<Integer> indexList = new ArrayList<>();

        for (String split : dashSplit) {
            String[] colonSplit = split.split(":");
            Preconditions.checkState(colonSplit.length <= 2);
            if(colonSplit.length==1) {
            	indexList.add(Integer.parseInt(colonSplit[0]));
            }
            else {
                int first = Integer.parseInt(colonSplit[0]);
                int last = Integer.parseInt(colonSplit[1]);
                Preconditions.checkState(first != last);
                if (first < last)
                    for (int i=first; i<=last; i++)
                    	indexList.add(i);
                else
                    for (int i=first; i>=last; i--)
                    	indexList.add(i);
            }
         }
        return indexList;
	}
	


	public static void main(String[] args) {
		
//		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-3.0.1/";	
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-3.1-maint/";
		DeformationModelEnum defModel = DeformationModelEnum.ELLIOT;
//		DeformationModelEnum defModel = DeformationModelEnum.GEO;
		ArrayList<FaultSystemSolution> fssList = getFaultSystemSolutionList(nshmModelDirPath, defModel);
		FaultSystemSolution bigFSS = fssList.get(fssList.size()-1);
		double[] partRates = bigFSS.calcParticRateForAllSects(0, 10);
		for(int i=0;i<partRates.length;i++)
			if(partRates[i]==0.0)
				System.out.println("Zero rate for sect "+i+"\t"+bigFSS.getRupSet().getFaultSectionData(i).getName());

		for(FaultSystemSolution fss:fssList) {
			for(int r=0; r<fss.getRupSet().getNumRuptures();r++) {
				if(fss.getRupSet().getAreaForRup(r) == 0) {
					System.out.println("Area zero for r="+r+"\t"+fss.getInfoString());
				}
			}
		}
		//		System.exit(0);

		
//		NshmErf erf = AK_FaultZones_creator.getNshmERF(nshmModelDirPath);
//		getBigFSS(nshmModelDirPath,defModel,erf);

//		forPeter();
//		
//		double length = 34.983*1e3;
//		double ddw = 15*1e3;
//		NSHM23_ScalingRelationships sh09_Mod = NSHM23_ScalingRelationships.WIDTH_LIMITED; // this is rake independent
//		double mag = sh09_Mod.getMag(ddw*length, length, ddw, ddw, Double.NaN);
//		System.out.println("\tmag="+(float)(mag));

		
//	    WC1994_MagLengthRelationship wcMagLength = new WC1994_MagLengthRelationship();
//	    double mag = 7.15;
//	    System.out.println(wcMagLength.getMedianLength(mag));
//	    System.out.println(wcMagLength.getMedianLength(mag,0.0));
//	    System.out.println(wcMagLength.getMedianLength(mag,-90));
//	    System.out.println(wcMagLength.getMedianLength(mag,90));
//
//	    double length = 34.983;
//	    wcMagLength = new WC1994_MagLengthRelationship();
//	    System.out.println(wcMagLength.getMedianMag(length));
//	    System.out.println(wcMagLength.getMedianMag(length,0.0));
//	    System.out.println(wcMagLength.getMedianMag(length,-90));
//	    System.out.println(wcMagLength.getMedianMag(length,90));

		
//		System.exit(0);
		
//		System.out.println(wcMagLength.getMedianMag(131.583));
//		wcMagLength.setRake(90);
//		System.out.println(wcMagLength.getMedianMag(131.583));
//		wcMagLength.setRake(0);
//		System.out.println(wcMagLength.getMedianMag(131.583));
//		wcMagLength.setRake(-90);
//		System.out.println(wcMagLength.getMedianMag(131.583));
//		wcMagLength.setRake(Double.NaN);
//		System.out.println(wcMagLength.getMedianMag(131.583));
		
		
		
	}

}
