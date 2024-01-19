package scratch.ned.nshm23;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSource;
import scratch.ned.nshm23.CEUS_FaultZones_creator.NSHM_SrcZoneData;

public class AK_FaultZones_creator {
	
	final static boolean D=true;
	
	final static double upperSesDepth = 1;
	final static double lowerSeisDepth = 15;
	
	// Following two are from Peter Powers email on Jan 8, 2024

	// Map of zone names and moment rates in N-M/yr
	static final Map<String, Double> ZONE_RATE_MAP = Map.of(
	      "Coast Shear Zone", 4.38E+17,
	      "Cook Inlet", 2.18E+18,
	      "Kodiak Shelf", 5.96E+17,
	      "Northern Alaska Range", 9.45E+17,
	      "Pamplona", 1.96E+18);

	// Map of moment rate in each zone from M>6.5 gridded seismicity
	static final Map<String, Double> ZONE_GRID_RATE_MAP = Map.of(
	      "Coast Shear Zone", 4.05E+16,
	      "Cook Inlet", 1.23E+17,
	      "Kodiak Shelf", 6.58E+16,
	      "Northern Alaska Range", 4.45E+17,
	      "Pamplona", 4.38E+17);

	
	
	public static void getFaultSourceZoneList(String nshmModelDirPath, double duration, boolean lineSource) {
		
		NshmErf erf=null;
		if(D) erf=getNshmERF(nshmModelDirPath);
		
		String geojsonPath, csvPath;
		
		// Cook Inlet; id=4602
		double netBranchWt=1;
		geojsonPath = nshmModelDirPath+"active-crust/zone/Cook Inlet/cook-inlet.geojson";
		csvPath = nshmModelDirPath+"active-crust/zone/Cook Inlet/cook-inlet.csv";
		processSourceZone(geojsonPath, csvPath, netBranchWt, erf, duration, lineSource);
//		System.exit(-1);
		
		// Kodiak Shelf; id=4603
		netBranchWt=1;
		geojsonPath = nshmModelDirPath+"active-crust/zone/Kodiak Shelf/kodiak-shelf.geojson";
		csvPath = nshmModelDirPath+"active-crust/zone/Kodiak Shelf/kodiak-shelf.csv";
		processSourceZone(geojsonPath, csvPath, netBranchWt, erf, duration, lineSource);
	
		// Coast Shear Zone; id=4600
		netBranchWt=0.3333;
		geojsonPath = nshmModelDirPath+"active-crust/zone/Coast Shear Zone/zone/coast-shear-zone.geojson";
		csvPath = nshmModelDirPath+"active-crust/zone/Coast Shear Zone/zone/coast-shear-zone.csv";
		processSourceZone(geojsonPath, csvPath, netBranchWt, erf, duration, lineSource);
		
		// Northern Alaska Range; id=4604
		netBranchWt=0.3333;
		geojsonPath = nshmModelDirPath+"active-crust/zone/Northern Alaska Range/zone/northern-alaska-range.geojson";
		csvPath = nshmModelDirPath+"active-crust/zone/Northern Alaska Range/zone/northern-alaska-range.csv";
		processSourceZone(geojsonPath, csvPath, netBranchWt, erf, duration, lineSource);
		
		// Pamplona; id=4601
		netBranchWt=0.3333;
		geojsonPath = nshmModelDirPath+"active-crust/zone/Pamplona/zone/pamplona.geojson";
		csvPath = nshmModelDirPath+"active-crust/zone/Pamplona/zone/pamplona.csv";
		processSourceZone(geojsonPath, csvPath, netBranchWt, erf, duration, lineSource);

	}
	
	
	private static AreaSourceNSHM23 processSourceZone(String geojsonPath, String csvPath,
			double netBranchWt, NshmErf erf, double duration, boolean lineSource) {

		NSHM_SrcZoneData data = new NSHM_SrcZoneData(geojsonPath, csvPath);
		String nshmName = data.name;
//		double strike = data.strike;
		Region zonePolygon = data.region;
		int id = data.id;
		double testTotRateAtZeroMagBin = data.totalRate;
		
		FocalMechanism focalMech=null;
		if(id==4600) {// the only strike-slip fault; Coast Shear Zone
			focalMech = new FocalMechanism(data.strike,90d,0d); // vertical strike slip
		}
		else {
			focalMech = new FocalMechanism(data.strike,50d,90d); // "REVERSE"
		}

		
		double moRate = ZONE_RATE_MAP.get(data.name) - ZONE_GRID_RATE_MAP.get(data.name);
		IncrementalMagFreqDist mfd = data.magProbDist.deepClone();
		mfd.scaleToTotalMomentRate(moRate*netBranchWt);

		double bValCalc = mfd.compute_bValue();
		double computedIncrRateAtZeroMagBin = mfd.getY(0) * Math.pow(10, -bValCalc*mfd.getX(0))/netBranchWt;
		double meanBranchRate = mfd.getTotalIncrRate(); // this includes netBranchWt
		double meanRate = meanBranchRate/netBranchWt; // to compare to csv file

		if(D) {

			System.out.println(data.name+"\n\tid = "+data.id+"\n\tStrike = "+data.strike+
					"\n\tmoRate+"+(float)moRate+"; moRateTest+"+(float)mfd.getTotalMomentRate()+
					"\n\tminMag = "+mfd.getX(0)+
					"\n\tmaxMag = "+mfd.getMaxMagWithNonZeroRate()+
					"\n\tbValCalc="+ (float)bValCalc+
					"\n\tmeanRate="+meanRate+
					"\n\tcomputedIncrRateAtZeroMagBin="+computedIncrRateAtZeroMagBin+
					"\n\ttestTotRateAtZeroMagBin="+testTotRateAtZeroMagBin+
					"\n\n\tRatio="+ (computedIncrRateAtZeroMagBin/testTotRateAtZeroMagBin)+"\n");

//			System.out.println("\tMag Prob Distribution:");
//			IncrementalMagFreqDist mpd = data.magProbDist;
//			for(int i=0; i<mpd.size();i++) {
//				if(mpd.getY(i)>0.0) {
//					System.out.println("\t"+(float)mpd.getX(i)+"\t"+(float)mpd.getY(i));
//				}
//			}
//			System.out.println("\tZone Polygon:");
//			for(Location loc: zonePolygon.getBorder()) {
//				System.out.println("\t"+loc.getLatitude()+"\t"+loc.getLongitude());
//			}
//			System.out.println("\n");
		}
		
//		if(Math.abs(meanRate/testTotRateAtZeroMagBin - 1.0) > 0.02)
//			throw new RuntimeException("rate discrepancy more then 2% for "+nshmName+"; rate = "+meanRate+"; rate from Peter's file = "+testTotRateAtZeroMagBin);

		// test against ERF
		LocationList erfLocList = new LocationList();
		if(erf != null && D) {
			double erfRate=0;
		    for(int s=0;s<erf.getNumSources();s++) {
		    	NshmSource src = (NshmSource)erf.getSource(s);
		    	if(src.getNSHM_ID() == id) {
		    		// this only works while Peter treats these as point sources
			    	Location loc = src.getRupture(0).getRuptureSurface().getEvenlyDiscretizedLocation(0);
		    		erfLocList.add(loc);
			    	for(int r=0;r<src.getNumRuptures();r++)
			    		erfRate += src.getRupture(r).getMeanAnnualRate(1.0);
		    	}
		    }		
//		    double erfRateTest = Math.abs(1.0-erfRate/(testTotRateAtZeroMagBin*netBranchWt));
//		    if(erfRateTest>1e-4)
//		    	throw new RuntimeException("csvTotRate inconsistent with ERF rate for "+nshmName);
		    
			System.out.println("\terfRate:"+(float)erfRate+"\tTestRatio:"+(float)(meanBranchRate/erfRate)+"\n");
		}
		
		AreaSourceNSHM23 src = new AreaSourceNSHM23(id, nshmName, zonePolygon, 0.1, 
				mfd, upperSesDepth, lowerSeisDepth, focalMech, duration, lineSource);
		src.setTectonicRegionType(TectonicRegionType.ACTIVE_SHALLOW);
		
		// test source
		double srcRateTest = 0;
		for(int r=0;r<src.getNumRuptures();r++)
			srcRateTest+= src.getRupture(r).getMeanAnnualRate(duration);
		if(Math.abs(srcRateTest/meanBranchRate - 1.0)> 0.001)
			throw new RuntimeException("Problem with source rate");

		// test number of gridded-region points
		GriddedRegion gridRegion = src.getGriddedRegion();
		if(D && (gridRegion.getNumLocations() != data.locList.size())) {
			double ratio = (double)gridRegion.getNumLocations()/(double)data.locList.size();
			System.out.println("\tDiff number of locs: here="+gridRegion.getNumLocations()+
					"; nshm="+data.locList.size()+"\tRatio="+(float)ratio+"\n");
		}
		
		return src;
		
	}

	
	 /**
	  * This returns the ERF with timespan duration set to 1.0
	  * @param nshmModelDirPath
	  * @return
	  */
	public static NshmErf getNshmERF(String nshmModelDirPath) {
	    Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW);
	    NshmErf erf = new NshmErf(Path.of(nshmModelDirPath), trts, IncludeBackgroundOption.EXCLUDE);
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
	    return erf;
	}



	public static void main(String[] args) {
		
//		System.exit(0);

		String nshmModelDirPath= "/Users/field/nshm-haz_data/nshm-alaska-main_Jan10_2024/";
		double duration = 50;
		boolean lineSource = false;
		getFaultSourceZoneList(nshmModelDirPath, duration, lineSource);
	}

}
