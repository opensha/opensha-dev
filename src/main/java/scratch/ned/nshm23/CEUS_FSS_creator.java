package scratch.ned.nshm23;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Set;


import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.json.Feature;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSource;

/**
 * This class creates a list of fault system solutions (FSS) out of the 2023 CEUS fault-based sources.  All fault 
 * sources are in a single FSS except the two Cheraw sources, which are handled separately due to their having 
 * floating ruptures (and these do not have variable depth to top of rupture in the FSS implementations here).
 * 
 * Note that Peter's directory names change in between versions (e.g., nshm-conus-6.b.4 vs nshm-conus-6.0.0)
 * 
 * @author field
 *
 */
public class CEUS_FSS_creator {
	
	final static boolean D = true;  
	
	 private static final Gson GSON = new GsonBuilder().create();
	 
	 public enum FaultModelEnum {
		 BOTH, PREFERRED, ALTERNATE 
	 }
	 
	 /**
	  * This map defines the sources unique to the preferred fault model.
	  * The keys are the NSHM section IDs and the values are the branch wts.
	  * The first part of comment is why it's here and the second part is the source name
	  * @return
	  */
	 private static HashMap<Integer,Double> getUniqueSourcesForPreferredFaultModelMap() {
		 HashMap<Integer,Double> map = new HashMap<Integer,Double>();
		 map.put(2180,0.5); //	2023 update	Cheraw
		 map.put(3801,0.6); //	more recent;	Eastern Rift Margin (South) (crittenden-co)
		 map.put(3053,0.5); //	more recent;	SSCn (AxS-AxN: Axial (south), Axial (north))
		 map.put(3054,0.5); //	more recent;	SSCn (AxS-BL: Axial (south), Bootheel)
		 map.put(3060,0.5);	//	more recent;	New Madrid - SSCn (New Madrid, north)
		 map.put(3062,0.5);	//	more recent;	SSCn (NMN-L: north, extended)
		 map.put(3071,0.5);	//	more recent;	New Madrid - SSCn (Reelfoot, north)
		 map.put(3073,0.5);	//	more recent;	SSCn (RFT-L: Reelfoot, extended)
		 map.put(3610,0.5); //  more recent; 	Meers (USGS)
		 return map;
	 }
	 
	 /**
	  * This map defines the sources unique to the alternative fault model.
	  * The keys are the NSHM section IDs and the values are the branch wts
	  * The first part of comment is why it's here and the second part is the source name

	  * @return
	  */
	 private static HashMap<Integer,Double> getUniqueSourcesForAlternateFaultModelMap() {
		 HashMap<Integer,Double> map = new HashMap<Integer,Double>();
		 map.put(3000,0.5);  // outdated;	New Madrid - USGS (west, north)
		 map.put(3001,0.5);  // outdated;	New Madrid - USGS (west, center)
		 map.put(3002,0.5);  //	outdated;	New Madrid - USGS (west, south)
		 map.put(3003,0.5);  //	outdated;	New Madrid - USGS (west)
		 map.put(3010,0.5);  //	outdated;	New Madrid - USGS (mid-west, north)
		 map.put(3011,0.5);  //	outdated;	New Madrid - USGS (mid-west, center)
		 map.put(3012,0.5);  //	outdated;	New Madrid - USGS (mid-west, south)
		 map.put(3013,0.5);  // outdated;	New Madrid - USGS (mid-west)
		 map.put(3020,0.5);  // outdated;	New Madrid - USGS (center, north)
		 map.put(3021,0.5);  // outdated;	New Madrid - USGS (center, center)
		 map.put(3022,0.5);  // outdated;	New Madrid - USGS (center, south)
		 map.put(3023,0.5);  // outdated;	New Madrid - USGS (center)
		 map.put(3030,0.5);  // outdated;	New Madrid - USGS (mid-east, north)
		 map.put(3031,0.5);  // outdated;	New Madrid - USGS (mid-east, center)
		 map.put(3032,0.5);  // outdated;	New Madrid - USGS (mid-east, south)
		 map.put(3033,0.5);  // outdated;	New Madrid - USGS (mid-east)
		 map.put(3040,0.5);  // outdated;	New Madrid - USGS (east, north)
		 map.put(3041,0.5);  // outdated;	New Madrid - USGS (east, center)
		 map.put(3042,0.5);  // outdated;	New Madrid - USGS (east, south)
		 map.put(3043,0.5);  // outdated;	New Madrid - USGS (east)
		 map.put(2181,0.5);  // outdated;	Cheraw (SSCn)
		 map.put(3805,0.4);	//	lower wt;	Eastern Rift Margin (South) (meeman-shelby)
		 map.put(3611,0.5); //  outdated; 	Meers (SSCn, cluster-out)		doesn't extend to include DOLE paleo site
		 return map;
	 }


	
	 /**
	  * I first tried mapping ruptures to surfaces, but some mapped to both USGS and CUES SSCn surfaces 
	  * (which are on different logic trees).  I now do this from Peter's files.
	  */
	private static void getSurfacesForSources() {
	}
	
	private static ArrayList<GeoJSONFaultSection> getFaultSectionList(String nshmModelDirPath) {
		ArrayList<GeoJSONFaultSection> list = new ArrayList<GeoJSONFaultSection>();
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/CO/Cheraw/features/Cheraw (SSCn).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/CO/Cheraw/features/Cheraw.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/Commerce/features/commerce.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (Axial, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (Axial, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (Bootheel).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (Charleston Uplift).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (New Madrid, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (New Madrid, west).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (Reelfoot, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - SSCn (Reelfoot, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (center, center).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (center, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (center, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (east, center).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (east, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (east, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (mid-east, center).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (mid-east, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (mid-east, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (mid-west, center).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (mid-west, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (mid-west, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (west, center).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (west, north).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/features/New Madrid - USGS (west, south).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/OK/Meers/features/Meers.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/OK/Meers/features/Meers (SSCn).geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (North)/features/eastern-rift-margin-north.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/features/crittenden-county.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/features/eastern-rift-margin-south-extension.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/features/eastern-rift-margin-south.geojson"));
		list.add(getFaultSection(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/features/meeman-shelby.geojson"));
		return list;
	}
	
	
	static GeoJSONFaultSection getFaultSection(String filePathString) {
		Feature feature=null;
		try {
			feature = Feature.read(new File(filePathString));
		} catch (IOException e) {
			System.out.println("Problem with input file: "+filePathString);
			e.printStackTrace();
		}
		return GeoJSONFaultSection.fromNSHMP_HazFeature(feature);
	}
	
	private static boolean areSourceRupSurfacesIdentical(ProbEqkSource src) {
		if(src.getNumRuptures()==1)
			return true;
		RuptureSurface firstSurface = src.getRupture(0).getRuptureSurface();
		int numPtsFirst = firstSurface.getEvenlyDiscretizedNumLocs();
		LocationList firstLocList = firstSurface.getEvenlyDiscritizedListOfLocsOnSurface();
		for(int s=1;s<src.getNumRuptures();s++) {
			RuptureSurface surface = src.getRupture(s).getRuptureSurface();
			if(surface.getEvenlyDiscretizedNumLocs() != numPtsFirst)
				return false;
			LocationList locList = surface.getEvenlyDiscritizedListOfLocsOnSurface();
			if(firstLocList.size() != locList.size())
				return false;
			for(int l=0;l<firstLocList.size();l+=5)
				if(!firstLocList.get(l).equals(locList.get(l)))
					return false;
		}
		
		return true;
	}
	
	
//	public static boolean isPointSource(NshmSource src) {
//    	int numLocs = src.getRupture(0).getRuptureSurface().getEvenlyDiscretizedNumLocs();
//    	if(numLocs == 1) 
//    		return true;
//    	else
//    		return false;
//	}
	
	
	public static ArrayList<FaultSystemSolution> getFaultSystemSolutionList(String nshmModelDirPath, FaultModelEnum fltModel) {
		
		if(D) System.out.println("fltModel = "+fltModel);

		ArrayList<FaultSystemSolution> fssList = new ArrayList<FaultSystemSolution>();
		
		ArrayList<Integer> srcZoneID_List = CEUS_FaultZones_creator.getSourceZoneID_List();

		// get parent fault sections from Peter's features files
		ArrayList<GeoJSONFaultSection> faultSectionData = getFaultSectionList(nshmModelDirPath);
//		try {
//			System.out.println(faultSectionData.get(0).toFeature().toJSON());
//		} catch (IOException e1) {
//			e1.printStackTrace();
//		}

		// make parSectID_List
		ArrayList<Integer> parSectID_List = new ArrayList<Integer>(); // this is NSHM ID for each parent section
		HashMap<Integer,String> parSecNameFromtID_Map = new HashMap<Integer,String>();
		if(D) System.out.println("index\tsectID\trake");
		for(int s=0;s<faultSectionData.size();s++) {
			GeoJSONFaultSection sect = faultSectionData.get(s);
			if(!parSectID_List.contains(sect.getSectionId())) {
				parSectID_List.add(sect.getSectionId());
				parSecNameFromtID_Map.put(sect.getSectionId(),sect.getName());
				if(D) System.out.println(s+"\t"+sect.getSectionId()+"\t"+sect.getAveRake()+"\t"+sect.getName());
			}
			else
				throw new RuntimeException("section IDs are not unique; duplicate: "+sect.getSectionId());
		}
		if (D)System.out.println("parSectID_List.size()="+parSectID_List.size());
		
		// Read from Peter's rupture-set.json and and cluster-set.json files
		HashMap<Integer,int[]> srcFltSectsMap = new HashMap<Integer,int[]>(); // the fault section used by each source (same order as above)
		getSrcIDsAndFaultSectionsLists(srcFltSectsMap, nshmModelDirPath);
		Set<Integer> srcIDsList = srcFltSectsMap.keySet();  // a list of all the source IDs (no duplicates)
		
		// print section(s) for each source
		for(int src:srcIDsList) {
			int[] sectArray = srcFltSectsMap.get(src);
			if(sectArray.length == 1) {
				System.out.println("Src "+src+" only uses:\t"+parSecNameFromtID_Map.get(sectArray[0])+"\n");
			} else {
				System.out.println("Src "+src+" uses:");
				for(int sect:sectArray)
					System.out.println("\t"+parSecNameFromtID_Map.get(sect));
//				System.out.println("\n"+src);
			}	
		}
//		System.exit(0);

		// some tests
		if(D) System.out.println("number of unique source IDs: "+srcIDsList.size());
		// make sure all sections are used and none are missing (with respect to Peter's files)
		ArrayList<Integer> testParSectID_List = new ArrayList<Integer>(); // fill up with the par IDs from all the sources
		for(int srcID:srcFltSectsMap.keySet()) {
//			for(int i=0;i<srcIDsList.size();i++) {
			int[] sects = srcFltSectsMap.get(srcID);
			if(sects.length==1)
				if(sects[0] != srcID) // source ID = faultSection ID if only one fault section used
					throw new RuntimeException("problem");

			if(D)System.out.print("\n"+srcID+"\t");
			for(int sect:sects) {
				if (D)System.out.print(sect+", ");	
				if(!testParSectID_List.contains(sect))
					testParSectID_List.add(sect);
			}
		}
		if (D) System.out.print("\n");
    	if(parSectID_List.size() != testParSectID_List.size())
    		throw new RuntimeException("parSectID_List.size() != testParSectID_List.size()");
    	for(int id : testParSectID_List)
    		if(!parSectID_List.contains(id))
    			throw new RuntimeException("parSectID_List does not contain: "+id);
    	if(D)System.out.println("parSectID_List passed tests");
		
		
		// create the ERF
		NshmErf erf = getNshmERF(nshmModelDirPath); // this excludes gridded seismicity except for that for fault zones
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
		ArrayList<Integer> floaterSrcID_List = new ArrayList<Integer>();  // these will have a separate FSS

	    if (D) System.out.println("NSHM ERF NumSources: " + erf.getNumSources());
		ArrayList<Integer> testSrcIDsList = new ArrayList<Integer>();  // a list of all the source IDs (no duplicates)
	    int numSrcZonePtSources=0; // this is all the gridded seismicity sources for the fault source zones
	    ArrayList<Integer> tempSrcZoneIDs = new ArrayList<Integer>();
	    for(int s=0;s<erf.getNumSources();s++) {
	    	NshmSource src = (NshmSource)erf.getSource(s);
	    	if(srcZoneID_List.contains(src.getNSHM_ID())) { // skip fault zone gridded source zones 
	    		numSrcZonePtSources+=1;
	    		if(!tempSrcZoneIDs.contains(src.getNSHM_ID()))	    	
	    			tempSrcZoneIDs.add(src.getNSHM_ID());
	    		continue;
	    	}
	    	boolean noFloaters = areSourceRupSurfacesIdentical(src); // look for ruptures that have area less than the full fault
	    	if(!noFloaters) {
	    		floaterSrcID_List.add(src.getNSHM_ID());
	    	}
	    	boolean srcIdEqualsSectionId = parSectID_List.contains(src.getNSHM_ID());
	    	if(!testSrcIDsList.contains(src.getNSHM_ID()))
	    		testSrcIDsList.add(src.getNSHM_ID());
	    	if (D)System.out.println(s+"\t"+src.getNumRuptures()+"\t"+noFloaters+"\t"+srcIdEqualsSectionId+"\t"+src.getNSHM_ID()+"\t"+src.getName());
	    }

	    if(D) {
	    	System.out.println("numSrcZonePtSources="+numSrcZonePtSources+"\tnumFltSrces="+(erf.getNumSources()-numSrcZonePtSources));
	    	System.out.println("numSrcZones="+tempSrcZoneIDs.size());
	    	System.out.println("tempSrcZoneIDs:");
	    	for(Integer nameAndID: tempSrcZoneIDs)
	    		System.out.println("\t"+nameAndID);
	    	System.out.println("Floater Source IDs:");
		    for(int id:floaterSrcID_List)
		    	System.out.println("\t"+id);
	    }

	    // run more tests; first make sure ERF sources are consistent with Peter's files
    	if(testSrcIDsList.size() != srcIDsList.size())
    		throw new RuntimeException("testSrcIDsList.size() != srcIDsList.size()");
    	for(int id : testSrcIDsList)
    		if(!srcIDsList.contains(id))
    			throw new RuntimeException("srcIDsList does not contain: "+id);
    	if(D)System.out.println("srcIDsList passed tests");

    	// Compute ERF MFDs for each source
    	HashMap<Integer, SummedMagFreqDist> mfdForSrcIdMap = new HashMap<Integer, SummedMagFreqDist>();
    	HashMap<Integer, Double> rakeForSrcIdMap = new HashMap<Integer, Double>();
    	HashMap<Integer, String> nameForSrcIdMap = new HashMap<Integer, String>();
    	for(int id:srcIDsList) {
        	mfdForSrcIdMap.put(id, getBlankMFD());
    	}
	    for(int s=0;s<erf.getNumSources();s++) {
	    	NshmSource src = (NshmSource)erf.getSource(s);
	    	if(srcZoneID_List.contains(src.getNSHM_ID())) 
	    		continue;	
	    	Integer srcID = src.getNSHM_ID();
	    	if(nameForSrcIdMap.keySet().contains(srcID)) { // check for name change
	    		String firstName = nameForSrcIdMap.get(srcID);
	    		if(!firstName.equals(src.getName()))  // all pass this test
	    			System.out.println("WARNING: Source name change in ERF for ID="+srcID+";\t"+firstName+",\t"+nameForSrcIdMap.get(srcID));
	    	}
	    	else { // add name
	    		nameForSrcIdMap.put(srcID, src.getName());
	    	}
	    	SummedMagFreqDist mfd = mfdForSrcIdMap.get(srcID);
	    	for(int r=0;r<src.getNumRuptures();r++) {
	    		double mag = src.getRupture(r).getMag();
	    		int iMag = mfd.getClosestXIndex(mag);
	    		double rate = src.getRupture(r).getProbability();  // rate approx equal to prob
	    		mfd.add(iMag, rate); // this requires the exact x value (no tolerance)
	    		double rake = src.getRupture(r).getAveRake();
	    		if(rakeForSrcIdMap.containsKey(srcID)) {
	    			if(rakeForSrcIdMap.get(srcID) != rake)
	    				throw new RuntimeException("rake change within source");
	    			else
	    				rakeForSrcIdMap.put(srcID, rake);
	    		}
	    	}
	    }
	    // print min and max mag for each src mfd
	    if(D) {
		    for(int srcID:mfdForSrcIdMap.keySet()) {
		    	SummedMagFreqDist mfd = mfdForSrcIdMap.get(srcID);
		    	System.out.println(srcID+"\t"+(float)mfd.getMinMagWithNonZeroRate()+"\t"+(float)mfd.getMaxMagWithNonZeroRate()+"\trake="+
		    	rakeForSrcIdMap.get(srcID)+"\t"+nameForSrcIdMap.get(srcID));
		    }
	    }
	    
		HashMap<Integer,Double> prefSrcWtmap = getUniqueSourcesForPreferredFaultModelMap();
		HashMap<Integer,Double> altSrcWtmap = getUniqueSourcesForAlternateFaultModelMap();
		
		// write out sources for each branch
		if(D) {
			System.out.println("Preferred model unique sources (ID, Wt, Name):");
			for(int src:prefSrcWtmap.keySet())
				System.out.println("\t"+src+"\t"+prefSrcWtmap.get(src)+"\t"+nameForSrcIdMap.get(src));
			System.out.println("Alternative model unique sources (ID, Wt, Name):");
			for(int src:altSrcWtmap.keySet())
				System.out.println("\t"+src+"\t"+altSrcWtmap.get(src)+"\t"+nameForSrcIdMap.get(src));
			System.out.println("Sources in both models (ID, Name):");
			for(int src:srcFltSectsMap.keySet()) {
				if(prefSrcWtmap.keySet().contains(src) || altSrcWtmap.keySet().contains(src))
					continue; // skip
				System.out.println("\t"+src+"\t"+nameForSrcIdMap.get(src));
			}
			
		}
		
		
		

	    
	    // Make the floater FSSs
	    for(int srcID:floaterSrcID_List) {
	    	
	    	double rateWt = 1.0;
    		if (fltModel == FaultModelEnum.PREFERRED) {
	    		if(prefSrcWtmap.keySet().contains(srcID))  // keep and set rateWt
	    			rateWt = 1.0/prefSrcWtmap.get(srcID);
	    		else if (altSrcWtmap.keySet().contains(srcID)) // skip it (not on this branch)
	    			continue;
    		}
    		else if (fltModel == FaultModelEnum.ALTERNATE) {
	    		if(prefSrcWtmap.keySet().contains(srcID))  // skip it
	    			continue;
	    		else if (altSrcWtmap.keySet().contains(srcID)) // keep and set rateWt
	    			rateWt = 1.0/altSrcWtmap.get(srcID);;
    		}
    		// for case BOTH, this will continue with rateWt=1;
    		if(D)System.out.println("Floater for src "+srcID+" has weight of "+(float) rateWt+" for fault-model branch "+fltModel);
    		
    		
	    	// find fault section (inefficient)
	    	GeoJSONFaultSection fltSect=null;
	    	for(GeoJSONFaultSection sect:faultSectionData) {
	    		if(sect.getSectionId() == srcID) {
	    			fltSect = sect;
	    			break;
	    		}
	    	}

	    	FaultSystemSolution fss_floater = getFaultSystemSolution(rateWt, fltSect, erf);
		    
		    // test total MFD
	    	Boolean testPassed = true;
	    	SummedMagFreqDist fssTotalMFD = getBlankMFD();
	    	for(int rup=0;rup<fss_floater.getRupSet().getNumRuptures();rup++) {
	    		double rate = fss_floater.getRateForRup(rup);
	    		double mag = fss_floater.getRupSet().getMagForRup(rup);
	    		int iMag = fssTotalMFD.getClosestXIndex(mag);
	    		fssTotalMFD.add(iMag, rate);
	    	}
		    // compare
	    	SummedMagFreqDist totSrcMFD = mfdForSrcIdMap.get(srcID);
		    for(int i=0;i<fssTotalMFD.size();i++) {
		    	double val1 = fssTotalMFD.getY(i)/rateWt;
		    	double val2 = totSrcMFD.getY(i);
		    	if(val1 == 0.0) {
		    		if(val2 != 0)
		    			throw new RuntimeException("PROBLEM: zero in one but not the other at index "+i+"\n"+fssTotalMFD+"\n"+fssTotalMFD);
		    	}
		    	else {
		    		double ratio = val1/val2;
		    		if(ratio<0.99 || ratio > 1.01) {
		    			throw new RuntimeException("PROBLEM: >1% difference at index "+i+"\n"+fssTotalMFD+"\n"+totSrcMFD);
		    		}
		    	}
		    }
		    if (D) System.out.println("MFD test for the following floater FSS passed: "+fltSect.getName());

		    fssList.add(fss_floater);
	    }
	    	    
	    // now make the big fss, we need the following without floater sources and with fault section indices that start from zero
	    HashMap<Integer, SummedMagFreqDist> mfdForSrcIdMapSubset = new HashMap<Integer, SummedMagFreqDist>();
	    HashMap<Integer, ArrayList<Integer>> surfListForSrcIdMapSubset = new HashMap<Integer, ArrayList<Integer>>();
	    ArrayList<GeoJSONFaultSection> faultSectionDataSubset = new ArrayList<GeoJSONFaultSection>();
	    HashMap<Integer, Integer> newFltIndexMap = new HashMap<Integer, Integer>();
	    
	    // get a list of source rateWts and a list of fault sections used for the specified fault model branch
	    ArrayList<Integer> srcUsedList = new ArrayList<Integer> ();
	    HashMap<Integer, Double> srcRateWtMap = new HashMap<Integer, Double>();
	    for(int srcID:mfdForSrcIdMap.keySet()) {  // set defaults
	    	if(floaterSrcID_List.contains(srcID)) { // skip floater
	    		if (D) System.out.println("Skipping fault section/floater "+srcID);
	    		continue;
	    	}
	    	srcUsedList.add(srcID);
	    	srcRateWtMap.put(srcID, 1.0);
	    }
		if (fltModel == FaultModelEnum.PREFERRED) {
			for(int srcID:altSrcWtmap.keySet()) { // remove these
				if(floaterSrcID_List.contains(srcID)) continue;
				srcUsedList.remove(Integer.valueOf(srcID));
				srcRateWtMap.remove(Integer.valueOf(srcID));
			}
			for(int srcID:prefSrcWtmap.keySet()) { // change these rateWts
				if(floaterSrcID_List.contains(srcID)) continue;
				srcRateWtMap.remove(Integer.valueOf(srcID)); // remove
				srcRateWtMap.put(srcID, 1.0/prefSrcWtmap.get(srcID)); // put correct one back
			}
		}
		if (fltModel == FaultModelEnum.ALTERNATE) {
			for(int srcID:prefSrcWtmap.keySet()) { // remove these
				if(floaterSrcID_List.contains(srcID)) continue;
				srcUsedList.remove(Integer.valueOf(srcID));
				srcRateWtMap.remove(Integer.valueOf(srcID));
			}
			for(int srcID:altSrcWtmap.keySet()) { // change these rateWts
				if(floaterSrcID_List.contains(srcID)) continue;
				srcRateWtMap.remove(Integer.valueOf(srcID)); // remove
				srcRateWtMap.put(srcID, 1.0/altSrcWtmap.get(srcID)); // put correct one back
			}
		}
		
		if(D) {
			System.out.println("srcUsedList.size()="+srcUsedList.size()+"\t"+srcRateWtMap.size());
			System.out.println("srcRateWtMap:");
			for(int srcID:srcRateWtMap.keySet()) {
				System.out.println("\t"+srcID+"\t"+srcRateWtMap.get(srcID));
			}

		}
		
		// Now make a list of used fault sections
	    ArrayList<Integer> usedFaultSectID_List = new ArrayList<Integer>();
	    for(int srcID:srcUsedList) {
	    	for(int fltID:srcFltSectsMap.get(srcID))
	    		if(!usedFaultSectID_List.contains(fltID))
	    			usedFaultSectID_List.add(fltID);
	    }
	    
	    // get a new instance of fault sections list so we can override the IDs
	    ArrayList<GeoJSONFaultSection> duplicateFaultSectionList = getFaultSectionList(nshmModelDirPath);

	    // make list of fault sections
	    int newFltIndex=0;
	    for(GeoJSONFaultSection fltSection:duplicateFaultSectionList) {
	    	int sectID = fltSection.getSectionId();
	    	if(usedFaultSectID_List.contains(sectID)) { 
		    	fltSection.setParentSectionId(sectID);
		    	fltSection.setParentSectionName(parSecNameFromtID_Map.get(sectID));
		    	fltSection.setSectionId(newFltIndex);
		    	faultSectionDataSubset.add(fltSection);
		    	newFltIndexMap.put(sectID, newFltIndex);
		    	newFltIndex+=1;
	    	}
	    }
	    
	    // now make HashMaps needed sources and translated ids:
	    for(int srcID:srcUsedList) {
	    	int[] oldSectForSrcArray = srcFltSectsMap.get(srcID);
	    	ArrayList<Integer> newSectForSrcList = new ArrayList<Integer>();
	    	for( int oldSectID:oldSectForSrcArray)
	    		newSectForSrcList.add(newFltIndexMap.get(oldSectID));
	    	surfListForSrcIdMapSubset.put(srcID, newSectForSrcList);
	    	SummedMagFreqDist mfd = mfdForSrcIdMap.get(srcID);
	    	mfd.scale(srcRateWtMap.get(srcID)); // this permanently changes the mfd
	    	mfdForSrcIdMapSubset.put(srcID, mfd);
	    }
	    
	    FaultSystemSolution bigFSS = getFaultSystemSolution(mfdForSrcIdMapSubset, 
	    		surfListForSrcIdMapSubset, faultSectionDataSubset);
	    
	    
	    if(D) { // test participation mfds for each fault section
	    	
	    	// for FSS:
	    	double fssTotalMomentRate = 0;
	    	SummedMagFreqDist fssTotalMFD = getBlankMFD();
	    	HashMap<Integer, SummedMagFreqDist> sectMfdMapFSS = new HashMap<Integer, SummedMagFreqDist>();
	    	for(int rup=0;rup<bigFSS.getRupSet().getNumRuptures();rup++) {
	    		double rate = bigFSS.getRateForRup(rup);
	    		double mag = bigFSS.getRupSet().getMagForRup(rup);
	    		fssTotalMomentRate += MagUtils.magToMoment(mag)*rate;
	    		int iMag = fssTotalMFD.getClosestXIndex(mag);
	    		List<Integer> sectsForRupList = bigFSS.getRupSet().getSectionsIndicesForRup(rup);
	    		for(int i:sectsForRupList) {
	    			if(!sectMfdMapFSS.keySet().contains(i)) {
	    				sectMfdMapFSS.put(i, getBlankMFD());
	    			}
	    			sectMfdMapFSS.get(i).add(iMag, rate);
	    		}
	    		fssTotalMFD.add(iMag, rate);
	    	}
	    	
	    	// From ERF:
	    	double fssTotalMomentRate2 = 0;
	    	double origTotalMomentRate = 0;
		    double aveDeltaMag=0;
		    double numMag =0;
	    	SummedMagFreqDist erfTotalMFD = getBlankMFD();
	    	HashMap<Integer, SummedMagFreqDist> sectMfdMapERF = new HashMap<Integer, SummedMagFreqDist>();
	    	for(int id:newFltIndexMap.keySet()) {
	    		sectMfdMapERF.put(id, getBlankMFD());
	    	}
		    for(int s=0;s<erf.getNumSources();s++) {
		    	NshmSource src = (NshmSource)erf.getSource(s);
		    	if(srcZoneID_List.contains(src.getNSHM_ID())) 
		    		continue;	
		    	Integer srcID = src.getNSHM_ID();
		    	// skip if not used
		    	if(!srcUsedList.contains(srcID)) { 
		    		continue;
		    	}
		    	for(int r=0;r<src.getNumRuptures();r++) {
		    		double mag = src.getRupture(r).getMag();
		    		int iMag = erfTotalMFD.getClosestXIndex(mag);
		    		double rate = src.getRupture(r).getProbability()*srcRateWtMap.get(srcID);  // rate approx equal to prob
		    		erfTotalMFD.add(iMag, rate); // this requires the exact x value (no tolerance)
		    		origTotalMomentRate += MagUtils.magToMoment(mag)*rate;
		    		fssTotalMomentRate2 += MagUtils.magToMoment(erfTotalMFD.getX(iMag))*rate;
		    		aveDeltaMag += erfTotalMFD.getX(iMag)-mag;
		    		numMag += 1;
			    	for(int sectID:srcFltSectsMap.get(srcID)) {
				    	SummedMagFreqDist mfd = sectMfdMapERF.get(sectID);
				    	mfd.add(iMag, rate);
			    	}
		    	}
		    }
		    
		    // compare moments
		    double tempRatio = fssTotalMomentRate/origTotalMomentRate;
		    double testRatio = fssTotalMomentRate/fssTotalMomentRate2;
		    aveDeltaMag /= numMag;
		    if(D) System.out.println("FSS:\n\tfssTotalMomentRate="+(float)fssTotalMomentRate+
		    		"\n\torigTotalMomentRate="+(float)origTotalMomentRate+"\n\tratio="+tempRatio+
		    		"\n\ttestRatio="+(float)testRatio+"\n\taveDeltaMag="+(float)aveDeltaMag);
		    
		    // compare
		    for(int i=0;i<fssTotalMFD.size();i++) {
		    	double val1 = fssTotalMFD.getY(i);
		    	double val2 = erfTotalMFD.getY(i);
		    	if(val1 == 0.0) {
		    		if(val2 != 0)
		    			throw new RuntimeException("PROBLEM: zero in one but not the other at index "+i+"\n"+fssTotalMFD+"\n"+erfTotalMFD);
		    	}
		    	else {
		    		double ratio = val1/val2;
		    		if(ratio<0.99 || ratio > 1.01) {
		    			throw new RuntimeException("PROBLEM: >1% difference at index "+i+"\n"+fssTotalMFD+"\n"+erfTotalMFD);
		    		}
		    	}
		    }
		    
		    for(int sectID:sectMfdMapERF.keySet()) {
		    	SummedMagFreqDist mfd1 = sectMfdMapERF.get(sectID);
		    	SummedMagFreqDist mfd2 = sectMfdMapFSS.get(newFltIndexMap.get(sectID));
			    for(int i=0;i<mfd1.size();i++) {
			    	double val1 = mfd1.getY(i);
			    	double val2 = mfd2.getY(i);
			    	if(val1 == 0.0) {
			    		if(val2 != 0)
			    			throw new RuntimeException("PROBLEM: zero in one bu not the other at index "+i+"\n"+mfd1+"\n"+mfd2);
			    	}
			    	else {
			    		double ratio = val1/val2;
			    		if(ratio<0.99 || ratio > 1.01) {
			    			throw new RuntimeException("PROBLEM: >1% difference at index "+i+"\n"+mfd1+"\n"+mfd2);
			    		}
			    	}
			    }
		    }
		    System.out.println("ERF versus bigFSS tests passed!!");
	    }   

	    
		fssList.add(bigFSS);
		
		return fssList;
	}
	
	
	/**
	 * This creates a fault system solution for a floating rupture source/fault.  This assumes 
	 * the NSHM source ID is the same as the faultSection ID.  This does not account for any
	 * variations in the depth to top of rupture in the NSHM model.
	 * @param rateWt
	 * @param faultSection
	 * @param erf
	 * @return
	 */
	protected static FaultSystemSolution getFaultSystemSolution(double rateWt, 
			GeoJSONFaultSection faultSection, NshmErf erf) {
		
		if(D) System.out.println("working on floater FSS for "+faultSection.getName());

	
    	// find full rupture area (actually, this is the max rupture area in the ERf)
    	double fullRupAreaKmsq = 0;
	    for(int s=0;s<erf.getNumSources();s++) {
	    	NshmSource src = (NshmSource)erf.getSource(s);
	    	
// temp fix for Peters ID duplicates; exception test shows it's an 
// AK problem (never occurs from CEUS), at least with Peters old model; also looks fixed in AK 3.0.0 model?
if(src.getName().equals("Unnamed fault system source")) { 
//	throw new RuntimeException("got one here");
	continue;
}
	    	if(src.getNSHM_ID() == faultSection.getSectionId()) { 
		    	for(int r=0;r<src.getNumRuptures();r++) {
		    		double rupArea = src.getRupture(r).getRuptureSurface().getArea();
		    		if(fullRupAreaKmsq<rupArea)
		    			fullRupAreaKmsq=rupArea;
		    	}
	    	}
	    }
	    
	    // now compute full rup vs floater MFDs
	    double origTotalMoment = 0;
	    double fssTotalMoment = 0;
	    double aveDeltaMag=0;
	    double numMag =0;
	    ArrayList<Float> rupTopDepthList = new ArrayList<Float>();
    	SummedMagFreqDist mfd_full = getBlankMFD();
    	SummedMagFreqDist mfd_float = getBlankMFD();
	    for(int s=0;s<erf.getNumSources();s++) {
	    	NshmSource src = (NshmSource)erf.getSource(s);
	
// temp fix for Peters ID duplicates; exception test shows it's an 
// AK problem (never occurs from CEUS), at least with Peters old model; also looks fixed in AK 3.0.0 model?
if(src.getName().equals("Unnamed fault system source")) // temp fix for Peters ID duplicates
	continue;

	    	if(src.getNSHM_ID() == faultSection.getSectionId()) { 
		    	for(int r=0;r<src.getNumRuptures();r++) {
		    		double mag = src.getRupture(r).getMag();
		    		int iMag = mfd_full.getClosestXIndex(mag);
		    		double rate = src.getRupture(r).getProbability();  // rate approx equal to prob
		    		double rupAreaKmsq = src.getRupture(r).getRuptureSurface().getArea();
		    		origTotalMoment += MagUtils.magToMoment(mag)*rate;
		    		fssTotalMoment += MagUtils.magToMoment(mfd_full.getX(iMag))*rate;
		    		aveDeltaMag += mfd_full.getX(iMag)-mag;
		    		numMag+=1;
// System.out.println(mag+"\t"+rupArea);

		    		if(rupAreaKmsq > 0.99*fullRupAreaKmsq) {
		    			mfd_full.add(iMag, rate*rateWt);
		    		}
		    		else {
		    			mfd_float.add(iMag, rate*rateWt);
		    		}
		    		float rupTopDepth = (float)src.getRupture(r).getRuptureSurface().getAveRupTopDepth();
		    		if(!rupTopDepthList.contains(rupTopDepth))
		    			rupTopDepthList.add(rupTopDepth);
		    	}
	    	}
	    }
	    
// System.out.println(mfd_full+"\n\n"+mfd_float);

	    aveDeltaMag /= numMag;
	    double tempRatio = fssTotalMoment/origTotalMoment;
	    if(D) System.out.println("   FSS:\n\tfssTotalMoment="+(float)fssTotalMoment+
	    		"\n\torigTotalMoment="+(float)origTotalMoment+"\n\tratio="+tempRatio+
	    		"\n\taveDeltaMag="+(float)aveDeltaMag+"\n\tfullRupArea="+fullRupAreaKmsq);
	    
	    List<List<Integer>> sectionForRups = new ArrayList<>();
	    ArrayList<Double> magForRupList =new ArrayList<Double>();
	    ArrayList<Double> rateForRupList =new ArrayList<Double>();
	    ArrayList<Double> areaForRupList =new ArrayList<Double>();
	    ArrayList<Double> lengthForRupList =new ArrayList<Double>();
	    
		// subsection the fault section
		double ddw = faultSection.getOrigDownDipWidth();
		List<? extends FaultSection> subsectionList = faultSection.getSubSectionsList(ddw/4.0, 0);
		double subSectLen = subsectionList.get(0).getTraceLength();
		WC1994_MagLengthRelationship wcMagLength = new WC1994_MagLengthRelationship();
		double fullMag = wcMagLength.getMedianMag(faultSection.getTraceLength());
		double lengthOfBiggestFloater = wcMagLength.getMedianLength(mfd_float.getMaxMagWithNonZeroRate());
		if(lengthOfBiggestFloater>faultSection.getTraceLength()) {
			double maxFloatMag= mfd_float.getMaxMagWithNonZeroRate();
			throw new RuntimeException("lengthOfBiggestFloater>faultSection.getTraceLength() not supported"+
					"\nlengthOfBiggestFloater="+lengthOfBiggestFloater+"\nTraceLength="+faultSection.getTraceLength()+
					"\nmaxFloatMag="+maxFloatMag+ "\nFor "+faultSection.getName());
		}
		if(D) {
			System.out.println("\tParent section length: "+faultSection.getTraceLength());
			System.out.println("\tSub-section length: "+subSectLen);
			System.out.println("\tNumSub-sections: "+subsectionList.size());
			System.out.println("\tFull mag: "+fullMag);
			System.out.println("\tMin full rup mag: "+mfd_full.getMinMagWithNonZeroRate());
			System.out.println("\tMax full rup mag: "+mfd_full.getMaxMagWithNonZeroRate());
			System.out.println("\tMin floater mag: "+mfd_float.getMinMagWithNonZeroRate());
			System.out.println("\tMax floater mag: "+mfd_float.getMaxMagWithNonZeroRate());
			System.out.println("\tlengthOfBiggestFloater: "+lengthOfBiggestFloater);
			System.out.println("\trupTopDepths in ERF Not applied here): ");
			for(double rupTopDepth:rupTopDepthList)
				System.out.println("\t\t"+(float)rupTopDepth);
		}
		int rupIndex = 0;
		boolean fullFaultFloater = false;
		for(int i=0;i<mfd_float.size();i++) {
			double rate = mfd_float.getIncrRate(i);
			if(rate == 0.0)continue;
			double mag = mfd_float.getX(i);
			double length = wcMagLength.getMedianLength(mag);
			int numSect = (int)Math.round(length/subSectLen);
			if(numSect==subsectionList.size())
				fullFaultFloater = true;
			int numRups = subsectionList.size()-numSect+1;
			if (D) System.out.println("Floater for M "+(float)mag+"\t"+(float)length+"\t"+numSect+"\t"+numRups);
			int firstSect = 0;
			int lastSect = numSect-1;
			int testNumRups=0;
			while(lastSect < subsectionList.size()) {
				testNumRups+=1;
				if (D) System.out.println("\t"+rupIndex+"\t"+firstSect+"\t"+lastSect);
				ArrayList<Integer> sectIndexForRupList = new ArrayList<>();
				for(int s=firstSect; s<=lastSect; s++)
					sectIndexForRupList.add(s);
				sectionForRups.add(sectIndexForRupList);
			    magForRupList.add(mag);
			    rateForRupList.add(rate/numRups);
			    areaForRupList.add(subSectLen*numSect*faultSection.getOrigDownDipWidth());
			    lengthForRupList.add(subSectLen*numSect);

				firstSect+=1;
				lastSect+=1;
				rupIndex+=1;
			}
			if(testNumRups != numRups)
				throw new RuntimeException("testNumRups != numRups");
		}
		if (D) System.out.println("\n\tfullFaultFloater = "+fullFaultFloater+"\n");
		
		// add full fault ruptures, one for each mag in MFD (duplicate here if floater has full rup)
		double fullLength = faultSection.getTraceLength();
		double fullArea = faultSection.getTraceLength()*faultSection.getOrigDownDipWidth();
		ArrayList<Integer> sectIndexForFullRupList = new ArrayList<>();
		for(int s=0; s<subsectionList.size(); s++)
			sectIndexForFullRupList.add(s);
		for(int i=0;i<mfd_full.size();i++) {
			double rate = mfd_full.getIncrRate(i);
			if(rate == 0.0)continue;
			double mag = mfd_full.getX(i);
			sectionForRups.add(sectIndexForFullRupList);
		    magForRupList.add(mag);
		    rateForRupList.add(rate);
		    areaForRupList.add(fullArea);
		    lengthForRupList.add(fullLength);
		}


		int numRups = sectionForRups.size();
		double[] mags = new double[numRups];
		double[] rakes = new double[numRups];  // defaults to zero here
		double[] rupAreas = new double[numRups];
		double[] rupLengths = new double[numRups];
		double[] rupRates = new double[numRups];
		for(int r=0;r<numRups;r++) {
			mags[r] = magForRupList.get(r);
			rupAreas[r] = areaForRupList.get(r)*1e6; // convert from km to m squared
			rupLengths[r] = lengthForRupList.get(r)*1e3; // convert from km to m
			rupRates[r] = rateForRupList.get(r);
			rakes[r] = faultSection.getAveRake();
		}

	    FaultSystemRupSet rupSet = new FaultSystemRupSet(
	    		subsectionList,
				sectionForRups,
				mags,
				rakes,
				rupAreas,
				rupLengths); 
	    
	    FaultSystemSolution fss= new FaultSystemSolution(rupSet, rupRates);
	    fss.setInfoString("FSS for floater source ID="+faultSection+"; Name: "+faultSection.getName());
	    return fss;
	}

	
	
	
	/**
	 * This assumes none of the sources have floaters (no partial fault-section ruptures).  
	 * This also assumes fault sections have incremental IDs starting at 0.  Rupture lengths 
	 * and areas are computed from fault sections provided.
	 * @param mfdForSrcIdMap
	 * @param surfListForSrcIdMap
	 * @return
	 */
	public static FaultSystemSolution getFaultSystemSolution(HashMap<Integer, SummedMagFreqDist> mfdForSrcIdMap, 
			HashMap<Integer, ArrayList<Integer>> surfListForSrcIdMap, 
			ArrayList<GeoJSONFaultSection> faultSectionData) {

	    // compute number of ruptures from nonzero rates in MFDs
	    int numRups = 0;
	    for(int srcID:mfdForSrcIdMap.keySet()) {
	    	SummedMagFreqDist mfd = mfdForSrcIdMap.get(srcID);
	    	for(int i=0;i<mfd.size();i++)
	    		if(mfd.getY(i) >0)
	    			numRups += 1;		
	    }
	    
		double[] mags = new double[numRups];
		double[] rakes = new double[numRups];  // defaults to zero here
		double[] rupAreas = new double[numRups];
		double[] rupLengths = new double[numRups];
		double[] rupRates = new double[numRups];
	    List<List<Integer>> sectionForRups = new ArrayList<>();
	    
	    if(D) System.out.println("numRups="+numRups);
	    
	    int r = 0;
	    for(int srcID:mfdForSrcIdMap.keySet()) {
	    	ArrayList<Integer> sectForRupList = surfListForSrcIdMap.get(srcID);
	    	double lengthKm = 0;
	    	double area = 0;
	    	double rake = faultSectionData.get(sectForRupList.get(0)).getAveRake();  // get rake of first section
	    	for(int id:sectForRupList) {
	    		lengthKm += faultSectionData.get(id).getTraceLength();
	    		area += faultSectionData.get(id).getArea(false);
	    		// check that rake is constant across sections
	    		double tempRake = faultSectionData.get(id).getAveRake();
	    		if(tempRake != rake)
	    			if (D) System.out.println("WARNING: Rake changes among sections for srcID="+srcID
	    					+"\n\t"+tempRake+" for "+faultSectionData.get(id).getName() 
	    					+ "\n\t"+rake+" for "+faultSectionData.get(sectForRupList.get(0)).getName());
	    	}
	    	SummedMagFreqDist mfd = mfdForSrcIdMap.get(srcID);
	    	for(int i=0;i<mfd.size();i++) {
	    		double rate = mfd.getY(i);
	    		if(rate >0) {
	    			rupRates[r] = rate;
	    			mags[r] = mfd.getX(i);
	    			rupAreas[r] = area;
	    			rakes[r] = rake;
	    			rupLengths[r] = lengthKm*1e3; // convert to meters
	    			sectionForRups.add(sectForRupList);
	    			r += 1;		
	    		}
	    	}
	    }

	    
	    FaultSystemRupSet rupSet = new FaultSystemRupSet(
				faultSectionData,
				sectionForRups,
				mags,
				rakes,
				rupAreas,
				rupLengths); 
	    
	    FaultSystemSolution fss= new FaultSystemSolution(rupSet, rupRates);
	    fss.setInfoString("Big 2023 CEUS Fault System Solution");
	    return fss;
	}
	
	
	
	 /**
	  * This returns the ERF with timespan duration set to 1.0
	  * @param nshmModelDirPath
	  * @return
	  */
	private static NshmErf getNshmERF(String nshmModelDirPath) {
	    Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.STABLE_SHALLOW);
	    NshmErf erf = new NshmErf(Path.of(nshmModelDirPath), trts, IncludeBackgroundOption.EXCLUDE);
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
	    return erf;
	}
	
	
	private static void testForPeter(String nshmModelDirPath) {
	    Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.STABLE_SHALLOW);
	    NshmErf erf = new NshmErf(Path.of(nshmModelDirPath), trts, IncludeBackgroundOption.EXCLUDE);
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
	    int s=232; // Cheraw USGS source
	    NshmSource src = (NshmSource)erf.getSource(s); 
    	System.out.println(s+"\t"+src.getName()+"\t"+src.getNSHM_ID()+"\t"+src.getNumRuptures());
	    for(int r=0;r<src.getNumRuptures();r++) {
	    	double area = src.getRupture(r).getRuptureSurface().getArea();
	    	double ddw = src.getRupture(r).getRuptureSurface().getAveWidth();
	    	double tor = src.getRupture(r).getRuptureSurface().getAveRupTopDepth();
	    	double len = src.getRupture(r).getRuptureSurface().getAveLength();
		    System.out.println("\t"+r+"\t"+src.getRupture(r).getMag()+"\t"+(float)area+"\t"+(float)ddw+"\t"+(float)tor+"\t"+(float)len);
	    }
	}
	
	/**
	 * This is for parsing Peter's cluster-set.json files
	 * @param filePath
	 */
	public static void parseClusterSetFile(String filePath, HashMap<Integer,int[]> srcFltSectsMap) {
		 String ID = "id";
		 String NAME = "name";
		 String SECTIONS = "sections";
		 String RUPTURE_SETS = "rupture-sets";

		Path path = Paths.get(filePath);
		
	    JsonObject obj=null;
	    try (BufferedReader br = Files.newBufferedReader(path)) {
	    	obj = JsonParser.parseReader(br).getAsJsonObject();
	    } catch (IOException ioe) {
	    	throw new RuntimeException(ioe);
	    }
	    
//	    System.out.println("Cluster set: "+obj.get(ID).getAsInt()+"\t"+obj.get(NAME).getAsString());
	    
	    JsonArray rupSets = obj.get(RUPTURE_SETS).getAsJsonArray();
	    for (JsonElement rupSet : rupSets) {
	    	JsonObject rupSetObj = rupSet.getAsJsonObject();
		    int srcID = rupSetObj.get(ID).getAsInt();
		    int[] fltSectionIDs;

		    if (rupSetObj.has(SECTIONS)) {
			      fltSectionIDs = GSON.fromJson(rupSetObj.get(SECTIONS), int[].class);
		    } else {
		    	fltSectionIDs = new int[]{srcID}; // section ID same as src ID
		    }
		    
		    if(!srcFltSectsMap.keySet().contains(srcID)) { // no need for redundancies because src IDs have unique rup surface
		    	srcFltSectsMap.put(srcID, fltSectionIDs);
		    }

//				System.out.println("\t"+srcID);
	    }
	}

	
	
	/**
	 * This is for parsing Peter's rupture-set.json files
	 * @param filePath
	 */
	public static void parseRuptureSetFile(String filePath,HashMap<Integer,int[]> srcFltSectsMap) {
		
		 String ID = "id";
		 String NAME = "name";
		 String SECTIONS = "sections";

		Path path = Paths.get(filePath);
		
	    JsonObject obj=null;
	    try (BufferedReader br = Files.newBufferedReader(path)) {
	    	obj = JsonParser.parseReader(br).getAsJsonObject();
	    } catch (IOException ioe) {
	    	throw new RuntimeException(ioe);
	    }
	    
//	    System.out.println("Rupture set: "+obj.get(ID).getAsInt()+"\t"+obj.get(NAME).getAsString());
	    
	    int srcID = obj.get(ID).getAsInt();
	    
	    int[] fltSectionIDs;
	    if (obj.has(SECTIONS)) {
	      fltSectionIDs = GSON.fromJson(obj.get(SECTIONS), int[].class);
	    } else {
	    	fltSectionIDs = new int[]{srcID}; // section ID same as src ID
	    }
//	    for(int sect:fltSectionIDs)
//	    	System.out.println("\t"+sect);
//if(fltSectionIDs.length==1)
//	System.out.println("ONLY ONE SUBSURF: "+srcID+"\t"+fltSectionIDs[0]+"\t"+obj.get(NAME).getAsString());
	    
	    if(!srcFltSectsMap.keySet().contains(srcID)) { // no need for redundancies because src IDs have unique rup surface
	    	srcFltSectsMap.put(srcID, fltSectionIDs);
	    }

	}
	
	/**
	 * I found these files by listing for "rupture-set.json" and "cluster-set.json" at various directory depths
	 * @param srcIDsList
	 * @param srcFltSectsList
	 */
	private static void getSrcIDsAndFaultSectionsLists(HashMap<Integer,int[]> srcFltSectsMap, String nshmModelDirPath) {
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/CO/Cheraw/usgs/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/Commerce/2-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/Commerce/3-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/OK/Meers/usgs/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (North)/1-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (North)/2-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/CO/Cheraw/sscn/recurrence-rate/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/OK/Meers/sscn/cluster-in/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/OK/Meers/sscn/cluster-out/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/crittenden-co/2-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/crittenden-co/3-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/crittenden-co/4-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/meeman-shelby/2-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/meeman-shelby/3-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/TN/Eastern Rift Margin (South)/meeman-shelby/4-eq/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/CO/Cheraw/sscn/slip-rate/full-rupture/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/CO/Cheraw/sscn/slip-rate/partial-rupture/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-out/reelfoot-extended/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-out/reelfoot-short/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/center/cluster-out/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/east/cluster-out/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/mid-east/cluster-out/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/mid-west/cluster-out/rupture-set.json", srcFltSectsMap);
		parseRuptureSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/west/cluster-out/rupture-set.json", srcFltSectsMap);

		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsaxn-rftl-nmnl/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsaxn-rftl-nmns/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsaxn-rfts-nmnl/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsaxn-rfts-nmns/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsbl-rftl-nmnl/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsbl-rftl-nmns/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsbl-rfts-nmnl/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/sscn/cluster-in/axsbl-rfts-nmns/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/center/cluster-in/all/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/center/cluster-in/center-south/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/east/cluster-in/all/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/east/cluster-in/center-south/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/mid-east/cluster-in/all/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/mid-east/cluster-in/center-south/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/mid-west/cluster-in/all/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/mid-west/cluster-in/center-south/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/west/cluster-in/all/cluster-set.json", srcFltSectsMap);
		parseClusterSetFile(nshmModelDirPath+"stable-crust/fault/MO/New Madrid/usgs/west/cluster-in/center-south/cluster-set.json", srcFltSectsMap);
	}
	
	private static SummedMagFreqDist getBlankMFD() {
		return new SummedMagFreqDist(5.05,80,0.05);
	}

	
	public static void main(String[] args) {
		
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.1.2/";
		// previous version of above won't work because hard-coded files changed

		ArrayList<FaultSystemSolution> fssList = getFaultSystemSolutionList(nshmModelDirPath,FaultModelEnum.PREFERRED);
		for(FaultSystemSolution fss:fssList) {
			int s=0;
			for(FaultSection sect:fss.getRupSet().getFaultSectionDataList()) {
				String magString = "Mags: ";
				for(int r: fss.getRupSet().getRupturesForSection(s))
					magString += (float)fss.getRupSet().getMagForRup(r)+", ";
				System.out.println(sect.getSectionId()+"\t"+sect.getParentSectionId()+"\t"+sect.getName()+"\t"+magString);
				s+=1;
			}
		}
		
//	    WC1994_MagLengthRelationship wcMagLength = new WC1994_MagLengthRelationship();
//	    double mag = 7.15;
//	    System.out.println(wcMagLength.getMedianLength(mag));
//	    System.out.println(wcMagLength.getMedianLength(mag,0.0));
//	    System.out.println(wcMagLength.getMedianLength(mag,-90));
//	    System.out.println(wcMagLength.getMedianLength(mag,90));

		
		
//	    // write attributes of something
//	    try {
//	    	File file = new File("junkRightHere");
//	    	DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
//	    	String lineString = "bla bla";
//	    	out.writeChars(lineString);
//	    	out.close();
//	    } catch (IOException e) {
//	    	e.printStackTrace();
//	    }


		
	}

}
