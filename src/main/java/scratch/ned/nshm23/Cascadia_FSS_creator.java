package scratch.ned.nshm23;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.opensha.commons.eq.MagUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.ApproxEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.fault.surface.ApproxGriddedSurface;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSource;
import scratch.ned.nshm23.CEUS_FSS_creator.FaultModelEnum;

public class Cascadia_FSS_creator {
	
	private final static Boolean D = true;
	
	 public enum FaultModelEnum {
		 ALL(1.0), 
		 TOP(0.2), 
		 MIDDLE(0.5), 
		 BOTTOM(0.3) ;
		 private double weight;
		 FaultModelEnum(double weight) {
			 this.weight=weight;
		 }
		 double getWeight() {
			 return weight;
		 }
	 }

	
	 /**
	  * This returns the ERF with timespan duration set to 1.0
	  * @param nshmModelDirPath
	  * @return
	  */
	private static NshmErf getNshmERF(String nshmModelDirPath) {
	    Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.SUBDUCTION_INTERFACE);
	    NshmErf erf = new NshmErf(Path.of(nshmModelDirPath), trts, IncludeBackgroundOption.EXCLUDE);
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
	    return erf;
	}
	

	private static ArrayList<GeoJSONFaultSection> getFaultSectionList(String nshmModelDirPath, FaultModelEnum fltMod) {
		ArrayList<GeoJSONFaultSection> list = new ArrayList<GeoJSONFaultSection>();

		String[] nameArray = {"top","middle","bottom"}; // default case for "ALL"
		switch (fltMod) {
		case ALL:
			break;  
		case TOP:
			nameArray = new String[1];
			nameArray[0] = "top";
			break;  
		case MIDDLE:
			nameArray = new String[1];
			nameArray[0] = "middle";
			break;
		case BOTTOM:
			nameArray = new String[1];
			nameArray[0] = "bottom";
			break;
		default:
			throw new RuntimeException("Problem with fault model");
		}
		
		for(String fltModName:nameArray) {
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-1 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-2 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-3 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-4 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-5 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-6 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-7 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 2-1 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 2-2 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 2-3 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 3-1 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 3-2 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 3-3 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-1 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-2 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-3 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-4 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-5 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-6 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-7 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-8 ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 4-9 ("+fltModName+").geojson"));
		}
		return list;
	}

	
	
	/**
	 * I found these files by listing for "rupture-set.json" and "cluster-set.json" at various directory depths
	 * @param srcIDsList
	 * @param srcFltSectsList
	 */
	private static void getSrcIDsAndFaultSectionsLists(HashMap<Integer,int[]> srcFltSectsMap, String nshmModelDirPath, FaultModelEnum fltMod) {
		
		if(D) System.out.println(fltMod);
		String[] nameArray = {"top","middle","bottom"}; // default case for "ALL"
		switch (fltMod) {
		case ALL:
			break;  
		case TOP:
			nameArray = new String[1];
			nameArray[0] = "top";
			break;  
		case MIDDLE:
			nameArray = new String[1];
			nameArray[0] = "middle";
			break;
		case BOTTOM:
			nameArray = new String[1];
			nameArray[0] = "bottom";
			break;
		default:
			throw new RuntimeException("Problem with fault model");
		}
		
		for(String name:nameArray) {
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/full-rupture/cluster-out/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA12/B/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA12/C/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA12/D/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA12/northern/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA17/B/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA17/C'/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA17/C/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA17/D/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA17/E/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/segmented/GEA17/F/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/unsegmented/GEA12-A/scaled/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/partial-rupture/unsegmented/GEA12-B/scaled/rupture-set.json", srcFltSectsMap);

			CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/full-rupture/cluster-in/cluster-7a/cluster-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/full-rupture/cluster-in/cluster-7b/cluster-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/"+name+"/full-rupture/cluster-in/cluster-8/cluster-set.json", srcFltSectsMap);

		}
		
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/full-rupture/cluster-out/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA12/B/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA12/C/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA12/D/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA12/northern/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA17/B/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA17/C'/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA17/C/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA17/D/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA17/E/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/segmented/GEA17/F/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/unsegmented/GEA12-A/scaled/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/partial-rupture/unsegmented/GEA12-B/scaled/rupture-set.json", srcFltSectsMap);
//
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/full-rupture/cluster-out/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA12/B/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA12/C/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA12/D/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA12/northern/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA17/B/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA17/C'/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA17/C/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA17/D/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA17/E/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/segmented/GEA17/F/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/unsegmented/GEA12-A/scaled/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/partial-rupture/unsegmented/GEA12-B/scaled/rupture-set.json", srcFltSectsMap);
//
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/full-rupture/cluster-out/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA12/B/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA12/C/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA12/D/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA12/northern/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA17/B/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA17/C'/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA17/C/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA17/D/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA17/E/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/segmented/GEA17/F/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/unsegmented/GEA12-A/scaled/rupture-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/partial-rupture/unsegmented/GEA12-B/scaled/rupture-set.json", srcFltSectsMap);
//
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/full-rupture/cluster-in/cluster-7a/cluster-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/full-rupture/cluster-in/cluster-7b/cluster-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/bottom/full-rupture/cluster-in/cluster-8/cluster-set.json", srcFltSectsMap);
//
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/full-rupture/cluster-in/cluster-7a/cluster-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/full-rupture/cluster-in/cluster-7b/cluster-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/middle/full-rupture/cluster-in/cluster-8/cluster-set.json", srcFltSectsMap);
//
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/full-rupture/cluster-in/cluster-7a/cluster-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/full-rupture/cluster-in/cluster-7b/cluster-set.json", srcFltSectsMap);
//		CEUS_FSS_creator.parseClusterSetFile(nshmModelDirPath+"subduction/interface/Cascadia/top/full-rupture/cluster-in/cluster-8/cluster-set.json", srcFltSectsMap);
	}


	public static FaultSystemSolution getFaultSystemSolution(String nshmModelDirPath, FaultModelEnum fltModel) {
		
		// rate weight for specified fault model branch
		double rateWt = 1.0/fltModel.getWeight();
		if(D)System.out.println("fltModel="+fltModel+"\nrateWt="+(float)rateWt);
		
	    // get fault section list for given fault model
	//	HashMap<Integer,GeoJSONFaultSection> faultSectionMap;
	    ArrayList<GeoJSONFaultSection> faultSectionList = getFaultSectionList(nshmModelDirPath, fltModel);
		// make parSectID_List & write attributes
		ArrayList<Integer> parSectID_List = new ArrayList<Integer>(); // this is NSHM ID for each parent section
		if(D) System.out.println("index\tsectID\trake");
		for(int s=0;s<faultSectionList.size();s++) {
			GeoJSONFaultSection sect = faultSectionList.get(s);
			if(!parSectID_List.contains(sect.getSectionId())) {
				parSectID_List.add(sect.getSectionId());
				if(D) System.out.println(s+
						"\tID="+sect.getSectionId()+
						"\tRake="+sect.getAveRake()+
						"\tDip="+(float)sect.getAveDip()+
						"\tLowDep="+(float)sect.getAveLowerDepth()+
						"\tUpperDep="+(float)sect.getOrigAveUpperDepth()+
						"\tDDW="+(float)sect.getOrigDownDipWidth()+
						"\tLength"+(float)sect.getTraceLength());
			}
			else
				throw new RuntimeException("section IDs are not unique; duplicate: "+sect.getSectionId());
		}
		if (D)System.out.println("parSectID_List.size()="+parSectID_List.size());
		
		// Read from Peter's rupture-set.json and and cluster-set.json files
		HashMap<Integer,int[]> srcFltSectsMap = new HashMap<Integer,int[]>(); // the fault section used by each source (same order as above)
		getSrcIDsAndFaultSectionsLists(srcFltSectsMap, nshmModelDirPath, fltModel);
		Set<Integer> srcIDsList = srcFltSectsMap.keySet();  // a list of all the source IDs (no duplicates)

		// some tests
		if(D) System.out.println("number of unique source IDs: "+srcIDsList.size());
		// make sure all sections are used and none are missing (with respect to Peter's files)
		ArrayList<Integer> testParSectID_List = new ArrayList<Integer>(); // fill up with the par IDs from all the sources
		for(int srcID:srcIDsList) {
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

//		System.exit(0);
		
		
		
		NshmErf erf = getNshmERF(nshmModelDirPath);
		System.out.println("erf.getNumSources() = "+erf.getNumSources());
//		int numPtSrc=0;
		ArrayList<Integer> testSrcIDsList = new ArrayList<Integer>();
		for(int s=0;s<erf.getNumSources();s++) {
			NshmSource src = (NshmSource)erf.getSource(s);
			int srcID = src.getNSHM_ID();
			// Following no longer needed because there are no other sources in ERF
//			if(CEUS_FSS_creator.isPointSource(src))
//				numPtSrc+=1;
//			else {
				// skip if not on the selected branch
				if(!srcIDsList.contains(srcID))
					continue;
				System.out.println("\t"+s+"\t"+srcID+"\t"+src.getNumRuptures()+"\t"+src.getName()+"\t"+src.getRupture(0).getClass().toString());
				if(!testSrcIDsList.contains(srcID))
					testSrcIDsList.add(srcID);
				else
					System.out.println("Duplicate srcID = "+srcID);
//			}
		}
		System.out.println("srcID_List.size() = "+testSrcIDsList.size());
//		System.out.println("numPtSrc = "+numPtSrc);
		
	    // run more tests; first make sure ERF sources are consistent with Peter's files
    	if(testSrcIDsList.size() != srcIDsList.size())
    		throw new RuntimeException("testSrcIDsList.size() != srcIDsList.size()");
    	for(int id : testSrcIDsList)
    		if(!srcIDsList.contains(id))
    			throw new RuntimeException("srcIDsList does not contain: "+id);
    	if(D)System.out.println("srcIDsList passed tests");
    	
    	
    	// Compute MFDs for each source
    	HashMap<Integer, SummedMagFreqDist> mfdForSrcIdMap = new HashMap<Integer, SummedMagFreqDist>();
    	HashMap<Integer, Double> rakeForSrcIdMap = new HashMap<Integer, Double>();
    	HashMap<Integer, String> nameForSrcIdMap = new HashMap<Integer, String>();
    	for(int id:srcIDsList) {
        	mfdForSrcIdMap.put(id, getBlankMFD());
    	}
	    for(int s=0;s<erf.getNumSources();s++) {
	    	NshmSource src = (NshmSource)erf.getSource(s);
//	    	if(CEUS_FSS_creator.isPointSource(src)) 
//	    		continue;	
	    	Integer srcID = src.getNSHM_ID();
	    	// skip if not on the flt mod branch
			if(!srcIDsList.contains(srcID)) 
				continue;
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
	    		double rate = src.getRupture(r).getMeanAnnualRate(1.0)*rateWt;
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
	    	System.out.println("Min and max mag, and rate, for each source ID");
		    for(int srcID:mfdForSrcIdMap.keySet()) {
		    	SummedMagFreqDist mfd = mfdForSrcIdMap.get(srcID);
		    	System.out.println("\t"+srcID+"\t"+(float)mfd.getMinMagWithNonZeroRate()+"\t"+(float)mfd.getMaxMagWithNonZeroRate()+"\trake="+
		    	rakeForSrcIdMap.get(srcID)+"\t"+nameForSrcIdMap.get(srcID));
		    }
	    }
	    
	    // change fault IDs to start from 0
	    HashMap<Integer, Integer> newFltIndexMap = new HashMap<Integer, Integer>();
	    int newFltIndex=0;
	    for(GeoJSONFaultSection fltSection:faultSectionList) {
	    	int sectID = fltSection.getSectionId();
		    fltSection.setParentSectionId(sectID);
		    fltSection.setSectionId(newFltIndex);
		    newFltIndexMap.put(sectID, newFltIndex);
		    newFltIndex+=1;
	    }
	    
	    // now make HashMap with translated fault section IDs:
	    HashMap<Integer, ArrayList<Integer>> translatedSurfListForSrcIdMap = new HashMap<Integer, ArrayList<Integer>>();
	    for(int srcID:srcIDsList) {
	    	int[] oldSectForSrcArray = srcFltSectsMap.get(srcID);
	    	ArrayList<Integer> newSectForSrcList = new ArrayList<Integer>();
	    	for( int oldSectID:oldSectForSrcArray)
	    		newSectForSrcList.add(newFltIndexMap.get(oldSectID));
	    	translatedSurfListForSrcIdMap.put(srcID, newSectForSrcList);
	    }
	    
	    FaultSystemSolution fss = CEUS_FSS_creator.getFaultSystemSolution(mfdForSrcIdMap, 
	    		translatedSurfListForSrcIdMap, faultSectionList);
	    
	    if(D) { // test participation mfds for each fault section
	    	
	    	// for FSS:
	    	double fssTotalMomentRate = 0;
	    	SummedMagFreqDist fssTotalMFD = getBlankMFD();
	    	HashMap<Integer, SummedMagFreqDist> sectMfdMapFSS = new HashMap<Integer, SummedMagFreqDist>();
	    	for(int rup=0;rup<fss.getRupSet().getNumRuptures();rup++) {
	    		double rate = fss.getRateForRup(rup);
	    		double mag = fss.getRupSet().getMagForRup(rup);
	    		fssTotalMomentRate += MagUtils.magToMoment(mag)*rate;
	    		int iMag = fssTotalMFD.getClosestXIndex(mag);
	    		List<Integer> sectsForRupList = fss.getRupSet().getSectionsIndicesForRup(rup);
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
//		    	if(CEUS_FSS_creator.isPointSource(src)) 
//		    		continue;	
		    	Integer srcID = src.getNSHM_ID();
		    	// skip if not on the flt mod branch
				if(!srcIDsList.contains(srcID)) 
					continue;
		    	for(int r=0;r<src.getNumRuptures();r++) {
		    		double mag = src.getRupture(r).getMag();
		    		int iMag = erfTotalMFD.getClosestXIndex(mag);
		    		double rate = src.getRupture(r).getMeanAnnualRate(1.0)*rateWt;  // rate approx equal to prob
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
		    if(D) {
		    	System.out.println("FSS:\n\tfssTotalMomentRate="+(float)fssTotalMomentRate+
		    		"\n\torigTotalMomentRate="+(float)origTotalMomentRate+"\n\tratio="+tempRatio+
		    		"\n\ttestRatio="+(float)testRatio+"\n\taveDeltaMag="+(float)aveDeltaMag);
		    	System.out.println("FSS Cumulative MFD (mag, rate, RI):");
		    	for(double mag = fssTotalMFD.getMinMagWithNonZeroRate(); mag <= fssTotalMFD.getMaxMagWithNonZeroRate(); mag+=0.05)
			    	System.out.println("\t"+(float)mag+"\t"+(float)fssTotalMFD.getCumRate(mag)+"\t"+(float)(1/fssTotalMFD.getCumRate(mag)));

		    }
		    
		    
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
			    			throw new RuntimeException("PROBLEM: zero in one but not the other at index "+i+"\n"+mfd1+"\n"+mfd2);
			    	}
			    	else {
			    		double ratio = val1/val2;
			    		if(ratio<0.99 || ratio > 1.01) {
			    			throw new RuntimeException("PROBLEM: >1% difference at index "+i+"\n"+mfd1+"\n"+mfd2);
			    		}
			    	}
			    }
		    }
		    System.out.println("ERF versus fss tests passed!!");
	    }   


		return fss;
	}
	


	public static void main(String[] args) {
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.1.2/";
		
		
		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.ALL);
		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.TOP);
		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.MIDDLE);
		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.BOTTOM);
		
//		GeoJSONFaultSection sect = CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Cascadia/features/Cascadia 1-1 (middle).geojson");
//		System.out.println("OrigAveUpperDepth = "+sect.getOrigAveUpperDepth());
//		System.out.println("AveLowerDepth() = "+sect.getAveLowerDepth());
//		System.out.println("Trace:\n\n"+sect.getFaultTrace());
//		System.exit(0);

		// THIS IS HOW PETER CREATES THE SURFACE (FROM gov.usgs.earthquake.nshmp.model.InterfaceSource.InterfaceSource(Builder))
//	    upperTrace = feature.traces.get(0);
//	    upperTraceLength = upperTrace.length();
//	    lowerTrace = feature.traces.get(feature.traces.size() - 1);
//	    lowerTraceLength = lowerTrace.length();
//
//	    surface = new ApproxGriddedSurface(upperTrace, lowerTrace, config.surfaceSpacing);
		
//		THIS IS OURS (same constructor arguments):
//		ApproxEvenlyGriddedSurface surf;

		
		



	}
	
	private static SummedMagFreqDist getBlankMFD() {
		return new SummedMagFreqDist(6.05,80,0.05);
	}


}
