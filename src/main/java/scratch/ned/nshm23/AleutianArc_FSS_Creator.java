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

public class AleutianArc_FSS_Creator {
	
	private final static Boolean D = true;
	
	 public enum FaultModelEnum {
		 ALL(1.0), 
		 GEODETIC(0.167), // = 0.167 from Fig 9 of report
		 GEOLOGIC_NARROW(0.4165), // (1-0.167)/2
		 GEOLOGIC_WIDE(0.4165) ;  // 
//		 GEODETIC(0.0666), // = 0.2/3.0
//		 GEOLOGIC_NARROW(0.4667), // = (1-0.2/3.0)/2
//		 GEOLOGIC_WIDE(0.4667) ;  // = (1-0.2/3.0)/2
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

		String[] nameArray = {"geodetic","narrow","wide"}; // default case for "ALL"
		switch (fltMod) {
		case ALL:
			break;  
		case GEODETIC:
			nameArray = new String[1];
			nameArray[0] = "geodetic";
			break;  
		case GEOLOGIC_NARROW:
			nameArray = new String[1];
			nameArray[0] = "narrow";
			break;
		case GEOLOGIC_WIDE:
			nameArray = new String[1];
			nameArray[0] = "wide";
			break;
		default:
			throw new RuntimeException("Problem with fault model");
		}
		
		for(String fltModName:nameArray) {
			
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Adak ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Amchitka ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Andreanof ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Attu ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Barren Islands ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Fox Islands ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Kenai ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Kodiak ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Komandorski ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Prince William Sound ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Sanak ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Semidi ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Shumagin ("+fltModName+").geojson"));
			list.add(CEUS_FSS_creator.getFaultSection(nshmModelDirPath+"subduction/interface/Aleutian Arc/features/Yakataga ("+fltModName+").geojson"));
		}
		return list;
	}

	
	
	/**
	 * I found these files by listing for "rupture-set.json" at various directory depths
	 * @param srcIDsList
	 * @param srcFltSectsList
	 */
	private static void getSrcIDsAndFaultSectionsLists(HashMap<Integer,int[]> srcFltSectsMap, String nshmModelDirPath, FaultModelEnum fltMod) {
		
		if(fltMod == FaultModelEnum.GEOLOGIC_WIDE || fltMod == FaultModelEnum.ALL) {
//			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/unsegmented/wide/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/barren/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/kenai/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/wide/yakataga/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/adak-amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/amchitka-attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/andreanof-adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/attu-komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/barren/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/fox-andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/kenai-barren-kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/kodiak-semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/pws-kenai/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/wide/yakataga/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/adak-amchitka-attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/amchitka-attu-komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/andreanof-adak-amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/fox-andreanof-adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/kenai-barren-kodiak-semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/pws-kenai-barren-kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/wide/yakataga/rupture-set.json", srcFltSectsMap);
		}

		if(fltMod == FaultModelEnum.GEOLOGIC_NARROW || fltMod == FaultModelEnum.ALL) {
//			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/unsegmented/narrow/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/barren/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/kenai/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geologic/narrow/yakataga/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/adak-amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/amchitka-attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/andreanof-adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/attu-komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/barren/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/fox-andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/kenai-barren-kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/kodiak-semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/pws-kenai/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/2-sections/geologic/narrow/yakataga/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/adak-amchitka-attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/amchitka-attu-komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/andreanof-adak-amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/fox-andreanof-adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/kenai-barren-kodiak-semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/pws-kenai-barren-kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/3-sections/geologic/narrow/yakataga/rupture-set.json", srcFltSectsMap);
		}

		if(fltMod == FaultModelEnum.GEODETIC || fltMod == FaultModelEnum.ALL) {
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/adak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/amchitka/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/andreanof/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/attu/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/barren/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/fox/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/kenai/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/kodiak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/komandorski/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/pws/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/sanak/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/semidi/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/shumagin/rupture-set.json", srcFltSectsMap);
			CEUS_FSS_creator.parseRuptureSetFile(nshmModelDirPath+"subduction/interface/Aleutian Arc/1-section/geodetic/yakataga/rupture-set.json", srcFltSectsMap);
		}
	}


	public static FaultSystemSolution getFaultSystemSolution(String nshmModelDirPath, FaultModelEnum fltModel) {
		
		// rate weight for specified fault model branch
		double rateWt = 1.0/fltModel.getWeight();
		if(D)System.out.println("fltModel="+fltModel+"\nrateWt="+(float)rateWt);
		
	    // get fault section list for given fault model
	//	HashMap<Integer,GeoJSONFaultSection> faultSectionMap;
	    ArrayList<GeoJSONFaultSection> faultSectionList = getFaultSectionList(nshmModelDirPath, fltModel);
		// make parSectID_List & write attributes
		ArrayList<Integer> sectID_List = new ArrayList<Integer>(); // this is NSHM ID for each parent section
		if(D) System.out.println("index\tsectID\trake");
		for(int s=0;s<faultSectionList.size();s++) {
			GeoJSONFaultSection sect = faultSectionList.get(s);
			if(!sectID_List.contains(sect.getSectionId())) {
				sectID_List.add(sect.getSectionId());
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
		if (D)System.out.println("sectID_List.size()="+sectID_List.size());
		
		// Read from Peter's rupture-set.json files
		HashMap<Integer,int[]> srcFltSectsMap = new HashMap<Integer,int[]>(); // the fault section used by each source (same order as above)
		getSrcIDsAndFaultSectionsLists(srcFltSectsMap, nshmModelDirPath, fltModel);
		Set<Integer> srcIDsList = srcFltSectsMap.keySet();  // a list of all the source IDs (no duplicates)

		// some tests
		if(D) System.out.println("number of unique source IDs: "+srcIDsList.size());
		// make sure all sections are used and none are missing (with respect to Peter's files)
		ArrayList<Integer> testSectID_List = new ArrayList<Integer>(); // fill up with the par IDs from all the sources
		for(int srcID:srcIDsList) {
//			for(int i=0;i<srcIDsList.size();i++) {
			int[] sects = srcFltSectsMap.get(srcID);
			
// FOLLOWING RULE NOT APPLICABLE HERE
//			if(sects.length==1) {
//				if(sects[0] != srcID) // source ID = faultSection ID if only one fault section used
//					throw new RuntimeException("problem");
//			}

			if(D)System.out.print("\n"+srcID+"\t");
			for(int sect:sects) {
				if (D)System.out.print(sect+", ");	
				if(!testSectID_List.contains(sect))
					testSectID_List.add(sect);
			}
		}
		
		if (D) System.out.print("\n");
    	if(sectID_List.size() != testSectID_List.size()) {
    		System.out.println("parSectID_List:\n"+sectID_List+"\n"+sectID_List.size());
    		System.out.println("testSectID_List:\n"+testSectID_List+"\n"+testSectID_List.size());
    		for(int id:sectID_List) if(!testSectID_List.contains(id)) System.out.println("missing one: "+id);
    		throw new RuntimeException("parSectID_List.size() != testParSectID_List.size()");
    	}
    	for(int id : testSectID_List)
    		if(!sectID_List.contains(id))
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
    		// skip if not on the selected branch
    		if(!srcIDsList.contains(srcID))
    			continue;
    		System.out.println("\t"+s+"\t"+srcID+"\t"+src.getNumRuptures()+"\t"+src.getName()+"\t"+src.getRupture(0).getClass().toString());
    		if(!testSrcIDsList.contains(srcID))
    			testSrcIDsList.add(srcID);
    		else
    			System.out.println("Duplicate srcID = "+srcID);
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
		    	int newID = newFltIndexMap.get(sectID);
    			System.out.println("RI for sect "+newID+"\t"+(1.0/mfd2.getTotalIncrRate())+"\t"+faultSectionList.get(newID).getName());
		    }
		    System.out.println("ERF versus fss tests passed!!");
	    }   
		return fss;
	}
	


	public static void main(String[] args) {
//		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-main_Jan10_2024/"; // old & not longer works
		
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-3.0.1/";
//		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.ALL);
		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.GEODETIC);
//		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.GEOLOGIC_NARROW);
//		getFaultSystemSolution(nshmModelDirPath, FaultModelEnum.GEOLOGIC_WIDE);
		
		
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
