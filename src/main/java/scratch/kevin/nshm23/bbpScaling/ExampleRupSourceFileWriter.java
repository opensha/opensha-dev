package scratch.kevin.nshm23.bbpScaling;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;

public class ExampleRupSourceFileWriter {

	public static void main(String[] args) throws IOException {
		List<int[]> parentCombs = new ArrayList<>();
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		
		File mainOutputDir = new File("/home/kevin/CyberShake/nshm23_scaling");
		
		List<? extends FaultSection> allSubSects = rupSet.getFaultSectionDataList();
		
		Map<Integer, List<FaultSection>> groupedByParent = allSubSects.stream()
				.collect(java.util.stream.Collectors.groupingBy(FaultSection::getParentSectionId));
		
		List<String> sectAbbrevs = new ArrayList<>();
		List<Integer> sectIDs = new ArrayList<>();
		
		FaultSection singleSect = null;
		
//		File outputDir = new File(mainOutputDir, "southern_saf");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		
//		sectAbbrevs.add("BB");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Big Bend"));
//		singleSect = groupedByParent.get(sectIDs.get(0)).get(0);
//		
//		sectAbbrevs.add("MN");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Mojave", "north"));
//		
//		sectAbbrevs.add("MS");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Mojave", "south"));
//		
//		sectAbbrevs.add("SBN");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Bernardino", "north"));
//		
//		sectAbbrevs.add("NB");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "North", "Branch"));
//		
////		sectAbbrevs.add("SBS");
////		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Bernardino", "south"));
//		
////		sectAbbrevs.add("SGP");
////		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Gorgonio", "Garnet"));
//		
//		sectAbbrevs.add("CC");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Coachella"));
		
//		File outputDir = new File(mainOutputDir, "northridge");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		
//		sectAbbrevs.add("NR");
//		sectIDs.add(189);
//		singleSect = groupedByParent.get(sectIDs.get(0)).get(0);
		
		File outputDir = new File(mainOutputDir, "compton");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		sectAbbrevs.add("CMPT");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Compton"));
		singleSect = groupedByParent.get(sectIDs.get(0)).get(0);
		
		for (int i=0; i<sectAbbrevs.size(); i++) {
			String sectName = sectAbbrevs.get(i);
			int sectID = sectIDs.get(i);
			String fullName = groupedByParent.get(sectID).get(0).getParentSectionName();
			System.out.println(sectName+" is "+sectID+": "+fullName);
		}
		
		List<Integer> rupsToSearch = rupSet.getRupturesForParentSection(sectIDs.get(0));
		
		List<String> rupNames = new ArrayList<>();
		List<Integer> rupIndexes = new ArrayList<>();
		
		Set<Integer> curParentIDs = new HashSet<>();
		String curName = null;
		for (int i=0; i<sectAbbrevs.size(); i++) {
			String sectName = sectAbbrevs.get(i);
			int sectID = sectIDs.get(i);
			Preconditions.checkState(sectID >= 0);
			
			if (curName == null)
				curName = "";
			else
				curName += "_";
			
			curParentIDs.add(sectID);
			curName += sectName;
			
			double maxLen = 0d;
			int idForMaxLen = -1;
			
			for (int rupIndex : rupsToSearch) {
				List<FaultSection> sectsForRup = rupSet.getFaultSectionDataForRupture(rupIndex);
				boolean allMatch = true;
				Set<Integer> rupParentIDs = new HashSet<>(curParentIDs.size());
				for (FaultSection sect : sectsForRup) {
					if (!curParentIDs.contains(sect.getParentSectionId())) {
						allMatch = false;
						break;
					}
					rupParentIDs.add(sect.getParentSectionId());
				}
				if (allMatch && rupParentIDs.size() == curParentIDs.size()) {
					double length = rupSet.getLengthForRup(rupIndex) * 1e-3;
					if (length > maxLen) {
						maxLen = length;
						idForMaxLen = rupIndex;
					}
				}
			}
			
			Preconditions.checkState(idForMaxLen >= 0,
					"Didn't find any ruptures that use each of (and nothing else): %s\n\t%s", curName, curParentIDs);
			
			double avgMag = rupSet.getMagForRup(idForMaxLen);
			System.out.println("Longest rup found for "+curName+": "+(float)maxLen+" km, M"+(float)avgMag);
			
			if (singleSect != null && singleSect.getParentSectionId() == sectID) {
				// also add single fault ruptures
				for (int rupIndex : rupSet.getRupturesForSection(singleSect.getSectionId())) {
					if (rupIndex == idForMaxLen)
						continue;
					List<FaultSection> sectsForRup = rupSet.getFaultSectionDataForRupture(rupIndex);
					boolean allMatch = true;
					for (FaultSection sect : sectsForRup) {
						if (sect.getParentSectionId() != sectID) {
							allMatch = false;
							break;
						}
					}
					if (allMatch) {
						String rupName = curName+"_partial_"+sectsForRup.size()+"sects";
						double length = rupSet.getLengthForRup(rupIndex) * 1e-3;
						System.out.println("Also adding single-fault "+rupName+": "
								+(float)length+" km, M"+(float)rupSet.getMagForRup(rupIndex));
						rupNames.add(rupName);
						rupIndexes.add(rupIndex);
					}
				}
			}
			
			rupNames.add(curName);
			rupIndexes.add(idForMaxLen);
		}
		
		int fixedSeed = 123456789;
		double dWid = 0.1;
		double dLen = 0.1;
		double cornerFreq = 0.15;
		
		NSHM23_ScalingRelationships[] scales = {
				NSHM23_ScalingRelationships.AVERAGE,
				NSHM23_ScalingRelationships.LOGA_C4p1,
				NSHM23_ScalingRelationships.LOGA_C4p2,
				NSHM23_ScalingRelationships.LOGA_C4p3,
				NSHM23_ScalingRelationships.WIDTH_LIMITED,
				NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD,
				NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN
		};
		
		for (int r=0; r<rupIndexes.size(); r++) {
			int rupIndex = rupIndexes.get(r);
			String name = rupNames.get(r);
			
			BBP_PlanarSurface bbpSurf = planarEquivalentSurface(rupSet, rupIndex);
			
			double hypoAlong = 0.5*bbpSurf.getLength();
			double hypoDown = 0.5*bbpSurf.getWidth();
			
			// these are in SI units, not km
			double area = rupSet.getAreaForRup(rupIndex);
			double length = rupSet.getLengthForRup(rupIndex);
			double width = area/length;
			double totOrigArea = 0;
			for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex)) {
				totOrigArea += sect.getArea(false);	// sq-m
			}
			double origDDW = totOrigArea / length;
			double aveRake = rupSet.getAveRakeForRup(rupIndex);
//			mags[r] = scale.getMag(rupAreas[r], rupLengths[r], rupAreas[r]/rupLengths[r], origDDW, rakes[r]);
			
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Scaling Relationship", "NSHM23 Magnitude", "NSHM23 Average Slip (m)");
			
			for (NSHM23_ScalingRelationships scale : scales) {
				double mag = scale.getMag(area, length, width, origDDW, aveRake);
				double slip = scale.getAveSlip(area, length, width, origDDW, aveRake);
				csv.addLine(scale.getShortName(), (float)mag+"", (float)slip+"");
			}
			
			csv.writeToFile(new File(outputDir, name+"_scaling.csv"));
			
			BBP_SourceFile src = new BBP_SourceFile(bbpSurf, rupSet.getMagForRup(rupIndex), hypoAlong, hypoDown, dWid, dLen, cornerFreq, fixedSeed);
			src.writeToFile(new File(outputDir, name+".src"));
		}
	}
	
	public static BBP_PlanarSurface planarEquivalentSurface(FaultSystemRupSet rupSet, int rupIndex) {
		RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 1d);
		
		Location firstLoc = surf.getFirstLocOnUpperEdge();
		Location lastLoc = surf.getLastLocOnUpperEdge();
		
		LocationVector vector = LocationUtils.vector(firstLoc, lastLoc);
		double strike = vector.getAzimuth();
		double rake = rupSet.getAveRakeForRup(rupIndex);
		double dip = surf.getAveDip();
		
		FocalMechanism mech = new FocalMechanism(strike, dip, rake);
		
		double vectorLen = vector.getHorzDistance();
		
		LocationVector halfVector = new LocationVector(vector.getAzimuth(), 0.5*vectorLen, 0d);
		Location topCenter = LocationUtils.location(firstLoc, halfVector);
		topCenter = new Location(topCenter.getLatitude(), topCenter.getLongitude(), surf.getAveRupTopDepth());
		
		double length = surf.getAveLength();
		double width = surf.getAveWidth();
		
		return new BBP_PlanarSurface(topCenter, length, width, mech);
	}

}
