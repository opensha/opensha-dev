package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;

import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.SlipEnabledRupSet;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class PaleoVisibleDataWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		// for Devin McPhillips, request via e-mail 5/27/2021
		
		File outputDir = new File("/tmp/paleo_calcs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Map<FaultModels, U3FaultSystemSolution> fms = new HashMap<>();
		fms.put(FaultModels.FM3_1, U3FaultSystemIO.loadSol(
				new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip")));
		fms.put(FaultModels.FM3_2, U3FaultSystemIO.loadSol(
				new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_2_ucerf3.zip")));
		
		ScalingRelationships scale = ScalingRelationships.MEAN_UCERF3;
		
		for (FaultModels fm : fms.keySet()) {
			U3FaultSystemSolution sol = fms.get(fm);
			FaultSystemRupSet rupSet = sol.getRupSet();
			
			ArrayList<PaleoRateConstraint> constrs = UCERF3_PaleoRateConstraintFetcher.getConstraints(
					sol.getRupSet().getFaultSectionDataList());
			
			Map<Integer, Double> traceLengthCache = Maps.newConcurrentMap();
			
			for (PaleoRateConstraint constr : constrs) {
				CSVFile<String> csv = new CSVFile<>(true);
				
				csv.addLine("Rupture Index", "Magnitude", "Rate", "Ave Slip (m)", "x/L");
				
				for (int rupIndex : rupSet.getRupturesForSection(constr.getSectionIndex())) {
					double totOrigArea = 0d;
					for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex))
						totOrigArea += sect.getTraceLength()*1e3*sect.getOrigDownDipWidth()*1e3;
					double origDDW = totOrigArea/rupSet.getLengthForRup(rupIndex);
					double aveSlip = scale.getAveSlip(rupSet.getAreaForRup(rupIndex),
							rupSet.getLengthForRup(rupIndex), origDDW);
					
					double distanceAlong = UCERF3InversionInputGenerator.getDistanceAlongRupture(
							rupSet.getFaultSectionDataForRupture(rupIndex), constr.getSectionIndex(), traceLengthCache);
					
					csv.addLine(rupIndex+"", rupSet.getMagForRup(rupIndex)+"", sol.getRateForRup(rupIndex)+"",
							aveSlip+"", distanceAlong+"");
				}
				
				String siteName = constr.getPaleoSiteName().replaceAll("\\W+", "_");
				while (siteName.startsWith("_"))
					siteName = siteName.substring(1);
				while (siteName.endsWith("_"))
					siteName = siteName.substring(0, siteName.length()-1);
				File outputFile = new File(outputDir, fm.encodeChoiceString()+"_"+siteName+".csv");
				csv.writeToFile(outputFile);
			}
		}
	}

}
