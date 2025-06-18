package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupCartoonGenerator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class MuertosMagnitudeTests {
	
	public static void main(String[] args) throws IOException {
		List<LogicTreeLevel<? extends RupSetScalingRelationship>> scaleLevels =
				List.of(PRVI25_LogicTreeBranch.CRUSTAL_SCALE, PRVI25_LogicTreeBranch.SUB_SCALE);
		
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		FaultSystemSolution crustalSolWithMuertos = FaultSystemSolution.load(new File(
				"/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2025_05_27-prvi25_crustal_branches-dmSample10x-mue_as_crustal/"
				+ "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged.zip"));
		FaultSystemRupSet crustalRupSetWithMuertos = crustalSolWithMuertos.getRupSet();
		
		int muertosParentID = FaultSectionUtils.findParentSectionID(crustalRupSetWithMuertos.getFaultSectionDataList(), "Muertos");
		List<Integer> allMuertosRups = crustalRupSetWithMuertos.getRupturesForParentSection(muertosParentID);
		List<Integer> muertosSingleFaultRups = new ArrayList<>();
		double rateSingle = 0d;
		double rateTotal = 0d;
		ClusterRuptures cRups = crustalRupSetWithMuertos.requireModule(ClusterRuptures.class);
		int cartoonPrintMod = 20;
		double minMagAlwaysPlot = 7.5d;
		File outputDir = new File("/tmp/muertos_rups");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(), "Failed to create output dir: %s",
				outputDir.getAbsolutePath());
		for (int rupIndex : allMuertosRups) {
			rateTotal += crustalSolWithMuertos.getRateForRup(rupIndex);
			boolean single = true;
			for (FaultSection sect : crustalRupSetWithMuertos.getFaultSectionDataForRupture(rupIndex)) {
				if (sect.getParentSectionId() != muertosParentID) {
					single = false;
					break;
				}
			}
			if (single) {
				muertosSingleFaultRups.add(rupIndex);
				rateSingle += crustalSolWithMuertos.getRateForRup(rupIndex);
			}
		}
		System.out.println("Have "+allMuertosRups.size()+" muertos rups ("+(
				allMuertosRups.size()-muertosSingleFaultRups.size())+" multi-fault)");
		
		// sort by mag
		Map<Integer, Double> rupMags = new HashMap<>();
		for (int rupIndex : allMuertosRups)
			rupMags.put(rupIndex, crustalRupSetWithMuertos.getMagForRup(rupIndex));
		List<Integer> rupsMagSorted = ComparablePairing.getSortedData(rupMags);
		
		DecimalFormat magDF = new DecimalFormat("0.0");
		System.out.println("Smallest rupture: "+rupsMagSorted.get(0)+": "+magDF.format(rupMags.get(rupsMagSorted.get(0))));
		System.out.println(cRups.get(rupsMagSorted.get(0)).toString(true));
		
		for (int i=0; i<rupsMagSorted.size(); i++) {
			int rupIndex = rupsMagSorted.get(i);
			double mag = rupMags.get(rupIndex);
			if (i % cartoonPrintMod == 0 || mag <= minMagAlwaysPlot) {
				RupCartoonGenerator.plotRupture(outputDir, "rup_"+i, cRups.get(rupIndex),
						"Rupture "+rupIndex+", M"+magDF.format(mag), false, true);
			}
		}
		
		for (LogicTreeLevel<? extends RupSetScalingRelationship> level : scaleLevels) {
			System.out.println(level.getName());
			
			MinMaxAveTracker levelAllRupsTrack = new MinMaxAveTracker();
			MinMaxAveTracker levelSingleFaultTrack = new MinMaxAveTracker();
			
			for (RupSetScalingRelationship scale : level.getNodes()) {
				if (scale.getNodeWeight(null) == 0d)
					continue;
				FaultSystemRupSet rupSetForScale = FaultSystemRupSet.builder(crustalRupSetWithMuertos.getFaultSectionDataList(),
						crustalRupSetWithMuertos.getSectionIndicesForAllRups())
						.forScalingRelationship(scale).build();
				
				MinMaxAveTracker allRupsTrack = new MinMaxAveTracker();
				MinMaxAveTracker singleFaultTrack = new MinMaxAveTracker();
				
				for (int rupIndex : allMuertosRups)
					allRupsTrack.addValue(rupSetForScale.getMagForRup(rupIndex));
				for (int rupIndex : muertosSingleFaultRups)
					singleFaultTrack.addValue(rupSetForScale.getMagForRup(rupIndex));
				
				levelAllRupsTrack.addFrom(allRupsTrack);
				levelSingleFaultTrack.addFrom(singleFaultTrack);
				
				System.out.println("\t"+scale.getName());
				System.out.println("\t\tSingle-fault:\t["+magDF.format(singleFaultTrack.getMin())+", "+magDF.format(singleFaultTrack.getMax())+"]");
				System.out.println("\t\tAll ruptures:\t["+magDF.format(allRupsTrack.getMin())+", "+magDF.format(allRupsTrack.getMax())+"]");
			}
			
			System.out.println("\tOverall");
			System.out.println("\t\tSingle-fault:\t["+magDF.format(levelSingleFaultTrack.getMin())+", "+magDF.format(levelSingleFaultTrack.getMax())+"]");
			System.out.println("\t\tAll ruptures:\t["+magDF.format(levelAllRupsTrack.getMin())+", "+magDF.format(levelAllRupsTrack.getMax())+"]");
			
			System.out.println();
		}
		
		System.out.println("Single-fault total rate as crustal:\t"+(float)rateSingle+" ("+magDF.format(1d/rateSingle)+" yrs)");
		System.out.println("Overall total rate as crustal:\t"+(float)rateTotal+" ("+magDF.format(1d/rateTotal)+" yrs)");
	}

}
