package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class SectMFD_Plots {
	
	public static void main(String[] args) throws IOException {
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/sect_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution u3Sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip")); // must be 3.1 for Pitas
		
		FaultSystemSolution methodsSol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		
		FaultSystemSolution modelSol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR"
				+ "/results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
//		String faultName = "Pitas Point (lower, Montalvo)";
//		String[] faultSearch = { "Pitas", "Montalvo" };
//		String prefix = "pitas_point_montalvo";
//		Range xRange = new Range(6.5d, 8.6);
//		Range yRange = new Range(1e-8, 1e-2);
		
//		String faultName = "San Cayetano";
//		String[] faultSearch = { "Cayetano" };
//		String prefix = "san_cayetano";
//		Range xRange = new Range(6.5d, 8.6);
//		Range yRange = new Range(1e-8, 1e-2);
		
//		String faultName = "Death Valley (No)";
//		String[] faultSearch = { "Death", "Valley", "No" };
//		String prefix = "death_valley_no";
//		Range xRange = new Range(6d, 8.2);
//		Range yRange = new Range(1e-8, 1e-2);
		
		String faultName = "Death Valley (Black Mtns Frontal)";
		String[] faultSearch = { "Death", "Valley", "Black" };
		String prefix = "death_valley_black_mtns";
		Range xRange = new Range(6d, 8.2);
		Range yRange = new Range(1e-8, 1e-2);
		
		List<FaultSystemSolution> sols = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<Color> colors = new ArrayList<>();
		
		sols.add(u3Sol);
		names.add("UCERF3");
		colors.add(Color.BLUE);
		
		sols.add(methodsSol);
		names.add("NSHM23 Methodology");
		colors.add(Color.GREEN.darker());
		
		sols.add(modelSol);
		names.add("NSHM23 Model");
		colors.add(Color.RED);
		
		List<DiscretizedFunc> incrFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
		List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> cmlChars = new ArrayList<>();
		
		for (int s=0; s<sols.size(); s++) {
			FaultSystemSolution sol = sols.get(s);
			String name = names.get(s);
			Color color = colors.get(s);
			
			FaultSystemRupSet rupSet = sol.getRupSet();
			
			int parentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), faultSearch);
			
			Preconditions.checkState(parentID > 0);
			IncrementalMagFreqDist incrMFD = sol.calcParticipationMFD_forParentSect(
					parentID, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			incrMFD.setName(name);
			EvenlyDiscretizedFunc cmlMFD = incrMFD.getCumRateDistWithOffset();
			
			incrFuncs.add(incrMFD);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
			
			cmlFuncs.add(cmlMFD);
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
		}
		
		// add nshm distr
		int numCopyLater = incrFuncs.size();
		
		BranchParentSectParticMFDs branchMFDs = modelSol.requireModule(BranchParentSectParticMFDs.class);
		int parentID = FaultSectionUtils.findParentSectionID(modelSol.getRupSet().getFaultSectionDataList(), faultSearch);
		
		IncrementalMagFreqDist[] incrPercentiles = branchMFDs.calcIncrementalSectFractiles(parentID, SolMFDPlot.standardFractiles);
		
		Color transColor = new Color(255, 0, 0, 30);
		PlotCurveCharacterstics minMaxChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, transColor);
		
		for (IncrementalMagFreqDist bounds : SolMFDPlot.processIncrFractiles(incrPercentiles)) {
			incrFuncs.add(bounds);
			incrChars.add(minMaxChar);
		}
		
		EvenlyDiscretizedFunc[] cmlPercentiles = branchMFDs.calcCumulativeSectFractiles(parentID, SolMFDPlot.standardFractiles);
		for (UncertainArbDiscFunc cmlBounds : SolMFDPlot.processCmlFractiles(cmlPercentiles, xRange.getLowerBound())) {
			cmlFuncs.add(cmlBounds);
			cmlChars.add(minMaxChar);
		}
		
		// now copy old ones on top
		for (int i=0; i<numCopyLater; i++) {
			DiscretizedFunc incr = incrFuncs.get(i).deepClone();
			incr.setName(null);
			DiscretizedFunc cml = cmlFuncs.get(i).deepClone();
			cml.setName(null);
			
			incrFuncs.add(incr);
			incrChars.add(incrChars.get(i));
			cmlFuncs.add(cml);
			cmlChars.add(cmlChars.get(i));
		}
		
		PlotSpec incrSpec = new PlotSpec(incrFuncs, incrChars, faultName, "Magnitude", "Incremental Rate (/yr)");
		incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		PlotSpec cmlSpec = new PlotSpec(cmlFuncs, cmlChars, faultName, "Magnitude", "Cumulative Rate (/yr)");
		cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setAxisLabelFontSize(26);
		gp.setTickLabelFontSize(22);
		
		gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_incr", gp, 800, 750, true, true, false);
		
		gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, true, false);
	}

}
