package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateSegmentationConstraint.RateCombiner;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;

import com.google.common.base.Preconditions;

class SegResultPlots {
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures/crustal_sol");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/data/kevin/nshm23/batch_inversions/2024_12_05-prvi25_crustal_branches-dmSample5x/"
						+ "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged.zip"));
		
		ClusterRuptures cRups = sol.getRupSet().requireModule(ClusterRuptures.class);
		PlausibilityConfiguration config = sol.getRupSet().requireModule(PlausibilityConfiguration.class);
		ClusterConnectionStrategy connStrat = config.getConnectionStrategy();
		
		SegmentationCalculator.WRITE_PDFS = true;
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		// overwrite the segmentation choice so that the colors used in the dist-depend plot are uniform
		sol.requireModule(LogicTreeBranch.class).setValue(NSHM23_SegmentationModels.NONE);
		SegmentationCalculator segCalc = new SegmentationCalculator(sol, cRups.getAll(),
				connStrat, config.getDistAzCalc(), new double[] { 0d });
		
		segCalc.setLegendFontSize(26);
		
		segCalc.plotDistDependComparison(outputDir, "dist_dependence", true, null, " ");
		
		// now annotated map
		GeographicMapMaker mapMaker = new RupSetMapMaker(rupSet,
				SlipRateFigures.CRUSTAL_FAULT_MAP_REG);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(210, 210, 210)));
//		plotter.setJumpLineThickness(4f);
		
		CPT cpt = SegmentationCalculator.getConnectionFractCPT();
		
		String label = "Passthrough Rate";
		mapMaker.plotJumpScalars(segCalc.calcJumpPassthroughs(0, RateCombiner.MIN), cpt, label);
		
		mapMaker.plot(outputDir, "passthrough_map", " ", 800);
	}
}
