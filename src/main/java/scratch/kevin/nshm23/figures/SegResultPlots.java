package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateSegmentationConstraint.RateCombiner;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;

import com.google.common.base.Preconditions;

public class SegResultPlots {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/seg_results");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File modelDir = new File(invsDir, "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		FaultSystemSolution fullSol = FaultSystemSolution.load(
				new File(modelDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		FaultSystemSolution noneSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_None.zip"));
		
		FaultSystemSolution lowSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_LowSeg.zip"));
		
		FaultSystemSolution middleSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_MidSeg.zip"));
		
		FaultSystemSolution highSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_HighSeg.zip"));
		
		FaultSystemSolution classicSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_Classic.zip"));
		
		ClusterRuptures cRups = fullSol.getRupSet().requireModule(ClusterRuptures.class);
		PlausibilityConfiguration config = fullSol.getRupSet().requireModule(PlausibilityConfiguration.class);
		ClusterConnectionStrategy connStrat = config.getConnectionStrategy();
		
		FaultSystemSolution[] sols = {
				fullSol,
				noneSol,
				lowSol,
				middleSol,
				highSol,
				classicSol
		};
		String[] prefixes = {
				"full_ba",
				"none",
				"low",
				"middle",
				"high",
				"classic"
		};
		
		SegmentationCalculator.WRITE_PDFS = true;
		
		for (int i=0; i<sols.length; i++) {
			// overwrite the segmentation choice so that the colors used in the dist-depend plot are uniform
			sols[i].requireModule(LogicTreeBranch.class).setValue(NSHM23_SegmentationModels.NONE);
			SegmentationCalculator segCalc = new SegmentationCalculator(sols[i], cRups.getAll(),
					connStrat, config.getDistAzCalc(), new double[] { 0d });
			String prefix = prefixes[i];
			
			segCalc.setLegendFontSize(24);
			
			segCalc.plotConnectionFracts(outputDir, "passthrough_map_"+prefix, " ");
			
			segCalc.plotDistDependComparison(outputDir, "dist_dependence_"+prefix, true, null, " ");
			
			if (i == 0) {
				// write an empty plot for Ned
				segCalc = new SegmentationCalculator(sols[i], cRups.getAll(),
						connStrat, config.getDistAzCalc(), new double[] { 10d });
				segCalc.plotDistDependComparison(outputDir, "dist_dependence_empty", true, null, " ");
			}
		}
	}

}
