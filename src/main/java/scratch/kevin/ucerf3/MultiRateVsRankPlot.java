package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;

import scratch.UCERF3.simulatedAnnealing.completion.AnnealingProgress;

public class MultiRateVsRankPlot {

	public static void main(String[] args) throws IOException {
		List<File> solFiles = new ArrayList<>();
		List<String> solLabels = new ArrayList<>();
		List<FaultSystemSolution> sols = new ArrayList<>();
		
		File baseDir = new File("/data/kevin/markdown/inversions");
		File outputDir = new File("/tmp");
		
//		String suffix = "_2hr";

//		solFiles.add(new File(baseDir, "2021_07_30-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-8-avg_anneal-5m-noWL-run0/solution.zip"));
//		solLabels.add("Perturb Scalar [1e-2,1e-8]");
//		
//		solFiles.add(new File(baseDir, "2021_07_29-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-10-avg_anneal-5m-noWL-run0/solution.zip"));
//		solLabels.add("Perturb Scalar [1e-2,1e-10]");
//		
//		solFiles.add(new File(baseDir, "2021_07_28-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-12-avg_anneal-5m-noWL-run0/solution.zip"));
//		solLabels.add("Perturb Scalar [1e-2,1e-12]");
//		
//		solFiles.add(new File(baseDir, "2021_07_29-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-14-avg_anneal-5m-noWL-run0/solution.zip"));
//		solLabels.add("Perturb Scalar [1e-2,1e-14]");
//		
//		solFiles.add(new File(baseDir, "2021_07_28-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-16-avg_anneal-5m-noWL-run0/solution.zip"));
//		solLabels.add("Perturb Scalar [1e-2,1e-16]");
		
		String suffix = "_5hr";
		
		solFiles.add(new File(baseDir, "2021_08_03-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-10-avg_anneal_5m-noWL-5hr-run0/solution.zip"));
		solLabels.add("Perturb Scalar [1e-2,1e-10], 5hr");
		
		solFiles.add(new File(baseDir, "2021_08_03-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-12-avg_anneal_5m-noWL-5hr-run0/solution.zip"));
		solLabels.add("Perturb Scalar [1e-2,1e-12], 5hr");
		
		solFiles.add(new File(baseDir, "2021_08_04-coulomb-u3_ref-perturb_exp_scale_1e-2_to_1e-14-avg_anneal_5m-noWL-5hr-run0/solution.zip"));
		solLabels.add("Perturb Scalar [1e-2,1e-14], 5hr");
		
		for (File solFile : solFiles)
			sols.add(FaultSystemSolution.load(solFile));
		
		List<DiscretizedFunc> rateFuncs = new ArrayList<>();
		List<DiscretizedFunc> energyFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, sols.size()-1);
		
		for (int s=0; s<sols.size(); s++) {
			FaultSystemSolution sol = sols.get(s);
			
			List<Double> values = new ArrayList<>(sol.getRupSet().getNumRuptures());
			for (double val : sol.getRateForAllRups())
				values.add(val);
			Collections.sort(values);
			Collections.reverse(values);
			
			DiscretizedFunc rateVsRank = new ArbitrarilyDiscretizedFunc(solLabels.get(s));
			for (int r=0; r<values.size(); r++)
				rateVsRank.set((double)r, values.get(r));
			
			DiscretizedFunc energyVsTime = new ArbitrarilyDiscretizedFunc(solLabels.get(s));
			AnnealingProgress progress = sol.requireModule(AnnealingProgress.class);
			for (int i=0; i<progress.size(); i++) {
				long millis = progress.getTime(i);
				double secs = millis/1000d;
				double mins = secs/60d;
				energyVsTime.set(mins, progress.getEnergies(i)[0]);
			}
			
			Color color = cpt.getColor(s);
			
			rateFuncs.add(rateVsRank);
			energyFuncs.add(energyVsTime);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3, color));
		}
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		PlotSpec rateSpec = new PlotSpec(rateFuncs, chars, "Rate vs Rank", "Rank", "Annual Rate");
		rateSpec.setLegendInset(true);
		
		gp.drawGraphPanel(rateSpec, false, true, null, new Range(1e-16, 1e-2));
		
		PlotUtils.writePlots(outputDir, "rate_vs_rank"+suffix, gp, 1000, 850, true, false, false);
		
		PlotSpec energySpec = new PlotSpec(energyFuncs, chars, "Energy vs Time", "Annealing Time (miniutes)", "Total Energy");
		energySpec.setLegendInset(true);
		
		gp.drawGraphPanel(energySpec, false, false, null, new Range(0, 300));
		
		PlotUtils.writePlots(outputDir, "energy_vs_time"+suffix, gp, 1000, 850, true, false, false);
	}

}
