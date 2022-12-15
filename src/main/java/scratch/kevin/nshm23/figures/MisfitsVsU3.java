package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;

import com.google.common.base.Preconditions;

class MisfitsVsU3 {

	public static void main(String[] args) throws IOException {
		File u3File = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged_full_modules.zip");
		File newFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_06-nshm23_u3_hybrid_branches-no_paleo_slip-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_FM3_1_CoulombRupSet_NoClassic_branch_averaged.zip");
		String newName = "NSHM23 Methodology"; String prefixAdd = "_no_classic";
//				+ "results_FM3_1_CoulombRupSet_branch_averaged.zip");
//		String newName = "NSHM23 Methodology"; String prefixAdd = "";
//				+ "node_branch_averaged/PaleoUncert_EvenFitPaleo.zip");
//		String newName = "NSHM23 Methodology"; String prefixAdd = "_even_paleo";
		
		File outputDir = new File("/tmp/u3_misfit_hists");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution u3Sol = FaultSystemSolution.load(u3File);
		FaultSystemSolution newSol = FaultSystemSolution.load(newFile);
		
		// slip rates
		for (boolean pDiff : new boolean[] {false,true}) {
			double[] u3Ratios = calcSlipRatios(u3Sol, pDiff);
			double[] newRatios = calcSlipRatios(newSol, pDiff);
			
			HistogramFunction u3Hist;
			String xAxisLabel;
			String prefix;
			if (pDiff) {
				u3Hist = HistogramFunction.getEncompassingHistogram(-49, 49d, 2.5d);
				xAxisLabel = "Solution vs Target Slip Rate, % Difference";
				prefix = "slip_rate_pDiffs";
			} else {
				u3Hist = HistogramFunction.getEncompassingHistogram(0.01d, 1.99d, 0.05);
				xAxisLabel = "Solution / Target, Slip Rate Ratio";
				prefix = "slip_rate_ratios";
			}
			HistogramFunction newHist = new HistogramFunction(u3Hist.getMinX(), u3Hist.getMaxX(), u3Hist.size());
			
			u3Hist.setName("UCERF3");
			newHist.setName(newName);
			
			populateHist(u3Hist, u3Ratios);
			populateHist(newHist, newRatios);
			
			plotComparisonHists(outputDir, prefix+prefixAdd, "Slip Rate Misfits", xAxisLabel, "Subsection Count",
					u3Hist, newHist, pDiff ? 10d : 0.2d);
		}
		
		// paleo data
		double[] u3PaleoZ = calcPaleoZScores(u3Sol);
		double[] newPaleoZ = calcPaleoZScores(newSol);
		
		HistogramFunction u3ZHist = HistogramFunction.getEncompassingHistogram(-2.99d, 2.99d, 0.25);
		HistogramFunction newZHist = new HistogramFunction(u3ZHist.getMinX(), u3ZHist.getMaxX(), u3ZHist.size());
		
		u3ZHist.setName("UCERF3");
		newZHist.setName(newName);
		
		populateHist(u3ZHist, u3PaleoZ);
		populateHist(newZHist, newPaleoZ);
		
		System.out.println("UCERF3 z-score stats:");
		double avg = 0d;
		double absAvg = 0d;
		for (double z : u3PaleoZ) {
			avg += z;
			absAvg += Math.abs(z);
		}
		avg /= u3PaleoZ.length;
		absAvg /= u3PaleoZ.length;
		System.out.println("\tAvg: "+(float)avg);
		System.out.println("\tAvg Absolute: "+(float)absAvg);
		
		System.out.println(newName+" z-score stats:");
		avg = 0d;
		absAvg = 0d;
		for (double z : newPaleoZ) {
			avg += z;
			absAvg += Math.abs(z);
		}
		avg /= newPaleoZ.length;
		absAvg /= newPaleoZ.length;
		System.out.println("\tAvg: "+(float)avg);
		System.out.println("\tAvg Absolute: "+(float)absAvg);
		
		plotComparisonHists(outputDir, "paleo_z_scores"+prefixAdd, "Paleoseismic Data Misfits", "Paleoseismic Rate z-score",
				"Paleoseismic Constraint Count", u3ZHist, newZHist, 0.5);
	}
	
	private static double[] calcSlipRatios(FaultSystemSolution sol, boolean pDiff) {
		double[] ret = new double[sol.getRupSet().getNumSections()];
		
		SectSlipRates targets = sol.getRupSet().requireModule(SectSlipRates.class);
		SolutionSlipRates solSlips = sol.requireModule(SolutionSlipRates.class);
		
		for (int s=0; s<ret.length; s++) {
			double solVal = solSlips.get(s);
			double targetVal = targets.getSlipRate(s);
			
			if (pDiff)
				ret[s] = 100d*(solVal - targetVal)/targetVal;
			else
				ret[s] = solVal/targetVal;
		}
		
		return ret;
	}
	
	private static double[] calcPaleoZScores(FaultSystemSolution sol) {
		PaleoseismicConstraintData data = sol.getRupSet().requireModule(PaleoseismicConstraintData.class);
		
		List<? extends SectMappedUncertainDataConstraint> rateConstraints = data.getPaleoRateConstraints();
		
		double[] ret = new double[rateConstraints.size()];
		
		for (int i=0; i<ret.length; i++) {
			SectMappedUncertainDataConstraint rateConstraint = rateConstraints.get(i);
			
			double solRate = sol.calcTotPaleoVisibleRateForSect(rateConstraint.sectionIndex, data.getPaleoProbModel());
			
			ret[i] = (solRate - rateConstraint.bestEstimate)/rateConstraint.getPreferredStdDev();
		}
		
		return ret;
	}
	
	private static void populateHist(HistogramFunction hist, double[] values) {
		for (double v : values)
			hist.add(hist.getClosestXIndex(v), 1d);
	}
	
	private static void plotComparisonHists(File outputDir, String prefix, String title,
			String xAxisLabel, String yAxisLabel, HistogramFunction u3Hist, HistogramFunction newHist,
			Double customXTick) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(u3Hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		
		funcs.add(newHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		
		HistogramFunction overlap = new HistogramFunction(u3Hist.getMinX(), u3Hist.getMaxX(), u3Hist.size());
		overlap.setName("Overlap");
		for (int i=0; i<overlap.size(); i++)
			overlap.set(i, Math.min(u3Hist.getY(i), newHist.getY(i)));
		funcs.add(overlap);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, new Color(200, 100, 200)));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendInset(true);
//		spec.setLegendVisible(true);
		
		Range xRange = new Range(u3Hist.getMinX()-0.5*u3Hist.getDelta(), u3Hist.getMaxX()+0.5*u3Hist.getDelta());
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setTickLabelFontSize(24);
		gp.setAxisLabelFontSize(26);
		gp.setPlotLabelFontSize(30);
		
		gp.drawGraphPanel(spec, false, false, xRange, null);
		
		if (customXTick != null)
			PlotUtils.setXTick(gp, customXTick);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 850, 700, true, true, false);
	}

}
