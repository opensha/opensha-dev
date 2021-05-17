package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;

public class ERFRupLengthHistGen {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File outputDir = new File("/tmp");
		String prefix = "erf_rupture_lengths";
		
		double duration = 30;
		
		File probsDir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_erf_probs");
		Map<MagDependentAperiodicityOptions, ERF_ProbsZipFileReader> probsMap = Maps.newHashMap();
		probsMap.put(null, new ERF_ProbsZipFileReader(new File(probsDir, "probs_30yr_POISSON.zip")));
		probsMap.put(MagDependentAperiodicityOptions.LOW_VALUES, new ERF_ProbsZipFileReader(
				new File(probsDir, "probs_30yr_"+MagDependentAperiodicityOptions.LOW_VALUES.name()+".zip")));
		probsMap.put(MagDependentAperiodicityOptions.MID_VALUES, new ERF_ProbsZipFileReader(
				new File(probsDir, "probs_30yr_"+MagDependentAperiodicityOptions.MID_VALUES.name()+".zip")));
		probsMap.put(MagDependentAperiodicityOptions.HIGH_VALUES, new ERF_ProbsZipFileReader(
				new File(probsDir, "probs_30yr_"+MagDependentAperiodicityOptions.HIGH_VALUES.name()+".zip")));
		
		double min = 25d;
		int num = 25;
		double delta = 50d;
		
		Map<FaultModels, FaultSystemRupSet> rupSetsMap = new HashMap<>();
		File rupSetDir = new File("/home/kevin/workspace/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions");
		rupSetsMap.put(FaultModels.FM3_1, FaultSystemIO.loadRupSet(new File(rupSetDir,
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip")));
		rupSetsMap.put(FaultModels.FM3_2, FaultSystemIO.loadRupSet(new File(rupSetDir,
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip")));
		
		Collection<LogicTreeBranch> branches = CompoundFaultSystemSolution.fromZipFile(
				new File(rupSetDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip")).getBranches();
		
		Map<MagDependentAperiodicityOptions, Map<LogicTreeBranch, HistogramFunction>> u3Funcs = new HashMap<>();
		
		XY_DataSetList u3List = new XY_DataSetList();
		List<Double> u3Weights = new ArrayList<>();
		BranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		
		for (MagDependentAperiodicityOptions cov : probsMap.keySet()) {
			ERF_ProbsZipFileReader reader = probsMap.get(cov);
			Map<LogicTreeBranch, HistogramFunction> valsMap = u3Funcs.get(cov);
			if (valsMap == null) {
				valsMap = new HashMap<>();
				u3Funcs.put(cov, valsMap);
			}
			for (LogicTreeBranch branch : branches) {
				FaultSystemRupSet rupSet = rupSetsMap.get(branch.getValue(FaultModels.class));
				double[] probs = reader.getProbabilities(branch);
				Preconditions.checkState(probs.length == rupSet.getNumRuptures());
				HistogramFunction hist = new HistogramFunction(min, num, delta);
				for (int r=0; r<probs.length; r++) {
					double rate = -Math.log(1 - probs[r])/duration;
					double length = rupSet.getLengthForRup(r)/1000d;
					int xIndex = hist.getClosestXIndex(length);
					hist.add(xIndex, rate);
				}
				hist.normalizeBySumOfY_Vals();
				valsMap.put(branch, hist);
				u3List.add(hist);
				u3Weights.add(weightProv.getWeight(branch)*FaultSystemSolutionERF.getWeightForCOV(cov));
			}
		}
		
		FractileCurveCalculator u3FractileCalc = new FractileCurveCalculator(u3List, u3Weights);
		
		// now calc UCERF2
		MeanUCERF2 u2erf = new MeanUCERF2();
		u2erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, MeanUCERF2.PROB_MODEL_WGCEP_PREF_BLEND);
		u2erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
		u2erf.getTimeSpan().setStartTime(2014);
		u2erf.getTimeSpan().setDuration(duration);
		u2erf.updateForecast();
		
		HistogramFunction u2hist = new HistogramFunction(min, num, delta);
		for (ProbEqkSource source : u2erf) {
			for (ProbEqkRupture rup : source) {
				double rate = rup.getMeanAnnualRate(duration);
				double length = rup.getRuptureSurface().getAveLength();
				int xIndex = u2hist.getClosestXIndex(length);
				u2hist.add(xIndex, rate);
			}
		}
		u2hist.normalizeBySumOfY_Vals();
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		u2hist.setName("UCERF2 Mean");
		funcs.add(u2hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		
		AbstractXY_DataSet u3Mean = u3FractileCalc.getMeanCurve();
		u3Mean.setName("UCERF3 Mean");
		funcs.add(u3Mean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		AbstractXY_DataSet u3Min = u3FractileCalc.getMinimumCurve();
		u3Min.setName("UCERF3 Min");
		funcs.add(u3Min);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));

		AbstractXY_DataSet u3Max = u3FractileCalc.getMaximumCurve();
		u3Max.setName("UCERF3 Max");
		funcs.add(u3Max);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		String title = "Rupture Length Comparison";
		String xAxisLabel = "Rupture Length (km)";
		String yAxisLabel = "Fraction of Earthquakes";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		Range xRange = new Range(0, u2hist.getMaxX()+0.5*delta);
		Range yRange = new Range(1e-6, 1);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(20);
		plotPrefs.setAxisLabelFontSize(22);
		plotPrefs.setPlotLabelFontSize(24);
		plotPrefs.setLegendFontSize(22);
		plotPrefs.setBackgroundColor(Color.WHITE);
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
	}
	
	

}
