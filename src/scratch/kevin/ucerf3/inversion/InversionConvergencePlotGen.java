package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class InversionConvergencePlotGen {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ZipException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File solDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/" +
				"scratch/InversionSolutions/");
		File compoundSol = new File(solDir, "2013_01_14-stampede_3p2_production_runs_combined_" +
				"MEAN_COMPOUND_SOL_with_indv_runs.zip");
		ZipFile zip = new ZipFile(compoundSol);
		
		File outputDir = new File("/tmp/convergence");
		if (!outputDir.exists())
			outputDir.mkdir();
		
		File avgFile = new File(solDir, "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_" +
				"NoFix_SpatSeisU3_mean_sol.zip");
		AverageFaultSystemSolution avgSol = FaultSystemIO.loadAvgInvSol(avgFile);
		
		File avgSol0rates = new File(solDir, "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_" +
				"M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_run00.bin");
		File avgSol0NoMins = new File(solDir, "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_" +
				"M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_run00_noMinRates.bin");
		double[] waterlevels = loadWaterlevels(avgSol0rates, avgSol0NoMins);
		writeWaterlevelResults(outputDir, "ref_branch",
				loadNonWaterlevelsForAvg(avgSol, waterlevels), avgSol.getRupSet().getNumRuptures());
		
		// these are for branch averages
		// now non waterlevel plots
		InversionFaultSystemRupSet rupSet = FaultSystemIO.loadInvRupSet(
				new File(compoundSol.getParentFile(),
						"2013_01_14-stampede_3p2_production_runs_combined_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		File fm3_1_bins = new File(solDir,
				"2013_01_14-stampede_3p2_production_runs_combined_FM3_1_bins.zip");
		writeWaterlevelResults(outputDir, "branch_runs_fm3_1_run0",
				loadNonWaterlevelsForBranches(new ZipFile(fm3_1_bins), "run0"), rupSet.getNumRuptures());
		writeResults(loadBranchAverages(zip, FaultModels.FM3_1, rupSet.getNumRuptures()),
				outputDir, "mean_fm3_1", rupSet);
		writeResults(loadBranchRuns(zip, FaultModels.FM3_1, "_rates_0.bin"),
				outputDir, "branch_runs_fm3_1_run0", rupSet);
		writeResults(loadBranchRuns(zip, FaultModels.FM3_1, null),
				outputDir, "branch_runs_fm3_1_mean", rupSet);
		
		rupSet = FaultSystemIO.loadInvRupSet(
				new File(compoundSol.getParentFile(),
						"2013_01_14-stampede_3p2_production_runs_combined_FM3_2_MEAN_BRANCH_AVG_SOL.zip"));
		writeResults(loadBranchAverages(zip, FaultModels.FM3_2, rupSet.getNumRuptures()),
				outputDir, "mean_fm3_2", rupSet);
		writeResults(loadBranchRuns(zip, FaultModels.FM3_2, "_rates_0.bin"),
				outputDir, "branch_runs_fm3_2_run0", rupSet);
		writeResults(loadBranchRuns(zip, FaultModels.FM3_2, null),
				outputDir, "branch_runs_fm3_2_mean", rupSet);
		
		writeResults(avgSol.getRatesForAllSols(), outputDir, "ref_branch", avgSol.getRupSet());
	}
	
	private static List<double[]> loadBranchAverages(ZipFile file, FaultModels fm, int numRups)
			throws IOException {
		Map<String, List<ZipEntry>> entriesMap = Maps.newHashMap();
		
		for (ZipEntry entry : Collections.list(file.entries())) {
			String name = entry.getName();
			if (!name.contains("_rates_"))
				continue;
			
			if (!name.startsWith(fm.getShortName()))
				continue;
			
			String rates_sub = name.substring(name.indexOf("_rates_"));
			
			List<ZipEntry> list = entriesMap.get(rates_sub);
			if (list == null) {
				list = Lists.newArrayList();
				entriesMap.put(rates_sub, list);
			}
			
			list.add(entry);
//			list.add(MatrixIO.doubleArrayFromInputStream(file.getInputStream(entry), entry.getSize()));
		}
		
		List<double[]> retList = Lists.newArrayList();
		for (String key : entriesMap.keySet()) {
			List<ZipEntry> entries = entriesMap.get(key);
			int numSols = entries.size();
			
			System.out.println("Averaging "+numSols+" runs for: "+key);
			double[] meanRates = new double[numRups];
			double rateMult = 1d/(double)numSols;
			for (ZipEntry entry : entries) {
				double[] rates = MatrixIO.doubleArrayFromInputStream(
						file.getInputStream(entry), entry.getSize());
				for (int r=0; r<meanRates.length; r++)
					meanRates[r] += rates[r]*rateMult;
			}
			retList.add(meanRates);
		}
		
		entriesMap = null;
		System.gc();
		
		return retList;
	}
	
	private static List<double[]> loadBranchRuns(ZipFile file, FaultModels fm, String nameGrep)
			throws IOException {
		List<double[]> retList = Lists.newArrayList();
		for (ZipEntry entry : Collections.list(file.entries())) {
			String name = entry.getName();
			
			if (!name.startsWith(fm.getShortName()))
				continue;
			
			if (nameGrep == null) {
				// use mean
				if (!name.contains("_rates.bin"))
					continue;
			} else {
				if (!name.contains(nameGrep))
					continue;
			}
			
			retList.add(MatrixIO.doubleArrayFromInputStream(file.getInputStream(entry), entry.getSize()));
		}
		
		System.gc();
		
		return retList;
	}
	
	private static void writeResults(List<double[]> ratesList, File dir, String prefix, InversionFaultSystemRupSet rupSet)
			throws IOException {
		int numRates = ratesList.get(0).length;
		int numLists = ratesList.size();
		List<RateEntry> entries = Lists.newArrayList();
		
		for (int r=0; r<numRates; r++) {
			double[] rates = new double[numLists];
			for (int i=0; i<numLists; i++)
				rates[i] = ratesList.get(i)[r];
			entries.add(new RateEntry(r, rates));
		}
		
		EvenlyDiscretizedFunc normStdDevFunc = new EvenlyDiscretizedFunc(0d, numRates, 1d);
		normStdDevFunc.setName("Std Dev / Mean");
		EvenlyDiscretizedFunc normRateFunc = new EvenlyDiscretizedFunc(0d, numRates, 1d);
		normRateFunc.setName("Cumulative Rate (normalized)");
		EvenlyDiscretizedFunc normMoRateFunc = new EvenlyDiscretizedFunc(0d, numRates, 1d);
		normMoRateFunc.setName("Cumulative Moment Rate (normalized)");
		DefaultXY_DataSet normStdDevVsRateFunc = new DefaultXY_DataSet();
		normStdDevVsRateFunc.setName("Std Dev / Mean vs Rate");
		DefaultXY_DataSet normStdDevVsMoRateFunc = new DefaultXY_DataSet();
		normStdDevVsMoRateFunc.setName("Std Dev / Mean vs Moment Rate");

		DefaultXY_DataSet nVsRateFunc = new DefaultXY_DataSet();
		nVsRateFunc.setName("N for SDOM within 10% vs Rate");
		DefaultXY_DataSet nVsMoRateFunc = new DefaultXY_DataSet();
		nVsMoRateFunc.setName("N for SDOM within 10% vs Rate");
		
		Collections.sort(entries);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		List<String> header = Lists.newArrayList(
				"Index", "Mean", "Std Dev", "Std Dev / Mean", "Min", "Max");
		
		csv.addLine(header);
		
		double rateSum = 0;
		double moRateSum = 0;
		for (int i=0; i<numRates; i++) {
			RateEntry rate = entries.get(i);
			List<String> line = Lists.newArrayList();
			
			line.add(rate.index+"");
			line.add(rate.mean+"");
			line.add(rate.stdDev+"");
			line.add(rate.normStdDev+"");
			line.add(StatUtils.min(rate.vals)+"");
			line.add(StatUtils.max(rate.vals)+"");
			
			csv.addLine(line);
			
			normStdDevFunc.set(i, rate.normStdDev);
			rateSum += rate.mean;
			normRateFunc.set(i, rateSum);
			double moRate = rate.mean*FaultMomentCalc.getMoment(
					rupSet.getAreaForRup(rate.index), rupSet.getAveSlipForRup(rate.index));
			moRateSum += moRate;
			normMoRateFunc.set(i, moRateSum);
			if (rate.mean > 0) {
				normStdDevVsRateFunc.set(rate.normStdDev, rate.mean);
				normStdDevVsMoRateFunc.set(rate.normStdDev, moRate);
//				if (nVsRateFunc.getNumCollisions() % 100 == 0 && nVsRateFunc.getNumCollisions()>0)
//					System.err.println("WARNING: Collisions: "+nVsRateFunc.getNumCollisions());
				nVsRateFunc.set(rate.nFor10Percent, rate.mean);
				nVsMoRateFunc.set(rate.nFor10Percent, moRate);
			}
		}
		
		double rateNorm = normStdDevFunc.getMaxY() / normRateFunc.getMaxY();
		double moRateNorm = normStdDevFunc.getMaxY() / normMoRateFunc.getMaxY();
		normRateFunc.scale(rateNorm);
		normMoRateFunc.scale(moRateNorm);
		
		ArrayList<XY_DataSet> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		funcs.add(normStdDevFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		writePlot(dir, prefix, funcs, chars, "Std Dev / Mean", null, "Std Dev / Mean", false);
		funcs.add(normRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		funcs.add(normMoRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		writePlot(dir, prefix+"_with_rates", funcs, chars, "Std Dev / Mean",
				null, "Std Dev / Mean", false);
		
		funcs.clear();
		chars.clear();
		funcs.add(normStdDevVsRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.RED));
		writePlot(dir, prefix+"_rate_scatter", funcs, chars, "Std Dev / Mean vs Rate",
				"Std Dev / Mean", "Rate", true);
		funcs.clear();
		chars.clear();
		funcs.add(normStdDevVsMoRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.BLACK));
		writePlot(dir, prefix+"_mo_rate_scatter", funcs, chars, "Std Dev / Mean vs Moment Rate",
				"Std Dev / Mean", "Moment Rate", true);
		
		funcs.clear();
		chars.clear();
		funcs.add(nVsRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.RED));
		writePlot(dir, prefix+"_n_vs_rate_scatter", funcs, chars, "N for SDOM within 10% vs Rate",
				"N", "Rate", true);
		funcs.clear();
		chars.clear();
		funcs.add(nVsMoRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.BLACK));
		writePlot(dir, prefix+"_n_vs_mo_rate_scatter", funcs, chars, "N for SDOM within 10% vs Moment Rate",
				"N", "Moment Rate", true);
		
		File csvFile = new File(dir, prefix+"_rates.csv");
		csv.writeToFile(csvFile);
	}
	
	private static void writePlot(File dir, String prefix,
			List<? extends XY_DataSet> funcs, ArrayList<PlotCurveCharacterstics> chars,
			String title, String xAxis, String yAxis, boolean yLog) throws IOException {
		writePlot(dir, prefix, funcs, chars, title, xAxis, yAxis, yLog, null);
	}
	
	private static void writePlot(File dir, String prefix,
			List<? extends XY_DataSet> funcs, ArrayList<PlotCurveCharacterstics> chars,
			String title, String xAxis, String yAxis, boolean yLog, double[][] axis) throws IOException {
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setYLog(yLog);
		CommandLineInversionRunner.setFontSizes(gp);
		
		if (axis != null)
			gp.setUserBounds(axis[0][0], axis[0][1], axis[1][0], axis[1][1]);
		
		gp.drawGraphPanel(xAxis, yAxis, funcs, chars, title);
		
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(new File(dir, prefix+".pdf").getAbsolutePath());
		gp.saveAsPNG(new File(dir, prefix+".png").getAbsolutePath());
		if (funcs.get(0).size() < 2000)
			gp.saveAsTXT(new File(dir, prefix+".txt").getAbsolutePath());
	}
	
	private static class RateEntry implements Comparable<RateEntry> {
		int index;
		double[] vals;
		double mean;
		double stdDev;
		double normStdDev;
		double nFor10Percent;
		public RateEntry(int index, double[] vals) {
			this.index = index;
			this.vals = vals;
			mean = StatUtils.mean(vals);
			stdDev = Math.sqrt(StatUtils.variance(vals, mean));
			normStdDev = stdDev / mean;
			nFor10Percent = Math.pow(10d*stdDev/mean, 2);
		}
		@Override
		public int compareTo(RateEntry o) {
			return Double.compare(o.normStdDev, normStdDev);
		}
	}
	
	private static EvenlyDiscretizedFunc loadNonWaterlevelsForBranches(ZipFile file, String nameGrep)
			throws IOException {
		// first find all noMinRates entries
		List<ZipEntry> noMinsEntries = Lists.newArrayList();
		for (ZipEntry entry : Collections.list(file.entries())) {
			if (!entry.getName().contains(nameGrep))
				continue;
			if (entry.getName().endsWith("_noMinRates.bin"))
				noMinsEntries.add(entry);
		}
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(1, noMinsEntries.size(), 1d);
		
		boolean[] aboves = null;
		for (int i=0; i<noMinsEntries.size(); i++) {
			ZipEntry entry = noMinsEntries.get(i);
			double[] rates =
					MatrixIO.doubleArrayFromInputStream(file.getInputStream(entry), entry.getSize());
			
			if (aboves == null)
				aboves = new boolean[rates.length];
			
			int cnt = 0;
			for (int r=0; r<rates.length; r++) {
				if (aboves[r]) {
					cnt++;
				} else if (rates[r] > 0) {
					aboves[r] = true;
					cnt++;
				}
			}
			
			func.set(i, (double)cnt);
		}
		
		return func;
	}
	
	private static double[] loadWaterlevels(File ratesFile, File ratesNoMinFile) throws IOException {
		double[] rates = MatrixIO.doubleArrayFromFile(ratesFile);
		double[] rateNomin = MatrixIO.doubleArrayFromFile(ratesNoMinFile);
		
		double[] waterlevels = new double[rates.length];
		
		for (int r=0; r<rates.length; r++)
			waterlevels[r] = rates[r] - rateNomin[r];
		
		return waterlevels;
	}

	private static EvenlyDiscretizedFunc loadNonWaterlevelsForAvg(
			AverageFaultSystemSolution avgSol, double[] waterlevels) throws IOException {
//		double[] waterlevels = new double[avgSol.getNumRuptures()];
//		double[] run0NoMins = MatrixIO.doubleArrayFromFile(run0NoMinsFile);
//		double[] run0Rates = avgSol.getRates(0);
//		int run0Above = 0;
//		for (int r=0; r<avgSol.getNumRuptures(); r++) {
//			waterlevels[r] = run0Rates[r] - run0NoMins[r];
//			Preconditions.checkState(waterlevels[r] >= 0);
//			if (run0Rates[r] > run0NoMins[r])
//				run0Above++;
//		}
//		System.out.println("Run0 above: "+run0Above);
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(1, avgSol.getNumSolutions(), 1d);
		
		boolean[] aboves = null;
		for (int i=0; i<avgSol.getNumSolutions(); i++) {
			double[] rates = avgSol.getRates(i);
			
			if (aboves == null)
				aboves = new boolean[rates.length];
			
			int cnt = 0;
			for (int r=0; r<rates.length; r++) {
				if (aboves[r]) {
					cnt++;
				} else if (rates[r] > waterlevels[r]) {
					aboves[r] = true;
					cnt++;
				}
			}
			
			func.set(i, (double)cnt);
		}
		
		return func;
	}
	
	private static void writeWaterlevelResults(File dir, String prefix, EvenlyDiscretizedFunc func, int numRups)
			throws IOException {
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(func);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		double[][] axis = new double[2][2];
		axis[0][0] = 0d;
		axis[0][1] = func.getMaxX();
		axis[1][0] = 0d;
		axis[1][1] = numRups;
		
		writePlot(dir, prefix+"_above_waterlevel", funcs, chars, "Ruptures Above Waterlevel",
				"# Runs", "# Rups Above Waterlevel", false, axis);
	}

}
