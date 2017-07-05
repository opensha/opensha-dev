package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.collections.EnumerationUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sra.calc.parallel.ThreadedEALCalc;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.CalculationExceptionHandler;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.inversion.BatchPlotGen;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ManyRunCompilation {
	
	private static double[][] loadRates(File zipFile, int numRups, int numRuns) throws ZipException, IOException {
		
		double[][] rates = new double[numRups][numRuns];
		
		ZipFile zip = new ZipFile(zipFile);
		List<ZipEntry> entries = EnumerationUtils.toList(zip.entries());
		Preconditions.checkState(entries.size() == numRuns, "zip entry size is off! "+entries.size()+" != "+numRuns);
		
		for (int i=0; i<numRuns; i++) {
			double[] runRates = MatrixIO.doubleArrayFromInputStream(zip.getInputStream(entries.get(i)), numRups*8l);
			
			for (int r=0; r<numRups; r++) {
				rates[r][i] = runRates[r];
			}
		}
		
		return rates;
	}
	
	public static double[][] loadRates(File dir, int numRups) throws IOException {
		ArrayList<File> files = new ArrayList<File>();
		
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				continue;
			if (!file.getName().endsWith(".bin"))
				continue;
			if (file.getName().contains("partic") || file.getName().contains("std_dev"))
				continue;
			
			files.add(file);
		}
		
		int numRuns = files.size();
		
		double[][] rates = new double[numRups][numRuns];
		
		int meanProfiles = 20;
		
		int meanProfileMod = numRuns / meanProfiles;
		
		File ratePlotsDir = new File(dir, "ratePlots");
		ratePlotsDir.mkdir();
		File sortedRatePlotsDir = new File(dir, "sortedRatePlots");
		sortedRatePlotsDir.mkdir();
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		ArrayList<DiscretizedFunc> funcs = new ArrayList<DiscretizedFunc>();
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		String rankXAxisLabel = "Rank";
		String idXAxisLabel = "Rupture ID";
		
		int numDigits = new String(""+(numRuns-1)).length();
		
		String title = "Rupture Rates";
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setYLog(true);
		
		EvenlyDiscretizedFunc meanFunc = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		meanFunc.setName("Mean Rupture Rate");
		EvenlyDiscretizedFunc stdDevFunc = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		stdDevFunc.setName("Std Dev of Rupture Rate");
		EvenlyDiscretizedFunc stdDevOfMeanFunc = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		stdDevOfMeanFunc.setName("Std Dev Of Mean Rupture Rate");
		
		for (int i=0; i<numRuns; i++) {
			double[] runRates = MatrixIO.doubleArrayFromFile(files.get(i));
			Preconditions.checkState(runRates.length == numRups,
					"Rate file is wrong size: "+runRates.length+" != "+numRups
					+" ("+files.get(i).getName()+")");
			
			String name = i+"";
			while (name.length() < numDigits)
				name = "0"+name;
			
			File rankFile = new File(ratePlotsDir, "rates_"+name+".png");
			if (!rankFile.exists()) {
				for (int j=0; j<runRates.length; j++)
					func.set(j, runRates[j]);
				
				gp.drawGraphPanel(idXAxisLabel, "Rate", funcs, chars, title);
				gp.getChartPanel().setSize(1000, 800);
				gp.saveAsPNG(rankFile.getAbsolutePath());
			}
			
			rankFile = new File(sortedRatePlotsDir, "sorted_rates_"+name+".png");
			if (!rankFile.exists()) {
				double[] sortedRates = Arrays.copyOf(runRates, runRates.length);
				Arrays.sort(sortedRates);
				
				int cnt = 0;
				for (int j=sortedRates.length; --j>=0;)
					func.set(cnt++, sortedRates[j]);
				
				gp.drawGraphPanel(rankXAxisLabel, "Rate", funcs, chars, title);
				gp.getChartPanel().setSize(1000, 800);
				gp.saveAsPNG(rankFile.getAbsolutePath());
			}
			
			for (int r=0; r<numRups; r++) {
				try {
					rates[r][i] = runRates[r];
				} catch (RuntimeException e) {
					System.out.println("r: "+r+", i: "+i+", numRups: "+numRups+", numRuns: "+numRuns);
					throw e;
				}
			}
			
			if (i > 0 && (i % meanProfileMod == 0 || i == (numRuns - 1))) {
				File meanFile = new File(ratePlotsDir, "mean_rates_"+name+".png");
				
				if (!meanFile.exists()) {
					int n = i+1;
					for (int r=0; r<numRups; r++) {
						double[] rupRates = Arrays.copyOf(rates[r], n);
						double mean = StatUtils.mean(rupRates);
						double stdDev = Math.sqrt(StatUtils.variance(rupRates, mean));
						double stdDevOfMean = stdDev / Math.sqrt(n);
						
						meanFunc.set(r, mean);
						stdDevFunc.set(r, stdDev);
						stdDevOfMeanFunc.set(r, stdDevOfMean);
					}
					
					double[] sortedRates = Arrays.copyOf(runRates, runRates.length);
					Arrays.sort(sortedRates);
					
					int cnt = 0;
					for (int j=sortedRates.length; --j>=0;)
						func.set(cnt++, sortedRates[j]);
					
					ArrayList<EvenlyDiscretizedFunc> funcs2 = Lists.newArrayList(stdDevFunc, stdDevOfMeanFunc, meanFunc);
					ArrayList<PlotCurveCharacterstics> chars2 = Lists.newArrayList(
							new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE),
							new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN),
							new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
					
					gp.drawGraphPanel(idXAxisLabel, "Rate", funcs2, chars2, "Mean/Std Devs after "+n+" runs");
					gp.getChartPanel().setSize(1000, 800);
					gp.saveAsPNG(meanFile.getAbsolutePath());
				}
			}
		}
		
		return rates;
	}
	
	private static void generateStabilityPlot(double[][] rates, double[] meanRates, File dir, FaultSystemRupSet rupSet) throws IOException {
		int totRuns = rates[0].length;
		
		EvenlyDiscretizedFunc maxFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
		maxFunc.setName("Maximum Individual Residual vs Overall Mean");
		EvenlyDiscretizedFunc totFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
		totFunc.setName("Total Residuals vs Overall Mean");
		EvenlyDiscretizedFunc avgFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
		avgFunc.setName("Average Individual Residual vs Overall Mean");
		EvenlyDiscretizedFunc medFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
		medFunc.setName("Median Individual Residual vs Overall Mean");
		
		int numRups = rates.length;
		double[] runningMeans = new double[numRups];
		
		double maxResidual = 0;
		double totResidual = 0;
		double[] residuals = new double[numRups];
		for (int rup=0; rup<numRups; rup++) {
			double rate = rates[rup][0];
			runningMeans[rup] = rate;
			double residual = Math.abs(rate - meanRates[rup]) / meanRates[rup];
			totResidual += residual;
			if (residual > maxResidual)
				maxResidual = residual;
			residuals[rup] = residual;
		}
		double avgResidual = totResidual / (double)rates.length;
		
		maxFunc.set(0, maxResidual);
		totFunc.set(0, totResidual);
		avgFunc.set(0, avgResidual);
		Arrays.sort(residuals);
		medFunc.set(0, median(residuals));
		
		for (int runs=2; runs<=totRuns; runs++) {
			maxResidual = 0;
			totResidual = 0;
			for (int rup=0; rup<rates.length; rup++) {
				double rate = runningMeans[rup] * (double)(runs-1) + rates[rup][runs-1];
//				for (int r=0; r<runs; r++)
//					rate += rates[rup][r];
				rate /= (double)runs;
				runningMeans[rup] = rate;
				double residual = (Math.abs(rate - meanRates[rup]) / meanRates[rup]);
				totResidual += residual;
				if (residual > maxResidual)
					maxResidual = residual;
				residuals[rup] = residual;
				
				if (runs == totRuns) {
					double diff = Math.abs(rate - meanRates[rup]);
					if (diff > 1e-14)
						throw new IllegalStateException("final rates don't match orig mean! abs("
								+rate+"-"+meanRates[rup]+")="+diff);
				}
			}
			avgResidual = totResidual / (double)rates.length;
			
			maxFunc.set(runs-1, maxResidual);
			totFunc.set(runs-1, totResidual);
			avgFunc.set(runs-1, avgResidual);
			Arrays.sort(residuals);
			double median = median(residuals);
			medFunc.set(runs-1, median);
			
			if (runs % 20 == 0) {
				System.out.println(runs+" runs: max="+maxResidual+"\ttot="+totResidual
						+"\tavg="+avgResidual+"\tmed="+median);
			}
		}
		
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		
		ArrayList<EvenlyDiscretizedFunc> funcs = Lists.newArrayList(maxFunc, totFunc, avgFunc, medFunc);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setYLog(true);
		gp.drawGraphPanel("# Runs", "Normalized Residuals vs Mean", funcs, chars, "Rup Rate Residuals Vs Runs");
		File plotFile = new File(dir, "rup_rate_residuals");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
		
		int[] nVals = { 1, 10, 100, 1000, 10000, 100000 };
		
		ArrayList<String> header = Lists.newArrayList("Rup ID", "Mag", "# Sects", "# Parent Sects",
				"Mean Rate", "Std Dev");
		for (int n : nVals)
			header.add("Mean/SDOM n="+n);
		
		// now csv file
		CSVFile<String> rupRateCSV = new CSVFile<String>(true);
		rupRateCSV.addLine(header);
		
		for (int r=0; r<numRups; r++) {
			double mean = StatUtils.mean(rates[r]);
			double stdDev = Math.sqrt(StatUtils.variance(rates[r], mean));
			
			double mag = rupSet.getMagForRup(r);
			List<Integer> sects = rupSet.getSectionsIndicesForRup(r);
			int numSects = sects.size();
			HashSet<Integer> parents = new HashSet<Integer>();
			for (int sect : sects) {
				int parent = rupSet.getFaultSectionData(sect).getParentSectionId();
				if (!parents.contains(parent))
					parents.add(parent);
			}
			
			ArrayList<String> line = Lists.newArrayList();
			line.add(r+"");
			line.add(mag+"");
			line.add(numSects+"");
			line.add(parents.size()+"");
			line.add(mean+"");
			line.add(stdDev+"");
			for (int n : nVals)
				line.add((mean / (stdDev / Math.sqrt(n)))+"");
			
			rupRateCSV.addLine(line);
			
//			int numRows = rupRateCSV.getNumRows();
//			for (int row=1; row<=numRows; row++) {
//				if (row == numRows) {
//					rupRateCSV.addLine(line);
//				}
//				double otherMoverSD = Double.parseDouble(rupRateCSV.get(row, 6));
//				if (meanOverStdDev >= otherMoverSD || Double.isNaN(otherMoverSD)) {
//					rupRateCSV.addLine(row, line);
//					break;
//				}
//			}
		}
		
		rupRateCSV.sort(6, 1, new Comparator<String>() {
			
			@Override
			public int compare(String o1, String o2) {
				Double d1 = Double.parseDouble(o1);
				Double d2 = Double.parseDouble(o2);
				return Double.compare(d1, d2);
			}
		});
		
		File csvFile = new File(dir, "rupRates.csv");
		rupRateCSV.writeToFile(csvFile);
	}
	
	private static void generateStabilityParticPlot(double[][] rates, double[] meanRates,
			File dir, FaultSystemRupSet rupSet) throws IOException {
		int totRuns = rates[0].length;
		int numRups = rates.length;
		int numSects = rupSet.getNumSections();
		
		ArrayList<EvenlyDiscretizedFunc> maxNormFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> totNormFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> avgNormFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> medNormFuncs = Lists.newArrayList();
		
		ArrayList<EvenlyDiscretizedFunc> maxFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> totFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> avgFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> medFuncs = Lists.newArrayList();
		
		ArrayList<EvenlyDiscretizedFunc> lowestMeanOverStdDevOfMeanFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> avgMeanOverStdDevOfMeanFuncs = Lists.newArrayList();
		ArrayList<EvenlyDiscretizedFunc> medMeanOverStdDevOfMeanFuncs = Lists.newArrayList();
		
		ArrayList<ArbitrarilyDiscretizedFunc> theoreticalNFuncs = Lists.newArrayList();
		ArrayList<ArbitrarilyDiscretizedFunc> percentAvobeMeanOverSDOMFuncs = Lists.newArrayList();
		ArrayList<DefaultXY_DataSet> scatterMeanOverSDOMFuncs = Lists.newArrayList();
		
		double min = 6;
		double max = 8.5;
		double delta = 0.5;
//		double max = 9;
//		double delta = 1;
		
		ArrayList<double[]> ranges = Lists.newArrayList();
		List<double[]> meanParticRatesList = Lists.newArrayList();
		List<double[]> stdDevParticRatesList = Lists.newArrayList();
		for (double start=min; (start+delta)<=max; start+=delta) {
			ranges.add(toArray(start, start+delta));
			
			String str = " ("+start+" <= mag < "+(start+delta)+")";
			
			EvenlyDiscretizedFunc maxFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			maxFunc.setName("Maximum Normalized Individual Partic Residual vs Overall Mean"+str);
			EvenlyDiscretizedFunc totFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			totFunc.setName("Total Normalized Partic Residuals vs Overall Mean"+str);
			EvenlyDiscretizedFunc avgFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			avgFunc.setName("Average Normalized Individual Partic Residual vs Overall Mean"+str);
			EvenlyDiscretizedFunc medFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			medFunc.setName("Median Normalized Individual Partic Residual vs Overall Mean"+str);
			
			maxNormFuncs.add(maxFunc);
			totNormFuncs.add(totFunc);
			avgNormFuncs.add(avgFunc);
			medNormFuncs.add(medFunc);
			
			maxFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			maxFunc.setName("Maximum Individual Partic Residual vs Overall Mean"+str);
			totFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			totFunc.setName("Total Partic Residuals vs Overall Mean"+str);
			avgFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			avgFunc.setName("Average Individual Partic Residual vs Overall Mean"+str);
			medFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			medFunc.setName("Median Individual Partic Residual vs Overall Mean"+str);
			
			EvenlyDiscretizedFunc meanOverStdDevOfMeanFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			meanOverStdDevOfMeanFunc.setName("Lowest Mean / Std. Dev Of Mean"+str);
			lowestMeanOverStdDevOfMeanFuncs.add(meanOverStdDevOfMeanFunc);
			EvenlyDiscretizedFunc avgOverStdDevOfMeanFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			avgOverStdDevOfMeanFunc.setName("Avg Mean / Std. Dev Of Mean"+str);
			avgMeanOverStdDevOfMeanFuncs.add(avgOverStdDevOfMeanFunc);
			EvenlyDiscretizedFunc medOverStdDevOfMeanFunc = new EvenlyDiscretizedFunc(1d, totRuns, 1d);
			medOverStdDevOfMeanFunc.setName("Median Mean / Std. Dev Of Mean"+str);
			medMeanOverStdDevOfMeanFuncs.add(medOverStdDevOfMeanFunc);
			
			ArbitrarilyDiscretizedFunc theoreticalNFunc = new ArbitrarilyDiscretizedFunc();
			theoreticalNFunc.setName("Theoretical N to reach target Mean / Std. Dev Of Mean"+str);
			theoreticalNFuncs.add(theoreticalNFunc);
			
			ArbitrarilyDiscretizedFunc percentAvobeMeanOverSDOMFunc = new ArbitrarilyDiscretizedFunc();
			percentAvobeMeanOverSDOMFunc.setName("% sects above targen Mean / Std. Dev Of Mean"+str);
			percentAvobeMeanOverSDOMFuncs.add(percentAvobeMeanOverSDOMFunc);
			
			DefaultXY_DataSet scatterMeanOverSDOMFunc = new DefaultXY_DataSet();
			scatterMeanOverSDOMFunc.setName("Mean Rate vs Mean / Std. Dev Of Mean"+str);
			scatterMeanOverSDOMFuncs.add(scatterMeanOverSDOMFunc);
			
			maxFuncs.add(maxFunc);
			totFuncs.add(totFunc);
			avgFuncs.add(avgFunc);
			medFuncs.add(medFunc);
			
			meanParticRatesList.add(new double[numSects]);
			stdDevParticRatesList.add(new double[numSects]);
		}
		
		Map<Integer, List<List<Integer>>> particRupsMap = Maps.newHashMap();
		
		for (int sect=0; sect<numSects; sect++) {
			List<List<Integer>> particRups = particRupsMap.get(sect);
			if (particRups == null) {
				particRups = Lists.newArrayList();
				for (int i=0; i<ranges.size(); i++)
					particRups.add(new ArrayList<Integer>());
				particRupsMap.put(sect, particRups);
			}
			List<Integer> rupsForSect = rupSet.getRupturesForSection(sect);
			for (int rup : rupsForSect) {
				double mag = rupSet.getMagForRup(rup);
				for (int i=0; i<ranges.size(); i++) {
					double[] range = ranges.get(i);
					if (mag >= range[0] && mag < range[1]) {
						particRups.get(i).add(rup);
						meanParticRatesList.get(i)[sect] += meanRates[rup];
						break;
					}
				}
			}
		}
		
		for (int i=0; i<ranges.size(); i++) {
			double[] stdDevs = stdDevParticRatesList.get(i);
			double[] means = meanParticRatesList.get(i);
			
			for (int sect=0; sect<numSects; sect++) {
				double[] particRates = new double[totRuns];
				for (int run=0; run<totRuns; run++) {
					//					if (meanParticRates[sect] < 1e-5)
					//						continue;
					particRates[run] = 0;
					List<Integer> rups = particRupsMap.get(sect).get(i);
					for (int rup : rups)
						particRates[run] += rates[rup][run];
				}
				double mean = StatUtils.mean(particRates);
				stdDevs[sect] = Math.sqrt(StatUtils.variance(particRates, mean));
				Preconditions.checkState(DataUtils.getPercentDiff(means[sect], mean) < 0.01, "uh oh, mean is off!");
			}
		}
		
		double[] runningMeans = new double[numRups];
		
		for (int rup=0; rup<numRups; rup++) {
			double rate = rates[rup][0];
			runningMeans[rup] = rate;
		}
		
		for (int runs=1; runs<=totRuns; runs++) {
			if (runs > 1) {
				for (int rup=0; rup<rates.length; rup++) {
					double rate = runningMeans[rup] * (double)(runs-1) + rates[rup][runs-1];
					rate /= (double)runs;
					runningMeans[rup] = rate;
				}
			}
			
			for (int i=0; i<ranges.size(); i++) {
				double[] meanParticRates = meanParticRatesList.get(i);
				double[] stdDevParticRates = stdDevParticRatesList.get(i);
				
				double maxNormalizedResidual = 0;
				int maxNormalizedIndex = 0;
				double maxPartic = 0;
				double totNormalizedResidual = 0;
				
				double maxResidual = 0;
				int maxIndex = 0;
				double totResidual = 0;
				
				double[] meanOverStdDevOfMeans = new double[numSects];
				
				ArrayList<Double> normalizedResiduals = new ArrayList<Double>();
				ArrayList<Double> residuals = new ArrayList<Double>();
				
				for (int sect=0; sect<numSects; sect++) {
//					if (meanParticRates[sect] < 1e-5)
//						continue;
					double particRate = 0;
					List<Integer> rups = particRupsMap.get(sect).get(i);
					for (int rup : rups)
						particRate += runningMeans[rup];
					
					double residual = Math.abs(particRate - meanParticRates[sect]);
					double normalizedResidual;
					if (meanParticRates[sect] == 0) {
						if (particRate == 0)
							normalizedResidual = 0;
						else
							normalizedResidual = 1e3;
					} else {
						normalizedResidual = residual / meanParticRates[sect];
					}
					
					if (runs == totRuns) {
						double diff = Math.abs(particRate - meanParticRates[sect]);
						if (diff > 1e-14)
							throw new IllegalStateException("final rates don't match orig mean! abs("
									+particRate+"-"+meanParticRates[sect]+")="+diff);
					}
					
					totResidual += residual;
					if (residual > maxResidual) {
						maxResidual = residual;
						maxIndex = sect;
					}
					residuals.add(residual);
					
					totNormalizedResidual += normalizedResidual;
					if (normalizedResidual > maxNormalizedResidual) {
						maxNormalizedResidual = normalizedResidual;
						maxNormalizedIndex = sect;
						maxPartic = particRate;
					}
					normalizedResiduals.add(normalizedResidual);
					
					meanOverStdDevOfMeans[sect] = meanParticRates[sect] / (stdDevParticRates[sect] / Math.sqrt((double)runs));
				}
				double avgNormResidual = totNormalizedResidual / (double)normalizedResiduals.size();
				
				maxNormFuncs.get(i).set(runs-1, maxNormalizedResidual);
				totNormFuncs.get(i).set(runs-1, totNormalizedResidual);
				avgNormFuncs.get(i).set(runs-1, avgNormResidual);
				Collections.sort(normalizedResiduals);
				double median = median(normalizedResiduals);
				medNormFuncs.get(i).set(runs-1, median);
				
				double avgResidual = totResidual / (double)residuals.size();
				
				maxFuncs.get(i).set(runs-1, maxResidual);
				totFuncs.get(i).set(runs-1, totResidual);
				avgFuncs.get(i).set(runs-1, avgResidual);
				Collections.sort(residuals);
				median = median(residuals);
				medFuncs.get(i).set(runs-1, median);
				
				lowestMeanOverStdDevOfMeanFuncs.get(i).set(runs-1, StatUtils.min(meanOverStdDevOfMeans));
				avgMeanOverStdDevOfMeanFuncs.get(i).set(runs-1, StatUtils.mean(meanOverStdDevOfMeans));
				medMeanOverStdDevOfMeanFuncs.get(i).set(runs-1, median(meanOverStdDevOfMeans));
				
				if (runs % 20 == 0) {
					System.out.println("RANGE: "+ranges.get(i)[0]+" => "+ranges.get(i)[1]);
					System.out.println(runs+" runs: max="+maxNormalizedResidual+"\ttot="+totNormalizedResidual
							+"\tavg="+avgNormResidual+"\tmed="+median);
					System.out.println("MAX at "+maxNormalizedIndex+": abs("+maxPartic+"-"
							+meanParticRates[maxNormalizedIndex]+")="+Math.abs(maxPartic-meanParticRates[maxNormalizedIndex]));
				}
			}
		}
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, ranges.size()-1);
		
		ArrayList<EvenlyDiscretizedFunc> funcs = Lists.newArrayList();
		funcs.addAll(maxNormFuncs);
		funcs.addAll(totNormFuncs);
		funcs.addAll(avgNormFuncs);
		funcs.addAll(medNormFuncs);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> scatterChars = Lists.newArrayList();
		for (int j=0; j<maxNormFuncs.size(); j++)
			scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.X, 1f, cpt.getColor((float)j)));
		for (int j=0; j<maxNormFuncs.size(); j++)
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt.getColor((float)j)));
		for (int j=0; j<maxNormFuncs.size(); j++)
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, cpt.getColor((float)j)));
		for (int j=0; j<maxNormFuncs.size(); j++)
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, cpt.getColor((float)j)));
		for (int j=0; j<maxNormFuncs.size(); j++)
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, cpt.getColor((float)j)));
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setYLog(true);
		
		gp.drawGraphPanel("# Runs", "Normalized Residuals vs Mean", funcs, chars, "Normalized Partic Rate Residuals Vs Runs");
		File plotFile = new File(dir, "partic_rate_normalized_residuals");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
		
		funcs.clear();
		funcs.addAll(maxFuncs);
		funcs.addAll(totFuncs);
		funcs.addAll(avgFuncs);
		funcs.addAll(medFuncs);
		
		gp.drawGraphPanel("# Runs", "Residuals vs Mean", funcs, chars, "Partic Rate Residuals Vs Runs");
		plotFile = new File(dir, "partic_rate_residuals");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
		
		// move the dashed ones up a touch
		for (int i=0; i<lowestMeanOverStdDevOfMeanFuncs.size(); i++)
			chars.remove(lowestMeanOverStdDevOfMeanFuncs.size());
		
		funcs.clear();
		funcs.addAll(lowestMeanOverStdDevOfMeanFuncs);
//		funcs.addAll(avgMeanOverStdDevOfMeanFuncs);
//		funcs.addAll(medMeanOverStdDevOfMeanFuncs);
		
		gp.setYLog(false);
		gp.drawGraphPanel("N", "Mean / Std Dev Of Mean", funcs, chars,
				"Lowest Mean / Std Dev Of Mean vs # Runs");
		plotFile = new File(dir, "partic_mean_over_std_dev_of_mean");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
		
		for (int i=0; i<ranges.size(); i++) {
			double[] meanParticRates = meanParticRatesList.get(i);
			double[] stdDevParticRates = stdDevParticRatesList.get(i);
			
			for (double targetVal=1; targetVal<=100; targetVal += 0.5) {
				double worstN = 0;
				
				for (int sect=0; sect<numSects; sect++) {
					if (particRupsMap.get(sect).get(i).isEmpty())
						continue;
					double n = Math.pow((targetVal * stdDevParticRates[sect]) / meanParticRates[sect], 2);
					if (n > worstN)
						worstN = n;
				}
				
				theoreticalNFuncs.get(i).set(worstN, targetVal);
			}
		}
		
		gp.drawGraphPanel("N", "Target Mean / Std Dev Of Mean", theoreticalNFuncs, chars,
				"Theoretical N Needed to reach given Mean / Std Dev Of Mean");
		plotFile = new File(dir, "partic_target_mean_over_std_dev_of_mean");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
		
		double[] aboveTargets = { 20d, 50d, 100d };
		
		for (double aboveTarget : aboveTargets) {
			for (int i=0; i<ranges.size(); i++) {
				double[] meanParticRates = meanParticRatesList.get(i);
				double[] stdDevParticRates = stdDevParticRatesList.get(i);
				
				for (int runs=1; runs<=1000; runs++) {
					double above = 0;
					double num = 0;
					for (int sect=0; sect<numSects; sect++) {
						if (particRupsMap.get(sect).get(i).isEmpty())
							continue;
//						if (stdDevParticRates[sect] == 0 || meanParticRates[sect] == 0)
//							System.out.println("We have a zero!!!! (mean="+meanParticRates[sect]+", stdDev="+stdDevParticRates[sect]+")");
						double meanOverSDOM = meanParticRates[sect] / (stdDevParticRates[sect] / Math.sqrt((double)runs));
						if (meanOverSDOM >= aboveTarget || (meanParticRates[sect] == 0))
							above++;
						
						if (runs == 1)
							scatterMeanOverSDOMFuncs.get(i).set(meanParticRates[sect], meanOverSDOM);
						num++;
					}
					above = (above / num) * 100d;
//					System.out.println(runs+", "+above);
					percentAvobeMeanOverSDOMFuncs.get(i).set((double)runs, above);
				}
			}
			
			gp.drawGraphPanel("N", "% above Mean / Std Dev Of Mean of "+aboveTarget, percentAvobeMeanOverSDOMFuncs, chars,
					"% Sects above Mean/SDOM of "+aboveTarget);
			plotFile = new File(dir, "partic_above_"+(int)aboveTarget+"_mean_over_std_dev_of_mean");
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
			gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
			

			if (aboveTarget == aboveTargets[0]) {
				gp.setXLog(true);
				gp.setYLog(true);
				gp.setUserBounds(1e-10, 1e-2, 1e-1, 1e3);
				gp.drawGraphPanel("Mean Rate", "Mean / Std Dev", scatterMeanOverSDOMFuncs, scatterChars,
						"Mean Rate vs Mean/SDOM of N=1");
				plotFile = new File(dir, "partic_scatter_mean_over_std_dev_of_mean");
				gp.getChartPanel().setSize(1000, 800);
				gp.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
				gp.saveAsPNG(plotFile.getAbsolutePath()+".png");
				gp.setXLog(false);
				gp.setYLog(false);
				gp.setUserBounds(null, null);
			}
		}
		
//		System.exit(0);
	}
	
	private static void generateStabilityEALPlot(double[][] rates, double[] meanRates,
			File dir, FaultSystemRupSet rupSet) throws IOException, InterruptedException {
//		File portFile =new File("/home/kevin/OpenSHA/portfolio_lec/" +
//				"Porter (07 NOv 2011) Portfolio LEC CEA proxy portfolio.csv");
//		Portfolio portfolio = Portfolio.createPortfolio(portFile);
//		
//		int numThreads = Runtime.getRuntime().availableProcessors();
//		
//		ScalarIMR[] imrs = new ScalarIMR[numThreads];
//		List<Asset> assets = portfolio.getAssetList();
//		double maxSourceDistance = 200;
//		CalculationExceptionHandler handler = new CalculationExceptionHandler() {
//			
//			@Override
//			public void calculationException(String errorMessage) {
//				System.out.println("ERROR: "+errorMessage);
//				System.exit(1);
//			}
//		};
//		
//		Site[] sites = new Site[numThreads];
//		
//		for (int i=0; i<numThreads; i++) {
//			imrs[i] = new CB_2008_AttenRel(null);
//			imrs[i].setParamDefaults();
//			
//			sites[i] = new Site(new Location(34, -118));
//			for (Parameter<?> param : imrs[i].getSiteParams())
//				sites[i].addParameter(param);
//		}
//		
//		int numRuns = rates[0].length;
//		int numAssets = assets.size();
//		
//		FaultSystemSolution meanSol = new SimpleFaultSystemSolution(rupSet, meanRates);
//		FaultSystemSolutionPoissonERF meanERF = new FaultSystemSolutionPoissonERF(meanSol);
//		meanERF.updateForecast();
//		
//		ERF[] erfs = {meanERF};
//		
//		double[] meanEALs = calcEALs(assets, erfs, imrs, handler, maxSourceDistance);
//		
//		for (int runs=1; runs<=numRuns; runs++) {
////			FaultSystemSolution sol = new SimpleFaultSystemSolution(rupSet, rupRateSolution)
//		}
	}
	
	private static double[] calcEALs(List<Asset> assets, ERF[] erfs, ScalarIMR[] imrs,
			CalculationExceptionHandler handler, double maxSourceDistance) throws InterruptedException {
		ThreadedEALCalc calc = new ThreadedEALCalc(assets, erfs, imrs, handler, maxSourceDistance, null);
		
		int[] batch = new int[assets.size()];
		for (int i=0; i<batch.length; i++)
			batch[i] = i;
		
		double[] eals = calc.calculateBatch(batch);
		return eals;
	}
	
	private static class RateRecord implements Comparable<RateRecord> {
		
		double mean, min, max, stdDev, lower, upper;
		
		public RateRecord(double mean, double min, double max, double lower, double upper, double stdDev) {
			this.mean = mean;
			this.min = min;
			this.max = max;
			this.stdDev = stdDev;
			this.lower = lower;
			this.upper = upper;
		}

		@Override
		public int compareTo(RateRecord o) {
			return Double.compare(mean, o.mean);
		}
		
	}
	
	private static double median(double[] sorted) {
		if (sorted.length % 2 == 1)
			return sorted[(sorted.length+1)/2-1];
		else
		{
			double lower = sorted[sorted.length/2-1];
			double upper = sorted[sorted.length/2];

			return (lower + upper) * 0.5;
		}	
	}
	
	private static double median(List<Double> sorted) {
		if (sorted.size() % 2 == 1)
			return sorted.get((sorted.size()+1)/2-1);
		else
		{
			double lower = sorted.get(sorted.size()/2-1);
			double upper = sorted.get(sorted.size()/2);

			return (lower + upper) * 0.5;
		}	
	}

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws ZipException 
	 * @throws RuntimeException 
	 * @throws GMT_MapException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException, GMT_MapException, RuntimeException, InterruptedException {
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_03_28-unconstrained-run-like-crazy/results");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_23-reprod_test/results/char");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_23-reprod_test/results/gr");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_26-fm2-reprod_test/results");
		
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_29-fm2-a-priori-test/results/VarAPrioriZero_VarAPrioriWt100");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_29-fm2-a-priori-test/results/VarAPrioriZero_VarAPrioriWt1000");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_29-fm2-a-priori-test/results/VarNone_VarAPrioriWt100");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_29-fm2-a-priori-test/results/VarNone_VarAPrioriWt1000");
		
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_30-fm2-a-priori-test/results/VarAPrioriZero_VarAPrioriWt1000_VarWaterlevel0");
		
		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_05_02-fm2-cooling-tests/results");
		File rupSetFile = new File(dir, "rupSet.zip");
		InversionFaultSystemRupSet rupSet = FaultSystemIO.loadInvRupSet(rupSetFile);
		int numRups = rupSet.getNumRuptures();
		System.out.println("Loaded rupSet with "+numRups+" ruptures");
//		int numRuns = 460;
//		File zipFile = new File(dir, "FM3_1_GLpABM_MaEllB_DsrTap_DrEllB_Unconst_VarAseis0.2_VarOffAseis0.5_VarMFDMod1_VarNone_Manyruns.zip");
//		String prefix = zipFile.getName().replaceAll(".zip", "");
		double[][] rates = loadRates(dir, numRups);
		int numRuns = rates[0].length;
//		String prefix = "FM3_1_GLpABM_MaEllB_DsrTap_DrEllB_Unconst_VarAseis0.1_VarOffAseis0.5_VarMFDMod1_VarNone";
//		String prefix = "FM3_1_GLpABM_MaEllB_DsrTap_DrEllB_Char_VarAseis0.1_VarOffAseis0.5_VarMFDMod1_VarNone";
//		String prefix = "FM3_1_GLpABM_MaEllB_DsrTap_DrEllB_GR_VarAseis0.1_VarOffAseis0.5_VarMFDMod1_VarNone";
//		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarAseis0.1_VarOffAseis0.5_VarMFDMod1_VarNone";
		
//		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarAPrioriZero_VarAPrioriWt100";
//		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarAPrioriZero_VarAPrioriWt1000";
//		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarNone_VarAPrioriWt100";
//		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarNone_VarAPrioriWt1000";
		
//		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarAPrioriZero_VarAPrioriWt1000_VarWaterlevel0";
		
		String prefix = "FM2_1_UC2ALL_MaAvU2_DsrTap_DrAveU2_Char_VarAPrioriZero_VarAPrioriWt1000_VarWaterlevel0";
		System.out.println("Loaded rates!");
		
		RateRecord[] rateRecords = new RateRecord[numRups];
		
		int numZeros = 0;
		double numZerosPer = 0;
		
		int upperInd = (int)(numRuns * 0.975+0.5);
		if (upperInd >= numRuns)
			upperInd = numRuns-1;
		int lowerInd = (int)(numRuns * 0.025+0.5);
		
		double[] meanRates = new double[numRups];
		double[] medianRates = new double[numRups];
		for (int r=0; r<numRups; r++) {
			double[] rupRates = rates[r];
			double mean = StatUtils.mean(rupRates);
			for (double rate : rupRates)
				if (rate == 0)
					numZerosPer++;
			meanRates[r] = mean;
			if (mean == 0)
				numZeros++;
			double min = StatUtils.min(rupRates);
			double max = StatUtils.max(rupRates);
			double stdDev = Math.sqrt(StatUtils.variance(rupRates, mean));
			
			double[] sorted = Arrays.copyOf(rupRates, numRuns);
			Arrays.sort(sorted);
			
			medianRates[r] = median(sorted);
			
//			highRates[r] = mean + stdDev;
			double upper = sorted[upperInd];
//			lowRates[r] = mean - stdDev;
			double lower = sorted[lowerInd];
			rateRecords[r] = new RateRecord(mean, min, max, lower, upper, stdDev);
		}
		
		numZerosPer /= numRuns;
		System.out.println("num zeros: "+numZeros);
		System.out.println("avg zeros per run: "+numZerosPer);
		
		System.out.println("Generating residuals plot...");
		generateStabilityPlot(rates, meanRates, dir, rupSet);
		System.out.println("DONE");
		
		System.out.println("Generating partic residuals plot...");
		generateStabilityParticPlot(rates, meanRates, dir, rupSet);
		System.out.println("DONE");
		
//		System.out.println("Generating EAL residuals plot...");
//		generateStabilityEALPlot(rates, meanRates, dir, rupSet);
//		System.out.println("DONE");
		
		RateRecord[] sortedRecords = Arrays.copyOf(rateRecords, rateRecords.length);
		
		Arrays.sort(sortedRecords);
		
		EvenlyDiscretizedFunc meanFunc = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		EvenlyDiscretizedFunc stdDevFunc = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		
		EvenlyDiscretizedFunc meanFuncByID = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		EvenlyDiscretizedFunc stdDevFuncByID = new EvenlyDiscretizedFunc(0d, numRups, 1d);
		
		int cnt = 0;
		for (int i=sortedRecords.length; --i>=0;) {
			RateRecord rec = sortedRecords[i];
			meanFunc.set(cnt, rec.mean);
			stdDevFunc.set(cnt, rec.stdDev);
			
			meanFuncByID.set(cnt, rateRecords[cnt].mean);
			stdDevFuncByID.set(cnt, rateRecords[cnt].stdDev);
			
			cnt++;
		}
		
		ArrayList<DiscretizedFunc> funcs = new ArrayList<DiscretizedFunc>();
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<PlotCurveCharacterstics>();
		
		funcs.add(stdDevFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1, Color.GREEN));
		funcs.add(meanFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1, Color.BLACK));
		
		String title = "Rupture Rate Distribution";
		
//		GraphWindow gw = new GraphWindow(funcs, title, chars);
//		gw.setX_AxisLabel(xAxisLabel);
//		gw.setYLog(true);
//		gw.getGraphWindow().setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setYLog(true);
		gp.drawGraphPanel("Mean Rate Rank", "Rate", funcs, chars, title);
		File rankFile = new File(dir, prefix+"_rate_dist");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(rankFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(rankFile.getAbsolutePath()+".png");
		
		funcs.set(1, meanFuncByID);
		funcs.set(0, stdDevFuncByID);
		
		gp.drawGraphPanel("Rupture ID", "Rate", funcs, chars, title);
		rankFile = new File(dir, prefix+"_rate_id_dist");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(rankFile.getAbsolutePath()+".pdf");
		gp.saveAsPNG(rankFile.getAbsolutePath()+".png");
		
		InversionFaultSystemSolution meanSol = new InversionFaultSystemSolution(rupSet, meanRates);
		FaultSystemIO.writeSol(meanSol, new File(dir, prefix+"_mean_sol.zip"));
		CommandLineInversionRunner.writeMFDPlots(meanSol, dir, prefix);
		
		ArrayList<PaleoRateConstraint> paleoConstraints = CommandLineInversionRunner.getPaleoConstraints(
				meanSol.getRupSet().getFaultModel(), meanSol.getRupSet());
		CommandLineInversionRunner.writePaleoPlots(paleoConstraints, null, meanSol, dir, prefix+"_mean");
		
		BatchPlotGen.makeMapPlots(meanSol, dir, prefix+"_mean");
		
		// now make the std dev plots
		int numSects = rupSet.getNumSections();
		
		Region region = new CaliforniaRegions.RELM_TESTING();
//		for (int i=0; i<numRuns; i++) {
//			if (i % 10 == 0)
//				System.out.println("Getting slip rates for solution "+i+"/"+numRuns);
//			if (i % 25 == 0)
//				System.gc();
//			double[] myRates = new double[numRups];
//			for (int r=0; r<numRups; r++)
//				myRates[r] = rates[r][i];
//			FaultSystemSolution mySol = new SimpleFaultSystemSolution(rupSet, myRates);
//			mySol.copyCacheFrom(rupSet);
//			double[] mySlipRates = mySol.calcSlipRateForAllSects();
//			
//			for (int s=0; s<numSects;s++) {
//				slipRates[s][i] = mySlipRates[s];
//			}
//		}
		
		System.out.println("Making participation plots...");
		
		rupSet.getRupturesForSection(0); // this initializes the cache
		
		ArrayList<double[]> ranges = new ArrayList<double[]>();
		ranges.add(toArray(6, 7));
		ranges.add(toArray(7, 8));
		ranges.add(toArray(8, 10));
		ranges.add(toArray(6.7, 10));
		
		for (double[] range : ranges) {
			double magLow = range[0];
			double magHigh = range[1];
			
			System.out.println("Range: "+magLow+"=>"+magHigh);
			
			double[][] partRates = new double[numSects][numRuns];
			
			try {
				AverageFaultSystemSolution.calcThreaded(rates, partRates, true, magLow, magHigh, rupSet);
			} catch (InterruptedException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
			
//			for (int i=0; i<numRuns; i++) {
//				if (i % 10 == 0)
//					System.out.println("Getting partic rates for solution "+i+"/"+numRuns);
//				if (i % 25 == 0)
//					System.gc();
//				double[] myRates = new double[numRups];
//				for (int r=0; r<numRups; r++)
//					myRates[r] = rates[r][i];
//				FaultSystemSolution mySol = new SimpleFaultSystemSolution(rupSet, myRates);
//				mySol.copyCacheFrom(rupSet);
//				double[] myPartRates = mySol.calcParticRateForAllSects(magLow, magHigh);
//				
//				for (int s=0; s<numSects;s++) {
//					partRates[s][i] = myPartRates[s];
//				}
//			}
			
			FaultBasedMapGen.plotParticipationStdDevs(rupSet, partRates, region, dir, prefix, false, magLow, magHigh);
		}
		
		System.out.println("Making slip std dev plot");
		double[][] slipRates = new double[numSects][numRuns];
		
		try {
			AverageFaultSystemSolution.calcThreaded(rates, slipRates, false, 0, 0, rupSet);
		} catch (InterruptedException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		FaultBasedMapGen.plotSolutionSlipRateStdDevs(rupSet, slipRates, region, dir, prefix, false);
		System.exit(0);
	}
	
	private static double[] toArray(double... vals) {
		return vals;
	}

}
