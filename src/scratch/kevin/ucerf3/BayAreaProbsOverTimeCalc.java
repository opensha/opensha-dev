package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FSSRupsInRegionCache;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class BayAreaProbsOverTimeCalc {
	
	private static EvenlyDiscretizedFunc[] calc(FaultSystemSolution sol, Region reg, int startYear, int endYear, int yearDelta,
			double minMag, int numMag, double deltaMag, double duration) {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.setParameter(HistoricOpenIntervalParam.NAME,
				(double)(FaultSystemSolutionERF.START_TIME_DEFAULT-1875));
		erf.setParameter(BPTAveragingTypeParam.NAME,
				BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
		erf.getTimeSpan().setDuration(duration);
		
		FSSRupsInRegionCache rupInRegionCache = new FSSRupsInRegionCache();
		
		EvenlyDiscretizedFunc[] ret = new EvenlyDiscretizedFunc[numMag];
		double min = startYear;
		double delta = yearDelta;
		int num = (int)Math.floor((double)(endYear - startYear)/(delta) + 1);
		
		for (int i=0; i<ret.length; i++)
			ret[i] = new EvenlyDiscretizedFunc(min, num, delta);
		
		System.out.println("Beginning calculation");
		
		for (int i=ret[0].size(); --i>=0;) {
			int year = (int) ret[0].getX(i);
			
			System.out.println(year);
			
			erf.getTimeSpan().setStartTime(year);
			erf.updateForecast();
			
			// add a half mag delta so that the cumulative rates line up with original magnitudes
			EvenlyDiscretizedFunc mfd = ERF_Calculator.getParticipationMagFreqDistInRegion(
					erf, reg, minMag+0.5*deltaMag, numMag, deltaMag, true, rupInRegionCache).getCumRateDistWithOffset();
			EvenlyDiscretizedFunc probs = FaultSysSolutionERF_Calc.calcProbsFromSummedMFD(mfd, duration);
			Preconditions.checkState(probs.size() == ret.length);
			for (int m=0; m<probs.size(); m++)
				ret[m].set(i, probs.getY(m));
		}
		
		System.out.println("DONE");
		
		return ret;
	}
	
	private static class CalcThread implements Task {
		
		private FaultSystemSolution sol;
		private Region reg;
		private int startYear;
		private int endYear;
		private int yearDelta;
		private double minMag;
		private int numMag;
		private double deltaMag;
		private double duration;
		
		private EvenlyDiscretizedFunc[] result;

		public CalcThread(FaultSystemSolution sol, Region reg, int startYear, int endYear, int yearDelta, double minMag,
				int numMag, double deltaMag, double duration) {
			super();
			this.sol = sol;
			this.reg = reg;
			this.startYear = startYear;
			this.endYear = endYear;
			this.yearDelta = yearDelta;
			this.minMag = minMag;
			this.numMag = numMag;
			this.deltaMag = deltaMag;
			this.duration = duration;
		}

		@Override
		public void compute() {
			result = calc(sol, reg, startYear, endYear, yearDelta, minMag, numMag, deltaMag, duration);
		}
		
	}
	
	private static void writeCSV(EvenlyDiscretizedFunc[] avgFuncs, double[] mags, File outputFile) throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList("Year");
		for (double mag : mags)
			header.add("Prob M≥"+(float)mag);
		csv.addLine(header);
		for (int i=0; i<avgFuncs[0].size(); i++) {
			List<String> line = Lists.newArrayList("");
			for (int j=0; j<mags.length; j++)
				line.add("");
			csv.addLine(line);
		}
		
		for (int i=0; i<avgFuncs.length; i++) {
			EvenlyDiscretizedFunc avgFunc = avgFuncs[i];
			
			for (int y=0; y<avgFunc.size(); y++) {
				int year = (int)avgFunc.getX(y);
				int row = y+1;
				int col = i+1;
				if (i == 0)
					csv.set(row, 0, year+"");
				csv.set(row, col, avgFunc.getY(y)+"");
			}
		}
		
		csv.writeToFile(outputFile);
	}
	
	public static EvenlyDiscretizedFunc[] loadCSV(File csvFile) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		int startYear = Integer.parseInt(csv.get(1, 0));
		int nextYear = Integer.parseInt(csv.get(2, 0));
		double yearDelta = nextYear - startYear;
		int numYears = csv.getNumRows()-1;
		
		EvenlyDiscretizedFunc[] ret = new EvenlyDiscretizedFunc[csv.getNumCols()-1];
		
		for (int i=0; i<ret.length; i++)
			ret[i] = new EvenlyDiscretizedFunc(startYear, numYears, yearDelta);
		
		for (int row=1; row<csv.getNumRows(); row++) {
			double year = Double.parseDouble(csv.get(row, 0));
			
			int yearIndex = row-1;
			double yearTest = ret[0].getX(yearIndex);
			Preconditions.checkState(yearTest == year, "Year mismatch! %s != %s", yearTest, year);
			
			for (int col=1; col<csv.getNumCols(); col++) {
				int magIndex = col-1;
				ret[magIndex].set(yearIndex, Double.parseDouble(csv.get(row, col)));
			}
		}
		
		return ret;
	}

	public static void main(String[] args) throws IOException, DocumentException, InterruptedException {
		double minMag = 6.7d;
		int numMag = 14;
		double deltaMag = 0.1;
		EvenlyDiscretizedFunc magsFunc = new EvenlyDiscretizedFunc(minMag, numMag, deltaMag);
		double[] mags = new double[magsFunc.size()];
		for (int i=0; i<magsFunc.size(); i++)
			mags[i] = magsFunc.getX(i);
		
		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/bay_area_fact_sheet/");
		
		boolean laComparison = true;
		
		String regName;
		Region reg;
		if (laComparison) {
			reg = new CaliforniaRegions.LA_BOX();
			mainDir = new File(mainDir, "la_compare");
			regName = "LA";
		} else {
			reg = new CaliforniaRegions.SF_BOX();
			regName = "SF";
		}
		int startYear = 1850;
		int endYear = 2014;
		int yearDelta = 1;
//		int endYear = 2010;
//		int yearDelta = 20;
		double duration = 30;
		
		File outputDir = new File(mainDir, "probs_over_time_start"+startYear+"_every"+yearDelta+"_end"+endYear);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File csvFile = new File(outputDir, "table.csv");
		
		EvenlyDiscretizedFunc[] avgResults;
		if (csvFile.exists()) {
			System.out.println("Loading from CSV");
			avgResults = loadCSV(csvFile);
		} else {
			System.out.println("Calculating CSV");
			File fssDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions");

			FaultSystemSolution fm31Sol = FaultSystemIO.loadSol(new File(fssDir,
					"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
			FaultSystemSolution fm32Sol = FaultSystemIO.loadSol(new File(fssDir,
					"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip"));
			
//			EvenlyDiscretizedFunc[] results31 = calc(fm31Sol, reg, startYear, endYear, yearDelta, minMag, numMag, deltaMag, duration);
//			EvenlyDiscretizedFunc[] results32 = calc(fm32Sol, reg, startYear, endYear, yearDelta, minMag, numMag, deltaMag, duration);
			
			CalcThread calc31 = new CalcThread(fm31Sol, reg, startYear, endYear, yearDelta, minMag, numMag, deltaMag, duration);
			CalcThread calc32 = new CalcThread(fm32Sol, reg, startYear, endYear, yearDelta, minMag, numMag, deltaMag, duration);
			
			List<CalcThread> tasks = Lists.newArrayList();
			tasks.add(calc31);
			tasks.add(calc32);
			
			new ThreadedTaskComputer(tasks).computeThreaded(2);
			
			EvenlyDiscretizedFunc[] results31 = calc31.result;
			EvenlyDiscretizedFunc[] results32 = calc32.result;
			
			avgResults = new EvenlyDiscretizedFunc[results31.length];
			for (int i=0; i<results31.length; i++) {
				EvenlyDiscretizedFunc result31 = results31[i];
				EvenlyDiscretizedFunc result32 = results32[i];
				avgResults[i] = new EvenlyDiscretizedFunc(result31.getMinX(), result31.getMaxX(), result31.size());
				
				for (int y=0; y<result31.size(); y++) {
					double avgVal = result31.getY(y)*0.5 + result32.getY(y)*0.5;
					avgResults[i].set(y, avgVal);
				}
			}
			
			System.out.println("Writing CSV");
			writeCSV(avgResults, mags, csvFile);
		}
		
		System.out.println("Plotting");
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<DiscretizedFunc> gainFuncs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(magsFunc.getMinX(), magsFunc.getMaxX());
		for (int i=0; i<avgResults.length; i++) {
			EvenlyDiscretizedFunc avgFunc = avgResults[i];
			
			float mag = (float)magsFunc.getX(i);
			avgFunc.setName("M≥"+mag);
			
			EvenlyDiscretizedFunc gainFunc = new EvenlyDiscretizedFunc(avgFunc.getMinX(), avgFunc.getMaxX(), avgFunc.size());
			for (int j=0; j<avgFunc.size(); j++)
				gainFunc.set(j, avgFunc.getY(j)/avgFunc.getY(0));
			gainFunc.setName(avgFunc.getName());
			
			funcs.add(avgFunc);
			gainFuncs.add(gainFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, cpt.getColor(mag)));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, regName+" Probabilities Over Time", "Year", "Probability");
		spec.setLegendVisible(true);
		
		int width = 1000;
		int height = 800;
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		// log plot
		gp.drawGraphPanel(spec, false, true, new Range(startYear, endYear), new Range(1e-3,1));
		gp.getChartPanel().setSize(width, height);
		gp.saveAsPNG(new File(outputDir, "probs_log.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "probs_log.pdf").getAbsolutePath());
		// now linear
		gp.drawGraphPanel(spec, false, false, new Range(startYear, endYear), new Range(0,1));
		gp.getChartPanel().setSize(width, height);
		gp.saveAsPNG(new File(outputDir, "probs_linear.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "probs_linear.pdf").getAbsolutePath());
		
		PlotSpec gainSpec = new PlotSpec(gainFuncs, chars, regName+" Probability Gains Over Time",
				"Year", "Probability Gain (ref="+startYear+")");
		gainSpec.setLegendVisible(true);
		gp.drawGraphPanel(gainSpec, false, false, new Range(startYear, endYear), null);
		gp.getChartPanel().setSize(width, height);
		gp.saveAsPNG(new File(outputDir, "gain.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "gain.pdf").getAbsolutePath());
	}

}
