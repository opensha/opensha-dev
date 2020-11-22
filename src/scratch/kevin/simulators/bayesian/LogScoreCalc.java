package scratch.kevin.simulators.bayesian;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import com.google.common.base.Preconditions;

public class LogScoreCalc {
	
	private int numFaults;
	
	private FaultStateEventCalc trainingData;
	private FaultStateEventCalc testData;

	public LogScoreCalc(FaultStateEventCalc trainingData, FaultStateEventCalc testData) {
		Preconditions.checkState(trainingData.getNumFaults() == testData.getNumFaults());
		this.numFaults = trainingData.getNumFaults();
		this.trainingData = trainingData;
		this.testData = testData;
	}
	
	public double calcLogScore(boolean[] state, double alpha) {
		int numTest = testData.getCount(state);
		double forecastRate = trainingData.getLaplaceProb(state, alpha);
		return -numTest*Math.log(forecastRate);
	}
	
	public double calcTotalLogScore(double alpha) {
		boolean[] state = new boolean[numFaults];
		return calcTotalLogScore(alpha, state, 0);
	}
	
	private double calcTotalLogScore(double alpha, boolean[] state, int pos) {
		double sum = 0d;
		if (pos == state.length-1) {
			state[pos] = true;
			sum += calcLogScore(state, alpha);
			state[pos] = false;
			sum += calcLogScore(state, alpha);
		} else {
			state[pos] = true;
			sum += calcTotalLogScore(alpha, state, pos+1);
			state[pos] = false;
			sum += calcTotalLogScore(alpha, state, pos+1);
		}
		return sum;
	}
	
	public DiscretizedFunc calcScoresForAlpha(double minAlpha, double maxAlpha, int num, boolean logSpacing) {
		Preconditions.checkArgument(maxAlpha > minAlpha);
		if (logSpacing) {
			Preconditions.checkArgument(minAlpha > 0);
			EvenlyDiscretizedFunc logFunc = new EvenlyDiscretizedFunc(Math.log10(minAlpha), Math.log10(maxAlpha), num);
			for (int i=0; i<num; i++)
				logFunc.set(i, calcTotalLogScore(Math.pow(10, logFunc.getX(i))));
			ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : logFunc)
				func.set(Math.pow(10, pt.getX()), pt.getY());
			return func;
		} else {
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(minAlpha, maxAlpha, num);
			for (int i=0; i<num; i++)
				func.set(i, calcTotalLogScore(func.getX(i)));
			return func;
		}
	}
	
	public double minimizeAlpha(double minAlpha, double maxAlpha, int num, boolean logSpacing) {
		DiscretizedFunc alphaFunc = calcScoresForAlpha(minAlpha, maxAlpha, num, logSpacing);
		return minimizeAlpha(alphaFunc);
	}
	
	public double minimizeAlpha(DiscretizedFunc alphaFunc) {
		return minimizeAlphaRecursive(alphaFunc, 5, 10);
	}
	
	private double minimizeAlphaRecursive(DiscretizedFunc alphaFunc, int numLeft, int numEach) {
		double min = alphaFunc.getMinY();
		int closestX = -1;
		for (int i=0; i<alphaFunc.size(); i++) {
			if (alphaFunc.getY(i) == min) {
				closestX = i;
				break;
			}
		}
		Preconditions.checkState(closestX >= 0 && closestX < alphaFunc.size());
		if (closestX == 0 || closestX == alphaFunc.size()-1 || numLeft == 0)
			return alphaFunc.getX(closestX);
		// now narrow in on it
		double xBefore = alphaFunc.getX(closestX-1);
		double xAfter = alphaFunc.getX(closestX+1);
		DiscretizedFunc miniAlphaFunc = calcScoresForAlpha(xBefore, xAfter, 1000, xBefore > 0);
		return minimizeAlphaRecursive(miniAlphaFunc, numLeft-1, numEach);
	}
	
	public void plotScoreVsAlpha(File outputDir, String prefix, String title, double minAlpha, double maxAlpha, int num,
			boolean logX, boolean logY) throws IOException {
		plotScoreVsAlpha(outputDir, prefix, title, calcScoresForAlpha(minAlpha, maxAlpha, num, logX), logX, logY);
	}
	
	public void plotScoreVsAlpha(File outputDir, String prefix, String title, DiscretizedFunc alphaFunc,
			boolean logX, boolean logY) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = new Range(alphaFunc.getMinX(), alphaFunc.getMaxX());
		Range yRange;
		double funcMinY = alphaFunc.getMinY();
		double funcMaxY = alphaFunc.getMaxY();
		if (logY) {
			funcMinY = Math.log10(funcMinY);
			funcMaxY = Math.log10(funcMaxY);
		}
		double funcWidthY = funcMaxY-funcMinY;
		yRange = new Range(funcMinY - 0.05*funcWidthY, funcMaxY + 0.05*funcWidthY);
		System.out.println(yRange);
		if (logY) {
			yRange = new Range(Math.pow(10, yRange.getLowerBound()), Math.pow(10, yRange.getUpperBound()));
			System.out.println(yRange);
		}
		
		double minScore = alphaFunc.getMinY();
		DiscretizedFunc minFunc = new ArbitrarilyDiscretizedFunc();
		minFunc.set(xRange.getLowerBound(), minScore);
		minFunc.set(xRange.getUpperBound(), minScore);
		
		funcs.add(minFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, Color.GRAY));
		
		double bestAlpha = minimizeAlpha(alphaFunc);
		System.out.println("Best Alpha: "+bestAlpha);
		DefaultXY_DataSet bestAlphaFunc = new DefaultXY_DataSet();
		bestAlphaFunc.set(bestAlpha, yRange.getLowerBound());
		bestAlphaFunc.set(bestAlpha, yRange.getUpperBound());
		
		funcs.add(bestAlphaFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, Color.GRAY));
		
		funcs.add(alphaFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		List<XYAnnotation> anns = new ArrayList<>();
		
		String alphaLabel;
		if (bestAlpha > 100)
			alphaLabel = intDF.format(bestAlpha);
		else if (bestAlpha > 10)
			alphaLabel = oneDigitDF.format(bestAlpha);
		else if (bestAlpha > 0.1)
			alphaLabel = twoDigitsDF.format(bestAlpha);
		else
			alphaLabel = scientificDF.format(bestAlpha);
		XYTextAnnotation alphaAnn = new XYTextAnnotation("Best Alpha: "+alphaLabel, bestAlpha, minScore);
		alphaAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
		alphaAnn.setTextAnchor(TextAnchor.TOP_LEFT);
		anns.add(alphaAnn);
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Alpha", "Log Score");
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, logX, logY, xRange, yRange);
		gp.getChartPanel().setSize(800, 700);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	static final DecimalFormat intDF = new DecimalFormat("0");
	static final DecimalFormat oneDigitDF = new DecimalFormat("0.0");
	static final DecimalFormat twoDigitsDF = new DecimalFormat("0.00");
	static final DecimalFormat scientificDF = new DecimalFormat("0E00");

	public static void main(String[] args) throws IOException {
		File rsDir = new File("/data/kevin/simulators/catalogs/rundir2585_1myr");
		String prefix = "catalog_event_probs_9faults_m7_10yr";
		FaultStateEventCalc trainingData = FaultStateEventCalc.loadStatesCSV(new File(rsDir, prefix+"_first_half.csv"));
		FaultStateEventCalc testData = FaultStateEventCalc.loadStatesCSV(new File(rsDir, prefix+"_second_half.csv"));
		
		LogScoreCalc scoreCalc = new LogScoreCalc(trainingData, testData);
		DiscretizedFunc alphaFunc = scoreCalc.calcScoresForAlpha(1e-2, 1e10, 1000, true);
		scoreCalc.plotScoreVsAlpha(rsDir, prefix+"_alpha_vs_score", "RSQSim 1st Half to Predict 2nd Half", alphaFunc, true, false);
		double bestAlpha = scoreCalc.minimizeAlpha(alphaFunc);
		DiscretizedFunc alphaFuncTight = scoreCalc.calcScoresForAlpha(bestAlpha/3d, bestAlpha*3d, 1000, false);
		scoreCalc.plotScoreVsAlpha(rsDir, prefix+"_alpha_vs_score_tight", "RSQSim 1st Half to Predict 2nd Half", alphaFuncTight, false, false);
		
		// now compare with U3
		// use whole RSQSim catalog
		trainingData = FaultStateEventCalc.loadStatesCSV(new File(rsDir, prefix+"_first_half.csv"));
		File u3Dir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_catalogs/2017_04_10-long-catalogs/batch0/100000yr_run0");
		String u3Prefix = prefix.replaceAll("catalog_", "u3_");
		testData = FaultStateEventCalc.loadStatesCSV(new File(u3Dir, u3Prefix+".csv"));
		
		String u3OutPrefix = prefix.replaceAll("catalog_", "catalog_vs_u3_");
		scoreCalc = new LogScoreCalc(trainingData, testData);
		alphaFunc = scoreCalc.calcScoresForAlpha(1e-2, 1e10, 1000, true);
		scoreCalc.plotScoreVsAlpha(rsDir, u3OutPrefix+"_alpha_vs_score", "RSQSim Entire Catalog to Predict UCERF3", alphaFunc, true, false);
		bestAlpha = scoreCalc.minimizeAlpha(alphaFunc);
		alphaFuncTight = scoreCalc.calcScoresForAlpha(bestAlpha/3d, bestAlpha*3d, 1000, false);
		scoreCalc.plotScoreVsAlpha(rsDir, u3OutPrefix+"_alpha_vs_score_tight", "RSQSim Entire Catalog to Predict UCERF3", alphaFuncTight, false, false);
	}

}
