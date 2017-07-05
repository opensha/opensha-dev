package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.GraphWindow;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.UCERF2_ComparisonSolutionFetcher;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

import com.google.common.collect.Lists;

public class RupLengthHistogramGenerator {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions");
		File solFile = new File(dir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		
		// fraction of earthquakes
		// rupture length
		
		double min = 25d;
		int num = 25;
		double delta = 50d;
		int[] counts = new int[num+1];
		double[] rates = new double[num+1];
		double[] moRates = new double[num+1];
		
		double halfDelta = 0.5*delta;
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
//		double scaleToTotal = 246;
		double scaleToTotal = 1;
		boolean doU2 = false;
		
		FaultSystemSolution u2sol = null;
		if (doU2)
			u2sol = UCERF2_ComparisonSolutionFetcher.getUCERF2Solution(FaultModels.FM2_1);
		
		HistogramFunction hist = new HistogramFunction(min, num, delta);
		HistogramFunction u2hist = new HistogramFunction(min, num, delta);
		HistogramFunction asDisrHist = new HistogramFunction(min, num, delta);
		HistogramFunction asDisrU2Hist = new HistogramFunction(min, num, delta);
		HistogramFunction momentHist = new HistogramFunction(min, num, delta);
		HistogramFunction u2momentHist = new HistogramFunction(min, num, delta);
		double maxLen = hist.getMaxX()+halfDelta;
		
		double totRate = StatUtils.sum(sol.getRateForAllRups());
		double rateScale = 1;
		if (scaleToTotal != 1)
			rateScale = scaleToTotal / totRate;
		
		double maxLengthFound = 0d;
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			double rate = sol.getRateForRup(r);
			double moment = MagUtils.magToMoment(rupSet.getMagForRup(r));
			double scaledRate = rate*rateScale;
			double length = rupSet.getLengthForRup(r)/1000d; // m to km
			
			if (length > maxLengthFound)
				maxLengthFound = length;
			
			int index;
			
			if (length < maxLen) {
				hist.add(length, scaledRate);
				momentHist.add(length, moment*rate);
				asDisrHist.add(length, 1d);
				index = hist.getClosestXIndex(length);
			} else {
				index = counts.length-1;
			}
			
			counts[index]++;
			rates[index] += rate;
			moRates[index] += moment*rate;
		}
//		asDisrHist.scale(1d/totRate);
		
		if (doU2) {
			FaultSystemRupSet u2rupSet = u2sol.getRupSet();
			totRate = StatUtils.sum(u2sol.getRateForAllRups());
			if (scaleToTotal != 1)
				rateScale = scaleToTotal / totRate;
			
			for (int r=0; r<u2rupSet.getNumRuptures(); r++) {
				double rate = u2sol.getRateForRup(r);
				double moment = MagUtils.magToMoment(u2rupSet.getMagForRup(r));
				double scaledRate = rate*rateScale;
				double length = u2rupSet.getLengthForRup(r)/1000d; // m to km
				
				if (length < maxLen) {
					u2hist.add(length, scaledRate);
					u2momentHist.add(length, moment*rate);
					asDisrU2Hist.add(length, 1d);
				}
			}
//			asDisrU2Hist.scale(1d/totRate);
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(hist);
		if (doU2)
			funcs.add(u2hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 5f, Color.RED));
		
		String title;
		if (scaleToTotal != 1)
			title = "Lengh Histogram (scaled to "+(int)scaleToTotal+" ruptures)";
		else
			title = "Length Histogram";
		GraphWindow gw = new GraphWindow(funcs, title, chars);
		gw.getGraphWidget().getGraphPanel().setBackgroundColor(Color.WHITE);
		gw.setX_AxisLabel("Length (km)");
		if (scaleToTotal != 1)
			gw.setY_AxisLabel("Num");
		else
			gw.setY_AxisLabel("Rate (1/yr)");
		
		funcs = Lists.newArrayList();
		funcs.add(asDisrHist);
		if (doU2)
			funcs.add(asDisrU2Hist);
		gw = new GraphWindow(funcs, "Length Histogram (as discretized)", chars);
		gw.getGraphWidget().getGraphPanel().setBackgroundColor(Color.WHITE);
		gw.setX_AxisLabel("Length (km)");
		gw.setY_AxisLabel("Num Ruptures");
		
		funcs = Lists.newArrayList();
		funcs.add(momentHist);
		if (doU2)
			funcs.add(u2momentHist);
		new GraphWindow(funcs, "Moment Rate Histogram", chars);
		
		funcs = Lists.newArrayList();
		funcs.add(momentHist.getCumulativeDistFunction());
		if (doU2)
			funcs.add(u2momentHist.getCumulativeDistFunction());
		chars.clear();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 5f, Color.RED));
		new GraphWindow(funcs, "Cumulative Moment Rates", chars);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Bin Start", "Bin End", "# Ruptures", "Sum Rate", "Sum Moment Rate");
		for (int i=0; i<counts.length; i++) {
			if (counts[i] == 0)
				break;
			List<String> line = Lists.newArrayList();
			
			double start, end;
			if (i == hist.size()) {
				start = hist.getMaxX()+halfDelta;
				end = maxLengthFound;
			} else {
				start = hist.getX(i)-halfDelta;
				end = hist.getX(i)+halfDelta;
			}
			
			line.add((float)start+"");
			line.add((float)end+"");
			line.add(counts[i]+"");
			line.add((float)rates[i]+"");
			line.add((float)moRates[i]+"");
			
			csv.addLine(line);
		}
		csv.writeToFile(new File("/tmp/length_counts.csv"));
	}

}
