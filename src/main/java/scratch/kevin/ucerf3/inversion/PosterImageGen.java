package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.collect.Lists;

public class PosterImageGen {
	
	private static final int png_thumb_width = 400;
	private static final int png_thumb_height = 400;
	
	private static final int png_med_width = 700;
	private static final int png_med_height = 700;
	
	private static final int width = 1000;
	private static final int height = 1000;
	
	private static final boolean tableOnly = false;
	private static final boolean highQuality = true;
	
	private static ArrayList<ArrayList<String>> wikiTable = new ArrayList<ArrayList<String>>();
	static {
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(0).add("!Dataset");
		wikiTable.get(0).add("!Northern California (Well Constrained)<br>39,075 elements");
		wikiTable.get(0).add("!Northern California (Poorly Constrained)<br>39,075 elements");
		wikiTable.get(0).add("!All California (Well Constrained)<br>198,260 elements");
		wikiTable.get(0).add("!All California (Poorly Constrained)<br>198,260 elements");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(1).add("!Energy vs Time");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(2).add("!Avg Energy vs Time");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(3).add("!Serial Time vs Parallel Time");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(4).add("!Time Speedup vs Time");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(5).add("!Std. Dev. vs Time");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(6).add("!Improvement vs Energy");
		wikiTable.add(new ArrayList<String>());
		wikiTable.get(7).add("!Time Speedup vs Threads");
	}
//	private static final String opensha_files_url = "http://opensha.usc.edu/ftp/kmilner/ucerf3/dsa_poster/";
	private static final String opensha_files_url = "http://opensha.usc.edu/ftp/kmilner/ucerf3/2011agu/";
	
	private static void saveImages(GraphWindow gw, File dir, String fName) throws IOException {
		GraphWidget gp = gw.getGraphWidget();
		gp.setBackgroundColor(Color.WHITE);
		gp.setSize(width, height);
		gp.drawGraph();
		gp.setVisible(true);
		
		gp.validate();
		gp.repaint();
		
		gp.saveAsPDF(new File(dir, fName+".pdf").getAbsolutePath());
		gp.saveAsPNG(new File(dir, fName+".png").getAbsolutePath());

		gw.setPlotLabelFontSize(20);
		gw.setAxisLabelFontSize(18);
		gw.setTickLabelFontSize(14);
		gp.setSize(png_med_width, png_med_height);
		gp.drawGraph();
		gp.setVisible(true);
		
		gp.validate();
		gp.repaint();
		gp.saveAsPNG(new File(dir, fName+".medium.png").getAbsolutePath());

		gw.setPlotLabelFontSize(16);
		gw.setAxisLabelFontSize(14);
		gw.setTickLabelFontSize(10);
		gp.setSize(png_thumb_width, png_thumb_height);
		gp.drawGraph();
		gp.setVisible(true);

		
		gp.validate();
		gp.repaint();
		gp.saveAsPNG(new File(dir, fName+".small.png").getAbsolutePath());
	}
	
	private static String getImageTableLine(String fName) {
		return "|["+opensha_files_url+fName+".png "+opensha_files_url+fName+".small.png]";
	}
	
	private static void handleDir(File dir, Range timeRange, Range energyRange,
			Range stdDevRange, Range improvementRange)
	throws IOException {
		// Range speedupRange, 
		
		int myAvgNumX = avgNumX;
		if (!highQuality)
			myAvgNumX = avgNumX/2;
		int myTargetNum = targetPPM;
		if (!highQuality)
			myTargetNum = targetPPM/2;
		
		HashMap<String, GraphWindow> windows = null;
		if (!tableOnly)
			windows = ResultPlotter.generatePlots(null, dir, highlight, coolType, threads, nodes,
				includeStartSubZero, plotAvg, bundleDsaBySubs, bundleTsaBySubs,
				myAvgNumX, myTargetNum, false, plots, null, null);
		GraphWindow gw;
		
//		wikiTable.get(0).add("!"+dir.getName());
		String prefix = dir.getName()+"_";
		String fName;
		
		fName = prefix+"e_vs_t";
		wikiTable.get(1).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.energy_vs_time_title);
			gw.setX_AxisRange(timeRange.getLowerBound(), timeRange.getUpperBound());
			gw.setY_AxisRange(energyRange.getLowerBound(), energyRange.getUpperBound());
			saveImages(gw, dir, fName);
		}
		
		fName = prefix+"avg_e_vs_t";
		wikiTable.get(2).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.avg_energy_vs_time_title);
			gw.setX_AxisRange(timeRange.getLowerBound(), timeRange.getUpperBound());
			gw.setY_AxisRange(energyRange.getLowerBound(), energyRange.getUpperBound());
			saveImages(gw, dir, fName);
		}

		fName = prefix+"st_vs_pt";
		wikiTable.get(3).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.time_comparison_title);
			gw.setX_AxisRange(timeRange.getLowerBound(), timeRange.getUpperBound());
//			gw.setY_AxisRange(speedupRange.getLowerBound(), speedupRange.getUpperBound());
			saveImages(gw, dir, fName);
		}

		fName = prefix+"spd_vs_t";
		wikiTable.get(4).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.time_speedup_vs_time_title);
			gw.setX_AxisRange(timeRange.getLowerBound(), timeRange.getUpperBound());
//			gw.setY_AxisRange(speedupRange.getLowerBound(), speedupRange.getUpperBound());
			saveImages(gw, dir, fName);
		}
		
		fName = prefix+"std_dev_vs_t";
		wikiTable.get(5).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.std_dev_vs_time_title);
			gw.setX_AxisRange(timeRange.getLowerBound(), timeRange.getUpperBound());
			gw.setY_AxisRange(stdDevRange.getLowerBound(), stdDevRange.getUpperBound());
			saveImages(gw, dir, fName);
		}
		
		fName = prefix+"imp_vs_t";
		wikiTable.get(6).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.improvement_vs_time_title);
			gw.setX_AxisRange(timeRange.getLowerBound(), timeRange.getUpperBound());
			gw.setY_AxisRange(improvementRange.getLowerBound(), improvementRange.getUpperBound());
			saveImages(gw, dir, fName);
		}
		
		fName = prefix+"spd_vs_thrd";
		wikiTable.get(7).add(getImageTableLine(fName));
		if (!tableOnly) {
			gw = windows.get(ResultPlotter.speedup_vs_threads_title);
			List<DiscretizedFunc> funcs = gw.getGraphWidget().getPlotSpec().getPlotFunctionsOnly();
			DiscretizedFunc spdFunc = funcs.get(0);
			spd_vs_thds.add(spdFunc);
			if (spd_vs_thds_comp.isEmpty()) {
				spd_vs_thds_comp.add(funcs.get(1));
				spd_vs_thds_comp.add(funcs.get(2));
			}
			saveImages(gw, dir, fName);
		}
	}
	
	private static ArrayList<DiscretizedFunc> spd_vs_thds = new ArrayList<DiscretizedFunc>();
	private static ArrayList<DiscretizedFunc> spd_vs_thds_comp = new ArrayList<DiscretizedFunc>();
	
	private static void handleSpdCompares(File dir) throws IOException {
		if (spd_vs_thds.size() != 4)
			return;
		
		spd_vs_thds.addAll(spd_vs_thds_comp);
		
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<PlotCurveCharacterstics>();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 8f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, PlotSymbol.FILLED_CIRCLE, 0f, Color.BLUE));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, PlotSymbol.FILLED_CIRCLE, 0f, Color.GREEN));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, PlotSymbol.FILLED_CIRCLE, 0f, Color.RED));
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GREEN));
		
		GraphWindow gw = ResultPlotter.getGraphWindow(spd_vs_thds,
				ResultPlotter.speedup_vs_threads_title, chars, ResultPlotter.threads_label,
				ResultPlotter.time_speedup_label, false);
		
		saveImages(gw, dir, "spd_vs_threads_comp");
		gw.setX_AxisRange(0d, 80d);
		gw.setY_AxisRange(0d, 80d);
		saveImages(gw, dir, "spd_vs_threads_evencomp");
	}
	
	private static String coolType = null;
	private static int threads = -1;
	private static int nodes = -1;
	private static boolean includeStartSubZero = false;
	private static boolean plotAvg = true;
	private static boolean bundleDsaBySubs = false;
	private static boolean bundleTsaBySubs = false;
	private static String highlight = null;
	
	private static int avgNumX = 600;
	private static int targetPPM = 4;
	
	private static ArrayList<String> plots = new ArrayList<String>();
	static {
		plots.add(ResultPlotter.energy_vs_time_title);
		plots.add(ResultPlotter.avg_energy_vs_time_title);
		plots.add(ResultPlotter.std_dev_vs_time_title);
		plots.add(ResultPlotter.improvement_vs_time_title);
		plots.add(ResultPlotter.time_comparison_title);
		plots.add(ResultPlotter.time_speedup_vs_time_title);
		plots.add(ResultPlotter.speedup_vs_threads_title);
	}
	
	private static void printTable() {
		System.out.println("***** TABLE *****\n");
		System.out.println("{| border=\"1\"");
		
		for (int i=0; i<wikiTable.size(); i++) {
			if (i > 0)
				System.out.println("|-");
			ArrayList<String> row = wikiTable.get(i);
			for (String cell : row)
				System.out.println(cell);
		}
		
		System.out.println("|}");
		
		System.out.println();
		System.out.println("Speedup Vs Threads Comparisons");
		System.out.println(getImageTableLine("spd_vs_threads_comp").substring(1));
		System.out.println(getImageTableLine("spd_vs_threads_evencomp").substring(1));
	}
	
	public static void main(String args[]) throws IOException {
//		File main = new File("/home/kevin/OpenSHA/UCERF3/test_inversion/bench/poster");
		File main = new File("/home/kevin/OpenSHA/UCERF3/test_inversion/bench/agu");
		
		File ncalConst = new File(main, "ncal_constrained");
		File ncalUnconst = new File(main, "ncal_unconstrained");
		File stateConst = new File(main, "allcal_constrained");
		File stateUnonst = new File(main, "allcal_unconstrained");
		
		handleDir(ncalConst,
				// time range
				new Range(0, 120),
				// energy plot range
				new Range(30, 70),
//				// speedup range
//				new Range(0.5, 21),
				// std. dev. range
				new Range(0, 10),
				// % improvement range
				new Range(0, 10));
		
		handleDir(ncalUnconst,
				// time range
				new Range(0, 120),
				// energy plot range
				new Range(0, 20),
//				// speedup range
//				new Range(0.5, 10),
				// std. dev. range
				new Range(0, 10),
				// % improvement range
				new Range(0, 10));
		
		handleDir(stateConst,
				// time range
				new Range(0, 480),
				// energy plot range
				new Range(25, 400),
//				// speedup range
//				new Range(0.5, 6),
				// std. dev. range
				new Range(0, 10),
				// % improvement range
				new Range(0, 10));
		
		handleDir(stateUnonst,
				// time range
				new Range(0, 480),
				// energy plot range
				new Range(0, 300),
//				// speedup range
//				new Range(0.5, 6.5),
				// std. dev. range
				new Range(0, 10),
				// % improvement range
				new Range(0, 10));
		
		handleSpdCompares(main);
		writeSubIGraphs(main);
		
		printTable();
		
		System.exit(0);
	}
	
	private static void writeSubIGraphs(File dir) throws IOException {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		
		func.set(100d, 18.869238);
		func.set(500d, 7.675522);
		func.set(50d, 18.955439);
		func.set(200d, 13.790538);
		func.set(25d, 18.423607);
		
		ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<PlotCurveCharacterstics>();
		chars.add(new PlotCurveCharacterstics(
				PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 8f, Color.BLACK));
		
		GraphWindow gw = ResultPlotter.getGraphWindow(funcs,
				"Speedup vs nSubIterations (40 Threads)", chars, "nSubIterations",
				"Speedup", false);
		
		saveImages(gw, dir, "spd_vs_subs_comp");
	}

}
