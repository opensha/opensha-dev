package scratch.kevin.cybershake.etasCalcs;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class DurationLogPlotter {
	
	private static final long day_millis = 24l*60l*60l*1000l;
	
	private static final SimpleDateFormat df = MPJTaskCalculator.df;
	private static final DecimalFormat meanFormat = new DecimalFormat("0.00");
	
	private static PlotSpec buildPlot(File logFile) throws IOException {
		Map<Integer, Date> startDates = Maps.newHashMap();
		Map<Integer, Date> endDates = Maps.newHashMap();
		
		BufferedReader tis =
				new BufferedReader(new InputStreamReader(new FileInputStream(logFile)));
		String str = tis.readLine();
		while(str != null) {
			if (str.startsWith("[") && (str.contains("calculating") || str.contains("completed"))) {
				if (str.contains("batch")) {
					str = tis.readLine();
					continue;
				}
				// isolate date
				String dateStr = str.substring(1);
				dateStr = dateStr.substring(0, dateStr.indexOf(" "));
				Date date;
				try {
					date = df.parse(dateStr);
				} catch (ParseException e) {
					str = tis.readLine();
					continue;
				}
				
				// now get index
				String[] split = str.split(" ");
				try {
					Integer index = Integer.parseInt(split[split.length-1]);
					if (str.contains("calculating")) {
//						Preconditions.checkState(!startDates.containsKey(index), "Duplicate for "+index);
						startDates.put(index, date);
					} else {
//						Preconditions.checkState(!endDates.containsKey(index));
						endDates.put(index, date);
					}
				} catch (NumberFormatException e) {
					str = tis.readLine();
					continue;
				}
			}
			str = tis.readLine();
		}
		tis.close();
		
		System.out.println("Successfully parsed "+startDates.size()+" start dates and "+endDates.size()+" end dates");
		
		Map<Integer, Double> deltaSecsMap = Maps.newHashMap();
		
		for (Integer index : endDates.keySet()) {
			if (!startDates.containsKey(index))
				continue;
			Date start = startDates.get(index);
			Date end = endDates.get(index);
			
			long diffMillis = end.getTime() - start.getTime();
			while (diffMillis < 0)
				diffMillis += day_millis;
			double diffSecs = diffMillis / 1000d;
			
			deltaSecsMap.put(index, diffSecs);
		}
		
//		HistogramFunction hist = new HistogramFunction(0.5, 100, 1d);
		HistogramFunction hist = new HistogramFunction(0.5, 1000, 1d);
		
		double[] vals = new double[deltaSecsMap.size()];
		int cnt = 0;
		
		for (Integer index : deltaSecsMap.keySet()) {
			double delta = deltaSecsMap.get(index)/60d; // convert to minutes
			int histInd = hist.getClosestXIndex(delta);
			hist.add(histInd, 1d);
			
			vals[cnt++] = delta;
		}
		
		double mean = StatUtils.mean(vals);
		System.out.println("Mean: "+mean+" m, Median: "+DataUtils.median(vals)+" m");
		
		List<XY_DataSet> elems = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		hist.setName("Histogram");
		elems.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
//		DefaultXY_DataSet xy = new DefaultXY_DataSet();
//		xy.set(mean, 0d);
//		xy.set(mean, hist.getMaxY()*1.25);
//		xy.setName("Mean: "+(float)mean);
//		
//		elems.add(xy);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		PlotSpec histSpec = new PlotSpec(elems, chars, "Simulation Time Hist", "Duration (minutes)", "Number");
		
		List<XYTextAnnotation> anns = Lists.newArrayList();
//		XYTextAnnotation ann = new XYTextAnnotation("Mean: "+meanFormat.format(mean)+" m", 10, hist.getMaxY()*1.2);
//		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		XYTextAnnotation ann = new XYTextAnnotation("Mean: "+meanFormat.format(mean)+" m", 10, 1000);
		ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
		ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
		anns.add(ann);
		histSpec.setPlotAnnotations(anns);
		
		return histSpec;
	}

	public static void main(String[] args) throws IOException {
		
		File[] logFiles = {
//						new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2015_05_13-mojave_7/"
//								+ "2015_05_13-mojave_7.pbs.o5290046"),
//						new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2015_06_10-mojave_7/"
//								+ "2015_06_10-mojave_7.pbs.o5388265")
//						new File("/tmp/2015_11_04-spontaneous-1000yr-full_td-maxChar10.0-noApplyLTR.pbs.o5997982")
		};
		
//		File logFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2015_05_13-mojave_7/"
//				+ "2015_05_13-mojave_7.pbs.o5290046");
//		File logFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2015_06_10-mojave_7/"
//				+ "2015_06_10-mojave_7.pbs.o5388265");
//		File logFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2015_10_09-mojave_m7-full_td-maxChar10.0/"
//				+ "2015_10_09-mojave_m7-full_td-maxChar10.0.pbs.o5880378");
		
		Preconditions.checkState(logFiles.length > 0);
		
		List<PlotSpec> specs = Lists.newArrayList();
		
		List<Range> xRanges = Lists.newArrayList();
		
		Range yRange = new Range(0, 0);
		
		for (File logFile : logFiles) {
			PlotSpec spec = buildPlot(logFile);
			specs.add(spec);
			
			yRange = new Range(0d, Math.max(yRange.getUpperBound(), getMaxRange(spec, false).getUpperBound()*1.25));
			Range xRange = new Range(0d, getMaxRange(spec, true).getUpperBound());
			xRanges.add(xRange);
		}
		List<Range> yRanges = Lists.newArrayList(yRange);
		
		if (specs.size() == 1) {
			GraphWindow gw = new GraphWindow(specs.get(0));
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		} else {
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setBackgroundColor(Color.WHITE);
			
			gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
			
			gp.getChartPanel().setSize(200+600*specs.size(), 400);
			gp.saveAsPNG(new File("/tmp/plot_duration.png").getAbsolutePath());
		}
	}
	
	private static Range getMaxRange(PlotSpec spec, boolean x) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		
		for (PlotElement elem : spec.getPlotElems()) {
			if (elem instanceof XY_DataSet) {
				XY_DataSet xy = (XY_DataSet)elem;
				if (x) {
					track.addValue(xy.getMinX());
					track.addValue(xy.getMaxX());
				} else {
					track.addValue(xy.getMinY());
					track.addValue(xy.getMaxY());
				}
			}
		}
		
		return new Range(track.getMin(), track.getMax());
	}

}
