package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.SynchIdens.SynchFaults;
import scratch.kevin.simulators.momRateVariation.SimulatorMomRateVarCalc;
import scratch.kevin.simulators.momRateVariation.UCERF3ComparisonAnalysis;
import scratch.kevin.simulators.momRateVariation.UCERF3_ETASComparisons;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class RecurrencePlotGen {
	
	public enum DistanceMetric {
		L1_NORM {
			@Override
			public double calc(double[] state1, double[] state2) {
				Preconditions.checkArgument(state1.length == state2.length);
				double sum = 0d;
				for (int i=0; i<state1.length; i++)
					sum += Math.abs(state1[i] - state2[i]);
				return sum;
			}
		},
		L2_NORM {
			@Override
			public double calc(double[] state1, double[] state2) {
				Preconditions.checkArgument(state1.length == state2.length);
				double sum = 0d;
				for (int i=0; i<state1.length; i++)
					sum += Math.pow(state1[i] - state2[i], 2);
				return sum;
			}
		},
		LINFINITY_NORM {
			@Override
			public double calc(double[] state1, double[] state2) {
				Preconditions.checkArgument(state1.length == state2.length);
				double ret = 0d;
				for (int i=0; i<state1.length; i++) {
					double diff = Math.abs(state1[i] - state2[i]);
					if (diff > ret)
						ret = diff;
				}
				return ret;
			}
		};
		
		public abstract double calc(double[] state1, double[] state2);
	}
	
	public static BitSet[] calcBitSet(List<double[]> fullPath, DistanceMetric distCalc, double threshold) {
		System.out.println("Calculating for "+distCalc.name()+", thresh="+threshold);
		BitSet[] ret = new BitSet[fullPath.size()];
		for (int i=0; i<fullPath.size(); i++)
			ret[i] = new BitSet(fullPath.size());
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		int numBelow = 0;
		int tot = 0;
		
		boolean below;
		for (int i=0; i<fullPath.size(); i++) {
			for (int j=0; j<fullPath.size(); j++) {
				double val = distCalc.calc(fullPath.get(i), fullPath.get(j));
				track.addValue(val);
				below = val <= threshold;
				if (below) {
					numBelow++;
					ret[i].set(j);
				}
				tot++;
			}
		}
		
		double percent = 100d*numBelow/(double)tot;
		System.out.println(numBelow+"/"+tot+" below threshold of "+threshold+" ("+(float)percent+" %)");
		System.out.println(track);
		
		return ret;
	}
	
	public static double[][] calcDist(List<double[]> fullPath, DistanceMetric distCalc) {
		System.out.println("Calculating for "+distCalc.name());
		double[][] data = new double[fullPath.size()][fullPath.size()];
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		
		for (int i=0; i<fullPath.size(); i++) {
			for (int j=0; j<fullPath.size(); j++) {
				double val = distCalc.calc(fullPath.get(i), fullPath.get(j));
				track.addValue(val);
				data[i][j] = val;
			}
		}
		
		System.out.println(track);
		
		return data;
	}
	
	public static double[][] calcRotated(List<double[]> fullPath, DistanceMetric distCalc, int width) {
		Preconditions.checkArgument(width % 2 == 1, "Width must be odd so that there's a point on the axis");
		
		long len = fullPath.size()*width;
		System.out.println("Calculating rotated for "+distCalc.name()+", size: "+fullPath.size()+" x "+width+" = "+len);
		
		int numBefore = width / 2;
		
		double[][] data = new double[fullPath.size()][width];
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		
		for (int i=0; i<fullPath.size(); i++) {
			for (int n=0; n<width; n++) {
				int j = i - numBefore + n;
				double val;
				if (j < 0 || j >= fullPath.size()) {
					val = Double.NaN;
				} else {
					val = distCalc.calc(fullPath.get(i), fullPath.get(j));
					track.addValue(val);
				}
				data[i][n] = val;
			}
		}
		
		System.out.println(track);
		
		return data;
	}
	
	private static BitSet[] distToBitSet(double[][] data, double threshold) {
		BitSet[] ret = new BitSet[data.length];
		for (int i=0; i<data.length; i++)
			ret[i] = new BitSet(data.length);
		
		for (int i=0; i<data.length; i++)
			for (int j=0; j<data.length; j++)
				if (data[i][j] <= threshold)
					ret[i].set(j);
		
		return ret;
	}
	
	private static int width = 650;
	private static int height = 700;
	
//	private static CPT discreteCPT;
//	static {
//		discreteCPT = new CPT();
//		discreteCPT.setBelowMinColor(Color.WHITE);
//		discreteCPT.add(new CPTVal(0f, Color.WHITE, 0.5f, Color.WHITE));
//		discreteCPT.add(new CPTVal(0.5f, Color.BLACK, 1f, Color.BLACK));
//		discreteCPT.setAboveMaxColor(Color.BLACK);
//		discreteCPT.setNanColor(Color.GRAY);
//	}
	
	private static PlotPreferences prefs = XYZGraphPanel.getDefaultPrefs();
	static {
		prefs.setPlotLabelFontSize(50);
		prefs.setAxisLabelFontSize(40);
		prefs.setTickLabelFontSize(25);
	}
	
	private static Color trans_white = new Color(255, 255, 255, 170);
	private static Font annotation_font = new Font(Font.SANS_SERIF, Font.BOLD, 60);
	
	private static CPT getDiscreteCPT(double threshold) {
		double maxZ = threshold * 5d;
		CPT discreteCPT = new CPT();
		discreteCPT.setBelowMinColor(Color.BLACK);
		discreteCPT.add(new CPTVal(0f, Color.BLACK, (float)threshold, Color.BLACK));
		discreteCPT.add(new CPTVal((float)threshold, Color.WHITE, (float)maxZ, Color.WHITE));
		discreteCPT.setAboveMaxColor(Color.WHITE);
		discreteCPT.setNanColor(Color.GRAY);
		return discreteCPT;
	}
	
	public static void plotDiscrete(BitSet[] data, DistanceMetric distCalc, double threshold, File outputFile, int zoom,
			double distSpacing) 	throws IOException {
		
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(
				data.length, data.length, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
		for (int i=0; i<data.length; i++)
			for (int j=0; j<data.length; j++)
				if (data[i].get(j))
					xyz.set(i, j, 0d);
				else
					xyz.set(i, j, 1d);
		
		String title = "Recurrence Plot, "+distCalc.name()+", thresh="+(float)threshold;
		plotSquare(xyz, title, getDiscreteCPT(threshold), outputFile, zoom);
	}
	
	public static void plotContinuous(double[][] data, DistanceMetric distCalc, double maxZ, File outputFile, int zoom,
			double distSpacing) throws IOException {
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(data, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().reverse();
		if (maxZ > 0)
			cpt = cpt.rescale(0d, maxZ);
		else
			cpt = cpt.rescale(0d, xyz.getMaxZ());
		String title = "Recurrence Plot, "+distCalc.name();
		
		plotSquare(xyz, title, cpt, outputFile, zoom);
	}
	
	public static void plotHybrid(double[][] data, DistanceMetric distCalc, double threshold, double maxZ,
			File outputFile, int zoom, double distSpacing) throws IOException {
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(data, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
		CPT cpt = getHybridCPT(threshold, maxZ);
		String title = "Recurrence Plot, "+distCalc.name()+", thresh="+(float)threshold;
		
		plotSquare(xyz, title, cpt, outputFile, zoom);
	}
	
	private static void plotSquare(EvenlyDiscrXYZ_DataSet xyz, String title, CPT cpt, File outputFile, int zoom)
			throws IOException {
		XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, title, "Years", "Years", "");
		
		if (zoom <= 0)
			zoom = xyz.getNumX();
		double max = (zoom+0.5)*xyz.getGridSpacingX();
		XYZGraphPanel panel = new XYZGraphPanel();
		panel.drawPlot(spec, false, false, new Range(0, max), new Range(0, max));
		
		if (outputFile == null) {
			// display it
			XYZPlotWindow window = new XYZPlotWindow(panel);
			window.setSize(width, height);
			window.setDefaultCloseOperation(XYZPlotWindow.EXIT_ON_CLOSE);
		} else {
			// write plot
			panel.getChartPanel().setSize(width, height);
			panel.saveAsPNG(outputFile.getAbsolutePath());
		}
	}
	
	private static CPT getHybridCPT(double threshold, double maxZ) throws IOException {
		Preconditions.checkState(threshold > 0);
		Preconditions.checkState(maxZ > threshold);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().reverse();
		cpt = cpt.rescale(0d, maxZ);
		cpt.setNanColor(Color.GRAY);
		
		// saturate
		for (CPTVal v : cpt) {
			v.minColor = saturate(v.minColor);
			v.maxColor = saturate(v.maxColor);
		}
		cpt.setAboveMaxColor(cpt.getMaxColor());
		
		// now threshold
		Color cAtThresh = cpt.getColor((float)threshold);
		for (int i=cpt.size(); --i >= 0;)
			if (cpt.get(i).start <= threshold)
				cpt.remove(i);
		Preconditions.checkState(!cpt.isEmpty(), "Threshold and maxZ too close!");
		cpt.add(0, new CPTVal((float)threshold, cAtThresh, cpt.get(0).start, cpt.get(0).minColor));
		cpt.add(0, new CPTVal(0f, Color.BLACK, (float)threshold, Color.BLACK));
		
		return cpt;
	}
	
	private static Color saturate(Color c) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		
		int saturationSteps = 2;
		
		for (int i=0; i<saturationSteps; i++) {
			r = (int)(0.5d*(r + 255d)+0.5);
			g = (int)(0.5d*(g + 255d)+0.5);
			b = (int)(0.5d*(b + 255d)+0.5);
		}
		
		return new Color(r, g, b);
	}
	
	private static boolean plot_rotated_preserve_min = true;
	
	public static void plotContinuousRotated(double[][] data, DistanceMetric distCalc, int widthEach,
			double maxZ, double distSpacing, File outputFile)
			throws IOException {
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().reverse();
		if (maxZ <= 0) {
			maxZ = 0;
			for (double[] vals : data)
				maxZ = Math.max(maxZ, StatUtils.max(vals));
		}
		cpt = cpt.rescale(0d, maxZ);
		cpt.setNanColor(Color.GRAY);
		
		plotRotated(data, distCalc, cpt, widthEach, distSpacing, outputFile);
	}
	
	public static void plotDiscreteRotated(double[][] data, DistanceMetric distCalc, int widthEach,
			double threshold, double distSpacing, File outputFile) throws IOException {
		plotRotated(data, distCalc, getDiscreteCPT(threshold), widthEach, distSpacing, outputFile);
	}
	
	public static void plotHybridRotated(double[][] data, DistanceMetric distCalc, int widthEach,
			double threshold, double maxZ, double distSpacing, File outputFile) throws IOException {
		plotRotated(data, distCalc, getHybridCPT(threshold, maxZ), widthEach, distSpacing, outputFile);
	}
	
	private static void plotRotated(double[][] data, DistanceMetric distCalc, CPT cpt, int widthEach,
			double distSpacing, File outputFile) throws IOException {
		List<double[][]> datas = Lists.newArrayList();
		datas.add(data);
		plotRotated(datas, null, distCalc, cpt, widthEach, distSpacing, outputFile);
	}
	
	private static XYZGraphPanel plotRotated(List<double[][]> datas, List<String> labels, DistanceMetric distCalc,
			CPT cpt, int widthEach, double distSpacing, File outputFile) throws IOException {
		return plotRotated(datas, labels, distCalc, cpt, widthEach, distSpacing, outputFile, null, null);
	}
	
	private static final Color[] rup_colors = { Color.RED.darker(), Color.BLUE.darker(),
			Color.GREEN.darker(), Color.BLACK, Color.ORANGE.darker() };
	
	private static XYZGraphPanel plotRotated(List<double[][]> datas, List<String> labels, DistanceMetric distCalc,
			CPT cpt, int widthEach, double distSpacing, File outputFile,
			List<List<int[]>> paths, List<String> idenNames) throws IOException {
		// now rotate to pixel space
		int shiftDist = widthEach / 2;
		int overlapEach = 1; // 0.5 on each end\
		
		Preconditions.checkState(widthEach % 2 == 1);
		Preconditions.checkState(widthEach >= 3);
		
		double cptTickUnit;
		if (cpt.getMaxValue() <= 1f)
			cptTickUnit = 0.1;
		else if (cpt.getMaxValue() <= 3f)
			cptTickUnit = 0.25;
		else
			cptTickUnit = (float)(Math.round(100d*cpt.getMaxValue()/100d)/10d);
		
		List<Range> yRanges = Lists.newArrayList();
		Range xRange = null;
		int plotWidth = 0;
		List<XYZPlotSpec> specs = Lists.newArrayList();
		
		for (int d=0; d<datas.size(); d++) {
			double[][] data = datas.get(d);
			int jWidth = data[0].length;
			
			int totalWidth = (data.length + data[0].length/2)*(widthEach - overlapEach) + widthEach;
			int totalHeight = (widthEach - shiftDist) + shiftDist*jWidth;
			
			double gridSpacingX = distSpacing*0.5d/shiftDist;
			double gridSpacingY = distSpacing/shiftDist;
			double deltaY = distSpacing*(jWidth+1d); // +1 to cover the height of each pixel, half on top and half on bottom
			double minY = -0.5*deltaY;
			
			double xOffset = distSpacing*0.25*(jWidth-1d);
//			double minX = -xOffset - 0.5;
			double minX = -xOffset;
			
			short[][] counts = new short[totalWidth][totalHeight];
			EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(totalWidth, totalHeight, minX, minY, gridSpacingX, gridSpacingY);
			// initialize to dummy val
			if (plot_rotated_preserve_min)
				for (int x=0; x<totalWidth; x++)
					for (int y=0; y<totalHeight; y++)
						xyz.set(x, y, Double.POSITIVE_INFINITY);
			
			for (int i=0; i<data.length; i++) {
				for (int j=0; j<data[i].length; j++) {
					// find center of this rotated pixel
					
					int negJ = (jWidth-1)-j;
					int y = shiftDist + shiftDist*negJ;
					
					int x = shiftDist + (widthEach-overlapEach)*i + shiftDist*j;
					
					double val = data[i][j];
					
					for (int yDist=0; yDist<=widthEach; yDist++) {
						int[] yAdds;
						if (yDist == 0)
							yAdds = new int[] {0};
						else
							yAdds = new int[] { yDist, -yDist};
						
						for (int yAdd : yAdds) {
							int xDelta = shiftDist - yDist;
							for (int xAdd=-xDelta; xAdd<=xDelta; xAdd++) {
								setRotated(xyz, counts, x+xAdd, y+yAdd, val);
							}
						}
					}
				}
			}
			
			// now post process
			for (int x=0; x<totalWidth; x++) {
				for (int y=0; y<totalHeight; y++) {
					if (counts[x][y] == 0)
						xyz.set(x, y, Double.NaN);
					else {
						if (counts[x][y] > 1 && !plot_rotated_preserve_min)
							xyz.set(x, y, xyz.get(x, y)/(double)counts[x][y]);
//						if (threshold > 0) {
//							// discrete plot
//							double newVal;
//							if (xyz.get(x, y) <= threshold)
//								newVal = 1d;
//							else
//								newVal = 0d;
//							xyz.set(x, y, newVal);
//						}
					}
				}
			}
			
			String title = "Recurrence Plot, "+distCalc.name();
			XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, title, "Years", "Δ Years", "");
			
			yRanges.add(new Range(xyz.getMinY(), xyz.getMaxY()));
			int myPlotWidth = 100 + (int)(1.75d*(height-140d)*data.length/data[0].length);
			if (myPlotWidth > plotWidth)
				plotWidth = myPlotWidth;
			Range myXRange = new Range(-0.5d, xyz.getMaxX()-xOffset);
			if (xRange == null)
				xRange = myXRange;
			else
				xRange = new Range(Math.min(xRange.getLowerBound(), myXRange.getLowerBound()),
						Math.max(xRange.getUpperBound(), myXRange.getUpperBound()));
			
			spec.setCPTTickUnit(cptTickUnit);
			
			specs.add(spec);
			
			if (paths != null && paths.size() > d) {
				// add actual ruptures
				List<int[]> path = paths.get(d);
				Preconditions.checkState(idenNames.size() == path.get(0).length);
				List<ArbitrarilyDiscretizedFunc> rupFuncs = Lists.newArrayList();
				List<PlotCurveCharacterstics> chars = Lists.newArrayList();
				for (int i=0; i<idenNames.size(); i++) {
					rupFuncs.add(new ArbitrarilyDiscretizedFunc());
					chars.add(new PlotCurveCharacterstics(
							PlotSymbol.FILLED_SQUARE, 5f, rup_colors[i % rup_colors.length]));
				}
				
				double ySpacing = distSpacing*2.2;
				double startY = xyz.getMinY() - ySpacing;
				double endY = startY - ySpacing*(idenNames.size()+0.5);
				
				for (int i=0; i<path.size(); i++) {
					int[] state = path.get(i);
					double x = i*distSpacing + 0.5*distSpacing;
					for (int j=0; j<state.length; j++) {
						if (state[j] == 0) {
							double y = startY - ySpacing*j;
							rupFuncs.get(j).set(x, y);
						}
					}
				}
				yRanges.set(yRanges.size()-1, new Range(endY, xyz.getMaxY()));
				spec.setXYElems(rupFuncs);
				spec.setXYChars(chars);
			}
		}
		
		if (labels != null) {
			for (int i=0; i<specs.size(); i++) {
				String label = labels.get(i);
				
				double x = xRange.getLowerBound() + 0.05*(xRange.getUpperBound()-xRange.getLowerBound());
				Range yRange = yRanges.get(i);
				double y = yRange.getLowerBound() + 0.9*(yRange.getUpperBound()-yRange.getLowerBound());
				
				XYTextAnnotation ann = new XYTextAnnotation(label, x, y);
				ann.setFont(annotation_font);
				ann.setTextAnchor(TextAnchor.TOP_LEFT);
				ann.setPaint(Color.BLACK);
				ann.setBackgroundPaint(trans_white);
				
				List<XYTextAnnotation> annotations = Lists.newArrayList(ann);
				specs.get(i).setPlotAnnotations(annotations);
			}
		}
		
//		if (zoom <= 0)
//			zoom = data.length;
		
		XYZGraphPanel panel = new XYZGraphPanel(prefs);
		panel.drawPlot(specs, false, false, Lists.newArrayList(xRange), yRanges);
		
		if (outputFile == null) {
			// display it
			XYZPlotWindow window = new XYZPlotWindow(panel);
			window.setSize(width, (height-140)/4+140);
			window.setDefaultCloseOperation(XYZPlotWindow.EXIT_ON_CLOSE);
		} else {
			// write plot
			int myHeight = height;
			if (specs.size() > 0)
				myHeight = (height-200)*specs.size()+500;
			panel.getChartPanel().setSize(plotWidth, myHeight);
			panel.saveAsPNG(outputFile.getAbsolutePath());
//			if (plotWidth < 5000) {
//				String path = outputFile.getAbsolutePath().replaceAll(".png", "")+".pdf";
//				panel.saveAsPDF(path);
//			}
		}
		return panel;
	}
	
	private static void setRotated(EvenlyDiscrXYZ_DataSet xyz, short[][] counts, int x, int y, double val) {
		if (Double.isNaN(val))
			return;
		counts[x][y]++;
		if (plot_rotated_preserve_min)
			val = Math.min(val, xyz.get(x, y));
		else
			val = val + xyz.get(x, y);
		xyz.set(x, y, val);
	}
	
	public enum CalcMetric {
		RR("Recurrence Rate", 0.005),
		DET("Determinism", 0.005),
		LAM("Laminarity", 0.005),
		DIAG_LEN("Diagonal Length", 0.5),
		ENTROPY("Entropy", 0.05);
		
		private String name;
		private double histDelta;
		private CalcMetric(String name, double histDelta) {
			this.name = name;
			this.histDelta = histDelta;
		}
		
		@Override
		public String toString() {
			return name;
		}
	}
	
	private static double[] calcMetrics(DistanceMetric distMetric, double threshold, CalcMetric[] calcMetrics,
			List<int[]> fullPath) {
		List<double[]> normPath = calcNormalizedPath(fullPath, calcMeanRIs(fullPath));
		double[][] data = calcRotated(normPath, distMetric, calc_metrics_width);
		
		return calcMetrics(threshold, calcMetrics, data);
	}
	
	private static double[] calcMetrics(double threshold, CalcMetric[] calcMetrics,
			double[][] data) {
		int middleIndex = data[0].length/2 + 1;
		
//		// now we need to get outside the main diagonal
//		int skipStates = 10;
		
		HistogramFunction skipStatesHist = new HistogramFunction(0d, 100, 1d);
		rowLoop:
		for (double[] row : data) {
			int mySkipStates = 0;
			for (int i=middleIndex; --i>=0;) {
				double val = row[i];
				if (Double.isNaN(val))
					continue rowLoop;
				if (row[i] > threshold)
					break;
				mySkipStates = middleIndex - i;
			}
			Preconditions.checkState(mySkipStates <= skipStatesHist.getMaxX());
			skipStatesHist.add(mySkipStates, 1d);
		}
		
		// choose skip states based on histogram
		skipStatesHist.normalizeBySumOfY_Vals();
		int skipStates;
		int maxSkipStates = 15;
		double sum = 0d;
		for (skipStates=0; skipStates<skipStatesHist.size() && skipStates<maxSkipStates; skipStates++) {
			sum += skipStatesHist.getY(skipStates);
			if (sum > 0.9)
				break;
		}
		
//		if (data.length > 20000) {
//			// it's RSQSim, lets debug
//			
//			// now plot the histogram
//			List<XY_DataSet> funcs = Lists.newArrayList();
//			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
//			
//			funcs.add(skipStatesHist);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
//			
//			funcs.add(new DefaultXY_DataSet(new double[] {skipStates, skipStates}, new double[] {0d, 1d}));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//			
//			PlotSpec spec = new PlotSpec(funcs, chars, "Skip States for thresh="+threshold, "States", "Fract");
//			new GraphWindow(spec);
//		}
		
		double[] ret = new double[calcMetrics.length];
		
		int maxIndex = middleIndex - skipStates;
		
		for (int m=0; m<calcMetrics.length; m++) {
			CalcMetric calcMetric = calcMetrics[m];
			
			switch (calcMetric) {
			case RR:
				ret[m] = calcRR(data, maxIndex, threshold);
				break;
			
			case DET:
				ret[m] = calcDET(data, maxIndex, threshold);
				break;
				
			case LAM:
				ret[m] = calcLAM(data, maxIndex, threshold);
				break;
				
			case DIAG_LEN:
				ret[m] = calcDIAG_LEN(data, maxIndex, threshold);
				break;
				
			case ENTROPY:
				ret[m] = calcENTROPY(data, maxIndex, threshold);
				break;
				
			default:
				throw new IllegalStateException("Metric not implemented: "+calcMetric);
			}
		}
		
		return ret;
	}
	
	private static double calcRR(double[][] data, int maxIndex, double threshold) {
		long total = 0;
		long below = 0;
		
		rowLoop:
			for (double[] row : data) {
				int myNumBelow = 0;
				// we don't want the full row, just one half not including the origin
				// also don't want any rows with NaN's (start or end of plot) as that would bias results
				for (int i=0; i<maxIndex; i++) {
					double val = row[i];
					if (Double.isNaN(val))
						continue rowLoop;
					if (val <= threshold)
						myNumBelow++;
					//						if (i == (middleIndex - skipStates - 1))
					//							Preconditions.checkState(val > threshold);
				}
				total += maxIndex;
				below += myNumBelow;
			}

		return (double)below/(double)total;
	}
	
	private static double calcDET(double[][] data, int maxIndex, double threshold) {
		long total = 0;
		long onDiag = 0;
		for (int r=0; r<data.length; r++) {
			double[] row = data[r];
			// we don't want the full row, just one half not including the origin
			for (int i=0; i<maxIndex; i++) {
				double val = row[i];
				if (Double.isNaN(val))
					continue;
				if (val <= threshold) {
					if ((r > 0  && data[r-1][i] <= threshold)
							|| (r < data.length-1 && data[r+1][i] <= threshold))
						onDiag++;
					total++;
				}
			}
		}

		return (double)onDiag/(double)total;
	}
	
	private static double calcLAM(double[][] data, int maxIndex, double threshold) {
		if (data.length <= 1)
			return Double.NaN;
		long total = 0;
		long onDiag = 0;
		for (int r=0; r<data.length; r++) {
			double[] row = data[r];
			// we don't want the full row, just one half not including the origin
			for (int i=0; i<maxIndex; i++) {
				double val = row[i];
				if (Double.isNaN(val))
					continue;
				if (val <= threshold) {
					if ((i > 0  && data[i-1][i] <= threshold)
							|| (i < row.length-1 && data[i+1][i] <= threshold))
						onDiag++;
					total++;
				}
			}
		}

		return (double)onDiag/(double)total;
	}
	
	private static double calcDIAG_LEN(double[][] data, int maxIndex, double threshold) {
		// we'll ignore any diagonals that start or end at the beginning
		
		Map<Integer, Integer> lenCount = doCalcDiagLengths(data, maxIndex, threshold);
		
		if (lenCount.isEmpty())
			return 0d;
		
		double mean = 0d;
		int totCount = 0;
		
		for (Integer len : lenCount.keySet()) {
			int count = lenCount.get(len);
			mean += len*count;
			totCount += count;
		}
		
		mean /= (double)totCount;
		
		return mean;
	}
	
	private static Map<Integer, Integer> doCalcDiagLengths(double[][] data, int maxIndex, double threshold) {
		// we'll ignore any diagonals that start or end at the beginning
		
		Map<Integer, Integer> lenCounts = Maps.newHashMap();
		
		for (int i=0; i<maxIndex; i++) {
			int curStart = -2;
			for (int r=0; r<data.length; r++) {
				double val = data[r][i];
				if (Double.isNaN(val)) {
					curStart = -2;
					continue;
				}
				if (val <= threshold) {
					if (curStart == -1)
						// start of a new diagonal
						curStart = r;
				} else {
					if (curStart >= 0) {
						// end of a diagonal
						int len = r - curStart;
						if (!lenCounts.containsKey(len))
							lenCounts.put(len, 1);
						else
							lenCounts.put(len, lenCounts.get(len)+1);
					}
					curStart = -1;
				}
			}
		}
		
		return lenCounts;
	}

	private static double calcENTROPY(double[][] data, int maxIndex, double threshold) {
		
		Map<Integer, Integer> lenCount = doCalcDiagLengths(data, maxIndex, threshold);
		
		if (lenCount.isEmpty())
			return 0d;
		
		int totCount = 0;
		for (int len : lenCount.keySet()) {
			int count = lenCount.get(len);
			totCount += count;
		}
		
		double sum = 0;
		
		for (int len : lenCount.keySet()) {
			int count = lenCount.get(len);
			double p = (double)count/(double)totCount;
			sum += p*Math.log(p);
		}

		return -sum;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/home/kevin/Simulators/recurrence_plots");
//		File outputDir = new File("/home/kevin/Simulators/recurrence_plots/noSkip");
		
		List<SynchFaults[]> faultSets = Lists.newArrayList();
		
		faultSets.add(new SynchFaults[] {SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CARRIZO});
		faultSets.add(new SynchFaults[] {SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO});
		faultSets.add(new SynchFaults[] {SynchFaults.SAF_MOJAVE, SynchFaults.SAF_CARRIZO, SynchFaults.SAF_COACHELLA,
				SynchFaults.SAN_JACINTO});
		
		/* UCERF3 Comparisons */
		File ucerf3MainDir = new File("/home/kevin/Simulators/time_series/ucerf3_compare");
		
//		MagDependentAperiodicityOptions[] ucerf3Comparisons = {MagDependentAperiodicityOptions.MID_VALUES};
//		String[] ucerf3DirNames = {"2015_07_30-MID_VALUES"};
//		MagDependentAperiodicityOptions[] ucerf3ForCombo = ucerf3Comparisons;
		
		MagDependentAperiodicityOptions[] ucerf3Comparisons = {MagDependentAperiodicityOptions.LOW_VALUES,
				MagDependentAperiodicityOptions.MID_VALUES, MagDependentAperiodicityOptions.HIGH_VALUES, null};
		String[] ucerf3DirNames = {"2015_08_05-LOW_VALUES", "2015_07_30-MID_VALUES", "2015_08_05-HIGH_VALUES", null};
		MagDependentAperiodicityOptions[] ucerf3ForCombo = {MagDependentAperiodicityOptions.MID_VALUES};
		
//		MagDependentAperiodicityOptions[] ucerf3Comparisons = null;
//		String[] ucerf3DirNames = null;
//		MagDependentAperiodicityOptions[] ucerf3ForCombo = ucerf3Comparisons;
		
		File[] etasCatalogs = {
//				new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//						+ "2015_11_30-spontaneous-1000yr-FelzerParams-mc20-full_td-noApplyLTR/results_m4.bin")
//						+ "2016_01_05-spontaneous-10000yr-mc10-applyGrGridded-full_td-noApplyLTR/results_m4.bin")
//						+ "2016_02_04-spontaneous-10000yr-full_td-subSeisSupraNucl-gridSeisCorr/results_m4.bin")
		};
		
		
		HashSet<MagDependentAperiodicityOptions> ucerf3ForComboSet = new HashSet<MagDependentAperiodicityOptions>();
		if (ucerf3ForCombo != null)
			for (MagDependentAperiodicityOptions cov : ucerf3ForCombo)
				ucerf3ForComboSet.add(cov);
		
		double distSpacing = 10d;
		boolean normalize = true;
		
		boolean plotRecurrence = true;
		boolean doMomRateComparison = false;
		boolean doCalcMetrics = false;
		boolean plotSpecialPosterFigs = false;
		EvenlyDiscretizedFunc threshFuncXVals = new EvenlyDiscretizedFunc(0.1, 10, 0.1);
		
		boolean[] poissons = { false, true };
		
//		CalcMetric[] calcMetrics = CalcMetric.values();
		CalcMetric[] calcMetrics = { CalcMetric.RR };
		
		List<Color> colors = Lists.newArrayList();
		
		for (boolean poisson : poissons) {
			if (poisson)
				colors.add(Color.BLUE);
			else
				colors.add(Color.RED);
		}
		if (ucerf3Comparisons != null) {
			for (MagDependentAperiodicityOptions cov : ucerf3Comparisons) {
				if (cov == null)
					colors.add(Color.CYAN);
				else switch (cov) {
				case LOW_VALUES:
					colors.add(Color.GREEN);
					break;
				case MID_VALUES:
					colors.add(Color.GREEN.darker());
					break;
				case HIGH_VALUES:
					colors.add(Color.GREEN.darker().darker());
					break;

				default:
					throw new IllegalStateException("Unknown cov: "+cov);
				}
			}
		}
		if (etasCatalogs != null)
			for (int i=0; i<etasCatalogs.length; i++)
				colors.add(Color.ORANGE);
		
		SimAnalysisCatLoader loader = new SimAnalysisCatLoader(true, null, true);
		List<? extends SimulatorEvent> events = loader.getEvents();
		List<SimulatorElement> elems = loader.getElements();
		
		List<List<List<SimulatorEvent>>> ucerf3Catalogs = null;
		FaultSystemSolution ucerf3Sol = null;
		Map<Integer, SimulatorElement> ucerf3Elems = null;
		int ucerf3StartYear = 2014;
		int numUCERF3Catalogs = 5;
		int totNumUCERF3 = 0;
		if (ucerf3Comparisons == null)
			ucerf3Comparisons = new MagDependentAperiodicityOptions[0];
		totNumUCERF3 = ucerf3Comparisons.length;
		if (etasCatalogs != null && etasCatalogs.length > 0)
			totNumUCERF3 += etasCatalogs.length;
		if (totNumUCERF3 > 0) {
			ucerf3Catalogs = Lists.newArrayList();
			ucerf3Sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/"
					+ "UCERF3/data/scratch/InversionSolutions/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
			ucerf3Elems = UCERF3ComparisonAnalysis.loadElements(ucerf3Sol.getRupSet());
			for (int i=0; i<totNumUCERF3; i++) {
				if (i >= ucerf3Comparisons.length) {
					// etas
					File catFile = etasCatalogs[i-ucerf3Comparisons.length];
					List<List<ETAS_EqkRupture>> etasCats = ETAS_CatalogIO.loadCatalogsBinary(catFile, 6d);
					List<List<SimulatorEvent>> myCatalogs = UCERF3_ETASComparisons.loadUCERF3EtasCatalogs(
							etasCats, ucerf3Sol, null, ucerf3Elems);
					
					ucerf3Catalogs.add(myCatalogs);
				} else if (ucerf3Comparisons[i] == null) {
					// poisson
					List<List<SimulatorEvent>> myCatalogs = Lists.newArrayList();
					int num;
					if (doCalcMetrics)
						num = 80;
					else
						num = numUCERF3Catalogs;
					
					int len = 10000;
					
					File ucerf3Dir = new File(ucerf3MainDir, ".poisson_cache");
					
					for (int j=0; j<num; j++) {
						File binFile = new File(ucerf3Dir, "len_"+len+"_cat_"+j+".bin");
						List<SimulatorEvent> u3Events;
						if (binFile.exists()) {
							u3Events = UCERF3ComparisonAnalysis.loadUCERF3CatalogBinary(ucerf3Sol, ucerf3Elems, binFile);
						} else {
							u3Events = UCERF3ComparisonAnalysis.calcPoissonCatalog(ucerf3Sol, ucerf3Elems, len);
							UCERF3ComparisonAnalysis.writeUCERF3CatalogBinary(u3Events, binFile);
						}
						myCatalogs.add(u3Events);
					}
					ucerf3Catalogs.add(myCatalogs);
				} else {
					File ucerf3Dir = new File(ucerf3MainDir, ucerf3DirNames[i]);
					List<List<SimulatorEvent>> myCatalogs = UCERF3ComparisonAnalysis.loadUCERF3Catalogs(
							ucerf3Dir, ucerf3Sol, null, ucerf3Elems, ucerf3StartYear);
					
					ucerf3Catalogs.add(myCatalogs);
				}
			}
		}
		
		List<DistanceMetric> metrics = Lists.newArrayList();
		List<Double> thresholds = Lists.newArrayList();
		
		populateDefaultMetrics(metrics, thresholds, normalize);
		
		for (SynchFaults[] faults : faultSets) {
			List<RuptureIdentifier> rupIdens = SynchIdens.getIndividualFaults(7, 10d, faults);
			List<String> idenNames = Lists.newArrayList();
			for (RuptureIdentifier iden : rupIdens)
				idenNames.add(iden.getName());
			
			String faultNames = Joiner.on("_").join(Lists.newArrayList(faults));
			
			// for combined plots
			List<List<int[]>> comboPlotPaths = Lists.newArrayList();
			List<List<? extends SimulatorEvent>> comboEventLists = null;
			List<Double> comboStartTimesForComparison = null;
			if (doMomRateComparison) {
				comboEventLists = Lists.newArrayList();
				comboStartTimesForComparison = Lists.newArrayList();
			}
			List<String> comboPlotNames = Lists.newArrayList();
			
			// for threshold histogram
			Map<CalcMetric, List<double[]>> threshFracts = null;
			Map<CalcMetric, List<EvenlyDiscretizedFunc>> threshFuncs = null;
			Map<CalcMetric, List<EvenlyDiscretizedFunc>> metricTimeSeries = null;
			List<String> threshNames = null;
			if (doCalcMetrics) {
				threshFracts = Maps.newHashMap();
				threshFuncs = Maps.newHashMap();
				metricTimeSeries = Maps.newHashMap();
				for (CalcMetric metric : calcMetrics) {
					threshFracts.put(metric, new ArrayList<double[]>());
					threshFuncs.put(metric, new ArrayList<EvenlyDiscretizedFunc>());
					metricTimeSeries.put(metric, new ArrayList<EvenlyDiscretizedFunc>());
				}
				threshNames = Lists.newArrayList();
			}
			
			for (boolean poisson : poissons) {
				File subDir;
				if (poisson)
					subDir = new File(outputDir, "rsqsim_poisson");
				else
					subDir = new File(outputDir, "rsqsim_standard");
				Preconditions.checkState(subDir.exists() || subDir.mkdir());
				
				List<? extends SimulatorEvent> myEvents;
				
				if (poisson) {
					List<SimulatorEvent> trueRandom = Lists.newArrayList();
					double start = events.get(0).getTime();
					double end = events.get(events.size()-1).getTime();
					int id = 0;
					for (SimulatorEvent e : events)
						trueRandom.add(e.cloneNewTime((end-start)*Math.random()+start, id++));
					Collections.sort(trueRandom);
					myEvents = trueRandom;
				} else {
					myEvents = events;
				}
				
				File mySubDir = new File(subDir, faultNames);
				Preconditions.checkState(mySubDir.exists() || mySubDir.mkdir());
				
				List<List<? extends SimulatorEvent>> eventsForStatesList = null;
				if (doMomRateComparison)
					eventsForStatesList = Lists.newArrayList();
				
				List<int[]> fullPath = MarkovChainBuilder.getStatesPath(
						distSpacing, MarkovChainBuilder.getMatchesLists(myEvents, rupIdens), 0d, eventsForStatesList);
				// skip the first 10,000 years while things are getting ramped up
				int pathStartIndex = 1000;
//				int pathStartIndex = 0;
				fullPath = fullPath.subList(pathStartIndex, fullPath.size());
					
//				System.out.println("Path length: "+fullPath.size());
				
				if (plotRecurrence)
					plotRecurrence(mySubDir, fullPath, distSpacing, normalize, metrics, thresholds, rupIdens);
				
				if (!poisson && plotSpecialPosterFigs) {
					// plot standard for poster
					File outputFile = new File(mySubDir, "orig_style.png");
					List<int[]> myPath = fullPath.subList(0, 500);
					List<double[]> normPath = calcNormalizedPath(myPath, calcMeanRIs(myPath));
					DistanceMetric metric = metrics.get(0);
					double threshold = thresholds.get(0);
					BitSet[] data = calcBitSet(normPath, metric, threshold);
					plotDiscrete(data, metric, threshold, outputFile, -1, distSpacing);
					data = null;
					double[][] rotated = calcRotated(normPath, metric, rotated_width);
					outputFile = new File(mySubDir, "orig_style_rotated.png");
					plotDiscreteRotated(rotated, metric, rotated_pixel_width,
							threshold, distSpacing, outputFile);
					outputFile = new File(mySubDir, "orig_style_rotated_hybrid.png");
					plotHybridRotated(rotated, metric, rotated_pixel_width,
							threshold, 5d*threshold, distSpacing, outputFile);
					
					// now a huuuuuge one for the poster
					myPath = fullPath.subList(0, 10000);
					normPath = calcNormalizedPath(myPath, calcMeanRIs(myPath));
					rotated = calcRotated(normPath, metric, 3*(rotated_width-1)+1);
					outputFile = new File(mySubDir, "huge_for_poster.png");
					plotHybridRotated(rotated, metric, rotated_pixel_width,
							threshold, 5d*threshold, distSpacing, outputFile);
				}
				
				String name;
				if (poisson)
					name = "RSQSim Poisson";
				else
					name = "RSQSim";
				comboPlotPaths.add(fullPath);
				comboPlotNames.add(name);
				
				if (doMomRateComparison) {
					eventsForStatesList = eventsForStatesList.subList(pathStartIndex, eventsForStatesList.size());
					double startTime = getStartYearForMomRateComparison(distSpacing, myEvents, eventsForStatesList);
					comboEventLists.add(myEvents);
					comboStartTimesForComparison.add(startTime);
				}
				
				if (doCalcMetrics) {
					List<double[]> normPath = calcNormalizedPath(fullPath, calcMeanRIs(fullPath));
					double[][] data = calcRotated(normPath, metrics.get(0), calc_metrics_width);
					
					double[] result = calcMetrics(thresholds.get(0), calcMetrics, data);
					List<double[]> funcResults = Lists.newArrayList();
					for (int t=0; t<threshFuncXVals.size(); t++)
						funcResults.add(calcMetrics(threshFuncXVals.getX(t), calcMetrics, data));
					for (int m=0; m<calcMetrics.length; m++) {
						threshFracts.get(calcMetrics[m]).add(new double[] {result[m]});
						EvenlyDiscretizedFunc threshFunc = new EvenlyDiscretizedFunc(
								threshFuncXVals.getMinX(), threshFuncXVals.size(), threshFuncXVals.getDelta());
						for (int t=0; t<threshFuncXVals.size(); t++)
							threshFunc.set(t, funcResults.get(t)[m]);
						threshFuncs.get(calcMetrics[m]).add(threshFunc);
					}
					
					// now time series
					EvenlyDiscretizedFunc[] timeSeries = new EvenlyDiscretizedFunc[calcMetrics.length];
					for (int m=0; m<calcMetrics.length; m++) {
						timeSeries[m] = new EvenlyDiscretizedFunc(0d, data.length, distSpacing);
						metricTimeSeries.get(calcMetrics[m]).add(timeSeries[m]);
					}
					for (int i=0; i<data.length; i++) {
						double[][] subData = {data[i]};
						result = calcMetrics(thresholds.get(0), calcMetrics, subData);
						for (int m=0; m<calcMetrics.length; m++)
							timeSeries[m].set(i, result[m]);
					}
					
					threshNames.add(name);
				}
			}
			
			if (ucerf3Catalogs != null) {
				List<RuptureIdentifier> ucerf3Idens = UCERF3ComparisonAnalysis.buildUCERF3_EquivIdens(
						rupIdens, elems, ucerf3Elems, ucerf3Sol.getRupSet());
				for (int i=0; i<ucerf3Catalogs.size(); i++) {
					List<List<SimulatorEvent>> catalogs = ucerf3Catalogs.get(i);
					
					MagDependentAperiodicityOptions cov;
					if (i < ucerf3Comparisons.length)
						cov = ucerf3Comparisons[i];
					else
						cov = null;
					String name;
					String comboName;
					if (cov == null) {
						if (i >= ucerf3Comparisons.length)
							name = "UCERF3-ETAS";
						else
							name = "UCERF3 Poisson";
						comboName = name;
					} else {
						name = "UCERF3 "+cov.name().replaceAll("_", " ")
							.replaceAll("VALUES", "COV");
						comboName = "UCERF3-TD";
					}
					
					File subDir;
					if (cov == null)
						if (i >= ucerf3Comparisons.length)
							subDir = new File(outputDir, "ucerf3_etas");
						else
							subDir = new File(outputDir, "ucerf3_poisson");
					else
						subDir = new File(outputDir, "ucerf3_"+cov.name());
					Preconditions.checkState(subDir.exists() || subDir.mkdir());
					
					subDir = new File(subDir, faultNames);
					if (subDir.exists())
						FileUtils.deleteRecursive(subDir);
					Preconditions.checkState(subDir.mkdir());
					
					List<SimulatorEvent> preEvents = null;
					int skipStates = 0;
					
					if (cov != null || i>= ucerf3Comparisons.length) {
						preEvents = UCERF3ComparisonAnalysis.getFakePreEvents(
								ucerf3Sol.getRupSet(), ucerf3Idens, ucerf3Elems, ucerf3StartYear);
						double oiBefore = preEvents.get(preEvents.size()-1).getTimeInYears();
						Preconditions.checkState(oiBefore < 0);
						skipStates = (int)(-oiBefore/distSpacing + 0.5);
						System.out.println("Fake event OI established at "+oiBefore
								+" years, skipping "+skipStates+" states");
					}
					
					List<List<int[]>> catalogPaths = Lists.newArrayList();
					List<String> catalogNames = Lists.newArrayList();
					List<List<? extends SimulatorEvent>> catalogEventLists = null;
					List<Double> catalogStartTimesForComparison = null;
					if (doMomRateComparison) {
						catalogEventLists = Lists.newArrayList();
						catalogStartTimesForComparison = Lists.newArrayList();
					}
					
					if (doCalcMetrics) {
						List<double[]> results = Lists.newArrayList();
						double catRate = 1d/catalogs.size();
						EvenlyDiscretizedFunc[] threshFunc = new EvenlyDiscretizedFunc[calcMetrics.length];
						for (int m=0; m<calcMetrics.length; m++)
								threshFunc[m] = new EvenlyDiscretizedFunc(
										threshFuncXVals.getMinX(), threshFuncXVals.size(), threshFuncXVals.getDelta());
						double[][] data = null;
						for (List<SimulatorEvent> catalog : catalogs) {
							catalog = Lists.newArrayList(catalog);
							if (preEvents != null)
								catalog.addAll(0, preEvents);
							List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing,
									MarkovChainBuilder.getMatchesLists(catalog, ucerf3Idens), 0d, null);
							List<double[]> normPath = calcNormalizedPath(fullPath, calcMeanRIs(fullPath));
							double[][] myData = calcRotated(normPath, metrics.get(0), calc_metrics_width);
							results.add(calcMetrics(thresholds.get(0), calcMetrics, myData));
							for (int t=0; t<threshFuncXVals.size(); t++)  {
								double[] ret = calcMetrics(threshFuncXVals.getX(t), calcMetrics, myData);
								for (int m=0; m<calcMetrics.length; m++)
									threshFunc[m].add(t, catRate*ret[m]);
							}
							if (data == null)
								data = myData;
						}
						for (int m=0; m<calcMetrics.length; m++) {
							double[] fracts = new double[results.size()];
							for (int c=0; c<fracts.length; c++)
								fracts[c] = results.get(c)[m];
							threshFracts.get(calcMetrics[m]).add(fracts);
							threshFuncs.get(calcMetrics[m]).add(threshFunc[m]);
						}
						
						// now time series
						EvenlyDiscretizedFunc[] timeSeries = new EvenlyDiscretizedFunc[calcMetrics.length];
						for (int m=0; m<calcMetrics.length; m++) {
							timeSeries[m] = new EvenlyDiscretizedFunc(0d, data.length, distSpacing);
							metricTimeSeries.get(calcMetrics[m]).add(timeSeries[m]);
						}
						for (int j=0; j<data.length; j++) {
							double[][] subData = {data[j]};
							double[] result = calcMetrics(thresholds.get(0), calcMetrics, subData);
							for (int m=0; m<calcMetrics.length; m++)
								timeSeries[m].set(i, result[m]);
						}
						
						threshNames.add(name);
					}
					
					// now sort by duration
					List<Double> durations = Lists.newArrayList();
					for (List<SimulatorEvent> catalog : catalogs)
						durations.add(UCERF3ComparisonAnalysis.getDurationYears(catalog));
					
					// sort by duration decreasing
					catalogs = ComparablePairing.getSortedData(durations, catalogs);
					Collections.reverse(catalogs);
					
					if (catalogs.size() > numUCERF3Catalogs)
						catalogs = catalogs.subList(0, numUCERF3Catalogs);
					
					for (int j=0; j<catalogs.size(); j++) {
						List<SimulatorEvent> catalog = catalogs.get(j);
//						double origCatLen = UCERF3ComparisonAnalysis.getDurationYears(catalog);
						if (preEvents != null) {
							catalog.addAll(0, preEvents);
//							double newCatLen = catalog.get(catalog.size()-1).getTimeInYears() - catalog.get(0).getTimeInYears();
//							double addedYears = newCatLen - origCatLen;
//							int approxAddedStates = (int)(addedYears/distSpacing + 0.5);
//							System.out.println("Catalog "+j+" lenght: "+origCatLen);
//							System.out.println("Length with pre events: "+newCatLen+", added "
//									+addedYears+" ("+approxAddedStates+" states)");
						}
						
						List<List<? extends SimulatorEvent>> eventsForStatesList = null;
						if (doMomRateComparison) {
							eventsForStatesList = Lists.newArrayList();
						}
						
						List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing,
								MarkovChainBuilder.getMatchesLists(catalog, ucerf3Idens), 0d, eventsForStatesList);
						fullPath = fullPath.subList(skipStates, fullPath.size());
						
						catalogPaths.add(fullPath);
						catalogNames.add(name+" Catalog "+j);
						
						// == here includes the first ETAS catalog if present 
						if (j == 0 &&
								(i < ucerf3Comparisons.length && ucerf3ForComboSet.contains(cov)
										|| i == ucerf3Comparisons.length)) {
//						if (j == 0 && ucerf3ForComboSet.contains(cov)) {
							comboPlotPaths.add(fullPath);
							comboPlotNames.add(comboName);
						}
						
						if (doMomRateComparison) {
							eventsForStatesList = eventsForStatesList.subList(skipStates, eventsForStatesList.size());
							double startTime = getStartYearForMomRateComparison(distSpacing, catalog, eventsForStatesList);
							catalogEventLists.add(catalog);
							catalogStartTimesForComparison.add(startTime);
							if (j == 0) {
								comboEventLists.add(catalog);
								comboStartTimesForComparison.add(startTime);
							}
						}
						
						File mySubDir;
						if (catalogs.size() > 1) {
							mySubDir = new File(subDir, "catalog"+j);
							Preconditions.checkState(mySubDir.exists() || mySubDir.mkdir());
						} else {
							mySubDir = subDir;
						}
						
//						fullPath = fullPath.subList(0, maxStates);
						
						if (plotRecurrence)
							plotRecurrence(mySubDir, fullPath, distSpacing, normalize,
								metrics, thresholds, rupIdens);
					}
					// now all catalogs combined
					if (catalogPaths.size() > 1 && plotRecurrence) {
						for (int m=0; m<metrics.size(); m++) {
							if (m == 0)
								plotMultiRecurrence(distSpacing, normalize, metrics.get(m), thresholds.get(m), catalogPaths,
									catalogNames, catalogEventLists, catalogStartTimesForComparison, subDir, idenNames);
							else
								plotMultiRecurrence(distSpacing, normalize, metrics.get(m), thresholds.get(m), catalogPaths,
										catalogNames, null, null, subDir, idenNames);
						}
					}
				}
			}
			
			// now combo plots
			File comboDir = new File(outputDir, "combined_plots");
			Preconditions.checkState(comboDir.exists() || comboDir.mkdir());
			comboDir = new File(comboDir, faultNames);
			Preconditions.checkState(comboDir.exists() || comboDir.mkdir());
			
			for (int i=0; i<metrics.size() && plotRecurrence; i++) {
				if (i == 0)
					plotMultiRecurrence(distSpacing, normalize, metrics.get(i), thresholds.get(i), comboPlotPaths,
						comboPlotNames, comboEventLists, comboStartTimesForComparison, comboDir, idenNames);
				else
					plotMultiRecurrence(distSpacing, normalize, metrics.get(i), thresholds.get(i), comboPlotPaths,
							comboPlotNames, null, null, comboDir, idenNames);
			}
			
			// now thresh fract plots
			if (doCalcMetrics) {
				for (CalcMetric calcMetric : calcMetrics) {
					plotMetricHists(metrics.get(0), thresholds.get(0), calcMetric,
							threshFracts.get(calcMetric), threshNames, colors, comboDir);
					plotMetricFuncs(metrics.get(0), calcMetric, threshFuncs.get(calcMetric),
							threshNames, colors, comboDir);
					plotMetricTimeSeries(metrics.get(0), calcMetric, metricTimeSeries.get(calcMetric), thresholds.get(0),
							threshNames, colors, comboDir);
				}
			}
		}
	}
	
	private static void plotMetricHists(DistanceMetric distMetric, double threshold, CalcMetric calcMetric,
			List<double[]> fracts, List<String> names, List<Color> colors, File outputDir) throws IOException {
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		double histDelta = calcMetric.histDelta;
		
		double maxX = 0;
		
		// first hists
		for (int i=fracts.size(); --i>=0;) {
			double[] fract = fracts.get(i);
			Color color = colors.get(i);
			
			if (fract.length > 1) {
				// do a histogram
				HistogramFunction hist = HistogramFunction.getEncompassingHistogram(
						StatUtils.min(fract)-histDelta, StatUtils.max(fract)+histDelta, histDelta);
				
				double weight = 1d/fract.length;
				for (double val : fract)
					if (Doubles.isFinite(val))
						hist.add(val, weight);
				
				maxX = Math.max(maxX, hist.getMaxX()+0.5*histDelta);
				
				funcs.add(hist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
			}
		}
		
		DecimalFormat df = new DecimalFormat("0.0000");
		
		// now mean lines
		for (int i=fracts.size(); --i>=0;) {
			double[] fract = fracts.get(i);
			String name = names.get(i);
			Color color = colors.get(i);
			
			double mean = StatUtils.mean(fract);
			maxX = Math.max(maxX, mean*1.1);
			XY_DataSet meanXY = new DefaultXY_DataSet();
			meanXY.setName(name+" ("+df.format(mean)+")");
			meanXY.set(mean, 0d);
			meanXY.set(mean, 0.5);
			meanXY.set(mean, 1d);
			funcs.add(meanXY);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
		}
		
		String title = calcMetric.toString()+", "+distMetric.name()+" ≤ "+(float)threshold;
		PlotSpec spec = new PlotSpec(funcs, chars, title, calcMetric.toString(), "Density");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, new Range(0, maxX), new Range(0, 1));
		gp.getChartPanel().setSize(800, 700);
		String prefix = distMetric.name()+"_"+calcMetric.name()+"_hist";
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	private static void plotMetricFuncs(DistanceMetric metric, CalcMetric calcMetric, List<EvenlyDiscretizedFunc> funcs,
			List<String> names, List<Color> colors, File outputDir) throws IOException {
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		for (int i=0; i<funcs.size(); i++) {
			funcs.get(i).setName(names.get(i));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors.get(i)));
		}
		
		Collections.reverse(funcs);
		Collections.reverse(chars);
		
		String title = calcMetric.toString()+", "+metric.name()+" ≤ Threshold";
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Threshold", calcMetric.toString());
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		EvenlyDiscretizedFunc func0 = funcs.get(0);
		gp.drawGraphPanel(spec, false, false, new Range(func0.getMinX(), func0.getMaxX()), null);
		gp.getChartPanel().setSize(800, 700);
		String prefix = metric.name()+"_"+calcMetric.name()+"_func";
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	private static void plotMetricTimeSeries(DistanceMetric metric, CalcMetric calcMetric, List<EvenlyDiscretizedFunc> funcs,
			double threshold, List<String> names, List<Color> colors, File outputDir) throws IOException {
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		for (int i=0; i<funcs.size(); i++) {
			funcs.get(i).setName(names.get(i));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors.get(i)));
		}
		
		Collections.reverse(funcs);
		Collections.reverse(chars);
		
		String title = calcMetric.toString()+" Time Series, "+metric.name()+" ≤ "+(float)threshold;
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Years", calcMetric.toString());
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		EvenlyDiscretizedFunc func0 = funcs.get(0);
		gp.drawGraphPanel(spec, false, false, new Range(func0.getMinX(), func0.getMaxX()), null);
		gp.getChartPanel().setSize(800, 700);
		String prefix = metric.name()+"_"+calcMetric.name()+"_ts";
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
	}
	
	/**
	 * Aligns catalog to recurrence plot
	 * @return
	 */
	private static double getStartYearForMomRateComparison(double distSpacing,
			List<? extends SimulatorEvent> catalog, List<List<? extends SimulatorEvent>> eventsForStatesList) {
		// first we need to detect the start time. this is <= to the first event time in eventsForStatesList
		
		// start by assuming that it is the time of the first event
		// events for each state are already sorted
		double startTime = Double.NaN;
		for (int i=0; i<eventsForStatesList.size(); i++) {
			if (eventsForStatesList.get(i).isEmpty())
				continue;
			startTime = eventsForStatesList.get(i).get(0).getTimeInYears();
			// if this isn't the first bin (first bin was empty), adjust
			startTime -= distSpacing*i;
			break;
		}
		Preconditions.checkState(!Double.isNaN(startTime));
		
		for (int i=1; i<eventsForStatesList.size(); i++) {
			if (eventsForStatesList.get(i).isEmpty())
				continue;
			double offset = i*distSpacing;
			double windowStart = startTime + offset;
			SimulatorEvent firstEvent = eventsForStatesList.get(i).get(0);
			double t = firstEvent.getTimeInYears();
			if (t < windowStart) {
				// our start time is off, adjust
				windowStart = t;
				startTime = windowStart - offset;
			}
		}
		
		// now validate
		for (int i=1; i<eventsForStatesList.size(); i++) {
			if (eventsForStatesList.get(i).isEmpty())
				continue;
			double offset = i*distSpacing;
			double windowStart = startTime + offset;
			double windowEnd = windowStart + distSpacing;
			for (SimulatorEvent e : eventsForStatesList.get(i)) {
				double t = e.getTimeInYears();
				Preconditions.checkState(t >= windowStart && t <= windowEnd,
						"Bad window detection at state %s. t=%s not in [%s %s]", i, t, windowStart, windowEnd);
			}
		}
		
		return startTime;
	}

	private static void plotMultiRecurrence(double distSpacing, boolean normalize, DistanceMetric metric,
			double threshold, List<List<int[]>> comboPlotPaths, List<String> comboPlotNames,
			List<List<? extends SimulatorEvent>> catalogs, List<Double> startTimesForComparison, File comboDir,
			List<String> idenNames) throws IOException {
		int startIndex = 0;
		int rotatedEndIndex = 500;
		for (List<int[]> fullPath : comboPlotPaths)
			if (fullPath.size() < rotatedEndIndex)
				rotatedEndIndex = fullPath.size();
		int plotStatesLen = rotatedEndIndex - startIndex;
		String threshStr = "thresh";
		if (normalize)
			threshStr += "Norm";
		double maxZ = threshold*5d;
		threshStr += (float)threshold;
		
		CPT cpt = getHybridCPT(threshold, maxZ);
		
		List<double[][]> datas = Lists.newArrayList();
		
		for (int j=0; j<comboPlotPaths.size(); j++) {
			List<int[]> fullPath = comboPlotPaths.get(j);
			
			List<double[]> rotatedSubPath = calcNormalizedPath(
					fullPath.subList(startIndex, rotatedEndIndex), calcMeanRIs(fullPath));
			
			// now rotated
			double[][] rotatedData = calcRotated(rotatedSubPath, metric, rotated_width);
			datas.add(rotatedData);
		}
		
		String name = "rp_"+metric.name()+"_hybrid_"+threshStr+"_rotated.png";
		XYZGraphPanel xyzGP = plotRotated(datas, comboPlotNames, metric, cpt, rotated_pixel_width, distSpacing,
				new File(comboDir, name), comboPlotPaths, idenNames);
		
		if (catalogs != null) {
			int windowLen = 25;
			double[] taper = SimulatorMomRateVarCalc.buildHanningTaper(windowLen);
			double halfDelta = distSpacing*0.5;
			List<PlotSpec> specs = Lists.newArrayList();
			Range xRange = null;
			List<Range> yRanges = Lists.newArrayList();
			Region region = new CaliforniaRegions.RELM_SOCAL();
			for (int i=0; i<comboPlotPaths.size(); i++) {
				double startTimeForComparison = startTimesForComparison.get(i);
				// the moment rate binning puts all events that happen before the given year in that bin
				// so construct year away using bin edges instead of bin starts (startTime)
				startTimeForComparison += distSpacing;
				double[] years = new double[plotStatesLen];
				for (int j=0; j<plotStatesLen; j++)
					years[j] = startTimeForComparison + j*distSpacing;
				// filter to So Cal
				RegionIden iden = new RegionIden(region);
				List<? extends SimulatorEvent> events = iden.getMatches(catalogs.get(i));
				double[] momRates = SimulatorMomRateVarCalc.calcTaperedMomRates(events, years, taper);
				
				EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(halfDelta, years.length, distSpacing);
				if (xRange == null)
					xRange = new Range(0, func.getMaxX()+halfDelta);
				for (int j=0; j<func.size(); j++)
					func.set(j, momRates[j]);
				yRanges.add(new Range(func.getMinY()*0.9, func.getMaxY()*1.1));
				List<DiscretizedFunc> funcs = Lists.newArrayList();
				List<PlotCurveCharacterstics> chars = Lists.newArrayList();
				funcs.add(func);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
				PlotSpec spec = new PlotSpec(funcs, chars, "Moment Rate Variation", "Years", "Moment Rate");
				double x = xRange.getLowerBound() + 0.05*(xRange.getUpperBound()-xRange.getLowerBound());
				Range yRange = yRanges.get(i);
				double y = yRange.getLowerBound() + 0.9*(yRange.getUpperBound()-yRange.getLowerBound());
				XYTextAnnotation ann = new XYTextAnnotation(comboPlotNames.get(i), x, y);
				ann.setFont(annotation_font);
				ann.setTextAnchor(TextAnchor.TOP_LEFT);
				ann.setPaint(Color.BLACK);
				ann.setBackgroundPaint(trans_white);
				
				List<XYTextAnnotation> annotations = Lists.newArrayList(ann);
				spec.setPlotAnnotations(annotations);
				
				specs.add(spec);
			}
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(prefs);
			gp.setBackgroundColor(Color.WHITE);
			
			List<Range> xRanges = Lists.newArrayList(xRange);
			gp.drawGraphPanel(specs, false, true, xRanges, yRanges);
			gp.getChartPanel().setSize(xyzGP.getChartPanel().getWidth(), (int)(0.6*xyzGP.getChartPanel().getHeight()));
			gp.saveAsPNG(new File(comboDir, "mom_rate_comparison.png").getAbsolutePath());
			gp.saveAsPDF(new File(comboDir, "mom_rate_comparison.pdf").getAbsolutePath());
		}
	}
	
	public static List<double[]> toDoublePath(List<int[]> fullPath) {
		List<double[]> doubles = Lists.newArrayList();
		
		for (int[] state : fullPath) {
			double[] newState = new double[state.length];
			for (int i=0; i<state.length; i++)
				newState[i] = state[i];
			doubles.add(newState);
		}
		
		return doubles;
	}
	
	public static double[] calcMeanRIs(List<int[]> path) {
		List<List<Double>> ris = Lists.newArrayList();
		
		for (int i=0; i<path.get(0).length; i++)
			ris.add(new ArrayList<Double>());
		
		for (int i=1; i<path.size(); i++) {
			int[] state = path.get(i);
			for (int j=0; j<state.length; j++) {
				if (state[j] == 0) {
					// just ruptured, add previous state
					ris.get(j).add((double)path.get(i-1)[j]);
				}
			}
		}
		
		double[] meanRIs = new double[ris.size()];
		
		for (int i=0; i<meanRIs.length; i++)
			meanRIs[i] = StatUtils.mean(Doubles.toArray(ris.get(i)));
		
		return meanRIs;
	}
	
	public static List<double[]> calcNormalizedPath(List<int[]> path, double[] meanRIs) {
		List<double[]> ret = Lists.newArrayList();
		
		for (int[] state : path) {
			double[] doubles = new double[state.length];
			
			for (int i=0; i<state.length; i++)
				doubles[i] = state[i]/meanRIs[i];
			
			ret.add(doubles);
		}
		
		return ret;
	}
	
	private static void populateDefaultMetrics(List<DistanceMetric> metrics, List<Double> thresholds, boolean normalize) {
		// units here are un normalized, will be converted assuming meanRI = 100;
		
		metrics.add(DistanceMetric.L1_NORM);
		thresholds.add(5d);
		
		metrics.add(DistanceMetric.LINFINITY_NORM);
		thresholds.add(2d);

		//				metrics.add(DistanceMetric.L2_NORM);
		//				thresholds.add(100d);
		
		if (normalize) {
			for (int i=0; i<thresholds.size(); i++)
				thresholds.set(i, 0.1*thresholds.get(i));
		}
	}
	
	private static final int rotated_width = 151;
	private static final int calc_metrics_width = 5*(rotated_width-1)+1;
	private static final int rotated_pixel_width = 9;

	public static void plotRecurrence(File outputDir, List<int[]> fullPath, double distSpacing,
			boolean normalize) throws IOException {
		List<DistanceMetric> metrics = Lists.newArrayList();
		List<Double> thresholds = Lists.newArrayList();
		
		populateDefaultMetrics(metrics, thresholds, normalize);
		
		plotRecurrence(outputDir, fullPath, distSpacing, normalize, metrics, thresholds, null);
	}

	public static void plotRecurrence(File outputDir, List<int[]> fullPath, double distSpacing,
			boolean normalize, List<DistanceMetric> metrics, List<Double> thresholds,
			List<RuptureIdentifier> idens) throws IOException {
//		File rotatedDir = new File(outputDir, "rotated");
//		Preconditions.checkState(rotatedDir.exists() || rotatedDir.mkdir());
		
		int startIndex = 0;
//		int squareEndIndex = 1000;
//		int squareZoomLevel = 300;
		int rotatedEndIndex = 3000;
		
//		if (squareEndIndex >= fullPath.size())
//			squareEndIndex = fullPath.size()-1;
		if (rotatedEndIndex >= fullPath.size())
			rotatedEndIndex = fullPath.size()-1;
		
//		int rotatedWidth = 51;
		
//		List<int[]> squareSubPath = fullPath.subList(startIndex, squareEndIndex);
//		List<int[]> rotatedSubPath = fullPath.subList(startIndex, rotatedEndIndex);
		
//		List<double[]> squareSubPath;
		List<double[]> rotatedSubPath;
		
		if (normalize) {
			double[] ris = calcMeanRIs(fullPath);
//			squareSubPath = calcNormalizedPath(fullPath.subList(startIndex, squareEndIndex), ris);
			rotatedSubPath = calcNormalizedPath(fullPath.subList(startIndex, rotatedEndIndex), ris);
		} else {
//			squareSubPath = toDoublePath(fullPath.subList(startIndex, squareEndIndex));
			rotatedSubPath = toDoublePath(fullPath.subList(startIndex, rotatedEndIndex));
		}
		
		for (int i=0; i<metrics.size(); i++) {
			DistanceMetric metric = metrics.get(i);
			double threshold = thresholds.get(i);
			
//			double[][] dists = calcDist(squareSubPath, metric);
//			BitSet[] bitSets = distToBitSet(dists, threshold);
			
			String threshStr = "thresh";
			if (normalize)
				threshStr += "Norm";
			threshStr += (float)threshold;
			
			double maxZ = threshold*5d;
			
//			for (int zoom : new int[] { -1, squareZoomLevel }) {
//				
//				String name = "rp_"+metric.name()+"_"+threshStr;
//				if (zoom > 0)
//					name += "_zoom"+zoom;
//				name += ".png";
//				plotDiscrete(bitSets, metric, threshold, new File(outputDir, name), zoom, distSpacing);
//				
//				name = "rp_"+metric.name()+"_continuous";
//				if (zoom > 0)
//					name += "_zoom"+zoom;
//				name += ".png";
//				plotContinuous(dists, metric, maxZ, new File(outputDir, name), zoom, distSpacing);
//				
//				String name = "rp_"+metric.name()+"_hybird_"+threshStr;
//				if (zoom > 0)
//					name += "_zoom"+zoom;
//				name += ".png";
//				plotHybrid(dists, metric, threshold, maxZ, new File(outputDir, name), zoom, distSpacing);
//			}
//			
//			dists = null;
//			bitSets = null;
			
			// now rotated
			double[][] rotatedData = calcRotated(rotatedSubPath, metric, rotated_width);
			
//			String name = "rp_"+metric.name()+"_"+threshStr+"_rotated.png";
//			plotDiscreteRotated(rotatedData, metric, rotatedPixelWidth, threshold, distSpacing, new File(rotatedDir, name));
//			
//			name = "rp_"+metric.name()+"_continuous_rotated.png";
//			plotContinuousRotated(rotatedData, metric, rotatedPixelWidth, maxZ, distSpacing, new File(rotatedDir, name));
			
			String name = "rp_"+metric.name()+"_hybrid_"+threshStr+"_rotated.png";
			plotHybridRotated(rotatedData, metric, rotated_pixel_width, threshold, maxZ, distSpacing, new File(outputDir, name));
			
			// now with annotations
			if (idens != null) {
				List<String> idenNames = Lists.newArrayList();
				for (RuptureIdentifier iden : idens)
					idenNames.add(iden.getName());
				CPT cpt = getHybridCPT(threshold, maxZ);
				List<List<int[]>> paths = Lists.newArrayList();
				paths.add(fullPath);
				List<double[][]> datas = Lists.newArrayList();
				datas.add(rotatedData);
				name = "rp_"+metric.name()+"_hybrid_"+threshStr+"_rotated_with_events.png";
				plotRotated(datas, null, metric, cpt, rotated_pixel_width, distSpacing,
						new File(outputDir, name), paths, idenNames);
			}
		}
	}

}
