package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Shaw_2009_ModifiedMagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagAreaRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class MagAreaScalingPlot extends AbstractPlot {
	
	private boolean slip;
	private Map<FaultStyle, DefaultXY_DataSet> scatterMap;

	private MinMaxAveTracker magTrack = new MinMaxAveTracker();
	private MinMaxAveTracker slipTrack = new MinMaxAveTracker();
	private MinMaxAveTracker areaTrack = new MinMaxAveTracker();
	
	private static final int max_scatter_points = 1000000;
	
	private static double[] csv_fractiles = { 0d, 0.025, 0.16, 0.5, 0.84, 0.975, 1d };
	private static double csv_mag_delta = 0.1;
	
	private static Map<FaultStyle, List<MagAreaRelationship>> comparisons;
	
	static {
		List<FaultStyle> styles = new ArrayList<>();
		for (FaultStyle style : FaultStyle.values())
			styles.add(style);
		styles.add(null);
		
		comparisons = new HashMap<>();
		
		for (FaultStyle style : styles) {
			List<MagAreaRelationship> rels = new ArrayList<>();
			
			rels.add(new WC1994_MagAreaRelationship());
			rels.add(new Ellsworth_B_WG02_MagAreaRel());
			rels.add(new HanksBakun2002_MagAreaRel());
			rels.add(new Shaw_2009_ModifiedMagAreaRel());
			
			double rake = Double.NaN;
			if (style != null) {
				switch (style) {
				case STRIKE_SLIP:
					rake = 0d;
					break;
				case REVERSE:
					rake = 90;
					break;
				case NORMAL:
					rake = -90d;
					break;

				default:
					break;
				}
			}
			
			for (MagAreaRelationship rel : rels)
				rel.setRake(rake);
			
			comparisons.put(style, rels);
		}
	}
	
	public MagAreaScalingPlot(boolean slip) {
		scatterMap = new HashMap<>();
		for (FaultStyle style : FaultStyle.values())
			scatterMap.put(style, new DefaultXY_DataSet());
		this.slip = slip;
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (!Doubles.isFinite(mag))
			return;
		double areaM2 = e.getArea(); // m^2
		double areaKM2 = areaM2/1e6;
		
		FaultStyle style = RSQSimUtils.calcFaultStyle(e, 15, 0.1);
		DefaultXY_DataSet scatter = scatterMap.get(style);
		
		if (slip) {
			double moment = 0d;
			for (EventRecord r : e)
				moment += r.getMoment();
			double meanSlip = FaultMomentCalc.getSlip(areaM2, moment);
			slipTrack.addValue(meanSlip);
			scatter.set(areaKM2, meanSlip);
		} else {
			scatter.set(areaKM2, mag);
		}
		
		areaTrack.addValue(areaKM2);
		magTrack.addValue(mag);
	}

	@Override
	public void finalizePlot() throws IOException {
		DefaultXY_DataSet masterScatter = new DefaultXY_DataSet();
		int width = getPlotWidth();
		int height = getPlotHeight();
		Range xRange, yRange;
		if (slip) {
			xRange = new Range(0d, areaTrack.getMax()*1.1);
			yRange = new Range(Math.floor(slipTrack.getMin()), Math.ceil(slipTrack.getMax()));
		} else {
			xRange = calcEncompassingLog10Range(areaTrack.getMin(), areaTrack.getMax());
			for (List<MagAreaRelationship> comps : comparisons.values()) {
				for (MagAreaRelationship comp : comps) {
					magTrack.addValue(comp.getMedianMag(areaTrack.getMin()));
					magTrack.addValue(comp.getMedianMag(areaTrack.getMax()));
				}
			}
			yRange = new Range(Math.floor(magTrack.getMin()), Math.ceil(magTrack.getMax()));
		}
		for (FaultStyle style : scatterMap.keySet()) {
			DefaultXY_DataSet scatter = scatterMap.get(style);
			if (scatter.size() == 0)
				continue;
			for (Point2D pt : scatter)
				masterScatter.set(pt);
			writeScatterPlots(scatter, slip, style, getCatalogName(), getOutputDir(), getOutputPrefix()+"_"+style.name(),
					xRange, yRange, width, height);
		}
		writeScatterPlots(masterScatter, slip, null, getCatalogName(), getOutputDir(), getOutputPrefix(), xRange, yRange, width, height);
	}
	
	public static void writeScatterPlots(XY_DataSet scatter, boolean meanSlip, String name, File outputDir, String prefix,
			int plotWidth, int plotHeight) throws IOException {
		writeScatterPlots(scatter, meanSlip, null, name, outputDir, prefix, null, null, plotWidth, plotHeight);
	}
	
	public static void writeScatterPlots(XY_DataSet scatter, boolean meanSlip, FaultStyle faultStyle, String name, File outputDir, String prefix,
			Range xRange, Range yRange, int plotWidth, int plotHeight) throws IOException {
		// scatter plot
		
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		XY_DataSet plotScatter = scatter;
		if (scatter.size() > max_scatter_points) {
			System.out.println("Filtering M-A scatter from "+scatter.size()+" to ~"+max_scatter_points+" points");
//			plotScatter = new DefaultXY_DataSet();
//			Random r = new Random();
//			for (int i=0; i<max_scatter_points; i++)
//				plotScatter.set(scatter.get(r.nextInt(scatter.size())));
			plotScatter = downsampleByMag(scatter, false, max_scatter_points);
			System.out.println("Filter done (random mag-dependent sample): "+plotScatter.size());
		}
		funcs.add(plotScatter);
		plotScatter.setName(name);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		System.out.println("Scatter y range: "+scatter.getMinY()+" "+scatter.getMaxY());
		
		List<EvenlyDiscretizedFunc> compFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> compChars = new ArrayList<>();
		
		Color[] compColors = { Color.RED.darker(), Color.BLUE.darker(), Color.GREEN.darker(), Color.ORANGE.darker(), Color.MAGENTA.darker() };
		
		if (!meanSlip) {
			List<MagAreaRelationship> list = comparisons.get(faultStyle);
			for (int m=0; m<list.size(); m++) {
				MagAreaRelationship ma = list.get(m);
				EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), 1000);
				
				if (ma instanceof WC1994_MagAreaRelationship)
					func.setName("W-C 1994");
				else if (ma instanceof Ellsworth_B_WG02_MagAreaRel)
					func.setName("EllsworthB");
				else if (ma instanceof HanksBakun2002_MagAreaRel)
					func.setName("H-B 2002");
				else if (ma instanceof Shaw_2009_ModifiedMagAreaRel)
					func.setName("Shaw 2009 (Mod)");
				
				for (int i=0; i<func.size(); i++)
					func.set(i, ma.getMedianMag(func.getX(i)));
				
				compFuncs.add(func);
				compChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColors[m % compColors.length]));
			}
		}
		
		funcs.addAll(compFuncs);
		chars.addAll(compChars);
		
		String title, yAxisLabel;
		if (meanSlip) {
			title = "Slip-Area Scaling";
			yAxisLabel = "Mean Slip (m)";
		} else {
			title = "Mag-Area Scaling";
			yAxisLabel = "Magnitude";
		}
		String xAxisLabel = "Area (km^2)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		double minY, maxY;
		if (meanSlip) {
			minY = Math.floor(scatter.getMinY());
			maxY = Math.ceil(scatter.getMaxY());
			if (xRange == null)
				xRange = new Range(0d, scatter.getMaxX()*1.1);
		} else {
			double minX = scatter.getMinX();
			double maxX = scatter.getMaxX();
			minY = scatter.getMinY();
			maxY = scatter.getMaxY();
			for (DiscretizedFunc func : compFuncs) {
				minX = Math.min(minX, func.getMinX());
				maxX = Math.max(maxX, func.getMaxX());
				minY = Math.min(minY, func.getMinY());
				maxY = Math.max(maxY, func.getMaxY());
			}
			minY = Math.floor(minY);
			maxY = Math.ceil(maxY);
			if (xRange == null)
				xRange = calcEncompassingLog10Range(minX, maxX);
		}
		if (yRange == null)
			yRange = new Range(minY, maxY);
		
		HeadlessGraphPanel gp = buildGraphPanel();
		gp.drawGraphPanel(plot, !meanSlip, false, xRange, yRange);
		gp.getChartPanel().setSize(plotWidth, plotHeight);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		
		int nx = 51;
		int ny = 51;
		
		double minX, maxX;
		if (meanSlip) {
			minX = xRange.getLowerBound();
			maxX = xRange.getUpperBound();
		} else {
			minX = Math.log10(xRange.getLowerBound());
			maxX = Math.log10(xRange.getUpperBound());
		}
		double gridSpacingX = (maxX - minX)/(nx-1);
		
		System.out.println("XYZ minX="+minX+", maxX="+maxX+", spacing="+gridSpacingX);
		
		double gridSpacingY = (maxY - minY)/(ny-1);
		
		// XYZ plot (2D hist)
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(nx, ny, minX, minY, gridSpacingX, gridSpacingY);
		
		for (Point2D pt : scatter) {
			int index;
			if (meanSlip)
				index = xyz.indexOf(pt.getX(), pt.getY());
			else
				index = xyz.indexOf(Math.log10(pt.getX()), pt.getY());
			if (index < 0 || index >= xyz.size())
				throw new IllegalStateException("Scatter point not in XYZ range. x: "
							+pt.getX()+" ["+xyz.getMinX()+" "+xyz.getMaxX()
						+"], y: "+pt.getY()+" ["+xyz.getMinY()+" "+xyz.getMaxY()+"]");
			xyz.set(index, xyz.get(index)+1);
		}
		// convert to density
		for (int i=0; i<xyz.size(); i++) {
			// convert to density
			Point2D pt = xyz.getPoint(i);
			double x = pt.getX();
			double binWidth;
			if (meanSlip)
				binWidth = gridSpacingX;
			else
				binWidth = Math.pow(10, x + 0.5*gridSpacingX) - Math.pow(10, x - 0.5*gridSpacingX);
			double binHeight = gridSpacingY;
			double area = binWidth * binHeight;
			xyz.set(i, xyz.get(i)*area);
		}
		xyz.scale(1d/xyz.getSumZ());
		
		// set all zero to NaN so that it will plot white
		for (int i=0; i<xyz.size(); i++) {
			if (xyz.get(i) == 0)
				xyz.set(i, Double.NaN);
		}
		xyz.log10();
		
		double minZ = Double.POSITIVE_INFINITY;
		double maxZ = Double.NEGATIVE_INFINITY;
		for (int i=0; i<xyz.size(); i++) {
			double val = xyz.get(i);
			if (!Doubles.isFinite(val))
				continue;
			if (val < minZ)
				minZ = val;
			if (val > maxZ)
				maxZ = val;
		}
		
		System.out.println("MinZ: "+minZ);
		System.out.println("MaxZ: "+maxZ);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		if ((float)minZ == (float)maxZ)
			cpt = cpt.rescale(minZ, minZ*2);
		else if (!Doubles.isFinite(minZ))
			cpt = cpt.rescale(0d, 1d);
		else
			cpt = cpt.rescale(minZ, maxZ);
		cpt.setNanColor(new Color(255, 255, 255, 0));
		
		String zAxisLabel = "Log10(Density)";
		
		if (!meanSlip)
			xAxisLabel = "Log10 "+xAxisLabel;
		XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, title, xAxisLabel, yAxisLabel, zAxisLabel);
		if (!meanSlip) {
			// add W-C
			funcs = Lists.newArrayList();
			chars = Lists.newArrayList(compChars);
			
			for (DiscretizedFunc func : compFuncs) {
				ArbitrarilyDiscretizedFunc logFunc = new ArbitrarilyDiscretizedFunc();
				for (Point2D pt : func)
					logFunc.set(Math.log10(pt.getX()), pt.getY());
				funcs.add(logFunc);
			}
			
			xyzSpec.setXYElems(funcs);
			xyzSpec.setXYChars(chars);
		}
		
		XYZGraphPanel xyzGP = buildXYZGraphPanel();
		xyzGP.drawPlot(xyzSpec, false, false, new Range(minX-0.5*gridSpacingX, maxX+0.5*gridSpacingX),
				new Range(Math.max(0d, yRange.getLowerBound()-0.5*gridSpacingY), yRange.getUpperBound()+0.5*gridSpacingY));
		// write plot
		xyzGP.getChartPanel().setSize(plotWidth, plotHeight);
		xyzGP.saveAsPNG(new File(outputDir, prefix+"_hist2D.png").getAbsolutePath());
		xyzGP.saveAsPDF(new File(outputDir, prefix+"_hist2D.pdf").getAbsolutePath());
		
		// now write CSV
		if (!meanSlip) {
			IncrementalMagFreqDist magFunc = MFDPlot.buildIncrementalFunc(scatter.getMinY(), csv_mag_delta);
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<String>();
			header.add("Magnitude");
			header.add("Mean");
			header.add("Standard Deviation");
			for (double f : csv_fractiles)
				header.add((float)f+" fractile");
			csv.addLine(header);
			
			for (int i=0; i<magFunc.size(); i++) {
				double mag = magFunc.getX(i);
				double minMag = mag - 0.5*csv_mag_delta;
				double maxMag = mag + 0.5*csv_mag_delta;
				
				List<Double> areasForBin = new ArrayList<>();
				
				for (Point2D pt : scatter)
					if (pt.getY() >= minMag && pt.getY() < maxMag)
						areasForBin.add(pt.getX());
				
				double[] areasArray = Doubles.toArray(areasForBin);
				List<String> line = new ArrayList<>();
				line.add((float)mag+"");
				line.add(StatUtils.mean(areasArray)+"");
				line.add(Math.sqrt(StatUtils.variance(areasArray))+"");
				for (double f : csv_fractiles)
					if (f == 0d)
						line.add(StatUtils.min(areasArray)+"");
					else
						line.add(StatUtils.percentile(areasArray, f*100d)+"");
				csv.addLine(line);
			}
			
			csv.writeToFile(new File(outputDir, prefix+".csv"));
		}
	}
	
	public static void plotMultiMagArea(Collection<CSVFile<String>> csvs, CSVFile<String> baselineCSV, String baselineName,
			boolean variationFractiles, boolean median, double[] fractiles, File outputDir, String prefix) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		String percentileString = null;
		if (fractiles != null) {
			List<String> fractileStrings = new ArrayList<>();
			for (int i=0; i<fractiles.length; i++)
				fractileStrings.add((float)(fractiles[i]*100d)+"");
			percentileString = Joiner.on(",").join(fractileStrings)+" %-ile";
			if (fractileStrings.size() > 1)
				percentileString += "s";
		}
		
		if (baselineCSV != null) {
			List<DiscretizedFunc> baselineFuncs = loadCSV(baselineCSV, median, fractiles, false);
			
			if (baselineName == null)
				baselineName = "";
			else
				baselineName += " ";
			if (median)
				baselineName += "Median";
			else
				baselineName += "Mean";
			
			baselineFuncs.get(0).setName(baselineName);
			funcs.add(baselineFuncs.get(0));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			
			if (fractiles != null) {
				for (int i=0; i<fractiles.length; i++) {
					DiscretizedFunc func = baselineFuncs.get(i+1);
					if (i == 0)
						func.setName(percentileString);
					funcs.add(func);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.BLACK));
				}
			}
		}
		
		boolean first = true;
		for (CSVFile<String> csv : csvs) {
			if (csv == baselineCSV)
				continue;
			
			List<DiscretizedFunc> csvFuncs = loadCSV(csv, median, fractiles, false);
			
			if (first) {
				if (baselineCSV != null) {
					csvFuncs.get(0).setName("Variations");
				} else {
					if (median)
						csvFuncs.get(0).setName("Median");
					else
						csvFuncs.get(0).setName("Mean");
				}
			}
			funcs.add(csvFuncs.get(0));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
			
			if (fractiles != null && variationFractiles) {
				for (int i=0; i<fractiles.length; i++) {
					DiscretizedFunc func = csvFuncs.get(i+1);
					if (first && i == 0)
						func.setName(percentileString);
					funcs.add(func);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(128, 128, 128, 60)));
				}
			}
			first = false;
		}
		
		String title = "Mag-Area Scaling";
		String xAxisLabel = "Area (km^2)";
		String yAxisLabel = "Magnitude";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0d;
		double minX = Double.POSITIVE_INFINITY;
		double maxX = 0d;
		for (DiscretizedFunc func : funcs) {
			minY = Math.min(minY, AbstractPlot.minNonZero(func));
			maxY = Math.max(maxY, func.getMaxY());
			minX = Math.min(minX, func.getMinX());
			for (Point2D pt : func)
				if (pt.getY() > 0)
					maxX = Math.max(maxX, pt.getX());
		}
		Range yRange = new Range(minY*0.9, maxY*1.1);
		Range xRange = new Range(minX*0.9, maxX*1.1);
		
		HeadlessGraphPanel gp = buildGraphPanel();
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		gp.drawGraphPanel(plot, true, false, xRange, yRange);
		gp.getChartPanel().setSize(650, 600);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	static List<DiscretizedFunc> loadCSV(CSVFile<String> csv, boolean median, double[] fractiles, boolean magX) {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		DiscretizedFunc centralFunc = new ArbitrarilyDiscretizedFunc();
		funcs.add(centralFunc);
		DiscretizedFunc[] fractileFuncs = null;
		int[] fractileCols = null;
		if (fractiles != null) {
			fractileFuncs = new DiscretizedFunc[fractiles.length];
			fractileCols = new int[fractiles.length];
			for (int i=0; i<fractiles.length; i++) {
				fractileFuncs[i] = new ArbitrarilyDiscretizedFunc();
				funcs.add(fractileFuncs[i]);
				for (int col=3; col<csv.getNumCols(); col++) {
					if (csv.get(0, col).startsWith((float)fractiles[i]+""))
						fractileCols[i] = col;
				}
				Preconditions.checkState(fractileCols[i] > 0);
			}
		}
		int centralCol;
		if (median) {
			centralCol = -1;
			for (int col=3; col<csv.getNumCols(); col++) {
				if (csv.get(0, col).startsWith("0.5 "))
					centralCol = col;
			}
			Preconditions.checkState(centralCol > 0);
		} else {
			centralCol = 1;
		}
		
		for (int row=1; row<csv.getNumRows(); row++) {
			double mag = Double.parseDouble(csv.get(row, 0));
			double centralVal = Double.parseDouble(csv.get(row, centralCol));
			if (Double.isNaN(centralVal))
				continue;
			if (magX)
				centralFunc.set(mag, centralVal);
			else
				centralFunc.set(centralVal, mag);
			if (fractileFuncs != null)
				for (int i=0; i<fractileFuncs.length; i++)
					if (magX)
						fractileFuncs[i].set(mag, Double.parseDouble(csv.get(row, fractileCols[i])));
					else
						fractileFuncs[i].set(Double.parseDouble(csv.get(row, fractileCols[i])), mag);
		}
		
		return funcs;
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

}
