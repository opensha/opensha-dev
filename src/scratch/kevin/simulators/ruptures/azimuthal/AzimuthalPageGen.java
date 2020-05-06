package scratch.kevin.simulators.ruptures.azimuthal;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.geometry.euclidean.twod.Line;
import org.apache.commons.math3.geometry.euclidean.twod.Segment;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

public class AzimuthalPageGen {
	
	private AzimuthalSiteConfig<Integer> config;
	private AzimuthalZipLoader loader;
	private double[] periods;
	
	private EvenlyDiscrXYZ_DataSet refXYZ;
	
	// everything stored in ln space
	private Map<Integer, EvenlyDiscrXYZ_DataSet[]> eventXYZs;
	private Map<Integer, EvenlyDiscrXYZ_DataSet[]> eventResidualXYZs;

	public AzimuthalPageGen(AzimuthalSiteConfig<Integer> config, AzimuthalZipLoader loader,
			double[] periods) throws IOException {
		this.config = config;
		this.loader = loader;
		this.periods = periods;
		
		refXYZ = config.getGC2XYZ();
		
		eventXYZs = new HashMap<>();
		eventResidualXYZs = new HashMap<>();
		
		List<Integer> rups = config.getRuptures();
		
		System.out.println("Building shakemaps...");
		for (int i=0; i<rups.size(); i++) {
			int id = rups.get(i);
			EvenlyDiscrXYZ_DataSet[] xyzs = new EvenlyDiscrXYZ_DataSet[periods.length];
			
			for (int p=0; p<periods.length; p++)
				xyzs[p] = refXYZ.copy();
			
			for (int siteIndex=0; siteIndex<refXYZ.size(); siteIndex++) {
				DiscretizedFunc rd50 = loader.readRotD50(siteIndex, id);
				for (int p=0; p<periods.length; p++)
					xyzs[p].set(siteIndex, Math.log(rd50.getInterpolatedY(periods[p])));
			}
			
			eventXYZs.put(id, xyzs);
			
			eventResidualXYZs.put(id, calcResidualXYZ(id, xyzs));
		}
		System.out.println("done");
	}
	
	private EvenlyDiscrXYZ_DataSet resample(EvenlyDiscrXYZ_DataSet xyz, int numPer) {
		EvenlyDiscrXYZ_DataSet interp = new EvenlyDiscrXYZ_DataSet((xyz.getNumX()-1)*numPer+1, (xyz.getNumY()-1)*numPer+1,
				xyz.getMinX(), xyz.getMinY(), xyz.getGridSpacingX()/(double)numPer);
//		System.out.println(xyz.getNumX()+" x "+xyz.getNumY());
//		System.out.println(xyz.getMinX()+" "+xyz.getMaxX()+" "+xyz.getMinY()+" "+xyz.getMaxY());
//		System.out.println(interp.getNumX()+" x "+interp.getNumY());
//		System.out.println(interp.getMinX()+" "+interp.getMaxX()+" "+interp.getMinY()+" "+interp.getMaxY());
		for (int i=0; i<interp.size(); i++) {
			Point2D pt = interp.getPoint(i);
			interp.set(i, xyz.bilinearInterpolation(pt.getX(), pt.getY()));
		}
		
		return interp;
	}
	
	private EvenlyDiscrXYZ_DataSet[] calcResidualXYZ(Integer rupture, EvenlyDiscrXYZ_DataSet[] xyzs) {
		EvenlyDiscrXYZ_DataSet[] ret = new EvenlyDiscrXYZ_DataSet[xyzs.length];
		
		Point2D[] rect = config.getRectangle(rupture);
		
		// figure out the maximum distance we can go in all directions
		double maxDist = refXYZ.getMaxY() - rect[1].getY();
		maxDist = Math.min(maxDist, -refXYZ.getMinY());
		maxDist = Math.min(maxDist, refXYZ.getMaxX() - rect[2].getX());
//		maxDist -= 0.5*refXYZ.getGridSpacingX();
		
		int numDist = (int)(2d * maxDist / refXYZ.getGridSpacingX());
		
		EvenlyDiscrXYZ_DataSet dists = null;
		
		for (int p=0; p<xyzs.length; p++) {
			// resample it to be really high resolution
			EvenlyDiscrXYZ_DataSet interp = resample(xyzs[p], 10);
			if (dists == null) {
				dists = interp.copy();
				for (int i=0; i<dists.size(); i++)
					dists.set(i, distToRect(dists.getPoint(i), rect));
			}
			
			EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(0d, maxDist, numDist+1);
			// now shift by half a bin
			distFunc = new EvenlyDiscretizedFunc(0.5*distFunc.getDelta(), numDist, distFunc.getDelta());
			List<List<Double>> binnedGMs = new ArrayList<>();
			for (int i=0; i<numDist; i++)
				binnedGMs.add(new ArrayList<>());
			
			for (int i=0; i<interp.size(); i++) {
				double dist = dists.get(i);
				if (dist > maxDist)
					continue;
				int distIndex = distFunc.getClosestXIndex(dist);
				binnedGMs.get(distIndex).add(interp.get(i));
			}
//			if (p == 0 && rupture == config.getRuptures().get(0))
//				for (int i=0; i<numDist; i++)
//					System.out.println(i+". x="+(float)distFunc.getX(i)+", "+binnedGMs.get(i).size()+" values");
			
			// now build dist func
			for (int i=0; i<numDist; i++) {
				List<Double> gms = binnedGMs.get(i);
				Preconditions.checkState(gms.size() >= 10, "Have %s ground motions for bin at distance %s. MaxDist=%s",
						gms.size(), distFunc.getX(i), maxDist);
				double median = DataUtils.median(Doubles.toArray(gms));
				distFunc.set(i, median);
			}
			
			// now calculate residuals
			ret[p] = refXYZ.copy();
			for (int i=0; i<ret[p].size(); i++) {
				double dist = distToRect(ret[p].getPoint(i), rect);
				if (dist >= maxDist) {
					ret[p].set(i, Double.NaN);
				} else {
					double median;
					if (dist > distFunc.getMaxX())
						median = distFunc.getY(distFunc.size()-1);
					else if (dist < distFunc.getMinX())
						median = distFunc.getY(0);
					else
						median = distFunc.getInterpolatedY(dist);
					ret[p].set(i, xyzs[p].get(i) - median);
				}
			}
		}
		
		return ret;
	}
	
	private double distToRect(Point2D pt, Point2D[] rect) {
		double x = pt.getX();
		double y = pt.getY();
		
		double minX = rect[0].getX();
		double maxX = rect[2].getX();
		double minY = rect[0].getY();
		double maxY = rect[2].getY();
		
		// check if it contains it
		if ((float)x >= (float)minX && (float)x <= (float)maxX && (float)y >= (float)minY && (float)y <= (float)maxY)
			return 0d;
		double minDist = Double.POSITIVE_INFINITY;
		
		// TODO be better
		int numEach = 10;
		for (int i=0; i<4; i++) {
			Point2D p1 = rect[i];
			Point2D p2 = rect[(i+1) % 4];
			for (int j=0; j<numEach+1; j++) {
				double fract = (double)j/(double)numEach;
				double myX = fract*p1.getX() + (1d-fract)*p2.getX();
				double myY = fract*p1.getY() + (1d-fract)*p2.getY();
				minDist = Math.min(minDist, pt.distance(myX, myY));
			}
		}
		
//		Vector2D v = new Vector2D(pt.getX(), pt.getY());
//		
//		for (int i=0; i<4; i++) {
//			Point2D p1 = rect[i];
//			Point2D p2 = rect[(i+1) % 4];
////			minDist = Double.min(minDist, pt.get)
//			Vector2D v1 = new Vector2D(p1.getX(), p1.getY());
//			Vector2D v2 = new Vector2D(p2.getX(), p2.getY());
//			Line line = new Line(v1, v2, 0.1d);
//			Segment seg = new Segment(v1, v2, line);
//			minDist = Math.min(minDist, seg.distance(v));
//		}
		Preconditions.checkState(Double.isFinite(minDist), "Bad distance calculation: %s for Point: %s", minDist, pt);
		return minDist;
	}

	private static XY_DataSet dippingXY(double length, double width) {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		xy.set(0d, 0d);
		xy.set(0d, length);
		xy.set(width, length);
		xy.set(width, 0d);
		xy.set(0d, 0d);
		return xy;
	}
	
	private static XY_DataSet ssXY(double length, double whiskerWidth) {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		double tickLen = 0.5*whiskerWidth;
		xy.set(-tickLen, 0d);
		xy.set(tickLen, 0d);
		xy.set(0d, 0d);
		xy.set(0d, length);
		xy.set(-tickLen, length);
		xy.set(tickLen, length);
		return xy;
	}
	
	private void calcXYZs(List<Integer> rups, EvenlyDiscrXYZ_DataSet[] medXYZs,
			EvenlyDiscrXYZ_DataSet[] stdXYZs, EvenlyDiscrXYZ_DataSet[] resXYZs) {
		for (int p=0; p<periods.length; p++) {
			medXYZs[p] = refXYZ.copy();
			stdXYZs[p] = refXYZ.copy();
			resXYZs[p] = refXYZ.copy();
		}
		
		for (int s=0; s<refXYZ.size(); s++) {
			for (int p=0; p<periods.length; p++) {
				double[] vals = new double[rups.size()];
				for (int i=0; i<vals.length; i++)
					vals[i] = eventXYZs.get(rups.get(i))[p].get(s);
				medXYZs[p].set(s, DataUtils.median(vals));
				stdXYZs[p].set(s, Math.sqrt(StatUtils.variance(vals)));
				// now residuals
				vals = new double[rups.size()];
				boolean nan = false;
				for (int i=0; i<vals.length; i++) {
					vals[i] = eventResidualXYZs.get(rups.get(i))[p].get(s);
					nan = nan || Double.isNaN(vals[i]);
				}
				if (nan)
					resXYZs[p].set(s, Double.NaN);
				else
					resXYZs[p].set(s, DataUtils.median(vals));
			}
		}
	}
	
	
	private File[] plotXYZs(File outputDir, String prefix, EvenlyDiscrXYZ_DataSet[] xyzs,
			List<Integer> rups, CPT cpt, String zLabel, boolean rescaleCPT) throws IOException {
		File[] ret = new File[periods.length];
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(20);
		plotPrefs.setAxisLabelFontSize(22);
		plotPrefs.setPlotLabelFontSize(24);
		plotPrefs.setLegendFontSize(20);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		double spacing = refXYZ.getGridSpacingX();
		double minX = refXYZ.getMinX() - 0.5*spacing;
		double maxX = refXYZ.getMaxX() + 0.5*spacing;
		double minY = refXYZ.getMinY() - 0.5*spacing;
		double maxY = refXYZ.getMaxY() + 0.5*spacing;
		
		double spanX = maxX - minX;
		double spanY = maxY - minY;
		
		double[] lengths = new double[rups.size()];
		double[] horzWidths = new double[rups.size()];
		
		for (int i=0; i<rups.size(); i++) {
			int id = rups.get(i);
			Point2D[] rect = config.getRectangle(id);
			lengths[i] = rect[1].getY();
			horzWidths[i] = rect[2].getX();
		}
		double minLen = StatUtils.min(lengths);
		double maxLen = StatUtils.max(lengths);
		double medLen = DataUtils.median(lengths);
		double minHorzWid = StatUtils.min(horzWidths);
		double maxHorzWid = StatUtils.max(horzWidths);
		double medHorzWid = DataUtils.median(horzWidths);

		Color medianColor = Color.CYAN.darker();
		Color extremaColor = new Color(medianColor.getRed(),
				medianColor.getGreen(), medianColor.getBlue(), 127);
		
		PlotCurveCharacterstics extremaChar =
				new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, extremaColor);
		PlotCurveCharacterstics medianChar =
				new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, medianColor);
		
		List<XY_DataSet> traceFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> traceChars = new ArrayList<>();
		if (medHorzWid > 0.1*medLen) {
			// dipping
			if ((float)minLen < (float)maxLen) {
				traceFuncs.add(dippingXY(minLen, minHorzWid));
				traceChars.add(extremaChar);
				traceFuncs.add(dippingXY(maxLen, maxHorzWid));
				traceChars.add(extremaChar);
			}
			traceFuncs.add(dippingXY(medLen, medHorzWid));
			traceChars.add(medianChar);
		} else {
			double whiskerWidth = 0.05*spanX;
			if ((float)minLen < (float)maxLen) {
				traceFuncs.add(ssXY(minLen, whiskerWidth));
				traceChars.add(extremaChar);
				traceFuncs.add(ssXY(maxLen, whiskerWidth));
				traceChars.add(extremaChar);
			}
			traceFuncs.add(ssXY(medLen, whiskerWidth));
			traceChars.add(medianChar);
		}
		
		int width = 800;
		int height = (int)(width*spanY/spanX);
		
		for (int p=0; p<periods.length; p++) {
			CPT myCPT;
			if (rescaleCPT) {
//				double max = Math.ceil(xyzs[p].getMaxZ());
//				double min = Math.floor(xyzs[p].getMinZ());
				double max = xyzs[p].getMaxZ();
				double min = xyzs[p].getMinZ();
				myCPT = cpt.rescale(min, max);
			} else {
				myCPT = cpt;
			}
			XYZPlotSpec spec = new XYZPlotSpec(xyzs[p], myCPT, " ", "RX (km)", "RY (km)",
					optionalDigitDF.format(periods[p])+"s "+zLabel);
			spec.setXYElems(traceFuncs);
			spec.setXYChars(traceChars);
			spec.setCPTPosition(RectangleEdge.BOTTOM);
			
			XYZGraphPanel gp = new XYZGraphPanel(plotPrefs);
			
			gp.drawPlot(spec, false, false, new Range(minX, maxX), new Range(minY, maxY));
//			gp.getYAxis().setStandardTickUnits(tus);
//			gp.getXAxis().setStandardTickUnits(tus);
			gp.getChartPanel().setSize(width, height);
			
			String myPrefix = prefix+"_"+optionalDigitDF.format(periods[p])+"s";
			
			ret[p] = new File(outputDir, myPrefix+".png");
			gp.saveAsPNG(ret[p].getAbsolutePath());
		}
		return ret;
	}
	
	private List<Integer> getRupsByFractDAS(double min, double max) {
		List<Integer> rups = new ArrayList<>();
		
		for (int rup : config.getRuptures()) {
			double len = config.getLength(rup);
			double das = config.getHypocenterDAS(rup);
			double fract = das/len;
			if (fract >= min && fract <= max)
				rups.add(rup);
		}
		
		return rups;
	}
	
	private List<Integer> getRupsByFractDDW(double min, double max) {
		List<Integer> rups = new ArrayList<>();
		
		for (int rup : config.getRuptures()) {
			double width = config.getWidth(rup);
			double ddw = config.getHypocenterDDW(rup);
			double fract = ddw/width;
			if (fract >= min && fract <= max)
				rups.add(rup);
		}
		
		return rups;
	}
	
	private void plotHypocenterDistribution(File resourcesDir, String prefix, List<Integer> rups)
			throws IOException {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		List<Double> lengths = new ArrayList<>();
		List<Double> widths = new ArrayList<>();
		for (Integer rupture : rups) {
			double hypoDAS = config.getHypocenterDAS(rupture);
			double hypoDDW = config.getHypocenterDDW(rupture);
			
			double len = config.getLength(rupture);
			hypoDAS /= len;
			double width = config.getWidth(rupture);
			hypoDDW /= width;
			
			lengths.add(len);
			widths.add(width);
			
			xy.set(hypoDAS, hypoDDW);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet outline = new DefaultXY_DataSet();
		outline.set(0d, 0d);
		outline.set(1d, 0d);
		outline.set(1d, 1d);
		outline.set(0d, 1d);
		outline.set(0d, 0d);
		
		funcs.add(outline);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Hypocenter Distribution",
				"Fractional Distance Along Strike", "Fractional Distance Down Dip");
		
		double length = DataUtils.median(Doubles.toArray(lengths));
		double width = DataUtils.median(Doubles.toArray(widths));
		
//		System.out.println("Len: "+length);
//		System.out.println("Width: "+width);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		Range xRange = new Range(-0.05, 1.05);
		Range yRange = new Range(-0.05, 1.05);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.setyAxisInverted(true);
		
		File file = new File(resourcesDir, prefix);
		int imgHeight = 400;
		int imgWidth = (int)((double)imgHeight*length/width);
		gp.getChartPanel().setSize(imgWidth, imgHeight);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
//		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private List<String> buildRupInfoLines(File resourcesDir, String prefix, List<Integer> rups, String hypocenterDesc) throws IOException {
		List<String> lines = new ArrayList<>();
		
		String line = "These plots show ground motions across "+rups.size()+" ruptures";
		if (hypocenterDesc != null)
			line += " with "+hypocenterDesc;
		line += ".";
		List<Double> lengths = new ArrayList<>();
		List<Double> widths = new ArrayList<>();
		for (Integer rup : rups) {
			lengths.add(config.getLength(rup));
			widths.add(config.getWidth(rup));
		}
		double[] lenArray = Doubles.toArray(lengths);
		double[] widArray = Doubles.toArray(widths);
		
		double minLen = StatUtils.min(lenArray);
		double maxLen = StatUtils.max(lenArray);
		double medLen = DataUtils.median(lenArray);
		
		double minWid = StatUtils.min(widArray);
		double maxWid = StatUtils.max(widArray);
		double medWid = DataUtils.median(widArray);
		
		if ((float)minLen < (float)maxLen)
			line += " Each surface is sized slightely differently. The lengh range is ["
					+(float)minLen+","+(float)maxLen+"] with a median length of "+medLen+" (km). The width range is ["
					+(float)minWid+","+(float)maxWid+"] with a median width of "+medWid+" (km).";
		else
			line += " Each surface is sized identically, with length="+(float)medLen+" (km) and width="+(float)medWid+" (km).";
		line += " The hypocenter distribution is shown below:";
		
		lines.add(line);
		
		plotHypocenterDistribution(resourcesDir, prefix, rups);
		lines.add("");
		lines.add("![hypo dist]("+resourcesDir.getName()+"/"+prefix+".png)");
		
		if (hypocenterDesc == null) {
			// add general metadata
			lines.add("");
			line = "Each column represents a defferent spectral period. The first row gives the log median "
					+ "ground motion (an average shakemap). The second row shows residuals at each point relative "
					+ "to the median ground motion at that distance (Rjb). Residuals are calculated individually for each "
					+ "ruprture, and then averaged spatially. The bottom row shows the standard deviation of log "
					+ "ground motions. ";
			if ((float)minLen < (float)maxLen) {
				if (config.getScenario().getDip() < 90d) {
					 line += "The meidan surface outline is drawn with a dark cyan line, and the extents of all surfaces "
					 		+ "light cyan lines.";
				} else {
					line += "The meidan length is annotated with a dark cyan line and ticks, and the extents of all surfaces "
					 		+ "light cyan lines and ticks.";
				}
			} else {
				if (config.getScenario().getDip() < 90d) {
					line += "The surface outline is drawn with a dark cyan line.";
				} else {
					line += "The surface length is annotated with a dark cyan line and ticks.";
				}
			}
			lines.add(line);
		}
		
		return lines;
	}
	
	public void generatePage(File outputDir, String modelName, List<String> headerLines)
			throws IOException {
		List<String> lines = new ArrayList<>();
		
		Scenario scenario = config.getScenario();
		lines.add("# "+modelName+", "+scenario.getName()+" Spatial Distributions");
		lines.add("");
		
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			lines.add("");
		}

		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Full Spatial Distributions");
		lines.add(topLink); lines.add("");
		CPT gmCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(-4, 0);
		CPT residualCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1, 1);
		residualCPT.setNanColor(Color.LIGHT_GRAY);
		CPT stdCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0d, 1d);
		
		lines.addAll(buildRupInfoLines(resourcesDir, "full_hypos", config.getRuptures(), null));
		lines.add("");
		
		TableBuilder table = buildPlotTable(resourcesDir, gmCPT, stdCPT, residualCPT, config.getRuptures(), "full");
		
		lines.addAll(table.build());
		lines.add("");
		
		boolean[] downDips;
		if (config.getScenario().getDip() < 90d)
			downDips = new boolean[] { false, true };
		else
			downDips = new boolean[] { false };
		
		for (boolean downDip : downDips) {
			if (downDips.length > 1) {
				if (downDip)
					lines.add("## Down-Dip Hypocenter Spatial Distributions");
				else
					lines.add("## Along-Strike Hypocenter Spatial Distributions");
				lines.add(topLink); lines.add("");
			}
			for (int i=0; i<3; i++) {
				double min = (double)i/3d;
				double max = (i+1d)/3d;
				String name, prefix, desc;
				List<Integer> rups;
				if (downDip) {
					if (i == 0) {
						name = "Top Third Down-Dip";
						prefix = "top_third_dip_hypos";
						desc = "hypocenters in the top third of the rupture down-dip (left third in map view)";
					} else if (i == 1) {
						name = "Center Third Down-Dip";
						prefix = "center_third_dip_hypos";
						desc = "hypocenters in the center third of the rupture down-dip";
					} else {
						name = "Bottom Third Down-Dip";
						prefix = "bottom_third_dip_hypos";
						desc = "hypocenters in the last third of the rupture down-tip (right third in map view)";
					}
					rups = getRupsByFractDDW(min, max);
				} else {
					if (i == 0) {
						name = "First Third Along-Strike";
						prefix = "first_third_strike_hypos";
						desc = "hypocenters in the first third of the rupture along-strike (bottom third in map view)";
					} else if (i == 1) {
						name = "Center Third Along-Strike";
						prefix = "center_third_strike_hypos";
						desc = "hypocenters in the center third of the rupture along-strike";
					} else {
						name = "Last Third Along-Strike";
						prefix = "last_third_strike_hypos";
						desc = "hypocenters in the last third of the rupture along-strike (top third in map view)";
					}
					rups = getRupsByFractDAS(min, max);
				}
				if (rups.size() < 2) {
					System.out.println("Skipping "+name+", no matching rups");
					continue;
				}
				System.out.println("Found "+rups.size()+" rups for "+name);
				
				if (downDips.length > 1)
					lines.add("### "+name+" Hypocenter Spatial Distributions");
				else
					lines.add("## "+name+" Hypocenter Spatial Distributions");
				lines.add(topLink); lines.add("");
				
				lines.addAll(buildRupInfoLines(resourcesDir, prefix, rups, desc));
				lines.add("");
				
				table = buildPlotTable(resourcesDir, gmCPT, stdCPT, residualCPT, rups, prefix);
				
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}

	TableBuilder buildPlotTable(File resourcesDir, CPT gmCPT, CPT stdCPT, CPT residualCPT,
			List<Integer> rups, String prefix)
			throws IOException {
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("");
		for (double period : periods)
			table.addColumn(optionalDigitDF.format(period)+" s");
		table.finalizeLine();
		EvenlyDiscrXYZ_DataSet[] medXYZs = new EvenlyDiscrXYZ_DataSet[periods.length];
		EvenlyDiscrXYZ_DataSet[] resXYZs = new EvenlyDiscrXYZ_DataSet[periods.length];
		EvenlyDiscrXYZ_DataSet[] stdXYZs = new EvenlyDiscrXYZ_DataSet[periods.length];
		calcXYZs(rups, medXYZs, stdXYZs, resXYZs);
		File[] meanPlots = plotXYZs(resourcesDir, prefix+"_mean", medXYZs, rups,
				gmCPT, "Ln(Median) RD50", true);
		File[] stdPlots = plotXYZs(resourcesDir, prefix+"_std_dev", stdXYZs, rups,
				stdCPT, "Standard Deviation", false);
		File[] resPlots = plotXYZs(resourcesDir, prefix+"_residual", resXYZs, rups,
				residualCPT, "Residual", false);
		table.initNewLine();
		table.addColumn("**Ln(Median)**");
		for (int p=0; p<periods.length; p++)
			table.addColumn("![Plot](resources/"+meanPlots[p].getName()+")");
		table.finalizeLine();
		table.initNewLine();
		table.addColumn("**Residuals**");
		for (int p=0; p<periods.length; p++)
			table.addColumn("![Plot](resources/"+resPlots[p].getName()+")");
		table.finalizeLine();
		table.initNewLine();
		table.addColumn("**Std. Dev.**");
		for (int p=0; p<periods.length; p++)
			table.addColumn("![Plot](resources/"+stdPlots[p].getName()+")");
		table.finalizeLine();
		return table;
	}
	
	protected static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");

	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException, DocumentException {
		File gitDir = new File("/home/kevin/git/rsqsim-analysis");
		File mainCatalogsDir = new File(gitDir, "catalogs");
		File mainBBPDir = new File("/data/kevin/bbp/parallel");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_4983.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_4860_10X.instance();
//		RSQSimCatalog catalog = null;
		
		double[] periods = {2d, 3d, 5d, 7.5d, 10d};
		
		File[] subDirs = mainBBPDir.listFiles();
		Arrays.sort(subDirs, new FileNameComparator());
		
		Map<Scenario, File> bbpDirs = new HashMap<>();
		
		for (File subDir : subDirs) {
			String name = subDir.getName();
			if (!name.contains("-azimuthal"))
				continue;
			if (catalog == null && !name.contains("-gp"))
				continue;
			if (catalog != null && !name.contains(catalog.getCatalogDir().getName()))
				continue;
			if (!(new File(subDir, "results_rotD.zip").exists()))
				continue;
			for (Scenario scenario : Scenario.values()) {
				File rupsFile = new File(subDir, scenario.getPrefix()+"_rups.csv");
				if (!rupsFile.exists())
					continue;
				bbpDirs.put(scenario, subDir);
			}
		}
		Preconditions.checkState(!bbpDirs.isEmpty());
		
		File catalogOutputDir = null;
		if (catalog != null) {
			catalogOutputDir = new File(mainCatalogsDir, catalog.getCatalogDir().getName());
			Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		}
		
		for (Scenario scenario : Scenario.values()) {
			if (!bbpDirs.containsKey(scenario))
				 continue;
			File bbpDir = bbpDirs.get(scenario);
			
			System.out.println("Doing "+scenario.getName());
			System.out.println("BBP dir: "+bbpDir.getAbsolutePath());
			
			VelocityModel vm = RSQSimBBP_Config.detectVM(bbpDir);
			
			File scenarioParentDir;
			String modelName;
			
			if (catalog == null) {
				scenarioParentDir = new File(gitDir, "gp_comparisons/"+bbpDir.getName());
				Preconditions.checkState(scenarioParentDir.exists() || scenarioParentDir.mkdir());
				modelName = "Graves & Pitarka (2016)";
			} else {
				File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
				Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
				scenarioParentDir = vmDir;
				modelName = catalog.getName();
			}
			
			File rupsFile = new File(bbpDir, scenario.getPrefix()+"_rups.csv");
			File locsFile = new File(bbpDir, scenario.getPrefix()+"_locs.csv");
			
			AzimuthalSiteConfig<Integer> config = AzimuthalSiteConfig.loadGeneric(scenario, locsFile, rupsFile);
			
			AzimuthalZipLoader loader = new AzimuthalZipLoader(new File(bbpDir, "results_rotD.zip"),
					config.getGC2XYZ().size(), scenario.getPrefix());
			
			AzimuthalPageGen pageGen = new AzimuthalPageGen(config, loader, periods);
			
			List<String> headerLines = new ArrayList<>();
			
			if (catalog != null) {
				headerLines.add("[RSQSim Catalog Details](../#"
						+MarkdownUtils.getAnchorName(catalog.getName())+")");
			}
			
			File scenarioDir = new File(scenarioParentDir, "azimuthal_"+scenario.getPrefix());
			
			pageGen.generatePage(scenarioDir, modelName, headerLines);
		}
		
		if (catalog != null) {
			catalog.writeMarkdownSummary(catalogOutputDir, true, false);
			RSQSimCatalog.writeCatalogsIndex(mainCatalogsDir);
		}
	}

}
