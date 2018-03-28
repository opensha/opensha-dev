package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
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
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.WarningParameter;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.MultiIMR_Averaged_AttenRel;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class ResidualScatterPlot {
	
	private int plotWidth = 900;
	private boolean writePDF = true;
	private int maxScatterPoints = 100000;
	private int numRegressionBins = 10;
	
	private XY_DataSet xy;
	private String xAxisLabel;
	private boolean logX;
	private String residualLabel;
	private String title;
	
	private Range xRange;
	private Range yRange;
	private EvenlyDiscretizedFunc xHistBins;
	
	double residualMean;
	double regressionIntercept;
	double regressionLinearSlope;
	private DiscretizedFunc regressionLinear;
	private EvenlyDiscretizedFunc regressionBins;
	private DiscretizedFunc[] regressionBinned;
	
	private double dataSigma;
	private double[] dataSigmaBinned;
	
	private EvenlyDiscretizedFunc gmpeSigma;
	
	public ResidualScatterPlot(XY_DataSet xy, String xAxisLabel, boolean logX, double gridSpacingX, String residualLabel, String title) {
		this(xy, xAxisLabel, logX, gridSpacingX, residualLabel, title, null, null);
	}
	
	public ResidualScatterPlot(XY_DataSet xy, String xAxisLabel, boolean logX, double gridSpacingX, String residualLabel, String title,
			ScalarIMR gmpeForSigma, String gmpeParameterName) {
		this.xy = xy;
		this.xAxisLabel = xAxisLabel;
		this.logX = logX;
		this.residualLabel = residualLabel;
		this.title = title;
		
		residualMean = 0d;
		for (Point2D pt : xy)
			residualMean += pt.getY();
		residualMean /= xy.size();
		
		yRange = getYRange(xy);
		
		double minX, maxX;
		if (logX) {
			minX = Math.log10(xy.getMinX());
			maxX = Math.log10(xy.getMaxX());
		} else {
			minX = xy.getMinX();
			maxX = xy.getMaxX();
		}
		int nx;
		if (Double.isNaN(gridSpacingX) || gridSpacingX <= 0) {
			nx = 51;
			gridSpacingX = (maxX - minX)/(nx-1);
		} else {
			double halfDelta = 0.5*gridSpacingX;
			double numBinsAwayFromZero = Math.floor(minX / gridSpacingX);
			double myMinX = numBinsAwayFromZero * gridSpacingX;
			// handle edge cases
			if (minX < myMinX-halfDelta)
				myMinX -= gridSpacingX;
			else if (minX > myMinX+halfDelta)
				myMinX += gridSpacingX;
//			Preconditions.checkState(minX <= myMinX + halfDelta && minX >= myMinX - halfDelta);
			nx = (int)Math.ceil((maxX - myMinX)/gridSpacingX + 1);
			minX = myMinX;
			System.out.println("Minimum bin edge adjusted to: "+minX);
		}
		
		while (xHistBins == null || xHistBins.getX(xHistBins.size()-1) - 0.5*gridSpacingX > maxX) {
			if (xHistBins == null) {
				if (nx > 1)
					nx--;
				else
					break;
			}
			xHistBins = new EvenlyDiscretizedFunc(minX + 0.5*gridSpacingX, nx, gridSpacingX);
		}
		maxX = xHistBins.getMaxX() + 0.5*gridSpacingX;
		
		if (logX)
			// xRange always in linear space
			xRange = new Range(Math.pow(10, minX), Math.pow(10, maxX));
		else
			xRange = new Range(minX, maxX);
		
		SimpleRegression linearRegression = getRegression(xy, logX, xRange);
		regressionIntercept = linearRegression.getIntercept();
		regressionLinearSlope = linearRegression.getSlope();
		System.out.println("Linear regression yInt: "+linearRegression.getIntercept());
		System.out.println("Linear regression slope: "+linearRegression.getSlope());
		regressionLinear = getRegressionFit(linearRegression, xy, logX, xRange);
		dataSigma = calcSigma(regressionLinear, xy);
		
		regressionBinned = new ArbitrarilyDiscretizedFunc[numRegressionBins];
		dataSigmaBinned = new double[numRegressionBins];
		double deltaX = (maxX - minX) / (numRegressionBins - 1);
		regressionBins = new EvenlyDiscretizedFunc(minX+0.5*deltaX, numRegressionBins, deltaX);
		System.out.println("Regression bins:\n\t"+Joiner.on("\t").join(asFloatIt(regressionBins.getXValuesIterator())));
		for (int i=0; i<numRegressionBins; i++) {
			double binCenter = regressionBins.getX(i);
			Range regressRange = new Range(binCenter - 0.5*deltaX, binCenter + 0.5*deltaX);
			if (logX)
				regressRange = new Range(Math.pow(10, regressRange.getLowerBound()), Math.pow(10, regressRange.getUpperBound()));
			regressionBinned[i] = getRegressionFit(xy, logX, regressRange);
			if (regressionBinned[i] != null)
				dataSigmaBinned[i] = calcSigma(regressionBinned[i], xy);
		}
		
		if (gmpeForSigma != null && gmpeParameterName != null) {
			gmpeSigma = new EvenlyDiscretizedFunc(regressionBins.getMinX(), regressionBins.getMaxX(), numRegressionBins);
			for (int i=0; i<numRegressionBins; i++) {
				double x = gmpeSigma.getX(i);
				if (logX)
					x = Math.pow(10, x);
				if (gmpeForSigma instanceof MultiIMR_Averaged_AttenRel) {
					((MultiIMR_Averaged_AttenRel)gmpeForSigma).setParameterInIMRs(gmpeParameterName, x);
				} else {
					Parameter<Double> gmpeParameter = gmpeForSigma.getParameter(gmpeParameterName);
					if (gmpeParameter instanceof WarningParameter)
						((WarningParameter<Double>)gmpeParameter).setValueIgnoreWarning(x);
					else
						gmpeParameter.setValue(x);
				}
				gmpeSigma.set(i, gmpeForSigma.getStdDev());
			}
			System.out.println(gmpeForSigma.getShortName()+" Sigma function:\n\t"
					+Joiner.on("\t").join(asFloatIt(gmpeSigma.getYValuesIterator())));
		}
	}
	
	private Iterator<Float> asFloatIt(final Iterator<Double> it) {
		return new Iterator<Float>() {

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Float next() {
				return it.next().floatValue();
			}
		};
	}
	
	public void plotScatter(File outputDir, String prefix) throws IOException {
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		XY_DataSet xy = this.xy;
		boolean isFiltered = false;
		if (xy.size() > maxScatterPoints) {
			isFiltered = true;
			System.out.println("Filtering scatter points from "+xy.size()+" to ~"+maxScatterPoints);
			// use fixed seed for reproducibility of downsampled plots
			Random r = new Random(xy.size());
			double rand = (double)maxScatterPoints/(double)xy.size();
			XY_DataSet filtered = new DefaultXY_DataSet();
			for (Point2D pt : xy)
				if (r.nextDouble() < rand)
					filtered.set(pt);
			System.out.println("\tNew size: "+filtered.size());
			xy = filtered;
		}
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, new Color(255, 80, 80)));
		
		DefaultXY_DataSet zeroLine = new DefaultXY_DataSet();
		zeroLine.set(xRange.getLowerBound(), 0d);
		zeroLine.set(xRange.getCentralValue(), 0d);
		zeroLine.set(xRange.getUpperBound(), 0d);
		funcs.add(zeroLine);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.DARK_GRAY));
		
		addResidualFuncs(funcs, chars, true);
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, residualLabel);
		spec.setLegendVisible(true);
		List<XYTextAnnotation> anns = buildAnnotations(false);
		if (isFiltered) {
			// add filtered note
			Font font = new Font(Font.SANS_SERIF, Font.BOLD, 20);
			
			double percent = 100d*xy.size()/this.xy.size();
			
			XYTextAnnotation ann = new XYTextAnnotation(
					" ("+maxScatterPoints+"/"+this.xy.size()+" points shown, "+annDF.format(percent)+" %)",
					xRange.getLowerBound(), 0.95*yRange.getLowerBound());
			ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
			ann.setFont(font);
			anns.add(ann);
		}
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setLegendFontSize(18);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, logX, false, xRange, yRange);
		gp.getYAxis().setStandardTickUnits(getYTickUnits());
		gp.getChartPanel().setSize(plotWidth, plotWidth-100);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		if (writePDF)
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	public void plot2DHist(File outputDir, String prefix) throws IOException {
		int nx = xHistBins.size();
		double gridSpacingX = xHistBins.getDelta();
		int ny = nx;
		double minY = yRange.getLowerBound();
		double maxY = yRange.getUpperBound();
		double gridSpacingY = (maxY - minY)/(ny-1);
		
		// XYZ plot (2D hist)
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(nx, ny, xHistBins.getMinX(), minY+0.5*gridSpacingY, gridSpacingX, gridSpacingY);
		
		for (Point2D pt : xy) {
			if (!xRange.contains(pt.getX()) || !yRange.contains(pt.getY()))
				continue;
			int index;
			if (logX)
				index = xyz.indexOf(Math.log10(pt.getX()), pt.getY());
			else
				index = xyz.indexOf(pt.getX(), pt.getY());
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
			if (logX)
				binWidth = Math.pow(10, x + 0.5*gridSpacingX) - Math.pow(10, x - 0.5*gridSpacingX);
			else
				binWidth = gridSpacingX;
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
		cpt.setNanColor(Color.WHITE);
		
		String zAxisLabel = "Log10(Density)";
		
		if (logX)
			xAxisLabel = "Log10 "+xAxisLabel;
		XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, title, xAxisLabel, residualLabel, zAxisLabel);
		ArrayList<XY_DataSet> funcs = new ArrayList<>();
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = this.xRange;
		if (logX)
			xRange = new Range(Math.log10(xRange.getLowerBound()), Math.log10(xRange.getUpperBound()));
		
		DefaultXY_DataSet zeroLine = new DefaultXY_DataSet();
		zeroLine.set(xRange.getLowerBound(), 0d);
		zeroLine.set(xRange.getCentralValue(), 0d);
		zeroLine.set(xRange.getUpperBound(), 0d);
		funcs.add(0, zeroLine);
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.DARK_GRAY));
		
		addResidualFuncs(funcs, chars, false);
		
		// add custom y axis guides
		for (int i=(int)Math.round(yRange.getLowerBound()); i<=(int)Math.round(yRange.getUpperBound()); i++) {
			if (i == 0)
				// already have one at zero
				continue;
			DefaultXY_DataSet line = new DefaultXY_DataSet();
			line.set(xRange.getLowerBound(), (double)i);
			line.set(xRange.getUpperBound(), (double)i);
			funcs.add(0, line);
			chars.add(0, new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.LIGHT_GRAY));
		}
		
		xyzSpec.setXYElems(funcs);
		xyzSpec.setXYChars(chars);
		xyzSpec.setPlotAnnotations(buildAnnotations(true));
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
		System.out.println("Plotting with xRange: "+xRange);
		xyzGP.drawPlot(xyzSpec, false, false, xRange, yRange);
		// write plot
		xyzGP.getYAxis().setStandardTickUnits(getYTickUnits());
		xyzGP.getChartPanel().setSize(plotWidth, plotWidth-100);
		xyzGP.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		if (writePDF)
			xyzGP.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	private List<XYTextAnnotation> buildAnnotations(boolean doLog) {
		List<XYTextAnnotation> anns = new ArrayList<>();
		
		double yPos = yRange.getUpperBound() * 0.95;
		double yDelta = yPos*0.12;
		
		double xMax = xRange.getUpperBound();
		double xMin = xRange.getLowerBound();
		if (logX) {
			xMax = Math.log10(xMax);
			xMin = Math.log10(xMin);
		}
		double xPos = xMin + 0.95*(xMax - xMin);
		if (logX && !doLog)
			xPos = Math.pow(10, xPos);
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		Color backgroundPaint = new Color(255, 255, 255, 255/2);
		
		XYTextAnnotation ann = new XYTextAnnotation("Mean Residual: "+annDF.format(residualMean), xPos, yPos);
		ann.setBackgroundPaint(backgroundPaint);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(font);
		anns.add(ann);
		
		yPos -= yDelta;
		ann = new XYTextAnnotation("Mean σ: "+annDF.format(dataSigma), xPos, yPos);
		ann.setBackgroundPaint(backgroundPaint);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(font);
		anns.add(ann);
		
		if (gmpeSigma != null) {
			double meanGMPE = 0d;
			for (Point2D pt : gmpeSigma)
				meanGMPE += pt.getY();
			meanGMPE /= gmpeSigma.size();
			yPos -= yDelta;
			ann = new XYTextAnnotation("GMPE Mean σ: "+annDF.format(meanGMPE), xPos, yPos);
			ann.setBackgroundPaint(backgroundPaint);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			ann.setFont(font);
			anns.add(ann);
		}
		
		// now regression top left
		yPos = yRange.getUpperBound() * 0.95;
		xPos = xMin + 0.05*(xMax - xMin);
		if (logX && !doLog)
			xPos = Math.pow(10, xPos);

		ann = new XYTextAnnotation("Linear LSQ Fit", xPos, yPos);
		ann.setBackgroundPaint(backgroundPaint);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(font);
		anns.add(ann);

		yPos -= yDelta;
		ann = new XYTextAnnotation("Intercept: "+regressionDF.format(regressionIntercept), xPos, yPos);
		ann.setBackgroundPaint(backgroundPaint);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(font);
		anns.add(ann);
		
		yPos -= yDelta;
		ann = new XYTextAnnotation("Slope: "+regressionDF.format(regressionLinearSlope), xPos, yPos);
		ann.setBackgroundPaint(backgroundPaint);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(font);
		anns.add(ann);
		
		return anns;
	}
	
	private static final DecimalFormat annDF = new DecimalFormat("0.00");
	private static final DecimalFormat regressionDF = new DecimalFormat("0.00E0");
	
	private void addResidualFuncs(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars, boolean scatter) {
		Color linearColor = Color.BLACK;
		Color binnedColor = Color.BLACK;
		Color gmpeColor = new Color(35, 72, 132);
		
		float residualThickness = 4f;
		float plusMinusThickness = 2f;
		
		// linear scatter
		regressionLinear.setName("Linear LS Fit");
		if (!scatter && logX)
			funcs.add(getInLogX(regressionLinear));
		else
			funcs.add(regressionLinear);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, residualThickness, linearColor));
		
		DiscretizedFunc binnedPoints = new ArbitrarilyDiscretizedFunc();
		binnedPoints.setName("Binned LS Fits");
		
		funcs.add(binnedPoints);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, residualThickness*4, binnedColor));
		
		boolean first = true;
		for (int i=0; i<regressionBinned.length; i++) {
			DiscretizedFunc regression = regressionBinned[i];
			if (regression == null)
				continue;
			
			double binStart = regression.getMinX();
			double binEnd = regression.getMaxX();
			double binMiddle = regressionBins.getX(i);
			if (logX)
				binMiddle = Math.pow(10, binMiddle);
			double binY = regressionBinned[i].getInterpolatedY(binMiddle);
			
			double binFractEachSide = 0.2;
			double tickStart, tickEnd;
			if (logX) {
				double s = Math.log10(binStart);
				double e = Math.log10(binEnd);
				double m = Math.log10(binMiddle);
				double delta = s - e;
				tickStart = Math.pow(10, m - binFractEachSide*delta);
				tickEnd = Math.pow(10, m + binFractEachSide*delta);
			} else {
				double delta = binEnd - binStart;
				tickStart = binMiddle - binFractEachSide*delta;
				tickEnd = binMiddle + binFractEachSide*delta;
			}
			
			if (!scatter && logX) {
				binMiddle = Math.log10(binMiddle);
				tickStart = Math.log10(tickStart);
				tickEnd = Math.log10(tickEnd);
			}
			
			binnedPoints.set(binMiddle, binY);
			addSigmaFuncs(funcs, chars, binMiddle, binY, tickStart, tickEnd, dataSigmaBinned[i],
					true, PlotLineType.SOLID, plusMinusThickness, binnedColor);
			if (first)
				funcs.get(funcs.size()-1).setName("± σ");
			
			if (gmpeSigma != null) {
				double sigma = gmpeSigma.getY(i);
				addSigmaFuncs(funcs, chars, binMiddle, binY, tickStart, tickEnd, sigma,
						false, PlotLineType.DOTTED, plusMinusThickness, gmpeColor);
				if (first)
					funcs.get(funcs.size()-1).setName("± GMPE-σ");
			}
			first = false;
		}
	}
	
	private void addSigmaFuncs(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars,
			double centerX, double centerY, double startX, double endX, double sigma,
			boolean vertLine, PlotLineType lineType, float thickness, Color color) {
		ArbitrarilyDiscretizedFunc plusSigma = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc minusSigma = new ArbitrarilyDiscretizedFunc();
		
		plusSigma.set(startX, centerY + sigma);
		plusSigma.set(endX, centerY + sigma);
		
		minusSigma.set(startX, centerY - sigma);
		minusSigma.set(endX, centerY - sigma);
		
		PlotCurveCharacterstics plotChar = new PlotCurveCharacterstics(lineType, thickness, color);
		
		if (vertLine) {
			XY_DataSet vert = new DefaultXY_DataSet();
			vert.set(centerX, centerY - sigma);
			vert.set(centerX, centerY + sigma);
			
			funcs.add(vert);
			chars.add(plotChar);
		}
		
		funcs.add(plusSigma);
		chars.add(plotChar);
		
		funcs.add(minusSigma);
		chars.add(plotChar);
	}
	
	private DiscretizedFunc getInLogX(DiscretizedFunc xy) {
		ArbitrarilyDiscretizedFunc logXY = new ArbitrarilyDiscretizedFunc();
		logXY.setName(xy.getName());
		for (Point2D pt : xy)
			logXY.set(Math.log10(pt.getX()), pt.getY());
		return logXY;
	}
	
	private static Range getYRange(XY_DataSet xy) {
//		double yExtent = Math.max(Math.abs(xy.getMaxY()), Math.abs(xy.getMinY()));
//		if (yExtent < 0.5)
//			yExtent = yExtent * 1.2;
//		else
//			yExtent = Math.ceil(yExtent);
//		return new Range(-yExtent, yExtent);
		return new Range(-3, 3);
	}
	
	private static TickUnits getYTickUnits() {
		TickUnits tus = new TickUnits();
		TickUnit tu = new NumberTickUnit(0.5);
		tus.add(tu);
		return tus;
	}
	
	private static SimpleRegression getRegression(XY_DataSet xy, boolean logX, Range xRange) {
		// regression
		SimpleRegression regression = new SimpleRegression();
		for (Point2D pt : xy) {
			if (!xRange.contains(pt.getX()))
				continue;
			if (logX)
				pt = new Point2D.Double(Math.log10(pt.getX()), pt.getY());
			regression.addData(pt.getX(), pt.getY());
		}
		if (regression.getN() == 0l)
			return null;
		return regression;
	}
	
	private static ArbitrarilyDiscretizedFunc getRegressionFit(XY_DataSet xy, boolean logX, Range xRange) {
		SimpleRegression regression = getRegression(xy, logX, xRange);
		if (regression == null)
			return null;
		return getRegressionFit(regression, xy, logX, xRange);
	}
	
	private static ArbitrarilyDiscretizedFunc getRegressionFit(SimpleRegression regression, XY_DataSet xy, boolean logX, Range xRange) {
		double b = regression.getIntercept();
		double m = regression.getSlope();
		
		ArbitrarilyDiscretizedFunc fit = new ArbitrarilyDiscretizedFunc();
		// use one to one for x values
		double[] regressXVals = { xRange.getLowerBound(), xRange.getCentralValue(), xRange.getUpperBound() };
		
		if (Double.isNaN(b) || Double.isNaN(m)) {
			// only 1 x value
			double myX = Double.NaN;
			SummaryStatistics stats = new SummaryStatistics();
			for (Point2D pt : xy) {
				if (xRange.contains(pt.getX())) {
					Preconditions.checkState((float)pt.getX() == (float)myX || Double.isNaN(myX));
					stats.addValue(pt.getY());
					myX = pt.getX();
				}
			}
			Preconditions.checkState(!Double.isNaN(myX));
			// mean with 1 x-value minimizes mean squared error
			double meanY = stats.getMean();
			
			for (double x : regressXVals)
				fit.set(x, meanY);
		} else {
			for (int i=0; i<regressXVals.length; i++) {
				double origX = regressXVals[i];
				double x;
				if (logX)
					x = Math.log10(origX);
				else
					x = origX;
				double y = m*x + b;
				fit.set(origX, y);
			}
		}
		
		return fit;
	}
	
	private double calcSigma(DiscretizedFunc regression, XY_DataSet xy) {
		double minX = regression.getMinX();
		double maxX = regression.getMaxX();
		if (logX)
			regression = getInLogX(regression);
		SummaryStatistics stats = new SummaryStatistics();
		for (Point2D pt : xy) {
			if (pt.getX() >= minX && pt.getX() <= maxX) {
				double regressVal;
				if (logX)
					regressVal = regression.getInterpolatedY(Math.log10(pt.getX()));
				else
					regressVal = regression.getInterpolatedY(pt.getX());
				stats.addValue(pt.getY() - regressVal);
			}
		}
		return stats.getStandardDeviation();
	}

	public void setPlotWidth(int plotWidth) {
		this.plotWidth = plotWidth;
	}

	public void setWritePDF(boolean writePDF) {
		this.writePDF = writePDF;
	}

	public void setMaxScatterPoints(int maxScatterPoints) {
		this.maxScatterPoints = maxScatterPoints;
	}
	
	public static <E> void plotPeriodDependentSigma(File outputDir, String prefix,
			Map<Site, List<RuptureComparison<E>>> siteCompsMap, SimulationRotDProvider<E> simProv, double... periods) throws IOException {
		List<Variance[]> siteResidualVariances = new ArrayList<>();
		Map<RuptureComparison<E>, Variance[][]> rupResidualVariances = new HashMap<>();
		Variance[] totalVars = new Variance[periods.length];
		for (int p=0; p<periods.length; p++)
			totalVars[p] = new Variance();
		int rupSimCount = 0;
		for (Site site : siteCompsMap.keySet()) {
			Variance[] siteVars = new Variance[periods.length];
			for (int p=0; p<periods.length; p++)
				siteVars[p] = new Variance();
			siteResidualVariances.add(siteVars);
			for (RuptureComparison<E> comp : siteCompsMap.get(site)) {
				double[] gmpeVals = new double[periods.length];
				for (int p=0; p<periods.length; p++)
					gmpeVals[p] = comp.getLogMean(site, periods[p]);
				E rupture = comp.getRupture();
				int numSims = simProv.getNumSimulations(site, rupture);
				Variance[][] rupVars = rupResidualVariances.get(comp);
				if (rupVars == null) {
					rupSimCount += numSims;
					rupVars = new Variance[numSims][periods.length];
					for (int i=0; i<numSims; i++)
						for (int p=0; p<periods.length; p++)
							rupVars[i][p] = new Variance();
					rupResidualVariances.put(comp, rupVars);
				}
				Preconditions.checkState(numSims > 0);
				for (int i=0; i<numSims; i++) {
					DiscretizedFunc simRD50 = simProv.getRotD50(site, rupture, i);
					for (int p=0; p<periods.length; p++) {
						double simVal = Math.log(simRD50.getY(periods[p]));
						double residual = simVal - gmpeVals[p];
						rupVars[i][p].increment(residual);
						siteVars[p].increment(residual);
						totalVars[p].increment(residual);
					}
				}
			}
		}
		
		// mean between event variances, tau^2
		double[] tauSq = new double[periods.length];
		for (int p=0; p<periods.length; p++) {
			double[] allVars = new double[siteResidualVariances.size()];
			for (int i=0; i<allVars.length; i++)
				allVars[i] = siteResidualVariances.get(i)[p].getResult();
			tauSq[p] = StatUtils.mean(allVars);
		}
		
		// mean within event variances, phi^2
		double[] phiSq = new double[periods.length];
		for (int p=0; p<periods.length; p++) {
			double[] allVars = new double[rupSimCount];
			int index = 0;
			for (RuptureComparison<E> comp : rupResidualVariances.keySet()) {
				Variance[][] rupVars = rupResidualVariances.get(comp);
				for (int i=0; i<rupVars.length; i++)
					allVars[index++] = rupVars[i][p].getResult();
			}
			Preconditions.checkState(index == allVars.length, "Bad index? %s != %s", index, allVars.length);
			phiSq[p] = StatUtils.mean(allVars);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		List<XYAnnotation> anns = new ArrayList<>();
		
		double barWidth = 0.35;
		double centerBuffer = 0.06;
		Color tauColor = Color.RED.darker();
		Color phiColor = Color.BLUE.darker();
		Font valFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		Color periodColor = Color.BLACK;
		Font periodFont = new Font(Font.SANS_SERIF, Font.BOLD, 30);
		PlotCurveCharacterstics sumSigmaChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK);
		PlotCurveCharacterstics sigmaChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY);
		
		Range xRange = new Range(0d, periods.length);
		Range yRange = new Range(0d, 1d);
		
		double centerHalfBuffer = 0.5*centerBuffer;
		
		for (int p=0; p<periods.length; p++) {
			double tau = Math.sqrt(tauSq[p]);
			double phi = Math.sqrt(phiSq[p]);
			double sigma = Math.sqrt(tauSq[p] + phiSq[p]);
			double totalSigma = Math.sqrt(totalVars[p].getResult());
			
			double curCenter = 0.5d + p;
			
			double tauX0 = curCenter - centerHalfBuffer - barWidth;
			addBoxAnns(anns, tau, tauColor, "τ", valFont, tauX0, barWidth);
			
			double phiX0 = curCenter + centerHalfBuffer;
			addBoxAnns(anns, phi, phiColor, "ϕ", valFont, phiX0, barWidth);
			
			DefaultXY_DataSet sumSigmaXY = new DefaultXY_DataSet();
			sumSigmaXY.set(tauX0, sigma);
			sumSigmaXY.set(phiX0 + barWidth, sigma);
			funcs.add(sumSigmaXY);
			chars.add(sumSigmaChar);
			
			DefaultXY_DataSet sigmaXY = new DefaultXY_DataSet();
			sigmaXY.set(tauX0, totalSigma);
			sigmaXY.set(phiX0 + barWidth, totalSigma);
			funcs.add(0, sigmaXY);
			chars.add(0, sigmaChar);
			
			XYTextAnnotation sigmaAnn = new XYTextAnnotation("√(τ²+ϕ²)="+annDF.format(sigma), curCenter, sigma);
			sigmaAnn.setFont(valFont);
			sigmaAnn.setTextAnchor(TextAnchor.BOTTOM_CENTER);
			sigmaAnn.setPaint(sumSigmaChar.getColor());
			anns.add(sigmaAnn);
			
			String periodStr;
			if (periods[p] == Math.round(periods[p]))
				periodStr = (int)periods[p]+"s";
			else
				periodStr = (float)periods[p]+"s";
			XYTextAnnotation periodAnn = new XYTextAnnotation(periodStr, curCenter, yRange.getUpperBound()*0.95);
			periodAnn.setFont(periodFont);
			periodAnn.setTextAnchor(TextAnchor.TOP_CENTER);
			periodAnn.setPaint(periodColor);
			anns.add(periodAnn);
			
			System.out.println("T="+periodStr+"\tτ="+annDF.format(tau)+"\tϕ="+annDF.format(phi)
				+"\tsqrt(τ*τ+ϕ*ϕ)="+annDF.format(sigma)+"\tσ="+annDF.format(totalSigma));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Sigma Components", "Period", "Standard Deviations");
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setLegendFontSize(18);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getXAxis().setTickLabelsVisible(false);
		gp.getPlot().setDomainGridlinesVisible(false);
		gp.getChartPanel().setSize(700, 550);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	private static void addBoxAnns(List<XYAnnotation> anns, double value, Color color, String symbol, Font font, double x0, double width) {
		XYBoxAnnotation box = new XYBoxAnnotation(x0, 0, x0+width, value, null, null, color);
		anns.add(box);
		double boxCenter = x0 + 0.5*width;
		XYTextAnnotation text = new XYTextAnnotation(symbol+"="+annDF.format(value), boxCenter, value);
		text.setFont(font);
		text.setTextAnchor(TextAnchor.BOTTOM_CENTER);
		text.setPaint(color);
		anns.add(text);
	}

}
