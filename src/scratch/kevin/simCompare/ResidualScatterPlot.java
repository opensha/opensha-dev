package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
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
	double regressionLinearSlope;
	private DiscretizedFunc regressionLinear;
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
		regressionLinearSlope = linearRegression.getSlope();
		regressionLinear = getRegressionFit(linearRegression, xy, logX, xRange);
		dataSigma = calcSigma(regressionLinear, xy);
		
		regressionBinned = new ArbitrarilyDiscretizedFunc[numRegressionBins];
		dataSigmaBinned = new double[numRegressionBins];
		double deltaX = (maxX - minX) / (numRegressionBins - 1);
		EvenlyDiscretizedFunc regressionBins = new EvenlyDiscretizedFunc(minX+0.5*deltaX, numRegressionBins, deltaX);
		System.out.println("Regression bins:\n"+regressionBins);
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
			System.out.println("GMPE Sigma function:\n"+gmpeSigma);
		}
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
		// add W-C
		ArrayList<XY_DataSet> funcs = new ArrayList<>();
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		addResidualFuncs(funcs, chars, false);
		
		Range xRange = this.xRange;
		if (logX)
			xRange = new Range(Math.log10(xRange.getLowerBound()), Math.log10(xRange.getUpperBound()));
		
		xyzSpec.setXYElems(funcs);
		xyzSpec.setXYChars(chars);
		xyzSpec.setPlotAnnotations(buildAnnotations(true));
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
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
		if (logX && doLog) {
			xMax = Math.log10(xMax);
			xMin = Math.log10(xMin);
		}
		double xPos = xMin + 0.95*(xMax - xMin);
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		
		XYTextAnnotation ann = new XYTextAnnotation("Mean Residual: "+annDF.format(residualMean), xPos, yPos);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(font);
		anns.add(ann);
		
		yPos -= yDelta;
		ann = new XYTextAnnotation("Linear SQ Fit slope: "+annDF.format(regressionLinearSlope), xPos, yPos);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(font);
		anns.add(ann);
		
		yPos -= yDelta;
		ann = new XYTextAnnotation("Mean σ: "+annDF.format(dataSigma), xPos, yPos);
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
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			ann.setFont(font);
			anns.add(ann);
		}
		
		return anns;
	}
	
	private static final DecimalFormat annDF = new DecimalFormat("0.00");
	
	private void addResidualFuncs(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars, boolean scatter) {
		Color linearColor = Color.BLACK;
		Color binnedColor = linearColor;
		Color gmpeColor = new Color(35, 72, 132);
		
		float residualThickness = 4f;
		float plusMinusThickness = 2f;
		
		// linear scatter
		regressionLinear.setName("Linear LS Fit");
		if (!scatter && logX)
			funcs.add(getInLogX(regressionLinear));
		else
			funcs.add(regressionLinear);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, residualThickness, linearColor));
		addSigmaFuncs(funcs, chars, regressionLinear, dataSigma, PlotLineType.DASHED, plusMinusThickness, linearColor, !scatter);
		funcs.get(funcs.size()-1).setName("± σ");
		
		if (gmpeSigma != null) {
			boolean first = true;
			// linear first
			
			double x0 = regressionLinear.getX(0);
			double y0 = regressionLinear.getY(0);
			
			if (logX)
				x0 = Math.log10(x0);
			
			for (int i=0; i<gmpeSigma.size(); i++) {
				double x = gmpeSigma.getX(i);
				double start = x - 0.5*gmpeSigma.getDelta();
				double end = x + 0.5*gmpeSigma.getDelta();
				double regressionLinearSlope = this.regressionLinearSlope;
				double startY;
				double endY;
				startY = (start - x0)*regressionLinearSlope + y0;
				endY = (end - x0)*regressionLinearSlope + y0;
				DiscretizedFunc tempRegression = new ArbitrarilyDiscretizedFunc();
				if (logX) {
					// gmpe sigma func is in log space, back to linear
					start = Math.pow(10, start);
					end = Math.pow(10, end);
				}
				tempRegression.set(start, startY);
				tempRegression.set(end, endY);
				double sigma = gmpeSigma.getY(i);
				addSigmaFuncs(funcs, chars, tempRegression, sigma, PlotLineType.DASHED, plusMinusThickness, gmpeColor, !scatter);
				if (first) {
					funcs.get(funcs.size()-1).setName("± GMPE-σ");
					first = false;
				}
			}
		}
		
		boolean first = true;
		for (int i=0; i<regressionBinned.length; i++) {
			DiscretizedFunc regression = regressionBinned[i];
			if (regression == null)
				continue;
			if (first)
				regressionBinned[i].setName("Binned LS Fits");
			if (!scatter && logX)
				funcs.add(getInLogX(regression));
			else
				funcs.add(regression);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, residualThickness, binnedColor));
			addSigmaFuncs(funcs, chars, regression, dataSigmaBinned[i], PlotLineType.DOTTED, plusMinusThickness, binnedColor, !scatter);
			if (first) {
				funcs.get(funcs.size()-1).setName("± σ");
				first = false;
			}
		}
		
		if (gmpeSigma != null) {
			// now binned
			first = true;
			for (int i=0; i<regressionBinned.length; i++) {
				if (regressionBinned[i] == null)
					continue;
				double sigma = gmpeSigma.getY(i);
				addSigmaFuncs(funcs, chars, regressionBinned[i], sigma, PlotLineType.DOTTED, plusMinusThickness, gmpeColor, !scatter);
				if (first) {
					funcs.get(funcs.size()-1).setName("± GMPE-σ");
					first = false;
				}
			}
		}
	}
	
	private void addSigmaFuncs(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars,
			DiscretizedFunc regression, double sigma, PlotLineType lineType, float thickness, Color color, boolean leaveInLogSpace) {
		if (logX)
			regression = getInLogX(regression);
		ArbitrarilyDiscretizedFunc plusSigma = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc minusSigma = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : regression) {
			double x = pt.getX();
			if (logX && !leaveInLogSpace)
				x = Math.pow(10, x);
			double y = pt.getY();
			plusSigma.set(x, y+sigma);
			minusSigma.set(x, y-sigma);
		}
		
		funcs.add(plusSigma);
		chars.add(new PlotCurveCharacterstics(lineType, thickness, color));
		
		funcs.add(minusSigma);
		chars.add(new PlotCurveCharacterstics(lineType, thickness, color));
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

}
