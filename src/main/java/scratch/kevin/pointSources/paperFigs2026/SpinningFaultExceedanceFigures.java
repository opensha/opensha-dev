package scratch.kevin.pointSources.paperFigs2026;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Leonard2010_MagLengthRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.calc.RuptureExceedProbCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RectangularSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import net.mahdilamb.colormap.Colors;

public class SpinningFaultExceedanceFigures {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("C:\\Users\\kmilner\\Downloads");
		
		double mag = 7d;
		FocalMech mech = FocalMech.STRIKE_SLIP;
		MagLengthRelationship ml = Leonard2010_MagLengthRelationship.STRIKE_SLIP;
		String prefix = "m7_ss";
		String label = "M7, Strike-Slip";
		
//		double mag = 7d;
//		FocalMech mech = FocalMech.REVERSE;
//		MagLengthRelationship ml = Leonard2010_MagLengthRelationship.DIP_SLIP;
//		String prefix = "m7_rev";
//		String label = "M7, Reverse";

		double period = 0d;
//		double[] distances = {25d};
		double[] distances = {10d, 20d, 40d};
		Range imlRange = new Range(1e-2, 3e0);

		double rake = mech.rake();
		double dip = mech.dip();
		double upperDepth = 1d;
		double lowerDepth = 14d;
		double length = ml.getMedianLength(mag);
		
		System.out.println("Length is: "+length);
		ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
		int numCenteredCalcSurfs = 360;
		int numPlotSurfs = 36;
		EvenlyDiscretizedFunc log10XVals =  new EvenlyDiscretizedFunc(Math.log10(imlRange.getLowerBound()),
				Math.log10(imlRange.getUpperBound()), 100);
		DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<log10XVals.size(); i++) {
			xVals.set(Math.pow(10, log10XVals.getX(i)), 0d);
			logXVals.set(Math.log(xVals.getX(i)), 0d);
		}
		
		List<PointSourceDistanceCorrection> corrs = new ArrayList<>();
		List<String> corrLabels = new ArrayList<>();
		List<Color> corrColors = new ArrayList<>();
		
		corrs.add(PointSourceDistanceCorrections.AVERAGE_SPINNING.get());
		corrLabels.add("Single Distance Correction");
		corrColors.add(Colors.tab_orange);
		
//		corrs.add(PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG.get());
//		corrLabels.add("Proposed Distance Correction");
//		corrColors.add(Colors.tab_purple);
		
		Location gridLoc = new Location(0d, 0d);
		
		EqkRupture rup = new EqkRupture(mag, rake, null, null);
		
		DecimalFormat oDF = new DecimalFormat("0.#");
		
		String xAxisLabel;
		String perPrefix;
		if (period == 0d) {
			perPrefix = "pga";
			gmm.setIntensityMeasure(PGA_Param.NAME);
			xAxisLabel = "Peak Ground Acceleration (g)";
		} else {
			perPrefix = oDF.format(period)+"s";
			gmm.setIntensityMeasure(SA_Param.NAME);
			SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
			xAxisLabel = oDF.format(period)+"s Spectral Acceleration (g)";
		}
		
		PointSurfaceBuilder surfBuilder = new PointSurfaceBuilder(gridLoc)
				.magnitude(mag).dip(dip).upperDepth(upperDepth).lowerDepth(lowerDepth).length(length);
		
		surfBuilder.strike(0d);
		surfBuilder.fractionalDAS(0.5d);
		surfBuilder.fractionalHypocentralDepth(0.5d);
		PointSurface ptSurf = surfBuilder.buildPointSurface();
		RectangularSurface[] centeredCalcSurfs = surfBuilder.buildRandRectSurfaces(numCenteredCalcSurfs);
		RectangularSurface[] plotSurfs = surfBuilder.buildRandRectSurfaces(numPlotSurfs);
		
		surfBuilder.sampleDASs().sampleHypocentralDepths();
		RectangularSurface[] uncenteredCalcSurfs = surfBuilder.buildRandRectSurfaces(numCenteredCalcSurfs*10);
		
		List<PlotSpec> plots = new ArrayList<>();
		for (double distance : distances) {
			System.out.println("Calculating for "+distance+" km");
			Location siteLoc = LocationUtils.location(gridLoc, 0d, distance);
			
			Site site = new Site(siteLoc);
			site.addParameterList(gmm.getSiteParams());
			
			gmm.setSite(site);
			
			DiscretizedFunc[] centeredExceedProbs = calcExceedProbs(rup, gmm, xVals, centeredCalcSurfs);
			DiscretizedFunc[] plotExceedProbs = calcExceedProbs(rup, gmm, xVals, plotSurfs);
			DiscretizedFunc[] uncenteredExceedProbs = calcExceedProbs(rup, gmm, xVals, uncenteredCalcSurfs);
			
			DiscretizedFunc minFunc = min(uncenteredExceedProbs);
			DiscretizedFunc maxFunc = max(uncenteredExceedProbs);
			DefaultXY_DataSet minMax = new DefaultXY_DataSet();
			
			minMax.set(imlRange.getLowerBound(), 1d);
			for (int i=0; i<minFunc.size(); i++) {
				double x = minFunc.getX(i);
				double y = minFunc.getY(i);
				if (imlRange.contains(x))
					minMax.set(x, y);
			}
			
			minMax.set(imlRange.getUpperBound(), 0d);
			for (int i=maxFunc.size(); --i>=0;) {
				double x = maxFunc.getX(i);
				double y = maxFunc.getY(i);
				if (imlRange.contains(x))
					minMax.set(x, y);
			}
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DefaultXY_DataSet fakeMinMax = new DefaultXY_DataSet();
			fakeMinMax.set(1e-10, 2);
			fakeMinMax.setName("Uncentered Range");
			funcs.add(fakeMinMax);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 4f, new Color(0, 0, 0, 40)));
			funcs.add(minMax);
			chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, new Color(0, 0, 0, 40)));
			
//			if (mech == FocalMech.STRIKE_SLIP) {
				for (int i=0; i<plotSurfs.length; i++) {
					DiscretizedFunc func = plotExceedProbs[i];
					if (i == 0)
						func.setName("Centered Surfaces");
					funcs.add(func);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
				}
//			} else {
//				for (boolean hw : new boolean[] {true,false}) {
//					boolean first = true;
//					for (int i=0; i<plotSurfs.length; i++) {
//						boolean isHW = plotSurfs[i].getDistanceX(siteLoc) >= 0d;
//						if (isHW != hw)
//							continue;
//						
//						DiscretizedFunc func = plotExceedProbs[i];
//						funcs.add(func);
//						if (hw) {
//							if (first)
//								func.setName("Centered Hanging-Wall Surfaces");
//							chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.DARK_GRAY));
//						} else {
//							if (first)
//								func.setName("Centered Footwall Surfaces");
//							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
//						}
//						first = false;
//					}
//				}
//			}
			
			float avgThickness = 4f;
			
			DiscretizedFunc centeredAvg = average(centeredExceedProbs);
			centeredAvg.setName("Centered Average");
			funcs.add(centeredAvg);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_blue));
			
			DiscretizedFunc uncenteredAvg = average(uncenteredExceedProbs);
			uncenteredAvg.setName("Uncentered Average");
			funcs.add(uncenteredAvg);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Color.BLACK));
			
			for (int c=0; c<corrs.size(); c++) {
				PointSourceDistanceCorrection corr = corrs.get(c);
				
				// calculate distance corrected
				rup.setRuptureSurface(ptSurf.getDistancedProtected(corr, TectonicRegionType.ACTIVE_SHALLOW, mag));
				RuptureExceedProbCalculator.calcExceedanceProbabilities(gmm, rup, logXVals);
				DiscretizedFunc corrExceedProbs = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<logXVals.size(); i++)
					corrExceedProbs.set(xVals.getX(i), logXVals.getY(i));
				corrExceedProbs.setName(corrLabels.get(c));
				funcs.add(corrExceedProbs);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, avgThickness, corrColors.get(c)));
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, label, xAxisLabel, "Exceedance Probability");
			if (plots.isEmpty())
				plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			
			double annX = xVals.getX((int)(xVals.size()*0.95));
			double annY = 0.95;
			XYTextAnnotation distAnn = new XYTextAnnotation(oDF.format(distance)+" km", annX, annY);
			distAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 34));
			distAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			plot.addPlotAnnotation(distAnn);
			
			plots.add(plot);
		}
		
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<plots.size(); i++)
			yRanges.add(new Range(0d, 1d));
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plots, true, false, List.of(imlRange), yRanges);
		
		PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix, gp, 1000, 200+300*plots.size(), true, true, false);
	}
	
	private static DiscretizedFunc[] calcExceedProbs(EqkRupture rup, ScalarIMR gmm,
			DiscretizedFunc xVals, RectangularSurface[] surfs) {
		System.out.println("Calculating exceedance probabilities for "+surfs.length+" surfaces");
		double[] linearX = new double[xVals.size()];
		double[] logX = new double[xVals.size()];
		for (int i=0; i<logX.length; i++) {
			linearX[i] = xVals.getX(i);
			logX[i] = Math.log(linearX[i]);
		}
		DiscretizedFunc logXVals = new LightFixedXFunc(logX, new double[logX.length]);
		
		DiscretizedFunc[] ret = new DiscretizedFunc[surfs.length];
		for (int r=0; r<surfs.length; r++) {
			rup.setRuptureSurface(surfs[r]);
			gmm.setEqkRupture(rup);
			gmm.getExceedProbabilities(logXVals);
			LightFixedXFunc probs = new LightFixedXFunc(linearX, new double[linearX.length]);
			for (int i=0; i<linearX.length; i++)
				probs.set(i, logXVals.getY(i));
			ret[r] = probs;
		}
		return ret;
	}
	
	private static DiscretizedFunc average(DiscretizedFunc[] funcs) {
		double scalar = 1d/funcs.length;
		double[] xVals = null;
		double[] yVals = null;
		
		for (DiscretizedFunc func : funcs) {
			if (xVals == null) {
				xVals = new double[func.size()];
				yVals = new double[func.size()];
				for (int i=0; i<func.size(); i++)
					xVals[i] = func.getX(i);
			}
			for (int i=0; i<func.size(); i++)
				yVals[i] = Math.fma(scalar, func.getY(i), yVals[i]);
		}
		return new LightFixedXFunc(xVals, yVals);
	}
	
	private static DiscretizedFunc min(DiscretizedFunc[] funcs) {
		double minVal = Double.POSITIVE_INFINITY;
		DiscretizedFunc minFunc = null;
		for (DiscretizedFunc func : funcs) {
			double val = func.getFirstInterpolatedX(0.5);
			if (val < minVal) {
				minVal = val;
				minFunc = func;
			}
		}
		return minFunc;
	}
	
	private static DiscretizedFunc max(DiscretizedFunc[] funcs) {
		double maxVal = Double.NEGATIVE_INFINITY;
		DiscretizedFunc maxFunc = null;
		for (DiscretizedFunc func : funcs) {
			double val = func.getFirstInterpolatedX(0.5);
			if (val > maxVal) {
				maxVal = val;
				maxFunc = func;
			}
		}
		return maxFunc;
	}

}
