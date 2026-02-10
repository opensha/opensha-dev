package scratch.kevin.pointSources.paperFigs2026;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Function;
import java.util.function.Supplier;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.netlib.util.intW;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.WeightedValue;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
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
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.RuptureExceedProbCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RectangularSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.SurfaceDistances;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection.FractileBin;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ErgodicIMR;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.latex.LaTeXUtils;

public class SpinningFaultExceedanceFigures {
	
	enum Sortables {
		RJB("Rjb", D -> D.getDistanceJB(), Colors.tab_blue),
		RRUP("Rrup", D -> D.getDistanceRup(), Colors.tab_green),
		RRUP_PLUS_RJB("Rjb + Rrup", D -> D.getDistanceJB()+D.getDistanceRup(), Color.DARK_GRAY);
		
		private String label;
		private Function<SurfaceDistances, Double> function;
		private Color color;
		private Color transColor;

		private Sortables(String Label, Function<SurfaceDistances, Double> function, Color color) {
			label = Label;
			this.function = function;
			this.color = color;
			this.transColor = new Color(color.getRed(), color.getGreen(), color.getBlue(), 80);
		}
	};

	public static void main(String[] args) throws IOException {
		File outputDir = new File(ConstantsAndSettings.FIGURES_DIR, "spinning_exceedance_probs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
//		FocalMech mech = FocalMech.STRIKE_SLIP;
		FocalMech mech = FocalMech.REVERSE;
		
		double mag;
		MagLengthRelationship ml;
		String prefix, geomLabel, imLabel, texLabel;
		
		switch (mech) {
		case STRIKE_SLIP:
			mag = 7.05d;
			ml = new WC1994_MagLengthRelationship();
//			ml = Leonard2010_MagLengthRelationship.STRIKE_SLIP;
			prefix = "m7_ss";
//			geomLabel = "M7.05, Vertical";
			imLabel = "M7.05, Strike-Slip";
			geomLabel = imLabel;
			texLabel = "MSevenSS";
			break;
		case REVERSE:
			mag = 7.05d;
			ml = new WC1994_MagLengthRelationship();
//			ml = Leonard2010_MagLengthRelationship.DIP_SLIP;
			prefix = "m7_rev";
//			geomLabel = "M7.05, Dipping";
			imLabel = "M7.05, Reverse";
			geomLabel = imLabel;
			texLabel = "MSevenRev";
			break;

		default:
			throw new IllegalStateException("Unexpected mech: "+mech);
		}
		
		Sortables[] sortables = Sortables.values();
		Sortables SORT_QUANTITY = Sortables.RRUP_PLUS_RJB;

		double period = 0d;
//		double[] distances = {25d};
//		double[] distances = {0d, 10d, 20d, 40d};
//		String[] texDistPrefixes = {"Zero", "Ten", "Twenty", "Forty"};
//		double[] distances = {10d, 20d, 40d};
//		String[] texDistPrefixes = {"Ten", "Twenty", "Forty"};
		double[] distances = {10d, 25d, 50d};
		String[] texDistPrefixes = {"Ten", "TwentyFive", "Fifty"};
		
		Preconditions.checkState(texDistPrefixes.length == distances.length);
		Range imlRange = new Range(3e-2, 3e0);
		Range medianIMLRange1 = new Range(1e-1, 1e-0);
		Range medianIMLRange2= new Range(5e-2, 5e-1);
		Range[] medianIMLRanges = new Range[distances.length];
		for (int d=0; d<distances.length; d++) {
			medianIMLRanges[d] = distances[d] > 30d ? medianIMLRange2 : medianIMLRange1;
		}

		double rake = mech.rake();
		double dip = mech.dip();
		double upperDepth = 1d;
		double lowerDepth = 14d;
		double length = ml.getMedianLength(mag);
		
		System.out.println("Length is: "+length);
		Supplier<ErgodicIMR> gmmRef = () -> {
			ErgodicIMR gmm = (ErgodicIMR)AttenRelRef.USGS_NSHM23_ACTIVE.get();
			gmm.getOtherParams().getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
			gmm.getOtherParams().getParameter(SigmaTruncLevelParam.NAME).setValue(3d);
			if (period == 0d) {
				gmm.setIntensityMeasure(PGA_Param.NAME);
			} else {
				gmm.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
			}
			return gmm;
		};
		
		int numCenteredCalcSurfs = 36000;
		int numUncenteredCalcSurfs = numCenteredCalcSurfs*100;
		int numPlotSurfs = 36;
		int outlineMapSurfCount = 360*2;
		EvenlyDiscretizedFunc log10DiscrXVals =  new EvenlyDiscretizedFunc(Math.log10(imlRange.getLowerBound()),
				Math.log10(imlRange.getUpperBound()), 100);
		DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<log10DiscrXVals.size(); i++) {
			xVals.set(Math.pow(10, log10DiscrXVals.getX(i)), 0d);
			logXVals.set(Math.log(xVals.getX(i)), 0d);
		}
		
		FileWriter texFW = new FileWriter(new File(outputDir, prefix+".tex"));
		
		List<PointSourceDistanceCorrection> corrs = new ArrayList<>();
		List<String> corrLabels = new ArrayList<>();
		List<PlotCurveCharacterstics> corrChars = new ArrayList<>();
		
		PointSourceDistanceCorrection nshm23Corr = PointSourceDistanceCorrections.NSHM_2013.get();
		String nshm23Name = "NSHM23 as-published";
		PointSourceDistanceCorrection averageCorr = PointSourceDistanceCorrections.AVERAGE_SPINNING_CENTERED.get();
		String averageName = "Average centered";
		DistanceDistributionCorrection proposedCorr = (DistanceDistributionCorrection)
				PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST.get();
		
		float avgThickness = 4f;
		
		corrs.add(nshm23Corr);
		corrLabels.add(nshm23Name);
		corrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_brown));
		
		corrs.add(averageCorr);
		corrLabels.add(averageName+" correction");
		corrChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, avgThickness, Colors.tab_orange));
		
		corrs.add(PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST.get());
		corrLabels.add("Proposed correction");
		corrChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, avgThickness, Colors.tab_purple));
		
		Location gridLoc = new Location(0d, 0d);
		
		System.out.println("Grid location: "+gridLoc);
		
		EqkRupture rup = new EqkRupture(mag, rake, null, null);
		
		String imlAxisLabel;
		String perPrefix;
		if (period == 0d) {
			perPrefix = "pga";
			imlAxisLabel = "Peak Ground Acceleration (g)";
		} else {
			perPrefix = oDF.format(period)+"s";
			imlAxisLabel = oDF.format(period)+"s Spectral Acceleration (g)";
		}
		
		PointSurfaceBuilder surfBuilder = new PointSurfaceBuilder(gridLoc)
				.magnitude(mag).dip(dip).upperDepth(upperDepth).lowerDepth(lowerDepth).length(length);
		
		surfBuilder.strike(0d);
		surfBuilder.fractionalDAS(0.5d);
		surfBuilder.fractionalHypocentralDepth(0.5d);
		PointSurface ptSurf = surfBuilder.buildPointSurface();
		
		texFW.write("% Rupture properties\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"Length", oDF.format(length))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"Dip", oDF.format(dip))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"Ztor", oDF.format(upperDepth))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"Zbot", oDF.format(lowerDepth))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"DDW", oDF.format(ptSurf.getAveWidth()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"HorzWidth", oDF.format(ptSurf.getAveHorizontalWidth()))+"\n");
		texFW.write("\n");
		
		RectangularSurface[] centeredCalcSurfs = surfBuilder.buildRandRectSurfaces(numCenteredCalcSurfs);
		RectangularSurface[] plotSurfs = surfBuilder.buildRandRectSurfaces(numPlotSurfs);
		RectangularSurface[] centeredMapOutlineSurfs = surfBuilder.buildRandRectSurfaces(outlineMapSurfCount);
		
		surfBuilder.sampleDASs().sampleHypocentralDepths();
		RectangularSurface[] uncenteredCalcSurfs = surfBuilder.buildRandRectSurfaces(numUncenteredCalcSurfs);
		
		surfBuilder.das(0d).fractionalHypocentralDepth(0.5);
		RectangularSurface[] extremeUncenteredSurfs = surfBuilder.buildRandRectSurfaces(outlineMapSurfCount);
		
		// examples
		surfBuilder.fractionalDAS(0.5).fractionalHypocentralDepth(0.5);
		RectangularSurface[] centeredExampleSurfs;
		if (mech == FocalMech.STRIKE_SLIP) {
			centeredExampleSurfs = surfBuilder.buildRandRectSurfaces(9);
		} else {
			centeredExampleSurfs = new RectangularSurface[] {
					surfBuilder.strike(-45d).buildRectSurface(),
					surfBuilder.strike(45d).buildRectSurface(),
			};
		}
		
		surfBuilder.das(10d).hypocentralDepth(5d).strike(-15);
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"UncenteredExampleDAS", "10")+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(texLabel+"UncenteredExampleZHyp", "5")+"\n");
		RectangularSurface uncenteredExampleSurf = surfBuilder.buildRectSurface();
		
		ErgodicIMR gmm0 = gmmRef.get();
		
		List<PlotSpec> plots = new ArrayList<>();
		for (int d=0; d<distances.length; d++) {
			double distance = distances[d];
			System.out.println("Calculating for "+distance+" km");
			Location siteLoc = LocationUtils.location(gridLoc, 0d, distance);
			
			System.out.println("\tSite location: "+siteLoc);
			
			Site site = new Site(siteLoc);
			site.addParameterList(gmm0.getSiteParams());
			
			gmm0.setSite(site);
			
			DiscretizedFunc[] centeredExceedProbs = calcExceedProbs(site, rup, gmmRef, xVals, centeredCalcSurfs);
			DiscretizedFunc[] plotExceedProbs = calcExceedProbs(site, rup, gmmRef, xVals, plotSurfs);
			DiscretizedFunc[] uncenteredExceedProbs = calcExceedProbs(site, rup, gmmRef, xVals, uncenteredCalcSurfs);
			
			DiscretizedFunc minFunc = min(uncenteredExceedProbs);
			DiscretizedFunc maxFunc = max(uncenteredExceedProbs);
			DefaultXY_DataSet minMax = new DefaultXY_DataSet();
			
			minMax.set(imlRange.getLowerBound(), 1d);
			for (int i=0; i<minFunc.size(); i++) {
				double x = minFunc.getX(i);
				double y = minFunc.getY(i);
				if (imlRange.contains(x) || Precision.equals(x, imlRange.getLowerBound(), 1e-8)
						 || Precision.equals(x, imlRange.getUpperBound(), 1e-8))
					minMax.set(x, y);
			}
			
			minMax.set(imlRange.getUpperBound(), 0d);
			for (int i=maxFunc.size(); --i>=0;) {
				double x = maxFunc.getX(i);
				double y = maxFunc.getY(i);
				if (imlRange.contains(x) || Precision.equals(x, imlRange.getLowerBound(), 1e-8)
						 || Precision.equals(x, imlRange.getUpperBound(), 1e-8))
					minMax.set(x, y);
			}
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			XY_DataSet fakeMinMax = new DefaultXY_DataSet();
			fakeMinMax.set(1e-10, 2000);
			fakeMinMax.setName("Uncentered range");
			Color polyColor = new Color(0, 0, 0, 40);
			funcs.add(fakeMinMax);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 4f, polyColor));
			funcs.add(minMax);
			chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, polyColor));
			
//			if (mech == FocalMech.STRIKE_SLIP) {
				for (int i=0; i<plotSurfs.length; i++) {
					DiscretizedFunc func = plotExceedProbs[i];
					if (i == 0)
						func.setName("Centered individual");
					funcs.add(func);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Colors.tab_lightblue));
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
			
			DiscretizedFunc centeredAvg = average(centeredExceedProbs);
			centeredAvg.setName("Centered average");
			funcs.add(centeredAvg);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_blue));
			
			DiscretizedFunc uncenteredAvg = average(uncenteredExceedProbs);
			uncenteredAvg.setName("Uncentered average");
			funcs.add(uncenteredAvg);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Color.BLACK));
			
			for (int c=0; c<corrs.size(); c++) {
				PointSourceDistanceCorrection corr = corrs.get(c);
				
				WeightedList<SurfaceDistances> corrDists = corr.getCorrectedDistances(
						siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance);
				System.out.println("Distances for "+corr);
				for (WeightedValue<SurfaceDistances> corrDist : corrDists)
					System.out.println("\tweight="+(float)corrDist.weight+"; "+corrDist.value);
				
				// calculate distance corrected
				rup.setRuptureSurface(ptSurf.getDistancedProtected(corr, TectonicRegionType.ACTIVE_SHALLOW, mag));
				RuptureExceedProbCalculator.calcExceedanceProbabilities(gmm0, rup, logXVals);
				DiscretizedFunc corrExceedProbs = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<logXVals.size(); i++)
					corrExceedProbs.set(xVals.getX(i), logXVals.getY(i));
				corrExceedProbs.setName(corrLabels.get(c));
				funcs.add(corrExceedProbs);
				chars.add(corrChars.get(c));
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, imLabel, imlAxisLabel, "Exceedance Probability");
			
			double annX = xVals.getX((int)(xVals.size()*0.97));
			double annY = 0.97;
			Font distFont = new Font(Font.SANS_SERIF, Font.BOLD, 30);
			XYTextAnnotation ann = new XYTextAnnotation(oDF.format(distance)+" km", annX, annY);
			ann.setFont(distFont);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			plot.addPlotAnnotation(ann);

			WeightedList<SurfaceDistances> tmpNshm23Dists = nshm23Corr.getCorrectedDistances(
					siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance);
			SurfaceDistances nshmDistsHW, nshmDistsFW;
			if (tmpNshm23Dists.size() == 1) {
				nshmDistsFW = tmpNshm23Dists.getValue(0);
				nshmDistsHW = null;
			} else {
				nshmDistsHW = null;
				nshmDistsFW = null;
				for (int i=0; i<tmpNshm23Dists.size(); i++) {
					SurfaceDistances dist = tmpNshm23Dists.getValue(i);
					if (dist.getDistanceX() >= 0)
						nshmDistsHW = dist;
					else
						nshmDistsFW = dist;
				}
				Preconditions.checkNotNull(nshmDistsHW);
				Preconditions.checkNotNull(nshmDistsFW);
			}
//			SurfaceDistances nshmDists = nshm23Corr.getCorrectedDistances(
//					siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance).getValue(0);
			WeightedList<SurfaceDistances> averageDists = averageCorr.getCorrectedDistances(
					siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance);
			WeightedList<SurfaceDistances> proposedDists = proposedCorr.getCorrectedDistances(
					siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance);
			
			plots.add(plot);
			
			// Rjb vs Rrup figure
			SurfaceDistances[] centeredDists = calcDists(siteLoc, centeredCalcSurfs);
			DefaultXY_DataSet centeredRupVsJB = calcRrupVsRjb(centeredDists);
			DefaultXY_DataSet centeredRupVsJB_fw = calcRrupVsRjb(centeredDists, false);
			DefaultXY_DataSet centeredRupVsJB_hw = calcRrupVsRjb(centeredDists, true);
			SurfaceDistances[] uncenteredDists = calcDists(siteLoc, uncenteredCalcSurfs);
			DefaultXY_DataSet uncenteredRupVsJB = calcRrupVsRjb(uncenteredDists);
			
//			DefaultXY_DataSet uncenteredPoly = ConvexHullCalc.calcConvexHull(uncenteredRupVsJB);
			
			double minDist = uncenteredRupVsJB.getMinX();
			double maxDist = uncenteredRupVsJB.getMaxY();
			if (minDist < 5)
				minDist = 0;
			else
				minDist = Math.floor(minDist/5d)*5d;
			maxDist = Math.ceil(maxDist/5d)*5d;
			Range distRange = new Range(minDist, maxDist);
			
			DefaultXY_DataSet uncenteredPoly = calcEnvelopePoly(uncenteredRupVsJB, 100);
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();

//			funcs.add(uncenteredRupVsJB);
//			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 1f, Color.BLACK));
//			funcs.add(centeredRupVsJB);
//			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Colors.tab_blue));
			
//			funcs.add(uncenteredPoly);
			funcs.add(fakeMinMax);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 4f, polyColor));
			funcs.add(uncenteredPoly);
			chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, polyColor));
//			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.GRAY));
			
			// failed experiment: contours
//			if (mech != FocalMech.STRIKE_SLIP) {
//				Path2D poly2d = buildPolygon(uncenteredPoly);
//				EvenlyDiscretizedFunc distGridding = new EvenlyDiscretizedFunc(distRange.getLowerBound(), distRange.getUpperBound(), 50);
//				EvenlyDiscrXYZ_DataSet uncenteredInterpXYZ = new EvenlyDiscrXYZ_DataSet(
//						distGridding.size(), distGridding.size(), distGridding.getMinX(), distGridding.getMinX(),
//						distGridding.getDelta());
//				// mask outside of our polygon
//				for (int i=0; i<uncenteredInterpXYZ.size(); i++)
//					if (!poly2d.contains(uncenteredInterpXYZ.getPoint(i)))
//						uncenteredInterpXYZ.set(i, Double.NaN);
//				// interpolate
//				System.out.println("Interpolating Rrup vs Rjb data");
//				ScatterInterpolator.interpolateXYDensity(uncenteredRupVsJB, uncenteredInterpXYZ, 1, 1d, 0);
//				// set NaNs to zero
//				for (int i=0; i<uncenteredInterpXYZ.size(); i++)
//					if (Double.isNaN(uncenteredInterpXYZ.get(i)))
//						uncenteredInterpXYZ.set(i, 0d);
//				uncenteredInterpXYZ.scale(1d/uncenteredInterpXYZ.getSumZ());
//				System.out.println("Interpolated data:");
//				for (int y=uncenteredInterpXYZ.getNumY(); --y>=0;) {
//					for (int x=0; x<uncenteredInterpXYZ.getNumX(); x++)
//						System.out.print(" "+(float)uncenteredInterpXYZ.get(x, y));
//					System.out.println();
//				}
//				
////				EvenlyDiscretizedFunc interpLevels = HistogramFunction.getEncompassingHistogram(
////						0.01, uncenteredInterpXYZ.getMaxZ()*0.9, uncenteredInterpXYZ.getMaxZ()*0.9/5d);
////				double[] levels = new double[interpLevels.size()];
////				for (int i=0; i<levels.length; i++)
////					levels[i] = interpLevels.getX(i);
////				
////				System.out.println("Contour levels:");
////				for (double level : levels)
////					System.out.print(" "+(float)level);
////				System.out.println();
////				double[] levels = { 0.01d, 0.05d, 0.1d };
//				EvenlyDiscretizedFunc invCDF = new EvenlyDiscretizedFunc(-10d, 0d, 100); // log10 units
//				for (int i=0; i<uncenteredInterpXYZ.size(); i++) {
//					double z = uncenteredInterpXYZ.get(i);
//					if (z > 0) {
//						double logZ = Math.log10(z);
//						for (int j=0; j<invCDF.size(); j++) {
//							if (logZ >= invCDF.getX(j))
//								invCDF.add(j, z);
//							else
//								break;
//						}
//					}
//				}
//				for (int i=0; i<invCDF.size(); i++)
//					invCDF.set(i, 1d - invCDF.getY(i));
//				System.out.println("Inv CDF:\n"+invCDF);
//				double[] massFracts = {0.05, 0.1, 0.2, 0.5};
//				List<Double> levels = new ArrayList<>();
//				System.out.println("Contour levels:");
//				for (double massFract : massFracts) {
//					if (massFract < invCDF.getY(0) || massFract > invCDF.getY(invCDF.size()-1))
//						continue;
//					double level = Math.pow(10, invCDF.getFirstInterpolatedX(massFract));
//					System.out.println("\t"+(float)massFract+" -> "+(float)level);
//					levels.add(level);
//				}
//				
//				if (!levels.isEmpty()) {
//					System.out.println("Contouring Rrup vs Rjb data");
//					Map<Double, List<XY_DataSet>> contours = XYZContourGenerator.contours(uncenteredInterpXYZ, Doubles.toArray(levels), false);
//					List<Double> contourLevels = new ArrayList<>(contours.keySet());
//					Collections.sort(contourLevels);
//					PlotCurveCharacterstics contourChar = new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, transColor(polyColor, 40));
//					int contourCount = 0;
//					for (double level : contourLevels) {
//						for (XY_DataSet contour : contours.get(level)) {
//							if (contourCount == 0)
//								System.out.println("Contour0:\n"+contour);
//							contourCount++;
//							funcs.add(contour);
//							chars.add(contourChar);
//						}
//					}
//					System.out.println("Done with "+contourCount+" contours");
//				}
//			}
			
			List<XY_DataSet> rupVsJBXYs = new ArrayList<>();
			List<PlotCurveCharacterstics> rupVsJBChars = new ArrayList<>();
			
			PlotCurveCharacterstics nshm23FWSymbol = new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 6f, Colors.tab_lightbrown);
			PlotCurveCharacterstics nshm23HWSymbol = new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 6f, Colors.tab_brown);
			
			if (mech != FocalMech.STRIKE_SLIP) {
				centeredRupVsJB_fw = thin(centeredRupVsJB_fw, 0.1);
				centeredRupVsJB_fw.setName("Centered, FW");
				funcs.add(centeredRupVsJB_fw);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_lightblue));
				centeredRupVsJB_hw = thin(reverse(centeredRupVsJB_hw), 0.1);
				centeredRupVsJB_hw.setName("Centered, HW");
				funcs.add(centeredRupVsJB_hw);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_blue));
				
				DefaultXY_DataSet nshm23RupVsJB = new DefaultXY_DataSet();
				nshm23RupVsJB.set(nshmDistsFW.getDistanceJB(), nshmDistsFW.getDistanceRup());
				nshm23RupVsJB.setName(nshm23Name+", FW");
				funcs.add(nshm23RupVsJB);
				chars.add(nshm23FWSymbol);
//				rupVsJBXYs.add(funcs.get(funcs.size()-1));
//				rupVsJBChars.add(chars.get(chars.size()-1));
				
				nshm23RupVsJB = new DefaultXY_DataSet();
				nshm23RupVsJB.set(nshmDistsHW.getDistanceJB(), nshmDistsHW.getDistanceRup());
				nshm23RupVsJB.setName(nshm23Name+", HW");
				funcs.add(nshm23RupVsJB);
				chars.add(nshm23HWSymbol);
				rupVsJBXYs.add(funcs.get(funcs.size()-1));
				rupVsJBChars.add(chars.get(chars.size()-1));
			} else {
				centeredRupVsJB.setName("Centered");
				funcs.add(centeredRupVsJB);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_blue));
				
				DefaultXY_DataSet nshm23RupVsJB = new DefaultXY_DataSet();
				nshm23RupVsJB.set(nshmDistsFW.getDistanceJB(), nshmDistsFW.getDistanceRup());
				nshm23RupVsJB.setName(nshm23Name);
				funcs.add(nshm23RupVsJB);
				chars.add(nshm23HWSymbol);
				rupVsJBXYs.add(funcs.get(funcs.size()-1));
				rupVsJBChars.add(chars.get(chars.size()-1));
			}
			
			PlotCurveCharacterstics averageDistHWSymbol = new PlotCurveCharacterstics(PlotSymbol.FILLED_DIAMOND, 5f, Colors.tab_orange);
			PlotCurveCharacterstics averageDistFWSymbol = new PlotCurveCharacterstics(PlotSymbol.FILLED_DIAMOND, 5f, Colors.tab_lightorange);
			for (WeightedValue<SurfaceDistances> dist : averageDists) {
				DefaultXY_DataSet avgRupVsJB = new DefaultXY_DataSet();
				avgRupVsJB.set(dist.value.getDistanceJB(), dist.value.getDistanceRup());
				PlotCurveCharacterstics pChar;
				if (averageDists.size() == 1) {
					avgRupVsJB.setName(averageName);
					pChar = averageDistHWSymbol;
				} else if (dist.value.getDistanceX() >= 0d) {
					avgRupVsJB.setName(averageName+", HW");
					pChar = averageDistHWSymbol;
				} else {
					avgRupVsJB.setName(averageName+", FW");
					pChar = averageDistFWSymbol;
				}
				funcs.add(avgRupVsJB);
				chars.add(pChar);
				rupVsJBXYs.add(funcs.get(funcs.size()-1));
				rupVsJBChars.add(chars.get(chars.size()-1));
			}
			DefaultXY_DataSet proposedRupVsJB_hw = new DefaultXY_DataSet();
			DefaultXY_DataSet proposedRupVsJB = new DefaultXY_DataSet();
			DefaultXY_DataSet proposedRupVsJB_fw = new DefaultXY_DataSet();
			for (WeightedValue<SurfaceDistances> dist : proposedDists) {
				if (dist.value.getDistanceX() >= 0d)
					proposedRupVsJB_hw.set(dist.value.getDistanceJB(), dist.value.getDistanceRup());
				else
					proposedRupVsJB_fw.set(dist.value.getDistanceJB(), dist.value.getDistanceRup());
				proposedRupVsJB.set(dist.value.getDistanceJB(), dist.value.getDistanceRup());
			}
			PlotCurveCharacterstics proposedDistHWSymbol = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, Color.BLACK);
			PlotCurveCharacterstics proposedDistFWSymbol = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, new Color(80, 80, 80));
			if (proposedRupVsJB_hw.size() > 0 && proposedRupVsJB_fw.size() > 0) {
				proposedRupVsJB_fw.setName("Proposed, FW");
				funcs.add(proposedRupVsJB_fw);
				chars.add(proposedDistFWSymbol);
				rupVsJBXYs.add(funcs.get(funcs.size()-1));
				rupVsJBChars.add(chars.get(chars.size()-1));
				proposedRupVsJB_hw.setName("Proposed, HW");
				funcs.add(proposedRupVsJB_hw);
				chars.add(proposedDistHWSymbol);
				rupVsJBXYs.add(funcs.get(funcs.size()-1));
				rupVsJBChars.add(chars.get(chars.size()-1));
			} else {
				proposedRupVsJB.setName("Proposed");
				funcs.add(proposedRupVsJB);
				chars.add(proposedDistHWSymbol);
				rupVsJBXYs.add(funcs.get(funcs.size()-1));
				rupVsJBChars.add(chars.get(chars.size()-1));
			}
			
			plot = new PlotSpec(funcs, chars, geomLabel+", "+oDF.format(distance)+" km", "Rjb (km)", "Rrup (km)");
			plot.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
			
			annX = distRange.getLowerBound() + 0.025*distRange.getLength();
			annY = distRange.getLowerBound() + 0.975*distRange.getLength();
			double annDeltaY = distRange.getLength()*0.04;
			
			List<String> distAnns = new ArrayList<>();
			distAnns.add("Length="+oDF.format(length)+", Ztor="+oDF.format(upperDepth)
					+", Dip="+oDF.format(mech.dip())+", DDW="+oDF.format(ptSurf.getAveWidth()));
			String distTexPrefix = texLabel+texDistPrefixes[d];
			texFW.write("% "+oDF.format(distance)+"km surface distances\n");
			texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMJBCenteredPercentile",
					percentileStr(calcPercentiles(centeredDists, nshmDistsFW, true, null)))+"\n");
			texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMJBUncenteredPercentile",
					percentileStr(calcPercentiles(uncenteredDists, nshmDistsFW, true, null)))+"\n");
			if (mech != FocalMech.STRIKE_SLIP) {
				distAnns.add("NSHM23 Rjb="+oneDF.format(nshmDistsFW.getDistanceJB())
						+", RrupFW="+oneDF.format(nshmDistsFW.getDistanceRup())
						+", RrupHW="+oneDF.format(nshmDistsHW.getDistanceRup()));
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMJB", oneDF.format(nshmDistsFW.getDistanceJB()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupFW", oneDF.format(nshmDistsFW.getDistanceRup()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupHW", oneDF.format(nshmDistsHW.getDistanceRup()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupFWCenteredPercentile",
						percentileStr(calcPercentiles(centeredDists, nshmDistsFW, false, false)))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupFWUncenteredPercentile",
						percentileStr(calcPercentiles(uncenteredDists, nshmDistsFW, false, false)))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupHWCenteredPercentile",
						percentileStr(calcPercentiles(centeredDists, nshmDistsHW, false, true)))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupHWUncenteredPercentile",
						percentileStr(calcPercentiles(uncenteredDists, nshmDistsHW, false, true)))+"\n");
			} else {
				distAnns.add("NSHM23 Rjb="+oneDF.format(nshmDistsFW.getDistanceJB())
						+", Rrup="+oneDF.format(nshmDistsFW.getDistanceRup()));
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMJB", oneDF.format(nshmDistsFW.getDistanceJB()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRup", oneDF.format(nshmDistsFW.getDistanceRup()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupCenteredPercentile",
						percentileStr(calcPercentiles(centeredDists, nshmDistsFW, false, null)))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRupUncenteredPercentile",
						percentileStr(calcPercentiles(uncenteredDists, nshmDistsFW, false, null)))+"\n");
			}
			
			MinMaxAveTracker jbCenteredRange = calcJBRange(centeredDists);
			MinMaxAveTracker rupCenteredRange = calcRupRange(centeredDists);
			MinMaxAveTracker jbUncenteredRange = calcJBRange(uncenteredDists);
			MinMaxAveTracker rupUncenteredRange = calcRupRange(uncenteredDists);
			
			rangeTexDefine(texFW, distTexPrefix+"CenteredJB", jbCenteredRange);
			rangeTexDefine(texFW, distTexPrefix+"CenteredRup", rupCenteredRange);
			rangeTexDefine(texFW, distTexPrefix+"UncenteredJB", jbUncenteredRange);
			rangeTexDefine(texFW, distTexPrefix+"UncenteredRup", rupUncenteredRange);
			if (mech == FocalMech.STRIKE_SLIP) {
				distAnns.add("Centered: avg. Rjb="+oneDF.format(jbCenteredRange.getAverage())+", Rrup="+oneDF.format(rupCenteredRange.getAverage()));
//				distAnns.add("  Rjb ∈ ["+oneDF.format(jbRange.getMin())+", "+oneDF.format(jbRange.getMax())+"]");
//				distAnns.add("  Rrup ∈ ["+oneDF.format(rupRange.getMin())+", "+oneDF.format(rupRange.getMax())+"]");
				distAnns.add("Uncentered: avg. Rjb="+oneDF.format(jbUncenteredRange.getAverage())+", Rrup="+oneDF.format(rupUncenteredRange.getAverage()));
				distAnns.add("  Rjb ∈ ["+oneDF.format(jbUncenteredRange.getMin())+", "+oneDF.format(jbUncenteredRange.getMax())+"]");
				distAnns.add("  Rrup ∈ ["+oneDF.format(rupUncenteredRange.getMin())+", "+oneDF.format(rupUncenteredRange.getMax())+"]");
			} else {
				for (boolean centered : new boolean[] {true, false}) {
					SurfaceDistances[] myDists = centered ? centeredDists : uncenteredDists;
					MinMaxAveTracker myHWJBRange = calcJBRange(myDists, true);
					MinMaxAveTracker myHWRupRange = calcRupRange(myDists, true);
					MinMaxAveTracker myFWJBRange = calcJBRange(myDists, false);
					MinMaxAveTracker myFWRupRange = calcRupRange(myDists, false);
					String myPrefix = distTexPrefix + (centered ? "Centered" : "Uncentered");
					rangeTexDefine(texFW, myPrefix+"HWJB", myHWJBRange);
					rangeTexDefine(texFW, myPrefix+"HWRup", myHWRupRange);
					rangeTexDefine(texFW, myPrefix+"FWJB", myFWJBRange);
					rangeTexDefine(texFW, myPrefix+"FWRup", myFWRupRange);
					double hwFract = calcHWFract(myDists);
					texFW.write(LaTeXUtils.defineValueCommand(myPrefix+"HWFract", twoDF.format(hwFract))+"\n");
					if (centered)
						distAnns.add("Centered:");
					else
						distAnns.add("Uncentered:");
					distAnns.add("  FW avg. Rjb="+oneDF.format(myFWJBRange.getAverage())+", Rrup="+oneDF.format(myFWRupRange.getAverage()));
					distAnns.add("  HW avg. Rjb="+oneDF.format(myHWJBRange.getAverage())+", Rrup="+oneDF.format(myHWRupRange.getAverage()));
					if (!centered) {
						distAnns.add("  Rjb ∈ ["+oneDF.format(jbUncenteredRange.getMin())+", "+oneDF.format(jbUncenteredRange.getMax())+"]");
						distAnns.add("  Rrup ∈ ["+oneDF.format(rupUncenteredRange.getMin())+", "+oneDF.format(rupUncenteredRange.getMax())+"]");
						distAnns.add("  HW fract: "+twoDF.format(hwFract));
					}
				}
				
			}
			texFW.write("\n");
			
			Font labelFont = new Font(Font.SANS_SERIF, Font.PLAIN, 20);
			for (String str : distAnns) {
				ann = new XYTextAnnotation(str, annX, annY);
				ann.setFont(labelFont);
				ann.setTextAnchor(TextAnchor.TOP_LEFT);
				plot.addPlotAnnotation(ann);
				annY -= annDeltaY;
			}
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plot, false, false, distRange, distRange);
			
			double tick;
			double histBin;
			if (distRange.getLength() > 40) {
				tick = 5;
				histBin = 1;
			} else if (distRange.getLength() > 15) {
				tick = 2;
				histBin = 0.5;
			} else {
				tick = 1;
				histBin = 0.25;
			}
			PlotUtils.setXTick(gp, tick);
			PlotUtils.setYTick(gp, tick);
			
//			PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix+"_"+oDF.format(distance)+"km_rrup_vs_jb",
//					gp, 800, 700, true, true, false);
			
			// now add Rjb hist below
			EvenlyDiscretizedFunc histBins = HistogramFunction.getEncompassingHistogram(0.001, maxDist, histBin);
			
			HistogramFunction centeredHist = calcJBHist(centeredDists, histBins);
			HistogramFunction uncenteredHist = calcJBHist(uncenteredDists, histBins);
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			funcs.add(uncenteredHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.LIGHT_GRAY));
			
			funcs.add(centeredHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_blue));
			
			Range histRange = new Range(0d, 1.05d);
			double lineUpper = 1d;
			
//			for (WeightedValue<SurfaceDistances> dist : proposedDists) {
//				DefaultXY_DataSet line = new DefaultXY_DataSet();
//				line.set(dist.value.getDistanceJB(), 0d);
//				line.set(dist.value.getDistanceJB(), histRange.getUpperBound());
//				funcs.add(line);
//				if (mech == FocalMech.STRIKE_SLIP || dist.value.getDistanceX() >= 0d)
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//				else
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, new Color(80, 80, 80)));
//			}
//			
//			DefaultXY_DataSet line = new DefaultXY_DataSet();
//			line.set(nshmDistsFW.getDistanceJB(), 0d);
//			line.set(nshmDistsFW.getDistanceJB(), lineUpper);
//			funcs.add(line);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_brown));
//			
//			double fwJB = calcJBRange(centeredDists, false).getAverage();
//			line = new DefaultXY_DataSet();
//			line.set(fwJB, 0d);
//			line.set(fwJB, lineUpper);
//			funcs.add(line);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_lightorange));
//			
//			double hwJB = calcJBRange(centeredDists, true).getAverage();
//			line = new DefaultXY_DataSet();
//			line.set(hwJB, 0d);
//			line.set(hwJB, lineUpper);
//			funcs.add(line);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_orange));
			
			for (int i=0; i<rupVsJBXYs.size(); i++) {
				XY_DataSet xy = rupVsJBXYs.get(i);
				PlotCurveCharacterstics pChar = rupVsJBChars.get(i);
				
				float thickness = pChar.getSymbol() == PlotSymbol.FILLED_CIRCLE ? 2f : avgThickness;
				
				for (Point2D pt : xy) {
					DefaultXY_DataSet line = new DefaultXY_DataSet();
					line.set(pt.getX(), 0d);
					line.set(pt.getX(), lineUpper);
					funcs.add(line);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, pChar.getColor()));
				}
			}
			
			for (int i=0; i<rupVsJBXYs.size(); i++) {
				XY_DataSet xy = rupVsJBXYs.get(i);
				PlotCurveCharacterstics pChar = rupVsJBChars.get(i);
				
				DefaultXY_DataSet modXY = new DefaultXY_DataSet();
				for (Point2D pt : xy)
					modXY.set(pt.getX(), lineUpper);
				
				funcs.add(modXY);
				chars.add(pChar);
			}
			
			PlotSpec histPlot = new PlotSpec(funcs, chars, null, plot.getXAxisLabel(), " ");
			
			gp.drawGraphPanel(List.of(plot, histPlot), false, false,
					List.of(distRange), List.of(distRange, histRange));
			
			PlotUtils.setSubPlotWeights(gp, 8, 2);
			
			CombinedDomainXYPlot combPlot = (CombinedDomainXYPlot)gp.getPlot();
			
			XYPlot subPlot = ((XYPlot)combPlot.getSubplots().get(1));
			subPlot.getRangeAxis().setTickLabelsVisible(false);
			subPlot.setRangeGridlinesVisible(false);
			
			PlotUtils.setXTick(gp, tick);
			PlotUtils.setYTick(gp, tick);
			
			PlotUtils.writePlots(outputDir, prefix+"_"+oDF.format(distance)+"km_rrup_vs_jb",
					gp, 800, 900, true, true, false);
			
			// mow median GM vs comparable
			double[] medians = new double[uncenteredExceedProbs.length];
			for (int i=0; i<medians.length; i++)
				medians[i] = uncenteredExceedProbs[i].getFirstInterpolatedX(0.5d);
			
			Boolean[] hwFlags;
			if (mech == FocalMech.STRIKE_SLIP || proposedRupVsJB_hw.size() == 0)
				hwFlags = new Boolean[] {null};
			else
				hwFlags = new Boolean[] {false, true};
			
			for (Boolean hwFlag : hwFlags) {
				if (hwFlag != null && hwFlag == false && distance == 0)
					continue;
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				double sortMin = Double.POSITIVE_INFINITY;
				double sortMax = 0d;
				List<double[]> sortScalarsList = new ArrayList<>();
				List<double[]> sortIMLsList = new ArrayList<>();
				WeightedList<FractileBin> fractiles = proposedCorr.getFractiles();
				List<double[]> sortBinLowerSortables = new ArrayList<>();
				List<double[]> sortBinUpperSortables = new ArrayList<>();
				double[] binLowerSortables = null;
				double[] binUpperSortables = null;
				
				WeightedList<SurfaceDistances> fractileDists;
				if (hwFlag == null) {
					fractileDists = proposedDists;
				} else {
					fractileDists = new WeightedList<>();
					for (WeightedValue<SurfaceDistances> val : proposedDists)
						if (hwFlag == val.value.getDistanceX() >= 0d)
							fractileDists.add(val);
				}
				
				for (Sortables sortable : sortables) {
					double[] sortValues, sortIMLs;
					SurfaceDistances[] sortDists;
					if (hwFlag == null) {
						sortDists = uncenteredDists;
						sortValues = new double[uncenteredDists.length];
						for (int i=0; i<sortValues.length; i++)
							sortValues[i] = sortable.function.apply(uncenteredDists[i]);
						sortIMLs = medians;
					} else {
						List<Double> sorts = new ArrayList<>();
						List<Double> imls = new ArrayList<>();
						List<SurfaceDistances> dists = new ArrayList<>();
						for (int i=0; i<uncenteredDists.length; i++) {
							if (hwFlag == uncenteredDists[i].getDistanceX() >= 0d) {
								dists.add(uncenteredDists[i]);
								sorts.add(sortable.function.apply(uncenteredDists[i]));
								imls.add(medians[i]);
							}
						}
						sortDists = dists.toArray(new SurfaceDistances[0]);
						sortValues = Doubles.toArray(sorts);
						sortIMLs = Doubles.toArray(imls);
					}
					sortScalarsList.add(sortValues);
					sortIMLsList.add(sortIMLs);
					sortMin = Math.min(sortMin, StatUtils.min(sortValues));
					sortMax = Math.max(sortMax, StatUtils.max(sortValues));
					DefaultXY_DataSet cloud = new DefaultXY_DataSet(sortValues, sortIMLs);
					
					DiscretizedFunc average = calcBinnedAverage(cloud, 50);
					
					average.setName(sortable.label);
					
					if (distance > 0d) {
						DefaultXY_DataSet polygon = calcEnvelopePoly(cloud, 50);
						funcs.add(polygon);
						chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, sortable.transColor));
					}
					
					funcs.add(average);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, sortable.color));
					
					double[] unityWeights = new double[sortValues.length];
					for (int i=0; i<unityWeights.length; i++)
						unityWeights[i] = 1d;
					LightFixedXFunc sortNCDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(sortValues, unityWeights);
					double[] myBinLowerSortables = new double[fractiles.size()];
					double[] myBinUpperSortables = new double[fractiles.size()];
					for (int f=0; f<fractiles.size(); f++) {
						FractileBin bin = fractiles.getValue(f);

						double lowerX = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(sortNCDF, bin.minimum);
						double upperX = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(sortNCDF, bin.maximum);
						myBinLowerSortables[f] = lowerX;
						myBinUpperSortables[f] = upperX;

//						double lowerY = f == 0 ? average.getY(0) : average.getInterpolatedY(lowerX);
//						double upperY = f == fractiles.size()-1 ? average.getY(average.size()-1) : average.getInterpolatedY(upperX);
//
//						if (lowerX > 0d) {
//							double tickUp = Math.pow(10, Math.log10(lowerY)+0.05);
//							double tickDown = Math.pow(10, Math.log10(lowerY)-0.05);
//							funcs.add(0, new DefaultXY_DataSet(lowerX, tickDown, lowerX, tickUp));
//							chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, sortable.color));
//						}
//
//						double tickUp = Math.pow(10, Math.log10(upperY)+0.05);
//						double tickDown = Math.pow(10, Math.log10(upperY)-0.05);
//						funcs.add(0, new DefaultXY_DataSet(upperX, tickDown, upperX, tickUp));
//						chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, sortable.color));
					}

					sortBinLowerSortables.add(myBinLowerSortables);
					sortBinUpperSortables.add(myBinUpperSortables);
					if (sortable == SORT_QUANTITY) {
						binLowerSortables = myBinLowerSortables;
						binUpperSortables = myBinUpperSortables;
					}

					// calculate and add manually
					int[] counts = new int[fractiles.size()];
					double[] sumJBs = new double[fractiles.size()];
					double[] sumRups = new double[fractiles.size()];
					double[] sumXs = new double[fractiles.size()];
					List<List<Double>> binIMLs = new ArrayList<>(fractiles.size());
					for (int f=0; f<fractiles.size(); f++)
						binIMLs.add(new ArrayList<>());

					for (int i=0; i<sortValues.length; i++) {
						for (int f=0; f<myBinLowerSortables.length; f++) {
							if (f == fractiles.size()-1 || sortValues[i] <= myBinUpperSortables[f]) {
								counts[f]++;
								sumJBs[f] += sortDists[i].getDistanceJB();
								sumRups[f] += sortDists[i].getDistanceRup();
								sumXs[f] += sortDists[i].getDistanceX();
								binIMLs.get(f).add(sortIMLs[i]);
								break;
							}
						}
					}

					DefaultXY_DataSet scatter = new DefaultXY_DataSet();
					DefaultXY_DataSet trueAvgScatter = new DefaultXY_DataSet();
					funcs.add(scatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, sortable.color));
					funcs.add(scatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, sortable.color.darker().darker()));
//					funcs.add(trueAvgScatter);
//					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, Color.BLACK));
					for (int f=0; f<fractiles.size(); f++) {
						if (counts[f] > 0) {
							double lower = myBinLowerSortables[f];
							double upper = myBinUpperSortables[f];
							double middle = 0.5*(lower + upper);
							double rJB = sumJBs[f]/(double)counts[f];
							double rRup = sumRups[f]/(double)counts[f];
							double rX = sumXs[f]/(double)counts[f];
							gmm0.setPropagationEffectParams(new SurfaceDistances.Precomputed(siteLoc, rRup, rJB, rX));
							scatter.set(middle, Math.exp(gmm0.getMean()));
							List<Double> myIMLs = binIMLs.get(f);
							DoubleSummaryStatistics imlStats = myIMLs.stream().mapToDouble(v->v).summaryStatistics();
							trueAvgScatter.set(middle, imlStats.getAverage());
							
//							funcs.add(new DefaultXY_DataSet(lower, imlStats.getMin(), lower, imlStats.getMax()));
//							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sortable.color));
//							funcs.add(new DefaultXY_DataSet(upper, imlStats.getMin(), upper, imlStats.getMax()));
//							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sortable.color));
//							funcs.add(new DefaultXY_DataSet(lower, imlStats.getAverage(), upper, imlStats.getAverage()));
//							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sortable.color));
						}
					}
				}
				
				Range medianIMLRange = medianIMLRanges[d];
				
				// add corr lines and bins
				PlotCurveCharacterstics binChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.BLACK);
//				XY_DataSet proposedGMs = new DefaultXY_DataSet();
				for (int f=0; f<fractiles.size(); f++) {
					double lower = binLowerSortables[f];
					double upper = binUpperSortables[f];
					if (f == 0) {
						DefaultXY_DataSet line = new DefaultXY_DataSet();
						line.set(lower, medianIMLRange.getLowerBound());
						line.set(lower, medianIMLRange.getUpperBound());
						funcs.add(line);
						chars.add(binChar);
					}
					
					DefaultXY_DataSet line = new DefaultXY_DataSet();
					line.set(upper, medianIMLRange.getLowerBound());
					line.set(upper, medianIMLRange.getUpperBound());
					funcs.add(line);
					chars.add(binChar);
					
//					gmm0.setPropagationEffectParams(fractileDists.getValue(f));
//					
//					double iml = Math.exp(gmm0.getMean());
					
//					line = new DefaultXY_DataSet();
//					line.set(lower, iml);
//					line.set(upper, iml);
//					if (f == 0)
//						line.setName("Proposed ");
//					funcs.add(line);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
//					
//					proposedGMs.set(0.5*(upper+lower), iml);
				}
////				proposedGMs.setName("Proposed");
//				funcs.add(proposedGMs);
//				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, Sortables.RRUP_PLUS_RJB.color));
////				proposedGMs = proposedGMs.deepClone();
////				proposedGMs.setName(null);
//				funcs.add(proposedGMs);
//				chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, Color.BLACK));
				
				if (sortMin < 5)
					sortMin = 0;
				else
					sortMin = Math.floor(sortMin/5d)*5d;
				sortMax = Math.ceil(sortMax/5d)*5d;
				Range sortRange = new Range(sortMin, sortMax);
				
				String title = imLabel+", "+oDF.format(distance)+" km";
				if (hwFlag != null) {
					if (hwFlag)
						title += ", Hanging Wall";
					else
						title += ", Footwall";
				}
				PlotSpec sortPlot = new PlotSpec(funcs, chars, title, "Sorting Quantity (ξ)", "Median "+imlAxisLabel);
				sortPlot.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				/*
				 * Variance plot
				 */
				
				double maxVar = 0d;
				DiscretizedFunc[] varFuncs = new DiscretizedFunc[sortables.length];
				for (int s=0; s<sortables.length; s++) {
					double[] scalars = sortScalarsList.get(s);
					double minScalar = StatUtils.min(scalars);
					double maxScalar = StatUtils.max(scalars);
					
					double myDelta = Math.min((maxScalar-minScalar)/50d, 1d);
					EvenlyDiscretizedFunc varBins = HistogramFunction.getEncompassingHistogram(minScalar+0.25*myDelta, maxScalar-0.25*myDelta, myDelta);
					
					List<Double> zeroValues = new ArrayList<>();
					
					double[] imls =sortIMLsList.get(s);
					List<List<Double>> binValues = new ArrayList<>(varBins.size());
					for (int i=0; i<varBins.size(); i++)
						binValues.add(new ArrayList<>());
					Preconditions.checkState(scalars.length == imls.length);
					for (int i=0; i<scalars.length; i++) {
						double value = scalars[i];
						double iml = imls[i];
						if (value == 0d) {
							zeroValues.add(iml);
						} else {
							int index = varBins.getClosestXIndex(value);
							binValues.get(index).add(iml);
						}
					}
					
					ArbitrarilyDiscretizedFunc varFunc = new ArbitrarilyDiscretizedFunc();
					
					for (int i=-1; i<varBins.size(); i++) {
						double x;
						double[] binIMLs;
						if (i < 0) {
							x = 0d;
							binIMLs = Doubles.toArray(zeroValues);
						} else {
							x = varBins.getX(i);
							binIMLs = Doubles.toArray(binValues.get(i));
						}
						if (binIMLs.length == 0)
							continue;
						double var;
						if (binIMLs.length == 1)
							var = 0d;
						else
							var = StatUtils.variance(binIMLs);
						double mean = StatUtils.mean(binIMLs);
						double sd = Math.sqrt(var);
						double cov = sd/mean;
						double plotVal = cov;
						if (i <= 0) {
							System.out.println(sortables[s]+" "+(float)distance+" km, hwFlag="+hwFlag+" x="+(float)x+" var="
									+(float)var+", sd="+(float)sd+", cov="+(float)cov);
//							System.out.print("\tValues: ");
//							for (int j=0; j<binIMLs.length && j<20; j++) {
//								if (j > 0)
//									System.out.print(", ");
//								System.out.print((float)binIMLs[j]);
//							}
//							System.out.println();
						}
						Preconditions.checkState(Double.isFinite(var), "Bad var=%s for %s", var, binIMLs);
						varFunc.set(x, plotVal);
						maxVar = Math.max(maxVar, plotVal);
					}
					if (minScalar > 0d)
						varFunc.set(minScalar, varFunc.getY(0));
					varFunc.set(maxScalar, varFunc.getY(varFunc.size()-1));
					
					System.out.print("\tValues: ");
					for (int i=0; i<varFunc.size(); i++) {
						if (i>0)
							System.out.print(", ");
						System.out.print((float)varFunc.getY(i));
					}
					System.out.println();
					
//					System.out.println(sortables[s]+" Variance Func:\n"+varFunc+"\n\n");
					
					varFuncs[s] = varFunc;
				}
				
				maxVar = Math.max(0.1, Math.ceil((maxVar*1.05)*20d)/20d);
				Range varRange = new Range(0d, maxVar);
				
				List<XY_DataSet> varPlotFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> varChars = new ArrayList<>();
				
				double varTick = varRange.getLength()/30d;
				for (int s=0; s<sortables.length; s++) {
					varPlotFuncs.add(varFuncs[s]);
					varChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, sortables[s].color));
					
//					double[] lowers = sortBinLowerSortables.get(s);
//					double[] uppers = sortBinUpperSortables.get(s);
//					
//					for (int f=0; f<lowers.length; f++) {
//						double y = f == 0 ? varFuncs[s].getY(0) : varFuncs[s].getInterpolatedY(lowers[f]);
//						varPlotFuncs.add(new DefaultXY_DataSet(lowers[f], y-varTick, lowers[f], y+varTick));
//						varChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sortables[s].color));
//					}
//					double x = uppers[uppers.length-1];
//					double y = varFuncs[s].getY(varFuncs[s].size()-1);
//					varPlotFuncs.add(new DefaultXY_DataSet(x, y-varTick, x, y+varTick));
//					varChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sortables[s].color));
				}
				
				PlotSpec varPlot = new PlotSpec(varPlotFuncs, varChars, null, sortPlot.getXAxisLabel(), "COV");
				
				/*
				 * Bin histogram plot
				 */
				
				EvenlyDiscretizedFunc sortHistBins = HistogramFunction.getEncompassingHistogram(
//						sortMin+1e-3, sortMax, sortMax > 50 ? 2d : 1d);
						sortMin+1e-3, sortMax, sortMax > 50 ? 1d : 0.5d);
				
				List<XY_DataSet> histFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> histChars = new ArrayList<>();
				
				for (int s=0; s<sortables.length; s++) {
					EvenlyDiscretizedFunc hist = new EvenlyDiscretizedFunc(sortHistBins.getMinX(), sortHistBins.getMaxX(), sortHistBins.size());
					double[] scalars = sortScalarsList.get(s);
					for (double value : scalars) {
						int index = hist.getClosestXIndex(value);
						hist.add(index, 1d);
					}
//					hist.scale(1d/hist.getMaxY());
					hist.scale(1d/hist.calcSumOfY_Vals());
					
					histFuncs.add(hist);
					histChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, sortables[s].transColor));
				}
				
				// mow rescale to max of any
				double maxHistY = 0d;
				for (XY_DataSet xy : histFuncs)
					maxHistY = Math.max(maxHistY, xy.getMaxY());
				for (XY_DataSet xy : histFuncs)
					((DiscretizedFunc)xy).scale(1d/maxHistY);
				
				PlotSpec sortHistPlot = new PlotSpec(histFuncs, histChars, null, sortPlot.getXAxisLabel(), "Density");
				
				// add fractile annotations
				double fLine0 = 1.2;
				double fLine1 = fLine0+0.2;
				
//				for (int s=0; s<sortables.length; s++) {
//					if (sortables[s] != SORT_QUANTITY) {
//						double[] lowers = sortBinLowerSortables.get(s);
//						double[] uppers = sortBinUpperSortables.get(s);
//						PlotCurveCharacterstics oBinChar = new PlotCurveCharacterstics(binChar.getLineType(), binChar.getLineWidth(), sortables[s].color);
//						histFuncs.add(new DefaultXY_DataSet(lowers[0], fLine0, lowers[0], fLine1));
//						histChars.add(oBinChar);
//						for (double upper : uppers) {
//							histFuncs.add(new DefaultXY_DataSet(upper, fLine0, upper, fLine1));
//							histChars.add(oBinChar);
//						}
//					}
//				}
				
				char sub1 = '₁';
//				Font weightFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
//				Font percentileFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
				Font weightFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
				Font percentileFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
				for (int f=0; f<fractiles.size(); f++) {
					FractileBin bin = fractiles.getValue(f);
					double weight = fractiles.getWeight(f);
					
					double lower = binLowerSortables[f];
					double upper = binUpperSortables[f];
					
					if (f == 0) {
						DefaultXY_DataSet line = new DefaultXY_DataSet();
						line.set(lower, fLine0);
						line.set(lower, fLine1);
						histFuncs.add(line);
						histChars.add(binChar);
						line = new DefaultXY_DataSet();
						line.set(lower, varRange.getLowerBound());
						line.set(lower, varRange.getUpperBound());
						varPlotFuncs.add(line);
						varChars.add(binChar);
						
						double lowerFractOfRange = (lower-sortRange.getLowerBound())/sortRange.getLength();
						boolean snap = lowerFractOfRange < 0.02;
						ann = new XYTextAnnotation(oDF.format(bin.minimum*100d)+"%", snap ? sortRange.getLowerBound() : lower, fLine0);
						ann.setFont(percentileFont);
						ann.setTextAnchor(snap ? TextAnchor.TOP_LEFT : TextAnchor.TOP_CENTER);
						sortHistPlot.addPlotAnnotation(ann);
					}
					
					histFuncs.add(new DefaultXY_DataSet(upper, fLine0, upper, fLine1));
					histChars.add(binChar);
					varPlotFuncs.add(new DefaultXY_DataSet(upper, varRange.getLowerBound(), upper, varRange.getUpperBound()));
					varChars.add(binChar);
					
					double upperFractOfRange = (upper-sortRange.getLowerBound())/sortRange.getLength();
					boolean snap = upperFractOfRange > 0.97;
					ann = new XYTextAnnotation(oDF.format(bin.maximum*100d)+"%", snap ? sortRange.getUpperBound() : upper, fLine0);
					ann.setFont(percentileFont);
					if (snap)
						ann.setTextAnchor(TextAnchor.TOP_RIGHT);
					else
						ann.setTextAnchor(TextAnchor.TOP_CENTER);
					sortHistPlot.addPlotAnnotation(ann);
					
					int digits = ((float)weight+"").length();
					double fractOfRange = (upper - lower)/sortRange.getLength();
//					boolean include = (digits == 3 && fractOfRange > 0.07) || fractOfRange > 0.09;
//					boolean includeShort = (digits == 3 && fractOfRange > 0.04) || fractOfRange > 0.03;
					boolean include = (digits == 3 && fractOfRange > 0.06) || fractOfRange > 0.07;
					boolean includeShort = (digits == 3 && fractOfRange > 0.04) || fractOfRange > 0.03;
					char subNum = (char)(sub1+f);
					String weightTxt = "w"+subNum+"="+(float)weight;
					if (!include && includeShort) {
						include = true;
						weightTxt = (float)weight+"";
					}
					if (include) {
						double centerX = 0.5*(upper + lower);
						double centerY = 0.5*(fLine0 + fLine1);
						ann = new XYTextAnnotation(weightTxt, centerX, centerY);
//						ann = new XYTextAnnotation(oDF.format(weight*100d)+"%", centerX, centerY);
						ann.setFont(weightFont);
						ann.setTextAnchor(TextAnchor.CENTER);
						ann.setBackgroundPaint(new Color(255, 255, 255, 80));
						sortHistPlot.addPlotAnnotation(ann);
					}
				}
				
//				gp.drawGraphPanel(sortPlot, false, true, sortRange, medianIMLRange);
//				
//				PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix+"_"+oDF.format(distance)+"km_sorting",
//						gp, 800, 700, true, true, false);
				
				gp.drawGraphPanel(List.of(sortPlot, sortHistPlot), List.of(false), List.of(true, false),
						List.of(sortRange), List.of(medianIMLRange, new Range(0d, fLine1)));
				
				PlotUtils.setSubPlotWeights(gp, 6, 2);
				
				combPlot = (CombinedDomainXYPlot)gp.getPlot();
				
				subPlot = ((XYPlot)combPlot.getSubplots().get(1));
				subPlot.getRangeAxis().setTickLabelsVisible(false);
				subPlot.setRangeGridlinesVisible(false);
				
				String sortPrefix = prefix+"_"+perPrefix+"_"+oDF.format(distance)+"km_sorting";
				if (hwFlag != null) {
					if (hwFlag)
						sortPrefix += "_hw";
					else
						sortPrefix += "_fw";
				}
				PlotUtils.writePlots(outputDir, sortPrefix, gp, 1000, 900, true, true, false);
				
				gp.drawGraphPanel(List.of(sortPlot, varPlot, sortHistPlot), List.of(false), List.of(true, false, false),
						List.of(sortRange), List.of(medianIMLRange, varRange, new Range(0d, fLine1)));
				
				PlotUtils.setSubPlotWeights(gp, 6, 2, 2);
				
				combPlot = (CombinedDomainXYPlot)gp.getPlot();
				
				subPlot = ((XYPlot)combPlot.getSubplots().get(2));
				subPlot.getRangeAxis().setTickLabelsVisible(false);
				subPlot.setRangeGridlinesVisible(false);
				
				sortPrefix = prefix+"_"+perPrefix+"_"+oDF.format(distance)+"km_sorting_var";
				if (hwFlag != null) {
					if (hwFlag)
						sortPrefix += "_hw";
					else
						sortPrefix += "_fw";
				}
				PlotUtils.writePlots(outputDir, sortPrefix, gp, 1000, 1000, true, true, false);
			}
			
			// map plot
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			double minY = -0.25*length;
			double maxY = Math.max(distance, length);
//			double maxX = Math.max(0.25*maxY, ptSurf.getAveHorizontalWidth()/2d);
			double maxX = 0.5*maxY;
			double buffer = Math.max(1d, 0.05*maxY);
			double halfBuffer = 0.5*buffer;
			Range xRange = new Range(-(maxX+halfBuffer), maxX+halfBuffer);
			Range yRange = new Range(minY-buffer, maxY+buffer);
			
			DefaultXY_DataSet gridLocXY = new DefaultXY_DataSet();
			gridLocXY.set(projectLoc(gridLoc, gridLoc));
			gridLocXY.setName("Grid source location");
			funcs.add(gridLocXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, Color.BLACK));
			
			DefaultXY_DataSet siteLocXY = new DefaultXY_DataSet();
			siteLocXY.set(projectLoc(gridLoc, siteLoc));
			siteLocXY.setName("Site location ("+oDF.format(distance)+" km away)");
			funcs.add(siteLocXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_INV_TRIANGLE, 6f, Colors.tab_red));
			
			// average Rjb annotations
			double centeredAvgJB = jbCenteredRange.getAverage();
			double uncenteredAvgJB = jbUncenteredRange.getAverage();
			
			List<XYAnnotation> anns = new ArrayList<>();
			if (centeredAvgJB > 10d) {
				double siteToGridAz = LocationUtils.azimuthRad(siteLoc, gridLoc);
				Point2D sitePt = siteLocXY.get(0);
				double jbDeltaX = 0.04*xRange.getLength();
				
				for (boolean centered : new boolean[] {true,false}) {
					double avgJB = centered ? centeredAvgJB : uncenteredAvgJB;
					Point2D jbPt = projectLoc(gridLoc, LocationUtils.location(siteLoc, siteToGridAz, avgJB));
					Preconditions.checkState((float)sitePt.getX() == (float)jbPt.getX(), "This assumes site and grid loc have same x");
					
					double x = sitePt.getX();
					double offsetX = x + jbDeltaX;
					if (!centered)
						offsetX += jbDeltaX;
					
					double annOffset = 0.2*jbDeltaX;
					annX = offsetX + annOffset;
					TextAnchor anchor = TextAnchor.BASELINE_CENTER;
					annY = 0.5*(sitePt.getY() + jbPt.getY());
//					TextAnchor anchor = TextAnchor.BASELINE_RIGHT;
//					annY = jbPt.getY();
					
//					String text = centered ? "Centered" : "Uncentered";
//					text += " avg. Rjb="+oDF.format(avgJB);
					String text = oDF.format(avgJB)+" km";
					
					XYTextAnnotation jbAnn = new XYTextAnnotation(text, annX, annY);
					jbAnn.setTextAnchor(anchor);
					jbAnn.setRotationAnchor(anchor);
					jbAnn.setRotationAngle(0.5*Math.PI);
					jbAnn.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 16));
					jbAnn.setBackgroundPaint(new Color(255, 255, 255, 60));
					anns.add(jbAnn);
					
					DefaultXY_DataSet jbXY = new DefaultXY_DataSet();
					jbXY.set(sitePt);
					jbXY.set(offsetX, sitePt.getY());
					jbXY.set(offsetX, jbPt.getY());
					jbXY.set(jbPt);
					funcs.add(jbXY);
					chars.add(new PlotCurveCharacterstics(centered ? PlotLineType.SOLID : PlotLineType.SHORT_DASHED, 1f, Color.BLACK));
				}
			}
			
			int exampleInsertIndex = funcs.size();

			PlotCurveCharacterstics centeredTraceChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_blue);
			if (mech == FocalMech.STRIKE_SLIP) {
				boolean first = true;
				for (RectangularSurface surf : centeredExampleSurfs) {

					DefaultXY_DataSet xy = new DefaultXY_DataSet();
					for (Point2D pt : projectLocList(gridLoc, surf.getUpperEdge(), false))
						xy.set(pt);
					if (first) {
						xy.setName("Example centered surfaces");
						first = false;
					}

					funcs.add(xy);
					chars.add(centeredTraceChar);
				}
			} else {
				boolean firstFW = true;
				boolean firstHW = true;
				PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(
						PlotLineType.SOLID, 1f, centeredTraceChar.getColor());
//						PlotLineType.SOLID, 1f, Color.DARK_GRAY);
				
				PlotCurveCharacterstics fwTraceChar = new PlotCurveCharacterstics(
						centeredTraceChar.getLineType(), centeredTraceChar.getLineWidth(), Colors.tab_lightblue);
				PlotCurveCharacterstics fwOutlineChar = new PlotCurveCharacterstics(
						outlineChar.getLineType(), outlineChar.getLineWidth(), fwTraceChar.getColor());
//						outlineChar.getLineType(), outlineChar.getLineWidth(), Color.DARK_GRAY);
				boolean[] hws = new boolean[centeredExampleSurfs.length];
				int numHW = 0;
				int numFW = 0;
				for (int s=0; s<centeredExampleSurfs.length; s++) {
					hws[s] = centeredExampleSurfs[s].getDistanceX(siteLoc) >= 0d;
					if (hws[s])
						numHW++;
					else
						numFW++;
				}
				// traces frist
				for (int s=0; s<centeredExampleSurfs.length; s++) {
					DefaultXY_DataSet traceXY = new DefaultXY_DataSet();
					for (Point2D pt : projectLocList(gridLoc, centeredExampleSurfs[s].getUpperEdge(), false))
						traceXY.set(pt);
					funcs.add(traceXY);
					if (hws[s]) {
						// hanging wall
						if (firstHW) {
							if (numHW > 1)
								traceXY.setName("Example centered HW surfaces");
							else
								traceXY.setName("Example centered HW surface");
						}
						firstHW = false;
						chars.add(centeredTraceChar);
					} else {
						// foot wall
						if (numFW > 1)
							traceXY.setName("Example centered FW surfaces");
						else
							traceXY.setName("Example centered FW surface");
						firstFW = false;
						chars.add(fwTraceChar);
					}
				}
				// now outlines
				for (int s=0; s<centeredExampleSurfs.length; s++) {
					DefaultXY_DataSet outlineXY = new DefaultXY_DataSet();
					for (Point2D pt : projectLocList(gridLoc, centeredExampleSurfs[s].getPerimeter(), true))
						outlineXY.set(pt);
					funcs.add(outlineXY);
					if (hws[s]) {
						// hanging wall
						firstHW = false;
						chars.add(outlineChar);
//						addTransOutlineFuncs(funcs, chars, gridLoc, mapPlotSurfs[s], outlineChar.getColor(), outlineChar.getLineWidth(), minOutlineAlpha);
					} else {
						// foot wall
						firstFW = false;
						chars.add(fwOutlineChar);
//						addTransOutlineFuncs(funcs, chars, gridLoc, mapPlotSurfs[s], fwOutlineChar.getColor(), fwOutlineChar.getLineWidth(), minOutlineAlpha);
					}
				}
			}
			
			// example uncentered
			DefaultXY_DataSet traceXY = new DefaultXY_DataSet();
			for (Point2D pt : projectLocList(gridLoc, uncenteredExampleSurf.getUpperEdge(), false))
				traceXY.set(pt);
			
			PlotCurveCharacterstics uncenteredExampleChar = new PlotCurveCharacterstics(
					PlotLineType.SOLID, centeredTraceChar.getLineWidth(), Color.DARK_GRAY);
			funcs.add(exampleInsertIndex, traceXY);
			chars.add(exampleInsertIndex, uncenteredExampleChar);
			
			// now add label here
			DefaultXY_DataSet fakeTrace = new DefaultXY_DataSet();
			fakeTrace.set(-1000d, -100d);
			fakeTrace.setName("Example uncentered surface");
			funcs.add(fakeTrace);
			chars.add(uncenteredExampleChar);
			
			if (mech != FocalMech.STRIKE_SLIP) {
				DefaultXY_DataSet outlineXY = new DefaultXY_DataSet();
				for (Point2D pt : projectLocList(gridLoc, uncenteredExampleSurf.getPerimeter(), true))
					outlineXY.set(pt);
				
				funcs.add(exampleInsertIndex, outlineXY);
				chars.add(exampleInsertIndex, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
			}
			
			// spin the centered surfaces
			Color centeredHWColor = overWhite(transColor(Colors.tab_blue, 80));
			Color centeredFWColor = overWhite(transColor(Colors.tab_lightblue, 80));
			if (mech == FocalMech.STRIKE_SLIP) {
				// simple, no hw/fw
				DefaultXY_DataSet xyCircle = new DefaultXY_DataSet();
				for (RectangularSurface surf : centeredMapOutlineSurfs)
					xyCircle.set(projectLoc(gridLoc, surf.getLastLocOnUpperEdge()));
				DefaultXY_DataSet fakeCircleXY = new DefaultXY_DataSet();
				fakeCircleXY.setName("Centered range");
				fakeCircleXY.set(-100d, -100d);
				funcs.add(fakeCircleXY);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, centeredHWColor));
				funcs.add(xyCircle);
				chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, centeredHWColor));
			} else {
				boolean[] hws = new boolean[centeredMapOutlineSurfs.length];
				for (int i=0; i<centeredMapOutlineSurfs.length; i++)
					hws[i] = centeredMapOutlineSurfs[i].getDistanceX(siteLoc) >= 0d;
				List<List<RectangularSurface>> hwBundles = new ArrayList<>();
				List<List<RectangularSurface>> fwBundles = new ArrayList<>();
				boolean curHW = false;
				List<RectangularSurface> curBundle = null;
				for (int i=0; i<centeredMapOutlineSurfs.length; i++) {
					if (curBundle == null || curHW != hws[i]) {
						if (curBundle != null) {
							// close out the previous one with this one to fill the gap
							curBundle.add(centeredMapOutlineSurfs[i]);
							if (curHW)
								hwBundles.add(curBundle);
							else
								fwBundles.add(curBundle);
						}
						curBundle = new ArrayList<>();
						curHW = hws[i];
					}
					curBundle.add(centeredMapOutlineSurfs[i]);
				}
				if (curHW == hws[0] && curBundle.size() < centeredMapOutlineSurfs.length) {
					// we can stitch the last one with the first one
					if (curHW)
						hwBundles.get(0).addAll(0, curBundle);
					else
						fwBundles.get(0).addAll(0, curBundle);
				} else {
					// still stitch to the first one to close the gap
					curBundle.add(centeredMapOutlineSurfs[0]);
					if (curHW) {
						hwBundles.add(curBundle);
					} else {
						fwBundles.add(curBundle);
					}
				}
				for (boolean hw : new boolean[] {true,false}) {
					List<List<RectangularSurface>> bundles = hw ? hwBundles : fwBundles;
					Color color = hw ? centeredHWColor : centeredFWColor;
					
					if (bundles.isEmpty())
						continue;
					
					for (int i=0; i<bundles.size(); i++) {
						DefaultXY_DataSet xyCircle = new DefaultXY_DataSet();
						for (RectangularSurface surf : bundles.get(i)) {
							// use middle because we're showing the strike range, not the upper edge range
							double horzWidth = surf.getAveHorizontalWidth();
							double dipDirRad = Math.toRadians(surf.getAveDipDirection());
							Location traceEndLocLoc = surf.getLastLocOnUpperEdge();
							Location centerEndLoc = LocationUtils.location(traceEndLocLoc, dipDirRad, 0.5*horzWidth);
							// but need to push it further out to match how far away the trace is
							double traceDist = LocationUtils.horzDistance(gridLoc, traceEndLocLoc);
							double centerDist = LocationUtils.horzDistance(gridLoc, centerEndLoc);
							if (traceDist > centerDist) {
								double centerAz = LocationUtils.azimuthRad(gridLoc, centerEndLoc);
								centerEndLoc = LocationUtils.location(centerEndLoc, centerAz, traceDist-centerDist);
							}
							xyCircle.set(projectLoc(gridLoc, centerEndLoc));
//							xyCircle.set(projectLoc(gridLoc, surf.getLastLocOnUpperEdge()));
						}
						xyCircle.set(projectLoc(gridLoc, gridLoc));
						xyCircle.set(xyCircle.get(0));
						if (i == 0) {
							DefaultXY_DataSet fakeCircleXY = new DefaultXY_DataSet();
							if (hw)
								fakeCircleXY.setName("Centered HW strike range");
							else
								fakeCircleXY.setName("Centered FW strike range");
							fakeCircleXY.set(-100d, -100d);
							funcs.add(fakeCircleXY);
							chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, color));
						}
						funcs.add(xyCircle);
						chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, color));
					}
				}
			}
			
			// now spin the most extreme uncentered surfaces
			DefaultXY_DataSet xyCircle = new DefaultXY_DataSet();
			for (RectangularSurface surf : extremeUncenteredSurfs)
				xyCircle.set(projectLoc(gridLoc, surf.getLastLocOnUpperEdge()));
			DefaultXY_DataSet fakeCircleXY = new DefaultXY_DataSet();
			fakeCircleXY.setName("Uncentered range");
			fakeCircleXY.set(-100d, -100d);
			funcs.add(fakeCircleXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, polyColor));
			funcs.add(xyCircle);
			chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, polyColor));
			
			plot = new PlotSpec(funcs, chars, geomLabel, "X (km)", "Y (km)");
			plot.setLegendInset(true);
			plot.setPlotAnnotations(anns);
			
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			
			gp.drawGraphPanel(plot, false, false, xRange, yRange);
			
			PlotUtils.setXTick(gp, 5d);
			PlotUtils.setYTick(gp, 5d);
			
			PlotUtils.writePlots(outputDir, prefix+"_"+oDF.format(distance)+"km_map", gp, 800, false, true, true, false);
		}
		
		plots.get(0).setLegendInset(RectangleAnchor.BOTTOM_LEFT);
//		plots.get(plots.size()-1).setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
		
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<plots.size(); i++)
			yRanges.add(new Range(0d, 1d));
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plots, true, false, List.of(imlRange), yRanges);
		
		PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix, gp, 1200, 200+300*plots.size(), true, true, false);
		
		texFW.close();
	}

	private static DecimalFormat oDF = new DecimalFormat("0.#");
	private static DecimalFormat oneDF = new DecimalFormat("0.0");
	private static DecimalFormat twoDF = new DecimalFormat("0.00");
	
	private static SurfaceDistances[] calcDists(Location siteLoc, RectangularSurface[] surfs) {
		List<CompletableFuture<SurfaceDistances>> futures = new ArrayList<>(surfs.length);
		for (RectangularSurface surf : surfs)
			futures.add(CompletableFuture.supplyAsync(
					() -> {return surf.calcDistances(siteLoc);}));
		
		SurfaceDistances[] dists = new SurfaceDistances[surfs.length];
		for (int i=0; i<surfs.length; i++)
			dists[i] = futures.get(i).join();
		return dists;
	}
	
	private static DefaultXY_DataSet calcRrupVsRjb(SurfaceDistances[] dists) {
		return calcRrupVsRjb(dists, null);
	}
	
	private static DefaultXY_DataSet calcRrupVsRjb(SurfaceDistances[] dists, Boolean hwFlag) {
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (SurfaceDistances dist : dists) {
			if (hwFlag != null && hwFlag != dist.getDistanceX()>=0d)
				continue;
			ret.set(dist.getDistanceJB(), dist.getDistanceRup());
		}
		return ret;
	}
	
	private static DefaultXY_DataSet reverse(DefaultXY_DataSet func) {
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (int i=func.size(); --i>=0;)
			ret.set(func.get(i));
		return ret;
	}
	
	/**
	 * Plot artifacts can occur if the functions are too dense; this thins out intermediate nearby points,
	 * preserving the first, last, and most extreme values.
	 * 
	 * @param func
	 * @param cartesianDistance minimum distance between points in the thinned out function
	 * @return thinned function
	 */
	private static DefaultXY_DataSet thin(DefaultXY_DataSet func, double cartesianDistance) {
		double minDistSq = cartesianDistance*cartesianDistance;
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		ret.set(func.get(0));
		double prevX = func.getX(0);
		double prevY = func.getY(0);
		double lastX = func.getX(func.size()-1);
		double lastY = func.getY(func.size()-1);
		double minX = func.getMinX();
		double minY = func.getMaxY();
		double maxX = func.getMinX();
		double maxY = func.getMaxY();
		int numForLastCheck = Integer.max(1, (func.size()/100));
		for (int i=1; i<func.size(); i++) {
			double x = func.getX(i);
			double y = func.getY(i);
			double prevDistSq = (x-prevX)*(x-prevX) + (y-prevY)*(y-prevY);
			boolean prevCheck = prevDistSq > minDistSq;
			double lastDistSq = (x-lastX)*(x-lastX) + (y-lastY)*(y-lastY);
			boolean lastCheck = func.size()-i > numForLastCheck || lastDistSq < minDistSq;
			boolean isExtreme = x == minX && y == minY || x == maxX && y == maxY;
			if (prevCheck && lastCheck || isExtreme) {
				prevX = x;
				prevY = y;
				ret.set(x, y);
			}
		}
		ret.set(lastX, lastY);
		return ret;
	}
	
	private static DiscretizedFunc[] calcExceedProbs(Site site, EqkRupture rup, Supplier<? extends ScalarIMR> gmmRef,
			DiscretizedFunc xVals, RectangularSurface[] surfs) {
		System.out.println("Calculating exceedance probabilities for "+surfs.length+" surfaces");
		double[] linearX = new double[xVals.size()];
		double[] logX = new double[xVals.size()];
		for (int i=0; i<logX.length; i++) {
			linearX[i] = xVals.getX(i);
			logX[i] = Math.log(linearX[i]);
		}
		
		ArrayDeque<ScalarIMR> gmms = new ArrayDeque<>();
		
		List<CompletableFuture<DiscretizedFunc>> futures = new ArrayList<>();
		for (RuptureSurface surf : surfs) {
			futures.add(CompletableFuture.supplyAsync(() -> {
				ScalarIMR gmm = null;
				synchronized (gmms) {
					if (!gmms.isEmpty())
						gmm = gmms.pop();
				}
				if (gmm == null) {
					gmm = gmmRef.get();
					gmm.setSite(site);
				}
				EqkRupture myRup = new EqkRupture(rup.getMag(), rup.getAveRake(), surf, null);
				gmm.setEqkRupture(myRup);
				LightFixedXFunc logFunc = new LightFixedXFunc(logX, new double[linearX.length]);
				gmm.getExceedProbabilities(logFunc);
				synchronized (gmms) {
					gmms.push(gmm);
				}
				return new LightFixedXFunc(linearX, logFunc.getYVals());
			}));
		}
		
		DiscretizedFunc[] ret = new DiscretizedFunc[surfs.length];
		for (int r=0; r<surfs.length; r++)
			ret[r] = futures.get(r).join();
		return ret;
	}
	
	private static HistogramFunction calcJBHist(SurfaceDistances[] dists, EvenlyDiscretizedFunc histBins) {
		HistogramFunction hist = new HistogramFunction(histBins.getMinX(), histBins.getMaxX(), histBins.size());
		for (SurfaceDistances dist : dists)
			hist.add(hist.getClosestXIndex(dist.getDistanceJB()), 1d);
		hist.scale(1d/hist.getMaxY());
		return hist;
	}
	
	private static MinMaxAveTracker calcJBRange(SurfaceDistances[] dists) {
		return calcJBRange(dists, null);
	}
	
	private static MinMaxAveTracker calcRupRange(SurfaceDistances[] dists) {
		return calcRupRange(dists, null);
	}
	
	private static MinMaxAveTracker calcJBRange(SurfaceDistances[] dists, Boolean hwFlag) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (SurfaceDistances dist : dists)
			if (hwFlag == null || hwFlag == dist.getDistanceX()>=0d)
				track.addValue(dist.getDistanceJB());
		return track;
	}
	
	private static MinMaxAveTracker calcRupRange(SurfaceDistances[] dists, Boolean hwFlag) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (SurfaceDistances dist : dists)
			if (hwFlag == null || hwFlag == dist.getDistanceX()>=0d)
				track.addValue(dist.getDistanceRup());
		return track;
	}
	
	private static double calcHWFract(SurfaceDistances[] dists) {
		int numHW = 0;
		for (SurfaceDistances dist : dists)
			if (dist.getDistanceX() >= 0)
				numHW++;
		return (double)numHW/(double)dists.length;
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
	
	private static DiscretizedFunc calcBinnedAverage(XY_DataSet scatter, int numBins) {
		if ((float)scatter.getMinX() == (float)scatter.getMaxX()) {
			double sum = 0d;
			for (Point2D pt : scatter)
				sum += pt.getY();
			ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
			ret.set(scatter.getMinX(), sum/(double)scatter.size());
			return ret;
		}
		EvenlyDiscretizedFunc binning = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), numBins);
		int[] binCounts = new int[binning.size()];
		double[] binSums = new double[binning.size()];
		
		for (Point2D pt : scatter) {
			int index = binning.getClosestXIndex(pt.getX());
			binCounts[index]++;
			binSums[index] += pt.getY();
		}
		
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<binCounts.length; i++) {
			if (binCounts[i] > 0)
				ret.set(binning.getX(i), binSums[i]/(double)binCounts[i]);
		}
		
		return ret;
	}
		
	
	private static DefaultXY_DataSet calcEnvelopePoly(XY_DataSet scatter, int numBins) {
		EvenlyDiscretizedFunc polyDiscr = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), numBins);
		double[] binMins = new double[polyDiscr.size()];
		for (int i=0; i<binMins.length; i++)
			binMins[i] = Double.POSITIVE_INFINITY;
		double[] binMaxs = new double[polyDiscr.size()];
		for (Point2D pt : scatter) {
			int index = polyDiscr.getClosestXIndex(pt.getX());
			binMins[index] = Math.min(binMins[index], pt.getY());
			binMaxs[index] = Math.max(binMaxs[index], pt.getY());
		}
		DefaultXY_DataSet poly = new DefaultXY_DataSet();
		for (int i=0; i<binMins.length; i++)
			if (Double.isFinite(binMins[i]))
				poly.set(polyDiscr.getX(i), binMins[i]);
		for (int i=binMaxs.length; --i>=0;)
			if (Double.isFinite(binMins[i]))
				poly.set(polyDiscr.getX(i), binMaxs[i]);
		poly.set(poly.get(0));
		return poly;
	}
	
	private static Point2D projectLoc(Location refLoc, Location loc) {
		if (LocationUtils.areSimilar(refLoc, loc))
			return new Point2D.Double(0d, 0d);
		Location refUp = new Location(refLoc.getLatitude()+1d, refLoc.getLongitude());
		Location refDown = new Location(refLoc.getLatitude()-1d, refLoc.getLongitude());
		double x = LocationUtils.distanceToLine(refDown, refUp, loc);
		
		Location refLeft = new Location(refLoc.getLatitude(), refLoc.getLongitude()-1d);
		Location refRight = new Location(refLoc.getLatitude(), refLoc.getLongitude()+1d);
		double y = LocationUtils.distanceToLine(refRight, refLeft, loc);
		return new Point2D.Double(x, y);
	}
	
	private static List<Point2D> projectLocList(Location refLoc, List<Location> locs, boolean connect) {
		List<Point2D> ret = new ArrayList<>(connect ? locs.size()+1 : locs.size());
		
		for (Location loc : locs)
			ret.add(projectLoc(refLoc, loc));
		if (connect)
			ret.add(ret.get(0));
		
		return ret;
	}
	
	private static void addTransOutlineFuncs(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars,
			Location refLoc, RectangularSurface surf, Color surfaceColor, float thickness, int minAlpha) {
		// draw forward and reverse edges
		int steps = 20;
		Color transColor = new Color(surfaceColor.getRed(), surfaceColor.getGreen(), surfaceColor.getBlue(), minAlpha);
		CPT transCPT = new CPT(0d,  (double)steps-1d, surfaceColor, transColor);
		
		double horzWidth = surf.getAveHorizontalWidth();
		double widthEach = horzWidth/steps;
		double dipDirRad = Math.toRadians(surf.getAveDipDirection());
		for (boolean forward : new boolean[] {false,true}) {
			Location top = forward ? surf.getLastLocOnUpperEdge() : surf.getFirstLocOnUpperEdge();
			for (int step=0; step<steps; step++) {
				Location bottom = LocationUtils.location(top, dipDirRad, widthEach);
				
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				xy.set(projectLoc(refLoc, top));
				xy.set(projectLoc(refLoc, bottom));
				
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, transCPT.getColor((double)step)));
				
				top = bottom;
			}
		}
		
		// now add bottom edge
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		xy.set(projectLoc(refLoc, LocationUtils.location(surf.getFirstLocOnUpperEdge(), dipDirRad, horzWidth)));
		xy.set(projectLoc(refLoc, LocationUtils.location(surf.getLastLocOnUpperEdge(), dipDirRad, horzWidth)));
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, transColor));
	}
	
	private static Color transColor(Color color, int alpha) {
		return new Color(color.getRed(), color.getGreen(), color.getBlue(), alpha);
	}
	
	public static Color overWhite(Color rgba) {
		float a = rgba.getAlpha() / 255f;

		int r = Math.round(rgba.getRed()   * a + 255f * (1f - a));
		int g = Math.round(rgba.getGreen() * a + 255f * (1f - a));
		int b = Math.round(rgba.getBlue()  * a + 255f * (1f - a));

		return new Color(r, g, b);
	}
	
	private static void rangeTexDefine(FileWriter texFW, String prefix, MinMaxAveTracker track) throws IOException {
		texFW.write(LaTeXUtils.defineValueCommand(prefix+"Avg", oneDF.format(track.getAverage()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(prefix+"Min", oneDF.format(track.getMin()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand(prefix+"Max", oneDF.format(track.getMax()))+"\n");
	}
	
	public static Path2D buildPolygon(XY_DataSet xy) {
		Path2D polygon = new Path2D.Double();
		boolean firstPoint = true;

		for (Point2D pt : xy) {
			if (firstPoint) {
				polygon.moveTo(pt.getX(), pt.getY());
				firstPoint = false;
			} else {
				polygon.lineTo(pt.getX(), pt.getY());
			}
		}

		// Ensure the polygon is closed
		if (xy.size() > 2) {
			polygon.closePath();
		}
		return polygon;
	}
	
	private static double calcPercentiles(SurfaceDistances[] dists, SurfaceDistances target, boolean rJB, Boolean hwFlag) {
		List<Double> values = new ArrayList<>();
		for (SurfaceDistances dist : dists) {
			if (hwFlag != null && hwFlag != dist.getDistanceX() >= 0)
				continue;
			values.add(rJB ? dist.getDistanceJB() : dist.getDistanceRup());
		}
		
		LightFixedXFunc nCDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(values, null);
		double val = rJB ? target.getDistanceJB() : target.getDistanceRup();
		if (val < nCDF.getMinX())
			return 0d;
		else if (val > nCDF.getMaxX())
			return 1d;
		else
			return nCDF.getInterpolatedY(val);
	}
	
	private static String percentileStr(double fractile) {
		return LaTeXUtils.numberAsOrdinal((int)Math.round(fractile*100d));
	}

}
