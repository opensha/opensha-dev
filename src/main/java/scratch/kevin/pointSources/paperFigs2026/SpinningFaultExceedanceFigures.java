package scratch.kevin.pointSources.paperFigs2026;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Function;
import java.util.function.Supplier;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Leonard2010_MagLengthRelationship;
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
		String prefix, label, texLabel;
		
		switch (mech) {
		case STRIKE_SLIP:
			mag = 7.05d;
			ml = Leonard2010_MagLengthRelationship.STRIKE_SLIP;
			prefix = "m7_ss";
			label = "M7.05, Strike-Slip";
			texLabel = "MSevenSS";
			break;
		case REVERSE:
			mag = 7.05d;
			ml = Leonard2010_MagLengthRelationship.DIP_SLIP;
			prefix = "m7_rev";
			label = "M7.05, Reverse";
			texLabel = "MSevenRev";
			break;

		default:
			throw new IllegalStateException("Unexpected mech: "+mech);
		}

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
		Range medianIMLRange = new Range(6e-2, 8e-1);

		double rake = mech.rake();
		double dip = mech.dip();
		double upperDepth = 1d;
		double lowerDepth = 14d;
		double length = ml.getMedianLength(mag);
		
		System.out.println("Length is: "+length);
		Supplier<ScalarIMR> gmmRef = () -> {
			ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
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
		List<Color> corrColors = new ArrayList<>();
		
		PointSourceDistanceCorrection nshm23Corr = PointSourceDistanceCorrections.NSHM_2013.get();
		String nshm23Name = "NSHM23 As Published";
		PointSourceDistanceCorrection averageCorr = PointSourceDistanceCorrections.AVERAGE_SPINNING_CENTERED.get();
		String averageName = "Average Centered";
		DistanceDistributionCorrection proposedCorr = (DistanceDistributionCorrection)
				PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST.get();
		
		corrs.add(nshm23Corr);
		corrLabels.add(nshm23Name);
		corrColors.add(Colors.tab_brown);
		
		corrs.add(averageCorr);
		corrLabels.add(averageName+" Correction");
		corrColors.add(Colors.tab_orange);
		
		corrs.add(PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST.get());
		corrLabels.add("Proposed Correction");
		corrColors.add(Colors.tab_purple);
		
		Location gridLoc = new Location(0d, 0d);
		
		EqkRupture rup = new EqkRupture(mag, rake, null, null);

		DecimalFormat oDF = new DecimalFormat("0.#");
		DecimalFormat oneDF = new DecimalFormat("0.0");
		DecimalFormat twoDF = new DecimalFormat("0.00");
		
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
		RectangularSurface[] centeredCalcSurfs = surfBuilder.buildRandRectSurfaces(numCenteredCalcSurfs);
		RectangularSurface[] plotSurfs = surfBuilder.buildRandRectSurfaces(numPlotSurfs);
		
		surfBuilder.sampleDASs().sampleHypocentralDepths();
		RectangularSurface[] uncenteredCalcSurfs = surfBuilder.buildRandRectSurfaces(numUncenteredCalcSurfs);
		
		ScalarIMR gmm0 = gmmRef.get();
		
		List<PlotSpec> plots = new ArrayList<>();
		for (int d=0; d<distances.length; d++) {
			double distance = distances[d];
			System.out.println("Calculating for "+distance+" km");
			Location siteLoc = LocationUtils.location(gridLoc, 0d, distance);
			
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
			fakeMinMax.setName("Uncentered Range");
			Color polyColor = new Color(0, 0, 0, 40);
			funcs.add(fakeMinMax);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 4f, polyColor));
			funcs.add(minMax);
			chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, polyColor));
			
//			if (mech == FocalMech.STRIKE_SLIP) {
				for (int i=0; i<plotSurfs.length; i++) {
					DiscretizedFunc func = plotExceedProbs[i];
					if (i == 0)
						func.setName("Centered Individual");
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
				
				WeightedList<SurfaceDistances> corrDists = corr.getCorrectedDistances(
						siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance);
				System.out.println("Distances for "+corr);
				for (WeightedValue<SurfaceDistances> corrDist : corrDists)
					System.out.println("\tweight="+oDF.format(corrDist.weight)+"; "+corrDist.value);
				
				// calculate distance corrected
				rup.setRuptureSurface(ptSurf.getDistancedProtected(corr, TectonicRegionType.ACTIVE_SHALLOW, mag));
				RuptureExceedProbCalculator.calcExceedanceProbabilities(gmm0, rup, logXVals);
				DiscretizedFunc corrExceedProbs = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<logXVals.size(); i++)
					corrExceedProbs.set(xVals.getX(i), logXVals.getY(i));
				corrExceedProbs.setName(corrLabels.get(c));
				funcs.add(corrExceedProbs);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, avgThickness, corrColors.get(c)));
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, label, imlAxisLabel, "Exceedance Probability");
			
			double annX = xVals.getX((int)(xVals.size()*0.97));
			double annY = 0.97;
			Font distFont = new Font(Font.SANS_SERIF, Font.BOLD, 30);
			XYTextAnnotation ann = new XYTextAnnotation(oDF.format(distance)+" km", annX, annY);
			ann.setFont(distFont);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			plot.addPlotAnnotation(ann);

			SurfaceDistances nshmDists = nshm23Corr.getCorrectedDistances(
					siteLoc, ptSurf, TectonicRegionType.ACTIVE_SHALLOW, mag, distance).getValue(0);
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
			
			if (mech != FocalMech.STRIKE_SLIP) {
				centeredRupVsJB_fw = thin(centeredRupVsJB_fw, 0.1);
				centeredRupVsJB_fw.setName("Centered, FW");
				funcs.add(centeredRupVsJB_fw);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_lightblue));
				centeredRupVsJB_hw = thin(reverse(centeredRupVsJB_hw), 0.1);
				centeredRupVsJB_hw.setName("Centered, HW");
				funcs.add(centeredRupVsJB_hw);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_blue));
			} else {
				centeredRupVsJB.setName("Centered");
				funcs.add(centeredRupVsJB);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_blue));
			}
			
			DefaultXY_DataSet nshm23RupVsJB = new DefaultXY_DataSet();
			nshm23RupVsJB.set(nshmDists.getDistanceJB(), nshmDists.getDistanceRup());
			nshm23RupVsJB.setName(nshm23Name);
			funcs.add(nshm23RupVsJB);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 4f, Colors.tab_brown));
			
			for (WeightedValue<SurfaceDistances> dist : averageDists) {
				DefaultXY_DataSet avgRupVsJB = new DefaultXY_DataSet();
				avgRupVsJB.set(dist.value.getDistanceJB(), dist.value.getDistanceRup());
				Color color;
				if (averageDists.size() == 1) {
					avgRupVsJB.setName(averageName);
					color = Colors.tab_orange;
				} else if (dist.value.getDistanceX() >= 0d) {
					avgRupVsJB.setName(averageName+", HW");
					color = Colors.tab_orange;
				} else {
					avgRupVsJB.setName(averageName+", FW");
					color = Colors.tab_lightorange;
				}
				funcs.add(avgRupVsJB);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 4f, color));
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
			if (proposedRupVsJB_hw.size() > 0 && proposedRupVsJB_fw.size() > 0) {
				proposedRupVsJB_fw.setName("Proposed, FW");
				funcs.add(proposedRupVsJB_fw);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 4f, new Color(80, 80, 80)));
				proposedRupVsJB_hw.setName("Proposed, HW");
				funcs.add(proposedRupVsJB_hw);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 4f, Color.BLACK));
			} else {
				proposedRupVsJB.setName("Proposed");
				funcs.add(proposedRupVsJB);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 4f, Color.BLACK));
			}
			
			plot = new PlotSpec(funcs, chars, label+", "+oDF.format(distance)+" km", "Rjb (km)", "Rrup (km)");
			plot.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
			
			annX = distRange.getLowerBound() + 0.025*distRange.getLength();
			annY = distRange.getLowerBound() + 0.975*distRange.getLength();
			double annDeltaY = distRange.getLength()*0.04;
			
			List<String> distAnns = new ArrayList<>();
			distAnns.add("Length="+oDF.format(length)+", Ztor="+oDF.format(upperDepth)
					+", Dip="+oDF.format(mech.dip())+", DDW="+oDF.format(ptSurf.getAveWidth()));
			distAnns.add("NSHM23 Rjb="+oneDF.format(nshmDists.getDistanceJB())
					+", Rrup="+oneDF.format(nshmDists.getDistanceRup()));
			
			String distTexPrefix = texLabel+texDistPrefixes[d];
			texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMJB", oneDF.format(nshmDists.getDistanceJB()))+"\n");
			texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"NSHMRup", oneDF.format(nshmDists.getDistanceRup()))+"\n");
			
			if (mech == FocalMech.STRIKE_SLIP) {
				MinMaxAveTracker jbRange = calcJBRange(centeredDists);
				MinMaxAveTracker rupRange = calcRupRange(centeredDists);
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"CenteredJB", oneDF.format(jbRange.getAverage()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"CenteredRup", oneDF.format(rupRange.getAverage()))+"\n");
				distAnns.add("Centered: avg. Rjb="+oneDF.format(jbRange.getAverage())+", Rrup="+oneDF.format(rupRange.getAverage()));
//				distAnns.add("  Rjb ∈ ["+oneDF.format(jbRange.getMin())+", "+oneDF.format(jbRange.getMax())+"]");
//				distAnns.add("  Rrup ∈ ["+oneDF.format(rupRange.getMin())+", "+oneDF.format(rupRange.getMax())+"]");
				jbRange = calcJBRange(uncenteredDists);
				rupRange = calcRupRange(uncenteredDists);
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredJB", oneDF.format(jbRange.getAverage()))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredRup", oneDF.format(rupRange.getAverage()))+"\n");
				distAnns.add("Uncentered: avg. Rjb="+oneDF.format(jbRange.getAverage())+", Rrup="+oneDF.format(rupRange.getAverage()));
				distAnns.add("  Rjb ∈ ["+oneDF.format(jbRange.getMin())+", "+oneDF.format(jbRange.getMax())+"]");
				distAnns.add("  Rrup ∈ ["+oneDF.format(rupRange.getMin())+", "+oneDF.format(rupRange.getMax())+"]");
			} else {
				MinMaxAveTracker jbRange = calcJBRange(centeredDists);
				MinMaxAveTracker rupRange = calcRupRange(centeredDists);
				double hwJB = calcJBRange(centeredDists, true).getAverage();
				double hwRup = calcRupRange(centeredDists, true).getAverage();
				double fwJB = calcJBRange(centeredDists, false).getAverage();
				double fwRup = calcRupRange(centeredDists, false).getAverage();
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"CenteredFWJB", oneDF.format(fwJB))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"CenteredFWRup", oneDF.format(fwRup))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"CenteredHWJB", oneDF.format(hwJB))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"CenteredHWRup", oneDF.format(hwRup))+"\n");
				distAnns.add("Centered:");
				distAnns.add("  FW avg. Rjb="+oneDF.format(fwJB)+", Rrup="+oneDF.format(fwRup));
				distAnns.add("  HW avg. Rjb="+oneDF.format(hwJB)+", Rrup="+oneDF.format(hwRup));
//				distAnns.add("  Rjb ∈ ["+oneDF.format(jbRange.getMin())+", "+oneDF.format(jbRange.getMax())+"]");
//				distAnns.add("  Rrup ∈ ["+oneDF.format(rupRange.getMin())+", "+oneDF.format(rupRange.getMax())+"]");
				jbRange = calcJBRange(uncenteredDists);
				rupRange = calcRupRange(uncenteredDists);
				hwJB = calcJBRange(uncenteredDists, true).getAverage();
				hwRup = calcRupRange(uncenteredDists, true).getAverage();
				fwJB = calcJBRange(uncenteredDists, false).getAverage();
				fwRup = calcRupRange(uncenteredDists, false).getAverage();
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredFWJB", oneDF.format(fwJB))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredFWRup", oneDF.format(fwRup))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredHWJB", oneDF.format(hwJB))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredHWRup", oneDF.format(hwRup))+"\n");
				distAnns.add("Uncentered:");
				distAnns.add("  FW avg. Rjb="+oneDF.format(fwJB)+", Rrup="+oneDF.format(fwRup));
				distAnns.add("  HW avg. Rjb="+oneDF.format(hwJB)+", Rrup="+oneDF.format(hwRup));
				distAnns.add("  Rjb ∈ ["+oneDF.format(jbRange.getMin())+", "+oneDF.format(jbRange.getMax())+"]");
				distAnns.add("  Rrup ∈ ["+oneDF.format(rupRange.getMin())+", "+oneDF.format(rupRange.getMax())+"]");
				double hwFract = calcHWFract(uncenteredDists);
				distAnns.add("  HW fract: "+twoDF.format(hwFract));
				texFW.write(LaTeXUtils.defineValueCommand(distTexPrefix+"UncenteredHWFract", twoDF.format(hwFract))+"\n");
			}
			
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
			
			for (WeightedValue<SurfaceDistances> dist : proposedDists) {
				DefaultXY_DataSet line = new DefaultXY_DataSet();
				line.set(dist.value.getDistanceJB(), histRange.getLowerBound());
				line.set(dist.value.getDistanceJB(), histRange.getUpperBound());
				funcs.add(line);
				if (mech == FocalMech.STRIKE_SLIP || dist.value.getDistanceX() >= 0d)
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				else
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, new Color(80, 80, 80)));
			}
			
			DefaultXY_DataSet line = new DefaultXY_DataSet();
			line.set(nshmDists.getDistanceJB(), histRange.getLowerBound());
			line.set(nshmDists.getDistanceJB(), histRange.getUpperBound());
			funcs.add(line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_brown));
			
			double fwJB = calcJBRange(centeredDists, false).getAverage();
			line = new DefaultXY_DataSet();
			line.set(fwJB, histRange.getLowerBound());
			line.set(fwJB, histRange.getUpperBound());
			funcs.add(line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_lightorange));
			
			double hwJB = calcJBRange(centeredDists, true).getAverage();
			line = new DefaultXY_DataSet();
			line.set(hwJB, histRange.getLowerBound());
			line.set(hwJB, histRange.getUpperBound());
			funcs.add(line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, avgThickness, Colors.tab_orange));
			
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
			
			PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix+"_"+oDF.format(distance)+"km_rrup_vs_jb",
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
				Sortables[] sortables = Sortables.values();
				WeightedList<FractileBin> fractiles = proposedCorr.getFractiles();
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
					
					if (sortable !=Sortables.RRUP_PLUS_RJB) {
						double[] unityWeights = new double[sortValues.length];
						for (int i=0; i<unityWeights.length; i++)
							unityWeights[i] = 1d;
						LightFixedXFunc sortNCDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(sortValues, unityWeights);
						double[] myBinLowerSortables = new double[fractiles.size()];
						double[] myBinUpperSortables = new double[fractiles.size()];
						for (int f=0; f<fractiles.size(); f++) {
							FractileBin bin = fractiles.getValue(f);
							
							myBinLowerSortables[f] = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(sortNCDF, bin.minimum);
							myBinUpperSortables[f] = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(sortNCDF, bin.maximum);
						}
						
						// calculate and add manually
						int[] counts = new int[fractiles.size()];
						double[] sumJBs = new double[fractiles.size()];
						double[] sumRups = new double[fractiles.size()];
						double[] sumXs = new double[fractiles.size()];
						
						for (int i=0; i<sortValues.length; i++) {
							for (int f=0; f<myBinLowerSortables.length; f++) {
								if (sortValues[i] <= myBinUpperSortables[f]) {
									counts[f]++;
									sumJBs[f] += sortDists[i].getDistanceJB();
									sumRups[f] += sortDists[i].getDistanceRup();
									sumXs[f] += sortDists[i].getDistanceX();
									break;
								}
							}
						}
						
						DefaultXY_DataSet scatter = new DefaultXY_DataSet();
						for (int f=0; f<fractiles.size(); f++) {
							double middle = 0.5*(myBinLowerSortables[f] + myBinUpperSortables[f]);
							double rJB = sumJBs[f]/(double)counts[f];
							double rRup = sumRups[f]/(double)counts[f];
							double rX = sumXs[f]/(double)counts[f];
							gmm0.setPropagationEffectParams(new SurfaceDistances.Precomputed(siteLoc, rRup, rJB, rX));
							scatter.set(middle, Math.exp(gmm0.getMean()));
						}
						funcs.add(scatter);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, sortable.color));
						funcs.add(scatter);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, sortable.color.darker().darker()));
					}
				}
				
				WeightedList<SurfaceDistances> fractileDists;
				if (hwFlag == null) {
					fractileDists = proposedDists;
				} else {
					fractileDists = new WeightedList<>();
					for (WeightedValue<SurfaceDistances> val : proposedDists)
						if (hwFlag == val.value.getDistanceX() >= 0d)
							fractileDists.add(val);
				}
				Preconditions.checkState(fractiles.size() == fractileDists.size(), "%s != %s", fractiles.size(), fractileDists.size());
				double[] finalSortables = sortScalarsList.get(sortScalarsList.size()-1);
				double[] unityWeights = new double[finalSortables.length];
				for (int i=0; i<unityWeights.length; i++)
					unityWeights[i] = 1d;
				LightFixedXFunc sortNCDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(finalSortables, unityWeights);
				double[] binLowerSortables = new double[fractiles.size()];
				double[] binUpperSortables = new double[fractiles.size()];
				for (int f=0; f<fractiles.size(); f++) {
					FractileBin bin = fractiles.getValue(f);
					
					binLowerSortables[f] = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(sortNCDF, bin.minimum);
					binUpperSortables[f] = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(sortNCDF, bin.maximum);
				}
				
				// add corr lines and bins
				PlotCurveCharacterstics binChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.BLACK);
				XY_DataSet proposedGMs = new DefaultXY_DataSet();
				for (int f=0; f<fractiles.size(); f++) {
					double lower = binLowerSortables[f];
					double upper = binUpperSortables[f];
					if (f == 0) {
						line = new DefaultXY_DataSet();
						line.set(lower, medianIMLRange.getLowerBound());
						line.set(lower, medianIMLRange.getUpperBound());
						funcs.add(line);
						chars.add(binChar);
					}
					
					line = new DefaultXY_DataSet();
					line.set(upper, medianIMLRange.getLowerBound());
					line.set(upper, medianIMLRange.getUpperBound());
					funcs.add(line);
					chars.add(binChar);
					
					gmm0.setPropagationEffectParams(fractileDists.getValue(f));
					
					double iml = Math.exp(gmm0.getMean());
					
//					line = new DefaultXY_DataSet();
//					line.set(lower, iml);
//					line.set(upper, iml);
//					if (f == 0)
//						line.setName("Proposed ");
//					funcs.add(line);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
					
					proposedGMs.set(0.5*(upper+lower), iml);
				}
//				proposedGMs.setName("Proposed");
				funcs.add(proposedGMs);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, Sortables.RRUP_PLUS_RJB.color));
//				proposedGMs = proposedGMs.deepClone();
//				proposedGMs.setName(null);
				funcs.add(proposedGMs);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, Color.BLACK));
				
				if (sortMin < 5)
					sortMin = 0;
				else
					sortMin = Math.floor(sortMin/5d)*5d;
				sortMax = Math.ceil(sortMax/5d)*5d;
				Range sortRange = new Range(sortMin, sortMax);
				
				String title = plot.getTitle();
				if (hwFlag != null) {
					if (hwFlag)
						title += ", Hanging-Wall";
					else
						title += ", Foot-Wall";
				}
				PlotSpec sortPlot = new PlotSpec(funcs, chars, title, "Sorting Quantity", "Median "+imlAxisLabel);
				sortPlot.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				EvenlyDiscretizedFunc sortHistBins = HistogramFunction.getEncompassingHistogram(sortMin+1e-3, sortMax, sortMax > 50 ? 2d : 1d);
				
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				for (int s=0; s<sortables.length; s++) {
					HistogramFunction hist = new HistogramFunction(sortHistBins.getMinX(), sortHistBins.getMaxX(), sortHistBins.size());
					double[] scalars = sortScalarsList.get(s);
					for (double value : scalars) {
						int index = hist.getClosestXIndex(value);
						hist.add(index, 1d);
					}
					hist.scale(1d/hist.getMaxY());
					
					funcs.add(hist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, sortables[s].transColor));
				}
				
				PlotSpec sortHistPlot = new PlotSpec(funcs, chars, null, sortPlot.getXAxisLabel(), " ");
				
				// add fractile annotations
				double fLine0 = 1.2;
				double fLine1 = fLine0+0.2;
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
						line = new DefaultXY_DataSet();
						line.set(lower, fLine0);
						line.set(lower, fLine1);
						funcs.add(line);
						chars.add(binChar);
						
						double lowerFractOfRange = (lower-sortRange.getLowerBound())/sortRange.getLength();
						boolean snap = lowerFractOfRange < 0.02;
						ann = new XYTextAnnotation(oDF.format(bin.minimum*100d)+"%", snap ? sortRange.getLowerBound() : lower, fLine0);
						ann.setFont(percentileFont);
						ann.setTextAnchor(snap ? TextAnchor.TOP_LEFT : TextAnchor.TOP_CENTER);
						sortHistPlot.addPlotAnnotation(ann);
					}
					
					line = new DefaultXY_DataSet();
					line.set(upper, fLine0);
					line.set(upper, fLine1);
					funcs.add(line);
					chars.add(binChar);
					
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
					boolean include = (digits == 3 && fractOfRange > 0.07) || fractOfRange > 0.09;
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
						sortHistPlot.addPlotAnnotation(ann);
					}
				}
				
//				gp.drawGraphPanel(sortPlot, false, true, sortRange, medianIMLRange);
//				
//				PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix+"_"+oDF.format(distance)+"km_sorting",
//						gp, 800, 700, true, true, false);
				
				gp.drawGraphPanel(List.of(sortPlot, sortHistPlot), List.of(false), List.of(true, false),
						List.of(sortRange), List.of(medianIMLRange, new Range(0d, fLine1)));
				
				PlotUtils.setSubPlotWeights(gp, 8, 2);
				
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
				PlotUtils.writePlots(outputDir, sortPrefix, gp, 800, 900, true, true, false);
			}
		}
		
		plots.get(0).setLegendInset(RectangleAnchor.BOTTOM_LEFT);
//		plots.get(plots.size()-1).setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
		
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<plots.size(); i++)
			yRanges.add(new Range(0d, 1d));
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plots, true, false, List.of(imlRange), yRanges);
		
		PlotUtils.writePlots(outputDir, prefix+"_"+perPrefix, gp, 1000, 200+300*plots.size(), true, true, false);
		
		texFW.close();
	}
	
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
	
	private static DiscretizedFunc[] calcExceedProbs(Site site, EqkRupture rup, Supplier<ScalarIMR> gmmRef,
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

}
