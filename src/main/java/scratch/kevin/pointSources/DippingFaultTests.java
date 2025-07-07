package scratch.kevin.pointSources;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.WeightedValue;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
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
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.earthquake.util.GriddedFiniteRuptureSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.PointSurface.DistanceCorrected;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PointSourceDistanceCorrection;
import org.opensha.sha.faultSurface.utils.PointSourceDistanceCorrections;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.faultSurface.utils.RjbDistributionDistanceCorrection;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class DippingFaultTests {

	public static void main(String[] args) throws IOException {
		double mag = 6.55;
		FocalMech mech = FocalMech.REVERSE;
		double[] dists = {0d, 1d, 3d, 5d, 10, 15d, 20d, 24d, 50d};
//		double[] dists = {5d};
//		double[] dists = {15d};
//		double[] dists = {50d};
		double[] periods = {0d, 0.2, 1d, 5d};
		
		double testJB = 41;
		
		boolean printPointRrups = dists.length == 1;
		
//		Boolean[] forceHWs = {null, false, true};
//		Boolean[] forceHWs = {true};
		Boolean[] forceHWs = {null};
		
		Location loc = new Location(0d, 0d);
		
		AttenRelRef gmmRef = AttenRelRef.USGS_NSHM23_ACTIVE;
		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
		
		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.FIVE_POINT_RJB_DIST_ALONG.get();
		boolean randDD = true;
		boolean randDAS = true;
		
//		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.FIVE_POINT_RJB_DIST.get();
////		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.TWENTY_POINT_RJB_DIST.get();
//		boolean randDD = false;
//		boolean randDAS = false;
		
		for (Boolean forceHW : forceHWs) {
			
			int numFinite = 1000;
			if (forceHW != null || randDD || randDAS)
				numFinite *= 5;
			
			GriddedFiniteRuptureSettings finiteSettings = GriddedFiniteRuptureSettings.DEFAULT_CROSSHAIR.forStrike(0d)
					.forNumSurfaces(numFinite).forSampleAlongStrike(randDAS).forSampleDownDip(randDD);
			boolean trackHWRanges = !randDAS && !randDD;
			
			File outputDir = new File("/tmp/dipping_fault_tests");
			if (forceHW != null) {
				if (forceHW)
					outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_hw_only");
				else
					outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_fw_only");
			}
			Preconditions.checkState(outputDir.exists() || outputDir.mkdirs(),
					"Failed to create output directory: " + outputDir.getAbsolutePath());
			
			PointSurfaceBuilder builder = new PointSurfaceBuilder(loc)
					.random(12345l)
					.magnitude(mag);
			
			double dip = mech.dip();
			double dipRad = Math.toRadians(dip);
			double zTop = PointSourceNshm.SURF_BUILDER_DEFAULT.depthForMag(mag);
			double widthDD = PointSourceNshm.SURF_BUILDER_DEFAULT.calcWidth(mag, zTop, dipRad);
			double length = PointSourceNshm.SURF_BUILDER_DEFAULT.calcLength(mag);
			builder.upperDepthWidthAndDip(zTop, widthDD, dip);
			builder.fractionalDAS(0.5).fractionalHypocentralDepth(0.5d);
			builder.length(length);
			
			DecimalFormat oDF = new DecimalFormat("0.##");
			DecimalFormat pDF = new DecimalFormat("0.#%");
			
			PointSurface pointSurf = builder.buildPointSurface();
			WeightedList<? extends RuptureSurface> finiteSurfs = builder.build(BackgroundRupType.FINITE, finiteSettings);
			
			System.out.println("Fault dimensions: length= "+oDF.format(length)+", zTop="+oDF.format(zTop)+", dip="+oDF.format(dip)
					+", ddw="+oDF.format(widthDD)+", horzWidth="+oDF.format(widthDD*Math.cos(dipRad)));
			
			System.out.println("Built 1 point surfaces and "+finiteSurfs.size()+" finite surfaces");
			
			ScalarIMR gmm = gmmRef.get();
			gmm.setIntensityMeasure(PGA_Param.NAME);
			
			Range imRange = new Range(-5d, 2d);
			Range distRange = new Range(0d, StatUtils.max(dists)+10d);
			EvenlyDiscretizedFunc logIMfunc = new EvenlyDiscretizedFunc(imRange.getLowerBound(), imRange.getUpperBound(), 1000);
			EvenlyDiscretizedFunc distFunc = HistogramFunction.getEncompassingHistogram(distRange.getLowerBound(), distRange.getUpperBound(), 1d);
			
			List<List<PlotSpec>> perCombGMPlots = new ArrayList<>();
			for (int p=0; p<periods.length; p++)
				perCombGMPlots.add(new ArrayList<>());
			
			for (double dist : dists) {
				Location siteLoc = dist == 0d ? loc : LocationUtils.location(loc, 0d, dist);
				Site site = new Site(siteLoc);
				site.addParameterList(gmm.getSiteParams());
				gmm.setSite(site);
				
				WeightedList<DistanceCorrected> myPointSurfs = pointSurf.getForDistanceCorrection(siteLoc, distCorr, trt, mag);
				WeightedList<? extends RuptureSurface> myFiniteSurfs;
				if (forceHW == null) {
					myFiniteSurfs = finiteSurfs;
				} else {
					WeightedList<RuptureSurface> matchingFiniteSurfs = new WeightedList<>();
					for (WeightedValue<? extends RuptureSurface> surfVal : finiteSurfs) {
						boolean surfHW = surfVal.value.getDistanceX(siteLoc) >= 0d;
						if (surfHW == forceHW.booleanValue())
							matchingFiniteSurfs.add(surfVal.value, surfVal.weight);
					}
					if (matchingFiniteSurfs.isEmpty()) {
						System.out.println("ForceHW="+forceHW+" but no finite surfaces match at rEpi="+oDF.format(dist)+", skipping\n");
						continue;
					} else {
						matchingFiniteSurfs.normalize();
						System.out.println("ForceHW="+forceHW+", found "+matchingFiniteSurfs.size()
								+" matching finite surfaces ("+pDF.format((double)matchingFiniteSurfs.size()/(double)finiteSurfs.size())+")");
					}
					myFiniteSurfs = matchingFiniteSurfs;
					
					WeightedList<DistanceCorrected> matchingPointSurfs = new WeightedList<>();
					for (WeightedValue<DistanceCorrected> surfVal : myPointSurfs) {
						boolean surfHW = surfVal.value.getDistanceX(siteLoc) >= 0d;
						if (surfHW == forceHW.booleanValue())
							matchingPointSurfs.add(surfVal.value, surfVal.weight);
					}
					if (matchingPointSurfs.isEmpty()) {
						System.out.println("ForceHW="+forceHW+" but no point surfaces match at rEpi="+oDF.format(dist)+", skipping\n");
						continue;
					} else {
						matchingPointSurfs.normalize();
						System.out.println("ForceHW="+forceHW+", found "+matchingPointSurfs.size()
								+" matching point surfaces ("+pDF.format((double)matchingPointSurfs.size()/(double)myPointSurfs.size())+")");
					}
					myPointSurfs = matchingPointSurfs;
				}
				
				List<XY_DataSet> rJBFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> rJBChars = new ArrayList<>();
				List<XY_DataSet> rRupFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> rRupChars = new ArrayList<>();
				List<XY_DataSet> rRupVsJBFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> rRupVsJBChars = new ArrayList<>();
				
				List<List<XY_DataSet>> gmDistFuncs = new ArrayList<>();
				for (int p=0; p<periods.length; p++)
					gmDistFuncs.add(new ArrayList<>());
				List<PlotCurveCharacterstics> gmDistChars = new ArrayList<>();
				
				double minDist = Double.POSITIVE_INFINITY;
				double maxDist = 0d;
				double histMax = 0d;
				
				String hwAnnAdd = ", HW:";
				
				for (boolean finite : new boolean[] { true, false }) {
					WeightedList<? extends RuptureSurface> surfs = finite ? myFiniteSurfs : myPointSurfs;
					String type = finite ? "finite" : "point";

					System.out.println("Calculating " + type + " surface for rEpi=" + oDF.format(dist) + " km");
					
					EvenlyDiscretizedFunc rJBHist = new EvenlyDiscretizedFunc(distFunc.getMinX(), distFunc.size(), distFunc.getDelta());
					EvenlyDiscretizedFunc rRupHist = new EvenlyDiscretizedFunc(distFunc.getMinX(), distFunc.size(), distFunc.getDelta());
					DefaultXY_DataSet rRupVsJB_hw = new DefaultXY_DataSet();
					DefaultXY_DataSet rRupVsJB_fw = new DefaultXY_DataSet();
					
					MinMaxAveTracker rJBtrack = new MinMaxAveTracker();
					MinMaxAveTracker rRuptrack = new MinMaxAveTracker();
					double rJBAvg = 0d;
					double rRupAvg = 0d;
					double fractHW = 0d;
					List<Range> hwRanges = new ArrayList<>();
					double curFirstHWStrike = Double.NaN;
					double curLastHWStrike = Double.NaN;
					
					EvenlyDiscretizedFunc[] exceedFuncs = new EvenlyDiscretizedFunc[periods.length];
					for (int p=0; p<periods.length; p++)
						exceedFuncs[p] = new EvenlyDiscretizedFunc(logIMfunc.getMinX(), logIMfunc.size(), logIMfunc.getDelta());
					
					Preconditions.checkState(surfs.isNormalized());
					
					double closestTestJBDiff = Double.POSITIVE_INFINITY;
					double closestTestJB = Double.NaN;
					double closestRupToTestJB = Double.NaN;
					
					for (WeightedValue<? extends RuptureSurface> weightedSurf : surfs) {
						RuptureSurface surf = weightedSurf.value;
						double weight = weightedSurf.weight;
						EqkRupture rup = new EqkRupture(mag, mech.rake(), surf, null);
						gmm.setEqkRupture(rup);
						
						double rJB = surf.getDistanceJB(siteLoc);
						double rRup = surf.getDistanceRup(siteLoc);
						boolean hw = surf.getDistanceX(siteLoc) >= 0d;
						
						if (finite && hw && Double.isFinite(testJB)) {
							double diff = Math.abs(testJB - rJB);
							if (diff < closestTestJBDiff) {
								closestTestJB = rJB;
								closestTestJBDiff = diff;
								closestRupToTestJB = rRup;
							}
						}
						
						if (!finite && printPointRrups && hw)
							System.out.println("\tEstimated hanging-wall rRup="+(float)rRup+" for rJB="+(float)rJB);
						
						if (hw)
							rRupVsJB_hw.set(rJB, rRup);
						else
							rRupVsJB_fw.set(rJB, rRup);
						
						minDist = Math.min(minDist, Math.min(rJB, rRup));
						maxDist = Math.max(maxDist, Math.max(rJB, rRup));
						
						rJBHist.add(distFunc.getClosestXIndex(rJB), weight);
						rRupHist.add(distFunc.getClosestXIndex(rRup), weight);
						rJBtrack.addValue(rJB);
						rRuptrack.addValue(rRup);
						rJBAvg += rJB * weight;
						rRupAvg += rRup * weight;
						if (hw) {
							fractHW += weight;
							if (trackHWRanges && !(surf instanceof PointSurface)) {
								// finite
								double strike = surf.getUpperEdge().getAveStrike();
								if (Double.isNaN(curFirstHWStrike))
									curFirstHWStrike = strike;
								curLastHWStrike = strike;
							}
						} else {
							if (Double.isFinite(curFirstHWStrike))
								hwRanges.add(new Range(curFirstHWStrike, curLastHWStrike));
							curFirstHWStrike = Double.NaN;
							curLastHWStrike = Double.NaN;
						}

						for (int p=0; p<periods.length; p++) {
							double period = periods[p];
							if (period == 0d) {
								gmm.setIntensityMeasure(PGA_Param.NAME);
							} else {
								gmm.setIntensityMeasure(SA_Param.NAME);
								SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
							}
							
							gmm.getExceedProbabilities(logIMfunc);
							
							for (int i=0; i<logIMfunc.size(); i++)
								exceedFuncs[p].add(i, weight*logIMfunc.getY(i));
						}
					}
					System.out.println("\trJB stats:\t"+rJBtrack);
					System.out.println("\trJB weighted avg:\t"+(float)rJBAvg);
					System.out.println("\trRup stats:\t"+rRuptrack);
					System.out.println("\trRup weighted avg:\t"+(float)rRupAvg);
					System.out.println("\tFract HW:\t"+(float)fractHW);
					
					if (finite && Double.isFinite(testJB) && closestTestJBDiff < 1) {
						System.out.println("Found a close finite FW case to rJB="+(float)testJB);
						System.out.println("\tcloseJB="+(float)closestTestJB);
						System.out.println("\trRup="+(float)closestRupToTestJB);
						System.out.println("Doing pt src equiv");
						double ptRrup = RjbDistributionDistanceCorrection.getCorrDistRup(closestTestJB, dist, zTop,
								pointSurf.getAveRupBottomDepth(), dipRad, length, pointSurf.getAveHorizontalWidth(), false);
						System.out.println("\tPt rRup="+(float)ptRrup);
					}
					
					if (finite)
						hwAnnAdd += " F=";
					else
						hwAnnAdd += " Pt=";
					hwAnnAdd += pDF.format(fractHW);
					
					if (Double.isFinite(curFirstHWStrike))
						hwRanges.add(new Range(curFirstHWStrike, curLastHWStrike));
					if (!hwRanges.isEmpty()) {
						System.out.print("\tHW ranges:");
						for (Range hwRange : hwRanges)
							System.out.print("\t["+oDF.format(hwRange.getLowerBound())+", "+oDF.format(hwRange.getUpperBound())+"]");
						System.out.println();
					}
					
					DefaultXY_DataSet meanRjbLine = new DefaultXY_DataSet();
					meanRjbLine.set(rJBAvg, 0d);
					meanRjbLine.set(rJBAvg, 1d);
					
					DefaultXY_DataSet meanRrupLine = new DefaultXY_DataSet();
					meanRrupLine.set(rRupAvg, 0d);
					meanRrupLine.set(rRupAvg, 1d);
					
					Color histColor;
					Color lineColor;
					PlotSymbol scatterSymbol;
					Color scatterColorFW;
					Color scatterColorHW;
					PlotLineType gmLineType;
					String name;
					if (finite) {
						name = "Finite";
//						histColor = Colors.tab_lightorange;
						lineColor = Colors.tab_orange;
						histColor = new Color(lineColor.getRed(), lineColor.getGreen(), lineColor.getBlue(), 127);
						scatterColorHW = Colors.tab_orange;
						scatterColorFW = Colors.tab_lightorange;
						scatterSymbol = PlotSymbol.BOLD_CROSS;
						gmLineType = PlotLineType.SOLID;
					} else {
						name = "Point Source";
//						histColor = Colors.tab_lightblue;
						lineColor = Colors.tab_blue;
						histColor = new Color(lineColor.getRed(), lineColor.getGreen(), lineColor.getBlue(), 127);
						scatterColorHW = Colors.tab_blue;
						scatterColorFW = Colors.tab_lightblue;
						scatterSymbol = PlotSymbol.BOLD_X;
						gmLineType = PlotLineType.SHORT_DASHED;
					}
					rJBHist.setName(name+" Distribution");
					rRupHist.setName(name+" Distribution");
					meanRjbLine.setName(name+" Mean");
					meanRrupLine.setName(name+" Mean");
					rRupVsJB_fw.setName(name+" (FW)");
					rRupVsJB_hw.setName(name+" (HW)");
					if (!finite) {
						rRupVsJB_fw.setName(rRupVsJB_fw.size()+" "+rRupVsJB_fw.getName());
						rRupVsJB_hw.setName(rRupVsJB_hw.size()+" "+rRupVsJB_hw.getName());
					}
					exceedFuncs[0].setName(name);
					
					rJBFuncs.add(rJBHist);
					rJBChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, histColor));
					rJBFuncs.add(meanRjbLine);
					rJBChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, lineColor));
					
					rRupFuncs.add(rRupHist);
					rRupChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, histColor));
					rRupFuncs.add(meanRrupLine);
					rRupChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, lineColor));
					
					if (forceHW == null || !forceHW) {
						rRupVsJBFuncs.add(rRupVsJB_fw);
						rRupVsJBChars.add(new PlotCurveCharacterstics(scatterSymbol, 5f, scatterColorFW));
					}
					if (forceHW == null || forceHW) {
						rRupVsJBFuncs.add(rRupVsJB_hw);
						rRupVsJBChars.add(new PlotCurveCharacterstics(scatterSymbol, 5f, scatterColorHW));
					}
					
					for (int p=0; p<periods.length; p++)
						gmDistFuncs.get(p).add(exceedFuncs[p]);
					gmDistChars.add(new PlotCurveCharacterstics(gmLineType, 3f, lineColor));
					
					histMax = Math.max(histMax, Math.max(rJBHist.getMaxY(), rRupHist.getMaxY()));
				}
				System.out.println();
				
				Range plotDistRange = new Range(Math.max(0d, minDist-5), maxDist+5d);
				Range plotHistYRange = new Range(0d, Math.min(1d, histMax+0.1));
				Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
				
				PlotSpec rJBplot = new PlotSpec(rJBFuncs, rJBChars, null, "Distance (km)", "Fraction");
				rJBplot.setLegendVisible(true);
				XYTextAnnotation rJBAnn = new XYTextAnnotation("rJB"+hwAnnAdd,
						plotDistRange.getLowerBound()+0.9*plotDistRange.getLength(), 0.8*plotHistYRange.getUpperBound());
				rJBAnn.setFont(annFont);
				rJBAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
				rJBplot.addPlotAnnotation(rJBAnn);
				
				PlotSpec rRupPlot = new PlotSpec(rRupFuncs, rRupChars, null, "Distance (km)", "Fraction");
				rRupPlot.setLegendVisible(false);
				XYTextAnnotation rRupAnn = new XYTextAnnotation("rRup"+hwAnnAdd, rJBAnn.getX(), rJBAnn.getY());
				rRupAnn.setFont(annFont);
				rRupAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
				rRupPlot.addPlotAnnotation(rRupAnn);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(List.of(rJBplot, rRupPlot), false, false, List.of(plotDistRange), List.of(plotHistYRange, plotHistYRange));
				
				String distPrefix = oDF.format(dist)+"km";
				PlotUtils.writePlots(outputDir, distPrefix+"_distances", gp, 800, 800, true, false, false);
				
				DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
				oneToOne.set(plotDistRange.getLowerBound(), plotDistRange.getLowerBound());
				oneToOne.set(plotDistRange.getUpperBound(), plotDistRange.getUpperBound());
				
				rRupVsJBFuncs.add(0, oneToOne);
				rRupVsJBChars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
				
				PlotSpec rRupVsJBPlot = new PlotSpec(rRupVsJBFuncs, rRupVsJBChars, null, "Distance JB (km)", "Distance Rup (km)");
				rRupVsJBPlot.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
				
				gp.drawGraphPanel(rRupVsJBPlot, false, false, plotDistRange, plotDistRange);
				
				PlotUtils.writePlots(outputDir, distPrefix+"_distance_ratios", gp, 800, false, true, false, false);
				
				List<PlotSpec> gmDistPlots = new ArrayList<>();
				List<Range> gmDistYRanges = new ArrayList<>();
				for (int p=0; p<periods.length; p++) {
					PlotSpec plot = new PlotSpec(gmDistFuncs.get(p), gmDistChars, " ", "Ln(IML)", "P(> IML)");
					plot.setLegendVisible(p == 0);
					String text = periods[p] == 0d ? "PGA" : oDF.format(periods[p])+"s SA";
					XYTextAnnotation periodAnn = new XYTextAnnotation(text+hwAnnAdd, imRange.getLowerBound()+imRange.getLength()*0.9, 0.9);
					periodAnn.setFont(annFont);
					periodAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
					plot.addPlotAnnotation(periodAnn);
					gmDistPlots.add(plot);
					gmDistYRanges.add(new Range(0d, 1d));
					
					perCombGMPlots.get(p).add(plot);
				}
				
				gp.drawGraphPanel(gmDistPlots, false, false, List.of(imRange), gmDistYRanges);
				
				PlotUtils.writePlots(outputDir, distPrefix+"_gms", gp, 800, 800, true, false, false);
			}
			
			for (int p=0; p<periods.length; p++) {
				
			}
		}
	}

}
