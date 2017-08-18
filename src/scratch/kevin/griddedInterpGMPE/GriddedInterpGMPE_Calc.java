package scratch.kevin.griddedInterpGMPE;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import javax.swing.JFrame;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Maps;

import scratch.aftershockStatisticsETAS.ETAS_ShakingForecastCalc;

public class GriddedInterpGMPE_Calc {
	
	private ScalarIMR gmpe;
	
	private DistanceInterpolator distInterp;
	private MagnitudeInterpolator magInterp;
	private AbstractGMPEInterpolation<?>[] otherInterps;
	
	private List<AbstractGMPEInterpolation<?>> allInterps;
	
	private IncrementalMagFreqDist inputMFD;
	
	private NDimArrayCalc arrayCalc;
	private Map<FocalMech, double[]> means;
	private Map<FocalMech, double[]> stdDevs;
	private NDimensionalLinearInterpolation interpolator;
	
	public GriddedInterpGMPE_Calc(ScalarIMR gmpe, DistanceInterpolator distInterp,
			MagnitudeInterpolator magInterp, AbstractGMPEInterpolation<?>... otherInterps) {
		this.gmpe = gmpe;
		
		this.distInterp = distInterp;
		this.magInterp = magInterp;
		this.otherInterps = otherInterps;
		
		allInterps = new ArrayList<>();
		
		allInterps.add(magInterp);
		allInterps.add(distInterp);
		for (AbstractGMPEInterpolation<?> interp : otherInterps)
			allInterps.add(interp);
		
		// build input MFD with mag range
		inputMFD = new IncrementalMagFreqDist(magInterp.getValue(0),
				magInterp.getValue(magInterp.getNumBins()-1),magInterp.getNumBins());
		for (int i=0; i<inputMFD.size(); i++)
			inputMFD.set(i, 1d);
		
		int[] dimensions = new int[allInterps.size()];
		for (int i=0; i<dimensions.length; i++)
			dimensions[i] = allInterps.get(i).getNumBins();
		arrayCalc = new NDimArrayCalc(dimensions);
		interpolator = new NDimensionalLinearInterpolation(dimensions.length);
		
		means = Maps.newHashMap();
		stdDevs = Maps.newHashMap();
	}
	
	public void precalc(double[] depths) {
		Location loc = new Location(0d, 0d);
		
		Map<FocalMech, Double> mechWtMap = new HashMap<>();
		double wtEach = 1d/FocalMech.values().length;
		for (FocalMech mech : FocalMech.values())
			mechWtMap.put(mech, wtEach);
		PointSource13b source = new PointSource13b(loc, inputMFD, 1d, depths, mechWtMap);
		
		System.out.println("Precalculating GMPE for "+allInterps.size()+" dimensions, "+arrayCalc.rawArraySize()+" values");
		
		Set<EqkRupture> ruptures = new HashSet<>();
		for (ProbEqkRupture rup : source)
			ruptures.add((EqkRupture)rup.clone());
		
		Site site = new Site(loc);
		site.addParameterList(gmpe.getSiteParams());
		gmpe.setSite(site);
		gmpe.setEqkRupture(ruptures.iterator().next());
		
		numRupsTrack = new MinMaxAveTracker();
		
		for (FocalMech mech : FocalMech.values()) {
			double[] means = new double[arrayCalc.rawArraySize()];
			double[] stdDevs = new double[arrayCalc.rawArraySize()];
			
			HashSet<EqkRupture> mechRups = new HashSet<>();
			for (EqkRupture rup : ruptures) {
				if (FocalMechInterpolator.forRake(rup.getAveRake()) == mech)
					mechRups.add(rup);
			}
			
			precalcRecursive(new int[0], source, mechRups, means, stdDevs);
			this.means.put(mech, means);
			this.stdDevs.put(mech, stdDevs);
		}
		System.out.println("Done precalculating!");
		System.out.println("Num Rups tracker: "+numRupsTrack);
	}
	
	private MinMaxAveTracker numRupsTrack;
	
	private void precalcRecursive(int[] upstreamIndexes, PointSource13b source, Set<EqkRupture> curRuptures, double[] means, double[] stdDevs) {
		int curIndex = upstreamIndexes.length;
		AbstractGMPEInterpolation<?> interp = allInterps.get(curIndex);
		int[] indexes = Arrays.copyOf(upstreamIndexes, curIndex+1);
		for (int i=0; i<interp.getNumBins(); i++) {
			indexes[indexes.length-1] = i;
			interp.setGMPE_Params(gmpe, source, i);
			Set<EqkRupture> rups = interp.getViableRuptures(curRuptures, i);
			Preconditions.checkState(!rups.isEmpty(), "No viable ruptures for %s=%s", interp.getName(), interp.getValue(i));
			
			if (indexes.length == allInterps.size()) {
				// time to actually calculate
				numRupsTrack.addValue(rups.size());
				int arrayIndex = arrayCalc.getIndex(indexes);
				Preconditions.checkState(means[arrayIndex] == 0, "Duplicate calc?");
				Preconditions.checkState(stdDevs[arrayIndex] == 0, "Duplicate calc?");
				double wtEach = 1d/rups.size();
				for (EqkRupture rup : rups) {
					gmpe.setEqkRupture(rup);
					// TODO remove hardcoded
					double current = (Double)gmpe.getParameter(Vs30_Param.NAME).getValue();
					double expected = (Double)allInterps.get(allInterps.size()-1).getValue(indexes[indexes.length-1]);
					Preconditions.checkState(current == expected,
							"Vs30 issue! Is actually %s, expected %s", current, expected);
					means[arrayIndex] += wtEach*gmpe.getMean();
					stdDevs[arrayIndex] += wtEach*gmpe.getStdDev();
				}
			} else {
				// pass to the next level
				precalcRecursive(indexes, source, rups, means, stdDevs);
			}
		}
	}
	
	public DiscretizedFunc[] calc(List<ProbEqkSource> sources, List<Site> sites, DiscretizedFunc xVals) {
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1);
		
		DiscretizedFunc[] curves = new DiscretizedFunc[sites.size()];
		
		for (int i=0; i<sites.size(); i++) {
			curves[i] = new LightFixedXFunc(xVals);
			// initialize the hazard function to 1.0
			for (int j=0; j<curves[i].size(); j++)
				curves[i].set(j, 1d);
		}
		
		double[] indexes = new double[allInterps.size()];
		
		for (ProbEqkSource source : sources) {
			// calculate distances
			Preconditions.checkState(source.getSourceSurface() instanceof PointSurface, "Only point sources supported");
			Location sourceLoc = ((PointSurface)source.getSourceSurface()).getLocation();
			double[] distances = new double[sites.size()];
			for (int s=0; s<sites.size(); s++)
				distances[s] = LocationUtils.horzDistanceFast(sourceLoc, sites.get(s).getLocation());
			for (ProbEqkRupture rup : source) {
				// these 2 are common to all
				Preconditions.checkState(rup.getMag() >= magInterp.getMin() && rup.getMag() <= magInterp.getMax(),
						"Rup mag=%s outisde of range [%s %s]", rup.getMag(), magInterp.getMin(), magInterp.getMax());
				indexes[0] = magInterp.getInterpolatedBinIndex(rup.getMag());
				FocalMech mech = FocalMechInterpolator.forRake(rup.getAveRake());
				double[] means = this.means.get(mech);
				double[] stdDevs = this.stdDevs.get(mech);
				for (int s=0; s<sites.size(); s++) {
					Site site = sites.get(s);
					if (distances[s] > distInterp.getMax())
						continue;
//					System.out.println("dist: "+distances[s]);
					// set dist and any other params
					indexes[1] = distInterp.getInterpolatedBinIndex(distances[s]);
					for (int j=2; j<allInterps.size(); j++)
						indexes[j] = allInterps.get(j).detectInterpolatedBinIndex(gmpe, site);
					
					if (Math.random() < 0.01)
						System.out.println("Indexes: "+(float)indexes[0]+" "+(float)indexes[1]+" "+(float)indexes[2]);
					
					// multi-dimensional interpolation
					double mean = interpolator.interpolate(means, arrayCalc, indexes);
					double stdDev = interpolator.interpolate(stdDevs, arrayCalc, indexes);
//					gmpe.setSite(site);
//					gmpe.setEqkRupture(rup);
//					double mean = gmpe.getMean();
//					double stdDev = gmpe.getStdDev();
					
					if (debugMeanScatter != null) {
						gmpe.setSite(site);
						gmpe.setEqkRupture(rup);
						double mean2 = gmpe.getMean();
						double stdDev2 = gmpe.getStdDev();
						debugMeanScatter.set(mean2, mean);
						debugStdDevScatter.set(stdDev2, stdDev);
					}
					
					for (int i=0; i<curves[s].size(); i++) {
						double x = logXVals.getX(i);
						double y = AttenuationRelationship.getExceedProbability(mean, stdDev, x, null, null);
						curves[s].set(i,curves[s].getY(i)*Math.pow(1-rup.getProbability(), y));
					}
				}
			}
		}
		
		// convert to exceedance probabilities (currently non-exceedance)
		for (DiscretizedFunc curve : curves)
			for (int i=0; i<curve.size(); i++)
				curve.set(i, 1d - curve.getY(i));
		
		return curves;
	}
	
	private DefaultXY_DataSet debugMeanScatter = new DefaultXY_DataSet();
	private DefaultXY_DataSet debugStdDevScatter = new DefaultXY_DataSet();
	
	private DiscretizedFunc[] calcTraditional(ERF erf, List<Site> sites, DiscretizedFunc xVals) {
		HazardCurveCalculator calc = new HazardCurveCalculator();
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1);
		
		DiscretizedFunc[] curves = new DiscretizedFunc[sites.size()];
		
		for (int i=0; i<sites.size(); i++) {
			DiscretizedFunc logCurve = new LightFixedXFunc(logXVals);
			calc.getHazardCurve(logCurve, sites.get(i), gmpe, erf);
			curves[i] = new LightFixedXFunc(xVals);
			// initialize the hazard function to 1.0
			for (int j=0; j<curves[i].size(); j++)
				curves[i].set(j, logCurve.getY(j));
		}
		
		return curves;
	}
	
	public static void main(String[] args) throws IOException {
		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(PGA_Param.NAME);
		
		DistanceInterpolator distInterp = new DistanceInterpolator(0d, 200d, 200);
		
		double refMag = 5d;
		double maxMag = 8.5d;
		double magDelta = 0.1;
		int numMag = (int)((maxMag - refMag)/magDelta + 0.5) + 1;
		System.out.println("NumMag: "+numMag);
		MagnitudeInterpolator magInterp = new MagnitudeInterpolator(refMag, maxMag, numMag);
		
//		FocalMechInterpolator mechInterp = new FocalMechInterpolator();
		
		DoubleParameterInterpolator vs30Interp = new DoubleParameterInterpolator(
				Vs30_Param.NAME, 180, 760, 200, false, false); // matches Wald Allen range
		
		GriddedInterpGMPE_Calc calc = new GriddedInterpGMPE_Calc(gmpe, distInterp, magInterp, vs30Interp);
		double[] depths = { 7, 2 };
		calc.precalc(depths);
		
		GeoDataSet rateModel = ArbDiscrGeoDataSet.loadXYZFile("/tmp/rateMap.txt", true);
		double b = 1;
		
		System.out.println(gmpe.getIntensityMeasure().getName());
		
		WaldAllenGlobalVs30 vs30Provider = new WaldAllenGlobalVs30();
		vs30Provider.setActiveCoefficients();
		
		double durationYears = 30d/365d;
		
		Map<FocalMech, Double> mechWts = new HashMap<>();
		mechWts.put(FocalMech.STRIKE_SLIP, 0.5);
		mechWts.put(FocalMech.NORMAL, 0.25);
		mechWts.put(FocalMech.REVERSE, 0.25);
//		mechWts.put(FocalMech.STRIKE_SLIP, 0d);
//		mechWts.put(FocalMech.NORMAL, 0.5);
//		mechWts.put(FocalMech.REVERSE, 0.5);
		List<Map<FocalMech, Double>> mechWtsList = new ArrayList<>();
		mechWtsList.add(mechWts);
		
		ERF erf = new ETAS_ShakingForecastCalc.GriddedForecast(rateModel, refMag, maxMag, b, mechWtsList, depths, durationYears);
		erf.updateForecast();
		
		double calcSpacing = 1.0;
		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(rateModel.getMaxLat(), rateModel.getMaxLon()),
						new Location(rateModel.getMinLat(), rateModel.getMinLon())), calcSpacing, null);
		List<Site> sites = new ArrayList<>();
		List<Double> vs30s = vs30Provider.getValues(calcRegion.getNodeList());
		for (int i=0; i<calcRegion.getNodeCount(); i++) {
			Site site = new Site(calcRegion.locationForIndex(i));
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>) param.clone());
//			site.getParameter(Double.class, Vs30_Param.NAME).setValue(vs30s.get(i));
			site.getParameter(Double.class, Vs30_Param.NAME).setValue(325d);
			sites.add(site);
		}
		
		System.out.println(sites.size()+" sites");
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		
		System.out.println("Calculating interpolated");
		Stopwatch watch = Stopwatch.createStarted();
		DiscretizedFunc[] interpCurves = calc.calc(erf.getSourceList(), sites, xVals);
		watch.stop();
		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" seconds");
		
		System.out.println("Calculating traditional");
		watch.reset();
		watch.start();
		DiscretizedFunc[] traditionalCurves = calc.calcTraditional(erf, sites, xVals);
		watch.stop();
		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" seconds");
		
		GriddedGeoDataSet interpMap = ETAS_ShakingForecastCalc.extractMap(calcRegion, interpCurves, false, 0.5);
		GriddedGeoDataSet traditionalMap = ETAS_ShakingForecastCalc.extractMap(calcRegion, traditionalCurves, false, 0.5);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, Math.max(interpMap.getMaxZ(), traditionalMap.getMaxZ()));
		
		XYZGraphPanel xyzGP = new XYZGraphPanel();
		
		XYZPlotSpec interpSpec = new XYZPlotSpec(interpMap, cpt, "Spatial Forecast", "Longitude", "Latitude",
				gmpe.getIntensityMeasure().getName());
		XYZPlotSpec traditionalSpec = new XYZPlotSpec(traditionalMap, cpt, "Spatial Forecast", "Longitude", "Latitude",
				gmpe.getIntensityMeasure().getName());
		
		List<XYZPlotSpec> specs = new ArrayList<>();
		specs.add(interpSpec);
		specs.add(traditionalSpec);
		
		Range xRange = new Range(calcRegion.getMinGridLon()-0.5*calcRegion.getLonSpacing(),
				calcRegion.getMaxGridLon()+0.5*calcRegion.getLonSpacing());
		Range yRange = new Range(calcRegion.getMinGridLat()-0.5*calcRegion.getLatSpacing(),
				calcRegion.getMaxGridLat()+0.5*calcRegion.getLatSpacing());
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(xRange);
		List<Range> yRanges = new ArrayList<>();
		yRanges.add(yRange);
		yRanges.add(yRange);
		
		xyzGP.drawPlot(specs, false, false, xRanges, yRanges, null);
		
		JFrame frame = new JFrame("");
		frame.setContentPane(xyzGP);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	
		GraphWindow gw = null;
		
		int curveMod = interpCurves.length/20;
		if (curveMod < 1)
			curveMod = 1;
		
		for (int i=0; i<interpCurves.length; i++) {
			if (i % curveMod > 0)
				continue;
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			interpCurves[i].setName("Interpolated");
			funcs.add(interpCurves[i]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
			
			traditionalCurves[i].setName("Traditional");
			funcs.add(traditionalCurves[i]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Site "+i, gmpe.getIntensityMeasure().getName(), "Probability");
			spec.setLegendVisible(true);
			
			if (gw == null) {
				gw = new GraphWindow(spec);
				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
				gw.setVisible(true);
			} else {
				gw.addTab(spec);
			}
			gw.setXLog(true);
			gw.setYLog(true);
		}
		
		if (calc.debugMeanScatter != null) {
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(calc.debugMeanScatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Mean Scatter", "GMPE", "Interpolated");
			gw = new GraphWindow(spec);
			Range range = new Range(Math.min(calc.debugMeanScatter.getMinX(), calc.debugMeanScatter.getMinY()),
					Math.max(calc.debugMeanScatter.getMaxX(), calc.debugMeanScatter.getMaxY()));
			gw.setX_AxisRange(range);
			gw.setY_AxisRange(range);
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			funcs.add(calc.debugStdDevScatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
			
			spec = new PlotSpec(funcs, chars, "Std Dev Scatter", "GMPE", "Interpolated");
			gw = new GraphWindow(spec);
			range = new Range(Math.min(calc.debugStdDevScatter.getMinX(), calc.debugStdDevScatter.getMinY()),
					Math.max(calc.debugStdDevScatter.getMaxX(), calc.debugStdDevScatter.getMaxY()));
			gw.setX_AxisRange(range);
			gw.setY_AxisRange(range);
		}
	}

}
