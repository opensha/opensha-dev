package scratch.kevin.griddedInterpGMPE;

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

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
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

import scratch.aftershockStatisticsETAS.ETAS_ShakingForecastCalc;

public class GriddedInterpGMPE_Calc {
	
	private ScalarIMR gmpe;
	
	private DistanceInterpolator distInterp;
	private MagnitudeInterpolator magInterp;
	private FocalMechInterpolator mechInterp;
	private AbstractGMPEInterpolation<?>[] otherInterps;
	
	private List<AbstractGMPEInterpolation<?>> allInterps;
	
	private IncrementalMagFreqDist inputMFD;
	private Map<FocalMech, Double> mechWtMap;
	
	private NDimArrayCalc arrayCalc;
	private double[] means;
	private double[] stdDevs;
	private NDimensionalLinearInterpolation interpolator;
	
	public GriddedInterpGMPE_Calc(ScalarIMR gmpe, DistanceInterpolator distInterp,
			MagnitudeInterpolator magInterp, FocalMechInterpolator mechInterp,
			AbstractGMPEInterpolation<?>... otherInterps) {
		this.gmpe = gmpe;
		
		this.distInterp = distInterp;
		this.magInterp = magInterp;
		this.mechInterp = mechInterp;
		this.otherInterps = otherInterps;
		
		allInterps = new ArrayList<>();
		
		allInterps.add(magInterp);
		allInterps.add(distInterp);
		allInterps.add(mechInterp);
		for (AbstractGMPEInterpolation<?> interp : otherInterps)
			allInterps.add(interp);
		
		// build input MFD with mag range
		inputMFD = new IncrementalMagFreqDist(magInterp.getValue(0),
				magInterp.getValue(magInterp.getNumBins()-1),magInterp.getNumBins());
		for (int i=0; i<inputMFD.size(); i++)
			inputMFD.set(i, 1d);
		
		mechWtMap = new HashMap<>();
		double wtEach = 1d/(double)mechInterp.getNumBins();
		for (FocalMech mech : mechInterp)
			mechWtMap.put(mech, wtEach);
		
		int[] dimensions = new int[allInterps.size()];
		for (int i=0; i<dimensions.length; i++)
			dimensions[i] = allInterps.get(i).getNumBins();
		arrayCalc = new NDimArrayCalc(dimensions);
		interpolator = new NDimensionalLinearInterpolation(dimensions.length);
	}
	
	public void precalc(double[] depths) {
		Location loc = new Location(0d, 0d);
		
		PointSource13b source = new PointSource13b(loc, inputMFD, 1d, depths, mechWtMap);
		
		System.out.println("Precalculating GMPE for "+allInterps.size()+" dimensions, "+arrayCalc.rawArraySize()+" values");
		
		means = new double[arrayCalc.rawArraySize()];
		stdDevs = new double[arrayCalc.rawArraySize()];
		
		Set<EqkRupture> ruptures = new HashSet<>();
		for (ProbEqkRupture rup : source)
			ruptures.add((EqkRupture)rup.clone());
		
		Site site = new Site(loc);
		site.addParameterList(gmpe.getSiteParams());
		gmpe.setSite(site);
		gmpe.setEqkRupture(ruptures.iterator().next());
		
		numRupsTrack = new MinMaxAveTracker();
		calcRecursive(new int[0], source, ruptures);
		System.out.println("Done precalculating!");
		System.out.println("Num Rups tracker: "+numRupsTrack);
	}
	
	private MinMaxAveTracker numRupsTrack;
	
	private void calcRecursive(int[] upstreamIndexes, PointSource13b source, Set<EqkRupture> curRuptures) {
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
				double wtEach = 1d/rups.size();
				for (EqkRupture rup : rups) {
					gmpe.setEqkRupture(rup);
					means[arrayIndex] = wtEach*gmpe.getMean();
					stdDevs[arrayIndex] = wtEach*gmpe.getStdDev();
				}
			} else {
				// pass to the next level
				calcRecursive(indexes, source, rups);
			}
		}
	}
	
	public DiscretizedFunc[] calc(List<ProbEqkSource> sources, List<Site> sites, DiscretizedFunc xVals) {
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1);
		
		DiscretizedFunc[] curves = new DiscretizedFunc[sites.size()];
		
		DiscretizedFunc condProbFunc = new LightFixedXFunc(logXVals);
		
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
				indexes[2] = mechInterp.getInterpolatedBinIndex(FocalMechInterpolator.forRake(rup.getAveRake()));
				for (int s=0; s<sites.size(); s++) {
					Site site = sites.get(s);
					if (distances[s] > distInterp.getMax())
						continue;
//					System.out.println("dist: "+distances[s]);
					// set dist and any other params
					indexes[1] = distInterp.getInterpolatedBinIndex(distances[s]);
					for (int j=3; j<allInterps.size(); j++)
						indexes[j] = allInterps.get(j).detectInterpolatedBinIndex(gmpe, site);
					
					// multi-dimensional interpolation
					double mean = interpolator.interpolate(means, arrayCalc, indexes);
					double stdDev = interpolator.interpolate(stdDevs, arrayCalc, indexes);
					
					for (int i=0; i<curves[s].size(); i++) {
						double x = condProbFunc.getX(i);
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
		
		DistanceInterpolator distInterp = new DistanceInterpolator(0d, 200d, 20);
		
		double refMag = 5d;
		double maxMag = 8.5d;
		double magDelta = 0.1;
		int numMag = (int)((maxMag - refMag)/magDelta + 0.5) + 1;
		System.out.println("NumMag: "+numMag);
		MagnitudeInterpolator magInterp = new MagnitudeInterpolator(refMag, maxMag, numMag);
		
		FocalMechInterpolator mechInterp = new FocalMechInterpolator();
		
		DoubleParameterInterpolator vs30Interp = new DoubleParameterInterpolator(
				Vs30_Param.NAME, 150, 760, 20, false, false);
		
		GriddedInterpGMPE_Calc calc = new GriddedInterpGMPE_Calc(gmpe, distInterp, magInterp, mechInterp, vs30Interp);
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
		List<Map<FocalMech, Double>> mechWtsList = new ArrayList<>();
		mechWtsList.add(mechWts);
		
		ERF erf = new ETAS_ShakingForecastCalc.GriddedForecast(rateModel, refMag, maxMag, b, mechWtsList, depths, durationYears);
		erf.updateForecast();
		
		double calcSpacing = 0.5;
		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(rateModel.getMaxLat(), rateModel.getMaxLon()),
						new Location(rateModel.getMinLat(), rateModel.getMinLon())), calcSpacing, null);
		List<Site> sites = new ArrayList<>();
		List<Double> vs30s = vs30Provider.getValues(calcRegion.getNodeList());
		for (int i=0; i<calcRegion.getNodeCount(); i++) {
			Site site = new Site(calcRegion.locationForIndex(i));
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>) param.clone());
			site.getParameter(Double.class, Vs30_Param.NAME).setValue(vs30s.get(i));
			sites.add(site);
		}
		
		System.out.println(sites.size()+" sites");
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		
		System.out.println("Calculating interpolated");
		Stopwatch watch = Stopwatch.createStarted();
		calc.calc(erf.getSourceList(), sites, xVals);
		watch.stop();
		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" seconds");
		
		System.out.println("Calculating traditional");
		watch.reset();
		watch.start();
		calc.calcTraditional(erf, sites, xVals);
		watch.stop();
		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" seconds");
	}

}
