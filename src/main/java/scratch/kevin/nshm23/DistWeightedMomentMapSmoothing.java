package scratch.kevin.nshm23;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class DistWeightedMomentMapSmoothing {
	
	public static GriddedGeoDataSet distKernelSmooth(GriddedGeoDataSet momentMap, ScalarIMR gmpe, double refMag, double period,
			boolean normalize) {
		return distKernelSmooth(new GriddedGeoDataSet[] {momentMap}, gmpe, refMag, period, normalize)[0];
	}
	
	public static GriddedGeoDataSet[] distKernelSmooth(GriddedGeoDataSet[] momentMaps, ScalarIMR gmpe, double refMag, double period,
			boolean normalize) {
		GriddedGeoDataSet[] ret = new GriddedGeoDataSet[momentMaps.length];
		GriddedRegion gridReg = momentMaps[0].getRegion();
		for (int i=0; i<ret.length; i++) {
			ret[i] = new GriddedGeoDataSet(gridReg);
			if (i > 0)
				Preconditions.checkState(momentMaps[i].getRegion().equalsRegion(gridReg), "Regions are inconsistent");
		}
		
		gmpe.setParamDefaults();
		if (period == 0d) {
			gmpe.setIntensityMeasure(PGA_Param.NAME);
		} else if (period > 0d) {
			gmpe.setIntensityMeasure(SA_Param.NAME);
			SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
		} else {
			throw new IllegalStateException("Period must be >= 0");
		}
		
		Location siteLoc = new Location(0d, 0d);
		Site site = new Site(siteLoc);
		site.addParameterList(gmpe.getSiteParams());
		
		gmpe.setSite(site);
		
		double maxDist = LocationUtils.horzDistanceFast(new Location(gridReg.getMinGridLat(), gridReg.getMinGridLon()),
				new Location(gridReg.getMaxGridLat(), gridReg.getMaxGridLon()));
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(0d, maxDist, 1000);
		
		for (int i=0; i<distFunc.size(); i++) {
			double dist = distFunc.getX(i);
			
			Location rupLoc = dist == 0d ? siteLoc : LocationUtils.location(siteLoc, 0d, dist);
			PointSurface rupSurf = new PointSurface(rupLoc);
			rupSurf.setDistCorrMagAndType(Double.NaN, PtSrcDistCorr.Type.NONE);
			rupSurf.setAveDip(90d);
			rupSurf.setAveStrike(0d);
			
			EqkRupture rup = new EqkRupture(refMag, 0d, rupSurf, rupLoc);
			
			gmpe.setEqkRupture(rup);
			
			double val = Math.exp(gmpe.getMean());
			distFunc.set(i, val);
		}
		
		distFunc.scale(1d/distFunc.getY(0));
		
		int nodes = gridReg.getNodeCount();
		
		List<Future<?>> futures = new ArrayList<>();
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		for (int i=0; i<nodes; i++) {
			boolean skip = true;
			for (GriddedGeoDataSet map : momentMaps) {
				if (map.get(i) > 0) {
					skip = false;
					break;
				}
			}
			if (skip)
				continue;
			int srcIndex = i;
			Location l1 = gridReg.getLocation(i);
			futures.add(exec.submit(new Runnable() {
				
				@Override
				public void run() {
					double[] scalars = new double[nodes];
					for (int j=0; j<nodes; j++) {
						Location l2 = gridReg.getLocation(j);
						double dist = LocationUtils.horzDistanceFast(l1, l2);
						if (dist < maxDist) {
							scalars[j] = distFunc.getInterpolatedY(dist);
						}
					}
					for (int m=0; m<ret.length; m++) {
						GriddedGeoDataSet src = momentMaps[m];
						GriddedGeoDataSet dest = ret[m];
						double moment = src.get(srcIndex);
						if (moment > 0d) {
							synchronized(dest) {
								for (int j=0; j<nodes; j++)
									dest.set(j, dest.get(j) + moment*scalars[j]);
							}
						}
					}
				}
			}));
		}
		
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		exec.shutdown();
		
		if (normalize)
			for (int m=0; m<ret.length; m++)
				ret[m].scale(momentMaps[m].getSumZ()/ret[m].getSumZ());
		return ret;
	}
	
	public static GriddedGeoDataSet hazMapSmooth(GriddedGeoDataSet momentMap, AttenRelRef gmpeRef, double refMag,
			double period, ReturnPeriods rp) {
		AbstractERF erf = new MomentScaledFixedMagERF(momentMap, refMag);
		
		return hazMapSmooth(momentMap, gmpeRef, period, rp, erf);
	}
	
	public static GriddedGeoDataSet hazMapSmooth(GriddedGeoDataSet momentMap, AttenRelRef gmpeRef,
			IncrementalMagFreqDist mfdShape, double period, ReturnPeriods rp) {
		AbstractERF erf = new MomentScaledMFDShapeERF(momentMap, mfdShape);
		
		return hazMapSmooth(momentMap, gmpeRef, period, rp, erf);
	}
	
	private static GriddedGeoDataSet hazMapSmooth(GriddedGeoDataSet momentMap, AttenRelRef gmpeRef, double period,
			ReturnPeriods rp, AbstractERF erf) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(momentMap.getRegion());
		
		GriddedRegion gridReg = momentMap.getRegion();
		double maxDist = Math.min(500d, LocationUtils.horzDistanceFast(new Location(gridReg.getMinGridLat(), gridReg.getMinGridLon()),
				new Location(gridReg.getMaxGridLat(), gridReg.getMaxGridLon())));
		
		ArrayDeque<ScalarIMR> gmpeDeque = new ArrayDeque<>();
		ArrayDeque<HazardCurveCalculator> curveCalcDqeue = new ArrayDeque<>();
		
		ScalarIMR refGMPE = checkOutGMPE(gmpeDeque, gmpeRef, period);
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(refGMPE.getIntensityMeasure());
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 0d);
		
		checkInGMPE(gmpeDeque, refGMPE);
		
		List<Future<Double>> futures = new ArrayList<>();
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		System.out.println("Calculating proxy hazard map for "+ret.size()+" sites");
		
		for (int i=0; i<ret.size(); i++) {
			Site site = new Site(ret.getLocation(i));
			futures.add(exec.submit(new Callable<Double>() {

				@Override
				public Double call() throws Exception {
					ScalarIMR gmpe = checkOutGMPE(gmpeDeque, gmpeRef, period);
					site.addParameterList(gmpe.getSiteParams());
					HazardCurveCalculator calc = checkOutCurveCalc(curveCalcDqeue, maxDist);
					
					DiscretizedFunc logCurve = logXVals.deepClone();
					calc.getHazardCurve(logCurve, site, gmpe, erf);
					DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
					for (int i=0; i<logCurve.size(); i++)
						curve.set(xVals.getX(i), logCurve.getY(i));
					
					double curveLevel = rp.oneYearProb;
					double curveVal;
					if (curveLevel > curve.getMaxY())
						curveVal = 0d;
					else if (curveLevel < curve.getMinY())
						// saturated
						curveVal = curve.getMaxX();
					else
						curveVal = curve.getFirstInterpolatedX_inLogXLogYDomain(curveLevel);
//					if (ret.indexOf(new Location(32.2, -115)) == nodeIndex) {
//						System.out.println("Debug for "+site.getLocation()+", val="+curveVal);
//						System.out.println("Curve:\n"+curve);
//					}
					
					checkInGMPE(gmpeDeque, gmpe);
					checkInCurveCalc(curveCalcDqeue, calc);
					
					return curveVal;
				}
			}));
		}
		
		MinMaxAveTracker imlTrack = new MinMaxAveTracker();
		
		for (int i=0; i<ret.size(); i++) {
			try {
				double val = futures.get(i).get();
				ret.set(i, val);
				imlTrack.addValue(val);
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		System.out.println("Done. Proxy Map z-range: "+imlTrack);
		
		exec.shutdown();
		
		return ret;
	}
	
	private static class MomentScaledFixedMagERF extends AbstractERF {
		
		private List<ProbEqkSource> sources;
		
		public MomentScaledFixedMagERF(GriddedGeoDataSet momentMap, double refMag) {
			sources = new ArrayList<>();
			
			double magMoment = MagUtils.magToMoment(refMag);
			
			System.out.println("Building moment-based ERF for M"+(float)refMag+" which equates to "+(float)magMoment);
			
			MinMaxAveTracker rateTrack = new MinMaxAveTracker();
			MinMaxAveTracker probTrack = new MinMaxAveTracker();
			
			for (int i=0; i<momentMap.size(); i++) {
				double momentRate = momentMap.get(i);
				if (momentRate == 0d || !Double.isFinite(momentRate))
					continue;
				
				double rate = momentRate/magMoment;
				double prob = 1d - Math.exp(-rate);
				
				rateTrack.addValue(rate);
				probTrack.addValue(prob);
				
				Location rupLoc = momentMap.getLocation(i);
				PointSurface rupSurf = new PointSurface(rupLoc);
				rupSurf.setDistCorrMagAndType(Double.NaN, PtSrcDistCorr.Type.NONE);
				rupSurf.setAveDip(90d);
				rupSurf.setAveStrike(0d);
				
				ProbEqkRupture rup = new ProbEqkRupture(refMag, 0d, prob, rupSurf, rupLoc);
				sources.add(new SingleLocSource(rup, rupLoc));
			}
			
			System.out.println("\tRate range: "+rateTrack);
			System.out.println("\tProb range: "+probTrack);
		}

		@Override
		public int getNumSources() {
			return sources.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return sources.get(idx);
		}

		@Override
		public void updateForecast() {}

		@Override
		public String getName() {
			return null;
		}
		
	}
	
	private static class MomentScaledMFDShapeERF extends AbstractERF {
		
		private List<ProbEqkSource> sources;
		
		public MomentScaledMFDShapeERF(GriddedGeoDataSet momentMap, IncrementalMagFreqDist mfd) {
			sources = new ArrayList<>();
			
			System.out.println("Building MFD-shaped moment-based ERF");
			
			MinMaxAveTracker rateTrack = new MinMaxAveTracker();
			MinMaxAveTracker probTrack = new MinMaxAveTracker();
			
			for (int i=0; i<momentMap.size(); i++) {
				double momentRate = momentMap.get(i);
				if (momentRate == 0d || !Double.isFinite(momentRate))
					continue;
				
				Location rupLoc = momentMap.getLocation(i);
				PointSurface rupSurf = new PointSurface(rupLoc);
				rupSurf.setDistCorrMagAndType(Double.NaN, PtSrcDistCorr.Type.NONE);
				rupSurf.setAveDip(90d);
				rupSurf.setAveStrike(0d);
				
				IncrementalMagFreqDist myMFD = mfd.deepClone();
				myMFD.scaleToTotalMomentRate(momentRate);
				
				List<ProbEqkRupture> rups = new ArrayList<>();
				for (Point2D pt : myMFD) {
					double rate = pt.getY();
					double prob = 1d - Math.exp(-rate);
					
					rateTrack.addValue(rate);
					probTrack.addValue(prob);
					
					ProbEqkRupture rup = new ProbEqkRupture(pt.getX(), 0d, prob, rupSurf, rupLoc);
					rups.add(rup);
				}
				
				
				sources.add(new SingleLocSource(rups, rupLoc));
			}
			
			System.out.println("\tRate range: "+rateTrack);
			System.out.println("\tProb range: "+probTrack);
		}

		@Override
		public int getNumSources() {
			return sources.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return sources.get(idx);
		}

		@Override
		public void updateForecast() {}

		@Override
		public String getName() {
			return null;
		}
		
	}
	
	private static class SingleLocSource extends ProbEqkSource {
		
		private List<ProbEqkRupture> rups;
		private Location rupLoc;

		public SingleLocSource(ProbEqkRupture rup, Location rupLoc) {
			this(List.of(rup), rupLoc);
		}

		public SingleLocSource(List<ProbEqkRupture> rups, Location rupLoc) {
			this.rups = rups;
			this.rupLoc = rupLoc;
			isPoissonian = true;
		}

		@Override
		public LocationList getAllSourceLocs() {
			LocationList ret = new LocationList();
			ret.add(rupLoc);
			return ret;
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return rups.get(0).getRuptureSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return LocationUtils.horzDistanceFast(site.getLocation(), rupLoc);
		}

		@Override
		public int getNumRuptures() {
			return rups.size();
		}

		@Override
		public ProbEqkRupture getRupture(int nRupture) {
			return rups.get(nRupture);
		}
		
	}
	
	private static synchronized ScalarIMR checkOutGMPE(ArrayDeque<ScalarIMR> deque, AttenRelRef gmpeRef, double period) {
		ScalarIMR ret;
		if (deque.isEmpty()) {
			ret = gmpeRef.get();
			
			ret.setParamDefaults();
			if (period == 0d) {
				ret.setIntensityMeasure(PGA_Param.NAME);
			} else if (period > 0d) {
				ret.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(ret.getIntensityMeasure(), period);
			} else {
				throw new IllegalStateException("Period must be >= 0");
			}
		} else {
			ret = deque.pop();
		}
		return ret;
	}
	
	private static synchronized void checkInGMPE(ArrayDeque<ScalarIMR> deque, ScalarIMR gmpe) {
		deque.push(gmpe);
	}
	
	private static synchronized HazardCurveCalculator checkOutCurveCalc(ArrayDeque<HazardCurveCalculator> deque, double maxDist) {
		HazardCurveCalculator ret;
		if (deque.isEmpty()) {
			ret = new HazardCurveCalculator();
			ret.setMaxSourceDistance(maxDist);
		} else {
			ret = deque.pop();
		}
		return ret;
	}
	
	private static synchronized void checkInCurveCalc(ArrayDeque<HazardCurveCalculator> deque, HazardCurveCalculator curveCalc) {
		deque.push(curveCalc);
	}
	
	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		GriddedRegion gridReg = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.1d, GriddedRegion.ANCHOR_0_0);
		GriddedGeoDataSet moRates = SingleSiteHazardAndDataComparisonPageGen.calcFSSPartic(sol, gridReg, new double[0], null).momentRateMap;
		
		hazMapSmooth(moRates, AttenRelRef.ASK_2014, 7d, 0d, ReturnPeriods.TWO_IN_50);
	}

}
