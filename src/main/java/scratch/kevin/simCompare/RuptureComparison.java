package scratch.kevin.simCompare;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.opensha.commons.data.Site;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

public abstract class RuptureComparison<E> {
	
	private E rupture;
	
	public RuptureComparison(E rupture) {
		this.rupture = rupture;
	}
	
	public E getRupture() {
		return rupture;
	}
	
	public abstract Collection<Site> getApplicableSites();
	
	public abstract double getLogMean(Site site, IMT imt);
	
	public abstract double getStdDev(Site site, IMT imt);
	
	public abstract double getDistanceRup(Site site);
	
	public abstract double getDistanceJB(Site site);
	
	public abstract boolean hasSite(Site site);
	
	public abstract boolean isComputed(Site site, IMT imt);
	
	public abstract EqkRupture getGMPERupture();
	
	public abstract double getMagnitude();
	
	public abstract double getAnnualRate();
	
	/**
	 * Occurrence time in years for this rupture, or NaN if not time based
	 * @return
	 */
	public abstract double getRuptureTimeYears();
	
	/**
	 * @param comps
	 * @return array containing minTime, maxTime, or null if times are NaN
	 */
	public static double[] getRuptureTimeRange(Iterable<? extends RuptureComparison<?>> comps) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (RuptureComparison<?> comp : comps) {
			double time = comp.getRuptureTimeYears();
			if (!Double.isFinite(time))
				return null;
			min = Math.min(time, min);
			max = Math.max(time, max);
		}
		return new double[] { min, max };
	}
	
	public static abstract class Cached<Y> extends RuptureComparison<Y> {
		
		private Map<Site, Map<IMT, double[]>> gms;
//		private Table<Site, IMT, Double> logMeans;
//		private Table<Site, IMT, Double> stdDevs;
		private Map<Site, Double> distanceRups;
		private Map<Site, Double> distanceJBs;

		public Cached(Y rupture) {
			super(rupture);
			
			gms = new ConcurrentHashMap<>();
//			logMeans = HashBasedTable.create();
//			stdDevs = HashBasedTable.create();
			distanceRups = new HashMap<>();
			distanceJBs = new HashMap<>();
		}
		
		public void addResult(Site site, IMT imt, double logMean, double stdDev) {
			Map<IMT, double[]> siteGMs = gms.get(site);
			if (siteGMs == null) {
				synchronized (gms) {
					siteGMs = gms.get(site);
					if (siteGMs == null) {
						siteGMs = new ConcurrentHashMap<>();
						gms.put(site, siteGMs);
					}
				}
			}
			double[] vals = { logMean, stdDev };
			siteGMs.putIfAbsent(imt, vals);
		}
		
		public void setDistances(Site site, double distanceRup, double distanceJB) {
			distanceJBs.put(site, distanceJB);
			distanceRups.put(site, distanceRup);
		}
		
		@Override
		public double getLogMean(Site site, IMT imt) {
			return gms.get(site).get(imt)[0];
		}
		
		@Override
		public double getStdDev(Site site, IMT imt) {
			return gms.get(site).get(imt)[1];
		}
		
		@Override
		public double getDistanceRup(Site site) {
			return distanceRups.get(site);
		}
		
		@Override
		public double getDistanceJB(Site site) {
			return distanceJBs.get(site);
		}

		@Override
		public boolean hasSite(Site site) {
			return gms.containsKey(site);
		}
		
		@Override
		public boolean isComputed(Site site, IMT imt) {
			Map<IMT, double[]> siteGMs = gms.get(site);
			return siteGMs != null && siteGMs.containsKey(imt);
		}
		
		public Set<IMT> getIMTs(Site site) {
			Map<IMT, double[]> imtMap = gms.get(site);
			return imtMap == null ? null : imtMap.keySet();
		}
		
		public void calculate(ScalarIMR gmpe, IMT... imts) {
			calculate(gmpe, getApplicableSites(), imts);
		}
		
		public void calculate(ScalarIMR gmpe, Site site, IMT... imts) {
			List<Site> sites = new ArrayList<>();
			sites.add(site);
			calculate(gmpe, sites, imts);
		}
		
		public void calculate(ScalarIMR gmpe, Collection<Site> sites, IMT... imts) {
			EqkRupture rup = getGMPERupture();
			gmpe.setEqkRupture(rup);
			for (Site site : sites) {
				gmpe.setSite(site);
				RuptureSurface surf = rup.getRuptureSurface();
				double distanceRup = surf.getDistanceRup(site.getLocation());
				double distanceJB = surf.getDistanceJB(site.getLocation());
				setDistances(site, distanceRup, distanceJB);
				for (IMT imt : imts) {
					imt.setIMT(gmpe);
					double logMean = gmpe.getMean();
					double stdDev = gmpe.getStdDev();
					addResult(site, imt, logMean, stdDev);
				}
			}
		}
		
	}

}
