package scratch.kevin.simCompare;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.ScalarIMR;
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
	
	public abstract double getLogMean(Site site, double period);
	
	public abstract double getStdDev(Site site, double period);
	
	public abstract double getDistanceRup(Site site);
	
	public abstract double getDistanceJB(Site site);
	
	public abstract boolean hasSite(Site site);
	
	public abstract boolean isComputed(Site site, double period);
	
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
		
		private Table<Site, Double, Double> logMeans;
		private Table<Site, Double, Double> stdDevs;
		private Map<Site, Double> distanceRups;
		private Map<Site, Double> distanceJBs;

		public Cached(Y rupture) {
			super(rupture);
			
			logMeans = HashBasedTable.create();
			stdDevs = HashBasedTable.create();
			distanceRups = new HashMap<>();
			distanceJBs = new HashMap<>();
		}
		
		public void addResult(Site site, double period, double logMean, double stdDev) {
			logMeans.put(site, period, logMean);
			stdDevs.put(site, period, stdDev);
		}
		
		public void setDistances(Site site, double distanceRup, double distanceJB) {
			distanceJBs.put(site, distanceJB);
			distanceRups.put(site, distanceRup);
		}
		
		@Override
		public double getLogMean(Site site, double period) {
			return logMeans.get(site, period);
		}
		
		@Override
		public double getStdDev(Site site, double period) {
			return stdDevs.get(site, period);
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
			return logMeans.containsRow(site);
		}
		
		@Override
		public boolean isComputed(Site site, double period) {
			return logMeans.contains(site, period);
		}
		
		public Set<Double> getPeriods(Site site) {
			return logMeans.row(site).keySet();
		}
		
		public void calculate(ScalarIMR gmpe, double... periods) {
			calculate(gmpe, getApplicableSites(), periods);
		}
		
		public void calculate(ScalarIMR gmpe, Site site, double... periods) {
			List<Site> sites = new ArrayList<>();
			sites.add(site);
			calculate(gmpe, sites, periods);
		}
		
		public void calculate(ScalarIMR gmpe, Collection<Site> sites, double... periods) {
			EqkRupture rup = getGMPERupture();
			gmpe.setEqkRupture(rup);
			gmpe.setIntensityMeasure(SA_Param.NAME);
			for (Site site : sites) {
				gmpe.setSite(site);
				RuptureSurface surf = rup.getRuptureSurface();
				double distanceRup = surf.getDistanceRup(site.getLocation());
				double distanceJB = surf.getDistanceJB(site.getLocation());
				setDistances(site, distanceRup, distanceJB);
				for (double period : periods) {
					SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
					double logMean = gmpe.getMean();
					double stdDev = gmpe.getStdDev();
					addResult(site, period, logMean, stdDev);
				}
			}
		}
		
	}

}
