package scratch.kevin.simCompare;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.opensha.commons.data.Site;

public abstract class RuptureComparisonFilter<E> {
	
	public abstract boolean matches(RuptureComparison<E> comp, Site site);
	
	public <Y extends RuptureComparison<E>> List<Y> getMatches(Collection<Y> comps, Site site) {
		List<Y> matches = new ArrayList<>();
		for (Y comp : comps)
			if (matches(comp, site))
				matches.add(comp);
		return matches;
	}
	
	public static class DistJBFilter<E> extends RuptureComparisonFilter<E> {
		
		private double min, max;

		public DistJBFilter(double min, double max) {
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			if (!comp.hasSite(site))
				return false;
			double dist = comp.getDistanceJB(site);
			return dist >= min && dist <= max;
		}
		
	}
	
	public static class DistRupFilter<E> extends RuptureComparisonFilter<E> {
		
		private double min, max;

		public DistRupFilter(double min, double max) {
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			if (!comp.hasSite(site))
				return false;
			double dist = comp.getDistanceJB(site);
			return dist >= min && dist <= max;
		}
		
	}
	
	public static class HasSiteFilter<E> extends RuptureComparisonFilter<E> {

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			return comp.hasSite(site);
		}
		
	}
	
	public static class MatchesSiteFilter<E> extends RuptureComparisonFilter<E> {
		
		private Site site;

		public MatchesSiteFilter(Site site) {
			this.site = site;
		}

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			return this.site.equals(site) && comp.hasSite(site);
		}
		
	}
	
	public static class MagFilter<E> extends RuptureComparisonFilter<E> {
		
		private double min, max;

		public MagFilter(double min, double max) {
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			double mag = comp.getMagnitude();
			return mag >= min && mag <= max;
		}
		
	}
	
	public static class AcceptAllFilter<E> extends RuptureComparisonFilter<E> {

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			return true;
		}
		
	}

}
