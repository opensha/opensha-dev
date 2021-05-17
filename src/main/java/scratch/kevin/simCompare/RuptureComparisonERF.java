package scratch.kevin.simCompare;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.geo.LocationList;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

public class RuptureComparisonERF<E> extends AbstractERF {
	
	private List<? extends RuptureComparison<E>> comps;
	private List<CompProbEqkSource> sources;

	public RuptureComparisonERF(List<? extends RuptureComparison<E>> comps) {
		this.comps = comps;
		timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		timeSpan.setDuration(1d);
	}

	@Override
	public int getNumSources() {
		return comps.size();
	}

	@Override
	public ProbEqkSource getSource(int idx) {
		return sources.get(idx);
	}

	@Override
	public void updateForecast() {
		sources = new ArrayList<>();
		double duration = timeSpan.getDuration();
		for (RuptureComparison<E> comp : comps)
			sources.add(new CompProbEqkSource(comp, duration));
	}

	@Override
	public String getName() {
		return "Simulation Wrapper ERF";
	}
	
	public class CompProbEqkRupture extends ProbEqkRupture {
		
		private RuptureComparison<E> comp;
		
		public CompProbEqkRupture(RuptureComparison<E> comp, EqkRupture rup, double prob) {
			super(rup.getMag(), rup.getAveRake(), prob, rup.getRuptureSurface(), rup.getHypocenterLocation());
			this.comp = comp;
		}
		
		public RuptureComparison<E> getRuptureComparison() {
			return comp;
		}
		
	}
	
	private class CompProbEqkSource extends ProbEqkSource {
		
		private CompProbEqkRupture rup;
		
		public CompProbEqkSource(RuptureComparison<E> comp, double duration) {
			double prob = 1 - Math.exp(-comp.getAnnualRate()*duration);
			rup = new CompProbEqkRupture(comp, comp.getGMPERupture(), prob);
		}

		@Override
		public LocationList getAllSourceLocs() {
			return rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return rup.getRuptureSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return rup.getRuptureSurface().getDistanceJB(site.getLocation());
		}

		@Override
		public int getNumRuptures() {
			return 1;
		}

		@Override
		public CompProbEqkRupture getRupture(int nRupture) {
			Preconditions.checkState(nRupture == 0);
			return rup;
		}
		
	}

}
