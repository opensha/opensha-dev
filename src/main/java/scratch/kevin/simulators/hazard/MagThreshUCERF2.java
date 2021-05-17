package scratch.kevin.simulators.hazard;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.RuptureSurface;

public class MagThreshUCERF2 extends MeanUCERF2 {
	
	private DoubleParameter magThreshParam;
	
	public MagThreshUCERF2() {
		this(6.5);
	}
	
	public MagThreshUCERF2(double minMag) {
		super();
		magThreshParam = new DoubleParameter("Minimum Magnitude", minMag);
		
		this.adjustableParams.addParameter(magThreshParam);
		this.getAdjustableParameterList();
	}

	@Override
	public void parameterChange(ParameterChangeEvent event) {
		super.parameterChange(event);
		// it can be removed on a probability model change, which triggers a list rebuild
		if (!adjustableParams.containsParameter(magThreshParam))
			adjustableParams.addParameter(magThreshParam);
	}

	@Override
	public ProbEqkSource getSource(int iSource) {
		return new MagThreshSource(super.getSource(iSource), magThreshParam.getValue());
	}
	
	private class MagThreshSource extends ProbEqkSource {
		
		private ProbEqkSource source;
		private List<Integer> rupIndexes;

		public MagThreshSource(ProbEqkSource source, double minMag) {
			this.source = source;
			this.isPoissonian = source.isPoissonianSource();
			
			rupIndexes = new ArrayList<>();
			for (int r=0; r<source.getNumRuptures(); r++)
				if (source.getRupture(r).getMag() >= minMag)
					rupIndexes.add(r);
		}

		@Override
		public LocationList getAllSourceLocs() {
			return source.getAllSourceLocs();
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return source.getSourceSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return source.getMinDistance(site);
		}

		@Override
		public int getNumRuptures() {
			return rupIndexes.size();
		}

		@Override
		public ProbEqkRupture getRupture(int nRupture) {
			return source.getRupture(rupIndexes.get(nRupture));
		}
		
	}

	@Override
	public ArrayList<ProbEqkSource> getSourceList() {
		ArrayList<ProbEqkSource> sources = new ArrayList<>();
		for (int i=0; i<getNumSources(); i++)
			sources.add(getSource(i));
		return sources;
	}
	
	public static void main(String[] args) {
		MagThreshUCERF2 erf = new MagThreshUCERF2(6.876);
		erf.updateForecast();
		for (Parameter<?> param : erf.getAdjustableParameterList())
			System.out.println(param.getName()+": "+param.getValue());
		
		System.out.println();
	}

}
