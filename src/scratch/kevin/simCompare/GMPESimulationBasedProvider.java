package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

public class GMPESimulationBasedProvider<E> implements SimulationRotDProvider<E> {
	
	private SimulationRotDProvider<E> simProv;
	private Map<E, RuptureComparison<E>> compsMap;
	private String name;
	private double[] periods;
	
	private NormalDistribution stdNorm;
	
	private Table<Site, E, DiscretizedFunc[]> cache;
	
	public GMPESimulationBasedProvider(SimulationRotDProvider<E> simProv, List<? extends RuptureComparison<E>> comps,
			String name, double[] periods) {
		this.name = name;
		this.simProv = simProv;
		this.periods = periods;
		compsMap = new HashMap<>();
		for (RuptureComparison<E> comp : comps)
			compsMap.put(comp.getRupture(), comp);
		stdNorm = new NormalDistribution(new Well19937c(comps.size()), 0d, 1d);
		cache = HashBasedTable.create();
	}
	
	public void clearCache() {
		cache.clear();
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public synchronized DiscretizedFunc getRotD50(Site site, E rupture, int index) throws IOException {
		if (cache.contains(site, rupture))
			return cache.get(site, rupture)[index];
		DiscretizedFunc[] funcs = new DiscretizedFunc[getNumSimulations(site, rupture)];
		RuptureComparison<E> comp = compsMap.get(rupture);
		Preconditions.checkNotNull(comp, "Comp not found for rupture: %s", rupture);
		for (int i=0; i<funcs.length; i++) {
			double[] vals = new double[periods.length];
			for (int p=0; p<periods.length; p++) {
				double logMean = comp.getLogMean(site, periods[p]);
				double stdDev = comp.getStdDev(site, periods[p]);
				double sample = stdNorm.sample();
				vals[p] = Math.exp(logMean + stdDev*sample);
			}
			funcs[i] = new LightFixedXFunc(periods, vals);
		}
		cache.put(site, rupture, funcs);
		return funcs[index];
	}

	@Override
	public DiscretizedFunc getRotD100(Site site, E rupture, int index) throws IOException {
		return null;
	}

	@Override
	public DiscretizedFunc[] getRotD(Site site, E rupture, int index) throws IOException {
		return null;
	}

	@Override
	public DiscretizedFunc getRotDRatio(Site site, E rupture, int index) throws IOException {
		return null;
	}

	@Override
	public int getNumSimulations(Site site, E rupture) {
		return simProv.getNumSimulations(site, rupture);
	}

	@Override
	public Collection<E> getRupturesForSite(Site site) {
		return simProv.getRupturesForSite(site);
	}

	@Override
	public boolean hasRotD50() {
		return true;
	}

	@Override
	public boolean hasRotD100() {
		return false;
	}

	@Override
	public double getAnnualRate(E rupture) {
		return simProv.getAnnualRate(rupture);
	}

	@Override
	public double getMinimumCurvePlotRate() {
		return simProv.getMinimumCurvePlotRate();
	}

	@Override
	public double getMagnitude(E rupture) {
		return simProv.getMagnitude(rupture);
	}

}
