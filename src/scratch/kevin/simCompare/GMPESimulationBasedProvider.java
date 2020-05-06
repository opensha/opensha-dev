package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.imr.param.IntensityMeasureParams.DurationTimeInterval;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

public class GMPESimulationBasedProvider<E> implements SimulationRotDProvider<E> {
	
	private SimulationRotDProvider<E> simProv;
	private Map<E, RuptureComparison<E>> compsMap;
	private String name;
	private IMT[] imts;
	private IMT[] saIMTs;
	private double[] saPeriods;
	private boolean hasPGV = false;
	private boolean hasDur = false;
	
	private NormalDistribution stdNorm;

	private Table<Site, E, DiscretizedFunc[]> cache;
	private Table<Site, E, double[]> pgvCache;
	
	public GMPESimulationBasedProvider(SimulationRotDProvider<E> simProv, List<? extends RuptureComparison<E>> comps,
			String name, IMT[] imts) {
		this.name = name;
		this.simProv = simProv;
		this.imts = imts;
		List<Double> saPeriods = new ArrayList<>();
		List<IMT> saIMTs = new ArrayList<>();
		for (IMT imt : imts) {
			if (imt.getParamName().equals(SA_Param.NAME)) {
				saPeriods.add(imt.getPeriod());
				saIMTs.add(imt);
			} else if (imt == IMT.PGV) {
				hasPGV = true;
			}
		}
		this.saPeriods = Doubles.toArray(saPeriods);
		this.saIMTs = saIMTs.toArray(new IMT[0]);
		compsMap = new HashMap<>();
		for (RuptureComparison<E> comp : comps)
			compsMap.put(comp.getRupture(), comp);
		stdNorm = new NormalDistribution(new Well19937c(comps.size()), 0d, 1d);
		cache = HashBasedTable.create();
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
			double[] vals = new double[imts.length];
			for (int p=0; p<saPeriods.length; p++) {
				double logMean = comp.getLogMean(site, saIMTs[p]);
				double stdDev = comp.getStdDev(site, saIMTs[p]);
				double sample = stdNorm.sample();
				vals[p] = Math.exp(logMean + stdDev*sample);
			}
			funcs[i] = new LightFixedXFunc(saPeriods, vals);
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
	public double getPGV(Site site, E rupture, int index) throws IOException {
		if (pgvCache.contains(site, rupture))
			return pgvCache.get(site, rupture)[index];
		double[] pgvs = new double[getNumSimulations(site, rupture)];
		RuptureComparison<E> comp = compsMap.get(rupture);
		Preconditions.checkNotNull(comp, "Comp not found for rupture: %s", rupture);
		for (int i=0; i<pgvs.length; i++) {
			double logMean = comp.getLogMean(site, IMT.PGV);
			double stdDev = comp.getStdDev(site, IMT.PGV);
			double sample = stdNorm.sample();
			pgvs[i] = Math.exp(logMean + stdDev*sample);
		}
		pgvCache.put(site, rupture, pgvs);
		return pgvs[index];
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
	public boolean hasPGV() {
		return hasPGV;
	}

	@Override
	public double getAnnualRate(E rupture) {
		return simProv.getAnnualRate(rupture);
	}

	@Override
	public double getMinimumCurvePlotRate(Site site) {
		return simProv.getMinimumCurvePlotRate(site);
	}

	@Override
	public double getMagnitude(E rupture) {
		return simProv.getMagnitude(rupture);
	}

	@Override
	public Location getHypocenter(E rupture, int index) {
		return simProv.getHypocenter(rupture, index);
	}

	@Override
	public double getDuration(Site site, E rupture, DurationTimeInterval interval, int index) throws IOException {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public boolean hasDurations() {
		return hasDur;
	}

}
