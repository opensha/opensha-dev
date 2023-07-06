package scratch.kevin.ruptureDirection;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class PreferredDirectionERF extends AbstractERF {
	
	private FaultSystemSolutionERF erf;
	private int numHypos;
	private List<PreferredDirection> directions;
	
	private ConcurrentMap<Integer, ProbEqkSource> overrideSourceCache = new ConcurrentHashMap<>();

	public PreferredDirectionERF(FaultSystemSolutionERF erf, int numHypos, PreferredDirection... directions) {
		this(erf, numHypos, List.of(directions));
	}
	
	public PreferredDirectionERF(FaultSystemSolutionERF erf, int numHypos, List<PreferredDirection> directions) {
		this.erf = erf;
		this.numHypos = numHypos;
		this.directions = directions;
	}

	@Override
	public int getNumSources() {
		return erf.getNumSources();
	}
	
	public void cacheSourcesParallel() {
		cacheSourcesParallel(FaultSysTools.defaultNumThreads());
	}
	
	public void cacheSourcesParallel(int threads) {
		HashSet<Integer> ruptures = new HashSet<>();
		for (PreferredDirection direction : directions)
			ruptures.addAll(direction.getAffectedRuptures());
		
		ExecutorService exec = Executors.newFixedThreadPool(threads);
		
		List<Integer> sourceIDs = new ArrayList<>();
		List<Future<ProbEqkSource>> sourceFutures = new ArrayList<>();
		
		int numFaultSources = erf.getNumFaultSystemSources();
		for (int sourceID=0; sourceID<numFaultSources; sourceID++) {
			int rupIndex = erf.getFltSysRupIndexForSource(sourceID);
			if (ruptures.contains(rupIndex)) {
				sourceIDs.add(sourceID);
				ProbEqkSource origSource = erf.getSource(sourceID);
				sourceFutures.add(exec.submit(new Callable<ProbEqkSource>() {

					@Override
					public ProbEqkSource call() throws Exception {
						return buildModifiedSource(origSource, rupIndex);
					}
				}));
			}
		}
		
		System.out.println("Building hypocenters for "+sourceFutures.size()+" source futures with "+threads+" threads");
		int printMod = 100;
		DecimalFormat twoDigits = new DecimalFormat("0.00");
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<sourceFutures.size(); i++) {
			if (i % printMod == 0) {
				double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
				System.out.println("Building for rupture "+i+"/"+sourceFutures.size()+" (rate="+twoDigits.format((double)i/secs)+" /s)");
			}
			if (i >= printMod*10 && printMod < 1000)
				printMod *= 10;
			int sourceID = sourceIDs.get(i);
			try {
				ProbEqkSource source = sourceFutures.get(i).get();
				overrideSourceCache.put(sourceID, source);
			} catch (InterruptedException | ExecutionException e) {
				exec.shutdown();
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		watch.stop();
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Done building hypocenters for "+sourceFutures.size()+" ruptures; took "+twoDigits.format(secs)+" s");
		
		exec.shutdown();
	}

	@Override
	public ProbEqkSource getSource(int sourceID) {
		// see if already cached
		ProbEqkSource override = overrideSourceCache.get(sourceID);
		if (override != null)
			return override;
		
		// see if it's a fault system source
		if (sourceID < erf.getNumFaultSystemSources()) {
			ProbEqkSource origSource = erf.getSource(sourceID);
			int rupIndex = erf.getFltSysRupIndexForSource(sourceID);
			ProbEqkSource modSource = buildModifiedSource(origSource, rupIndex);
			if (modSource != null) {
				overrideSourceCache.putIfAbsent(sourceID, modSource);
				return modSource;
			}
		}
		
		return erf.getSource(sourceID);
	}
	
	private ProbEqkSource buildModifiedSource(ProbEqkSource origSource, int rupIndex) {
		List<PreferredDirection> matches = null;
		for (PreferredDirection direction : directions) {
			if (direction.affectsRupture(rupIndex)) {
				if (matches == null)
					matches = new ArrayList<>();
				matches.add(direction);
			}
		}
		if (matches != null) {
			// this rupture is governed by preferred directions

			List<List<Location>> rupHypos = new ArrayList<>();
			List<Double> rateFracts = new ArrayList<>();
			double sumRateFracts = 0d;
			for (PreferredDirection direction : matches) {
				List<Location> hypos = direction.drawPreferredDirectionHypos(rupIndex, numHypos);
				if (hypos == null) {
					System.err.println("Warning: couldn't build hypocenters for "+direction.getName()+" and rupture "+rupIndex);
					continue;
				}
				rupHypos.add(hypos);
				double rateFract = direction.getFractNuclOnTargetSects(rupIndex)*direction.getTargetFractInPreferred();
				Preconditions.checkState(rateFract > 0, "Direction %s reports that it affects %s, but rate fract is %s",
						direction.getName(), rupIndex, rateFract);
				rateFracts.add(rateFract);
				sumRateFracts += rateFract;
			}
			Preconditions.checkState(sumRateFracts <= 1d,
					"Sum rate fract exceeds 1? Directions must overlap: %s", sumRateFracts);

			double duration = getTimeSpan().getDuration();

			List<ProbEqkRupture> subRuptures = new ArrayList<>();
			for (ProbEqkRupture origRup : origSource) {
				double origRate = origRup.getMeanAnnualRate(duration);

				if (sumRateFracts < 1d) {
					// add a version without a hypocenter to represent the fraction that should be uniform
					double leftoverRate = origRate * (1d - sumRateFracts);
					subRuptures.add(cloneNewRateHypo(origRup, leftoverRate, duration, null));
				}
				// now add in preferred direction hypocenters
				for (int i=0; i<rupHypos.size(); i++) {
					List<Location> hypos = rupHypos.get(i);
					double rateFract = rateFracts.get(i);

					double rateEach = origRate*rateFract/(double)hypos.size();
					for (Location hypo : hypos)
						subRuptures.add(cloneNewRateHypo(origRup, rateEach, duration, hypo));
				}
			}

			ProbEqkSource source = new RupListSource(origSource, subRuptures);
			return source;
		}
		return null;
	}
	
	private static class RupListSource extends ProbEqkSource {
		
		private ProbEqkSource origSource;
		private List<ProbEqkRupture> rups;

		public RupListSource(ProbEqkSource origSource, List<ProbEqkRupture> rups) {
			this.origSource = origSource;
			this.rups = rups;
		}

		@Override
		public LocationList getAllSourceLocs() {
			return origSource.getAllSourceLocs();
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return origSource.getSourceSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return origSource.getMinDistance(site);
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
	
	private ProbEqkRupture cloneNewRateHypo(ProbEqkRupture origRup, double newRate, double duration, Location hypo) {
		double prob = 1-Math.exp(-newRate*duration);
		return new ProbEqkRupture(origRup.getMag(), origRup.getAveRake(), prob, origRup.getRuptureSurface(), hypo);
	}

	@Override
	public void updateForecast() {
		overrideSourceCache.clear();
		erf.updateForecast();
	}

	@Override
	public String getName() {
		return "Preferred Dirction ERF";
	}

	@Override
	public Parameter getParameter(String paramName) {
		return erf.getParameter(paramName);
	}

	@Override
	public TimeSpan getTimeSpan() {
		return erf.getTimeSpan();
	}

	@Override
	public ParameterList getAdjustableParameterList() {
		return erf.getAdjustableParameterList();
	}

}
