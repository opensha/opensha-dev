package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.util.Precision;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.TrueMeanRuptureMappings;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.util.GriddedFiniteRuptureSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Table;

public class TrueMeanERFDemo {

	public static void main(String[] args) throws IOException {
		File dataDir = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3-gridded_rebuild");
		
		System.out.println("Loading solutions");
		
		// download these files from:
		// https://data.opensha.org/ftp/kmilner/markdown/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3-gridded_rebuild/
		
		// true mean solution file that contains data and mappings for each unique on-fault rupture across each branch
		FaultSystemSolution trueMeanSol = FaultSystemSolution.load(
				new File(dataDir, "true_mean_solution.zip"));
		// downsampled model with 1k branches (fault and gridded
		SolutionLogicTree slt = SolutionLogicTree.load(
				new File(dataDir, "results_full_gridded_downsampled_1k.zip"));
		// this contains the branch-averaged gridded seismicity model, including tweaks to active/stable weighting
		// made in the revised model
		FaultSystemSolution baSol = FaultSystemSolution.load(
				new File(dataDir, "results_WUS_FM_v3_branch_averaged_gridded.zip"));
		
		// use the branch-averaged gridded seismicity model
		GridSourceList baGridList = baSol.requireModule(GridSourceList.class);
		trueMeanSol.setGridSourceProvider(baGridList);
		// mapping between grid location, TRT, and source index in the grid sources
		Table<Integer, TectonicRegionType, Integer> baLocToSourceTable = HashBasedTable.create();
		for (int s=0; s<baGridList.getNumSources(); s++) {
			int locIndex = baGridList.getLocationIndexForSource(s);
			TectonicRegionType trt = baGridList.tectonicRegionTypeForSourceIndex(s);
			baLocToSourceTable.put(locIndex, trt, s);
		}
		
		System.out.println("Building true mean ERF");
		BaseFaultSystemSolutionERF trueMeanERF = new BaseFaultSystemSolutionERF();
		trueMeanERF.setSolution(trueMeanSol);
		// make sure it uses the branch-specific magnitudes
		trueMeanERF.setParameter(UseRupMFDsParam.NAME, true);
		// make sure gridded seismicity is enabled
		trueMeanERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		// tell it to draw random finite ruptures for each gridded seismicity rupture
		// these will be reproducible, i.e., each time will draw the same one for each rupture
		GriddedSeismicitySettings gridSettings = GriddedSeismicitySettings.DEFAULT
				.forSurfaceType(BackgroundRupType.FINITE)
				.forFiniteRuptureSettings(GriddedFiniteRuptureSettings.DEFAULT_SINGLE)
				// ignore gridded ruptures below this magnitude (if they happen to be attached)
				.forMinimumMagnitude(5d);
		trueMeanERF.setGriddedSeismicitySettings(gridSettings);
		// this speeds up grid source calculations/mappings, but uses more memory; disable if you run into trouble
		trueMeanERF.setCacheGridSources(true);
		// make sure forecast duration is 1 year
		trueMeanERF.getTimeSpan().setDuration(1d);
		trueMeanERF.updateForecast();
		
		// do whatever calculation you want on the true mean ERF here
		// cache data for each source & rupture combination
		for (int sourceIndex=0; sourceIndex<trueMeanERF.getNumSources(); sourceIndex++) {
			ProbEqkSource source = trueMeanERF.getSource(sourceIndex);
			for (int rupIndex=0; rupIndex<source.getNumRuptures(); rupIndex++) {
				ProbEqkRupture rup = source.getRupture(rupIndex);
				// TODO: do calculation
				// store results indexed by [sourceIndex][rupIndex]
			}
		}
		
		// this stores the data for finding branch-specific ruptures in the true mean solution
		TrueMeanRuptureMappings branchMappings = trueMeanSol.getRupSet().requireModule(TrueMeanRuptureMappings.class);
		
		// this will suppress prints on each branch
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		// iterate over the solution logic tree
		LogicTree<?> logicTree = slt.getLogicTree();
		// this will be used to preload the next solution in parallel because most of the time in this loop is spent
		// just waiting on I/O to load the solutions for each branch
		CompletableFuture<FaultSystemSolution> preloadFuture = null;
		Stopwatch watch = Stopwatch.createStarted();
		for (int b=0; b<logicTree.size(); b++) {
			LogicTreeBranch<?> branch = logicTree.getBranch(b);
			// will usually be the same for all branches if downsampled, but a branch can be drawn multiple times which
			// increases its weight
			double branchWeight = logicTree.getBranchWeight(branch);
			
			System.out.println("Handling branch "+b+"/"+logicTree.size()+": "+branch+" (wt="+(float)branchWeight+")");
			System.out.println("\tLoading solution...");
			FaultSystemSolution branchSol;
			if (preloadFuture == null)
				branchSol = slt.forBranch(branch);
			else
				branchSol = preloadFuture.join();
			FaultSystemRupSet branchRupSet = branchSol.getRupSet();
			GridSourceList branchGridList = branchSol.requireModule(GridSourceList.class);
			System.out.println("\t\tSolution is ready");
			
			// start the next load in parallel
			int nextIndex = b+1;
			if (nextIndex < logicTree.size()) {
				preloadFuture = CompletableFuture.supplyAsync(()->{
					FaultSystemSolution futureSol;
					try {
						futureSol = slt.forBranch(logicTree.getBranch(nextIndex));
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					// have it proload the grid sources
					futureSol.getGridSourceProvider();
					return futureSol;
				});
			}
			
			System.out.println("\tValidating FSS rupture mappings...");
			
			// this gives, for each on-fault rupture in my branch, the corresponding index in the true mean solution
			int[] branchRupIndexes = branchMappings.getRuptureMappings(branch);
			
			int numMappedFSS = 0;
			double mappedFSSRate = 0d;
			
			// find mappings for each fault-based rupture
			for (int i=0; i<branchRupIndexes.length; i++) {
				// this is the rate for that rupture on this branch
				// you can convert to a Poisson probability
				double rupRate = branchSol.getRateForRup(i);
				
				int trueMeanIndex = branchRupIndexes[i];
				if (trueMeanIndex < 0) {
					// not mapped, must be zero rate
					Preconditions.checkState(rupRate == 0d);
					continue;
				} else if (rupRate == 0d) {
					// mapped, but still zero
					continue;
				}
				
				// source index in the true mean ERF for this rupture
				int trueMeanSourceIndex = trueMeanERF.getSrcIndexForFltSysRup(trueMeanIndex);
				ProbEqkSource trueMeanSource = trueMeanERF.getSource(trueMeanSourceIndex);
				// now we need to find the one with the correct magnitude
				int trueMeanRupIndexInSource = -1;
				double rupMag = branchRupSet.getMagForRup(i);
				for (int rupIndex=0; rupIndex<trueMeanSource.getNumRuptures(); rupIndex++) {
					ProbEqkRupture rup = trueMeanSource.getRupture(rupIndex);
					if (Precision.equals(rup.getMag(), rupMag)) {
						Preconditions.checkState(trueMeanRupIndexInSource < 0, "Duplicate match found?");
						trueMeanRupIndexInSource = rupIndex;
					}
				}
				
				Preconditions.checkState(trueMeanRupIndexInSource >= 0,
						"No matching rupture found for branch %s FSS rupture %s with rate=%s",
						b, i, rupRate);
				
				// now we know that the true mean ERF's rupture for this branch is at:
				// 		sourceIndex=trueMeanSourceIndex, ruptureIndex=trueMeanRupIndexInSource
				// you can pre-calculate and store data using those the ERF source and rupture indexes
				
				numMappedFSS++;
				mappedFSSRate += rupRate;
			}
			System.out.println("\t\tMapped "+numMappedFSS+"/"+branchRupIndexes.length+" FSS ruptures with total rate "+(float)mappedFSSRate);
			
			System.out.println("\tValidating gridded rupture mappings...");
			// now find mappings for each gridded seismicity rupture
			Preconditions.checkState(branchGridList.getNumLocations() == baGridList.getNumLocations(),
					"Branch %s grid list has %s locations, but BA has %s", b, branchGridList.getNumLocations(), baGridList.getNumLocations());
			Preconditions.checkState(branchGridList.getNumSources() <= baGridList.getNumSources(),
					"Branch %s grid list has %s sources, but BA has %s", b, branchGridList.getNumSources(), baGridList.getNumSources());
			int numMappedGrid = 0;
			double mappedGridRate = 0d;
			// this will be used to speed up rupture mapping; most sources will have the same rupture lists over and
			// over again, so doing the index source for each rupture is wasteful; this keeps track of the indexes
			// from the prior source and reuses it if they agree
			int[] prevMappings = null;
			for (int s=0; s<branchGridList.getNumSources(); s++) {
				int gridIndex = branchGridList.getLocationIndexForSource(s);
				TectonicRegionType trt = branchGridList.tectonicRegionTypeForSourceIndex(s);
				Preconditions.checkState(baLocToSourceTable.contains(gridIndex, trt));
				// true mean ERF source index is the number of FSS sources, plus source index in the grid list
				int trueMeanSourceIndex = trueMeanERF.getNumFaultSystemSources()+baLocToSourceTable.get(gridIndex, trt);
				ProbEqkSource source = trueMeanERF.getSource(trueMeanSourceIndex);
				Preconditions.checkState(source instanceof PointSource,
						"True mean source %s isn't a point source? %s", trueMeanSourceIndex, source.getClass());
				PointSource ptSource = (PointSource)source;
				int numSourceRups = ptSource.getNumRuptures();
				List<ProbEqkRupture> sourceRups = ptSource.getRuptureList();
				
				Location loc = branchGridList.getLocation(gridIndex);
				Preconditions.checkState(LocationUtils.areSimilar(loc, ptSource.getLocation()));
				List<GriddedRupture> gridRups = branchGridList.getRuptures(branchGridList.tectonicRegionTypeForSourceIndex(s), gridIndex);
				
				int[] myMappings = null;
				if (prevMappings == null || gridRups.size() != prevMappings.length) {
					// definitely not a match with the prior source, track our mappings in case they can be used in the
					// future;
					prevMappings = null;
					myMappings = new int[gridRups.size()];
				}
				for (int g=0; g<gridRups.size(); g++) {
					GriddedRupture gridRup = gridRups.get(g);
					if (gridRup.properties.magnitude < gridSettings.minimumMagnitude)
						continue;
					// find a matching rupture in the ERF source
					int trueMeanRupIndexInSource;
					if (prevMappings != null && prevMappings[g] < numSourceRups && isGridMatch(gridRup, sourceRups.get(prevMappings[g]))) {
						// already matched by a prior source, just use that
						trueMeanRupIndexInSource = prevMappings[g];
					} else {
						// need to search for it (slow)
						prevMappings = null;
						trueMeanRupIndexInSource = -1;
						for (int rupIndex=0; rupIndex<numSourceRups; rupIndex++) {
							ProbEqkRupture rup = sourceRups.get(rupIndex);
							if (isGridMatch(gridRup, rup)) {
								Preconditions.checkState(trueMeanRupIndexInSource < 0, "Duplicate match found?");
								trueMeanRupIndexInSource = rupIndex;
							}
						}
					}
					Preconditions.checkState(trueMeanRupIndexInSource >= 0,
							"No matching rupture found for branch %s grid source %s rupture: %s",
							b, s, gridRup);
					if (myMappings != null)
						myMappings[g] = trueMeanRupIndexInSource;
					
					// this is the rate for that rupture on this branch
					// you can convert to a Poisson probability
					double rupRate = gridRup.rate;
					
					// now we know that the true mean ERF's rupture for this branch is at:
					// 		sourceIndex=trueMeanSourceIndex, ruptureIndex=trueMeanRupIndexInSource
					// you can pre-calculate and store data using those the ERF source and rupture indexes
					
					numMappedGrid++;
					mappedGridRate += rupRate;
				}
				if (myMappings != null)
					prevMappings = myMappings;
			}
			System.out.println("\t\tMapped "+numMappedGrid+" gridded ruptures across "+baGridList.getNumSources()
					+" sources with total rate "+(float)mappedGridRate);
			System.out.println("\tDone; "+timeLeftStr(watch, b+1, logicTree.size())+" remaining");
			
		}
		System.out.println("DONE");
		watch.stop();
		System.out.println("Took "+timeStr(watch));
	}
	
	private static DecimalFormat twoDF = new DecimalFormat("0.00");
	
	private static String timeStr(Stopwatch watch) {
		return timeStr(watch.elapsed(TimeUnit.MILLISECONDS)/1000d);
	}
	
	private static String timeLeftStr(Stopwatch watch, int done, int tot) {
		double totSecs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		double secsPer = totSecs/(double)done;
		double secsLeft = secsPer * (tot-done);
		return timeStr(secsLeft);
	}
	
	private static String timeStr(double secs) {
		double mins = secs/60d;
		double hours = mins/60d;
		if (hours > 1.5)
			return twoDF.format(hours)+" h";
		else if (mins > 1.5)
			return twoDF.format(mins)+" m";
		return twoDF.format(secs)+" s";
	}
	
	private static boolean isGridMatch(GriddedRupture gridRup, ProbEqkRupture erfRup) {
		return Precision.equals(erfRup.getMag(), gridRup.properties.magnitude)
				&& Precision.equals(erfRup.getAveRake(), gridRup.properties.rake);
	}

}
