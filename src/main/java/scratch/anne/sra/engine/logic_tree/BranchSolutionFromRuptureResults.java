package scratch.anne.sra.engine.logic_tree;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
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
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.TrueMeanRuptureMappings;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Table;

import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.io.readers.RuptureLossLookupReader;
import scratch.anne.sra.io.readers.TrueMeanSetupFromRuptureResults;
import scratch.anne.sra.io.readers.TrueMeanSetupFromRuptureResults.TrueMeanSetup;
import scratch.anne.sra.io.writers.LogicTreeBranchEalSummaryExporter;
import scratch.anne.sra.io.writers.LogicTreeBranchRuptureContributionExporter;

public class BranchSolutionFromRuptureResults {

	public static void main(String[] args) throws IOException {

		// ------------------------------------------------------------
		// Everything through updateForecast() (formerly lines 1-112:
		// loading trueMeanSol, filtering/reattaching the grid source
		// provider, building meanLocToSourceTable, constructing and
		// configuring trueMeanERF) is now driven by a previously
		// exported RuptureResultsExporter CSV instead of a hardcoded
		// construction sequence.
		// ------------------------------------------------------------
		Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\BERM\\results\\EAL\\true_mean_solution");
		Path inputFolder = Path.of("RES1-W1-MC_Porter-vuln_vs-360");
		
		
		// Cap the number of branches processed, for quick tests. Set to Integer.MAX_VALUE to run the full tree.
//		int maxBranches = 10;
		int maxBranches = Integer.MAX_VALUE;

		// defaults to the folder name, matching PortfolioExpectedLossDualCalcByCSV's
		// own fallback. Override this if that run used an explicit file_tag.
		String fileTag = inputFolder.getFileName().toString();
		// String fileTag = "some_other_tag"; // uncomment/edit if file_tag was set explicitly

		Path exportCsv = baseFolder.resolve(inputFolder).resolve(fileTag + "_eal-per-rup.csv");

		System.out.println("Rebuilding true mean ERF from exported HazardParameters: " + exportCsv);
		TrueMeanSetup setup = TrueMeanSetupFromRuptureResults.rebuild(exportCsv);

		FaultSystemSolution trueMeanSol = setup.trueMeanSol();
		BaseFaultSystemSolutionERF trueMeanERF = setup.trueMeanERF();
		GridSourceList meanGridList = setup.meanGridList();
		Table<Integer, TectonicRegionType, Integer> meanLocToSourceTable = setup.meanLocToSourceTable();
		float gridMinMag = setup.gridMinMag();

		System.out.println("True mean ERF ready: "
				+ trueMeanERF.getNumSources() + " sources, "
				+ trueMeanERF.getNumFaultSystemSources() + " fault-based");

		// ------------------------------------------------------------
		// Fixed per-rupture loss lookup: read ONCE from the exported
		// true-mean per-rupture CSV (ruptureLoss column only -- see
		// RuptureLossLookupReader). This never changes across branches;
		// only each branch's RATES change, per confirmed requirements.
		// likelihood/lossContribution/relativeContribution from that
		// export are intentionally ignored: likelihood is about to be
		// replaced by each branch's own rate, and the other two are
		// derived from likelihood*ruptureLoss, so neither is meaningful
		// input here.
		// ------------------------------------------------------------
		Map<RuptureKey, Double> ruptureLossLookup = RuptureLossLookupReader.read(exportCsv);
		System.out.println("Loaded fixed loss lookup for " + ruptureLossLookup.size() + " ruptures");

		// downsampled model with 1k branches (fault and gridded) -- this
		// part of the demo still needs its own source separately, since
		// it isn't part of what gets built/exported per-branch.
		File dataDir = new File(
				"C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\"
				+ "my_scratch\\logic_tree_reweighting\\2024_02_02-nshm23_branches-WUS_FM_v3-gridded_rebuild"
		);
		SolutionLogicTree slt = SolutionLogicTree.load(
				new File(dataDir, "results_full_gridded_downsampled_1k.zip"));

		// this stores the data for finding branch-specific ruptures in the true mean solution
		TrueMeanRuptureMappings branchMappings = trueMeanSol.getRupSet().requireModule(TrueMeanRuptureMappings.class);

		// this will suppress prints on each branch
		ModuleContainer.VERBOSE_DEFAULT = false;

		// iterate over the solution logic tree
		LogicTree<?> logicTree = slt.getLogicTree();
		
		int numBranchesToRun = Math.min(maxBranches, logicTree.size());
		if (numBranchesToRun < logicTree.size())
			System.out.println("TEST MODE: only running " + numBranchesToRun + "/" + logicTree.size() + " branches");
		
		String fileSuffix = numBranchesToRun < logicTree.size() ? "_subset-n-" + numBranchesToRun : "";
		Path branchEALsCsv = baseFolder.resolve(inputFolder).resolve(
		        fileTag + "_branch-eals" + fileSuffix + ".csv");
		
		Path branchContributionsParquet = baseFolder.resolve(inputFolder).resolve(fileTag + "_branch-rupture-contributions" + fileSuffix + ".parquet");

		// per-branch EAL results, indexed by branch index b. Fine as a
		// plain array up to millions of branches (a few MB at most);
		// swap for incremental disk writes later if parallelizing across
		// processes rather than threads within one process.
		double[] branchEALs = new double[numBranchesToRun];
		// per-branch total mapped rates, for the sanity check printed
		// during each branch (sum of annual occurrence rates across all
		// mapped ruptures in that category -- not a probability, not the
		// same as branchWeight).
		double[] branchFssRates = new double[numBranchesToRun];
		double[] branchGriddedRates = new double[numBranchesToRun];

		// this will be used to preload the next solution in parallel because most of the time in this loop is spent
		// just waiting on I/O to load the solutions for each branch
		CompletableFuture<FaultSystemSolution> preloadFuture = null;
		Stopwatch watch = Stopwatch.createStarted();
		
		try (LogicTreeBranchRuptureContributionExporter contribExporter = new LogicTreeBranchRuptureContributionExporter(branchContributionsParquet)) {
			for (int b=0; b<numBranchesToRun; b++) {
				LogicTreeBranch<?> branch = logicTree.getBranch(b);
				// will usually be the same for all branches if downsampled, but a branch can be drawn multiple times which
				// increases its weight
				double branchWeight = logicTree.getBranchWeight(branch);
	
				System.out.println("Handling branch "+b+"/"+numBranchesToRun+": "+branch+" (wt="+(float)branchWeight+")");
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
				if (nextIndex < numBranchesToRun) {
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
	
				// per-branch EAL accumulator: local to this branch/thread only,
				// no shared mutable state -- safe if this loop is parallelized
				// later (see BranchEALCalculator's class-level javadoc).
				BranchEALCalculator branchCalc = new BranchEALCalculator(ruptureLossLookup);
	
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
					// accumulate this rupture's contribution to the branch's EAL,
					// using the FIXED loss recorded for this rupture in the
					// true-mean run (see BranchEALCalculator)
					branchCalc.accumulateFault(trueMeanSourceIndex, trueMeanRupIndexInSource, rupRate);
					
		            Double rupLoss = ruptureLossLookup.get(new RuptureKey(trueMeanSourceIndex, trueMeanRupIndexInSource));
		            if (rupLoss != null)
		                contribExporter.write(b, trueMeanSourceIndex, trueMeanRupIndexInSource, rupLoss, rupRate);

	
					numMappedFSS++;
					mappedFSSRate += rupRate;
				}
				System.out.println("\t\tMapped "+numMappedFSS+"/"+branchRupIndexes.length+" FSS ruptures with total rate "+(float)mappedFSSRate);
				branchFssRates[b] = mappedFSSRate;
	
				System.out.println("\tValidating gridded rupture mappings...");
				// now find mappings for each gridded seismicity rupture
				Preconditions.checkState(branchGridList.getNumLocations() == meanGridList.getNumLocations(),
						"Branch %s grid list has %s locations, but mean has %s", b, branchGridList.getNumLocations(), meanGridList.getNumLocations());
				Preconditions.checkState(branchGridList.getNumSources() <= meanGridList.getNumSources(),
						"Branch %s grid list has %s sources, but mean has %s", b, branchGridList.getNumSources(), meanGridList.getNumSources());
				int numMappedGrid = 0;
				double mappedGridRate = 0d;
				// this will be used to speed up rupture mapping; most sources will have the same rupture lists over and
				// over again, so doing the index search for each rupture is wasteful; this keeps track of the indexes
				// from the prior source and reuses it if they agree, but if anything varies, it'll start fresh.
				int[] prevMappings = null;
				int numSourceMappingsReused = 0;
				for (int s=0; s<branchGridList.getNumSources(); s++) {
					int gridIndex = branchGridList.getLocationIndexForSource(s);
					TectonicRegionType trt = branchGridList.tectonicRegionTypeForSourceIndex(s);
					Preconditions.checkState(meanLocToSourceTable.contains(gridIndex, trt));
					// true mean ERF source index is the number of FSS sources, plus source index in the grid list
					int trueMeanSourceIndex = trueMeanERF.getNumFaultSystemSources()+meanLocToSourceTable.get(gridIndex, trt);
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
						if (gridRup.properties.magnitude < gridMinMag)
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
						branchCalc.accumulateGridded(trueMeanSourceIndex, trueMeanRupIndexInSource, rupRate);
						
			            Double rupLoss = ruptureLossLookup.get(new RuptureKey(trueMeanSourceIndex, trueMeanRupIndexInSource));
			            if (rupLoss != null)
			                contribExporter.write(b, trueMeanSourceIndex, trueMeanRupIndexInSource, rupLoss, rupRate);

	
						numMappedGrid++;
						mappedGridRate += rupRate;
					}
					if (myMappings != null)
						prevMappings = myMappings;
					else if (prevMappings != null)
						numSourceMappingsReused++;
				}
				System.out.println("\t\tMapped "+numMappedGrid+" gridded ruptures across "+branchGridList.getNumSources()
						+" sources with total rate "+(float)mappedGridRate+" ("+numSourceMappingsReused+" source mappings reused)");
				branchGriddedRates[b] = mappedGridRate;
	
				// EAL for this branch, recalculated from this branch's rates
				// against the FIXED true-mean loss lookup -- cheap, no
				// GMM/vulnerability recomputation needed (confirmed: only
				// rates vary per branch, not losses).
				branchEALs[b] = branchCalc.getEAL();
				System.out.println("\tBranch "+b+" EAL: "+branchEALs[b]);
	
				System.out.println("\tDone with branch "+b+"/"+numBranchesToRun+"; "+timeStr(watch)+" elapsed, estimated "
						+timeLeftStr(watch, b+1, numBranchesToRun)+" remaining");
	
			}
		}
		System.out.println("DONE");
		watch.stop();
		System.out.println("Took "+timeStr(watch));

		List<LogicTreeBranchEalSummaryExporter.BranchEalRow> rows = new java.util.ArrayList<>(numBranchesToRun);
		for (int b = 0; b < numBranchesToRun; b++) {
		    rows.add(new LogicTreeBranchEalSummaryExporter.BranchEalRow(
		            b,
		            logicTree.getBranchWeight(logicTree.getBranch(b)),
		            branchEALs[b],
		            branchFssRates[b],
		            branchGriddedRates[b],
		            logicTree.getBranch(b)));
		}
		LogicTreeBranchEalSummaryExporter.write(rows, logicTree, null, branchEALsCsv);
		System.out.println("Wrote branch EALs to: " + branchEALsCsv);
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