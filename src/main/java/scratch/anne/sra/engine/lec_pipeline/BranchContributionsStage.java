package scratch.anne.sra.engine.lec_pipeline;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
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

import scratch.anne.sra.domain.hazard.HazardParametersMetadataBuilder;
import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.engine.logic_tree.BranchEALCalculator;
import scratch.anne.sra.io.readers.RuptureLossLookupReader;
import scratch.anne.sra.io.readers.TrueMeanSetupFromRuptureResults;
import scratch.anne.sra.io.readers.TrueMeanSetupFromRuptureResults.TrueMeanSetup;
import scratch.anne.sra.io.writers.LogicTreeBranchEalSummaryExporter;
import scratch.anne.sra.io.writers.LogicTreeBranchRuptureContributionExporter;

/**
 * Pipeline stage 2: per-branch EALs and per-branch rupture contributions,
 * rebuilt from stage 1's exported per-rupture EAL CSV.
 *
 * <p>Replicates {@code BranchSolutionFromRuptureResults}'s body. No
 * existing class is modified; this class only calls the same public
 * readers/calculators/writers the original script already used.</p>
 *
 * <p>Gated by existence: the output file names depend only on
 * {@code config.getMaxBranches()} (see
 * {@link PipelineConfig#subsetSuffixFor()}), NOT on the actual logic-tree
 * size, so existence is checked FIRST, directly from config, with no
 * SolutionLogicTree (or true-mean ERF) load at all. The SLT/ERF are only
 * loaded if the outputs are actually missing and this stage has real
 * work to do.</p>
 *
 * <p>When this stage does run with a finite {@code maxBranches}, it
 * requires {@code maxBranches < logicTree.size()}: a finite value is
 * read as a deliberate, explicit truncation request, and if the real
 * tree isn't actually larger than that, the request doesn't describe
 * reality -- this is treated as a configuration error (fail loudly) and
 * NOT silently reduced to "run everything" the way {@code Math.min}
 * would. If you want to run the full tree, set {@code maxBranches} to
 * {@code Integer.MAX_VALUE} explicitly rather than relying on a large
 * finite number happening to exceed the tree's actual size.</p>
 */
public final class BranchContributionsStage {

    private BranchContributionsStage() {}

    public static void run(PipelineConfig config) throws IOException {

        Path branchEALsCsv = config.getBranchEalsCsv();
        Path branchContributionsParquet = config.getBranchContributionsParquet();

        if (Files.exists(branchEALsCsv) && Files.exists(branchContributionsParquet)) {
            System.out.println("[Stage 2] Skipping -- already exist: " + branchEALsCsv + ", " + branchContributionsParquet);
            return;
        }

        System.out.println("[Stage 2] Not found -- running branch contribution calculation...");

        Path exportCsv = config.getPerRuptureEalCsv();

        System.out.println("[Stage 2] Rebuilding true mean ERF from exported HazardParameters: " + exportCsv);
        TrueMeanSetup setup = TrueMeanSetupFromRuptureResults.rebuild(exportCsv);

        FaultSystemSolution trueMeanSol = setup.trueMeanSol();
        BaseFaultSystemSolutionERF trueMeanERF = setup.trueMeanERF();
        GridSourceList meanGridList = setup.meanGridList();
        Table<Integer, TectonicRegionType, Integer> meanLocToSourceTable = setup.meanLocToSourceTable();
        float gridMinMag = setup.gridMinMag();

        System.out.println("[Stage 2] True mean ERF ready: "
                + trueMeanERF.getNumSources() + " sources, "
                + trueMeanERF.getNumFaultSystemSources() + " fault-based");

        // downsampled model (fault and gridded) -- this part of the demo
        // still needs its own source separately, since it isn't part of
        // what gets built/exported per-branch.
        SolutionLogicTree slt = SolutionLogicTree.load(config.getSolutionLogicTreeZip().toFile());

        LogicTree<?> logicTree = slt.getLogicTree();

        int maxBranches = config.getMaxBranches();
        if (maxBranches != Integer.MAX_VALUE) {
            Preconditions.checkState(maxBranches < logicTree.size(),
                    "maxBranches=%s is >= the logic tree's actual size (%s) -- a finite maxBranches is a "
                    + "deliberate truncation request and must genuinely truncate. Set maxBranches to "
                    + "Integer.MAX_VALUE if you intend to run the full tree.",
                    maxBranches, logicTree.size());
        }
        int numBranchesToRun = Math.min(maxBranches, logicTree.size());

        System.out.println("[Stage 2] Running " + numBranchesToRun + "/" + logicTree.size() + " branches...");

        // ------------------------------------------------------------
        // Fixed per-rupture loss lookup: read ONCE from stage 1's export
        // (ruptureLoss column only -- see RuptureLossLookupReader). This
        // never changes across branches; only each branch's RATES change.
        // ------------------------------------------------------------
        Map<RuptureKey, Double> ruptureLossLookup = RuptureLossLookupReader.read(exportCsv);
        System.out.println("[Stage 2] Loaded fixed loss lookup for " + ruptureLossLookup.size() + " ruptures");

        // this stores the data for finding branch-specific ruptures in the true mean solution
        TrueMeanRuptureMappings branchMappings = trueMeanSol.getRupSet().requireModule(TrueMeanRuptureMappings.class);

        // this will suppress prints on each branch
        ModuleContainer.VERBOSE_DEFAULT = false;

        if (numBranchesToRun < logicTree.size())
            System.out.println("TEST MODE: only running " + numBranchesToRun + "/" + logicTree.size() + " branches");

        // per-branch EAL results, indexed by branch index b.
        double[] branchEALs = new double[numBranchesToRun];
        double[] branchFssRates = new double[numBranchesToRun];
        double[] branchGriddedRates = new double[numBranchesToRun];

        // this will be used to preload the next solution in parallel because most of the time in this loop is spent
        // just waiting on I/O to load the solutions for each branch
        CompletableFuture<FaultSystemSolution> preloadFuture = null;
        Stopwatch watch = Stopwatch.createStarted();

        // numBranchesToRun/logicTree.size() is the pipeline-level maxBranches
        // truncation (already reflected in the filename's _subset-n-N
        // suffix, if any -- see PipelineConfig#subsetSuffixFor), but the
        // solution logic tree zip itself (here, potentially already a
        // DOWNSAMPLED subset of a much larger original tree, e.g. 1000
        // branches downsampled from 2 million) has no record anywhere
        // else -- record both, PLUS the full HazardParameters that
        // produced the true-mean ERF these branch ruptures are mapped
        // against (same field set RuptureResultsExporter already writes
        // into _eal-per-rup.csv -- see HazardParametersMetadataBuilder),
        // PLUS the portfolio/response-model provenance RuptureResultsExporter
        // also wrote into that same _eal-per-rup.csv (see
        // PortfolioMetadataBuilder/ResponseModelMetadataBuilder), recovered
        // here via setup.portfolioAndResponseModelMetadata() with no new
        // file I/O of its own (TrueMeanSetupFromRuptureResults.rebuild
        // already parsed the full metadata block to build setup),
        // directly in this Parquet's own file-level metadata, since this
        // is the one file every downstream consumer (StreamingBranchLecBuilder,
        // PseudoLossCurveStage, LecSummaryStage) reads -- see
        // LogicTreeBranchRuptureContributionExporter's javadoc.
        Map<String, String> contribMetadata = new java.util.HashMap<>(
                HazardParametersMetadataBuilder.build(setup.hazardParameters()));
        contribMetadata.putAll(setup.portfolioAndResponseModelMetadata());
        contribMetadata.put("numBranches", numBranchesToRun + "/" + logicTree.size());
        contribMetadata.put("solutionLogicTreeZip", config.getSolutionLogicTreeZip().getFileName().toString());
        contribMetadata.put("fileTag", config.getFileTag());

        try (LogicTreeBranchRuptureContributionExporter contribExporter =
                new LogicTreeBranchRuptureContributionExporter(branchContributionsParquet, contribMetadata)) {
            for (int b = 0; b < numBranchesToRun; b++) {
                LogicTreeBranch<?> branch = logicTree.getBranch(b);
                double branchWeight = logicTree.getBranchWeight(branch);

                System.out.println("Handling branch " + b + "/" + numBranchesToRun + ": " + branch + " (wt=" + (float) branchWeight + ")");
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
                int nextIndex = b + 1;
                if (nextIndex < numBranchesToRun) {
                    preloadFuture = CompletableFuture.supplyAsync(() -> {
                        FaultSystemSolution futureSol;
                        try {
                            futureSol = slt.forBranch(logicTree.getBranch(nextIndex));
                        } catch (IOException e) {
                            throw ExceptionUtils.asRuntimeException(e);
                        }
                        futureSol.getGridSourceProvider();
                        return futureSol;
                    });
                }

                BranchEALCalculator branchCalc = new BranchEALCalculator(ruptureLossLookup);

                System.out.println("\tValidating FSS rupture mappings...");

                int[] branchRupIndexes = branchMappings.getRuptureMappings(branch);

                int numMappedFSS = 0;
                double mappedFSSRate = 0d;

                for (int i = 0; i < branchRupIndexes.length; i++) {
                    double rupRate = branchSol.getRateForRup(i);

                    int trueMeanIndex = branchRupIndexes[i];
                    if (trueMeanIndex < 0) {
                        Preconditions.checkState(rupRate == 0d);
                        continue;
                    } else if (rupRate == 0d) {
                        continue;
                    }

                    int trueMeanSourceIndex = trueMeanERF.getSrcIndexForFltSysRup(trueMeanIndex);
                    ProbEqkSource trueMeanSource = trueMeanERF.getSource(trueMeanSourceIndex);
                    int trueMeanRupIndexInSource = -1;
                    double rupMag = branchRupSet.getMagForRup(i);
                    for (int rupIndex = 0; rupIndex < trueMeanSource.getNumRuptures(); rupIndex++) {
                        ProbEqkRupture rup = trueMeanSource.getRupture(rupIndex);
                        if (Precision.equals(rup.getMag(), rupMag)) {
                            Preconditions.checkState(trueMeanRupIndexInSource < 0, "Duplicate match found?");
                            trueMeanRupIndexInSource = rupIndex;
                        }
                    }

                    Preconditions.checkState(trueMeanRupIndexInSource >= 0,
                            "No matching rupture found for branch %s FSS rupture %s with rate=%s",
                            b, i, rupRate);

                    branchCalc.accumulateFault(trueMeanSourceIndex, trueMeanRupIndexInSource, rupRate);

                    Double rupLoss = ruptureLossLookup.get(new RuptureKey(trueMeanSourceIndex, trueMeanRupIndexInSource));
                    if (rupLoss != null)
                        contribExporter.write(b, trueMeanSourceIndex, trueMeanRupIndexInSource, rupLoss, rupRate);

                    numMappedFSS++;
                    mappedFSSRate += rupRate;
                }
                System.out.println("\t\tMapped " + numMappedFSS + "/" + branchRupIndexes.length + " FSS ruptures with total rate " + (float) mappedFSSRate);
                branchFssRates[b] = mappedFSSRate;

                System.out.println("\tValidating gridded rupture mappings...");
                Preconditions.checkState(branchGridList.getNumLocations() == meanGridList.getNumLocations(),
                        "Branch %s grid list has %s locations, but mean has %s", b, branchGridList.getNumLocations(), meanGridList.getNumLocations());
                Preconditions.checkState(branchGridList.getNumSources() <= meanGridList.getNumSources(),
                        "Branch %s grid list has %s sources, but mean has %s", b, branchGridList.getNumSources(), meanGridList.getNumSources());
                int numMappedGrid = 0;
                double mappedGridRate = 0d;
                int[] prevMappings = null;
                int numSourceMappingsReused = 0;
                for (int s = 0; s < branchGridList.getNumSources(); s++) {
                    int gridIndex = branchGridList.getLocationIndexForSource(s);
                    TectonicRegionType trt = branchGridList.tectonicRegionTypeForSourceIndex(s);
                    Preconditions.checkState(meanLocToSourceTable.contains(gridIndex, trt));
                    int trueMeanSourceIndex = trueMeanERF.getNumFaultSystemSources() + meanLocToSourceTable.get(gridIndex, trt);
                    ProbEqkSource source = trueMeanERF.getSource(trueMeanSourceIndex);
                    Preconditions.checkState(source instanceof PointSource,
                            "True mean source %s isn't a point source? %s", trueMeanSourceIndex, source.getClass());
                    PointSource ptSource = (PointSource) source;
                    int numSourceRups = ptSource.getNumRuptures();
                    List<ProbEqkRupture> sourceRups = ptSource.getRuptureList();

                    Location loc = branchGridList.getLocation(gridIndex);
                    Preconditions.checkState(LocationUtils.areSimilar(loc, ptSource.getLocation()));
                    List<GriddedRupture> gridRups = branchGridList.getRuptures(branchGridList.tectonicRegionTypeForSourceIndex(s), gridIndex);

                    int[] myMappings = null;
                    if (prevMappings == null || gridRups.size() != prevMappings.length) {
                        prevMappings = null;
                        myMappings = new int[gridRups.size()];
                    }
                    for (int g = 0; g < gridRups.size(); g++) {
                        GriddedRupture gridRup = gridRups.get(g);
                        if (gridRup.properties.magnitude < gridMinMag)
                            continue;
                        int trueMeanRupIndexInSource;
                        if (prevMappings != null && prevMappings[g] < numSourceRups && isGridMatch(gridRup, sourceRups.get(prevMappings[g]))) {
                            trueMeanRupIndexInSource = prevMappings[g];
                        } else {
                            prevMappings = null;
                            trueMeanRupIndexInSource = -1;
                            for (int rupIndex = 0; rupIndex < numSourceRups; rupIndex++) {
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

                        double rupRate = gridRup.rate;

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
                System.out.println("\t\tMapped " + numMappedGrid + " gridded ruptures across " + branchGridList.getNumSources()
                        + " sources with total rate " + (float) mappedGridRate + " (" + numSourceMappingsReused + " source mappings reused)");
                branchGriddedRates[b] = mappedGridRate;

                branchEALs[b] = branchCalc.getEAL();
                System.out.println("\tBranch " + b + " EAL: " + branchEALs[b]);

                System.out.println("\tDone with branch " + b + "/" + numBranchesToRun + "; " + timeStr(watch) + " elapsed, estimated "
                        + timeLeftStr(watch, b + 1, numBranchesToRun) + " remaining");
            }
        }
        System.out.println("[Stage 2] DONE");
        watch.stop();
        System.out.println("Took " + timeStr(watch));

        List<LogicTreeBranchEalSummaryExporter.BranchEalRow> rows = new ArrayList<>(numBranchesToRun);
        for (int b = 0; b < numBranchesToRun; b++) {
            rows.add(new LogicTreeBranchEalSummaryExporter.BranchEalRow(
                    b,
                    logicTree.getBranchWeight(logicTree.getBranch(b)),
                    branchEALs[b],
                    branchFssRates[b],
                    branchGriddedRates[b],
                    logicTree.getBranch(b)));
        }
        // numBranchesToRun/logicTree.size() is the pipeline-level maxBranches
        // truncation (already reflected in the filename's _subset-n-N
        // suffix, if any -- see PipelineConfig#subsetSuffixFor). Same
        // full provenance (numBranches, solutionLogicTreeZip, fileTag,
        // HazardParameters) already written into
        // branchContributionsParquet's own file-level metadata above --
        // reused here rather than rebuilt, so both outputs agree exactly
        String metadataComment = formatMetadataComment(contribMetadata);
        LogicTreeBranchEalSummaryExporter.write(rows, logicTree, metadataComment, branchEALsCsv);
        System.out.println("[Stage 2] Wrote branch EALs to: " + branchEALsCsv);
        System.out.println("[Stage 2] Wrote branch contributions to: " + branchContributionsParquet);
    }

    /**
     * Flattens a key-value map into one {@code key=val,key2=val2,...}
     * line, for the CSV comment line -- same shape as
     * {@code PseudoLossCurveStage}/{@code LecSummaryStage}'s own
     * identical helper.
     */
    private static String formatMetadataComment(Map<String, String> metadata) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, String> e : metadata.entrySet()) {
            if (sb.length() > 0)
                sb.append(',');
            sb.append(e.getKey()).append('=').append(e.getValue());
        }
        return sb.toString();
    }

    private static final java.text.DecimalFormat TWO_DF = new java.text.DecimalFormat("0.00");

    private static String timeStr(Stopwatch watch) {
        return timeStr(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d);
    }

    private static String timeLeftStr(Stopwatch watch, int done, int tot) {
        double totSecs = watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
        double secsPer = totSecs / (double) done;
        double secsLeft = secsPer * (tot - done);
        return timeStr(secsLeft);
    }

    private static String timeStr(double secs) {
        double mins = secs / 60d;
        double hours = mins / 60d;
        if (hours > 1.5)
            return TWO_DF.format(hours) + " h";
        else if (mins > 1.5)
            return TWO_DF.format(mins) + " m";
        return TWO_DF.format(secs) + " s";
    }

    private static boolean isGridMatch(GriddedRupture gridRup, ProbEqkRupture erfRup) {
        return Precision.equals(erfRup.getMag(), gridRup.properties.magnitude)
                && Precision.equals(erfRup.getAveRake(), gridRup.properties.rake);
    }
}
