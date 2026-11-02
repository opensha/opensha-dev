package scratch.anne.sra.engine.lec_pipeline;

import java.nio.file.Path;

import scratch.anne.sra.io.writers.AssetRiskForRuptureWriter;

/**
 * All top-level options for the 4-stage true-mean LEC pipeline (per-rupture
 * EALs -> branch contributions -> pseudo loss curves -> LEC summary), plus
 * every derived path used by any stage.
 *
 * <p>Centralizing path-naming here means the four stages (and the four
 * standalone scripts they replicate -- PortfolioExpectedLossDualCalcByCSV,
 * BranchSolutionFromRuptureResults, PseudoLossCurveExportExample,
 * LogicTreeLecSummary) never duplicate or risk disagreeing on a file name;
 * every stage asks this class for its input/output paths instead of
 * rebuilding them from fileTag/fileSuffix locally.</p>
 *
 * <p>Naming convention: {@code maxBranches} truncation (stage 2 onward)
 * uses the existing {@code _subset-n-N} suffix when less than the full
 * logic tree is run. The pseudo-loss-curve CSV's separate, smaller
 * display-only branch sample (stage 3) uses a distinct
 * {@code _sample-k-M} suffix, so a reader can tell "this run only
 * processed N of the branches" (subset) apart from "this run processed
 * everything, but the CSV only shows M of them for readability" (sample)
 * -- these are not the same kind of truncation and must not share a
 * suffix.</p>
 */
public final class PipelineConfig {

    // ------------------------------------------------------------------
    // Run location
    // ------------------------------------------------------------------
    private final Path baseFolder;
    private final Path inputFolder;

    // ------------------------------------------------------------------
    // Stage 1 (per-rupture EALs) options
    // ------------------------------------------------------------------
    private final boolean keepAssetLoss;
    private final double fRelative;
    private final double fAbsolute;

    // ------------------------------------------------------------------
    // Stage 2 (branch contributions) options
    // ------------------------------------------------------------------
    private final int maxBranches;
    private final Path solutionLogicTreeZip;

    // ------------------------------------------------------------------
    // Stage 3 (pseudo loss curves) options
    // ------------------------------------------------------------------
    /** Branch indices to include in the pseudo-loss-curve CSV (display-only; Parquet always covers every branch). */
    private final java.util.Set<Integer> csvSampleBranches;

    // ------------------------------------------------------------------
    // Stage 4 (LEC summary) options
    // ------------------------------------------------------------------
    /**
     * Placeholder only -- no COV-mode accumulation path exists yet (see
     * LecSummaryStage). Reserved so switching this later naturally
     * changes stage 4's output file names (via {@link #covModeTag()}),
     * triggering only stage 4 to re-run; stages 1-3 are entirely
     * COV-agnostic and never need to be rerun when this changes.
     */
    private final String covMode;

    /** Weighted percentiles to report across the branch ensemble, per level (EAL row and every likelihood-level row). */
    private final double[] percentilesOfInterest;

    public PipelineConfig(
            Path baseFolder,
            Path inputFolder,
            boolean keepAssetLoss,
            double fRelative,
            double fAbsolute,
            int maxBranches,
            Path solutionLogicTreeZip,
            java.util.Set<Integer> csvSampleBranches,
            String covMode,
            double[] percentilesOfInterest
    ) {
        this.baseFolder = baseFolder;
        this.inputFolder = inputFolder;
        this.keepAssetLoss = keepAssetLoss;
        this.fRelative = fRelative;
        this.fAbsolute = fAbsolute;
        this.maxBranches = maxBranches;
        this.solutionLogicTreeZip = solutionLogicTreeZip;
        this.csvSampleBranches = csvSampleBranches;
        this.covMode = covMode;
        this.percentilesOfInterest = percentilesOfInterest;
    }

    // ------------------------------------------------------------------
    // Plain option getters
    // ------------------------------------------------------------------

    public Path getBaseFolder() { return baseFolder; }
    public Path getInputFolder() { return inputFolder; }
    public boolean isKeepAssetLoss() { return keepAssetLoss; }
    public double getFRelative() { return fRelative; }
    public double getFAbsolute() { return fAbsolute; }
    public int getMaxBranches() { return maxBranches; }
    public Path getSolutionLogicTreeZip() { return solutionLogicTreeZip; }
    public java.util.Set<Integer> getCsvSampleBranches() { return csvSampleBranches; }
    public String getCovMode() { return covMode; }
    public double[] getPercentilesOfInterest() { return percentilesOfInterest; }

    // ------------------------------------------------------------------
    // Derived naming
    // ------------------------------------------------------------------

    public Path getRunFolder() {
        return baseFolder.resolve(inputFolder);
    }

    /** Defaults to the input folder's name, matching every original script's own fallback. */
    public String getFileTag() {
        return inputFolder.getFileName().toString();
    }

    /**
     * {@code _subset-n-N} whenever {@code maxBranches} is finite (i.e.
     * not {@code Integer.MAX_VALUE}), else empty -- computable from
     * config alone, with NO dependency on the actual logic-tree size.
     *
     * <p>A finite {@code maxBranches} is a deliberate, explicit request
     * to truncate the run; the resulting file name always reflects that
     * request as given, regardless of how large the real logic tree
     * turns out to be. In particular, this does NOT silently fall back
     * to "no truncation after all" if {@code maxBranches} happens to be
     * {@code >=} the tree's actual size -- see
     * {@link BranchContributionsStage}, which treats that situation as
     * a user-configuration error (an explicit truncation request that
     * didn't actually truncate anything) rather than a case to reduce
     * away quietly.</p>
     */
    public String subsetSuffixFor() {
        return maxBranches == Integer.MAX_VALUE ? "" : "_subset-n-" + maxBranches;
    }

    /** {@code _sample-k-M} for the pseudo-loss-curve CSV's own smaller display subset. */
    public String sampleSuffix() {
        return "_sample-k-" + csvSampleBranches.size();
    }

    /** {@code covMode} folded into stage 4's output names, so changing it re-triggers only stage 4. */
    public String covModeTag() {
        return "_cov-" + covMode;
    }

    // ------------------------------------------------------------------
    // Stage 1 output paths
    // ------------------------------------------------------------------

    public Path getEalAllAssetsCsv() {
        return getRunFolder().resolve(getFileTag() + "_eal-all-assets.csv");
    }

    public Path getAggregatedEalCsv() {
        return getRunFolder().resolve(getFileTag() + "_aggregated-eal.csv");
    }

    public Path getPerRuptureEalCsv() {
        return getRunFolder().resolve(getFileTag() + "_eal-per-rup.csv");
    }

    /**
     * The per-asset-per-rupture subdirectory NAME (not a full resolved
     * path) -- {@code AssetRiskForRuptureWriter} takes an output
     * directory and this tag separately, resolving
     * {@code outputDir.resolve(tag)} itself (see its constructor).
     * Delegates to {@link AssetRiskForRuptureWriter#taggedWithThresholds}
     * directly so two runs with different {@code fRelative}/
     * {@code fAbsolute} thresholds never collide on the same folder,
     * and the threshold used is visible from the folder name itself.
     */
    public String getAssetEalPerRupTag() {
        return AssetRiskForRuptureWriter.taggedWithThresholds(
                getFileTag() + "_asset-eal-per-rup", fRelative, fAbsolute);
    }

    // ------------------------------------------------------------------
    // Stage 2 output paths
    // ------------------------------------------------------------------

    public Path getBranchEalsCsv() {
        return getRunFolder().resolve(
                getFileTag() + "_branch-eals" + subsetSuffixFor() + ".csv");
    }

    public Path getBranchContributionsParquet() {
        return getRunFolder().resolve(
                getFileTag() + "_branch-rupture-contributions" + subsetSuffixFor() + ".parquet");
    }

    // ------------------------------------------------------------------
    // Stage 3 output paths
    // ------------------------------------------------------------------

    public Path getPseudoLossCurveCsv() {
        return getRunFolder().resolve(
                getFileTag() + "_branch-pseudo-loss-curve" + subsetSuffixFor()
                        + sampleSuffix() + ".csv");
    }

    public Path getPseudoLossCurveParquet() {
        return getRunFolder().resolve(
                getFileTag() + "_branch-pseudo-loss-curve" + subsetSuffixFor() + ".parquet");
    }

    // ------------------------------------------------------------------
    // Stage 4 output paths
    // ------------------------------------------------------------------

    public Path getBranchEalsWithLecCsv() {
        return getRunFolder().resolve(
                getFileTag() + "_branch-eals-with-lec" + subsetSuffixFor()
                        + covModeTag() + ".csv");
    }

    public Path getBranchPercentilesCsv() {
        return getRunFolder().resolve(
                getFileTag() + "_branch-percentiles" + subsetSuffixFor()
                        + covModeTag() + ".csv");
    }
}
