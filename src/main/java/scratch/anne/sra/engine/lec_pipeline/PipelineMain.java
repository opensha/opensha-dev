package scratch.anne.sra.engine.lec_pipeline;

import java.nio.file.Path;
import java.util.Set;

/**
 * Single entry point for the 4-stage true-mean LEC pipeline:
 *
 * <pre>
 *   1. PerRuptureEalStage        -- per-rupture EALs for the true-mean solution
 *   2. BranchContributionsStage  -- per-branch EALs + per-branch rupture contributions
 *   3. PseudoLossCurveStage      -- full discrete-rank pseudo loss curves
 *   4. LecSummaryStage           -- per-branch LEC summary + ensemble percentile summary
 * </pre>
 *
 * <p>Each stage checks whether its own output(s) already exist before
 * doing any work, so re-running this script after a partial or complete
 * prior run only recomputes what is actually missing. No existing class
 * is modified anywhere in this pipeline -- every stage is a thin
 * orchestration layer over the same public readers/writers/calculators
 * the original four standalone scripts
 * (PortfolioExpectedLossDualCalcByCSV, BranchSolutionFromRuptureResults,
 * PseudoLossCurveExportExample, LogicTreeLecSummary) already used.</p>
 */
public final class PipelineMain {

    // ====================================================================
    // OPTIONS -- everything you're likely to want to change lives here.
    // ====================================================================

    // --- Run location ---
    private static final Path BASE_FOLDER = Path.of(
            "C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\BERM\\results\\EAL\\true_mean_solution");
//    private static final Path INPUT_FOLDER = Path.of("RES1-W1-MC_Porter-vuln_vs-360");
    private static final Path INPUT_FOLDER = Path.of("short-portfolio");

    // --- Stage 1 (per-rupture EALs) ---
    /** Whether to also write per-asset-per-rupture loss detail (stage 1 only). */
    private static final boolean KEEP_ASSET_LOSS = true;
    private static final double F_RELATIVE = 0.01d;
    private static final double F_ABSOLUTE = 0d;

    // --- Stage 2 (branch contributions) ---
    /** Cap the number of branches processed, for quick tests. Set to Integer.MAX_VALUE to run the full tree. */
    private static final int MAX_BRANCHES = 10;
//     private static final int MAX_BRANCHES = Integer.MAX_VALUE;
    private static final Path SOLUTION_LOGIC_TREE_ZIP = Path.of(
            "C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\"
            + "my_scratch\\logic_tree_reweighting\\2024_02_02-nshm23_branches-WUS_FM_v3-gridded_rebuild",
            "results_full_gridded_downsampled_1k.zip");

    // --- Stage 3 (pseudo loss curves) ---
    /** Branches shown in the human-readable pseudo-loss-curve CSV (Parquet always covers every branch regardless). */
    private static final Set<Integer> CSV_SAMPLE_BRANCHES = Set.of(0, 1, 2, 3, 4, 5, 6, 7, 8, 9);

    // --- Stage 4 (LEC summary) ---
    /**
     * Placeholder only -- no COV-mode accumulation path exists yet (see
     * LecSummaryStage's javadoc). Changing this only changes stage 4's
     * output file names for now; it has no other effect until COV-mode
     * support is actually implemented.
     */
    private static final String COV_MODE = "discrete-rank";
    /** Weighted percentiles to report across the branch ensemble, per level (EAL row and every likelihood-level row). */
    private static final double[] PERCENTILES_OF_INTEREST = { 1, 2.5, 5, 10, 16, 50, 84, 90, 95, 97.5, 99 };

    // ====================================================================

    public static void main(String[] args) throws Exception {

        PipelineConfig config = new PipelineConfig(
                BASE_FOLDER,
                INPUT_FOLDER,
                KEEP_ASSET_LOSS,
                F_RELATIVE,
                F_ABSOLUTE,
                MAX_BRANCHES,
                SOLUTION_LOGIC_TREE_ZIP,
                CSV_SAMPLE_BRANCHES,
                COV_MODE,
                PERCENTILES_OF_INTEREST
        );

        PerRuptureEalStage.run(config);
        BranchContributionsStage.run(config);
        PseudoLossCurveStage.run(config);
        LecSummaryStage.run(config);

        System.out.println("Pipeline complete.");
    }
}
