package scratch.kevin.prvi25.figures;

import java.io.File;

public class PRVI_Paths {
	
	public static final File PAPER_DIR =new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper");
	public static final File FIGURES_DIR =new File(PAPER_DIR, "Figures");
	
	public static final File INV_DIR = new File("/data/kevin/nshm23/batch_inversions/");
	
	public static final File CRUSTAL_DIR = new File(INV_DIR, "2025_01_17-prvi25_crustal_branches-dmSample10x");
	public static final File CRUSTAL_SOL_GRIDDED = new File(CRUSTAL_DIR, "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip");
	public static final File CRUSTAL_SOL_SUPRA_ONLY = new File(CRUSTAL_DIR, "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged.zip");
	public static final File CRUSTAL_SLT = new File(CRUSTAL_DIR, "results.zip");
	
	public static final File SUBDUCTION_DIR = new File(INV_DIR, "2025_01_17-prvi25_subduction_branches");
	public static final File SUBDUCTION_SOL_LARGE = new File(SUBDUCTION_DIR, "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip");
	public static final File SUBDUCTION_SOL_SMALL = new File(SUBDUCTION_DIR, "results_PRVI_SUB_FM_SMALL_branch_averaged_gridded.zip");
	public static final File SUBDUCTION_SOLS_COMBINED = new File(SUBDUCTION_DIR, "results_PRVI_SUB_FMs_combined_branch_averaged_gridded.zip");
	public static final File SUBDUCTION_SLT = new File(SUBDUCTION_DIR, "results.zip");
	
	public static final File COMBINED_DIR = new File(INV_DIR, "2025_01_17-prvi25_crustal_subduction_combined_branches");
	public static final File COMBINED_SOL = new File(COMBINED_DIR, "combined_branch_averaged_solution.zip");

}
