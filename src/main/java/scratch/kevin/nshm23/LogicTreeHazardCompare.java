package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.ShawSegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSectConstraintModels;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;

public class LogicTreeHazardCompare {
	
	private static boolean TITLES = false;
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		boolean currentWeights = false;
		boolean compCurrentWeights = false;
		
//		File mainDir = new File(invDir, "2021_11_24-nshm23_draft_branches-FM3_1");
//		String mainName = "NSHM23 Draft";
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = new File(invDir, "2021_11_23-u3_branches-FM3_1-5h");
//		String compName = "UCERF3 Redo";
//		LogicTreeNode[] compSubsetNodes = null;
//		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3_redo");
////		File compDir = new File(invDir, "2021_11_30-u3_branches-orig_calcs-5h");
////		String compName = "UCERF3 As Published";
////		LogicTreeNode[] compSubsetNodes = null;
////		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3_as_published");

////		File mainDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip");
////		String mainName = "Draft-Model-Coulomb-Shaw-Seg-Jump-Prob";
//		File mainDir = new File(invDir, "2022_02_22-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-StrictEquivJumpProb-2000ip");
//		String mainName = "Seg-Strict-Equiv-Jump-Prob";
////		File mainDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-CappedRdst-2000ip");
////		String mainName = "Draft-Model-Coulomb-Shaw-Seg-Capped-Redist";
////		String mainName = "Max-Dist-Seg";
////		File mainDir = new File(invDir, "2022_02_08-nshm23_u3_hybrid_branches-seg_bin_dist_capped_distr-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
////		String mainName = "Draft-Model-Coulomb-Shaw-Seg";
//		LogicTreeNode[] subsetNodes = null;
////		LogicTreeNode[] subsetNodes = { SubSectConstraintModels.TOT_NUCL_RATE };
////		File compDir = null;
////		String compName = null;
////		File outputDir = new File(mainDir, "hazard_maps");
////		File compDir = new File(invDir, "2022_02_10-nshm23_u3_hybrid_branches-seg-no_adj-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
////		String compName = "No-Seg";
////		File outputDir = new File(mainDir, "hazard_maps_comp_no_seg");
//////		File compDir = new File(invDir, "2022_02_10-nshm23_u3_hybrid_branches-seg-capped_self_contained-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
////		File compDir = new File(invDir, "2022_02_10-nshm23_u3_hybrid_branches-seg-greedy_self_contained-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
////		String compName = "Self-Contained";
////		File outputDir = new File(mainDir, "hazard_maps_comp_self_contained");
//		File compDir = new File(invDir, "2022_02_11-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		String compName = "Strict-Cutoff-Seg";
//		File outputDir = new File(mainDir, "hazard_maps_comp_strict_cutoff");
////		File compDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-CappedRdst-2000ip");
////		String compName = "Seg-Capped-Redist";
////		File outputDir = new File(mainDir, "hazard_maps_comp_capped_redist");
////		File compDir = new File(invDir, "2022_02_08-nshm23_u3_hybrid_branches-seg_bin_dist_capped_distr-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
////		String compName = "Capped-Distribution-Seg";
////		File outputDir = new File(mainDir, "hazard_maps_comp_capped_distr_tot_rate");
////		File compDir = new File(invDir, "2022_01_28-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-SubB1-5000ip");
////		String compName = "5000 Iters/Rup";
////		File outputDir = new File(mainDir, "hazard_maps_comp_5000ip");
//		LogicTreeNode[] compSubsetNodes = null;
////		LogicTreeNode[] compSubsetNodes = { SubSectConstraintModels.TOT_NUCL_RATE };
		
		/*
		 * Sweep of draft models, late Feb 2022
		 */
//		File mainDir = new File(invDir, "2022_02_16-u3_branches-orig_calc_params-FM3_1");
//		String mainName = "UCERF3 2022 Reproduction";
		
//		File mainDir = new File(invDir, "2022_02_17-u3_branches-orig_calc_params-new_avg-converged-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Avg Converged";
		
//		File mainDir = new File(invDir, "2022_02_20-u3_branches-orig_calc_params-new_avg-converged-noWL-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Avg Converged No-WL";
		
//		File mainDir = new File(invDir, "2022_02_20-u3_branches-orig_calc_params-new_avg-converged-new_perturb-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Avg Converged New-Perturb";
		
//		File mainDir = new File(invDir, "2022_02_22-u3_branches-orig_calc_params-new_avg-converged-noWL-new_perturb-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Avg Converged No-WL New-Perturb";
		
//		File mainDir = new File(invDir, "2022_02_21-u3_branches-orig_calc_params-new_avg-converged-noWL-new_perturb-try_zero_often-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Avg Converged New-Perturb No-WL Try-Zero";
		
//		File mainDir = new File(invDir, "2022_02_22-u3_branches-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Anneal Params";
		
//		File mainDir = new File(invDir, "2022_02_22-u3_branches-thinned_0.1-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Anneal Params, Reduced RS";
		
//		File mainDir = new File(invDir, "2022_03_25-u3_branches-coulomb-FM3_1-2000ip");
//		String mainName = "UCERF3 New-Anneal Params, Coulomb RS";
		
//		File mainDir = new File(invDir, "2022_02_23-nshm23_u3_hybrid_branches-no_seg-FM3_1-U3RupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		String mainName = "NSHM23 Draft, U3 RS, No Seg";
		
//		File mainDir = new File(invDir, "2022_02_23-nshm23_u3_hybrid_branches-no_seg-FM3_1-U3RedRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		String mainName = "NSHM23 Draft, U3 Reduced RS, No Seg";
		
//		File mainDir = new File(invDir, "2022_02_23-nshm23_u3_hybrid_branches-no_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, No Seg";
		
//		File mainDir = new File(invDir, "2022_02_27-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, Strict Seg";
//		currentWeights = true;
		
//		File mainDir = new File(invDir, "2022_02_27-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip-branch-translated");
//		String mainName = "NSHM23 Draft, Coulomb RS, Strict Seg Translated";
		
//		File mainDir = new File(invDir, "2022_02_27-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip-branch-translated-min3km");
//		String mainName = "NSHM23 Draft, Coulomb RS, Strict Seg Translated >=3km";
		
//		File mainDir = new File(invDir, "2022_02_23-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-NoAdj-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, Seg No-Adj";
		
//		File mainDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, Seg Thresh-Avg";
		
//		File mainDir = new File(invDir, "2022_02_22-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-StrictEquivJumpProb-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, Seg Thresh-Avg @ Strict Bins";
		
//		File mainDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String mainName = "NSHM23 Draft, Coulomb RS, Seg Thresh-Avg 1km-Shift";
		
//		File mainDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-CappedRdst-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, Seg Capped-Redist";
		
//		File mainDir = new File(invDir, "2022_03_03-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-NuclMFD-SubB1-JumpProb-2000ip");
//		String mainName = "NSHM23 Draft, Coulomb RS, Seg Thresh-Avg 1km-Shift, Nucl-MFD";
		
//		File compDir = new File(invDir, "2022_03_24-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip");
//		String compName = "Previous Step";
//		File outputDir = new File(mainDir, "hazard_maps_comp_prev");

//		File mainDir = new File(invDir, "2022_04_27-u3_branches-scaleLowerDepth1.3-FM3_1-2000ip");
//		String mainName = "UCERF3, 1.3x Lower-Seis Depth";
		
//		File mainDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-cluster_specific_inversion-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String mainName = "NSHM23 Draft, Cluster-Specific Inversions";
		
//		File mainDir = new File(invDir, "2022_05_20-nshm23_u3_hybrid_branches-scaleLowerDepth1.3-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String mainName = "NSHM23 Draft, 1.3x Lower-Seis Depth";
		
//		File mainDir = new File(invDir, "2022_05_20-nshm23_u3_hybrid_branches-test_scale_rels-scaleLowerDepth1.3-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String mainName = "NSHM23 Draft, Test Scale Rels, 1.3x Lower-Seis Depth";
		
//		File mainDir = new File(invDir, "2022_05_20-nshm23_u3_hybrid_branches-test_scale_rels-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String mainName = "NSHM23 Draft, Test Scale Rels";
		
//		File mainDir = new File(invDir, "2022_07_23-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep");
		File mainDir = new File(invDir, "2022_07_25-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR-IncludeThruCreep");
		String mainName = "NSHM23 Draft";
		
		LogicTreeNode[] subsetNodes = null;
		LogicTreeNode[] compSubsetNodes = null;
//		LogicTreeNode[] compSubsetNodes = { FaultModels.FM3_1 };
//		File compDir = null;
//		String compName = null;
//		File outputDir = new File(mainDir, "hazard_maps");
//		File compDir = new File(invDir, "2022_07_23-nshm23_branches-no_seg-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-IncludeThruCreep");
//		String compName = "No Segmentation";
//		File outputDir = new File(mainDir, "hazard_maps_comp_no_seg");
		File compDir = new File(invDir, "2022_07_23-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep");
		String compName = "Prev Segmentation";
		File outputDir = new File(mainDir, "hazard_maps_comp_prev_seg");
//		File compDir = new File(invDir, "2021_11_30-u3_branches-orig_calcs-5h");
//		String compName = "UCERF3 As Published";
//		File outputDir = new File(mainDir, "hazard_maps_comp_ucerf3_as_published");
//		String compName = "UCERF3 As Published, Uniform";
//		File outputDir = new File(mainDir, "hazard_maps_comp_ucerf3_as_published_uniform");
//		LogicTreeNode[] compSubsetNodes = { FaultModels.FM3_1, SlipAlongRuptureModels.UNIFORM };
//		File compDir = new File(invDir, "2022_03_24-u3_branches-FM3_1-2000ip");
//		String compName = "UCERF3 New Anneal";
//		File outputDir = new File(mainDir, "hazard_maps_comp_ucerf3_new_anneal");
//		File compDir = new File(invDir, "2022_03_13-u3_branches-coulomb-bilateral-FM3_1-2000ip");
//		String compName = "UCERF3 Coulomb Bilateral";
//		File outputDir = new File(mainDir, "hazard_maps_comp_ucerf3_coulomb_bilateral");
//		String compName = "UCERF3 New Anneal, Uniform";
//		File outputDir = new File(mainDir, "hazard_maps_comp_ucerf3_new_anneal_uniform");
//		LogicTreeNode[] compSubsetNodes = { FaultModels.FM3_1, SlipAlongRuptureModels.UNIFORM };
//		File compDir = new File(invDir, "2022_02_27-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		String compName = "Strict Segmentation";
//		File outputDir = new File(mainDir, "hazard_maps_comp_strict_seg");
//		File compDir = new File(invDir, "");
//		String compName = "Threshold Avg Shift-1km";
//		String compName = "NSHM23 Draft";
//		File outputDir = new File(mainDir, "hazard_maps_comp_thresh_avg_shift_1km");
//		File compDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String compName = "NSHM23 Draft, Test Scale Rels";
//		File outputDir = new File(mainDir, "hazard_maps_comp_test_scale_rels");
//		File compDir = new File(invDir, "2022_05_20-nshm23_u3_hybrid_branches-test_scale_rels-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		String compName = "UCERF3 Full Segmented";
//		File outputDir = new File(mainDir, "hazard_maps_comp_fully_segmented");
		
		/**
		 * NSHM18 ingredients
		 */
//		File mainDir = new File(invDir, "2022_04_08-nshm18_branches-shift_seg_1km-NSHM18_WUS_NoCA-CoulombRupSet-DsrUni-TotNuclRate-NoRed-JumpProb-10000ip");
//		String mainName = "NSHM18 Ingred, Coulomb RS, Seg Thresh-Avg 1km-Shift";
////		File mainDir = new File(invDir, "2022_04_08-nshm18_branches-no_seg-NSHM18_WUS_NoCA-AzimuthalRupSet-DsrUni-TotNuclRate-NoRed-10000ip");
////		String mainName = "NSHM18 Ingred, Azimuthal RS, No Seg";
////		File mainDir = new File(invDir, "2022_04_08-nshm18_branches-no_seg-NSHM18_WUS_NoCA-FullSegRupSet-DsrUni-TotNuclRate-NoRed-10000ip");
////		String mainName = "NSHM18 Ingred, Full Segmented RS";
//		
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = new File(invDir, "2022_04_08-nshm18_branches-no_seg-NSHM18_WUS_NoCA-FullSegRupSet-DsrUni-TotNuclRate-NoRed-10000ip");
//		String compName = "Fully Segmented RS";
//		File outputDir = new File(mainDir, "hazard_maps_comp_full_segmented");
////		File compDir = new File(invDir, "2022_04_08-nshm18_branches-shift_seg_1km-NSHM18_WUS_NoCA-CoulombRupSet-DsrUni-TotNuclRate-NoRed-JumpProb-10000ip");
////		String compName = "Coulomb RS w/ Seg";
////		File outputDir = new File(mainDir, "hazard_maps_comp_coulomb");
//		LogicTreeNode[] compSubsetNodes = null;
		
		SolutionLogicTree solTree = SolutionLogicTree.load(new File(mainDir, "results.zip"));
		
		ReturnPeriods[] rps = ReturnPeriods.values();
		double[] periods = { 0d, 1d };
		double spacing = 0.1;
//		double spacing = 0.25;
		
		LogicTree<?> tree = solTree.getLogicTree();
		if (subsetNodes != null)
			tree = tree.matchingAll(subsetNodes);
		if (currentWeights)
			tree.setWeightProvider(new BranchWeightProvider.CurrentWeights());
		LogicTreeHazardCompare mapper = new LogicTreeHazardCompare(solTree, tree,
				new File(mainDir, "results_hazard.zip"), rps, periods, spacing);
		
		LogicTreeHazardCompare comp = null;
		if (compDir != null) {
			SolutionLogicTree compSolTree = SolutionLogicTree.load(new File(compDir, "results.zip"));
			LogicTree<?> compTree = compSolTree.getLogicTree();
			if (compSubsetNodes != null)
				compTree = compTree.matchingAll(compSubsetNodes);
			if (compCurrentWeights)
				compTree.setWeightProvider(new BranchWeightProvider.CurrentWeights());
			comp = new LogicTreeHazardCompare(compSolTree, compTree,
					new File(compDir, "results_hazard.zip"), rps, periods, spacing);
		}
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		mapper.buildReport(outputDir, mainName, comp, compName);
	}
	private ReturnPeriods[] rps;
	private double[] periods;
	
	private List<? extends LogicTreeBranch<?>> branches;
	private List<Double> weights;
	private double totWeight;
	private Table<ReturnPeriods, Double, GriddedGeoDataSet[]> maps;
	private CPT logCPT;
	private CPT spreadCPT;
	private CPT spreadDiffCPT;
	private CPT pDiffCPT;
	private CPT percentileCPT;
	
	private SolHazardMapCalc mapper;
	private SolutionLogicTree solLogicTree;

	public LogicTreeHazardCompare(SolutionLogicTree solLogicTree, File mapsZipFile,
			ReturnPeriods[] rps, double[] periods, double spacing) throws IOException {
		this(solLogicTree, solLogicTree.getLogicTree(), mapsZipFile, rps, periods, spacing);
	}

	public LogicTreeHazardCompare(SolutionLogicTree solLogicTree, LogicTree<?> tree, File mapsZipFile,
			ReturnPeriods[] rps, double[] periods, double spacing) throws IOException {
		this.solLogicTree = solLogicTree;
		this.rps = rps;
		this.periods = periods;

		logCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3d, 1d);
		spreadCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0, 1d);
		spreadDiffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1d, 1d);
		spreadDiffCPT.setNanColor(Color.GRAY);
		pDiffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-100d, 100d);
		pDiffCPT.setNanColor(Color.GRAY);
		percentileCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 100d);
		percentileCPT.setNanColor(Color.BLACK);
		
		ZipFile zip = new ZipFile(mapsZipFile);
		
		branches = tree.getBranches();
		weights = new ArrayList<>();
		
		maps = HashBasedTable.create();
		
		BranchWeightProvider weightProv = tree.getWeightProvider();
		
		GriddedRegion gridReg = null;
		for (int i=0; i<branches.size(); i++) {
			LogicTreeBranch<?> branch = branches.get(i);
			double weight = weightProv.getWeight(branch);
			Preconditions.checkState(weight > 0, "Bad weight=%s for branch %s, weightProv=%s",
					weight, branch, weightProv.getClass().getName());
			weights.add(weight);
			totWeight += weight;
			
			if (gridReg == null) {
				FaultSystemSolution sol = solLogicTree.forBranch(branch);
				Region region = ReportMetadata.detectRegion(sol);
				if (region == null)
					// just use bounding box
					region = RupSetMapMaker.buildBufferedRegion(sol.getRupSet().getFaultSectionDataList());
				gridReg = new GriddedRegion(region, spacing, GriddedRegion.ANCHOR_0_0);
				
				mapper = new SolHazardMapCalc(sol, null, gridReg, periods);
			}
			
			System.out.println("Processing maps for "+branch);
			
			String dirName = branch.buildFileName();
			
			for (ReturnPeriods rp : rps) {
				for (double period : periods) {
					
					ZipEntry entry = zip.getEntry(dirName+"/"+MPJ_LogicTreeHazardCalc.mapPrefix(period, rp)+".txt");
					
					GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
					BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
					String line = bRead.readLine();
					int index = 0;
					while (line != null) {
						line = line.trim();
						if (!line.startsWith("#")) {
							StringTokenizer tok = new StringTokenizer(line);
							double lon = Double.parseDouble(tok.nextToken());
							double lat = Double.parseDouble(tok.nextToken());
							double val = Double.parseDouble(tok.nextToken());
							Location loc = new Location(lat, lon);
							Preconditions.checkState(LocationUtils.areSimilar(loc, gridReg.getLocation(index)));
							xyz.set(index++, val);
						}
						line = bRead.readLine();
					}
					Preconditions.checkState(index == gridReg.getNodeCount());
					GriddedGeoDataSet[] rpPerMaps = maps.get(rp, period);
					if (rpPerMaps == null) {
						rpPerMaps = new GriddedGeoDataSet[branches.size()];
						maps.put(rp, period, rpPerMaps);
					}
					rpPerMaps[i] = xyz;
				}
			}
		}
		
		zip.close();
		
		System.out.println("Total weight: "+totWeight);
	}
	
	private GriddedGeoDataSet buildMean(GriddedGeoDataSet[] maps) {
		GriddedGeoDataSet avg = new GriddedGeoDataSet(maps[0].getRegion(), false);
		
		for (int i=0; i<avg.size(); i++) {
			double val = 0d;
			for (int j=0; j<maps.length; j++)
				val += maps[j].get(i)*weights.get(j);
			val /= totWeight;
			avg.set(i, val);
		}
		
		return avg;
	}
	
	private GriddedGeoDataSet buildMean(List<GriddedGeoDataSet> maps, List<Double> weights) {
		GriddedGeoDataSet avg = new GriddedGeoDataSet(maps.get(0).getRegion(), false);
		
		double totWeight = 0d;
		for (Double weight : weights)
			totWeight += weight;
		
		for (int i=0; i<avg.size(); i++) {
			double val = 0d;
			for (int j=0; j<maps.size(); j++)
				val += maps.get(j).get(i)*weights.get(j);
			val /= totWeight;
			avg.set(i, val);
		}
		
		return avg;
	}
	
	private GriddedGeoDataSet buildMin(GriddedGeoDataSet[] maps) {
		GriddedGeoDataSet min = new GriddedGeoDataSet(maps[0].getRegion(), false);
		
		for (int i=0; i<min.size(); i++) {
			double val = Double.POSITIVE_INFINITY;
			for (int j=0; j<maps.length; j++)
				val = Math.min(val, maps[j].get(i));
			min.set(i, val);
		}
		
		return min;
	}
	
	private GriddedGeoDataSet buildMax(GriddedGeoDataSet[] maps) {
		GriddedGeoDataSet max = new GriddedGeoDataSet(maps[0].getRegion(), false);
		
		for (int i=0; i<max.size(); i++) {
			double val = Double.NEGATIVE_INFINITY;
			for (int j=0; j<maps.length; j++)
				val = Math.max(val, maps[j].get(i));
			max.set(i, val);
		}
		
		return max;
	}
	
	private GriddedGeoDataSet buildPDiffFromRange(GriddedGeoDataSet min, GriddedGeoDataSet max, GriddedGeoDataSet comp) {
		GriddedGeoDataSet diff = new GriddedGeoDataSet(min.getRegion(), false);
		
		for (int i=0; i<max.size(); i++) {
			double compVal = comp.get(i);
			double minVal = min.get(i);
			double maxVal = max.get(i);
			double pDiff;
			if (compVal >= minVal && compVal <= maxVal)
				pDiff = 0d;
			else if (compVal < minVal)
				pDiff = 100d*(compVal-minVal)/minVal;
			else
				pDiff = 100d*(compVal-maxVal)/maxVal;
			diff.set(i, pDiff);
		}
		
		return diff;
	}
	
	private GriddedGeoDataSet buildSpread(GriddedGeoDataSet min, GriddedGeoDataSet max) {
		GriddedGeoDataSet diff = new GriddedGeoDataSet(min.getRegion(), false);
		
		for (int i=0; i<max.size(); i++) {
			double minVal = min.get(i);
			double maxVal = max.get(i);
			diff.set(i, maxVal - minVal);
		}
		
		return diff;
	}
	
	private GriddedGeoDataSet buildPDiff(GriddedGeoDataSet ref, GriddedGeoDataSet comp) {
		GriddedGeoDataSet diff = new GriddedGeoDataSet(ref.getRegion(), false);
		
		for (int i=0; i<ref.size(); i++) {
			double compVal = comp.get(i);
			double refVal = ref.get(i);
			double pDiff = 100d*(compVal-refVal)/refVal;
			diff.set(i, pDiff);
		}
		
		return diff;
	}
	
	private ArbDiscrEmpiricalDistFunc[] buildArbDiscrEmpiricalDists(GriddedGeoDataSet[] maps, List<Double> weights) {
		return buildArbDiscrEmpiricalDists(List.of(maps), weights);
	}
	
	private ArbDiscrEmpiricalDistFunc[] buildArbDiscrEmpiricalDists(List<GriddedGeoDataSet> maps, List<Double> weights) {
		Preconditions.checkState(maps.size() == weights.size());
		ArbDiscrEmpiricalDistFunc[] ret = new ArbDiscrEmpiricalDistFunc[maps.get(0).size()];
		
		double totWeight = 0d;
		for (double weight : weights)
			totWeight += weight;
		
		for (int i=0; i<ret.length; i++) {
			ret[i] = new ArbDiscrEmpiricalDistFunc();
			for (int j=0; j<maps.size(); j++) {
				double weight = weights.get(j)/totWeight;
				Preconditions.checkState(weight > 0d, "bad weight: %s", weight);
				GriddedGeoDataSet map = maps.get(j);
				
				double val = map.get(i);
				ret[i].set(val, weight);
			}
		}
		
		for (int i=0; i<ret.length; i++) {
			float sum = (float)ret[i].calcSumOfY_Vals();
			Preconditions.checkState(sum == 1f, "sum of distribution y-values is not 1: %s", sum);
		}
		
		return ret;
	}
	
	private GriddedGeoDataSet calcMapAtPercentile(ArbDiscrEmpiricalDistFunc[] cellDists, GriddedRegion gridReg,
			double percentile) {
		Preconditions.checkState(gridReg.getNodeCount() == cellDists.length);
		GriddedGeoDataSet ret = new GriddedGeoDataSet(gridReg, false);
		
		for (int i=0; i<ret.size(); i++) {
			double val = cellDists[i].getInterpolatedFractile(percentile/100d);
			ret.set(i, val);
		}
		
		return ret;
	}
	
	private GriddedGeoDataSet calcPercentileWithinDist(ArbDiscrEmpiricalDistFunc[] cellDists, GriddedGeoDataSet comp) {
		Preconditions.checkState(comp.size() == cellDists.length);
		GriddedGeoDataSet ret = new GriddedGeoDataSet(comp.getRegion(), false);
		
		for (int i=0; i<ret.size(); i++) {
			double compVal = comp.get(i);
			
			DiscretizedFunc ncd = cellDists[i].getNormalizedCumDist();
			
			double percentile;
			if (compVal < ncd.getMinX() || compVal > ncd.getMaxX()) {
				percentile = Double.NaN;
			} else {
				percentile = 100d * ncd.getInterpolatedY(compVal);
			}
			ret.set(i, percentile);
		}
		
		return ret;
	}
	
//	private GriddedGeoDataSet calcPercentile(GriddedGeoDataSet[] maps, GriddedGeoDataSet comp) {
//		return calcPercentile(List.of(maps), weights, comp);
//	}
//	 
//	private GriddedGeoDataSet calcPercentile(List<GriddedGeoDataSet> maps, List<Double> weights, GriddedGeoDataSet comp) {
//		Preconditions.checkState(maps.size() == weights.size());
//		GriddedGeoDataSet ret = new GriddedGeoDataSet(maps.get(0).getRegion(), false);
//		
//		double totWeight;
//		if (weights == this.weights) {
//			totWeight = this.totWeight;
//		} else {
//			totWeight = 0d;
//			for (double weight : weights)
//				totWeight += weight;
//		}
//		
//		for (int i=0; i<ret.size(); i++) {
//			double compVal = comp.get(i);
//			double weightAbove = 0d; // weight of tree thats above comp value
//			double weightEqual = 0d; // weight of tree thats equal to comp value
//			double weightBelow = 0d; // weight of tree thats below comp value
//			int numAbove = 0;
//			int numBelow = 0;
//			for (int j=0; j<maps.size(); j++) {
//				double val = maps.get(j).get(i);
//				double weight = weights.get(j);
//				if (((float) compVal == (float)val) || (Double.isNaN(compVal)) && Double.isNaN(val)) {
//					weightEqual += weight;
//				} else if (compVal < val) {
//					weightAbove += weight;
//					numAbove++;
//				} else if (compVal > val) {
//					numBelow++;
//					weightBelow += weight;
//				}
//			}
//			if (weightEqual != 0d) {
//				// redistribute any exactly equal to either side
//				weightAbove += 0.5*weightEqual;
//				weightBelow += 0.5*weightEqual;
//			}
//			// normalize by total weight
//			weightAbove /= totWeight;
//			weightBelow /= totWeight;
//			double percentile;
//			if (numAbove == maps.size() || numBelow == maps.size())
//				percentile = Double.NaN;
//			else
//				percentile = 100d*weightBelow;
////			if (numAbove == maps.size() || numBelow == maps.size()) {
////				percentile = Double.NaN;
////			} else if (weightAbove > weightBelow) {
////				// more of the distribution is above my value
////				// this means that the percentile is <50
////				percentile = 100d*weightBelow/totWeight;
////			} else {
////				// more of the distribution is below my value
////				// this means that the percentile is >50
////				percentile = 100d*weightAbove/totWeight;
////				percentile = 100d*(1d - (weightAbove/totWeight));
////			}
//			ret.set(i, percentile);
//		}
//		
//		return ret;
//	}
	
	public void buildReport(File outputDir, String name, LogicTreeHazardCompare comp, String compName) throws IOException {
		List<String> lines = new ArrayList<>();
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		lines.add("# "+name+" Hazard Maps");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "_[(top)](#table-of-contents)_";
		
		for (double period : periods) {
			String perLabel, perPrefix;
			if (period == 0d) {
				perLabel = "PGA (g)";
				perPrefix = "pga";
			} else {
				perLabel = (float)period+"s SA";
				perPrefix = (float)period+"s";
			}
			
			for (ReturnPeriods rp : rps) {
				String label = perLabel+", "+rp.label;
				String prefix = perPrefix+"_"+rp.name();
				
				System.out.println(label);
				
				lines.add("## "+label);
				lines.add(topLink);
				lines.add("");
				
				GriddedGeoDataSet[] maps = this.maps.get(rp, period);
				Preconditions.checkNotNull(maps);
				for (int i=0; i<maps.length; i++)
					Preconditions.checkNotNull(maps[i], "map %s is null", i);
				
				GriddedRegion region = maps[0].getRegion();
				
				ArbDiscrEmpiricalDistFunc[] mapArbDiscrs = buildArbDiscrEmpiricalDists(maps, weights);
				GriddedGeoDataSet mean = buildMean(maps);
				GriddedGeoDataSet median = calcMapAtPercentile(mapArbDiscrs, region, 50d);
				GriddedGeoDataSet max = buildMax(maps);
				GriddedGeoDataSet min = buildMin(maps);
				GriddedGeoDataSet spread = buildSpread(log10(min), log10(max));
				
				File hazardCSV = new File(resourcesDir, prefix+".csv");
				writeHazardCSV(hazardCSV, mean, median, min, max);
				
//				table.addLine(meanMinMaxSpreadMaps(mean, min, max, spread, name, label, prefix, resourcesDir));
				
				ArbDiscrEmpiricalDistFunc[] cmapArbDiscrs = null;
				GriddedGeoDataSet cmean = null;
				GriddedGeoDataSet cmedian = null;
				GriddedGeoDataSet cmin = null;
				GriddedGeoDataSet cmax = null;
				GriddedGeoDataSet cspread = null;
				
				if (comp != null) {
					GriddedGeoDataSet[] cmaps = comp.maps.get(rp, period);
					Preconditions.checkNotNull(cmaps);
					for (int i=0; i<cmaps.length; i++)
						Preconditions.checkNotNull(cmaps[i], "map %s is null", i);

					cmapArbDiscrs = buildArbDiscrEmpiricalDists(cmaps, comp.weights);
					cmean = comp.buildMean(cmaps);
					cmedian = calcMapAtPercentile(cmapArbDiscrs, region, 50d);
					cmax = comp.buildMax(cmaps);
					cmin = comp.buildMin(cmaps);
					cspread = comp.buildSpread(log10(cmin), log10(cmax));
//					table.addLine(meanMinMaxSpreadMaps(cmean, cmin, cmax, cspread, compName, label, prefix+"_comp", resourcesDir));
					
					File compHazardCSV = new File(resourcesDir, prefix+"_comp.csv");
					writeHazardCSV(compHazardCSV, cmean, cmedian, cmin, cmax);
					
					lines.add("Download Hazard CSVs: ["+hazardCSV.getName()+"]("+resourcesDir.getName()+"/"+hazardCSV.getName()
						+")  ["+compHazardCSV.getName()+"]("+resourcesDir.getName()+"/"+compHazardCSV.getName()+")");
				} else {
					lines.add("Download Hazard CSV: ["+hazardCSV.getName()+"]("+resourcesDir.getName()+"/"+hazardCSV.getName()+")");
				}
				lines.add("");
				
//				table.invert();
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				
				File meanMapFile = mapper.plotMap(resourcesDir, prefix+"_mean", log10(mean),
						logCPT, TITLES ? name : " ", "Log10 Weighted-Average, "+label);
				GriddedGeoDataSet.writeXYZFile(mean, new File(resourcesDir, prefix+"_mean.xyz"));
				File medianMapFile = mapper.plotMap(resourcesDir, prefix+"_median", log10(median),
						logCPT, TITLES ? name : " ", "Log10 Weighted-Median, "+label);
				GriddedGeoDataSet.writeXYZFile(median, new File(resourcesDir, prefix+"_median.xyz"));
				
				if (cmean == null) {
					// no comparison
					table.initNewLine();
					table.addColumn(MarkdownUtils.boldCentered("Weighted-Average"));
					table.addColumn(MarkdownUtils.boldCentered("Weighted-Median"));
					table.finalizeLine().initNewLine();
					table.addColumn("![Mean Map]("+resourcesDir.getName()+"/"+meanMapFile.getName()+")");
					table.addColumn("![Median Map]("+resourcesDir.getName()+"/"+medianMapFile.getName()+")");
					table.finalizeLine();
				} else {
					// comparison
					File cmeanMapFile = mapper.plotMap(resourcesDir, prefix+"_comp_mean", log10(cmean),
							logCPT, TITLES ? compName : " ", "Log10 Weighted-Average, "+label);
					GriddedGeoDataSet.writeXYZFile(cmean, new File(resourcesDir, prefix+"_comp_mean.xyz"));
					File cmedianMapFile = mapper.plotMap(resourcesDir, prefix+"_comp_median", log10(cmedian),
							logCPT, TITLES ? compName : " ", "Log10 Weighted-Median, "+label);
					GriddedGeoDataSet.writeXYZFile(cmedian, new File(resourcesDir, prefix+"_comp_median.xyz"));
					
					table.initNewLine();
					table.addColumn(MarkdownUtils.boldCentered("Primary Weighted-Average"));
					table.addColumn(MarkdownUtils.boldCentered("Comparison Weighted-Average"));
					table.finalizeLine().initNewLine();
					table.addColumn("![Mean Map]("+resourcesDir.getName()+"/"+meanMapFile.getName()+")");
					table.addColumn("![Mean Map]("+resourcesDir.getName()+"/"+cmeanMapFile.getName()+")");
					table.finalizeLine();
					
					addDiffInclExtremesLines(mean, name, min, max, cmean, compName, prefix+"_mean", resourcesDir, table, "Mean", label);
					
					table.initNewLine();
					table.addColumn(MarkdownUtils.boldCentered("Primary Weighted-Median"));
					table.addColumn(MarkdownUtils.boldCentered("Comparison Weighted-Median"));
					table.finalizeLine().initNewLine();
					table.addColumn("![Median Map]("+resourcesDir.getName()+"/"+medianMapFile.getName()+")");
					table.addColumn("![Median Map]("+resourcesDir.getName()+"/"+cmedianMapFile.getName()+")");
					table.finalizeLine();
					
					addDiffInclExtremesLines(median, name, min, max, cmedian, compName, prefix+"_median", resourcesDir, table, "Median", label);
					
					table.addLine(MarkdownUtils.boldCentered("Comparison Mean Percentile"),
							MarkdownUtils.boldCentered("Comparison Median Percentile"));
					GriddedGeoDataSet cMeanPercentile = calcPercentileWithinDist(mapArbDiscrs, cmean);
					GriddedGeoDataSet cMedianPercentile = calcPercentileWithinDist(mapArbDiscrs, cmedian);
					table.initNewLine();
					File map = mapper.plotMap(resourcesDir, prefix+"_comp_mean_percentile", cMeanPercentile,
							percentileCPT, TITLES ? name+" vs "+compName : " ", "Comparison Mean %-ile, "+label);
					table.addColumn("![Mean Percentile Map]("+resourcesDir.getName()+"/"+map.getName()+")");
					map = mapper.plotMap(resourcesDir, prefix+"_comp_median_percentile", cMedianPercentile,
							percentileCPT, TITLES ? name+" vs "+compName : " ", "Comparison Median %-ile, "+label);
					table.addColumn("![Median Percentile Map]("+resourcesDir.getName()+"/"+map.getName()+")");
//					GriddedGeoDataSet spreadDiff = new GriddedGeoDataSet(spread.getRegion(), false);
//					for (int i=0; i<spreadDiff.size(); i++)
//						spreadDiff.set(i, spread.get(i)-cspread.get(i));
//					map = mapper.plotMap(resourcesDir, prefix+"_comp_spread_diff", spreadDiff,
//							spreadDiffCPT, name+" vs "+compName, "Log10 Spread Difference, "+label);
//					table.addColumn("![Spread Diff]("+resourcesDir.getName()+"/"+map.getName()+")");
					table.finalizeLine();
				}
				
				lines.addAll(table.build());
				lines.add("");

				table = MarkdownUtils.tableBuilder();
				if (cmean == null) {
					// no comparison
					table.addLine(MarkdownUtils.boldCentered("Minimum"), MarkdownUtils.boldCentered("Maximum"),
							MarkdownUtils.boldCentered("Spread"));
					
					File minMapFile = mapper.plotMap(resourcesDir, prefix+"_min", log10(min),
							logCPT, TITLES ? name : " ", "Log10 Minimum, "+label);
					File maxMapFile = mapper.plotMap(resourcesDir, prefix+"_max", log10(max),
							logCPT, TITLES ? name : " ", "Log10 Maximum, "+label);
					File spreadMapFile = mapper.plotMap(resourcesDir, prefix+"_spread", spread,
							spreadCPT, TITLES ? name : " ", "Log10 Spread, "+label);
					
					table.initNewLine();
					table.addColumn("![Min Map]("+resourcesDir.getName()+"/"+minMapFile.getName()+")");
					table.addColumn("![Max Map]("+resourcesDir.getName()+"/"+maxMapFile.getName()+")");
					table.addColumn("![Spread Map]("+resourcesDir.getName()+"/"+spreadMapFile.getName()+")");
					table.finalizeLine();
				} else {
					// comparison
					addMapCompDiffLines(min, name, cmin, compName, prefix+"_min", resourcesDir, table,
							"Minimum", "Minimum "+label, logCPT, false);
					addMapCompDiffLines(max, name, cmax, compName, prefix+"_max", resourcesDir, table,
							"Maximum", "Maximum"+label, logCPT, false);
					addMapCompDiffLines(spread, name, cspread, compName, prefix+"_spread", resourcesDir, table,
							"Spread", "Log10 Spread "+label, spreadCPT, true);
				}
				
				lines.addAll(table.build());
				lines.add("");
				
				lines.add("### "+label+" Logic Tree Comparisons");
				lines.add(topLink); lines.add("");
				
				// plot mean percentile
				table = MarkdownUtils.tableBuilder();
				
				GriddedGeoDataSet meanPercentile = calcPercentileWithinDist(mapArbDiscrs, mean);
				File meanPercentileMap = mapper.plotMap(resourcesDir, prefix+"_mean_percentile",
						meanPercentile, percentileCPT, TITLES ? "Branch-Averaged Percentiles" : " ",
						"Branch Averaged %-ile, "+label);
//				GriddedGeoDataSet median = calc
				table.addLine("Branch Averaged Map", "Branch Averaged Percentiles");
				table.addLine("![BA map]("+resourcesDir.getName()+"/"+meanMapFile.getName()+")",
						"![BA percentiles]("+resourcesDir.getName()+"/"+meanPercentileMap.getName()+")");
				table.addLine(MarkdownUtils.boldCentered("Median Map"), MarkdownUtils.boldCentered("Mean vs Median"));
				GriddedGeoDataSet meanMedDiff = buildPDiff(median, mean);
				File meanMedDiffMap = mapper.plotMap(resourcesDir, prefix+"_mean_median_diff",
						meanMedDiff, pDiffCPT, TITLES ? name : " ", "Mean - Median, % Difference, "+label, true);
				table.addLine("![BA Median map]("+resourcesDir.getName()+"/"+medianMapFile.getName()+")",
						"![Median vs Mean]("+resourcesDir.getName()+"/"+meanMedDiffMap.getName()+")");
//				File meanMap = new File(outputDir, prefix+"_median.png");
				lines.add("Branched-average hazard can be dominated by outlier branches. This map shows the percentile"
						+ " at which the branch averaged map lies; areas far from the 50-th percentile are likely "
						+ "outlier-dominated.");
				lines.add("");
				lines.addAll(table.build());
				lines.add("");
				
				for (LogicTreeLevel<?> level : solLogicTree.getLogicTree().getLevels()) {
					HashMap<LogicTreeNode, List<GriddedGeoDataSet>> choiceMaps = new HashMap<>();
					HashMap<LogicTreeNode, List<Double>> choiceWeights = new HashMap<>();
					for (int i=0; i<branches.size(); i++) {
						LogicTreeBranch<?> branch = branches.get(i);
						LogicTreeNode choice = branch.getValue(level.getType());
						List<GriddedGeoDataSet> myChoiceMaps = choiceMaps.get(choice);
						if (myChoiceMaps == null) {
							myChoiceMaps = new ArrayList<>();
							choiceMaps.put(choice, myChoiceMaps);
							choiceWeights.put(choice, new ArrayList<>());
						}
						myChoiceMaps.add(maps[i]);
						choiceWeights.get(choice).add(weights.get(i));
					}
					if (choiceMaps.size() > 1) {
						lines.add("#### "+level.getName());
						lines.add(topLink); lines.add("");
						HashMap<LogicTreeNode, GriddedGeoDataSet> choiceMeans = new HashMap<>();
						for (LogicTreeNode choice : choiceMaps.keySet())
							choiceMeans.put(choice, buildMean(choiceMaps.get(choice), choiceWeights.get(choice)));
						
						table = MarkdownUtils.tableBuilder();
						table.initNewLine();
						table.addColumn("**Choice**").addColumn("**Vs Mean**");
						List<LogicTreeNode> choices = new ArrayList<>(choiceMaps.keySet());
						Collections.sort(choices, nodeNameCompare);
						for (LogicTreeNode choice : choices)
							table.addColumn("**Vs "+choice.getShortName()+"**");
						table.finalizeLine();
						
						File choicesCSV = new File(resourcesDir, prefix+"_"+level.getShortName().replaceAll("\\W+", "_")+".csv");
						writeChoiceHazardCSV(choicesCSV, mean, choices, choiceMeans);
						lines.add("Download Choice Hazard CSV: ["+choicesCSV.getName()+"]("+resourcesDir.getName()+"/"+choicesCSV.getName()+")");
						lines.add("");
						
						TableBuilder mapTable = MarkdownUtils.tableBuilder();
						mapTable.addLine("", "Choice Mean vs Full Mean", "Choice Percentile in Full Dist", 
								"Choice Percentile in Dist Without");
						
						MinMaxAveTracker runningDiffAvg = new MinMaxAveTracker();
						MinMaxAveTracker runningAbsDiffAvg = new MinMaxAveTracker();
						
						for (LogicTreeNode choice : choices) {
							table.initNewLine().addColumn("**"+choice.getShortName()+"**");
							
							GriddedGeoDataSet choiceMap = choiceMeans.get(choice);
							table.addColumn(mapPDiffStr(choiceMap, mean, null, null));
							
							for (LogicTreeNode oChoice : choices) {
								if (choice == oChoice)
									table.addColumn("");
								else
									table.addColumn(mapPDiffStr(choiceMap, choiceMeans.get(oChoice),
											runningDiffAvg, runningAbsDiffAvg));
							}
							
							table.finalizeLine();
							
							// now maps
							GriddedGeoDataSet pDiff = buildPDiff(mean, choiceMap);
							
							mapTable.initNewLine().addColumn("**"+choice.getShortName()+"**");
							File map = mapper.plotMap(resourcesDir, prefix+"_"+choice.getFilePrefix()+"_pDiff",
									pDiff, pDiffCPT, TITLES ? choice.getShortName()+" Comparison" : " ",
									choice.getShortName()+" - Mean, % Difference, "+label, true);
							mapTable.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
							GriddedGeoDataSet percentile = calcPercentileWithinDist(mapArbDiscrs, choiceMap);
							map = mapper.plotMap(resourcesDir, prefix+"_"+choice.getFilePrefix()+"_percentile",
									percentile, percentileCPT, TITLES ? choice.getShortName()+" Comparison" : " ",
									choice.getShortName()+" %-ile, "+label);
							mapTable.addColumn("![Percentile Map]("+resourcesDir.getName()+"/"+map.getName()+")");
							List<GriddedGeoDataSet> mapsWithout = new ArrayList<>();
							List<Double> weightsWithout = new ArrayList<>();
							for (int i=0; i<branches.size(); i++) {
								LogicTreeBranch<?> branch = branches.get(i);
								if (!branch.hasValue(choice)) {
									mapsWithout.add(maps[i]);
									weightsWithout.add(weights.get(i));
								}
							}
							Preconditions.checkState(!mapsWithout.isEmpty());
							ArbDiscrEmpiricalDistFunc[] mapsWithoutArbDiscrs = buildArbDiscrEmpiricalDists(mapsWithout, weightsWithout);
							GriddedGeoDataSet percentileWithout = calcPercentileWithinDist(mapsWithoutArbDiscrs, choiceMap);
							map = mapper.plotMap(resourcesDir, prefix+"_"+choice.getFilePrefix()+"_percentile_without",
									percentileWithout, percentileCPT, TITLES ? choice.getShortName()+" Comparison" : " ",
									choice.getShortName()+" %-ile, "+label);
							mapTable.addColumn("![Percentile Map]("+resourcesDir.getName()+"/"+map.getName()+")");
							mapTable.finalizeLine();
						}
						lines.add("");
						lines.add("Mean absolute difference: "+twoDigits.format(runningAbsDiffAvg.getAverage())+"%");
						lines.add("");
						lines.addAll(table.build());
						lines.add("");
						lines.addAll(mapTable.build());
						lines.add("");
					}
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private void writeHazardCSV(File outputFile, GriddedGeoDataSet mean, GriddedGeoDataSet median,
			GriddedGeoDataSet min, GriddedGeoDataSet max) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		csv.addLine("Location Index", "Latitutde", "Longitude", "Weighted Mean", "Weighted Median", "Min", "Max");
		
		GriddedRegion reg = mean.getRegion();
		for (int i=0; i<reg.getNodeCount(); i++) {
			Location loc = reg.getLocation(i);
			
			csv.addLine(i+"", (float)loc.getLatitude()+"", (float)loc.getLongitude()+"", mean.get(i)+"",
					median.get(i)+"", min.get(i)+"", max.get(i)+"");
		}
		
		csv.writeToFile(outputFile);
	}
	
	private void writeChoiceHazardCSV(File outputFile, GriddedGeoDataSet mean, List<LogicTreeNode> choices,
			Map<LogicTreeNode, GriddedGeoDataSet> choiceMeans) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("Location Index");
		header.add("Latitutde");
		header.add("Longitude");
		header.add("Weighted Mean");
		for (LogicTreeNode choice : choices)
			header.add(choice.getShortName());
		csv.addLine(header);
		
		GriddedRegion reg = mean.getRegion();
		for (int i=0; i<reg.getNodeCount(); i++) {
			Location loc = reg.getLocation(i);
			
			List<String> line = new ArrayList<>();
			line.add(i+"");
			line.add((float)loc.getLatitude()+"");
			line.add((float)loc.getLongitude()+"");
			line.add(mean.get(i)+"");
			for (LogicTreeNode choice : choices)
				line.add(choiceMeans.get(choice).get(i)+"");
			
			csv.addLine(line);
		}
		
		csv.writeToFile(outputFile);
	}
	
	private static String mapPDiffStr(GriddedGeoDataSet map, GriddedGeoDataSet ref,
			MinMaxAveTracker runningDiffAvg, MinMaxAveTracker runningAbsDiffAvg) {
		double mean = 0d;
		double meanAbs = 0d;
		for (int i=0; i<map.size(); i++) {
			double z1 = map.get(i);
			double z2 = ref.get(i);
			double pDiff = 100d*(z1-z2)/z2;
			if (z1 == z2)
				pDiff = 0d;
			if (Double.isFinite(pDiff)) {
				mean += pDiff;
				meanAbs += Math.abs(pDiff);
				if (runningDiffAvg != null)
					runningDiffAvg.addValue(pDiff);
				if (runningAbsDiffAvg != null)
					runningAbsDiffAvg.addValue(Math.abs(pDiff));
			}
		}
		mean /= (double)map.size();
		meanAbs /= (double)map.size();
		
		return "Mean: "+twoDigits.format(mean)+"%, Mean Abs: "+twoDigits.format(meanAbs)+"%";
	}
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	
	private static final Comparator<LogicTreeNode> nodeNameCompare = new Comparator<LogicTreeNode>() {

		@Override
		public int compare(LogicTreeNode o1, LogicTreeNode o2) {
			return o1.getShortName().compareTo(o2.getShortName());
		}
	};
	
	private static GriddedGeoDataSet log10(GriddedGeoDataSet map) {
		map = map.copy();
		map.log10();
		return map;
	}
	
	private void addDiffInclExtremesLines(GriddedGeoDataSet primary, String name, GriddedGeoDataSet min,
			GriddedGeoDataSet max, GriddedGeoDataSet comparison, String compName, String prefix, File resourcesDir,
			TableBuilder table, String type, String label) throws IOException {
		
		table.addLine(MarkdownUtils.boldCentered(type+" Difference"),
				MarkdownUtils.boldCentered("Comparison "+type+" Diff From Extremes"));
		
		table.initNewLine();
		
		GriddedGeoDataSet pDiff = buildPDiff(comparison, primary);
		File map = mapper.plotMap(resourcesDir, prefix+"_comp_pDiff", pDiff, pDiffCPT, TITLES ? name+" vs "+compName : " ",
				(TITLES ? "Primary - Comparison, " : "")+"% Difference, "+label, true);
		table.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		GriddedGeoDataSet pDiffFromRange = buildPDiffFromRange(min, max, comparison);
		map = mapper.plotMap(resourcesDir, prefix+"_comp_pDiff_range", pDiffFromRange, pDiffCPT, TITLES ? name+" vs "+compName : " ",
				"Comparison - Extremes, % Difference, "+label, true);
		table.addColumn("![Range Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		table.finalizeLine();
	}
	
	private void addMapCompDiffLines(GriddedGeoDataSet primary, String name, GriddedGeoDataSet comparison,
			String compName, String prefix, File resourcesDir, TableBuilder table, String type, String label,
			CPT cpt, boolean spread) throws IOException {

		CPT diffCPT;
		GriddedGeoDataSet diff;
		String diffLabel;
		if (spread) {
			diff = new GriddedGeoDataSet(primary.getRegion(), false);
			for (int i=0; i<diff.size(); i++)
				diff.set(i, primary.get(i)-comparison.get(i));
			diffCPT = spreadDiffCPT;
			diffLabel = "Primary - Comparison, "+label;
		} else {
			diff = buildPDiff(comparison, primary);;
			primary = log10(primary);
			comparison = log10(comparison);
			diffCPT = pDiffCPT;
			diffLabel = "Primary - Comparison, % Difference, "+label;
		}
		
		table.addLine(MarkdownUtils.boldCentered("Primary "+type), MarkdownUtils.boldCentered("Comparison "+type),
				MarkdownUtils.boldCentered(spread ? "Difference" : "% Difference"));
		
		table.initNewLine();
		File map = mapper.plotMap(resourcesDir, prefix, primary, cpt, name, (spread ? "" : "Log10 ")+label);
		table.addColumn("!["+type+"]("+resourcesDir.getName()+"/"+map.getName()+")");
		map = mapper.plotMap(resourcesDir, prefix+"_comp", comparison, cpt, name, (spread ? "" : "Log10 ")+label);
		table.addColumn("!["+type+"]("+resourcesDir.getName()+"/"+map.getName()+")");
		map = mapper.plotMap(resourcesDir, prefix+"_comp_pDiff", diff, diffCPT, name+" vs "+compName, diffLabel, !spread);
		table.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		table.finalizeLine();
	}

}
