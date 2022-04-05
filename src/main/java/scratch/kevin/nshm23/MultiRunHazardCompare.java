package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipFile;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class MultiRunHazardCompare {
	
	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		List<File> runDirs = new ArrayList<>();
		List<String> runNames = new ArrayList<>();
		
		List<File> fixedComps = new ArrayList<>();
		List<String> fixedCompNames = new ArrayList<>();
		List<String> fixedCompDirNames = new ArrayList<>();
		
		boolean replot = false;
		
		fixedComps.add(new File(mainDir, "2021_11_30-u3_branches-orig_calcs-5h"));
		fixedCompNames.add("UCERF3 As Published");
		fixedCompDirNames.add("hazard_maps_comp_u3_as_published");
		
		runDirs.add(new File(mainDir, "2022_03_25-u3_branches-orig_calc_params-FM3_1"));
		runNames.add("UCERF3 Reproduction");
		
		runDirs.add(new File(mainDir, "2022_03_25-u3_branches-orig_calc_params-new_avg-converged-FM3_1-2000ip"));
		runNames.add("UCERF3 Converged, NewAvg");
		
		runDirs.add(new File(mainDir, "2022_03_25-u3_branches-orig_calc_params-new_avg-converged-noWL-FM3_1-2000ip"));
		runNames.add("UCERF3 Converged, NewAvg, No WL");
		
		runDirs.add(new File(mainDir, "2022_03_25-u3_branches-orig_calc_params-new_avg-converged-noWL-new_perturb-FM3_1-2000ip"));
		runNames.add("UCERF3 Converged, NewAvg, No WL, Exp Perturb");
		
		runDirs.add(new File(mainDir, "2022_03_24-u3_branches-FM3_1-2000ip"));
		runNames.add("UCERF3 New Anneal");
		fixedComps.add(runDirs.get(runDirs.size()-1));
		fixedCompNames.add(runNames.get(runNames.size()-1));
		fixedCompDirNames.add("hazard_maps_comp_u3_new_anneal");
		
		runDirs.add(new File(mainDir, "2022_03_25-u3_branches-thinned_0.1-FM3_1-2000ip"));
		runNames.add("UCERF3 New Anneal Reduced RS");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-no_seg-FM3_1-U3RedRupSet-DsrUni-TotNuclRate-SubB1-2000ip"));
		runNames.add("NSHM23 Draft, U3 Reduced RS, No Seg");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-no_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, No Seg");
		
		runDirs.add(new File(mainDir, "2022_03_30-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-NoAdj-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Shift-1km Seg, No Adj");
		
		runDirs.add(new File(mainDir, "2022_03_24-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Thresh-Avg, Shift-1km");
		fixedComps.add(runDirs.get(runDirs.size()-1));
		fixedCompNames.add(runNames.get(runNames.size()-1));
		fixedCompDirNames.add("hazard_maps_comp_nshm23_draft");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-10000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Thresh-Avg, Shift-1km, Extra Converged");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-CappedRdst-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Capped Redistribution, Shift-1km");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-NuclMFD-SubB1-JumpProb-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Thresh-Avg, Shift-1km, Nucl MFD");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SysAvg-JumpProb-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Thresh-Avg, Shift-1km, Syst-Avg Sub-Seis Reduction");
		
		runDirs.add(new File(mainDir, "2022_03_31-nshm23_u3_hybrid_branches-no_scale_adj_mfds-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip"));
		runNames.add("NSHM23 Draft, Coulomb RS, Thresh-Avg, Shift-1km, No Scale Adj MFDs");
		
		runDirs.add(new File(mainDir, "2022_03_25-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip-branch-translated-min3km"));
		runNames.add("NSHM23 Draft, Coulomb RS, Strict Segmentation Min 3km");
		
		List<File> allDirs = new ArrayList<>();
		allDirs.addAll(runDirs);
		allDirs.addAll(fixedComps);
		System.out.println("Making sure all inputs exist...");
		for (File dir : allDirs) {
			System.out.println("\t"+dir.getName());
			Preconditions.checkState(dir.exists(), "%s doesn't exist", dir.getAbsoluteFile());
			File resultsFile = new File(dir, "results.zip");
			Preconditions.checkState(resultsFile.exists(), "%s doesn't exist", resultsFile.getAbsoluteFile());
			new ZipFile(resultsFile).close();
			File hazardFile = new File(dir, "results_hazard.zip");
			Preconditions.checkState(hazardFile.exists(), "%s doesn't exist", hazardFile.getAbsoluteFile());
			new ZipFile(hazardFile).close();
		}
		
		File indexDir = new File(mainDir, "2022_03_30-sweeps_index");
		Preconditions.checkState(indexDir.exists() || indexDir.mkdir());
		File resourcesDir = new File(indexDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<LogicTreeHazardCompare> fixedCompInstances = new ArrayList<>();
		for (int i=0; i<fixedComps.size(); i++)
			fixedCompInstances.add(loadHaz(fixedComps.get(i), fixedCompNames.get(i)));
		
		LogicTreeHazardCompare prev = null;
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.initNewLine();
		table.addColumn("Name").addColumn("Compared To Previous");
		for (String name : fixedCompNames)
			table.addColumn(name);
		table.finalizeLine();
		
		String targetImageName = "resources/pga_TWO_IN_50_mean_comp_pDiff.png";
		
		for (int i=0; i<runDirs.size(); i++) {
			File runDir = runDirs.get(i);
			String runName = runNames.get(i);
			
			LogicTreeHazardCompare haz = loadHaz(runDir, runName);
			
			table.initNewLine().addColumn("["+runName+"](../"+runDir.getName()+")");
			
			if (prev == null) {
				table.addColumn("_(N/A)_");
			} else {
				File outputDir = new File(runDir, "hazard_maps_comp_prev");
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				
				File srcFile = new File(outputDir, targetImageName);
				
				if (replot || !srcFile.exists())
					haz.buildReport(outputDir, runName, prev, runNames.get(i-1));
				
				File destFile = new File(resourcesDir, i+"_vs_"+(i-1)+".png");
				Files.copy(srcFile, destFile);
				
				table.addColumn("![Comparison]("+resourcesDir.getName()+"/"+destFile.getName()+")");
			}
			
			for (int j=0; j<fixedCompDirNames.size(); j++) {
				File outputDir = new File(runDir, fixedCompDirNames.get(j));
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				
				File srcFile = new File(outputDir, targetImageName);
				
				if (replot || !srcFile.exists())
					haz.buildReport(outputDir, runName, fixedCompInstances.get(j), fixedCompNames.get(j));
				
				File destFile = new File(resourcesDir, i+"_vs_comp_"+j+".png");
				Files.copy(srcFile, destFile);
				
				table.addColumn("![Comparison]("+resourcesDir.getName()+"/"+destFile.getName()+")");
			}
			
			table.finalizeLine();
			
			prev = haz;
		}
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Model Sweeps");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		MarkdownUtils.writeReadmeAndHTML(lines, indexDir);
	}
	
	private static final ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50 };
	private static final double[] periods = { 0d, 1d };
	
	private static final double spacing = 0.1;
	
	private static LogicTreeHazardCompare loadHaz(File dir, String name) throws IOException {
		SolutionLogicTree slt = SolutionLogicTree.load(new File(dir, "results.zip"));
		LogicTree<?> tree = slt.getLogicTree();
		boolean hasFM32 = false;
		for (LogicTreeBranch<?> branch : tree) {
			if (branch.hasValue(FaultModels.FM3_2)) {
				hasFM32 = true;
				break;
			}
		}
		if (hasFM32)
			tree = tree.matchingAll(FaultModels.FM3_1);
		return new LogicTreeHazardCompare(slt, tree, new File(dir, "results_hazard.zip"), rps, periods, spacing);
	}

}
