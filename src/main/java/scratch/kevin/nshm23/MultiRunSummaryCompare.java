package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.zip.ZipFile;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;

public class MultiRunSummaryCompare {
	
	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		List<File> runDirs = new ArrayList<>();
		List<String> runNames = new ArrayList<>();
		
		List<File> fixedComps = new ArrayList<>();
		List<String> fixedCompNames = new ArrayList<>();
		List<String> fixedCompDirNames = new ArrayList<>();
		List<PlotCurveCharacterstics> fixedCurveChars = new ArrayList<>();

		PlotCurveCharacterstics primaryChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK);
		PlotCurveCharacterstics prevChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY);
		PlotCurveCharacterstics dataChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.GREEN.darker());
		
		boolean replot = false;
		
//		File indexDir = new File(mainDir, "2022_06_10-sweeps_index");
//		
//		runDirs.add(new File(mainDir, "2022_03_24-u3_branches-FM3_1-2000ip"));
//		runNames.add("UCERF3 New Anneal");
//		fixedComps.add(runDirs.get(runDirs.size()-1));
//		fixedCompNames.add(runNames.get(runNames.size()-1));
//		fixedCompDirNames.add("hazard_maps_comp_u3_new_anneal");
//		fixedCurveChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLUE.darker()));
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_seg-full_sys_inv-FM3_1-U3RedRupSet-DsrUni-TotNuclRate-SubB1-IncludeThruCreep"));
//		runNames.add("Draft w/ U3 reduced RS, full-sys inv, no seg");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_seg-FM3_1-U3RedRupSet-DsrUni-TotNuclRate-SubB1-IncludeThruCreep"));
//		runNames.add("Draft w/ U3 reduced RS, cluster-spec inv, no seg");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_seg-FM3_1-Coulomb5kmRupSet-DsrUni-TotNuclRate-SubB1-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS (5km max jump dist), cluster-spec inv, no seg");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS, no seg");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-IncludeThruCreep-branch-translated-min3km"));
//		runNames.add("Draft w/ Coulomb RS, strict seg");
//		fixedComps.add(runDirs.get(runDirs.size()-1));
//		fixedCompNames.add(runNames.get(runNames.size()-1));
//		fixedCompDirNames.add("hazard_maps_comp_nshm23_draft_strict_min_3km");
//		fixedCurveChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA.darker()));
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-StrictEquivThreshAvg-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS, inequality seg, 2km shift, thresh avg @ strict bins");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-NoAdj-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS, inequality seg, 2km shift, no adjustment");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS, inequality seg, 2km shift, thresh avg");
//		fixedComps.add(runDirs.get(runDirs.size()-1));
//		fixedCompNames.add(runNames.get(runNames.size()-1));
//		fixedCompDirNames.add("hazard_maps_comp_nshm23_draft_thresh_2km");
//		fixedCurveChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.ORANGE.darker()));
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgRelGR-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS, inequality seg, 2km shift, rel-gr thresh avg");
//		
//		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
//		runNames.add("Draft w/ Coulomb RS, inequality seg, 2km shift, iterative rel-gr thresh avg");
//		fixedComps.add(runDirs.get(runDirs.size()-1));
//		fixedCompNames.add(runNames.get(runNames.size()-1));
//		fixedCompDirNames.add("hazard_maps_comp_nshm23_draft_iter_rel_gr_thresh_2km");
//		fixedCurveChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED.darker()));
		
		File indexDir = new File(mainDir, "2022_06_10-sensitivity-tests");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Draft Preferred Model");
		fixedComps.add(runDirs.get(runDirs.size()-1));
		fixedCompNames.add(runNames.get(runNames.size()-1));
		fixedCompDirNames.add("hazard_maps_comp_nshm23_draft_iter_rel_gr_thresh_2km");
		fixedCurveChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED.darker()));
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep-10000ip"));
		runNames.add("5x converged");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-full_sys_inv-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Full-Sys-Inv");
		fixedComps.add(runDirs.get(runDirs.size()-1));
		fixedCompNames.add(runNames.get(runNames.size()-1));
		fixedCompDirNames.add("hazard_maps_comp_full_sys_inv");
		fixedCurveChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE.darker()));
		
		runDirs.add(new File(mainDir, "2022_07_08-nshm23_u3_hybrid_branches-full_sys_inv-no_reweight_use_prev-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Full-Sys-Inv, no-reweight, use prev");
		
		runDirs.add(new File(mainDir, "2022_07_08-nshm23_u3_hybrid_branches-full_sys_inv-no_reweight_use_prev_avg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Full-Sys-Inv, no-reweight, use average");
		
		runDirs.add(new File(mainDir, "2022_07_11-nshm23_u3_hybrid_branches-full_sys_inv-no_reweight_use_orig-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Full-Sys-Inv, no-reweight, use initial");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_reweight_use_orig-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no-reweight, use initial");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_reweight_use_orig-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no-reweight, use initial");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_paleo_parkfield-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no paleo/parkfield");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_scale_adj_mfds-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no scale adj. mfds");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-no_mfd_sigma_data_adj-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no mfd sigma adjust");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-NoShift-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no seg shift");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift1km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("1km seg shift");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift3km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("3km seg shift");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR"));
		runNames.add("with & without creeping section");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-NuclMFD-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("nucleation MFD branch");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SysAvgSubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("sys-avg b=1 sub-seis reduction");
		
		runDirs.add(new File(mainDir, "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("no sub-seis reduction");
		
		runDirs.add(new File(mainDir, "2022_07_12-nshm23_u3_hybrid_branches-classic_sa-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Classic SA Cooling");
		
		runDirs.add(new File(mainDir, "2022_07_12-nshm23_u3_hybrid_branches-exp_perturb-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Exponential Perturb");
		
		runDirs.add(new File(mainDir, "2022_07_12-nshm23_u3_hybrid_branches-limit_zeros-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("Limit Zero Rates");
		
		runDirs.add(new File(mainDir, "2022_07_12-nshm23_u3_hybrid_branches-no_avg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep"));
		runNames.add("No-Averaging");
		
//		runDirs.add(new File(mainDir, ""));
//		runNames.add("");
		
		Preconditions.checkState(indexDir.exists() || indexDir.mkdir());
		
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
			
			table.initNewLine().addColumn("**["+runName+"](../"+runDir.getName()+")**");
			
			if (prev == null) {
				table.addColumn("_(N/A)_");
			} else {
				File outputDir = new File(runDir, "hazard_maps_comp_prev");
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				
				File prevNameFile = new File(outputDir, "prev_dir_name.txt");
				boolean replotPrevious = true;
				if (prevNameFile.exists()) {
					String prevName = Files.readLines(prevNameFile, Charset.defaultCharset()).get(0).trim();
					if (prevName.equals(runDir.getName()))
						replotPrevious = false;
				}
				
				File srcFile = new File(outputDir, targetImageName);
				
				if (replot || !srcFile.exists() || replotPrevious) {
					haz.buildReport(outputDir, runName, prev, runNames.get(i-1));
					// write out the dir name
					FileWriter fw = new FileWriter(prevNameFile);
					fw.write(runDir.getName()+"\n");
					fw.close();
				}
				
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
		int tocIndex = lines.size();
		String topLink = "_[(top)](#table-of-contents)_";
		lines.add("## Hazard");
		lines.add(topLink); lines.add("");
		lines.addAll(table.build());
		lines.add("");

		TableBuilder fullTable = MarkdownUtils.tableBuilder();
		fullTable.addLine("Name", "Inremental", "Cumulative");
		TableBuilder zoomTable = MarkdownUtils.tableBuilder();
		zoomTable.addLine("Name", "Inremental", "Cumulative");

		Range xRange = new Range(5d, 9d);
		Range yRange = new Range(1e-7, 1e1);
		Range xRangeZoom = new Range(6d, 8.5d);
		Range yRangeZoom = new Range(1e-4, 1e0);
		
		List<IncrementalMagFreqDist> fixedMFDs = new ArrayList<>();
		System.out.println("Loading comparison MFDs");
		for (int i=0; i<fixedComps.size(); i++) {
			IncrementalMagFreqDist mfd = loadMFD(loadBA_Sol(fixedComps.get(i)), xRange);
			mfd.setName(fixedCompNames.get(i));
			fixedMFDs.add(mfd);
		}
		
		List<InversionMisfitStats> runMisfitStats = new ArrayList<>();
		List<Integer> runNumNonZeros = new ArrayList<>();
		List<Double> runFractNonZeros = new ArrayList<>();
		
		IncrementalMagFreqDist prevMFD = null;
		for (int i=0; i<runDirs.size(); i++) {
			File runDir = runDirs.get(i);
			String runName = runNames.get(i);
			
			FaultSystemSolution sol = loadBA_Sol(runDir);
			
			runMisfitStats.add(sol.getModule(InversionMisfitStats.class));
			
			IncrementalMagFreqDist mfd = loadMFD(sol, xRange);
			mfd.setName(runName);
			
			int numNonZero = 0;
			for (double rate : sol.getRateForAllRups())
				if ((float)rate > 0f)
					numNonZero++;
			
			runNumNonZeros.add(numNonZero);
			runFractNonZeros.add((double)numNonZero/(double)sol.getRupSet().getNumRuptures());
			
			List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			double roundedMaxMag = mfd.getX(mfd.getClosestXIndex(sol.getRupSet().getMaxMag()));
			GutenbergRichterMagFreqDist[] dataMFDs = {
					new GutenbergRichterMagFreqDist(5.05, 9.95, 50),
					new GutenbergRichterMagFreqDist(5.05, 9.95, 50),
					new GutenbergRichterMagFreqDist(5.05, 9.95, 50)
			};
			dataMFDs[0].setAllButTotMoRate(5.05, roundedMaxMag, TotalMag5Rate.RATE_6p5.getRateMag5(), 1d);
			dataMFDs[1].setAllButTotMoRate(5.05, roundedMaxMag, TotalMag5Rate.RATE_7p9.getRateMag5(), 1d);
			dataMFDs[2].setAllButTotMoRate(5.05, roundedMaxMag, TotalMag5Rate.RATE_9p6.getRateMag5(), 1d);
			
			for (IncrementalMagFreqDist dataMFD : dataMFDs) {
				if (incrFuncs.isEmpty())
					dataMFD.setName("UCERF3 Data Constraints");
				else
					dataMFD.setName(null);
				incrFuncs.add(dataMFD);
				chars.add(dataChar);
			}
			
			incrFuncs.add(mfd);
			chars.add(primaryChar);
			
			if (prevMFD != null && !fixedCompNames.contains(prevMFD.getName())) {
				incrFuncs.add(prevMFD);
				chars.add(prevChar);
			}
			
			for (int j=0; j<fixedComps.size(); j++) {
				IncrementalMagFreqDist fixedMFD = fixedMFDs.get(j);
				if (!fixedMFD.getName().equals(runName)) {
					incrFuncs.add(fixedMFD);
					chars.add(fixedCurveChars.get(j));
				}
			}
			
			if (sol.getRupSet().hasModule(InversionTargetMFDs.class)) {
				InversionTargetMFDs targets = sol.getRupSet().getModule(InversionTargetMFDs.class);
				
				IncrementalMagFreqDist target = targets.getTotalOnFaultSupraSeisMFD();
				if (target != null) {
					target.setName("Target");
					Color color = SolMFDPlot.SUPRA_SEIS_TARGET_COLOR;
					incrFuncs.add(target);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
					
					if (target instanceof UncertainIncrMagFreqDist) {
						UncertainBoundedIncrMagFreqDist sigmaIncrBounds =
								((UncertainIncrMagFreqDist)target).estimateBounds(UncertaintyBoundType.ONE_SIGMA);
						sigmaIncrBounds.setName("± σ");
						
						incrFuncs.add(sigmaIncrBounds);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
								new Color(color.getRed(), color.getGreen(), color.getBlue(), 60)));
					}
				}
			}
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			for (boolean zoom : new boolean[] {false, true}) {
				Range myXRange, myYRange;
				String prefix = "mfds_"+i;
				if (zoom) {
					myXRange = xRangeZoom;
					myYRange = yRangeZoom;
					table = zoomTable;
					prefix += "_zoom";
				} else {
					myXRange = xRange;
					myYRange = yRange;
					table = fullTable;
				}
				
				PlotSpec spec = new PlotSpec(incrFuncs, chars, runName, "Magnitude", "Incremental Rate (/yr)");
				spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				gp.drawGraphPanel(spec, false, true, myXRange, myYRange);
				PlotUtils.writePlots(resourcesDir, prefix, gp, 900, 800, true, true, false);
				
				table.initNewLine().addColumn("**["+runName+"](../"+runDir.getName()+")**");
				
				table.addColumn("![MFD]("+resourcesDir.getName()+"/"+prefix+".png)");
				
				prefix += "_cml";
				
				List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
				for (int j=0; j<incrFuncs.size(); j++) {
					IncrementalMagFreqDist incr = incrFuncs.get(j);
					PlotCurveCharacterstics pChar = chars.get(j);
					//				if (incr instancoe)
					if (pChar.getLineType() == PlotLineType.SHADED_UNCERTAIN) {
						Preconditions.checkState(incr instanceof UncertainBoundedIncrMagFreqDist);
						UncertainBoundedIncrMagFreqDist bounded = (UncertainBoundedIncrMagFreqDist)incr;
						
						EvenlyDiscretizedFunc cumulative = bounded.getCumRateDistWithOffset();
						
						EvenlyDiscretizedFunc upperCumulative = bounded.getUpper().getCumRateDistWithOffset();
						EvenlyDiscretizedFunc lowerCumulative = bounded.getLower().getCumRateDistWithOffset();
						Preconditions.checkState(cumulative.size() == upperCumulative.size());
						for (int k=0; k<cumulative.size(); k++) {
							upperCumulative.set(k, Math.max(cumulative.getY(k), upperCumulative.getY(k)));
							lowerCumulative.set(k, Math.max(0, Math.min(cumulative.getY(k), lowerCumulative.getY(k))));
						}
						
						UncertainArbDiscFunc cmlBounded = new UncertainArbDiscFunc(cumulative, lowerCumulative, upperCumulative);
						cmlBounded.setName("± σ");
						cmlFuncs.add(cmlBounded);
					} else {
						cmlFuncs.add(incr.getCumRateDistWithOffset());
					}
				}
				spec = new PlotSpec(cmlFuncs, chars, runName, "Magnitude", "Cumulative Rate (/yr)");
				spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				gp.drawGraphPanel(spec, false, true, myXRange, myYRange);
				
				PlotUtils.writePlots(resourcesDir, prefix, gp, 900, 800, true, true, false);
				
				table.addColumn("![MFD]("+resourcesDir.getName()+"/"+prefix+".png)");
				
				table.finalizeLine();
			}
			
			prevMFD = mfd;
		}
		lines.add("## MFDs");
		lines.add(topLink); lines.add("");
		lines.addAll(fullTable.build());
		lines.add("");
		
		lines.add("## Zoomed MFDs");
		lines.add(topLink); lines.add("");
		lines.addAll(zoomTable.build());
		lines.add("");
		
		HashSet<String> uncertWtConstraints = new HashSet<>();
		for (InversionMisfitStats stats : runMisfitStats)
			if (stats != null)
				for (MisfitStats constrStats : stats.getStats())
					if (constrStats.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
						uncertWtConstraints.add(constrStats.range.shortName);
		
		if (!uncertWtConstraints.isEmpty()) {
			List<String> names = new ArrayList<>(uncertWtConstraints);
			Collections.sort(names);
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumn("Name");
			for (String name : names)
				table.addColumn(name);
			table.finalizeLine();
			
			DecimalFormat twoDigits = new DecimalFormat("0.00");
			
			for (int r=0; r<runNames.size(); r++) {
				String runName = runNames.get(r);
				
				InversionMisfitStats runStats = runMisfitStats.get(r);
				
				String bld = fixedCompNames.contains(runName) ? "**" : "";
				
				table.initNewLine().addColumn(bld+runName+bld);
				
				for (String constrName : names) {
					MisfitStats stats = null;
					if (runStats != null) {
						for (MisfitStats stat : runStats.getStats()) {
							if (stat.range.shortName.equals(constrName)) {
								stats = stat;
								break;
							}
						}
					}
					if (stats == null)
						table.addColumn("_(N/A)_");
					else
						table.addColumn(bld+twoDigits.format(stats.absMean)+bld);
				}
				table.finalizeLine();
			}
			
			lines.add("## Constraint MADs");
			lines.add(topLink); lines.add("");
			lines.addAll(table.build());
			lines.add("");
		}
		
		DecimalFormat percentDF = new DecimalFormat("0.00%");
		DecimalFormat countDF = new DecimalFormat("#");
		countDF.setGroupingUsed(true);
		countDF.setGroupingSize(3);
		
		table = MarkdownUtils.tableBuilder();
		table.addLine("Name", "Num Non-Zero Ruptures", "Percent Non-Zero Ruptures");
		
		for (int r=0; r<runNames.size(); r++) {
			String runName = runNames.get(r);
			
			String bld = fixedCompNames.contains(runName) ? "**" : "";
			
			table.initNewLine().addColumn(bld+runName+bld);
			table.addColumn(countDF.format(runNumNonZeros.get(r)));
			table.addColumn(percentDF.format(runFractNonZeros.get(r)));
			table.finalizeLine();
		}
		
		lines.add("## Non-Zero Ruptures");
		lines.add(topLink); lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");
		
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
	
	private static FaultSystemSolution loadBA_Sol(File dir) throws IOException {
		File baFile = null;
		for (File file : dir.listFiles()) {
			if (file.getName().endsWith("_branch_averaged.zip")) {
				Preconditions.checkState(baFile == null, "multiple branch averaged files found in %s", dir.getName());
				baFile = file;
			}
		}
		return FaultSystemSolution.load(baFile);
	}
	
	private static IncrementalMagFreqDist loadMFD(FaultSystemSolution sol, Range xRange) {
		IncrementalMagFreqDist defaultMFD = SolMFDPlot.initDefaultMFD(xRange.getLowerBound()+0.01, xRange.getUpperBound()-0.01);
		
		return sol.calcTotalNucleationMFD(defaultMFD.getMinX(), defaultMFD.getMaxX(), defaultMFD.getDelta());
	}

}
