package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

public class SegMFDsComparePageGen {
	
	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");

////		File primaryDir = new File(mainDir, "2022_02_15-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3/"
//////				+ "JumpProb"); String primaryName = "Jump-Prob";
//////				+ "RupProb"); String primaryName = "Rupture-Prob";
//////				+ "CappedRdst"); String primaryName = "Capped Redist";
////				+ "Greedy"); String primaryName = "Greedy";
//////				+ "StictEquivJumpProb"); String primaryName = "Strict-Equiv-Jump-Prob";
////		File primaryDir = new File(mainDir, "2022_02_16-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3/"
//////			+ "JumpProbGt1km"); String primaryName = "Jump-Prob>1km";
////			+ "RupDistrWorstJumpProb"); String primaryName = "Rup-Distr-Worst-Jump";
//		File primaryDir = new File(mainDir, "2022_02_16-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3_Shift1km/"
////				+ "JumpProb"); String primaryName = "Jump-Prob-Shift-1km";
//				+ "CappedRdst"); String primaryName = "Capped-Redist-Shift-1km";
//		FaultSystemSolution primarySol = FaultSystemSolution.load(
//				new File(primaryDir, "solution.zip"));
//		
//		File compDir = new File(mainDir, "2022_02_15-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3/"
//				+ "None"); String compName = "None";
//		File outputDir = new File(primaryDir, "sect_targets_vs_none");
////		File compDir = new File(mainDir, "2022_02_15-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3/"
////				+ "JumpProb"); String compName = "Jump-Prob";
////		File outputDir = new File(primaryDir, "sect_targets_vs_jump_prob");
//////				+ "CappedRdst"); String compName = "Capped Redist";
//////		File outputDir = new File(primaryDir, "sect_targets_vs_capped_redist");
//		FaultSystemSolution compSol = FaultSystemSolution.load(
//				new File(compDir, "solution.zip"));
		
//		File compDir = new File(mainDir, "2022_02_14-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-U3_ZENG-Shaw09Mod-DsrUni-TotNuclRate-SubB1-SupraB0.8-2000ip");
//		String compName = "Strict-Seg";
//		FaultSystemSolution compSol = FaultSystemSolution.load(
//				new File(compDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
//		File outputDir = new File(primaryDir, "sect_targets_vs_strict_seg");

//		File primaryDir = new File(mainDir, "2022_05_27-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg");
//		String primaryName = "Thresh-Avg";
		File primaryDir = new File(mainDir, "2022_06_03-nshm23_u3_hybrid_branches-cluster_specific_inversion-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgRelGR-IncludeThruCreep");
		String primaryName = "Rel-GR-Thresh-Avg";
		FaultSystemSolution primarySol = FaultSystemSolution.load(new File(primaryDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		
//		File compDir = new File(mainDir, "2022_05_25-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-branch-translated-min3km");
//		String compName = "Strict-Seg";
//		File outputDir = new File(primaryDir, "sect_targets_vs_strict_seg");
//		File compDir = new File(mainDir, "2022_06_06-nshm23_u3_hybrid_branches-cluster_specific_inversion-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-CappedRdst-IncludeThruCreep");
//		String compName = "Capped-Redistribution";
//		File outputDir = new File(primaryDir, "sect_targets_vs_capped_redist");
//		File compDir = new File(mainDir, "2022_05_27-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg");
//		String compName = "Thresh-Avg";
//		File outputDir = new File(primaryDir, "sect_targets_vs_thresh_avg");
		File compDir = new File(mainDir, "2022_06_07-nshm23_u3_hybrid_branches-cluster_specific_inversion-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep");
		String compName = "Rel-GR-Thresh-Avg-Iters";
		File outputDir = new File(primaryDir, "sect_targets_vs_rel_gr_thresh_avg_iters");
		FaultSystemSolution compSol = FaultSystemSolution.load(new File(compDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		
		FaultSystemRupSet rupSet = primarySol.getRupSet();
		Preconditions.checkState(rupSet.isEquivalentTo(compSol.getRupSet()));
		
		InversionTargetMFDs primaryTargets = rupSet.requireModule(InversionTargetMFDs.class);
		InversionTargetMFDs compTargets = compSol.getRupSet().requireModule(InversionTargetMFDs.class);
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(rupSet, ReportMetadata.detectRegion(rupSet));
		
		IncrementalMagFreqDist defaultMFD = SolMFDPlot.initDefaultMFD(rupSet.getMinMag(), rupSet.getMaxMag());
		Range mfdXRange = new Range(defaultMFD.getMinX()-0.5*defaultMFD.getDelta(),
				defaultMFD.getMaxX()+0.5*defaultMFD.getDelta());
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Section Target MFDs Comparison: "+primaryName+" vs "+compName);
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "_[(top)](#table-of-contents)_";
		
		lines.add("## Regional Total MFDs");
		lines.add(topLink); lines.add("");
		
		List<Region> regions = new ArrayList<>();
		regions.add(null);
		if (rupSet.hasModule(RegionsOfInterest.class))
			regions.addAll(rupSet.requireModule(RegionsOfInterest.class).getRegions());
		
		for (Region region : regions) {
			String title, prefix;
			IncrementalMagFreqDist primaryTarget = null, compTarget = null;
			if (region == null) {
				title = "Total";
				prefix = "reg_mfds_total";
				primaryTarget = primaryTargets.getTotalOnFaultSupraSeisMFD();
				compTarget = compTargets.getTotalOnFaultSupraSeisMFD();
			} else {
				Preconditions.checkState(region.getName() != null);
				title = region.getName();
				prefix = "reg_mfds_"+region.getName().replaceAll("\\W+", "_");
				
				double[] fractsInside = region == null ?
						null : rupSet.getFractSectsInsideRegion(region, MFDInversionConstraint.MFD_FRACT_IN_REGION_TRACE_ONLY);
				
				List<? extends IncrementalMagFreqDist> primarySupraTargets = primaryTargets.getOnFaultSupraSeisNucleationMFDs();
				List<? extends IncrementalMagFreqDist> compSupraTargets = compTargets.getOnFaultSupraSeisNucleationMFDs();
				
				for (int s=0; s<primarySupraTargets.size(); s++) {
					IncrementalMagFreqDist primaryMFD = primarySupraTargets.get(s);
					IncrementalMagFreqDist compMFD = compSupraTargets.get(s);
					if (primaryTarget == null) {
						primaryTarget = new SummedMagFreqDist(primaryMFD.getMinX(), primaryMFD.size(), primaryMFD.getDelta());
						compTarget = new SummedMagFreqDist(compMFD.getMinX(), compMFD.size(), compMFD.getDelta());
					}
					if (fractsInside != null) {
						if (fractsInside[s] == 0d) {
							continue;
						} else if (fractsInside[s] < 1) {
							primaryMFD = primaryMFD.deepClone();
							primaryMFD.scale(fractsInside[s]);
							compMFD = compMFD.deepClone();
							compMFD.scale(fractsInside[s]);
						}
					}
					((SummedMagFreqDist)primaryTarget).addIncrementalMagFreqDist(primaryMFD);
					((SummedMagFreqDist)compTarget).addIncrementalMagFreqDist(compMFD);
				}
			}
			
			IncrementalMagFreqDist primarySolMFD = primarySol.calcNucleationMFD_forRegion(
					region, primaryTarget.getMinX(), primaryTarget.getMaxX(), primaryTarget.getDelta(),
					MFDInversionConstraint.MFD_FRACT_IN_REGION_TRACE_ONLY);
			IncrementalMagFreqDist compSolMFD = compSol.calcNucleationMFD_forRegion(
					region, compTarget.getMinX(), compTarget.getMaxX(), compTarget.getDelta(),
					MFDInversionConstraint.MFD_FRACT_IN_REGION_TRACE_ONLY);
			
			primaryTarget.setName(primaryName+" Target");
			primarySolMFD.setName(primaryName+" Solution");
			compTarget.setName(compName+" Target");
			compSolMFD.setName(compName+" Solution");
			
			writeMFDPlot(resourcesDir, prefix, title, primaryTarget, primarySolMFD, compTarget, compSolMFD, mfdXRange);
			
			if (regions.size() > 1) {
				lines.add("**"+title+"**");
				lines.add("");
			}
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			table.addLine("![Incremental]("+resourcesDir.getName()+"/"+prefix+".png)",
					"![Incremental]("+resourcesDir.getName()+"/"+prefix+"_cumulative.png)");
			
			lines.addAll(table.build());
		}
		
		double[] primaryTargetRates = calcTargetNuclRates(primaryTargets);
		double[] primarySolRates = primarySol.calcNucleationRateForAllSects(0d, Double.POSITIVE_INFINITY);
		double[] compTargetRates = calcTargetNuclRates(compTargets);
		double[] compSolRates = compSol.calcNucleationRateForAllSects(0d, Double.POSITIVE_INFINITY);
		
		double maxRate = 0d;
		double minRate = Double.POSITIVE_INFINITY;
		double[][] allRates = {primaryTargetRates, primarySolRates, compTargetRates, compSolRates};
		for (double[] rates : allRates) {
			for (double rate : rates) {
				maxRate = Math.max(maxRate, rate);
				if (rate > 0d)
					minRate = Math.min(minRate, rate);
			}
		}
		maxRate = Math.pow(10, Math.ceil(Math.log10(maxRate)));
		minRate = Math.pow(10, Math.floor(Math.log10(minRate)));
		CPT nuclRateCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(Math.log10(minRate), Math.log10(maxRate));
		CPT pDiffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-100d, 100d);
		
		String prefix;
		
		lines.add("## Section Nucleation Rates");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		// will invert this table
		for (boolean primary : new boolean[] {true, false}) {
			String name = primary ? primaryName : compName;
			double[] targets = primary ? primaryTargetRates : compTargetRates;
			double[] sols = primary ? primarySolRates : compSolRates;
			table.initNewLine();
			table.addColumn(MarkdownUtils.boldCentered(name+" Target Nucleation Rates"));
			mapMaker.plotSectScalars(log10(targets), nuclRateCPT, "Log10 Target Nucleation Rates");
			prefix = (primary ? "primary" : "comp")+"_target_nucl_rates";
			mapMaker.plot(resourcesDir, prefix, name);
			table.addColumn("![Target Nucleation Rates]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			table.addColumn(MarkdownUtils.boldCentered(name+" Solution Nucleation Rates"));
			mapMaker.plotSectScalars(log10(sols), nuclRateCPT, "Log10 Solution Nucleation Rates");
			prefix = (primary ? "primary" : "comp")+"_sol_nucl_rates";
			mapMaker.plot(resourcesDir, prefix, name);
			table.addColumn("![Solution Nucleation Rates]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			table.addColumn(MarkdownUtils.boldCentered(name+" Solution - Target, % Diffs"));
			mapMaker.plotSectScalars(pDiff(sols, targets), pDiffCPT, "Solution - Target, % Difference");
			prefix = (primary ? "primary" : "comp")+"_sol_target_diff_nucl_rates";
			mapMaker.plot(resourcesDir, prefix, name);
			table.addColumn("![Difference]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			File scatterPlot = logScatter(resourcesDir, (primary ? "primary" : "comp")+"_sol_target_nucl_rates_scatter",
					name+" Nucleation Rate Scatter", "Target Nucleation Rate", "Solution Nucleation Rate", targets, sols);
			table.addColumn("![Scatter]("+resourcesDir.getName()+"/"+scatterPlot.getName()+")");
			
			table.finalizeLine();
		}
		table.invert();
		lines.addAll(table.build());
		
		lines.add("## Section Nucleation Comparison");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		
		table.addLine(primaryName+" vs "+compName+", Targets", primaryName+" vs "+compName+", Solution");
		
		table.initNewLine();
		
		double[] targetPDiffs = pDiff(primaryTargetRates, compTargetRates);
		double[] solPDiffs = pDiff(primarySolRates, compSolRates);
		mapMaker.plotSectScalars(targetPDiffs, pDiffCPT, primaryName+" - "+compName+", % Difference");
		prefix = "target_diff_nucl_rates";
		mapMaker.plot(resourcesDir, prefix, "Target Rate Comparison");
		table.addColumn("![Difference]("+resourcesDir.getName()+"/"+prefix+".png)");
		mapMaker.plotSectScalars(solPDiffs, pDiffCPT, primaryName+" - "+compName+", % Difference");
		prefix = "sol_diff_nucl_rates";
		mapMaker.plot(resourcesDir, prefix, "Solution Rate Comparison");
		table.addColumn("![Difference]("+resourcesDir.getName()+"/"+prefix+".png)");
		
		table.finalizeLine();
		
		table.initNewLine();
		File scatterPlot = logScatter(resourcesDir, "target_nucl_rates_scatter", "Target Nucleation Rate Scatter",
				primaryName, compName, primaryTargetRates, compTargetRates);
		table.addColumn("![Scatter]("+resourcesDir.getName()+"/"+scatterPlot.getName()+")");
		scatterPlot = logScatter(resourcesDir, "sol_nucl_rates_scatter", "Solution Nucleation Rate Scatter",
				primaryName, compName, primaryTargetRates, compTargetRates);
		table.addColumn("![Scatter]("+resourcesDir.getName()+"/"+scatterPlot.getName()+")");
		table.finalizeLine();
		
		lines.addAll(table.build());
		
		lines.add("## Biggest Section Differences");
		lines.add(topLink); lines.add("");
		
		Map<Integer, Double> sectAbsDiffsMap = new HashMap<>();
		for (int s=0; s<rupSet.getNumSections(); s++) {
			double maxDiff = Math.max(Math.abs(targetPDiffs[s]), Math.abs(solPDiffs[s]));
			sectAbsDiffsMap.put(s, maxDiff);
		}
		List<Integer> sorted = ComparablePairing.getSortedData(sectAbsDiffsMap);
		// sort from biggest to smallest
		Collections.reverse(sorted);
		int maxNum = 20;
		boolean onePerParent = true;
		HashSet<Integer> prevParents = new HashSet<>();
		
		if (onePerParent) {
			lines.add("_Note: only one section per parent is included below._");
			lines.add("");
		}
		if (primarySol.hasModule(RupMFDsModule.class) || compSol.hasModule(RupMFDsModule.class)) {
			lines.add("_Note: rupture MFDs are present, which means the one or both of the solutions are branch-averaged. "
					+ "In that case, columns below for 'Single Fault' and 'Controlling Jump Distance(s)' only apply to "
					+ "the branch-averaged mean magnitude._");
			lines.add("");
		}
		
		if (!rupSet.hasModule(ClusterRuptures.class))
			rupSet.addModule(ClusterRuptures.singleStranged(rupSet));
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		
		SolutionSlipRates primarySolSlips = primarySol.getModule(SolutionSlipRates.class);
		if (primarySolSlips == null)
			primarySolSlips = SolutionSlipRates.calc(primarySol,
					primarySol.getRupSet().requireModule(AveSlipModule.class),
					primarySol.getRupSet().requireModule(SlipAlongRuptureModel.class));
		SolutionSlipRates compSolSlips = compSol.getModule(SolutionSlipRates.class);
		if (compSolSlips == null)
			compSolSlips = SolutionSlipRates.calc(compSol,
					compSol.getRupSet().requireModule(AveSlipModule.class),
					compSol.getRupSet().requireModule(SlipAlongRuptureModel.class));
		
		int num = 0;
		for (int sectIndex : sorted) {
			FaultSection sect = rupSet.getFaultSectionData(sectIndex);
			if (onePerParent) {
				int parentID = sect.getParentSectionId();
				if (prevParents.contains(parentID))
					continue;
				prevParents.add(parentID);
			}
			
			lines.add("### "+sectIndex+". "+sect.getSectionName());
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("", primaryName, compName, "% Difference");
			table.initNewLine();
			table.addColumn("Target Rate");
			table.addColumn((float)primaryTargetRates[sectIndex]);
			table.addColumn((float)compTargetRates[sectIndex]);
			table.addColumn((float)targetPDiffs[sectIndex]+" %");
			table.finalizeLine();
			table.initNewLine();
			table.addColumn("Solution Rate");
			table.addColumn((float)primarySolRates[sectIndex]);
			table.addColumn((float)compSolRates[sectIndex]);
			table.addColumn((float)solPDiffs[sectIndex]+" %");
			table.finalizeLine();
			double primaryTargetSlipRate = primarySol.getRupSet().getSlipRateForSection(sectIndex);
			double compTargetSlipRate = primarySol.getRupSet().getSlipRateForSection(sectIndex);
			double primarySolSlipRate = primarySolSlips.get(tocIndex);
			double compSolSlipRate = compSolSlips.get(tocIndex);
			table.initNewLine();
			table.addColumn("Target Slip Rate");
			table.addColumn((float)primaryTargetSlipRate);
			table.addColumn((float)compTargetSlipRate);
			table.addColumn((float)(100d*(primaryTargetSlipRate-compTargetSlipRate)/compTargetSlipRate)+" %");
			table.finalizeLine();
			table.initNewLine();
			table.addColumn("Solution Slip Rate");
			table.addColumn((float)primarySolSlipRate);
			table.addColumn((float)compSolSlipRate);
			table.addColumn((float)(100d*(primarySolSlipRate-compSolSlipRate)/compSolSlipRate)+" %");
			table.finalizeLine();
			table.initNewLine();
			table.addColumn("Solution Slip Misfits");
			table.addColumn((float)(100d*(primarySolSlipRate-primaryTargetSlipRate)/primaryTargetSlipRate)+" %");
			table.addColumn((float)(100d*(compSolSlipRate-compTargetSlipRate)/compTargetSlipRate)+" %");
			table.addColumn("");
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			IncrementalMagFreqDist primaryTargetMFD = primaryTargets.getOnFaultSupraSeisNucleationMFDs().get(sectIndex);
			IncrementalMagFreqDist primarySolMFD = primarySol.calcNucleationMFD_forSect(sectIndex,
					primaryTargetMFD.getMinX(), primaryTargetMFD.getMaxX(), primaryTargetMFD.size());
			IncrementalMagFreqDist compTargetMFD = compTargets.getOnFaultSupraSeisNucleationMFDs().get(sectIndex);
			IncrementalMagFreqDist compSolMFD = compSol.calcNucleationMFD_forSect(sectIndex,
					compTargetMFD.getMinX(), compTargetMFD.getMaxX(), compTargetMFD.size());
			
			primaryTargetMFD.setName(primaryName+" Target");
			primarySolMFD.setName(primaryName+" Solution");
			compTargetMFD.setName(compName+" Target");
			compSolMFD.setName(compName+" Solution");
			
			String title = sect.getSectionName();
			prefix = sect.getSectionName().replaceAll("\\W+", "_");
			
			int minMagIndex = Integer.MAX_VALUE;
			int maxMagIndex = 0;
			
			double roundMin = Double.POSITIVE_INFINITY, roundMax = Double.NEGATIVE_INFINITY;
			for (IncrementalMagFreqDist mfd : new IncrementalMagFreqDist[] {
					primaryTargetMFD, primarySolMFD, compTargetMFD, compSolMFD}) {
				for (int i=0; i<mfd.size(); i++) {
					if (mfd.getY(i) > 0) {
						roundMin = Math.min(mfd.getX(i), roundMin);
						roundMax = Math.max(mfd.getX(i), roundMax);
						minMagIndex = Integer.min(minMagIndex, i);
						maxMagIndex = Integer.max(maxMagIndex, i);
					}
				}
			}
			
			roundMin = 0.5*Math.floor(2d*roundMin);
			roundMax = 0.5*Math.ceil(2d*roundMax);
			Range myXRange = new Range(roundMin, roundMax);
			
			writeMFDPlot(resourcesDir, prefix, title, primaryTargetMFD, primarySolMFD, compTargetMFD, compSolMFD, myXRange);
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("![Incremental]("+resourcesDir.getName()+"/"+prefix+".png)",
					"![Incremental]("+resourcesDir.getName()+"/"+prefix+"_cumulative.png)");
			lines.addAll(table.build());
			lines.add("");
			
			List<HashSet<Jump>> magBinJumps = new ArrayList<>();
			for (int i=0; i<primaryTargetMFD.size(); i++)
				magBinJumps.add(new HashSet<>());
			boolean[] singleFaults = new boolean[primaryTargetMFD.size()];
			EvenlyDiscretizedFunc binMinJumps = new EvenlyDiscretizedFunc(primaryTargetMFD.getMinX(),
					primaryTargetMFD.size(), primaryTargetMFD.getDelta());
			for (int i=0; i<binMinJumps.size(); i++)
				binMinJumps.set(i, Double.POSITIVE_INFINITY);
			for (int rupIndex : rupSet.getRupturesForSection(sectIndex)) {
				int magIndex = primaryTargetMFD.getClosestXIndex(rupSet.getMagForRup(rupIndex));
				Jump controlling = null;
				for (Jump jump : cRups.get(rupIndex).getJumpsIterable()) {
					if (jump.fromSection.getSectionId() > jump.toSection.getSectionId())
						jump = jump.reverse();
					if (controlling == null || controlling.distance < jump.distance
							|| (controlling.distance == jump.distance && Jump.id_comparator.compare(jump, controlling) < 0))
						controlling = jump;
				}
				if (controlling == null) {
					binMinJumps.set(magIndex, 0d);
					singleFaults[magIndex] = true;
				} else {
					binMinJumps.set(magIndex, Math.min(binMinJumps.getY(magIndex), controlling.distance));
					HashSet<Jump> curs = magBinJumps.get(magIndex);
					
					if (curs.contains(controlling)) {
						// keep the shorter version
						Jump prev = null;
						for (Jump jump : curs) {
							if (jump.equals(controlling)) {
								prev = jump;
								break;
							}
						}
						if (prev.distance > controlling.distance) {
							// replace it
							curs.remove(prev);
							curs.add(controlling);
						}
					} else {
						curs.add(controlling);
					}
				}
			}
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("Magnitude", primaryName+" Target", "Solution", compName+" Target", "Solution",
					"Single Fault?", "Controlling Jump Distance(s)");
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			PlotCurveCharacterstics distChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK);
			DefaultXY_DataSet curLine = new DefaultXY_DataSet();
			double binDelta = 0.5*binMinJumps.getDelta();
			for (int i=minMagIndex; i<=maxMagIndex; i++) {
				double x = binMinJumps.getX(i);
				double y = binMinJumps.getY(i);
				if (!Double.isFinite(y)) {
					if (curLine.size() > 0) {
						curLine.set(x-binDelta, 0);
						funcs.add(curLine);
						chars.add(distChar);
						curLine = new DefaultXY_DataSet();
					}
				} else {
					if (curLine.size() == 0)
						curLine.set(x-binDelta, 0d);
					curLine.set(x-binDelta, y);
					curLine.set(x+binDelta, y);
				}
			}
			if (curLine.size() > 0) {
				curLine.set(curLine.getMaxX(), 0);
				funcs.add(curLine);
				chars.add(distChar);
			}
			
			PlotSpec magDistSpec = new PlotSpec(funcs, chars, sect.getName(),
					"Magnitude", "Smallest Controlling Jump Distance (km)");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			// convert empty bins to NaNs
			for (int i=0; i<binMinJumps.size(); i++)
				if (Double.isInfinite(binMinJumps.getY(i)))
					binMinJumps.set(i, Double.NaN);
			
			Range xRange = new Range(primaryTargetMFD.getX(minMagIndex)-0.5*primaryTargetMFD.getDelta(),
					primaryTargetMFD.getX(maxMagIndex)+0.5*primaryTargetMFD.getDelta());
			gp.drawGraphPanel(magDistSpec, false, false, xRange, null);
			
			prefix += "_mag_jumps";
			PlotUtils.writePlots(resourcesDir, prefix, gp, 900, 650, true, false, false);
			
			lines.add("![Distance Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			for (int i=minMagIndex; i<=maxMagIndex; i++) {
				table.initNewLine();
				
				if (i < primaryTargetMFD.size()) {
					table.addColumn((float)primaryTargetMFD.getX(i));
					table.addColumn((float)primaryTargetMFD.getY(i));
					table.addColumn((float)primarySolMFD.getY(i));
				} else {
					table.addColumn((float)compTargetMFD.getX(i));
					table.addColumn("_N/A_");
					table.addColumn("_N/A_");
				}
				if (i < compTargetMFD.size()) {
					table.addColumn((float)compTargetMFD.getY(i));
					table.addColumn((float)compSolMFD.getY(i));
				} else {
					table.addColumn("_N/A_");
					table.addColumn("_N/A_");
				}
				if (i >= primaryTargetMFD.size()) {
					table.addColumn("_N/A_");
					table.addColumn("_N/A_");
				} else {
					if (singleFaults[i])
						table.addColumn("**Yes**");
					else
						table.addColumn("No");
					List<Jump> jumps = new ArrayList<>(magBinJumps.get(i));
					if (jumps.isEmpty()) {
						table.addColumn("_(none)_");
					} else {
						jumps.sort(Jump.dist_comparator);
						String jumpStr = null;
						for (Jump jump : jumps) {
							if (jumpStr == null)
								jumpStr = "";
							else
								jumpStr += ", ";
							jumpStr += twoDigits.format(jump.distance);
						}
						jumpStr += " km";
						table.addColumn(jumpStr);
					}
					
				}
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
			
			num++;
			if (num == maxNum)
				break;
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final Color primaryColor = Color.RED;
	private static final Color compColor = Color.BLUE;
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	
	private static double[] calcTargetNuclRates(InversionTargetMFDs mfds) {
		List<? extends IncrementalMagFreqDist> supras = mfds.getOnFaultSupraSeisNucleationMFDs();
		
		double[] ret = new double[supras.size()];
		for (int s=0; s<ret.length; s++)
			ret[s] = supras.get(s).calcSumOfY_Vals();
		
		return ret;
	}
	
	private static double[] log10(double[] vals) {
		double[] ret = new double[vals.length];
		for (int i=0; i<ret.length; i++)
			ret[i] = Math.log10(vals[i]);
		return ret;
	}
	
	private static double[] pDiff(double[] primary, double[] comparison) {
		double[] ret = new double[primary.length];
		for (int i=0; i<ret.length; i++) {
			double z1 = primary[i];
			double z2 = comparison[i];
			double val;
			if (z1 == 0d && z2 == 0d)
				val = 0d;
			else if (z2 == 0d)
				val = Double.POSITIVE_INFINITY;
			else
				val = 100d*(z1-z2)/z2;
			ret[i] = val;
		}
		return ret;
	}
	
	private static void writeMFDPlot(File resourcesDir, String prefix, String title,
			IncrementalMagFreqDist primaryTargetMFD, IncrementalMagFreqDist primarySolMFD,
			IncrementalMagFreqDist compTargetMFD, IncrementalMagFreqDist compSolMFD,
			Range mfdXRange) throws IOException {
		List<IncrementalMagFreqDist> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		funcs.add(primaryTargetMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 4f, primaryColor));
		funcs.add(primarySolMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, primaryColor));
		funcs.add(compTargetMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 4f, compColor));
		funcs.add(compSolMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, compColor));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Magnitude", "Incremental Nucleation Rate (1/yr)");
		spec.setLegendInset(true);
		
		List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
		for (IncrementalMagFreqDist mfd : funcs)
			cmlFuncs.add(mfd.getCumRateDistWithOffset());
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0;
		for (DiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 1e-10 && mfdXRange.contains(pt.getX())) {
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		for (DiscretizedFunc func : cmlFuncs) {
			for (Point2D pt : func) {
				if (pt.getY() > 1e-10 && mfdXRange.contains(pt.getX())) {
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minY))), Math.pow(10, Math.ceil(Math.log10(maxY))));
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, true, mfdXRange, yRange);
		
		PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 700, true, true, false);
		
		
		spec = new PlotSpec(cmlFuncs, chars, title, "Magnitude", "Cumulative Nucleation Rate (1/yr)");
		spec.setLegendInset(true);
		
		gp.drawGraphPanel(spec, false, true, mfdXRange, yRange);
		
		PlotUtils.writePlots(resourcesDir, prefix+"_cumulative", gp, 800, 700, true, true, false);
	}
	
	private static File logScatter(File resourcesDir, String prefix, String title, String xAxisName, String yAxisName,
			double[] xData, double[] yData) throws IOException {
		DefaultXY_DataSet scatter = new DefaultXY_DataSet(xData, yData);
		
		double minNonZero = Double.POSITIVE_INFINITY;
		for (Point2D pt : scatter) {
			if (pt.getX() > 0)
				minNonZero = Math.min(minNonZero, pt.getX());
			if (pt.getY() > 0)
				minNonZero = Math.min(minNonZero, pt.getY());
		}
		
		double max = Math.max(scatter.getMaxX(), scatter.getMaxY());
		
		Range range = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))),
				Math.pow(10, Math.ceil(Math.log10(max))));
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisName, yAxisName);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, true, true, range, range);
		
		PlotUtils.writePlots(resourcesDir, prefix, gp, 800, false, true, false, false);
		
		return new File(resourcesDir, prefix+".png");
	}

}
