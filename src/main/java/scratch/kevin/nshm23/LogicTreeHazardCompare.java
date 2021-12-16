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
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
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
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.nshm23.logicTree.SegmentationModel;
import scratch.nshm23.logicTree.SubSectConstraintModel;

public class LogicTreeHazardCompare {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
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
		
//		File mainDir = new File(invDir, "2021_11_30-u3_branches-orig_calcs-5h");
//		String mainName = "UCERF3 As Published";
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = null;
//		String compName = null;
//		LogicTreeNode[] compSubsetNodes = null;
//		File outputDir = new File(mainDir, "hazard_maps");
		
//		File mainDir = new File(invDir, "2021_11_30-nshm23_draft_branches-FM3_1-FaultSpec");
//		String mainName = "NSHM23 Draft";
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = new File(invDir, "2021_11_23-u3_branches-FM3_1-5h");
//		String compName = "UCERF3 Redo";
//		LogicTreeNode[] compSubsetNodes = null;
//		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3_redo");
		
//		File mainDir = new File(invDir, "2021_11_24-nshm23_draft_branches-FM3_1");
//		String mainName = "SectTotRate";
//		LogicTreeNode[] subsetNodes = {SubSectConstraintModel.TOT_NUCL_RATE};
//		File compDir = new File(invDir, "2021_11_24-nshm23_draft_branches-FM3_1");
//		String compName = "SectMFD";
//		LogicTreeNode[] compSubsetNodes = {SubSectConstraintModel.NUCL_MFD};
//		File outputDir = new File(mainDir, "hazard_maps_sect_rate_vs_mfd");
		
//		File mainDir = new File(invDir, "2021_11_24-nshm23_draft_branches-FM3_1");
//		String mainName = "NSHM23 Draft SectTotRate";
//		LogicTreeNode[] subsetNodes = {SubSectConstraintModel.TOT_NUCL_RATE};
//		File outputDir = new File(mainDir, "hazard_maps_sect_rate_vs_ucerf3_redo");
////		String mainName = "NSHM23 Draft SectMFD";
////		LogicTreeNode[] subsetNodes = {SubSectConstraintModel.NUCL_MFD};
////		File outputDir = new File(mainDir, "hazard_maps_sect_mfd_vs_ucerf3_redo");
//		File compDir = new File(invDir, "2021_11_23-u3_branches-FM3_1-5h");
//		String compName = "UCERF3 Redo";
//		LogicTreeNode[] compSubsetNodes = null;
		
//		File mainDir = new File(invDir, "2021_12_01-nshm23_draft_branches-no_paleo-no_parkfield-FM3_1-SysAvg");
//		String mainName = "No Paleo/Park, Single MFD NSHM23 Draft";
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = new File(invDir, "2021_12_01-u3_branches-no_paleo-no_parkfield-single_mfd_reg-FM3_1-5h");
//		String compName = "UCERF3";
//		LogicTreeNode[] compSubsetNodes = null;
//		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3_no_paleo_park");
		
//		File mainDir = new File(invDir, "2021_12_03-nshm23_draft_branches-no_paleo-no_parkfield-FM3_1-FaultSpec");
//		String mainName = "No Paleo/Park, Single MFD NSHM23 Draft";
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = new File(invDir, "2021_12_01-u3_branches-no_paleo-no_parkfield-single_mfd_reg-FM3_1-5h");
//		String compName = "UCERF3";
//		LogicTreeNode[] compSubsetNodes = null;
//		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3_no_paleo_park");
		
//		File mainDir = new File(invDir, "2021_12_12-nshm23_draft_branches-coulomb-FM3_1-ZENGBB-Shaw09Mod-FaultSpec");
//		String mainName = "CoulombRS, Zeng/ShawMod, Segmentation";
//		LogicTreeNode[] subsetNodes = null;
////		File compDir = new File(invDir, "2021_12_11-u3_branches-coulomb-FM3_1-ZENGBB-Shaw09Mod-5h");
////		String compName = "CoulombRS, Zeng/ShawMod, UCERF3";
////		LogicTreeNode[] compSubsetNodes = null;
////		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3");
//		File compDir = mainDir;
//		String compName = mainName.replace("Segmentation", "No Segmentation");
//		LogicTreeNode[] compSubsetNodes = { SegmentationModel.NONE };
//		File outputDir = new File(mainDir, "hazard_maps_vs_no_seg");
		
//		File mainDir = new File(invDir, "2021_12_08-nshm23_draft_branches-FM3_1-TotNuclRate-SubB1-2h");
//		String mainName = "NSHM23 Draft, Sect-Nucl-Rate, Sub-Seis b=1";
//		LogicTreeNode[] subsetNodes = null;
//		File compDir = new File(invDir, "2021_11_23-u3_branches-FM3_1-5h");
//		String compName = "UCERF3";
//		LogicTreeNode[] compSubsetNodes = null;
//		File outputDir = new File(mainDir, "hazard_maps_vs_ucerf3");
		
		File mainDir = new File(invDir, "2021_12_15-nshm23_draft_branches-coulomb-ineq-FM3_1-ZENGBB-Shaw09Mod-TotNuclRate-SubB1");
		String mainName = "NSHM23 Draft With Segmentation";
		LogicTreeNode[] subsetNodes = null;
		File compDir = null;
		String compName = null;
		LogicTreeNode[] compSubsetNodes = null;
		File outputDir = new File(mainDir, "hazard_maps");
		
		SolutionLogicTree solTree = SolutionLogicTree.load(new File(mainDir, "results.zip"));
		
		ReturnPeriods[] rps = ReturnPeriods.values();
		double[] periods = { 0d, 1d };
		double spacing = 0.1;
		
		LogicTreeHazardCompare mapper = new LogicTreeHazardCompare(solTree,
				subsetNodes == null ? solTree.getLogicTree() : solTree.getLogicTree().matchingAll(subsetNodes),
				new File(mainDir, "results_hazard.zip"), rps, periods, spacing);
		
		LogicTreeHazardCompare comp = null;
		if (compDir != null) {
			SolutionLogicTree compSolTree = SolutionLogicTree.load(new File(compDir, "results.zip"));
			comp = new LogicTreeHazardCompare(compSolTree,
					compSubsetNodes == null ? compSolTree.getLogicTree() : compSolTree.getLogicTree().matchingAll(compSubsetNodes),
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
		pDiffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-100d, 100d);
		percentileCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 100d);
		percentileCPT.setNanColor(Color.BLACK);
		
		ZipFile zip = new ZipFile(mapsZipFile);
		
		branches = tree.getBranches();
		weights = new ArrayList<>();
		
		maps = HashBasedTable.create();
		
		GriddedRegion gridReg = null;
		for (int i=0; i<branches.size(); i++) {
			LogicTreeBranch<?> branch = branches.get(i);
			double weight = branch.getBranchWeight();
			weights.add(weight);
			totWeight += weight;
			
			if (gridReg == null) {
				FaultSystemSolution sol = solLogicTree.forBranch(branch);
				Region region = ReportMetadata.detectRegion(sol);
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
	
	private GriddedGeoDataSet calcPercentile(GriddedGeoDataSet[] maps, GriddedGeoDataSet comp) {
		return calcPercentile(List.of(maps), weights, comp);
	}
	
	private GriddedGeoDataSet calcPercentile(List<GriddedGeoDataSet> maps, List<Double> weights, GriddedGeoDataSet comp) {
		Preconditions.checkState(maps.size() == weights.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(maps.get(0).getRegion(), false);
		
		double totWeight;
		if (weights == this.weights) {
			totWeight = this.totWeight;
		} else {
			totWeight = 0d;
			for (double weight : weights)
				totWeight += weight;
		}
		
		for (int i=0; i<ret.size(); i++) {
			double compVal = comp.get(i);
			double weightAbove = 0d; // weight of tree thats above comp value
			double weightEqual = 0d; // weight of tree thats equal to comp value
			double weightBelow = 0d; // weight of tree thats below comp value
			int numAbove = 0;
			int numBelow = 0;
			for (int j=0; j<maps.size(); j++) {
				double val = maps.get(j).get(i);
				double weight = weights.get(j);
				if (((float) compVal == (float)val) || (Double.isNaN(compVal)) && Double.isNaN(val)) {
					weightEqual += weight;
				} else if (compVal < val) {
					weightAbove += weight;
					numAbove++;
				} else if (compVal > val) {
					numBelow++;
					weightBelow += weight;
				}
			}
			if (weightEqual != 0d) {
				// redistribute any exactly equal to either side
				weightAbove += 0.5*weightEqual;
				weightBelow += 0.5*weightEqual;
			}
			// normalize by total weight
			weightAbove /= totWeight;
			weightBelow /= totWeight;
			double percentile;
			if (numAbove == maps.size() || numBelow == maps.size())
				percentile = Double.NaN;
			else
				percentile = 100d*weightBelow;
//			if (numAbove == maps.size() || numBelow == maps.size()) {
//				percentile = Double.NaN;
//			} else if (weightAbove > weightBelow) {
//				// more of the distribution is above my value
//				// this means that the percentile is <50
//				percentile = 100d*weightBelow/totWeight;
//			} else {
//				// more of the distribution is below my value
//				// this means that the percentile is >50
//				percentile = 100d*weightAbove/totWeight;
//				percentile = 100d*(1d - (weightAbove/totWeight));
//			}
			ret.set(i, percentile);
		}
		
		return ret;
	}
	
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
				
				lines.add("## "+label);
				lines.add(topLink);
				lines.add("");
				
				GriddedGeoDataSet[] maps = this.maps.get(rp, period);
				Preconditions.checkNotNull(maps);
				for (int i=0; i<maps.length; i++)
					Preconditions.checkNotNull(maps[i], "map %s is null", i);
				
				GriddedGeoDataSet mean = buildMean(maps);
				GriddedGeoDataSet max = buildMax(maps);
				GriddedGeoDataSet min = buildMin(maps);
				GriddedGeoDataSet spread = buildSpread(log10(min), log10(max));
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine(meanMinMaxSpreadMaps(mean, min, max, spread, name, label, prefix, resourcesDir));
				
				GriddedGeoDataSet cmean = null;
				GriddedGeoDataSet cmin = null;
				GriddedGeoDataSet cmax = null;
				GriddedGeoDataSet cspread = null;
				
				if (comp != null) {
					GriddedGeoDataSet[] cmaps = comp.maps.get(rp, period);
					Preconditions.checkNotNull(cmaps);
					for (int i=0; i<cmaps.length; i++)
						Preconditions.checkNotNull(cmaps[i], "map %s is null", i);
					
					cmean = comp.buildMean(cmaps);
					cmax = comp.buildMax(cmaps);
					cmin = comp.buildMin(cmaps);
					cspread = comp.buildSpread(log10(cmin), log10(cmax));
					table.addLine(meanMinMaxSpreadMaps(cmean, cmin, cmax, cspread, compName, label, prefix+"_comp", resourcesDir));
				}
				
				table.invert();
				
				if (cmean != null) {
					// add min vs min and max vs max comparisons
					table.addLine(MarkdownUtils.boldCentered("Min vs Min Comparison"),
							MarkdownUtils.boldCentered("Max vs Max Comparison"));
					table.initNewLine();
					GriddedGeoDataSet pDiff = buildPDiff(min, cmin);
					File map = mapper.plotMap(resourcesDir, prefix+"_comp_min_pDiff", pDiff, pDiffCPT, name+" vs "+compName,
							"Min Comparison, % Difference, "+label, true);
					table.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
					pDiff = buildPDiff(max, cmax);
					map = mapper.plotMap(resourcesDir, prefix+"_comp_max_pDiff", pDiff, pDiffCPT, name+" vs "+compName,
							"Max Comparison, % Difference, "+label, true);
					table.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
					table.finalizeLine();
				}
				
				lines.addAll(table.build());
				
				if (cmean != null) {
					table = MarkdownUtils.tableBuilder();
					
					table.addLine(MarkdownUtils.boldCentered("Mean-Map Difference"),
							MarkdownUtils.boldCentered("Comparison Diff From Extremes"));
					
					GriddedGeoDataSet pDiff = buildPDiff(mean, cmean);
					table.initNewLine();
					File map = mapper.plotMap(resourcesDir, prefix+"_comp_pDiff", pDiff, pDiffCPT, name+" vs "+compName,
							"Comparison - Primary, % Difference, "+label, true);
					table.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
					GriddedGeoDataSet pDiffFromRange = buildPDiffFromRange(min, max, cmean);
					map = mapper.plotMap(resourcesDir, prefix+"_comp_pDiff_range", pDiffFromRange, pDiffCPT, name+" vs "+compName,
							"Comparison - Extremes, % Difference, "+label, true);
					table.addColumn("![Range Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
					table.finalizeLine();
					
					table.addLine(MarkdownUtils.boldCentered("Comparison Mean Percentile"),
							MarkdownUtils.boldCentered("Spread Difference"));
					GriddedGeoDataSet cMeanPercentile = calcPercentile(maps, cmean);
					table.initNewLine();
					map = mapper.plotMap(resourcesDir, prefix+"_comp_percentile", cMeanPercentile,
							percentileCPT, name+" vs "+compName, "Comparison %-ile, "+label);
					table.addColumn("![Percentile Map]("+resourcesDir.getName()+"/"+map.getName()+")");
					GriddedGeoDataSet spreadDiff = new GriddedGeoDataSet(spread.getRegion(), false);
					for (int i=0; i<spreadDiff.size(); i++)
						spreadDiff.set(i, spread.get(i)-cspread.get(i));
					map = mapper.plotMap(resourcesDir, prefix+"_comp_spread_diff", spreadDiff,
							spreadDiffCPT, name+" vs "+compName, "Log10 Spread Difference, "+label);
					table.addColumn("![Spread Diff]("+resourcesDir.getName()+"/"+map.getName()+")");
					table.finalizeLine();
					
					lines.add("");
					lines.addAll(table.build());
				}
				lines.add("");
				lines.add("### "+label+" Logic Tree Comparisons");
				lines.add(topLink); lines.add("");
				
				// plot mean percentile
				table = MarkdownUtils.tableBuilder();
				
				File meanMap = new File(outputDir, prefix+"_mean.png");
				GriddedGeoDataSet meanPercentile = calcPercentile(maps, mean);
				File meanPercentileMap = mapper.plotMap(resourcesDir, prefix+"_mean_percentile",
						meanPercentile, percentileCPT, "Branch-Averaged Percentiles",
						"Branch Averaged %-ile, "+label);
				table.addLine("Branch Averaged Map", "Branch Averaged Percentiles");
				table.addLine("![BA map]("+resourcesDir.getName()+"/"+meanMap.getName()+")",
						"![BA percentiles]("+resourcesDir.getName()+"/"+meanPercentileMap.getName()+")");
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
									pDiff, pDiffCPT, choice.getShortName()+" Comparison",
									choice.getShortName()+" - Mean, % Difference, "+label, true);
							mapTable.addColumn("![Difference Map]("+resourcesDir.getName()+"/"+map.getName()+")");
							GriddedGeoDataSet percentile = calcPercentile(maps, choiceMap);
							map = mapper.plotMap(resourcesDir, prefix+"_"+choice.getFilePrefix()+"_percentile",
									percentile, percentileCPT, choice.getShortName()+" Comparison",
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
							GriddedGeoDataSet percentileWithout = calcPercentile(mapsWithout, weightsWithout, choiceMap);
							map = mapper.plotMap(resourcesDir, prefix+"_"+choice.getFilePrefix()+"_percentile_without",
									percentileWithout, percentileCPT, choice.getShortName()+" Comparison",
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
	
	private List<String> meanMinMaxSpreadMaps(GriddedGeoDataSet mean, GriddedGeoDataSet min, GriddedGeoDataSet max,
			GriddedGeoDataSet spread, String title, String label, String prefix, File resourcesDir) throws IOException {
		List<String> ret = new ArrayList<>();
		
		ret.add(MarkdownUtils.boldCentered("Weighted-Average"));
		File map = mapper.plotMap(resourcesDir, prefix+"_mean", log10(mean), logCPT, title, "Log10 Weighted-Average, "+label);
		ret.add("![Mean Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		
		ret.add(MarkdownUtils.boldCentered("Spread"));
		map = mapper.plotMap(resourcesDir, prefix+"_spread", spread, spreadCPT, title, "Log10 Spread, "+label);
		ret.add("![Spread Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		
		ret.add(MarkdownUtils.boldCentered("Minimum"));
		map = mapper.plotMap(resourcesDir, prefix+"_min", log10(min), logCPT, title, "Log10 Min, "+label);
		ret.add("![Min Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		
		ret.add(MarkdownUtils.boldCentered("Maximum"));
		map = mapper.plotMap(resourcesDir, prefix+"_max", log10(max), logCPT, title, "Log10 Max, "+label);
		ret.add("![Max Map]("+resourcesDir.getName()+"/"+map.getName()+")");
		
		return ret;
	}

}
