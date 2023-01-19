package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels.MinisectionSlipRecord;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class OutlierReplacementPageGen {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/markdown/nshm23-misc");
		double[] ycs = { 2d, 3.5d, 5d };
		
		NSHM23_DeformationModels.OUTLIER_SUB_USE_BOUND = false;
		NSHM23_DeformationModels.OUTLIER_SUB_LOG = false;
		File reportDir = new File(mainDir, "outlier_detection_linear");
		
//		NSHM23_DeformationModels.OUTLIER_SUB_USE_BOUND = false;
//		NSHM23_DeformationModels.OUTLIER_SUB_LOG = true;
//		File reportDir = new File(mainDir, "outlier_detection_log");
		
		Preconditions.checkState(reportDir.exists() || reportDir.mkdir());
		File resourcesDir = new File(reportDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Deformation Model Outlier Replacement "
				+(NSHM23_DeformationModels.OUTLIER_SUB_LOG ? "(log-domain)" : "(linear)"));
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		
		NSHM23_DeformationModels[] defModels = {
				NSHM23_DeformationModels.GEOLOGIC,
				NSHM23_DeformationModels.EVANS,
				NSHM23_DeformationModels.POLLITZ,
				NSHM23_DeformationModels.SHEN_BIRD,
				NSHM23_DeformationModels.ZENG
		};
		
		List<Map<Integer, List<MinisectionSlipRecord>>> dmOrigMiniMaps = new ArrayList<>();
		
		for (NSHM23_DeformationModels dm : defModels)
			dmOrigMiniMaps.add(dm.getMinisections(fm));
		
		Region mapReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		NSHM23_DeformationModels.ORIGINAL_WEIGHTS = false;
		Map<Integer, List<MinisectionSlipRecord>> reweightedAvgMiniMap = NSHM23_DeformationModels.AVERAGE.getMinisections(fm);
		NSHM23_DeformationModels.ORIGINAL_WEIGHTS = true;
		Map<Integer, List<MinisectionSlipRecord>> origAvgMiniMap = NSHM23_DeformationModels.AVERAGE.getMinisections(fm);
		Map<Integer, List<MinisectionSlipRecord>> medianMiniMap = NSHM23_DeformationModels.MEDIAN.getMinisections(fm);
		
		Map<Integer, FaultSection> geolSects = fm.getFaultSectionIDMap();
		
		List<String> stateNames = new ArrayList<>();
		List<Region> stateRegions = new ArrayList<>();
		stateNames.add("Full Model");
		stateRegions.add(null);
		for (NSHM23_SingleStates state : NSHM23_SingleStates.values()) {
			stateNames.add(state.getName());
			stateRegions.add(state.loadRegion());
		}
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		
		CPT ratioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1, 1d);
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().rescale(-5d, 5d);
		
		for (double yc : ycs) {
			lines.add("## Critical Value: Yc="+oDF.format(yc));
			lines.add(topLink); lines.add("");
			
			String ycPrefix = "yc_"+oDF.format(yc).replace('.', 'p');
			String ycLabel = "Yc="+oDF.format(yc);
			
			List<File> mapPlots = new ArrayList<>();
			// these will be overall, above, below
			List<int[]> counts = new ArrayList<>();
			List<double[]> fracts = new ArrayList<>();
			List<double[]> momFracts = new ArrayList<>();
			List<Double> momChanges = new ArrayList<>();
			List<File> scatterToOrigPlots = new ArrayList<>();
			List<File> scatterToMedianPlots = new ArrayList<>();
			List<Integer> origNumZeros = new ArrayList<>();
			List<Integer> filteredNumZeros = new ArrayList<>();
			List<Double> origFractZeros = new ArrayList<>();
			List<Double> filteredFractZeros = new ArrayList<>();
			
			List<double[]> stateOrigMomentsList = new ArrayList<>();
			List<double[]> stateFilteredMomentsList = new ArrayList<>();
			
			NSHM23_DeformationModels.ORIGINAL_WEIGHTS = true;
			List<Map<Integer, List<MinisectionSlipRecord>>> dmFilteredList = new ArrayList<>();
			List<Double> dmWeightsList = new ArrayList<>();
			
			for (int d=0; d<defModels.length; d++) {
				Map<Integer, List<MinisectionSlipRecord>> origMinis = dmOrigMiniMaps.get(d);
				List<Integer> toRemove = new ArrayList<>();
				for (Integer id : origMinis.keySet())
					if (!geolSects.containsKey(id))
						toRemove.add(id);
				for (Integer id : toRemove)
					origMinis.remove(id);
				Map<Integer, List<MinisectionSlipRecord>> filtered = NSHM23_DeformationModels.applyOutlierSubstitution(
						fm, origMinis, yc, NSHM23_DeformationModels.OUTLIER_SUB_LOG);
				
				dmFilteredList.add(filtered);
				dmWeightsList.add(defModels[d].getNodeWeight(null));
				
				int totNum = 0;
				int numFiltered = 0;
				int numAbove = 0;
				int numBelow = 0;
				double totOrigMoment = 0d;
				double totFinalMoment = 0d;
				double momentFiltered = 0d;
				double momentFilteredAbove = 0d;
				double momentFilteredBelow = 0d;
				int origZeros = 0;
				int filteredZeros = 0;
				
				double[] stateOrigMoments = new double[stateRegions.size()];
				stateOrigMomentsList.add(stateOrigMoments);
				double[] stateFilteredMoments = new double[stateRegions.size()];
				stateFilteredMomentsList.add(stateFilteredMoments);
				
				List<FaultSection> origSects = new ArrayList<>();
				List<FaultSection> filteredSects = new ArrayList<>();

				List<Color> filterColors = new ArrayList<>();
				List<Double> filterColorComps = new ArrayList<>();
				
				for (Integer id : origMinis.keySet()) {
					List<MinisectionSlipRecord> origRecs = origMinis.get(id);
					List<MinisectionSlipRecord> filteredRecs = filtered.get(id);
					
					FaultSection fullSect = geolSects.get(id);
					if (fullSect == null)
						// extra fault not in final model
						continue;
					double ddw = fullSect.getOrigDownDipWidth();
					
					for (int i=0; i<origRecs.size(); i++) {
						MinisectionSlipRecord origRec = origRecs.get(i);
						MinisectionSlipRecord filteredRec = filteredRecs.get(i);
						
						double len = LocationUtils.horzDistanceFast(origRec.startLoc, origRec.endLoc);
						
						double area = len*ddw;
						double areaM = area*1e6;
						
						double myOrigMoment = FaultMomentCalc.getMoment(areaM, origRec.slipRate);
						double myFilteredMoment = FaultMomentCalc.getMoment(areaM, filteredRec.slipRate);
						totOrigMoment += myOrigMoment;
						totFinalMoment += myFilteredMoment;
						
						origSects.add(miniToSect(origRec, fullSect, origSects.size(), i));
						filteredSects.add(miniToSect(filteredRec, fullSect, filteredSects.size(), i));
						
						for (int r=0; r<stateRegions.size(); r++) {
							Region reg = stateRegions.get(r);
							
							if (reg == null || reg.contains(origRec.startLoc) || reg.contains(origRec.endLoc)) {
								stateOrigMoments[r] += myOrigMoment;
								stateFilteredMoments[r] += myFilteredMoment;
							}
						}
						
						if (origRec.slipRate == 0d)
							origZeros++;
						if (filteredRec.slipRate == 0d)
							filteredZeros++;
						
						totNum++;
						if ((float)origRec.slipRate != (float)filteredRec.slipRate) {
							numFiltered++;
							momentFiltered += myFilteredMoment;
							if (origRec.slipRate > filteredRec.slipRate) {
								numAbove++;
								momentFilteredAbove += myFilteredMoment;
								filterColors.add(Color.RED.darker());
								filterColorComps.add(1d);
							} else {
								numBelow++;
								momentFilteredBelow += myFilteredMoment;
								filterColors.add(Color.BLUE.darker());
								filterColorComps.add(1d);
							}
						} else {
							filterColors.add(Color.GRAY);
							filterColorComps.add(0d);
						}
					}
				}
				
				RupSetMapMaker mapMaker = new RupSetMapMaker(origSects, mapReg);
				
				mapMaker.plotSectColors(filterColors, null, null, filterColorComps);
				
				String dmName = defModels[d].getShortName();
				String prefix = ycPrefix+"_"+defModels[d].name();
				String mapPrefix = prefix+"_replacement_map";
				
				mapMaker.plot(resourcesDir, mapPrefix, dmName+" Replacement Map, "+ycLabel);
				
				mapPlots.add(new File(resourcesDir, mapPrefix+".png"));
				counts.add(new int[] {numFiltered, numAbove, numBelow});
				fracts.add(new double[] {(double)numFiltered/(double)totNum,
						(double)numAbove/(double)totNum,
						(double)numBelow/(double)totNum});
				momFracts.add(new double[] {(double)momentFiltered/(double)totFinalMoment,
						(double)momentFilteredAbove/(double)totFinalMoment,
						(double)momentFilteredBelow/(double)totFinalMoment});
				momChanges.add((totFinalMoment-totOrigMoment)/totOrigMoment);
				scatterToOrigPlots.add(scatterPlot(filtered, "Filtered, "+ycLabel, origMinis, "Original",
						resourcesDir, prefix+"_orig_scatter", Color.BLUE.darker(), dmName+", "+ycLabel+" Compared to Original"));
				scatterToMedianPlots.add(scatterPlot(filtered, "Filtered, "+ycLabel, medianMiniMap, "Median",
						resourcesDir, prefix+"_median_scatter", Color.GREEN.darker(), dmName+", "+ycLabel+" Compared to Median"));
				
				origNumZeros.add(origZeros);
				filteredNumZeros.add(filteredZeros);
				origFractZeros.add((double)origZeros/(double)totNum);
				filteredFractZeros.add((double)filteredZeros/(double)totNum);
			}
			
			// summary table
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			table.initNewLine().addColumn("");
			for (NSHM23_DeformationModels dm : defModels)
				table.addColumn(dm.getShortName());
			table.finalizeLine();
			
			table.initNewLine().addColumn("Total Num Filtered");
			for (int i=0; i<counts.size(); i++)
				table.addColumn(counts.get(i)[0]+" ("+pDF.format(fracts.get(i)[0])+")");
			table.finalizeLine();
			
			table.initNewLine().addColumn("Above "+ycLabel);
			for (int i=0; i<counts.size(); i++)
				table.addColumn(counts.get(i)[1]+" ("+pDF.format(fracts.get(i)[1])+")");
			table.finalizeLine();
			
			table.initNewLine().addColumn("Below "+ycLabel);
			for (int i=0; i<counts.size(); i++)
				table.addColumn(counts.get(i)[2]+" ("+pDF.format(fracts.get(i)[2])+")");
			table.finalizeLine();
			
			table.initNewLine().addColumn("Moment-weighted Filter %");
			for (double[] fract : momFracts)
				table.addColumn(pDF.format(fract[0]));
			table.finalizeLine();
			
			table.initNewLine().addColumn("Filtered Moment Change");
			for (double change : momChanges)
				table.addColumn((change >= 0d ? "+" : "")+pDF.format(change));
			table.finalizeLine();
			
			table.initNewLine().addColumn("Original Num Zeros");
			for (int i=0; i<counts.size(); i++)
				table.addColumn(origNumZeros.get(i)+" ("+pDF.format(origFractZeros.get(i))+")");
			table.finalizeLine();
			
			table.initNewLine().addColumn("Filtered Num Zeros");
			for (int i=0; i<counts.size(); i++)
				table.addColumn(filteredNumZeros.get(i)+" ("+pDF.format(filteredFractZeros.get(i))+")");
			table.finalizeLine();
			
			lines.add("### Summary Table, Yc="+oDF.format(yc)+"");
			lines.add(topLink); lines.add("");
			
			lines.addAll(table.build());
			lines.add("");
			
			// dm plots table
			table = MarkdownUtils.tableBuilder();
			
			table.initNewLine();
			for (NSHM23_DeformationModels dm : defModels)
				table.addColumn(dm.getShortName());
			table.finalizeLine();
			
			table.initNewLine();
			for (File mapPlot : mapPlots)
				table.addColumn("![Replacement Map]("+resourcesDir.getName()+"/"+mapPlot.getName()+")");
			table.finalizeLine();
			
			table.initNewLine();
			for (File mapPlot : scatterToOrigPlots)
				table.addColumn("![Scatter Plot]("+resourcesDir.getName()+"/"+mapPlot.getName()+")");
			table.finalizeLine();
			
			table.initNewLine();
			for (File mapPlot : scatterToMedianPlots)
				table.addColumn("![Scatter Plot]("+resourcesDir.getName()+"/"+mapPlot.getName()+")");
			table.finalizeLine();
			
			lines.add("### DM-Specific Plots, Yc="+oDF.format(yc)+"");
			lines.add(topLink); lines.add("");
			lines.add("The following table shows slip rate replacments in map view, with those above "+ycLabel+" in "
					+ "red, those below in blue, and those unchanged in gray. Then, slip rate scatters are shown "
					+ "for each model relative to the unfiltered version, and relative to the median model.");
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
			
			lines.add("### State-By-State Moment Analysis, Yc="+oDF.format(yc)+"");
			lines.add(topLink); lines.add("");
			
			lines.add("This table shows how the moment in the filtered model compares to the original defomation model, "
					+ "broken down by individual states.");
			lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			
			table.initNewLine().addColumn("State");
			for (NSHM23_DeformationModels dm : defModels)
				table.addColumn(dm.getShortName());
			table.finalizeLine();
			
			for (int i=0; i<stateNames.size(); i++) {
				table.initNewLine().addColumn("__"+stateNames.get(i)+"__");
				for (int d=0; d<defModels.length; d++) {
					double origMoment = stateOrigMomentsList.get(d)[i];
					double filteredMoment = stateFilteredMomentsList.get(d)[i];
					
					double change = (filteredMoment - origMoment)/origMoment;
					table.addColumn((change >= 0 ? "+" : "")+pDF.format(change));
				}
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
			
			lines.add("### Average Filtered Model Comparisons, Yc="+oDF.format(yc)+"");
			lines.add(topLink); lines.add("");
			
			lines.add("These maps compare the average model (after filtering) to models without filtering");
			lines.add("");
			
			// average models
			Map<Integer, List<MinisectionSlipRecord>> filteredAvgMinis = NSHM23_DeformationModels.averageMinisections(
					dmFilteredList, dmWeightsList);
			
			List<Double> ratiosToOrigAvg = new ArrayList<>();
			List<Double> diffsToOrigAvg = new ArrayList<>();
			
			List<Double> ratiosToReweightAvg = new ArrayList<>();
			List<Double> diffsToReweightAvg = new ArrayList<>();
			
			List<Double> ratiosToMedian = new ArrayList<>();
			List<Double> diffsToMedian = new ArrayList<>();
			
			List<FaultSection> sects = new ArrayList<>();
			for (Integer id : filteredAvgMinis.keySet()) {
				FaultSection origSect = geolSects.get(id);
				
				List<MinisectionSlipRecord> filteredRecs = filteredAvgMinis.get(id);
				List<MinisectionSlipRecord> origRecs = origAvgMiniMap.get(id);
				List<MinisectionSlipRecord> reweightRecs = reweightedAvgMiniMap.get(id);
				List<MinisectionSlipRecord> medianRecs = medianMiniMap.get(id);
				
				for (int i=0; i<filteredRecs.size(); i++) {
					double filtered = filteredRecs.get(i).slipRate;
					double orig = origRecs.get(i).slipRate;
					double reweight = reweightRecs.get(i).slipRate;
					double median = medianRecs.get(i).slipRate;
					
					ratiosToOrigAvg.add(filtered/orig);
					diffsToOrigAvg.add(filtered - orig);
					
					ratiosToReweightAvg.add(filtered/reweight);
					diffsToReweightAvg.add(filtered - reweight);
					
					ratiosToMedian.add(filtered/median);
					diffsToMedian.add(filtered - median);
					
					sects.add(miniToSect(filteredRecs.get(i), origSect, sects.size(), i));
				}
			}
			
			RupSetMapMaker mapMaker = new RupSetMapMaker(sects, mapReg);
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine("vs Original Average", "Vs Unfilterd Reweighted Average", "Vs Median");
			
			// ratios
			table.initNewLine();
			
			mapMaker.plotSectScalars(toLog10(ratiosToOrigAvg), ratioCPT, "Log10 (Filtered-Avg / Original-Avg)");
			String prefix = ycPrefix+"_ratio_vs_orig_avg";
			mapMaker.plot(resourcesDir, prefix, "Ratio of Filtered-Avg to Original-Avg");
			table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			mapMaker.plotSectScalars(toLog10(ratiosToReweightAvg), ratioCPT, "Log10 (Filtered-Avg / Unfiltered-Reweighted-Avg)");
			prefix = ycPrefix+"_ratio_vs_reweighted_avg";
			mapMaker.plot(resourcesDir, prefix, "Ratio of Filtered-Avg to Unfiltered-Reweighted-Avg");
			table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			mapMaker.plotSectScalars(toLog10(ratiosToMedian), ratioCPT, "Log10 (Filtered-Avg / Median)");
			prefix = ycPrefix+"_ratio_vs_median";
			mapMaker.plot(resourcesDir, prefix, "Ratio of Filtered-Avg to Median");
			table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			table.finalizeLine();
			
			// diffs
			table.initNewLine();

			mapMaker.plotSectScalars(diffsToOrigAvg, diffCPT, "Filtered-Avg - Original-Avg (mm/yr)");
			prefix = ycPrefix+"_diff_vs_orig_avg";
			mapMaker.plot(resourcesDir, prefix, "Diff between Filtered-Avg and Original-Avg");
			table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

			mapMaker.plotSectScalars(diffsToReweightAvg, diffCPT, "Filtered-Avg - Unfiltered-Reweighted-Avg (mm/yr)");
			prefix = ycPrefix+"_diff_vs_reweighted_avg";
			mapMaker.plot(resourcesDir, prefix, "Diff between Filtered-Avg and Unfiltered-Reweighted-Avg");
			table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

			mapMaker.plotSectScalars(diffsToMedian, diffCPT, "Filtered-Avg - Median (mm/yr)");
			prefix = ycPrefix+"_diff_vs_median";
			mapMaker.plot(resourcesDir, prefix, "Diff between Filtered-Avg and Median");
			table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

			table.finalizeLine();
			
			// scatter
			table.initNewLine();
			
			prefix = ycPrefix+"_scatter_vs_orig_avg";
			File plot = scatterPlot(filteredAvgMinis, "Filtered Average, "+yc, origAvgMiniMap, "Original Average",
					resourcesDir, prefix, Color.BLUE.darker(), "Filtered-Avg vs Original-Avg");
			table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
			
			prefix = ycPrefix+"_scatter_vs_reweighted_avg";
			plot = scatterPlot(filteredAvgMinis, "Filtered Average, "+yc, reweightedAvgMiniMap, "Unfiltered-Reweighted Average",
					resourcesDir, prefix, Color.RED.darker(), "Filtered-Avg vs Unfiltered-Reweighted-Avg");
			table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
			
			prefix = ycPrefix+"_scatter_vs_median";
			plot = scatterPlot(filteredAvgMinis, "Filtered Average, "+yc, medianMiniMap, "Median",
					resourcesDir, prefix, Color.GREEN.darker(), "Filtered-Avg vs Median");
			table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
			
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		lines.add("## Unfiltered Comparison");
		lines.add(topLink); lines.add("");
		
		lines.add("This section shows how the unfiltered models compare to each other.");
		lines.add("");
		
		List<Double> ratiosReweightToOrig = new ArrayList<>();
		List<Double> diffsReweightToOrig = new ArrayList<>();
		
		List<Double> ratiosOrigToMedian = new ArrayList<>();
		List<Double> diffsOrigToMedian = new ArrayList<>();
		
		List<Double> ratiosReweightToMedian = new ArrayList<>();
		List<Double> diffsReweightToMedian = new ArrayList<>();
		
		List<FaultSection> sects = new ArrayList<>();
		for (Integer id : reweightedAvgMiniMap.keySet()) {
			FaultSection origSect = geolSects.get(id);
			
			List<MinisectionSlipRecord> origRecs = origAvgMiniMap.get(id);
			List<MinisectionSlipRecord> reweightRecs = reweightedAvgMiniMap.get(id);
			List<MinisectionSlipRecord> medianRecs = medianMiniMap.get(id);
			
			for (int i=0; i<origRecs.size(); i++) {
				double orig = origRecs.get(i).slipRate;
				double reweight = reweightRecs.get(i).slipRate;
				double median = medianRecs.get(i).slipRate;
				
				ratiosReweightToOrig.add(reweight/orig);
				diffsReweightToOrig.add(reweight - orig);
				
				ratiosReweightToMedian.add(reweight/median);
				diffsReweightToMedian.add(reweight - median);
				
				ratiosOrigToMedian.add(orig/median);
				diffsOrigToMedian.add(orig - median);
				
				sects.add(miniToSect(origRecs.get(i), origSect, sects.size(), i));
			}
		}
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(sects, mapReg);
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("Unfiltered Reweighted vs Original", "Unfiltered Reweighted vs Median", "Original vs Median");
		
		// ratios
		table.initNewLine();
		
		mapMaker.plotSectScalars(toLog10(ratiosReweightToOrig), ratioCPT, "Log10 (Unfiltered-Reweighted-Avg / Original-Avg)");
		String prefix = "ratio_reweight_vs_orig";
		mapMaker.plot(resourcesDir, prefix, "Ratio of Unfiltered-Reweighted-Avg to Original-Avg");
		table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
		
		mapMaker.plotSectScalars(toLog10(ratiosReweightToMedian), ratioCPT, "Log10 (Unfiltered-Reweighted-Avg / Median)");
		prefix = "ratio_reweight_vs_median";
		mapMaker.plot(resourcesDir, prefix, "Ratio of Unfiltered-Reweighted-Avg to Median");
		table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
		
		mapMaker.plotSectScalars(toLog10(ratiosOrigToMedian), ratioCPT, "Log10 (Original-Avg / Median)");
		prefix = "ratio_origt_vs_median";
		mapMaker.plot(resourcesDir, prefix, "Ratio of Original-Avg to Median");
		table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
		
		table.finalizeLine();
		
		// diffs
		table.initNewLine();

		mapMaker.plotSectScalars(diffsReweightToOrig, diffCPT, "Unfiltered-Reweighted-Avg - Original-Avg (mm/yr)");
		prefix = "diff_reweight_vs_orig";
		mapMaker.plot(resourcesDir, prefix, "Diff between Unfiltered-Reweighted-Avg and Original-Avg");
		table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

		mapMaker.plotSectScalars(diffsReweightToMedian, diffCPT, "Unfiltered-Reweighted-Avg - Median (mm/yr)");
		prefix = "diff_reweight_vs_median";
		mapMaker.plot(resourcesDir, prefix, "Diff between Unfiltered-Reweighted-Avg and Median");
		table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

		mapMaker.plotSectScalars(diffsOrigToMedian, diffCPT, "Original-Avg - Median (mm/yr)");
		prefix = "diff_origt_vs_median";
		mapMaker.plot(resourcesDir, prefix, "Diff between Original-Avg and Median");
		table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

		table.finalizeLine();
		
		// scatter
		table.initNewLine();
		
		prefix = "scatter_reweight_vs_orig";
		File plot = scatterPlot(reweightedAvgMiniMap, "Unfiltered-Reweighted-Avg", origAvgMiniMap, "Original Average",
				resourcesDir, prefix, Color.BLUE.darker(), "Unfiltered-Reweighted-Avg vs Original-Avg");
		table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
		
		prefix = "scatter_reweight_vs_median";
		plot = scatterPlot(reweightedAvgMiniMap, "Unfiltered-Reweighted-Avg", medianMiniMap, "Median",
				resourcesDir, prefix, Color.GREEN.darker(), "Unfiltered-Reweighted-Avg vs Median");
		table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
		
		prefix = "scatter_origt_vs_median";
		plot = scatterPlot(origAvgMiniMap, "Original-Avg", medianMiniMap, "Median",
				resourcesDir, prefix, Color.GREEN.darker(), "Original-Avg vs Median");
		table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
		
		table.finalizeLine();
		
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Individual Model Comparisons (without replacement)");
		lines.add(topLink); lines.add("");
		
		for (int d=0; d<defModels.length; d++) {
			NSHM23_DeformationModels dm = defModels[d];
			
			lines.add("### "+dm.getShortName());
			lines.add(topLink); lines.add("");
			
			String dmPrefix = "dm_"+dm.name();
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine(dm.getShortName()+" Vs Average", dm.getShortName()+" vs Reweighted Average",
					dm.getShortName()+" vs Median");
			
			Map<Integer, List<MinisectionSlipRecord>> dmMinis = dmOrigMiniMaps.get(d);
			
			List<Double> ratiosToOrig = new ArrayList<>();
			List<Double> diffsToOrig = new ArrayList<>();
			
			List<Double> ratiosToReweight = new ArrayList<>();
			List<Double> diffsToReweight = new ArrayList<>();
			
			List<Double> ratiosToMedian = new ArrayList<>();
			List<Double> diffsToMedian = new ArrayList<>();
			
			sects = new ArrayList<>();
			for (Integer id : reweightedAvgMiniMap.keySet()) {
				FaultSection origSect = geolSects.get(id);
				
				List<MinisectionSlipRecord> dmRecs = dmMinis.get(id);
				List<MinisectionSlipRecord> origRecs = origAvgMiniMap.get(id);
				List<MinisectionSlipRecord> reweightRecs = reweightedAvgMiniMap.get(id);
				List<MinisectionSlipRecord> medianRecs = medianMiniMap.get(id);
				
				for (int i=0; i<origRecs.size(); i++) {
					double dmVal = dmRecs.get(i).slipRate;
					double orig = origRecs.get(i).slipRate;
					double reweight = reweightRecs.get(i).slipRate;
					double median = medianRecs.get(i).slipRate;
					
					ratiosToOrig.add(dmVal/orig);
					diffsToOrig.add(dmVal - orig);
					
					ratiosToReweight.add(dmVal/reweight);
					diffsToReweight.add(dmVal - reweight);
					
					ratiosToMedian.add(dmVal/median);
					diffsToMedian.add(dmVal - median);
					
					sects.add(miniToSect(dmRecs.get(i), origSect, sects.size(), i));
				}
			}
			
			// ratios
			table.initNewLine();
			
			mapMaker.plotSectScalars(toLog10(ratiosToOrig), ratioCPT, "Log10 ("+dm.getShortName()+" / Original-Avg)");
			prefix = dmPrefix+"_ratio_vs_orig_avg";
			mapMaker.plot(resourcesDir, prefix, "Ratio of "+dm.getShortName()+" to Original-Avg");
			table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			mapMaker.plotSectScalars(toLog10(ratiosToReweight), ratioCPT, "Log10 ("+dm.getShortName()+" / Reweighted-Avg)");
			prefix = dmPrefix+"_ratio_vs_reweight_avg";
			mapMaker.plot(resourcesDir, prefix, "Ratio of "+dm.getShortName()+" to Reweighted-Avg");
			table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			mapMaker.plotSectScalars(toLog10(ratiosToMedian), ratioCPT, "Log10 ("+dm.getShortName()+" / Median)");
			prefix = dmPrefix+"_ratio_vs_median";
			mapMaker.plot(resourcesDir, prefix, "Ratio of "+dm.getShortName()+" to Median");
			table.addColumn("![Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			table.finalizeLine();
			
			// diffs
			table.initNewLine();
			
			mapMaker.plotSectScalars(toLog10(diffsToOrig), diffCPT, dm.getShortName()+" - Original-Avg");
			prefix = dmPrefix+"_diff_vs_orig_avg";
			mapMaker.plot(resourcesDir, prefix, "Diff of "+dm.getShortName()+" and Original-Avg");
			table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			mapMaker.plotSectScalars(toLog10(diffsToReweight), diffCPT, dm.getShortName()+" - Reweighted-Avg");
			prefix = dmPrefix+"_diff_vs_reweight_avg";
			mapMaker.plot(resourcesDir, prefix, "Diff of "+dm.getShortName()+" and Reweighted-Avg");
			table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

			mapMaker.plotSectScalars(toLog10(diffsToMedian), diffCPT, dm.getShortName()+" - Median");
			prefix = dmPrefix+"_diff_vs_median";
			mapMaker.plot(resourcesDir, prefix, "Diff of "+dm.getShortName()+" and Median");
			table.addColumn("![diff]("+resourcesDir.getName()+"/"+prefix+".png)");

			table.finalizeLine();
			
			// scatter
			table.initNewLine();
			
			prefix = dmPrefix+"_scatter_vs_orig_avg";
			plot = scatterPlot(dmMinis, dm.getShortName(), origAvgMiniMap, "Original Average",
					resourcesDir, prefix, Color.BLUE.darker(), dm.getShortName()+" vs Original-Avg");
			table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
			
			prefix = dmPrefix+"_scatter_vs_reweight_avg";
			plot = scatterPlot(dmMinis, dm.getShortName(), reweightedAvgMiniMap, "Reweighted Average",
					resourcesDir, prefix, Color.RED.darker(), dm.getShortName()+" vs Reweighted-Avg");
			table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
			
			prefix = dmPrefix+"_scatter_vs_median";
			plot = scatterPlot(dmMinis, dm.getShortName(), medianMiniMap, "Median",
					resourcesDir, prefix, Color.GREEN.darker(), dm.getShortName()+" vs Median");
			table.addColumn("![scatter]("+resourcesDir.getName()+"/"+plot.getName()+")");
			
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, reportDir);
	}
	
	private static FaultSection miniToSect(MinisectionSlipRecord mini, FaultSection refSect, int index, int minisectionIndex) {
		FaultSectionPrefData sect = new FaultSectionPrefData();
		sect.setSectionId(index);
		sect.setSectionName(refSect.getSectionName()+" Minisection "+minisectionIndex);
		sect.setParentSectionId(refSect.getSectionId());
		sect.setParentSectionName(refSect.getSectionName());
		
		FaultTrace trace = new FaultTrace(sect.getSectionName());
		trace.add(mini.startLoc);
		trace.add(mini.endLoc);
		sect.setFaultTrace(trace);
		
		sect.setAveUpperDepth(refSect.getOrigAveUpperDepth());
		sect.setAveLowerDepth(refSect.getAveLowerDepth());
		sect.setAveDip(refSect.getAveDip());
		sect.setAveRake(refSect.getAveRake());
		sect.setAveSlipRate(mini.slipRate);
		
		return sect;
	}
	
	private static File scatterPlot(Map<Integer, List<MinisectionSlipRecord>> filteredMinis, String filteredName,
			Map<Integer, List<MinisectionSlipRecord>> refMinis, String refName, File outputDir, String prefix,
			Color color, String title) throws IOException {
		Range range = new Range(1e-3, 1e2);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		for (Integer id : filteredMinis.keySet()) {
			if (!refMinis.containsKey(id))
				continue;
			List<MinisectionSlipRecord> filteredRecs = filteredMinis.get(id);
			List<MinisectionSlipRecord> refRecs = refMinis.get(id);
			
			for (int i=0; i<filteredRecs.size(); i++)
				xy.set(Double.max(refRecs.get(i).slipRate, range.getLowerBound()), Double.max(filteredRecs.get(i).slipRate, range.getLowerBound()));
		}
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, color));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, refName, filteredName);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, true, true, range, range);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 600, false, true, true, false);
		
		return new File(outputDir, prefix+".png");
	}
	
	private static List<Double> toLog10(List<Double> vals) {
		List<Double> ret = new ArrayList<>(vals.size());
		for (double val : vals)
			ret.add(Math.log10(val));
		return ret;
	}

}
