package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.CubedGriddedRegion;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.nshm23.MomentRateCompNSHM18;

public class WUS_HazardChangePageGen {

	public static void main(String[] args) throws IOException {
		
		String entryName = "mean_map_pga_TWO_IN_50.txt";
		String wrapperEntryName = "map_pga_TWO_IN_50.txt";
		String hazLabel = "PGA, 2% in 50 yrs";
		String dirPrefix = "pga_2in50";
		
		// load hazard maps
		File nshm23HazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		GriddedGeoDataSet nshm23Hazard = CA_HazardChangeFigures.loadXYZ(nshm23HazardFile, entryName);
		
		File nshm23GridHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_only_0.1deg.zip");
		GriddedGeoDataSet nshm23GridHazard = CA_HazardChangeFigures.loadXYZ(nshm23GridHazardFile, entryName);
		Preconditions.checkState(nshm23GridHazard.size() == nshm23Hazard.size());
//		GriddedGeoDataSet nshm23GridHazard = null;
		
		File nshm18_23gridHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_04-nshm18-grid_src_from_23-hazard-ask2014-0.1deg-noSub/results_hazard.zip");
		GriddedGeoDataSet nshm18_23gridHazard = CA_HazardChangeFigures.loadXYZ(nshm18_23gridHazardFile, wrapperEntryName);
		Preconditions.checkState(nshm18_23gridHazard.size() == nshm23Hazard.size());
		
		File nshm18HazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_17-nshm18-hazard-ask2014-0.1deg-noSub/results_hazard.zip");
		GriddedGeoDataSet nshm18Hazard = CA_HazardChangeFigures.loadXYZ(nshm18HazardFile, wrapperEntryName);
		Preconditions.checkState(nshm18Hazard.size() == nshm23Hazard.size());
		
		File nshm18GridHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_18-nshm18-hazard-ask2014-0.1deg-noSub-griddedOnly/results_hazard.zip");
		GriddedGeoDataSet nshm18GridHazard = CA_HazardChangeFigures.loadXYZ(nshm18GridHazardFile, wrapperEntryName);
		Preconditions.checkState(nshm18GridHazard.size() == nshm23Hazard.size());
//		GriddedGeoDataSet nshm18GridHazard = null;
		
		// load moment rates from NSHM18 and NSHM23 deformation models
		MomentRateCompNSHM18.LINEAR_RAMP = true;
		MomentRateCompNSHM18.GEO_ONLY = false;
		
		GriddedRegion gridReg = nshm23Hazard.getRegion();
		
		CubedGriddedRegion cgr = new CubedGriddedRegion(gridReg);
		GriddedGeoDataSet momentRates23 = MomentRateCompNSHM18.getMomentRatesNSHM23(gridReg, cgr);
		List<FaultSection> sects18 = new ArrayList<>();
		GriddedGeoDataSet momentRates18 = MomentRateCompNSHM18.getMomentRatesNSHM18(gridReg, cgr, sects18);
		
		System.out.println("Total NSHM18 moment rate: "+(float)momentRates18.getSumZ());
		System.out.println("Total NSHM23 moment rate: "+(float)momentRates23.getSumZ());
		
		GriddedGeoDataSet momentRatio = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<momentRatio.size(); i++) {
			double rate23 = momentRates23.get(i);
			double rate18 = momentRates18.get(i);
			
			if (rate23 == 0d && rate18 == 0d)
				momentRatio.set(i, Double.NaN);
			else
				momentRatio.set(i, rate23/rate18);
		}
		
		double[] momentThresholds = { Double.POSITIVE_INFINITY, 10d, 5d, 2d, 1.5d };
		
		GriddedGeoDataSet faultHazardRatio = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet faultHazardPDiff = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet faultHazardDiff = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet fullHazardRatio = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet fullHazardPDiff = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet fullHazardDiff = new GriddedGeoDataSet(gridReg, false);
		
		boolean gridded = nshm23GridHazard != null && nshm18GridHazard != null;
		GriddedGeoDataSet griddedHazardRatio, griddedHazardPDiff, griddedHazardDiff;
		if (gridded) {
			griddedHazardRatio = new GriddedGeoDataSet(gridReg, false);
			griddedHazardPDiff = new GriddedGeoDataSet(gridReg, false);
			griddedHazardDiff = new GriddedGeoDataSet(gridReg, false);
		} else {
			griddedHazardRatio = null;
			griddedHazardPDiff = null;
			griddedHazardDiff = null;
		}
		
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			double full23 = nshm23Hazard.get(i);
			double fault18 = nshm18_23gridHazard.get(i);
			double full18 = nshm18Hazard.get(i);
			
			// these hold gridded seismicity constant (at '23)
			double faultRatio = full23/fault18;
			double faultDiff = full23 - fault18;
			
			faultHazardRatio.set(i, faultRatio);
			faultHazardPDiff.set(i, 100d*(full23 - fault18)/fault18);
			faultHazardDiff.set(i, faultDiff);
			
			double fullRatio = full23/full18;
			double fullDiff = full23 - full18;
			
			fullHazardRatio.set(i, fullRatio);
			fullHazardPDiff.set(i, 100d*(full23 - full18)/full18);
			fullHazardDiff.set(i, fullDiff);
			
			if (gridded) {
				double grid23 = nshm23GridHazard.get(i);
				double grid18 = nshm18GridHazard.get(i);
				
				// add fault contributions to highlight only gridded
				double faultPortion = Math.max(0, full23 - grid23);
				grid23 += faultPortion;
				grid18 += faultPortion;
				
				double griddedRatio = grid23 / grid18;
				double griddedDiff = grid23 - grid18;
				
				griddedHazardRatio.set(i, griddedRatio);
				griddedHazardPDiff.set(i, 100d*(grid23 - grid18)/grid18);
				griddedHazardDiff.set(i, griddedDiff);
			}
		}
		
		Color transparent = new Color(255, 255, 255, 0);
		
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		hazCPT.setNanColor(transparent);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(transparent);
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2d, 0.2d);
		diffCPT.setNanColor(transparent);
		
		List<FaultSection> sectsWithIndexes18 = new ArrayList<>();
		for (FaultSection sect : sects18) {
			sect = sect.clone();
			sect.setSectionId(sectsWithIndexes18.size());
			sectsWithIndexes18.add(sect);
		}
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v2), gridReg);
		RupSetMapMaker mapMaker18 = new RupSetMapMaker(sectsWithIndexes18, gridReg);
		
		for (RupSetMapMaker map : new RupSetMapMaker[] {mapMaker, mapMaker18}) {
			map.setWriteGeoJSON(false);
			map.setWritePDFs(true);
			map.setSectOutlineChar(null);
			map.setRegionOutlineChar(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, new Color(0, 0, 0, 180)));
			map.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(100, 100, 100, 127)));
		}
		
		File outputDir = new File(nshm23HazardFile.getParentFile(), "hazard_comparisons_nshm18_"+dirPrefix);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Hazard Comparisons with NSHM18");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Fault-Based Hazard Comparisons");
		lines.add(topLink); lines.add("");
		
		lines.add("This section compares NSHM23 with NSHM18, but holding the gridded seismicity model constant in "
				+ "order to focus on fault-based changes (due both to ingredient and methodological changes).");
		lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("NSHM23", "NSHM18 w/ NSHM23 Gridded");
		
		table.initNewLine();
		
		mapMaker.plotXYZData(asLog10(nshm23Hazard), hazCPT, "NSHM23, "+hazLabel+" (g)");
		mapMaker.plot(resourcesDir, "hazard_nshm23", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/hazard_nshm23.png)");
		
		mapMaker18.plotXYZData(asLog10(nshm18_23gridHazard), hazCPT, "NSHM18 (w/ '23 Grid), "+hazLabel+" (g)");
		mapMaker18.plot(resourcesDir, "hazard_nshm18_23grid", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/hazard_nshm18_23grid.png)");
		
		table.finalizeLine();
		
		table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
		
		table.initNewLine();
		
		mapMaker.plotXYZData(faultHazardPDiff, pDiffCPT, "NSHM23 vs NSHM18 (w/ '23 Grid), % Change, "+hazLabel);
		mapMaker.plot(resourcesDir, "fault_hazard_pDiff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/fault_hazard_pDiff.png)");
		
		mapMaker.plotXYZData(faultHazardDiff, diffCPT, "NSHM23 - NSHM18 (w/ '23 Grid), "+hazLabel+" (g)");
		mapMaker.plot(resourcesDir, "fault_hazard_diff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/fault_hazard_diff.png)");
		
		table.finalizeLine();
		
		lines.addAll(table.build());
		lines.add("");
		
		ChangeStats faultStats = new ChangeStats();
		for (int i=0; i<nshm23Hazard.size(); i++)
			faultStats.addValue(nshm23Hazard.get(i), nshm18_23gridHazard.get(i));
		
		table = MarkdownUtils.tableBuilder();
		table.addLine(faultStats.pDiffTableHeader());
		table.addLine(faultStats.pDiffTableLine());
		table.addLine(faultStats.diffTableHeader());
		table.addLine(faultStats.diffTableLine());
		
		lines.add("__Hazard Comparison Statistics:__");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("### Fault and Deformation Model Moment Changes");
		lines.add(topLink); lines.add("");
		
		lines.add("This section shows how fault-based deformation model moment changed between NSHM23 and NSHM18. We "
				+ "will use this map to mask hazard changes, highlighting only areas where changes are primarily due to "
				+ "methdological differences. In the ratio map, areas where moment exists in NSHM23 but not in NSHM18 "
				+ "are shown in yellow (i.e., a fault was added), and those with moment in NSHM18 but not in NSHM23 are "
				+ "shown in green (i.e., a fault was removed).");
		lines.add("");
		
		CPT logMomentRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
		logMomentRatioCPT.setAboveMaxColor(Color.YELLOW);
		logMomentRatioCPT.setBelowMinColor(Color.GREEN);
		logMomentRatioCPT.setNanColor(transparent);
		
		CPT momentDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-5e15, 5e15);
		momentDiffCPT.setNanColor(transparent);
		
		GriddedGeoDataSet logMomentRatio = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet momentDiff = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<logMomentRatio.size(); i++) {
			double moment23 = momentRates23.get(i);
			double moment18 = momentRates18.get(i);
			
			if (moment23 == 0d && moment18 == 0d) {
				logMomentRatio.set(i, Double.NaN);
				momentDiff.set(i, Double.NaN);
			} else {
				double ratio = momentRatio.get(i);
				double logRatio = Math.log10(ratio);
				if (Double.isFinite(logRatio)) {
					// neither one was zero
					// check to see if we should saturate it
					if ((float)logRatio > logMomentRatioCPT.getMaxValue())
						logRatio = (double)logMomentRatioCPT.getMaxValue();
					else if ((float)logRatio < logMomentRatioCPT.getMinValue())
						logRatio = (double)logMomentRatioCPT.getMinValue();
				}
				logMomentRatio.set(i, logRatio);
				
				momentDiff.set(i, moment23 - moment18);
			}
		}
		
		table = MarkdownUtils.tableBuilder();
		
		table.addLine("Ratio", "Difference");
		
		table.initNewLine();
		
		mapMaker.plotXYZData(logMomentRatio, logMomentRatioCPT, "Log10 (NSHM23 Fault Moment) / (NSHM18 Fault Moment)");
		mapMaker.plot(resourcesDir, "moment_ratio", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/moment_ratio.png)");
		
		mapMaker.plotXYZData(momentDiff, momentDiffCPT, "NSHM23 Fault Moment) - (NSHM18 Fault Moment)");
		mapMaker.plot(resourcesDir, "moment_diff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/moment_diff.png)");
		
		table.finalizeLine();
		
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("### Fault-Based Hazard Comparisons, Moment Chage Masked");
		lines.add(topLink); lines.add("");
		
		lines.add("This section masks hazard changes, only showing them near faults where NSHM23 fault-moment is within "
				+ "a given factor of NSHM18 fault-moment. Various thresholds are used, starting with one that only "
				+ "includes areas with faults in both models, and then the threshold decreases to only show areas where "
				+ "moment is similar in the two models. Areas where the sign of the hazard change differs from the "
				+ "sign of moment change are shown regardless of the threshold.");
		lines.add("");
		
		TableBuilder threshMapTable = MarkdownUtils.tableBuilder();
		threshMapTable.addLine("Ratio", "Difference");
		
		List<String> threshLabels = new ArrayList<>();
		List<ChangeStats> threshStats = new ArrayList<>();
		
		for (double threshold : momentThresholds) {
			String threshLabel, threshPrefix;
			if (Double.isFinite(threshold)) {
				threshLabel = oDF.format(threshold)+"x Moment Threshold";
				threshPrefix = "mom_thresh_"+oDF.format(threshold);
			} else {
				threshLabel = "Any Moment Threshold";
				threshPrefix = "mom_thresh_inf";
			}
			
			threshLabels.add(threshLabel);
			
			double lowerThresh = 1d/threshold;
			
			ChangeStats stats = new ChangeStats();
			threshStats.add(stats);
			GriddedGeoDataSet maskedPDiff = new GriddedGeoDataSet(gridReg, false);
			GriddedGeoDataSet maskedDiff = new GriddedGeoDataSet(gridReg, false);
			for (int i=0; i<nshm23Hazard.size(); i++) {
				double moRatio = momentRatio.get(i);
				boolean skip = Double.isNaN(moRatio) || (Double.isInfinite(threshold) && Double.isInfinite(moRatio));
				if (!skip) {
					// check the threshold
					boolean outside = moRatio > threshold || moRatio < lowerThresh;
					if (outside) {
						// make sure the sign is the same
						double hazDiff = faultHazardDiff.get(i);
						double moDiff = momentDiff.get(i);
						if (hazDiff > 0 && moDiff > 0)
							skip = true;
						else if (hazDiff < 0 && moDiff < 0)
							skip = true;
					}
				}
				if (skip) {
					maskedPDiff.set(i, Double.NaN);
					maskedDiff.set(i, Double.NaN);
				} else {
					stats.addValue(nshm23Hazard.get(i), nshm18_23gridHazard.get(i));
					maskedPDiff.set(i, faultHazardPDiff.get(i));
					maskedDiff.set(i, faultHazardDiff.get(i));
				}
			}
			
			threshMapTable.initNewLine();
			
			mapMaker.plotXYZData(maskedPDiff, pDiffCPT, "NSHM23 vs NSHM18 (w/ '23 Grid), % Change, "+hazLabel);
			mapMaker.plot(resourcesDir, "fault_hazard_pDiff_"+threshPrefix, threshLabel);
			threshMapTable.addColumn("![Map]("+resourcesDir.getName()+"/fault_hazard_pDiff_"+threshPrefix+".png)");
			
			mapMaker.plotXYZData(maskedDiff, diffCPT, "NSHM23 - NSHM18 (w/ '23 Grid), "+hazLabel+" (g)");
			mapMaker.plot(resourcesDir, "fault_hazard_diff_"+threshPrefix, threshLabel);
			threshMapTable.addColumn("![Map]("+resourcesDir.getName()+"/fault_hazard_diff_"+threshPrefix+".png)");
			
			threshMapTable.finalizeLine();
		}
		
		lines.addAll(threshMapTable.build());
		lines.add("");
		
		String withinLabel = MarkdownUtils.boldCentered("% Locs Within Moment Threshold");
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumns("Threshold", withinLabel);
		table.addColumns(threshStats.get(0).pDiffTableHeader());
		table.finalizeLine();
		for (int i=0; i<threshLabels.size(); i++) {
			table.initNewLine();
			table.addColumn(threshLabels.get(i));
			ChangeStats stats = threshStats.get(i);
			table.addColumn(pDF.format((double)stats.totalNum/(double)gridReg.getNodeCount()));
			table.addColumns(stats.pDiffTableLine());
			table.finalizeLine();
		}
		table.initNewLine();
		table.addColumns("Threshold", withinLabel);
		table.addColumns(threshStats.get(0).diffTableHeader());
		table.finalizeLine();
		for (int i=0; i<threshLabels.size(); i++) {
			table.initNewLine();
			table.addColumn(threshLabels.get(i));
			ChangeStats stats = threshStats.get(i);
			table.addColumn(pDF.format((double)stats.totalNum/(double)gridReg.getNodeCount()));
			table.addColumns(stats.diffTableLine());
			table.finalizeLine();
		}
		lines.add("__Near-Fault Hazard Comparison Statistics:__");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		ChangeStats griddedStats = null;
		if (gridded) {
			// gridded differences
			lines.add("## Gridded Seismicity Model Hazard Changes");
			lines.add(topLink); lines.add("");
			
			lines.add("These plots compare the gridded seismicity components of the models. For the ratio and "
					+ "difference maps, the contributions of NSHM23 fault sources is held constant (added to both "
					+ "models).");
			lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine("NSHM23, Gridded Only", "NSHM18, Gridded Only");
			
			table.initNewLine();
			
			mapMaker.plotXYZData(asLog10(nshm23GridHazard), hazCPT, "NSHM23, Gridded Only, "+hazLabel+" (g)");
			mapMaker.plot(resourcesDir, "hazard_gridded_nshm23", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/hazard_gridded_nshm23.png)");
			
			mapMaker18.plotXYZData(asLog10(nshm18GridHazard), hazCPT, "NSHM18, Gridded Only, "+hazLabel+" (g)");
			mapMaker18.plot(resourcesDir, "hazard_gridded_nshm18", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/hazard_gridded_nshm18.png)");
			
			table.finalizeLine();
			
			table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
			
			table.initNewLine();
			
			mapMaker.plotXYZData(griddedHazardPDiff, pDiffCPT, "NSHM23 vs NSHM18, Gridded Only, % Change, "+hazLabel);
			mapMaker.plot(resourcesDir, "gridded_hazard_pDiff", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_pDiff.png)");
			
			mapMaker.plotXYZData(griddedHazardDiff, diffCPT, "NSHM23 - NSHM18, Gridded Only, "+hazLabel+" (g)");
			mapMaker.plot(resourcesDir, "gridded_hazard_diff", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_diff.png)");
			
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			griddedStats = new ChangeStats();
			for (int i=0; i<nshm23Hazard.size(); i++)
				griddedStats.addValue(nshm23Hazard.get(i), nshm18Hazard.get(i));
			
			table = MarkdownUtils.tableBuilder();
			table.addLine(griddedStats.pDiffTableHeader());
			table.addLine(griddedStats.pDiffTableLine());
			table.addLine(griddedStats.diffTableHeader());
			table.addLine(griddedStats.diffTableLine());
			
			lines.add("__Hazard Comparison Statistics:__");
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
			
			// now mask outside of faults
			lines.add("### Gridded Seismicity Model Hazard Changes, Away From Faults");
			lines.add(topLink); lines.add("");
			
			lines.add("Here, we mask hazard changes only showing them away from faults in either model.");
			lines.add("");
			
			ChangeStats stats = new ChangeStats();
			GriddedGeoDataSet maskedPDiff = new GriddedGeoDataSet(gridReg, false);
			GriddedGeoDataSet maskedDiff = new GriddedGeoDataSet(gridReg, false);
			for (int i=0; i<nshm23Hazard.size(); i++) {
				double mo23 = momentRates23.get(i);
				double mo18 = momentRates18.get(i);
				if (mo23 > 0 || mo18 > 0) {
					maskedPDiff.set(i, Double.NaN);
					maskedDiff.set(i, Double.NaN);
				} else {
					double haz23 = nshm23Hazard.get(i);
					double grid23 = nshm23GridHazard.get(i);
					double fault23 = Math.max(0, haz23 - grid23);
					double grid18 = nshm18GridHazard.get(i);
					stats.addValue(haz23, grid18+fault23);
					maskedPDiff.set(i, griddedHazardPDiff.get(i));
					maskedDiff.set(i, griddedHazardDiff.get(i));
				}
			}
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("Ratio", "Difference");
			
			table.initNewLine();
			
			mapMaker.plotXYZData(maskedPDiff, pDiffCPT, "NSHM23 vs NSHM18, Gridded Only, % Change, "+hazLabel);
			mapMaker.plot(resourcesDir, "gridded_hazard_pDiff_away_from_faults", "Gridded Hazard Away From Faults");
			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_pDiff_away_from_faults.png)");
			
			mapMaker.plotXYZData(maskedDiff, diffCPT, "NSHM23 - NSHM18, Gridded Only, "+hazLabel+" (g)");
			mapMaker.plot(resourcesDir, "gridded_hazard_diff_away_from_faults", "Gridded Hazard Away From Faults");
			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_diff_away_from_faults.png)");
			
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			withinLabel = MarkdownUtils.boldCentered("% Locs Away From Faults");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn(withinLabel);
			table.addColumns(stats.pDiffTableHeader());
			table.finalizeLine().initNewLine();
			table.addColumn(pDF.format((double)stats.totalNum/(double)maskedPDiff.size()));
			table.addColumns(stats.pDiffTableLine());
			table.finalizeLine();
			table.initNewLine();
			table.addColumn(withinLabel);
			table.addColumns(stats.diffTableHeader());
			table.finalizeLine().initNewLine();
			table.addColumn(pDF.format((double)stats.totalNum/(double)maskedPDiff.size()));
			table.addColumns(stats.diffTableLine());
			table.finalizeLine();
			lines.add("__Gridded Hazard Comparison Statistics Away From Faults:__");
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
		}
		
		// full hazard difference
		lines.add("## Complete Model Hazard Changes");
		lines.add(topLink); lines.add("");
		
		lines.add("These plots compare the full models, including both fault and gridded seismicity changes (both "
				+ "ingredients and methodology).");
		lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		
		table.addLine("Full NSHM23", "Full NSHM18");
		
		table.initNewLine();
		
		// this is already plotted
		table.addColumn("![Map]("+resourcesDir.getName()+"/hazard_nshm23.png)");
		
		mapMaker18.plotXYZData(asLog10(nshm18Hazard), hazCPT, "NSHM18, "+hazLabel+" (g)");
		mapMaker18.plot(resourcesDir, "hazard_nshm18", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/hazard_nshm18.png)");
		
		table.finalizeLine();
		
		table.addLine(MarkdownUtils.boldCentered("Full Ratio"), MarkdownUtils.boldCentered("Full Difference"));
		
		table.initNewLine();
		
		mapMaker.plotXYZData(fullHazardPDiff, pDiffCPT, "NSHM23 vs NSHM18, % Change, "+hazLabel);
		mapMaker.plot(resourcesDir, "full_hazard_pDiff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/full_hazard_pDiff.png)");
		
		mapMaker.plotXYZData(fullHazardDiff, diffCPT, "NSHM23 - NSHM18, "+hazLabel+" (g)");
		mapMaker.plot(resourcesDir, "full_hazard_diff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/full_hazard_diff.png)");
		
		table.finalizeLine();
		
		// add already plotted fault-only comparison
		table.addLine(MarkdownUtils.boldCentered("Fault-Only Ratio"), MarkdownUtils.boldCentered("Fault-Only Difference"));
		
		table.initNewLine();
		table.addColumn("![Map]("+resourcesDir.getName()+"/fault_hazard_pDiff.png)");
		table.addColumn("![Map]("+resourcesDir.getName()+"/fault_hazard_diff.png)");
		table.finalizeLine();
		
		if (gridded) {
			// add already plotted gridded-only comparison
			table.addLine(MarkdownUtils.boldCentered("Gridded-Only Ratio"), MarkdownUtils.boldCentered("Gridded-Only Difference"));
			
			table.initNewLine();
			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_pDiff.png)");
			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_diff.png)");
			
			table.finalizeLine();
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		ChangeStats fullStats = new ChangeStats();
		for (int i=0; i<nshm23Hazard.size(); i++)
			fullStats.addValue(nshm23Hazard.get(i), nshm18Hazard.get(i));
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine().addColumn("");
		table.addColumns(fullStats.pDiffTableHeader()).finalizeLine();
		
		table.initNewLine().addColumn("Full Models");
		table.addColumns(fullStats.pDiffTableLine()).finalizeLine();
		table.initNewLine().addColumn("Faults Only");
		table.addColumns(faultStats.pDiffTableLine()).finalizeLine();
		if (gridded) {
			table.initNewLine().addColumn("Gridded Only");
			table.addColumns(griddedStats.pDiffTableLine()).finalizeLine();
		}
		
		table.initNewLine().addColumn("");
		table.addColumns(fullStats.diffTableHeader()).finalizeLine();
		
		table.initNewLine().addColumn("Full Models");
		table.addColumns(fullStats.diffTableLine()).finalizeLine();
		table.initNewLine().addColumn("Faults Only");
		table.addColumns(faultStats.diffTableLine()).finalizeLine();
		if (gridded) {
			table.initNewLine().addColumn("Gridded Only");
			table.addColumns(griddedStats.diffTableLine()).finalizeLine();
		}
		
		lines.add("__Hazard Comparison Statistics:__");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		if (gridded) {
			// add fault/grid attribution plots
			
			CPT cork = GMT_CPT_Files.DIVERGING_CORK_UNIFORM.instance();
			CPT broc = GMT_CPT_Files.DIVERGING_BROC_UNIFORM.instance();
			// they get too dark at the edges (hard to tell color), trim them
			cork = cork.rescale(-1d, 1d);
			broc = broc.rescale(-1d, 1d);
			CPT modCork = new CPT();
			CPT modBroc = new CPT();
			for (double v=-0.81; (float)v<=0.8f; v+=0.01) {
				float start = (float)(v-0.01);
				Color startColorCorc = cork.getColor(start);
				Color startColorBroc = broc.getColor(start);
				float end = (float)v;
				Color endColorCorc = cork.getColor(end);
				Color endColorBroc = broc.getColor(end);
				modCork.add(new CPTVal(start, startColorCorc, end, endColorCorc));
				modBroc.add(new CPTVal(start, startColorBroc, end, endColorBroc));
			}
			modCork.setBelowMinColor(modCork.getMinColor());
			modCork.setAboveMaxColor(modCork.getMaxColor());
			modBroc.setBelowMinColor(modBroc.getMinColor());
			modBroc.setAboveMaxColor(modBroc.getMaxColor());
			
			
			CPT pDiffAttributionCPT = CA_HazardChangeFigures.getAttributionCPT(modCork, 10d, 50d);
			CPT diffAttributionCPT = CA_HazardChangeFigures.getAttributionCPT(modBroc, 0.05, 0.2d);
			
			GriddedGeoDataSet pDiffAttribution = new GriddedGeoDataSet(gridReg, false);
			GriddedGeoDataSet diffAttribution = new GriddedGeoDataSet(gridReg, false);
			
			int totFaultAttributed = 0;
			int totGridAttributed = 0;
			
			int faultAbove10percent = 0;
			int gridAbove10percent = 0;
			
			int faultAbove0p05 = 0;
			int gridAbove0p05 = 0;
			
			for (int i=0; i<diffAttribution.size(); i++) {
				double fullPDiff = Math.abs(fullHazardPDiff.get(i));
				double fullDiff = Math.abs(fullHazardDiff.get(i));
				
				double faultDiff = Math.abs(faultHazardDiff.get(i));
				double gridDiff = Math.abs(griddedHazardDiff.get(i));
				
				double sign = faultDiff > gridDiff ? 1d : -1d;
				
				pDiffAttribution.set(i, sign*fullPDiff);
				diffAttribution.set(i, sign*fullDiff);
				
				if (faultDiff > gridDiff) {
					totFaultAttributed++;
					if (fullPDiff > 10d)
						faultAbove10percent++;
					if (fullDiff > 0.05)
						faultAbove0p05++;
				} else {
					totGridAttributed++;
					if (fullPDiff > 10d)
						gridAbove10percent++;
					if (fullDiff > 0.05)
						gridAbove0p05++;
				}
			}
			
			String pDiffLabel = "Gridded       ←       |Hazard % Change|       →        Faults";
			String diffLabel  = "Gridded       ←      |Hazard Difference|      →        Faults";
			
			// full hazard difference
			lines.add("### Fault vs Gridded, Hazard Change Attribution");
			lines.add(topLink); lines.add("");
			
			lines.add("This section shows hazard changes, colored by their primary contribution (fault or gridded) "
					+ "rather than their sign. Small changes are masked out.");
			lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine("Ratio Attribution", "Difference Attribution");
			
			mapMaker.plotXYZData(pDiffAttribution, pDiffAttributionCPT, pDiffLabel);
			mapMaker.plot(resourcesDir, "ratio_attribution", "Hazard Change Attribution");
			table.addColumn("![Map]("+resourcesDir.getName()+"/ratio_attribution.png)");
			
			mapMaker.plotXYZData(diffAttribution, diffAttributionCPT, diffLabel);
			mapMaker.plot(resourcesDir, "diff_attribution", "Hazard Change Attribution");
			table.addColumn("![Map]("+resourcesDir.getName()+"/diff_attribution.png)");
			
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine("", "Attributed to Faults", "Attributed to Gridded", "Sum");
			
			double numLocs = gridReg.getNodeCount();
			
			table.initNewLine().addColumn("Full Map");
			table.addColumn(pDF.format((double)totFaultAttributed/numLocs));
			table.addColumn(pDF.format((double)totGridAttributed/numLocs));
			table.addColumn(pDF.format(1d));
			table.finalizeLine();
			
			table.initNewLine().addColumn("> 10% Change");
			table.addColumn(pDF.format((double)faultAbove10percent/numLocs));
			table.addColumn(pDF.format((double)gridAbove10percent/numLocs));
			table.addColumn(pDF.format((double)(faultAbove10percent + gridAbove10percent)/numLocs));
			table.finalizeLine();
			
			table.initNewLine().addColumn("> 0.05g Change");
			table.addColumn(pDF.format((double)faultAbove0p05/numLocs));
			table.addColumn(pDF.format((double)gridAbove0p05/numLocs));
			table.addColumn(pDF.format((double)(faultAbove0p05 + gridAbove0p05)/numLocs));
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static GriddedGeoDataSet asLog10(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		xyz.log10();
		return xyz;
	}
	
	private static class ChangeStats {
		private int totalNum = 0;
		
		// PDiff
		private int within1Percent = 0;
		private int within5Percent = 0;
		private int within10Percent = 0;
		private double maxPDiff = Double.NEGATIVE_INFINITY;
		private double minPDiff = Double.POSITIVE_INFINITY;
		private double sumPDiff = 0d;
		
		// Diff
		private int within0p05G = 0;
		private int within0p1G = 0;
		private int within0p2G = 0;
		private double maxDiff = Double.NEGATIVE_INFINITY;
		private double minDiff = Double.POSITIVE_INFINITY;
		private double sumDiff = 0d;
		
		public void addValue(double newVal, double refVal) {
			if (!Double.isFinite(newVal) || !Double.isFinite(refVal))
				return;
			
			totalNum++;
			
			double pDiff = 100d*(newVal - refVal)/refVal;
			double diff = newVal - refVal;
			
			if (Math.abs(pDiff) < 1d)
				within1Percent++;
			if (Math.abs(pDiff) < 5d)
				within5Percent++;
			if (Math.abs(pDiff) < 10d)
				within10Percent++;
			maxPDiff = Math.max(maxPDiff, pDiff);
			minPDiff = Math.min(minPDiff, pDiff);
			sumPDiff += pDiff;
			
			if (Math.abs(diff) < 0.05d)
				within0p05G++;
			if (Math.abs(diff) < 0.1d)
				within0p1G++;
			if (Math.abs(diff) < 0.2d)
				within0p2G++;
			maxDiff = Math.max(maxDiff, diff);
			minDiff = Math.min(minDiff, diff);
			sumDiff += diff;
		}
		
		public List<String> pDiffTableHeader() {
			List<String> header = new ArrayList<>(List.of("Within 1%", "Within 5%", "Within 10%", "% Range", "Average % Change"));
			for (int i=0; i<header.size(); i++)
				header.set(i, MarkdownUtils.boldCentered(header.get(i)));
			return header;
		}
		
		public List<String> diffTableHeader() {
			List<String> header = new ArrayList<>(List.of("Within 0.05g", "Within 0.1g", "Within 0.2g", "Diff Range (g)", "Average Diff (g)"));
			for (int i=0; i<header.size(); i++)
				header.set(i, MarkdownUtils.boldCentered(header.get(i)));
			return header;
		}
		
		public List<String> pDiffTableLine() {
			return List.of(
					pDF.format((double)within1Percent/(double)totalNum),
					pDF.format((double)within5Percent/(double)totalNum),
					pDF.format((double)within10Percent/(double)totalNum),
					"["+twoDigits.format(minPDiff)+"%, "+twoDigits.format(maxPDiff)+"%]",
					twoDigits.format(sumPDiff/(double)totalNum)+"%"
					);
		}
		
		public List<String> diffTableLine() {
			return List.of(
					pDF.format((double)within0p05G/(double)totalNum),
					pDF.format((double)within0p1G/(double)totalNum),
					pDF.format((double)within0p2G/(double)totalNum),
					"["+twoDigits.format(minDiff)+", "+twoDigits.format(maxDiff)+"]",
					twoDigits.format(sumDiff/(double)totalNum)
					);
		}
	}

}
