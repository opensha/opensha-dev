package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.CubedGriddedRegion;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_FaultCubeAssociations;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.SpecialCases;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import scratch.kevin.nshm23.SimpleSmoothHazardMapCalc;
import scratch.kevin.nshm23.MomentRateCompNSHM18;
import scratch.kevin.nshm23.SingleSiteHazardAndDataComparisonPageGen;
import scratch.kevin.nshm23.SingleSiteHazardAndDataComparisonPageGen.RegionalParticipationResult;

public class WUS_HazardChangePageGen {

	public static void main(String[] args) throws IOException {
		
//		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		ReturnPeriods rp = ReturnPeriods.TEN_IN_50;
//		double period = 0d;
		double period = 1d;
		
		String entryName, wrapperEntryName, hazLabel, dirPrefix;
		if (period == 0d) {
			entryName = "mean_map_pga_"+rp.name()+".txt";
			wrapperEntryName = "map_pga_"+rp.name()+".txt";
			hazLabel = "PGA, "+rp.label;
			dirPrefix = "pga";
		} else {
			entryName = "mean_map_"+(float)period+"s_"+rp.name()+".txt";
			wrapperEntryName = "map_"+(float)period+"s_"+rp.name()+".txt";
			hazLabel = oDF.format(period)+"s SA, "+rp.label;
			dirPrefix = "sa_"+oDF.format(period)+"s";
		}
		
		dirPrefix += "_"+rp.name().toLowerCase().replace("two", "2").replace("ten", "10").replace("_", "");
		
		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;
		ScalarIMR gmpe = gmpeRef.get();
		String gmpeName = "ASK (2014)";
		boolean hasSubduction = false;
		
		// load hazard maps
		File nshm23HazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		GriddedGeoDataSet nshm23Hazard = CA_HazardChangeFigures.loadXYZ(nshm23HazardFile, entryName);
		
//		File outputDir = new File(nshm23HazardFile.getParentFile(), "hazard_comparisons_nshm18_"+dirPrefix);
		File outputDir = new File(nshm23HazardFile.getParent().replace("-ba_only", ""), "hazard_comparisons_nshm18_"+dirPrefix);
		
		File nshm23GridHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_only_0.1deg.zip");
		GriddedGeoDataSet nshm23GridHazard = CA_HazardChangeFigures.loadXYZ(nshm23GridHazardFile, entryName);
		Preconditions.checkState(nshm23GridHazard.size() == nshm23Hazard.size());
//		GriddedGeoDataSet nshm23GridHazard = null;
		
		FaultSystemSolution nshm23BASol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		File nshm18_23gridHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_27-nshm18-grid_src_from_23-hazard-ask2014-0.1deg-noSub/results_hazard.zip");
		GriddedGeoDataSet nshm18_23gridHazard = CA_HazardChangeFigures.loadXYZ(nshm18_23gridHazardFile, wrapperEntryName);
		Preconditions.checkState(nshm18_23gridHazard.size() == nshm23Hazard.size());
		
		File nshm18HazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_27-nshm18-hazard-ask2014-0.1deg-noSub/results_hazard.zip");
		GriddedGeoDataSet nshm18Hazard = CA_HazardChangeFigures.loadXYZ(nshm18HazardFile, wrapperEntryName);
		Preconditions.checkState(nshm18Hazard.size() == nshm23Hazard.size());
		
		File nshm18GridHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_27-nshm18-hazard-ask2014-0.1deg-noSub-griddedOnly/results_hazard.zip");
		GriddedGeoDataSet nshm18GridHazard = CA_HazardChangeFigures.loadXYZ(nshm18GridHazardFile, wrapperEntryName);
		Preconditions.checkState(nshm18GridHazard.size() == nshm23Hazard.size());
//		GriddedGeoDataSet nshm18GridHazard = null;
		
		boolean doNSHM18Ingredients = true;
		
		GriddedGeoDataSet nshm18IngredNewScaleHazard = null;
		GriddedGeoDataSet nshm18IngredNewScaleClassicHazard = null;
		GriddedGeoDataSet nshm18IngredWCHazard = null;
		GriddedGeoDataSet nshm18IngredWCClassicHazard = null;
		if (doNSHM18Ingredients) {
			File nshm18IngredNewScaleHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
					+ "2023_01_25-nshm18_branches-new_scale-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED"
					+ "-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded/results_hazard_include_0.1deg.zip");
			nshm18IngredNewScaleHazard = CA_HazardChangeFigures.loadXYZ(nshm18IngredNewScaleHazardFile, entryName);
			Preconditions.checkState(nshm18IngredNewScaleHazard.size() == nshm23Hazard.size());
			
			File nshm18IngredNewScaleClassicHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
					+ "2023_01_25-nshm18_branches-new_scale-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED"
					+ "-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded-classic_only/results_hazard_include_0.1deg.zip");
			nshm18IngredNewScaleClassicHazard = CA_HazardChangeFigures.loadXYZ(nshm18IngredNewScaleClassicHazardFile, entryName);
			Preconditions.checkState(nshm18IngredNewScaleClassicHazard.size() == nshm23Hazard.size());
			
			File nshm18IngredWCHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
					+ "2023_01_26-nshm18_branches-wc_94-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED"
					+ "-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded/results_hazard_include_0.1deg.zip");
			nshm18IngredWCHazard = CA_HazardChangeFigures.loadXYZ(nshm18IngredWCHazardFile, entryName);
			Preconditions.checkState(nshm18IngredWCHazard.size() == nshm23Hazard.size());
			
			File nshm18IngredWCClassicHazardFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
					+ "2023_01_26-nshm18_branches-wc_94-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED"
					+ "-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded-classic_only/results_hazard_include_0.1deg.zip");
			nshm18IngredWCClassicHazard = CA_HazardChangeFigures.loadXYZ(nshm18IngredWCClassicHazardFile, entryName);
			Preconditions.checkState(nshm18IngredWCClassicHazard.size() == nshm23Hazard.size());
		}
		
		boolean recalcHazard = false;
		
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
		
		double[] momentThresholdPercents = { Double.POSITIVE_INFINITY, 100d, 50d, 20d };
		
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
//		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		CPT pDiffCPT = CA_HazardChangeFigures.getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance(), 10d, 50d);
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
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		String cascadiaSentence = hasSubduction ? "" : " The Cascadia subduction zone implementation is excluded in these comparisons.";
		
		lines.add("# Hazard Comparisons with NSHM18");
		lines.add("");
		lines.add("This report compares hazard between the NSHM23 and NSHM18 WUS ERFs, attempting to tease apart the "
				+ "contributions from the fault and gridded seismicity models to overall hazard changes.");
		lines.add("");
		lines.add("All comparisons are done with a branch-averaged mean ERF model (as is used in the final NSHM), and a "
				+ "single ground motion model, "+gmpeName+", using uniform and default site conditions. All maps are for "
				+ hazLabel+"."+cascadiaSentence);
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Fault-Based Hazard Comparisons");
		lines.add(topLink); lines.add("");
		
		lines.add("This section compares NSHM23 with NSHM18, but holding the gridded seismicity model constant in "
				+ "order to focus on fault-based changes (due both to ingredient and methodological changes) in areas "
				+ "where hazard is dominated by modeled faults.");
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
		if (MomentRateCompNSHM18.LINEAR_RAMP) {
			lines.add("To smooth out minor geometry changes between the two models, moment is distributed around "
					+ "faults using a linear ramp in (3-D) up to "
					+ oDF.format(NSHM23_SingleRegionGridSourceProvider.DEFAULT_MAX_FAULT_NUCL_DIST)+" km; this is the "
					+ "same ramp used to carve out gridded seismicity near faults when constructing the NSHM23 gridded "
					+ "seismicity model.");
			lines.add("");
		}
		
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
		mapMaker.plot(resourcesDir, "moment_ratio", "Fault Moment Change, Ratio");
		table.addColumn("![Map]("+resourcesDir.getName()+"/moment_ratio.png)");
		
		mapMaker.plotXYZData(momentDiff, momentDiffCPT, "NSHM23 Fault Moment) - (NSHM18 Fault Moment)");
		mapMaker.plot(resourcesDir, "moment_diff", "Fault Moment Change, Difference");
		table.addColumn("![Map]("+resourcesDir.getName()+"/moment_diff.png)");
		
		table.finalizeLine();
		
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("### Fault-Based Hazard Comparisons, Moment Chage Masked");
		lines.add(topLink); lines.add("");
		
		lines.add("This section masks hazard changes, only showing them near faults where NSHM23 fault-moment is within "
				+ "a given factor of NSHM18 fault-moment. Hazard changes in areas that remain may still be largely due "
				+ "to moment changes, but this helps to narrow the search space for identification of areas where "
				+ "methodological changes may dominate hazard changes.");
		lines.add("");
		lines.add("Various thresholds are used, starting with one that only includes areas with faults in both models, "
				+ "and then the threshold decreases to only show areas where moment is similar in the two models. Areas "
				+ "where the sign of the hazard change differs from the sign of moment change are shown regardless of "
				+ "the threshold value.");
		lines.add("");
		
		TableBuilder threshMapTable = MarkdownUtils.tableBuilder();
		threshMapTable.addLine("Ratio", "Difference");
		
		List<String> threshLabels = new ArrayList<>();
		List<ChangeStats> threshStats = new ArrayList<>();
		
		for (double pThresh : momentThresholdPercents) {
			String threshLabel, threshPrefix;
			if (Double.isFinite(pThresh)) {
				threshLabel = "<"+oDF.format(pThresh)+"% Moment Change Threshold";
				threshPrefix = "mom_thresh_p"+oDF.format(pThresh);
			} else {
				threshLabel = "Areas w/ Moment in Both Models";
				threshPrefix = "mom_thresh_inf";
			}
			
			threshLabels.add(threshLabel);
			
			ChangeStats stats = new ChangeStats();
			threshStats.add(stats);
			GriddedGeoDataSet maskedPDiff = new GriddedGeoDataSet(gridReg, false);
			GriddedGeoDataSet maskedDiff = new GriddedGeoDataSet(gridReg, false);
			for (int i=0; i<nshm23Hazard.size(); i++) {
				double mom23 = momentRates23.get(i);
				double mom18 = momentRates18.get(i);
				boolean skip = Double.isNaN(mom23) || Double.isNaN(mom18) || mom23 == 0d || mom18 == 0d;
				if (!skip) {
					// check the threshold
					double pDiff = Math.abs(100d*(mom23 - mom18)/mom18);
					if (pDiff > pThresh) {
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
		
//		lines.add("### Fault-Based Hazard and Moment Change Ratios");
//		lines.add(topLink); lines.add("");
//		
//		CSVFile<String> csv = new CSVFile<>(false);
//		csv.addLine("Grid Cell Index", "Latitude", "Longitude",
//				"NSHM23 Fault Moment (N-m)", "NSHM18 Fault Moment (N-m)", "Moment Change Ratio",
//				"NSHM23 Fault Hazard, "+hazLabel+" (g)", "NSHM18 Fault Hazard, "+hazLabel+" (g)", "Fault Hazard Ratio",
//				"Fault(s) Associated");
//		
//		GriddedGeoDataSet hazMomRatioMap = new GriddedGeoDataSet(gridReg, false);
//		DefaultXY_DataSet hazMomScatter = new DefaultXY_DataSet();
//		
//		List<? extends FaultSection> nshm23Faults = NSHM23_DeformationModels.AVERAGE.build(NSHM23_FaultModels.NSHM23_v2);
//		FaultGridAssociations assoc23 = new NSHM23_FaultCubeAssociations(
//				nshm23Faults, cgr, NSHM23_SingleRegionGridSourceProvider.DEFAULT_MAX_FAULT_NUCL_DIST);
//		
//		for (int i=0; i<hazMomRatioMap.size(); i++) {
//			double mom23 = momentRates23.get(i);
//			double mom18 = momentRates18.get(i);
//			
//			if (mom23 > 0d && mom18 > 0d) {
//				List<String> line = new ArrayList<>();
//				line.add(i+"");
//				Location loc = gridReg.getLocation(i);
//				line.add((float)loc.getLatitude()+"");
//				line.add((float)loc.getLongitude()+"");
//				line.add((float)mom23+"");
//				line.add((float)mom18+"");
//				double momRatio = mom23/mom18;
//				line.add((float)momRatio+"");
//				double haz23 = nshm23Hazard.get(i);
//				double haz18 = nshm18_23gridHazard.get(i);
//				double hazRatio = haz23/haz18;
//				line.add((float)haz23+"");
//				line.add((float)haz18+"");
//				line.add((float)hazRatio+"");
//				
//				Map<Integer, Double> assocSects = assoc23.getSectionFracsOnNode(i);
//				HashSet<String> faultNames = new HashSet<>();
//				if (assocSects != null) {
//					for (Integer id : assocSects.keySet()) {
//						FaultSection sect = nshm23Faults.get(id);
//						if (!faultNames.contains(sect.getParentSectionName())) {
//							faultNames.add(sect.getParentSectionName());
//							line.add(sect.getParentSectionName());
//						}
//					}
//				}
//				csv.addLine(line);
//				
//				hazMomScatter.set(momRatio, hazRatio);
//				hazMomRatioMap.set(i, Math.log10(hazRatio)/Math.log10(momRatio));
//			} else {
//				hazMomRatioMap.set(i, Double.NaN);
//			}
//		}
//		
//		csv.writeToFile(new File(resourcesDir, "fault_moment_hazard_comparison.csv"));
//		
//		lines.add("This section shows the relationship between near-fault moment changes and hazard changes.");
//		lines.add("");
//		lines.add("The map shows a ratio of ratios in an attempt to identify areas where the Log10 hazard change ratio "
//				+ "differes in sign or scale from the Log10 moment change ratio:");
//		lines.add("");
//		lines.add("* Dark red areas: places where hazard change was greater than moment change, but they were in the "
//				+ "same direction (i.e., both increasing, or both decreasing)");
//		lines.add("* Dark blue areas: places where hazard change was greater than moment change, but they were in the "
//				+ "different directions (i.e., one increased and the other decreased)");
//		lines.add("* Lightly colored areas: places where hazard change was smaller than moment change");
//		lines.add("");
//		lines.add("");
//		lines.add("Download Data: [fault_moment_hazard_comparison.csv]("+resourcesDir.getName()
//			+"/fault_moment_hazard_comparison.csv)");
//		lines.add("");
//		
//		List<XY_DataSet> funcs = new ArrayList<>();
//		List<PlotCurveCharacterstics> chars = new ArrayList<>();
//		
//		Range range = new Range(1e-1, 1e1);
//		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
//		oneToOne.set(range.getLowerBound(), range.getLowerBound());
//		oneToOne.set(range.getUpperBound(), range.getUpperBound());
//		
//		funcs.add(oneToOne);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
//		
//		DefaultXY_DataSet boundedScatter = new DefaultXY_DataSet();
//		for (Point2D pt : hazMomScatter)
//			boundedScatter.set(Math.max(range.getLowerBound(), Math.min(range.getUpperBound(), pt.getX())),
//					Math.max(range.getLowerBound(), Math.min(range.getUpperBound(), pt.getY())));
//		
////		funcs.add(hazMomScatter);
//		funcs.add(boundedScatter);
//		chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, new Color(0, 0, 0, 60)));
//		
//		PlotSpec scatterSpec = new PlotSpec(funcs, chars, "Fault Hazard and Moment Comparison", "Moment Ratio", "Hazard Ratio");
//		
//		HeadlessGraphPanel gp = PlotUtils.initHeadless();
//		gp.drawGraphPanel(scatterSpec, true, true, range, range);
//		
//		table = MarkdownUtils.tableBuilder();
//		
//		table.addLine("Scatter Plot", "Map View");
//		
//		table.initNewLine();
//		PlotUtils.writePlots(resourcesDir, "haz_mom_ratios_scatter", gp, 800, false, true, false, false);
//		table.addColumn("![scatter plot]("+resourcesDir.getName()+"/haz_mom_ratios_scatter.png)");
//		
//		CPT hazMomRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-2d, 2d);
//		hazMomRatioCPT.setNanColor(transparent);
//		mapMaker.plotXYZData(hazMomRatioMap, hazMomRatioCPT, "Log10(Hazard Change Ratio) / Log10(Moment Change Ratio)");
//		
//		mapMaker.plot(resourcesDir, "haz_mom_ratios_map", "Fault Hazard and Moment Comparison");
//		table.addColumn("![map plot]("+resourcesDir.getName()+"/haz_mom_ratios_map.png)");
//		table.finalizeLine();
//		
//		lines.addAll(table.build());
//		lines.add("");
		
		lines.add("### NSHM18 Special Cases");
		lines.add(topLink); lines.add("");
		
		lines.add("NSHM18 includes a number of special cases that are not carried forward. These also affect hazard "
				+ "comparisons beyond the moment changes shown above, and are plotted below.");
		lines.add("");
		
		SpecialCases.plotSpecialCases(resourcesDir, "nshm18_special_cases", gridReg);
		lines.add("![Special Cases]("+resourcesDir.getName()+"/nshm18_special_cases.png)");
		lines.add("");
		lines.add(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON",
				resourcesDir.getName()+"/nshm18_special_cases.geojson"));
		lines.add("");
		
		if (doNSHM18Ingredients) {
			lines.add("### NSHM18 Ingredient Runs vs NSHM18");
			lines.add(topLink); lines.add("");
			
			lines.add("We also ran inversions using NSHM18 deformation models as input, but using updated methodologies. "
					+ "This helps to isolate regions where methodological changes (the inversion, as well as elmination "
					+ "of the special cases above) dominate hazard change. All comparisons hold the NSHM23 gridded "
					+ "seismicity model constant.");
			lines.add("");
			lines.add("We ran this test both using updated scaling relationships, as well as using Wells and "
					+ "Coppersmith (94); for the latter, results are masked in California as UCERF3 used "
					+ "their own scaling relationships. We also show results using the classic segmentation branches "
					+ "only, which should be most similar to NSHM18 outside of California.");
			lines.add("");
			
			GriddedGeoDataSet compMap = nshm18_23gridHazard;
			List<GriddedGeoDataSet> numerators = new ArrayList<>();
			List<String> longNames = new ArrayList<>();
			List<String> shortNames = new ArrayList<>();
			List<String> prefixes = new ArrayList<>();
			List<Boolean> maskCAs = new ArrayList<>();
			
			XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
			List<Region> caRegions = new ArrayList<>();
			for (int i=0; i<caOutlines.length; i++) {
				LocationList outline = new LocationList();
				for (Point2D pt : caOutlines[i])
					outline.add(new Location(pt.getY(), pt.getX()));
				caRegions.add(new Region(outline, BorderType.MERCATOR_LINEAR));
			}
			// now add some rectangles to mask out the offshort components
			Location lowLeftMask = new Location(30, -126);
			caRegions.add(new Region(lowLeftMask, new Location(32.718562, -114.719390))); // AZ/CA lower corner
			caRegions.add(new Region(lowLeftMask, new Location(36.15, -116.89)));
			caRegions.add(new Region(lowLeftMask, new Location(42, -120))); // OR/CA/NV corner
			
			numerators.add(nshm18IngredNewScaleHazard);
			longNames.add("NSHM18 Ingredients, NSHM23 Methodology (new scaling)");
			shortNames.add("Fault Methodologies");
			prefixes.add("nshm18_ingred");
			maskCAs.add(false);
			
			numerators.add(nshm18IngredWCHazard);
			longNames.add("NSHM18 Ingredients, NSHM23 Methodology (W-C '94)");
			shortNames.add("Fault Methodologies w/ WC '94");
			prefixes.add("nshm18_ingred_wc94");
			maskCAs.add(true);
			
			numerators.add(nshm18IngredNewScaleClassicHazard);
			longNames.add("NSHM18 Ingredients, NSHM23 Methodology (new scaling, classic only)");
			shortNames.add("Fault Methodologies (Classic Only)");
			prefixes.add("nshm18_ingred_classic");
			maskCAs.add(false);
			
			numerators.add(nshm18IngredWCClassicHazard);
			longNames.add("NSHM18 Ingredients, NSHM23 Methodology (W-C '94, classic only");
			shortNames.add("Fault Methodologies (Classic Only) w/ WC '94");
			prefixes.add("nshm18_ingred_wc94_classic");
			maskCAs.add(true);
			
			table = MarkdownUtils.tableBuilder();
			
			for (int i=0; i<numerators.size(); i++) {
				table.addLine(MarkdownUtils.boldCentered(longNames.get(i)), "");

//				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				GriddedGeoDataSet numerator = numerators.get(i);
				String shortName = shortNames.get(i);
				String prefix = prefixes.get(i);
				
				if (maskCAs.get(i)) {
					// apply mask
					numerator = numerator.copy();
					for (int j=0; j<numerator.size(); j++) {
						for (Region caReg : caRegions) {
							if (caReg.contains(numerator.getLocation(j))) {
								numerator.set(j, Double.NaN);
								break;
							}
						}
					}
				}
				
				table.initNewLine();
				
				mapMaker18.plotXYZData(mapPDiff(numerator, compMap), pDiffCPT, shortName+", % Change, "+hazLabel);
				mapMaker18.plot(resourcesDir, prefix+"_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_pDiff.png)");
				
				mapMaker18.plotXYZData(mapDiff(numerator, compMap), diffCPT, shortName+", Difference, "+hazLabel+" (g)");
				mapMaker18.plot(resourcesDir, prefix+"_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_diff.png)");
				
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		ChangeStats griddedStats = null;
		if (gridded) {
			// gridded differences
			lines.add("## Gridded Seismicity Model Hazard Changes");
			lines.add(topLink); lines.add("");
			
			lines.add("These plots compare the gridded seismicity components of the models. For the ratio and "
					+ "difference maps, the contributions of NSHM23 fault sources are held constant (added to both "
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
			
//			// now mask outside of faults
//			lines.add("### Gridded Seismicity Model Hazard Changes, Away From Faults");
//			lines.add(topLink); lines.add("");
//			
//			lines.add("Here, we mask gridded seismicity hazard changes only showing them away from faults in either model.");
//			lines.add("");
//			
//			ChangeStats stats = new ChangeStats();
//			GriddedGeoDataSet maskedPDiff = new GriddedGeoDataSet(gridReg, false);
//			GriddedGeoDataSet maskedDiff = new GriddedGeoDataSet(gridReg, false);
//			for (int i=0; i<nshm23Hazard.size(); i++) {
//				double mo23 = momentRates23.get(i);
//				double mo18 = momentRates18.get(i);
//				if (mo23 > 0 || mo18 > 0) {
//					maskedPDiff.set(i, Double.NaN);
//					maskedDiff.set(i, Double.NaN);
//				} else {
//					double haz23 = nshm23Hazard.get(i);
//					double grid23 = nshm23GridHazard.get(i);
//					double fault23 = Math.max(0, haz23 - grid23);
//					double grid18 = nshm18GridHazard.get(i);
//					stats.addValue(haz23, grid18+fault23);
//					maskedPDiff.set(i, griddedHazardPDiff.get(i));
//					maskedDiff.set(i, griddedHazardDiff.get(i));
//				}
//			}
//			
//			table = MarkdownUtils.tableBuilder();
//			table.addLine("Ratio", "Difference");
//			
//			table.initNewLine();
//			
//			mapMaker.plotXYZData(maskedPDiff, pDiffCPT, "NSHM23 vs NSHM18, Gridded Only, % Change, "+hazLabel);
//			mapMaker.plot(resourcesDir, "gridded_hazard_pDiff_away_from_faults", "Gridded Hazard Away From Faults");
//			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_pDiff_away_from_faults.png)");
//			
//			mapMaker.plotXYZData(maskedDiff, diffCPT, "NSHM23 - NSHM18, Gridded Only, "+hazLabel+" (g)");
//			mapMaker.plot(resourcesDir, "gridded_hazard_diff_away_from_faults", "Gridded Hazard Away From Faults");
//			table.addColumn("![Map]("+resourcesDir.getName()+"/gridded_hazard_diff_away_from_faults.png)");
//			
//			table.finalizeLine();
//			
//			lines.addAll(table.build());
//			lines.add("");
//			
//			withinLabel = MarkdownUtils.boldCentered("% Locs Away From Faults");
//			table = MarkdownUtils.tableBuilder();
//			table.initNewLine();
//			table.addColumn(withinLabel);
//			table.addColumns(stats.pDiffTableHeader());
//			table.finalizeLine().initNewLine();
//			table.addColumn(pDF.format((double)stats.totalNum/(double)maskedPDiff.size()));
//			table.addColumns(stats.pDiffTableLine());
//			table.finalizeLine();
//			table.initNewLine();
//			table.addColumn(withinLabel);
//			table.addColumns(stats.diffTableHeader());
//			table.finalizeLine().initNewLine();
//			table.addColumn(pDF.format((double)stats.totalNum/(double)maskedPDiff.size()));
//			table.addColumns(stats.diffTableLine());
//			table.finalizeLine();
//			lines.add("__Gridded Hazard Comparison Statistics Away From Faults:__");
//			lines.add("");
//			lines.addAll(table.build());
//			lines.add("");
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
			
			
			CPT pDiffAttributionCPT = CA_HazardChangeFigures.getCenterMaskedCPT(modCork, 10d, 50d);
			CPT diffAttributionCPT = CA_HazardChangeFigures.getCenterMaskedCPT(modBroc, 0.05, 0.2d);
			
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
			
			String pDiffLabel = "Gridded      ←        |Hazard % Change|        →       Faults";
			String diffLabel  = "Gridded      ←     |Hazard Difference (g)|     →       Faults";
			
			// full hazard difference
			lines.add("### Fault vs Gridded, Hazard Change Attribution");
			lines.add(topLink); lines.add("");
			
			lines.add("This section shows hazard changes, colored by the primary contributor to hazard change (faults "
					+ "or gridded seismicity) rather than the sign of that change. Small changes are masked out.");
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
		
		double[] magThresholds = {5d, 6.5d, 7.5};
		
		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.2.0");
		boolean subduction = false;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW, TectonicRegionType.STABLE_SHALLOW);
		if (subduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}

		HazardModel model = HazardModel.load(erfPath);
		NshmErf faultERF = new NshmErf(model, trts, IncludeBackgroundOption.EXCLUDE);
		System.out.println("NSHM Fault ERF size: " + faultERF.getNumSources());
		faultERF.getTimeSpan().setDuration(1.0);
		faultERF.updateForecast();
		NshmErf gridERF = new NshmErf(model, trts, IncludeBackgroundOption.ONLY);
		System.out.println("NSHM Grid ERF size: " + gridERF.getNumSources());
		gridERF.getTimeSpan().setDuration(1.0);
		gridERF.updateForecast();
		
		RegionalParticipationResult nshm18FaultPartic = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				faultERF, gridReg, magThresholds, null);
		RegionalParticipationResult nshm18GridPartic = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				gridERF, gridReg, magThresholds, null);

		GriddedGeoDataSet[] nshm18ParticRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			nshm18ParticRates[m] = sumMap(nshm18FaultPartic.particRateMaps[m], nshm18GridPartic.particRateMaps[m]);
		GriddedGeoDataSet[] nshm18NuclRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			nshm18NuclRates[m] = sumMap(nshm18FaultPartic.nuclRateMaps[m], nshm18GridPartic.nuclRateMaps[m]);
		GriddedGeoDataSet nshm18MomentRate = sumMap(nshm18FaultPartic.momentRateMap, nshm18GridPartic.momentRateMap);
		
		System.out.println("Building NSHM23 rate maps");
		RegionalParticipationResult nshm23FaultPartic = SingleSiteHazardAndDataComparisonPageGen.calcFSSFaultPartic(
				nshm23BASol, gridReg, magThresholds, null);
		RegionalParticipationResult nshm23GridPartic = SingleSiteHazardAndDataComparisonPageGen.calcFSSGriddedPartic(
				nshm23BASol, gridReg, magThresholds, null);
		
		GriddedGeoDataSet[] nshm23Rates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			nshm23Rates[m] = sumMap(nshm23FaultPartic.particRateMaps[m], nshm23GridPartic.particRateMaps[m]);
		GriddedGeoDataSet nshm23MomentRate = sumMap(nshm23FaultPartic.momentRateMap, nshm23GridPartic.momentRateMap);
		
		GriddedGeoDataSet[] ratePDiffs = new GriddedGeoDataSet[magThresholds.length];
		GriddedGeoDataSet[] rateDiffs = new GriddedGeoDataSet[magThresholds.length];
		
		for (int m=0; m<magThresholds.length; m++) {
			ratePDiffs[m] = new GriddedGeoDataSet(gridReg, false);
			rateDiffs[m] = new GriddedGeoDataSet(gridReg, false);
		}
		
		GriddedGeoDataSet erfMomentPDiff = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet erfMomentDiff = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			for (int m=0; m<magThresholds.length; m++) {
				ratePDiffs[m].set(i, 100d*(nshm23Rates[m].get(i) - nshm18ParticRates[m].get(i))/nshm18ParticRates[m].get(i));
				rateDiffs[m].set(i, nshm23Rates[m].get(i) - nshm18ParticRates[m].get(i));
			}
			erfMomentPDiff.set(i, 100d*(nshm23MomentRate.get(i) - nshm18MomentRate.get(i))/nshm18MomentRate.get(i));
			erfMomentDiff.set(i, nshm23MomentRate.get(i) - nshm18MomentRate.get(i));
		}
		
		CPT moCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(10d, 18d);
		moCPT.setNanColor(transparent);
		
		CPT moRateDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1e16, 1e16);
		moRateDiffCPT.setNanColor(transparent);
		
		CPT pDiff100CPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-100d, 100d);
		pDiff100CPT.setNanColor(transparent);
		
		lines.add("## Full Model Rate and Moment Rate Comparisons");
		lines.add(topLink); lines.add("");
		
		lines.add("### Full Model Rate Comparisons");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		
		for (int m=0; m<magThresholds.length; m++) {
			String magPrefix = "m"+oDF.format(magThresholds[m]);
			String magLabel = "M≥"+oDF.format(magThresholds[m]);
			table.addLine(MarkdownUtils.boldCentered("NSHM23, "+magLabel), MarkdownUtils.boldCentered("NSHM18, "+magLabel));
			
			double maxRate = 0d;
			double[] allRateDiffs = new double[gridReg.getNodeCount()];
			for (int i=0; i<allRateDiffs.length; i++) {
				double rate23 = nshm23Rates[m].get(i);
				double rate18 = nshm18ParticRates[m].get(i);
				maxRate = Math.max(maxRate, Math.max(rate23, rate18));
				allRateDiffs[i] = Math.abs(rate23-rate18);
				
				if (rate23 == 0d)
					nshm23Rates[m].set(i, Double.NaN);
				if (rate18 == 0d)
					nshm18ParticRates[m].set(i, Double.NaN);
			}
			double diffStat = StatUtils.percentile(allRateDiffs, 95d);
			double logMax = Math.ceil(Math.log10(maxRate));
			
			CPT nuclCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMax-6d, logMax);
			nuclCPT.setNanColor(transparent);
			
//			double diffMax = Math.pow(10, (int)(Math.log10(diffStat)+0.5));
			double diffMax = Math.pow(10, logMax-2);
			
			CPT rateDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-diffMax, diffMax);
			rateDiffCPT.setNanColor(transparent);
			
			table.initNewLine();
			
			mapMaker.plotXYZData(asLog10(nshm23Rates[m]), nuclCPT, "Participation Rate, "+magLabel);
			mapMaker.plot(resourcesDir, "nshm23_rate_"+magPrefix, " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/nshm23_rate_"+magPrefix+".png)");
			
			mapMaker18.plotXYZData(asLog10(nshm18ParticRates[m]), nuclCPT, "Participation Rate, "+magLabel);
			mapMaker18.plot(resourcesDir, "nshm18_rate_"+magPrefix, " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/nshm18_rate_"+magPrefix+".png)");
			
			table.finalizeLine();
			
			table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
			
			table.initNewLine();
			
			mapMaker.plotXYZData(ratePDiffs[m], pDiff100CPT, "% Change in Participation Rate, "+magLabel);
			mapMaker.plot(resourcesDir, "rate_pDiff_"+magPrefix, " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/rate_pDiff_"+magPrefix+".png)");
			
			mapMaker.plotXYZData(rateDiffs[m], rateDiffCPT, "Participation Rate Difference, "+magLabel);
			mapMaker.plot(resourcesDir, "rate_diff_"+magPrefix, " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/rate_diff_"+magPrefix+".png)");
			
			table.finalizeLine();
		}
		
		lines.addAll(table.build());
		
		lines.add("### Full Model Moment Rate Comparison");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		
		table.addLine(MarkdownUtils.boldCentered("NSHM23"), MarkdownUtils.boldCentered("NSHM18"));
		
		table.initNewLine();
		
		mapMaker.plotXYZData(asLog10(nshm23MomentRate), moCPT, "Log10 Moment Rate (N-m)");
		mapMaker.plot(resourcesDir, "full_nshm23_moment_rate", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/full_nshm23_moment_rate.png)");
		
		mapMaker18.plotXYZData(asLog10(nshm18MomentRate), moCPT, "Log10 Moment Rate (N-m)");
		mapMaker18.plot(resourcesDir, "full_nshm18_moment_rate", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/full_nshm18_moment_rate.png)");
		
		table.finalizeLine();
		
		table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
		
		table.initNewLine();
		
		mapMaker.plotXYZData(erfMomentPDiff, pDiff100CPT, "% Change in Moment Rate");
		mapMaker.plot(resourcesDir, "full_moment_rate_pDiff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/full_moment_rate_pDiff.png)");
		
		mapMaker.plotXYZData(erfMomentDiff, moRateDiffCPT, "Moment Rate Difference (N-m)");
		mapMaker.plot(resourcesDir, "full_moment_rate_diff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/full_moment_rate_diff.png)");
		
		table.finalizeLine();
		
		lines.addAll(table.build());
		
		for (boolean momentBased : new boolean[] { true, false}) {
			String rateName;
			String ratePrefix;
			
			if (momentBased) {
				rateName = "Moment Rate";
				ratePrefix = "moment_rate";
			} else {
				rateName = "M"+oDF.format(magThresholds[0])+" Rate";
				ratePrefix = "m"+oDF.format(magThresholds[0])+"_rate";
			}
			
			lines.add("## Smoothed "+rateName+" Hazard vs Model Hazard Change Comparisons");
			lines.add(topLink); lines.add("");
			
			System.out.println("Distance-smoothing "+rateName+" map...");
			// always include fault and gridded, but only variying the one we care about at a time
			
			GriddedGeoDataSet input18Map;
			GriddedGeoDataSet input18FaultMap;
			GriddedGeoDataSet input18GridMap;
			GriddedGeoDataSet input23Map;
			
			if (momentBased) {
				input18Map = nshm18MomentRate;
				input18FaultMap = sumMap(nshm18FaultPartic.momentRateMap, nshm23GridPartic.momentRateMap);
				input18GridMap = sumMap(nshm23FaultPartic.momentRateMap, nshm18GridPartic.momentRateMap);
				input23Map = nshm23MomentRate;
			} else {
				input18Map = nshm18NuclRates[0];
				input18FaultMap = sumMap(nshm18FaultPartic.nuclRateMaps[0], nshm23GridPartic.nuclRateMaps[0]);
				input18GridMap = sumMap(nshm23FaultPartic.nuclRateMaps[0], nshm18GridPartic.nuclRateMaps[0]);
				input23Map = sumMap(nshm23FaultPartic.nuclRateMaps[0], nshm23GridPartic.nuclRateMaps[0]);
			}
			
			GriddedGeoDataSet nshm23SmoothedTotalHazard;
			GriddedGeoDataSet nshm23SmoothedFaultHazard;
			GriddedGeoDataSet nshm23SmoothedGridHazard;
			GriddedGeoDataSet nshm18SmoothedTotalHazard;
			GriddedGeoDataSet nshm18SmoothedFaultHazard;
			GriddedGeoDataSet nshm18SmoothedGridHazard;
			
			GriddedGeoDataSet[] inputSmoothMaps = {
					input23Map, input18Map, input18FaultMap, input18GridMap};
			GriddedGeoDataSet[] smoothedMaps = new GriddedGeoDataSet[inputSmoothMaps.length];
			
			File[] mapFiles = {
					new File(resourcesDir, "hazard_"+ratePrefix+"_23.xyz"),
					new File(resourcesDir, "hazard_"+ratePrefix+"_18.xyz"),
					new File(resourcesDir, "hazard_"+ratePrefix+"_18fault.xyz"),
					new File(resourcesDir, "hazard_"+ratePrefix+"_18grid.xyz")
			};
			
			if (!recalcHazard) {
				// see if they're cached
				for (int i=0; i<mapFiles.length; i++) {
					if (mapFiles[i].exists()) {
						ArbDiscrGeoDataSet xyz = ArbDiscrGeoDataSet.loadXYZFile(mapFiles[i].getAbsolutePath(), false);
						if (xyz.size() == gridReg.getNodeCount()) {
							smoothedMaps[i] = new GriddedGeoDataSet(gridReg);
							for (int j=0; j<xyz.size(); j++)
								smoothedMaps[i].set(j, xyz.get(j));
						}
					}
				}
			}
			
			if (momentBased) {
				boolean doMFD = true;
//				GutenbergRichterMagFreqDist mfdShape = new GutenbergRichterMagFreqDist(6.5d, 5, 0.25, 1e16, 1d);
				GutenbergRichterMagFreqDist mfdShape = new GutenbergRichterMagFreqDist(6d, 7, 0.25, 1e16, 1d);
				double refMag = 7d;
				gmpe.setParamDefaults();
				
				String magDescription;
				if (doMFD) {
					for (int i=0; i<inputSmoothMaps.length; i++)
						if (smoothedMaps[i] == null)
							smoothedMaps[i] = SimpleSmoothHazardMapCalc.hazMapMomentSmooth(inputSmoothMaps[i], gmpeRef, mfdShape, period, rp);
					magDescription = "a G-R b="+oDF.format(mfdShape.get_bValue())+" MFD with magnitudes in the range ["
							+oDF.format(mfdShape.getMinX())+","+oDF.format(mfdShape.getMaxX())+"] that satisfies each cell's "
							+ "moment rate.";
				} else {
					for (int i=0; i<inputSmoothMaps.length; i++)
						if (smoothedMaps[i] == null)
							smoothedMaps[i] = SimpleSmoothHazardMapCalc.hazMapMomentSmooth(inputSmoothMaps[i], gmpeRef, refMag, period, rp);
					magDescription = "fixed M="+oDF.format(refMag)+" ruptures that satisfy each cell's moment rate.";
				}
				lines.add("Here, we calculate simple hazard maps by placing point sources at every grid location with "
						+ magDescription+" We then compare hazard changes from this simple model to the actual model hazard "
						+ "change. Areas where simplified model hazard change is similar to full model hazard change are "
						+ "likely dominated by ingredient (moment rate) changes, and those that differ may be affected by "
						+ "methodological changes (although this simplified comparison is not definitive).");
				lines.add("");
				lines.add("Note that these comparisons use final model moment rate maps, as opposed to the prior moment "
						+ "comparisons that used deformation model moment rates directly. Thus, any slip rate misfits will be "
						+ "incorporated into these comparisons.");
				lines.add("");
			} else {
				GutenbergRichterMagFreqDist mfdShape = new GutenbergRichterMagFreqDist(magThresholds[0], 6, 0.5, 1e16, 1d);
				gmpe.setParamDefaults();
				
				String magDescription;
				for (int i=0; i<inputSmoothMaps.length; i++)
					if (smoothedMaps[i] == null)
						smoothedMaps[i] = SimpleSmoothHazardMapCalc.hazMapRateSmooth(inputSmoothMaps[i], gmpeRef, mfdShape, period, rp);
				magDescription = "a G-R b="+oDF.format(mfdShape.get_bValue())+" MFD with magnitudes in the range ["
						+oDF.format(mfdShape.getMinX())+","+oDF.format(mfdShape.getMaxX())+"] that matches the total "
						+ "nucleation M"+oDF.format(magThresholds[0])+" rate in each cell.";
				lines.add("Here, we calculate simple hazard maps by placing point sources at every grid location with "
						+ magDescription+" We then compare hazard changes from this simple model to the actual model hazard "
						+ "change. Areas where simplified model hazard change is similar to full model hazard change are "
						+ "likely dominated by ingredient (total rate) changes, and those that differ may be affected by "
						+ "methodological changes (although this simplified comparison is not definitive).");
				lines.add("");
			}
			
			for (int i=0; i<mapFiles.length; i++)
				ArbDiscrGeoDataSet.writeXYZFile(smoothedMaps[i], mapFiles[i]);
			
			nshm23SmoothedTotalHazard = smoothedMaps[0];
			nshm23SmoothedFaultHazard = smoothedMaps[0];
			nshm23SmoothedGridHazard = smoothedMaps[0];
			nshm18SmoothedTotalHazard = smoothedMaps[1];
			nshm18SmoothedFaultHazard = smoothedMaps[2];
			nshm18SmoothedGridHazard = smoothedMaps[3];
			
			for (int n=0; n<3; n++) {
				String label, prefix;
				GriddedGeoDataSet nshm23Smoothed, nshm18Smoothed;
				GriddedGeoDataSet hazardPDiff, hazardRatio;
				
				GriddedGeoDataSet origHaz18;
				
				// TODO haz from full curves? wrapperCurvesEntryName
				
				if (n == 0) {
					label = "Full Model";
					prefix = "smoothed_full_"+ratePrefix;
					nshm23Smoothed = nshm23SmoothedTotalHazard;
					nshm18Smoothed = nshm18SmoothedTotalHazard;
					hazardPDiff = fullHazardPDiff;
					hazardRatio = fullHazardRatio;
					origHaz18 = nshm18Hazard;
				} else if (n == 1) {
					label = "Fault";
					prefix = "smoothed_fault_"+ratePrefix;
					nshm23Smoothed = nshm23SmoothedFaultHazard;
					nshm18Smoothed = nshm18SmoothedFaultHazard;
					hazardPDiff = faultHazardPDiff;
					hazardRatio = faultHazardRatio;
					origHaz18 = nshm18_23gridHazard;
					if (!momentBased)
						continue;
				} else {
					label = "Gridded";
					prefix = "smoothed_gridded_"+ratePrefix;
					nshm23Smoothed = nshm23SmoothedGridHazard;
					nshm18Smoothed = nshm18SmoothedGridHazard;
					hazardPDiff = griddedHazardPDiff;
					hazardRatio = griddedHazardRatio;
					origHaz18 = nshm18GridHazard;
				}
				
				GriddedGeoDataSet smoothedMomentPDiff = new GriddedGeoDataSet(gridReg);
				GriddedGeoDataSet smoothedMomentRatio = new GriddedGeoDataSet(gridReg);
				GriddedGeoDataSet smoothedMomentDiff = new GriddedGeoDataSet(gridReg);
				GriddedGeoDataSet methodologyPDiff = new GriddedGeoDataSet(gridReg);
				for (int i=0; i<gridReg.getNodeCount(); i++) {
					double moment23 = nshm23Smoothed.get(i);
					double moment18 = nshm18Smoothed.get(i);
					smoothedMomentPDiff.set(i, 100d*(moment23-moment18)/moment18);
					smoothedMomentRatio.set(i, moment23/moment18);
					smoothedMomentDiff.set(i, moment23-moment18);
					
					// estimate the change from methodological sources only
					double haz18 = origHaz18.get(i);
					double haz23 = nshm23Hazard.get(i);
					// estimate of what things would be like in 18 if we had 23 moment
					double modModHaz18 = haz18*moment23/moment18;
					methodologyPDiff.set(i, 100d*(haz23-modModHaz18)/modModHaz18);
				}
				
				table = MarkdownUtils.tableBuilder();
				
				String moLabel = "Simplified "+label+" "+rateName+" Hazard";
				String moShortLabel = "Simplified "+rateName+" Hazard";
				
				table.addLine(moLabel+" Ratio", label+" Hazard Ratio");
				
				table.initNewLine();
				
				mapMaker.plotXYZData(smoothedMomentPDiff, pDiffCPT, "% Change in "+moLabel);
				mapMaker.plot(resourcesDir, prefix+"_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_pDiff.png)");
				
				mapMaker.plotXYZData(hazardPDiff, pDiffCPT, "% Change in "+label+" Hazard");
				mapMaker.plot(resourcesDir, prefix+"_hazard_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_hazard_pDiff.png)");
				
				table.finalizeLine();
				
				GriddedGeoDataSet hazMomRatioMap = new GriddedGeoDataSet(gridReg, false);
				GriddedGeoDataSet hazMomDiffMap = new GriddedGeoDataSet(gridReg, false);
				DefaultXY_DataSet hazMomScatter = new DefaultXY_DataSet();
				
				for (int i=0; i<hazMomRatioMap.size(); i++) {
					double momPDiff = smoothedMomentPDiff.get(i);
					double hazPDiff = hazardPDiff.get(i);
					hazMomScatter.set(momPDiff, hazPDiff);
					hazMomRatioMap.set(i, Math.log10(hazardRatio.get(i))/Math.log10(smoothedMomentRatio.get(i)));
					hazMomDiffMap.set(i, hazPDiff - momPDiff);
				}
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				Range range = new Range(-100d, 100d);
				DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
				oneToOne.set(range.getLowerBound(), range.getLowerBound());
				oneToOne.set(range.getUpperBound(), range.getUpperBound());
				
				DefaultXY_DataSet boundedScatter = new DefaultXY_DataSet();
				for (Point2D pt : hazMomScatter)
					boundedScatter.set(Math.max(range.getLowerBound(), Math.min(range.getUpperBound(), pt.getX())),
							Math.max(range.getLowerBound(), Math.min(range.getUpperBound(), pt.getY())));
				
//				funcs.add(hazMomScatter);
				funcs.add(boundedScatter);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, new Color(0, 0, 0, 60)));
				
				funcs.add(oneToOne);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
				
				PlotSpec scatterSpec = new PlotSpec(funcs, chars, label+" Hazard and "+rateName+" Change Comparison",
						moLabel+" % Change", label+" Hazard % Change");
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				gp.drawGraphPanel(scatterSpec, false, false, range, range);
				
				table.addLine(MarkdownUtils.boldCentered("Scatter Plot"), MarkdownUtils.boldCentered("Hazard Change Comparison"));
				
				table.initNewLine();
				PlotUtils.writePlots(resourcesDir, prefix+"_haz_ratios_scatter", gp, 800, false, true, false, false);
				table.addColumn("![scatter plot]("+resourcesDir.getName()+"/"+prefix+"_haz_ratios_scatter.png)");
				
				CPT smoothHazDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse();
				smoothHazDiffCPT = CA_HazardChangeFigures.getCenterMaskedCPT(smoothHazDiffCPT, 10, 50);
				smoothHazDiffCPT.setNanColor(transparent);
				mapMaker.plotXYZData(hazMomDiffMap, smoothHazDiffCPT, "Hazard % Change - "+moShortLabel+" % Change");
				mapMaker.plot(resourcesDir, prefix+"_haz_diffs_map", label+" Hazard and "+rateName+" Change Comparison");
				table.addColumn("![map plot]("+resourcesDir.getName()+"/"+prefix+"_haz_diffs_map.png)");
				
				table.finalizeLine();
				
				lines.addAll(table.build());
				lines.add("");
			}
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
	
	static GriddedGeoDataSet asLog10(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		xyz.log10();
		return xyz;
	}
	
	static GriddedGeoDataSet sumMap(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, map1.get(i)+map2.get(i));
		return ret;
	}
	
	static GriddedGeoDataSet mapPDiff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, 100d*(map1.get(i)-map2.get(i))/map2.get(i));
		return ret;
	}
	
	static GriddedGeoDataSet mapDiff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, map1.get(i)-map2.get(i));
		return ret;
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
