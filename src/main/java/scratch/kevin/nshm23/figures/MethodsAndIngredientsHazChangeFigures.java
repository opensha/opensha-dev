package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.function.Consumer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.title.Title;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;

import com.google.common.base.Preconditions;

public class MethodsAndIngredientsHazChangeFigures {
	
	public static void main(String[] args) throws IOException {
		doCA();
//		doWUS();
	}
	
	public static void doCA() throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/u3_haz_change_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FileWriter logFW = new FileWriter(new File(outputDir, "log.txt"));
		
		FaultSystemRupSet rupSetU3 = FaultSystemRupSet.load(new File(
				"/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip"));
		FaultSystemRupSet rupSet23 = FaultSystemRupSet.load(new File(invsDir,
				"2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		File u3Haz31File = new File(invsDir,
				"2022_12_08-u3-FM3_1-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File u3Converged31HazFile = new File(invsDir,
				"2022_03_24-u3_branches-FM3_1-2000ip-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File u3HazFile = new File(invsDir,
				"2023_02_09-u3-both_fms-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File u3_23GridHazFile = new File(invsDir,
				"2023_09_21-u3-both_fms-ba_only-nshm23_gridded/"
				+ "results_hazard_include_0.1deg.zip");
		
		File methodsU3GridHazFile = new File(invsDir,
				"2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-u3_gridded/"
				+ "results_hazard_include_0.1deg.zip");
		
		File methods23GridHazFile = new File(invsDir,
				"2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded/"
				+ "results_hazard_include_0.1deg.zip");
		
		File modelHazFile = new File(invsDir,
				"2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		String entryName = "mean_map_pga_TWO_IN_50.txt";
		String hazLabel = "PGA, 2% in 50 yrs";
		
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		Region[] caRegions = new Region[caOutlines.length];
		for (int i=0; i<caOutlines.length; i++) {
			LocationList outline = new LocationList();
			for (Point2D pt : caOutlines[i])
				outline.add(new Location(pt.getY(), pt.getX()));
			caRegions[i] = new Region(outline, BorderType.MERCATOR_LINEAR);
		}

		GriddedGeoDataSet u3Map = mask(caRegions, loadXYZ(u3HazFile, entryName));
		GriddedGeoDataSet u3fm31Map = mask(caRegions, loadXYZ(u3Haz31File, entryName));
		GriddedGeoDataSet u3ConvergedMap = mask(caRegions, loadXYZ(u3Converged31HazFile, entryName));
		GriddedGeoDataSet u3_23GridMap = mask(caRegions, loadXYZ(u3_23GridHazFile, entryName));
		GriddedGeoDataSet methodsU3GridMap = mask(caRegions, loadXYZ(methodsU3GridHazFile, entryName));
		GriddedGeoDataSet methods23GridMap = mask(caRegions, loadXYZ(methods23GridHazFile, entryName));
		GriddedGeoDataSet modelMap = mask(caRegions, loadXYZ(modelHazFile, entryName));
		
		GriddedRegion refReg = u3Map.getRegion();
		GeographicMapMaker mapMakerU3 = new RupSetMapMaker(rupSetU3, refReg);
		GeographicMapMaker mapMaker23 = new RupSetMapMaker(rupSet23, refReg);
		
		for (GeographicMapMaker mapMaker : new GeographicMapMaker[] {mapMakerU3, mapMaker23}) {
			mapMaker.setSectOutlineChar(null);
			mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(100, 100, 100, 127)));
		}
		
		List<XYAnnotation> methodAnns = new ArrayList<>();
		
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		Stroke baseStroke = PlotLineType.SOLID.buildStroke(2f);
//		Color lineColor = new Color(0, 0, 0, 127);
		Color lineColor = Color.BLACK;
		
		double annX = -115.7;
		double annY = 37.75;
		double lineX = -117;
		double lineY = 36.8;
		XYTextAnnotation ann = new XYTextAnnotation("Death Valley (No)", annX, annY);
		ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
		ann.setFont(annFont);
		methodAnns.add(ann);
		methodAnns.add(new XYLineAnnotation(lineX, lineY, annX, annY-0.1, baseStroke, lineColor));
		
		annX = -121.5;
		annY = 35;
		lineX = -119.2;
		lineY = 34.5;
		ann = new XYTextAnnotation("San Cayetano", annX, annY);
		ann.setTextAnchor(TextAnchor.CENTER_RIGHT);
		ann.setFont(annFont);
		methodAnns.add(ann);
		methodAnns.add(new XYLineAnnotation(lineX, lineY, annX+0.1, annY-0.1, baseStroke, lineColor));
		
		List<XYAnnotation> ingredAnns = new ArrayList<>();
		
		annX = -114.8;
		annY = 36.6;
		lineX = -116.3;
		lineY = 35.4;
		ann = new XYTextAnnotation("E.C.S.Z.", annX, annY);
		ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
		ann.setFont(annFont);
		ingredAnns.add(ann);
		ingredAnns.add(new XYLineAnnotation(lineX, lineY, annX, annY-0.1, baseStroke, lineColor));
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(new Color(255, 255, 255, 0));
		plotHazardChange(outputDir, "u3_converged_vs_u3", mapMakerU3, pDiffCPT, refReg, u3ConvergedMap, u3fm31Map,
				"Converged UCERF3 vs UCERF3, % Change, "+hazLabel, logFW);
		plotHazardChange(outputDir, "methodology_vs_u3", mapMakerU3, pDiffCPT, refReg, methodsU3GridMap, u3Map,
				"NSHM23 Methodology vs UCERF3, % Change, "+hazLabel, logFW, methodAnns);
		plotHazardChange(outputDir, "ingredients", mapMaker23, pDiffCPT, refReg, modelMap, methods23GridMap,
				"NSHM23 vs UCERF3 Ingredients, % Change, "+hazLabel, logFW, ingredAnns);
		plotHazardChange(outputDir, "full_change", mapMaker23, pDiffCPT, refReg, modelMap, u3_23GridMap,
				"NSHM23 vs UCERF3, % Change, "+hazLabel, logFW);
		plotHazardChange(outputDir, "test_method_grid_change", mapMakerU3, pDiffCPT, refReg, methods23GridMap, methodsU3GridMap,
				"Methods, Grid23 vs GridU3, % Change, "+hazLabel, logFW);
		
		// attribution
		CPT attributionCPT = getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse(), 10d, 50d);
		
		GriddedGeoDataSet attXYZ = new GriddedGeoDataSet(u3Map.getRegion(), false);
		
		int numMethod = 0;
		int numIngredient = 0;
		int numWithin = 0;
		for (int i=0; i<attXYZ.size(); i++) {
			double u3Val = u3Map.get(i);
			Location loc = attXYZ.getLocation(i);
			
			double u3_23GridVal = u3_23GridMap.get(u3_23GridMap.indexOf(loc));
			double methodsU3GridVal = methodsU3GridMap.get(methodsU3GridMap.indexOf(loc));
			double methods23GridVal = methods23GridMap.get(methods23GridMap.indexOf(loc));
			double modelVal = modelMap.get(modelMap.indexOf(loc));
			
			double totalDiff = Math.abs(modelVal - u3_23GridVal);
			double modelDiff = Math.abs(modelVal - methods23GridVal);
			double methodsDiff = Math.abs(methodsU3GridVal - u3Val);
			
			double pDiff = 100d*totalDiff/u3Val;
			if (modelDiff < methodsDiff)
				// when methods dominate, go negative
				pDiff = -pDiff;
			attXYZ.set(i, pDiff);
			if (Double.isFinite(pDiff)) {
				if (Math.abs(pDiff) > 10d) {
					if (modelDiff > methodsDiff)
						numIngredient++;
					else
						numMethod++;
				} else {
					numWithin++;
				}
			}
		}
		
		int totNum = numMethod + numIngredient + numWithin;
		logPrint(logFW, "Attribution stats:");
		logPrint(logFW, "\tWithin 10%: "+pDF.format((double)numWithin/(double)totNum));
		logPrint(logFW, "\tMethodology: "+pDF.format((double)numMethod/(double)totNum));
		logPrint(logFW, "\tIngredients: "+pDF.format((double)numIngredient/(double)totNum));
		logPrint(logFW, "Attribution stats, of those exceeding 10%:");
		int numExceeding = totNum - numWithin;
		logPrint(logFW, "\tMethodology: "+pDF.format((double)numMethod/(double)numExceeding));
		logPrint(logFW, "\tIngredients: "+pDF.format((double)numIngredient/(double)numExceeding));
		
		mapMaker23.plotXYZData(attXYZ, attributionCPT, "Methodology     ←     |Hazard % Change|     →     Ingredients");
		PlotSpec plot = mapMaker23.buildPlot(" ");
		int width = mapMaker23.getDefaultPlotWidth();
		// plot both sides as positive
		DecimalFormat intDF = new DecimalFormat("0");
		NumberFormat format = new NumberFormat() {
			
			@Override
			public Number parse(String source, ParsePosition parsePosition) {
				return intDF.parse(source, parsePosition);
			}
			
			@Override
			public StringBuffer format(long number, StringBuffer toAppendTo, FieldPosition pos) {
				if (number < 0l)
					number = -number;
				if (number == -0l)
					number = 0l;
				return intDF.format(number, toAppendTo, pos);
			}
			
			@Override
			public StringBuffer format(double number, StringBuffer toAppendTo, FieldPosition pos) {
				return intDF.format(Math.abs(number), toAppendTo, pos);
			}
		};
		mapMaker23.plot(outputDir, "comb_hazard_attribution", plot, width, new Consumer<HeadlessGraphPanel>() {

			@Override
			public void accept(HeadlessGraphPanel gp) {
				for (Title subtitle : (List<Title>)gp.getPlot().getChart().getSubtitles()) {
					if (subtitle instanceof PaintScaleLegend) {
						NumberAxis axis = ((NumberAxis)((PaintScaleLegend)subtitle).getAxis());
						PlotUtils.setTick(axis, 10d);
						axis.setNumberFormatOverride(format);
					}
				}
			}
		});
		
		logFW.close();
	}
	
	public static void doWUS() throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/wus_haz_change_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FileWriter logFW = new FileWriter(new File(outputDir, "log.txt"));
		
		FaultSystemRupSet rupSet18 = FaultSystemRupSet.load(new File(invsDir,
				"2023_04_13-nshm18_branches-new_scale-u3_paleo-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM18_WUS_PlusU3_FM_3p1_CoulombRupSet_branch_averaged.zip"));
		FaultSystemRupSet rupSet23 = FaultSystemRupSet.load(new File(invsDir,
				"2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		File nshm18HazFile = new File(invsDir,
				"2023_01_27-nshm18-hazard-ask2014-0.1deg-noSub/"
				+ "results_hazard.zip");
		
		File nshm18GridOnlyHazFile = new File(invsDir,
				"2023_01_27-nshm18-hazard-ask2014-0.1deg-noSub-griddedOnly/"
				+ "results_hazard.zip");
		
		File nshm18_23GridHazFile = new File(invsDir,
				"2023_07_10-nshm18-grid_src_from_23-wus-hazard-ask2014-0.1deg-noSub/"
				+ "results_hazard.zip");
		
		File methods23GridHazFile = new File(invsDir,
				"2023_04_13-nshm18_branches-new_scale-u3_paleo-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded/"
				+ "results_hazard_include_0.1deg.zip");
		
		File modelHazFile = new File(invsDir,
				"2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File modelGridOnlyHazFile = new File(invsDir,
				"2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_only_0.1deg.zip");
		
		String nshmHazEntryName = "map_pga_TWO_IN_50.txt";
		String entryName = "mean_map_pga_TWO_IN_50.txt";
		String hazLabel = "PGA, 2% in 50 yrs";
		
		List<Region> wusRegionsList = new ArrayList<>();
		for (NSHM23_SingleStates state : NSHM23_SingleStates.values()) {
			XY_DataSet[] outlines = PoliticalBoundariesData.loadUSState(state.getStateName());
			for (int i=0; i<outlines.length; i++) {
				LocationList outline = new LocationList();
				for (Point2D pt : outlines[i])
					outline.add(new Location(pt.getY(), pt.getX()));
				wusRegionsList.add(new Region(outline, BorderType.MERCATOR_LINEAR));
			}
		}
		Region[] wusRegions = wusRegionsList.toArray(new Region[0]);

		GriddedGeoDataSet nshm18Map = mask(wusRegions, loadXYZ(nshm18HazFile, nshmHazEntryName));
		GriddedGeoDataSet nshm18GridOnlyMap = mask(wusRegions, loadXYZ(nshm18GridOnlyHazFile, nshmHazEntryName));
		GriddedGeoDataSet u3_23GridMap = mask(wusRegions, loadXYZ(nshm18_23GridHazFile, nshmHazEntryName));
		GriddedGeoDataSet methods23GridMap = mask(wusRegions, loadXYZ(methods23GridHazFile, entryName));
		GriddedGeoDataSet modelMap = mask(wusRegions, loadXYZ(modelHazFile, entryName));
		GriddedGeoDataSet modelGridOnlyMap = mask(wusRegions, loadXYZ(modelGridOnlyHazFile, entryName));
		
		GriddedRegion refReg = nshm18Map.getRegion();
		GeographicMapMaker mapMaker18 = new RupSetMapMaker(rupSet18, refReg);
		GeographicMapMaker mapMaker23 = new RupSetMapMaker(rupSet23, refReg);
		
		for (GeographicMapMaker mapMaker : new GeographicMapMaker[] {mapMaker18, mapMaker23}) {
			mapMaker.setSectOutlineChar(null);
			mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(100, 100, 100, 127)));
		}
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(new Color(255, 255, 255, 0));
		plotHazardChange(outputDir, "methodology", mapMaker18, pDiffCPT, refReg, methods23GridMap, u3_23GridMap,
				"NSHM23 Methodology vs NSHM18, % Change, "+hazLabel, logFW);
		plotHazardChange(outputDir, "ingredients", mapMaker23, pDiffCPT, refReg, modelMap, methods23GridMap,
				"NSHM23 vs NSHM18 Ingredients, % Change, "+hazLabel, logFW);
		plotHazardChange(outputDir, "full_change", mapMaker23, pDiffCPT, refReg, modelMap, u3_23GridMap,
				"NSHM23 vs NSHM18, % Change, "+hazLabel, logFW);
//		plotHazardChange(outputDir, "test_method_grid_change", mapMaker18, pDiffCPT, refReg, methods23GridMap, methodsU3GridMap,
//				"Methods, Grid23 vs GridU3, % Change, "+hazLabel);
		
		// attribution
		CPT attributionCPT = getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse(), 10d, 50d);
		CPT cork = GMT_CPT_Files.DIVERGING_CORK_UNIFORM.instance();
		// CORK gets too dark at the edges (hard to tell color), trim it
		cork = cork.rescale(-1d, 1d);
		CPT modCork = new CPT();
		for (double v=-0.81; (float)v<=0.8f; v+=0.01) {
			float start = (float)(v-0.01);
			Color startColorCorc = cork.getColor(start);
			float end = (float)v;
			Color endColorCorc = cork.getColor(end);
			modCork.add(new CPTVal(start, startColorCorc, end, endColorCorc));
		}
		modCork.setBelowMinColor(modCork.getMinColor());
		modCork.setAboveMaxColor(modCork.getMaxColor());
		CPT faultGridAttributionCPT = getCenterMaskedCPT(modCork, 10d, 50d);
		
		GriddedGeoDataSet attXYZ = new GriddedGeoDataSet(nshm18Map.getRegion(), false);
		GriddedGeoDataSet faultGridAttXYZ = new GriddedGeoDataSet(nshm18Map.getRegion(), false);
		GriddedGeoDataSet grid18PlusFaults23 = new GriddedGeoDataSet(nshm18Map.getRegion(), false);
		
		int numMethod = 0;
		int numIngredient = 0;
		int numWithin = 0;
		for (int i=0; i<attXYZ.size(); i++) {
//			double u3Val = nshm18Map.get(i);
			Location loc = attXYZ.getLocation(i);
			
			double u3_23GridVal = u3_23GridMap.get(u3_23GridMap.indexOf(loc));
//			double methodsU3GridVal = methodsU3GridMap.get(methodsU3GridMap.indexOf(loc));
			double methods23GridVal = methods23GridMap.get(methods23GridMap.indexOf(loc));
			double modelVal = modelMap.get(modelMap.indexOf(loc));
			
			double totalDiff = Math.abs(modelVal - u3_23GridVal);
			double modelDiff = Math.abs(modelVal - methods23GridVal);
			double methodsDiff = Math.abs(methods23GridVal - u3_23GridVal);
			
			double pDiff = 100d*totalDiff/u3_23GridVal;
			if (modelDiff < methodsDiff)
				// when methods dominate, go negative
				pDiff = -pDiff;
			attXYZ.set(i, pDiff);
			if (Double.isFinite(pDiff)) {
				if (Math.abs(pDiff) > 10d) {
					if (modelDiff > methodsDiff)
						numIngredient++;
					else
						numMethod++;
				} else {
					numWithin++;
				}
			}
			
			double grid18val = nshm18GridOnlyMap.get(nshm18GridOnlyMap.indexOf(loc));
			double grid23val = modelGridOnlyMap.get(modelGridOnlyMap.indexOf(loc));
			double faultAdd = Math.max(0d, modelVal - grid23val);
			double grid18PlusFaults23val = grid18val+faultAdd;
			grid18PlusFaults23.set(i, grid18PlusFaults23val);
			
			double full18val = nshm18Map.get(nshm18Map.indexOf(loc));
			
			double fullDiff = modelVal - full18val;
			
			double faultDiff = Math.abs(totalDiff);
			double gridDiff = Math.abs(u3_23GridVal - full18val);
			
			double sign = faultDiff > gridDiff ? 1d : -1d;
			
			double fullPDiff = Math.abs(100d*fullDiff/full18val);
			
			faultGridAttXYZ.set(i, sign*fullPDiff);
		}

		plotHazardChange(outputDir, "grid_change", mapMaker23, pDiffCPT, refReg, modelMap, grid18PlusFaults23,
				"NSHM23 vs NSHM18, Gridded, % Change, "+hazLabel, logFW);
		plotHazardChange(outputDir, "total_change", mapMaker23, pDiffCPT, refReg, modelMap, nshm18Map,
				"NSHM23 vs NSHM18, Total, % Change, "+hazLabel, logFW);
		
		int totNum = numMethod + numIngredient + numWithin;
		logPrint(logFW, "Attribution stats:");
		logPrint(logFW, "\tWithin 10%: "+pDF.format((double)numWithin/(double)totNum));
		logPrint(logFW, "\tMethodology: "+pDF.format((double)numMethod/(double)totNum));
		logPrint(logFW, "\tIngredients: "+pDF.format((double)numIngredient/(double)totNum));
		logPrint(logFW, "Attribution stats, of those exceeding 10%:");
		int numExceeding = totNum - numWithin;
		logPrint(logFW, "\tMethodology: "+pDF.format((double)numMethod/(double)numExceeding));
		logPrint(logFW, "\tIngredients: "+pDF.format((double)numIngredient/(double)numExceeding));
		
		mapMaker23.plotXYZData(attXYZ, attributionCPT, "Methodology     ←     |Hazard % Change|     →     Ingredients");
		PlotSpec plot = mapMaker23.buildPlot(" ");
		int width = mapMaker23.getDefaultPlotWidth();
		// plot both sides as positive
		DecimalFormat intDF = new DecimalFormat("0");
		NumberFormat format = new NumberFormat() {
			
			@Override
			public Number parse(String source, ParsePosition parsePosition) {
				return intDF.parse(source, parsePosition);
			}
			
			@Override
			public StringBuffer format(long number, StringBuffer toAppendTo, FieldPosition pos) {
				if (number < 0l)
					number = -number;
				if (number == -0l)
					number = 0l;
				return intDF.format(number, toAppendTo, pos);
			}
			
			@Override
			public StringBuffer format(double number, StringBuffer toAppendTo, FieldPosition pos) {
				return intDF.format(Math.abs(number), toAppendTo, pos);
			}
		};
		mapMaker23.plot(outputDir, "comb_hazard_attribution", plot, width, new Consumer<HeadlessGraphPanel>() {

			@Override
			public void accept(HeadlessGraphPanel gp) {
				for (Title subtitle : (List<Title>)gp.getPlot().getChart().getSubtitles()) {
					if (subtitle instanceof PaintScaleLegend) {
						NumberAxis axis = ((NumberAxis)((PaintScaleLegend)subtitle).getAxis());
						PlotUtils.setTick(axis, 10d);
						axis.setNumberFormatOverride(format);
					}
				}
			}
		});
		
		mapMaker23.plotXYZData(faultGridAttXYZ, faultGridAttributionCPT, "Gridded      ←        |Hazard % Change|        →       Faults");
		plot = mapMaker23.buildPlot(" ");
		mapMaker23.plot(outputDir, "fault_gridded_hazard_attribution", plot, width, new Consumer<HeadlessGraphPanel>() {

			@Override
			public void accept(HeadlessGraphPanel gp) {
				for (Title subtitle : (List<Title>)gp.getPlot().getChart().getSubtitles()) {
					if (subtitle instanceof PaintScaleLegend) {
						NumberAxis axis = ((NumberAxis)((PaintScaleLegend)subtitle).getAxis());
						PlotUtils.setTick(axis, 10d);
						axis.setNumberFormatOverride(format);
					}
				}
			}
		});
		
		logFW.close();
	}
	
	private static void logPrint(FileWriter fw, String line) throws IOException {
		System.out.println(line);
		fw.write(line+"\n");
	}
	
	static GriddedGeoDataSet loadXYZ(File zipFile, String entryName) throws IOException {
		System.out.println("Loading "+entryName+" from "+zipFile.getAbsolutePath());
		ZipFile zip = new ZipFile(zipFile);
		
		ZipEntry regEntry = zip.getEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(regEntry)));
		GriddedRegion gridReg = GriddedRegion.fromFeature(Feature.read(bRead));
		
		ZipEntry mapEntry = zip.getEntry(entryName);
		InputStream is = zip.getInputStream(mapEntry);
		Preconditions.checkNotNull(is, "IS is null for %s", entryName);
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(mapEntry)));
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
		
		zip.close();
		
		return xyz;
	}
	
	public static CPT getCenterMaskedCPT(CPT cptBasis, double maskRange, double fullRange) {
		CPT leftCPT = getHalfCPT(cptBasis.reverse(), maskRange, fullRange);
		CPT rightCPT = getHalfCPT(cptBasis, maskRange, fullRange);
		
		CPT attributionCPT = new CPT();
		for (int i=leftCPT.size(); --i>=0;) {
			CPTVal val = leftCPT.get(i);
			attributionCPT.add(new CPTVal(-val.end, val.maxColor, -val.start, val.minColor));
		}
		for (CPTVal val : rightCPT)
			attributionCPT.add(val);
		attributionCPT.setBelowMinColor(leftCPT.getAboveMaxColor());
		attributionCPT.setAboveMaxColor(rightCPT.getAboveMaxColor());
		attributionCPT.setNanColor(new Color(255, 255, 255, 0));
		
		return attributionCPT;
	}
	
	private static CPT getHalfCPT(CPT cpt, double maskRange, double fullRange) {
		cpt = cpt.rescale(-1d, 1d);
		CPT ret = new CPT();
		Color zeroColor = cpt.getColor(0f);
		ret.add(new CPTVal(0f, zeroColor, (float)maskRange, zeroColor));
		// now scale up to full range
		for (double x=0.01; x<1d; x+=0.01) {
			double mappedX = maskRange + (fullRange-maskRange)*x;
			
			double prevX = ret.getMaxValue();
			Color prevColor = ret.getMaxColor();
			
			ret.add(new CPTVal(prevX, prevColor, (float)mappedX, cpt.getColor((float)x)));
		}
		ret.add(new CPTVal(ret.getMaxValue(), ret.getMaxColor(), (float)fullRange, cpt.getMaxColor()));
		ret.setBelowMinColor(ret.getMinColor());
		ret.setAboveMaxColor(ret.getMaxColor());
		ret.setNanColor(Color.WHITE);
		return ret;
	}
	
	static GriddedGeoDataSet mask(Region[] includeRegions, GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		GriddedRegion gridReg = xyz.getRegion();
		double halfLatSpacing = gridReg.getLatSpacing()*0.5;
		double halfLonSpacing = gridReg.getLonSpacing()*0.5;
		for (int i=0; i<xyz.size(); i++) {
			boolean inside = false;
			Location center = xyz.getLocation(i);
			
			// first just check center
			for (Region reg : includeRegions) {
				if (reg.contains(center)) {
					inside = true;
					break;
				}
			}
			if (!inside) {
				// test full cell intersection
				Location upLeft = new Location(center.getLatitude()+halfLatSpacing, center.getLongitude()-halfLonSpacing);
				Location botRight = new Location(center.getLatitude()-halfLatSpacing, center.getLongitude()+halfLonSpacing);
				Region cell = new Region(upLeft, botRight);
				
				for (Region reg : includeRegions) {
					double dist = reg.distanceToLocation(center);
					if (dist < 30d && Region.intersect(cell, reg) != null) {
						inside = true;
						break;
					}
				}
			}
			if (!inside)
				xyz.set(i, Double.NaN);
		}
		return xyz;
	}
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static void plotHazardChange(File outputDir, String prefix, GeographicMapMaker mapMaker, CPT pDiffCPT,
			GriddedRegion refReg, GriddedGeoDataSet numerator, GriddedGeoDataSet denominator, String label,
			FileWriter logFW) throws IOException {
		plotHazardChange(outputDir, prefix, mapMaker, pDiffCPT, refReg, numerator, denominator, label, logFW, null);
	}
	
	private static void plotHazardChange(File outputDir, String prefix, GeographicMapMaker mapMaker, CPT pDiffCPT,
			GriddedRegion refReg, GriddedGeoDataSet numerator, GriddedGeoDataSet denominator, String label,
			FileWriter logFW, List<? extends XYAnnotation> anns) throws IOException {
		GriddedGeoDataSet pDiff = new GriddedGeoDataSet(refReg, false);
		
		MinMaxAveTracker meanAbsTrack = new MinMaxAveTracker();
		MinMaxAveTracker meanTrack = new MinMaxAveTracker();
		int numWithin10 = 0;
		int numWithin5 = 0;
		int numWithin1 = 0;
		int numValid = 0;
		for (int i=0; i<pDiff.size(); i++) {
			Location loc = pDiff.getLocation(i);
			
			double v1 = numerator.get(numerator.indexOf(loc));
			double v2 = denominator.get(denominator.indexOf(loc));
			
			if (Double.isFinite(v1) && Double.isFinite(v2)) {
				double val = 100d*(v1-v2)/v2;
				pDiff.set(i, val);
				
				numValid++;
				if (Math.abs(val) < 10d)
					numWithin10++;
				if (Math.abs(val) < 5d)
					numWithin5++;
				if (Math.abs(val) < 1d)
					numWithin1++;
				
				meanAbsTrack.addValue(Math.abs(val));
				meanTrack.addValue(val);
			} else {
				pDiff.set(i, Double.NaN);
			}
		}
		
		logPrint(logFW, "Plotting "+prefix+", "+label);
		logPrint(logFW, "\tRange: ["+twoDigits.format(meanTrack.getMin())+"%, "+twoDigits.format(meanTrack.getMax())+"%]");
		logPrint(logFW, "\tAverage: "+twoDigits.format(meanTrack.getAverage())+"%");
		logPrint(logFW, "\tAverage Absolute: "+twoDigits.format(meanAbsTrack.getAverage())+"%");
		logPrint(logFW, "\tWithin 1%: "+pDF.format((double)numWithin1/(double)numValid));
		logPrint(logFW, "\tWithin 5%: "+pDF.format((double)numWithin5/(double)numValid));
		logPrint(logFW, "\tWithin 10%: "+pDF.format((double)numWithin10/(double)numValid));
		
		
		mapMaker.plotXYZData(pDiff, pDiffCPT, label);
		
		if (anns != null)
			mapMaker.setAnnotations(anns);
		mapMaker.plot(outputDir, prefix, " ");
		mapMaker.clearAnnotations();
	}

}
