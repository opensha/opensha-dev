package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import scratch.kevin.nshm23.SimpleSmoothHazardMapCalc;
import scratch.kevin.nshm23.SingleSiteHazardAndDataComparisonPageGen;
import scratch.kevin.nshm23.SingleSiteHazardAndDataComparisonPageGen.RegionalParticipationResult;

public class WrapperComparisonPageGen {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/markdown/nshm23-misc/wrapper_rate_comparisons");
		File modelsDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models");
		
		String name1 = "NSHM23";
		String dirName1 = "nshm-conus-6.a.6";
		String name2 = "NSHM18";
		String dirName2 = "nshm-conus-5.3.0";
		File outputDir = new File(baseDir, "6.a.6_vs_5.3.0");
		Region reg = NSHM23_RegionLoader.loadFullConterminousUS();
		
		File extGriddedHazFile1 = new File(modelsDir, "ext_hazard_calcs/conus-2023-erf-6a6-grid/map_conus-2023-erf-6a6-GRID_vs760_PGA_02475yrs.gmt");
		File extGriddedHazFile2 = new File(modelsDir, "ext_hazard_calcs/conus-2018-530-GRID/map_conus-2018-530-GRID_vs760_PGA_02475yrs.gmt");
		
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		double period = 0d;
		String hazLabel = "PGA, "+rp.label;
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(List.of(), reg);
		mapMaker.setDefaultPlotWidth(1200);
		
		HazardModel model1 = HazardModel.load(new File(modelsDir, dirName1).toPath());
		HazardModel model2 = HazardModel.load(new File(modelsDir, dirName2).toPath());
		
		GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
		GriddedRegion smoothHazGridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Wrapped "+name1+" vs "+name2);
		lines.add("");
		
		lines.add("This comparison uses a wrapper to nshmp-haz to compare "+name1+" (tag: `"+dirName1
				+"`) to "+name2+" (tag: `"+dirName2+"`).");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW, TectonicRegionType.STABLE_SHALLOW);
		boolean subduction = false;
		if (subduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}
		
		double[] magThresholds = { 5d, 6d, 7d, 7.5d };
		
		NshmErf faultERF1 = new NshmErf(model1, trts, IncludeBackgroundOption.EXCLUDE);
		faultERF1.updateForecast();
		NshmErf gridERF1 = new NshmErf(model1, trts, IncludeBackgroundOption.ONLY);
		gridERF1.updateForecast();
		
		NshmErf faultERF2 = new NshmErf(model2, trts, IncludeBackgroundOption.EXCLUDE);
		faultERF2.updateForecast();
		NshmErf gridERF2 = new NshmErf(model2, trts, IncludeBackgroundOption.ONLY);
		gridERF2.updateForecast();
		
		RegionalParticipationResult faultPartic1 = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				faultERF1, gridReg, magThresholds, null);
		RegionalParticipationResult gridPartic1 = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				gridERF1, gridReg, magThresholds, null);
		RegionalParticipationResult faultPartic2 = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				faultERF2, gridReg, magThresholds, null);
		RegionalParticipationResult gridPartic2 = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				gridERF2, gridReg, magThresholds, null);
		
		IncludeBackgroundOption[] griddedOps = {IncludeBackgroundOption.EXCLUDE,
				IncludeBackgroundOption.ONLY, IncludeBackgroundOption.INCLUDE};
		
		Color transparent = new Color(255, 255, 255, 0);
		
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		hazCPT.setNanColor(transparent);
		
		CPT pDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance(), 10d, 50d);
		pDiffCPT.setNanColor(transparent);
		
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2d, 0.2d);
		diffCPT.setNanColor(transparent);
		
		CPT moCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(10d, 18d);
		moCPT.setNanColor(transparent);
		
		CPT moRateDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1e16, 1e16);
		moRateDiffCPT.setNanColor(transparent);
		
		CPT pDiff100CPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-100d, 100d);
		pDiff100CPT.setNanColor(transparent);
		
		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;
		
		boolean writeAsWeGo = !new File(outputDir, "index.html").exists();
		
		for (IncludeBackgroundOption griddedOp : griddedOps) {
			String label, prefix;
			GriddedGeoDataSet[] particRates1, particRates2;
			GriddedGeoDataSet[] nuclRates1, nuclRates2;
			GriddedGeoDataSet moRates1, moRates2;
			File extMap1 = null, extMap2 = null;
			if (griddedOp == IncludeBackgroundOption.EXCLUDE) {
				label = "Faults Only";
				
				particRates1 = faultPartic1.particRateMaps;
				particRates2 = faultPartic2.particRateMaps;
				nuclRates1 = faultPartic1.nuclRateMaps;
				nuclRates2 = faultPartic2.nuclRateMaps;
				moRates1 = faultPartic1.momentRateMap;
				moRates2 = faultPartic2.momentRateMap;
				
				prefix = "fault";
			} else if (griddedOp == IncludeBackgroundOption.ONLY) {
				label = "Gridded Seismicity Only";
				
				particRates1 = gridPartic1.particRateMaps;
				particRates2 = gridPartic2.particRateMaps;
				nuclRates1 = gridPartic1.nuclRateMaps;
				nuclRates2 = gridPartic2.nuclRateMaps;
				moRates1 = gridPartic1.momentRateMap;
				moRates2 = gridPartic2.momentRateMap;
				
				extMap1 = extGriddedHazFile1;
				extMap2 = extGriddedHazFile2;
				prefix = "grid";
			} else {
				Preconditions.checkState(griddedOp == IncludeBackgroundOption.INCLUDE);
				label = "Faults and Gridded Seismicity";

				particRates1 = new GriddedGeoDataSet[magThresholds.length];
				particRates2 = new GriddedGeoDataSet[magThresholds.length];
				nuclRates1 = new GriddedGeoDataSet[magThresholds.length];
				nuclRates2 = new GriddedGeoDataSet[magThresholds.length];
				for (int m=0; m<magThresholds.length; m++) {
					particRates1[m] = WUS_HazardChangePageGen.sumMap(faultPartic1.particRateMaps[m], gridPartic1.particRateMaps[m]);
					particRates2[m] = WUS_HazardChangePageGen.sumMap(faultPartic2.particRateMaps[m], gridPartic2.particRateMaps[m]);
					nuclRates1[m] = WUS_HazardChangePageGen.sumMap(faultPartic1.nuclRateMaps[m], gridPartic1.nuclRateMaps[m]);
					nuclRates2[m] = WUS_HazardChangePageGen.sumMap(faultPartic2.nuclRateMaps[m], gridPartic2.nuclRateMaps[m]);
				}
				
				moRates1 = WUS_HazardChangePageGen.sumMap(faultPartic1.momentRateMap, gridPartic1.momentRateMap);
				moRates2 = WUS_HazardChangePageGen.sumMap(faultPartic2.momentRateMap, gridPartic2.momentRateMap);
				prefix = "combined";
			}
			
			lines.add("## "+label);
			lines.add(topLink); lines.add("");
			
			boolean hazard = extMap1 != null && extMap2 != null;
			
			if (hazard) {
				// plot hazard
				
				lines.add("### Hazard Maps, "+label);
				lines.add(topLink); lines.add("");
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				
				table.addLine(name1, name2);
				
				GriddedGeoDataSet haz1 = loadExtMap(extGriddedHazFile1, gridReg);
				GriddedGeoDataSet haz2 = loadExtMap(extGriddedHazFile2, gridReg);
				
				table.initNewLine();
				mapMaker.plotXYZData(asLog10(haz1), hazCPT, name1+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_haz1", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_haz1.png)");
				mapMaker.plotXYZData(asLog10(haz2), hazCPT, name2+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_haz2", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_haz2.png)");
				table.finalizeLine();
				
				GriddedGeoDataSet pDiff = WUS_HazardChangePageGen.mapPDiff(haz1, haz2);
				GriddedGeoDataSet diff = WUS_HazardChangePageGen.mapDiff(haz1, haz2);
				
				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				table.initNewLine();
				mapMaker.plotXYZData(pDiff, pDiffCPT, name1+" vs "+name2+", % Change, "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_haz_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_haz_pDiff.png)");
				mapMaker.plotXYZData(diff, diffCPT, name1+" - "+name2+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_haz_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_haz_diff.png)");
				table.finalizeLine();
				
				lines.addAll(table.build());
				lines.add("");
			}
			
			// rate maps
			lines.add("### Rate Maps, "+label);
			lines.add(topLink); lines.add("");
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			for (int m=0; m<magThresholds.length; m++) {
				String magPrefix = "m"+oDF.format(magThresholds[m]);
				String magLabel = "M≥"+oDF.format(magThresholds[m]);
				table.addLine(MarkdownUtils.boldCentered("NSHM23, "+magLabel), MarkdownUtils.boldCentered("NSHM18, "+magLabel));
				
				GriddedGeoDataSet rates1 = zerosToNaNs(particRates1[m]);
				GriddedGeoDataSet rates2 = zerosToNaNs(particRates2[m]);
				
				double maxRate = 0d;
				for (int i=0; i<gridReg.getNodeCount(); i++) {
					double rate1 = rates1.get(i);
					if (Double.isFinite(rate1))
						maxRate = Math.max(maxRate, rate1);
					double rate2 = rates2.get(i);
					if (Double.isFinite(rate2))
						maxRate = Math.max(maxRate, rate2);
				}
				double logMax = maxRate > 0d ? Math.ceil(Math.log10(maxRate)) : -2d;
				
				CPT nuclCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMax-6d, logMax);
				nuclCPT.setNanColor(transparent);
				
//				double diffMax = Math.pow(10, (int)(Math.log10(diffStat)+0.5));
				double diffMax = Math.pow(10, logMax-2);
				
				CPT rateDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-diffMax, diffMax);
				rateDiffCPT.setNanColor(transparent);
				
				table.initNewLine();
				
				mapMaker.plotXYZData(asLog10(rates1), nuclCPT, "Participation Rate, "+magLabel);
				mapMaker.plot(resourcesDir, prefix+"_"+magPrefix+"_1", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_"+magPrefix+"_1.png)");
				
				mapMaker.plotXYZData(asLog10(rates2), nuclCPT, "Participation Rate, "+magLabel);
				mapMaker.plot(resourcesDir, prefix+"_"+magPrefix+"_2", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_"+magPrefix+"_2.png)");
				
				table.finalizeLine();
				
				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				table.initNewLine();
				
				GriddedGeoDataSet ratePDiffs = WUS_HazardChangePageGen.mapPDiff(rates1, rates2);
				GriddedGeoDataSet rateDiffs = WUS_HazardChangePageGen.mapDiff(rates1, rates2);
				
				mapMaker.plotXYZData(ratePDiffs, pDiff100CPT, "% Change in Participation Rate, "+magLabel);
				mapMaker.plot(resourcesDir, prefix+"_"+magPrefix+"_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_"+magPrefix+"_pDiff.png)");
				
				mapMaker.plotXYZData(rateDiffs, rateDiffCPT, "Participation Rate Difference, "+magLabel);
				mapMaker.plot(resourcesDir, prefix+"_"+magPrefix+"_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_"+magPrefix+"_diff.png)");
				
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
			
			if (hazard) {
				// rate-based simplified hazard maps
				lines.add("### Rate-Based Simplified Hazard Maps, "+label);
				lines.add(topLink); lines.add("");
				
				GutenbergRichterMagFreqDist mfdShape = new GutenbergRichterMagFreqDist(magThresholds[0], 6, 0.5, 1e16, 1d);
				
				GriddedGeoDataSet smoothHaz1 = SimpleSmoothHazardMapCalc.hazMapRateSmooth(
						nuclRates1[0], gmpeRef, mfdShape, period, rp, smoothHazGridReg);
				GriddedGeoDataSet smoothHaz2 = SimpleSmoothHazardMapCalc.hazMapRateSmooth(
						nuclRates2[0], gmpeRef, mfdShape, period, rp, smoothHazGridReg);
				
				table = MarkdownUtils.tableBuilder();
				
				table.addLine(name1, name2);
				
				table.initNewLine();
				mapMaker.plotXYZData(asLog10(smoothHaz1), hazCPT, name1+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_rate_haz1", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_rate_haz1.png)");
				mapMaker.plotXYZData(asLog10(smoothHaz2), hazCPT, name2+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_rate_haz2", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_rate_haz2.png)");
				table.finalizeLine();
				
				GriddedGeoDataSet pDiff = WUS_HazardChangePageGen.mapPDiff(smoothHaz1, smoothHaz2);
				GriddedGeoDataSet diff = WUS_HazardChangePageGen.mapDiff(smoothHaz1, smoothHaz2);
				
				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				table.initNewLine();
				mapMaker.plotXYZData(pDiff, pDiffCPT, name1+" vs "+name2+", % Change, "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_rate_haz_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_rate_haz_pDiff.png)");
				mapMaker.plotXYZData(diff, diffCPT, name1+" - "+name2+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_rate_haz_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_rate_haz_diff.png)");
				table.finalizeLine();
				
				lines.addAll(table.build());
				lines.add("");
			}
			
			// moment maps
			lines.add("### Moment Maps, "+label);
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine(MarkdownUtils.boldCentered(name1), MarkdownUtils.boldCentered(name2));
			
			table.initNewLine();
			
			mapMaker.plotXYZData(asLog10(moRates1), moCPT, "Log10 Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, prefix+"_moment_rate1", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate1.png)");
			
			mapMaker.plotXYZData(asLog10(moRates2), moCPT, "Log10 Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, prefix+"_moment_rate2", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate2.png)");
			
			table.finalizeLine();
			
			table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
			
			table.initNewLine();
			
			GriddedGeoDataSet moPDiff = WUS_HazardChangePageGen.mapPDiff(zerosToNaNs(moRates1), zerosToNaNs(moRates2));
			GriddedGeoDataSet moDiff = WUS_HazardChangePageGen.mapDiff(zerosToNaNs(moRates1), zerosToNaNs(moRates2));
			
			mapMaker.plotXYZData(moPDiff, pDiff100CPT, "% Change in Moment Rate");
			mapMaker.plot(resourcesDir, prefix+"_moment_rate_pDiff", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate_pDiff.png)");
			
			mapMaker.plotXYZData(moDiff, moRateDiffCPT, "Moment Rate Difference (N-m)");
			mapMaker.plot(resourcesDir, prefix+"_moment_rate_diff", " ");
			table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate_diff.png)");
			
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			if (hazard) {
				// rate-based simplified hazard maps
				lines.add("### Moment Rate-Based Simplified Hazard Maps, "+label);
				lines.add(topLink); lines.add("");
				
				GutenbergRichterMagFreqDist mfdShape = new GutenbergRichterMagFreqDist(6d, 4, 0.5, 1e16, 1d);
				
				GriddedGeoDataSet smoothHaz1 = SimpleSmoothHazardMapCalc.hazMapMomentSmooth(
						moRates1, gmpeRef, mfdShape, period, rp, smoothHazGridReg);
				GriddedGeoDataSet smoothHaz2 = SimpleSmoothHazardMapCalc.hazMapMomentSmooth(
						moRates2, gmpeRef, mfdShape, period, rp, smoothHazGridReg);
				
				table = MarkdownUtils.tableBuilder();
				
				table.addLine(name1, name2);
				
				table.initNewLine();
				mapMaker.plotXYZData(asLog10(smoothHaz1), hazCPT, name1+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_moment_rate_haz1", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate_haz1.png)");
				mapMaker.plotXYZData(asLog10(smoothHaz2), hazCPT, name2+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_moment_rate_haz2", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate_haz2.png)");
				table.finalizeLine();
				
				GriddedGeoDataSet pDiff = WUS_HazardChangePageGen.mapPDiff(smoothHaz1, smoothHaz2);
				GriddedGeoDataSet diff = WUS_HazardChangePageGen.mapDiff(smoothHaz1, smoothHaz2);
				
				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				table.initNewLine();
				mapMaker.plotXYZData(pDiff, pDiffCPT, name1+" vs "+name2+", % Change, "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_moment_rate_haz_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate_haz_pDiff.png)");
				mapMaker.plotXYZData(diff, diffCPT, name1+" - "+name2+", "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_moment_rate_haz_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_moment_rate_haz_diff.png)");
				table.finalizeLine();
				
				lines.addAll(table.build());
				lines.add("");
			}
			
			if (writeAsWeGo) {
				// write intermediate
				List<String> lines2 = new ArrayList<>(lines);
				// add TOC
				lines2.addAll(tocIndex, MarkdownUtils.buildTOC(lines2, 2));
				lines2.add(tocIndex, "## Table Of Contents");

				// write markdown
				MarkdownUtils.writeReadmeAndHTML(lines2, outputDir);
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static GriddedGeoDataSet loadExtMap(File file, GriddedRegion gridReg) throws IOException {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(gridReg);
		
		for (int i=0; i<ret.size(); i++)
			ret.set(i, Double.NaN);
		
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith("#"))
				continue;
			String[] split = line.split(",");
			Preconditions.checkState(split.length == 3);
			Location loc = new Location(Double.parseDouble(split[1]), Double.parseDouble(split[0]));
			double val = Double.parseDouble(split[2]);
			
			int index = gridReg.indexForLocation(loc);
			if (index >= 0)
				ret.set(index, val);
		}
		
		return ret;
	}
	
	static GriddedGeoDataSet asLog10(GriddedGeoDataSet xyz) {
		xyz = zerosToNaNs(xyz);
		xyz.log10();
		return xyz;
	}
	
	static GriddedGeoDataSet zerosToNaNs(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		for (int i=0; i<xyz.size(); i++)
			if (xyz.get(i) == 0d)
				xyz.set(i, Double.NaN);
		return xyz;
	}

}