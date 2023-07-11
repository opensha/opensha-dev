package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.CubedGriddedRegion;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_FaultCubeAssociations;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class MomentRateCompNSHM18 {
	
	public static boolean LINEAR_RAMP = true;
	public static boolean GEO_ONLY = false;

	public static void main(String[] args) throws IOException {
		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		region = new Region(new Location(region.getMinLat(), region.getMinLon()),
				new Location(region.getMaxLat(), region.getMaxLon()));
		GriddedRegion gridReg = new GriddedRegion(region, 0.1, GriddedRegion.ANCHOR_0_0);
		
		File outputDir = new File("/home/kevin/markdown/nshm23-misc/fault_mo_rates");
		if (GEO_ONLY)
			outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_geo");
		if (LINEAR_RAMP)
			outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_linear_ramp");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		CubedGriddedRegion cgr = new CubedGriddedRegion(gridReg);
		
		int width = 1000;
		
		// calculate NSHM18 first
		GriddedGeoDataSet nshm18 = new GriddedGeoDataSet(gridReg, false);
		
		// NSHM18 outside of CA
		NSHM18_DeformationModels dm18 = GEO_ONLY ? NSHM18_DeformationModels.GEOL : NSHM18_DeformationModels.BRANCH_AVERAGED;
		List<? extends FaultSection> nshm18sectsNoCA = dm18.build(NSHM18_FaultModels.NSHM18_WUS_NoCA);
		addFaultMoRates(nshm18, cgr, nshm18sectsNoCA, 1d);
		
		// UCERF3
		List<FaultSection> allSects18 = new ArrayList<>(nshm18sectsNoCA);
		DeformationModels u3DM = GEO_ONLY ? DeformationModels.GEOLOGIC : DeformationModels.MEAN_UCERF3;
		for (FaultModels u3fm : FaultModels.values()) {
			double fmWeight = u3fm.getNodeWeight(null);
			if (fmWeight == 0d)
				continue;
			
			List<? extends FaultSection> fullDM = u3DM.build(u3fm);
			
			List<FaultSection> modDM = new ArrayList<>();
			for (FaultSection sect : fullDM) {
				if (sect.getParentSectionId() == 721 || sect.getParentSectionId() == 719)
					continue;
				sect.setSectionId(modDM.size());
				modDM.add(sect);
			}
			
			addFaultMoRates(nshm18, cgr, modDM, fmWeight);
			allSects18.addAll(modDM);
		}
		
		// now NSHM23
		GriddedGeoDataSet nshm23 = new GriddedGeoDataSet(gridReg, false);
		
		NSHM23_DeformationModels dm23 = GEO_ONLY ? NSHM23_DeformationModels.GEOLOGIC : NSHM23_DeformationModels.AVERAGE;
		List<? extends FaultSection> nshm23sects = dm23.build(NSHM23_FaultModels.NSHM23_v2);
		
		addFaultMoRates(nshm23, cgr, nshm23sects, 1d);
		
		System.out.println("Total NSHM18 moment rate: "+(float)nshm18.getSumZ());
		System.out.println("Total NSHM23 moment rate: "+(float)nshm23.getSumZ());
		
		double minMoRate = Double.POSITIVE_INFINITY;
		double maxMoRate = Double.NEGATIVE_INFINITY;
		GriddedGeoDataSet diff = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet pDiff = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet ratio = new GriddedGeoDataSet(gridReg, false);
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Location Index", "Latitutde", "Longitude", "NSHM23 Moment Rate (N-m)", "NSHM18 Moment Rate (N-m)");
		
		double maxDiff = 0d;
		for (int i=0; i<nshm23.size(); i++) {
			double v1 = nshm23.get(i);
			double v2 = nshm18.get(i);
			if (v1 > 0)
				minMoRate = Math.min(minMoRate, v1);
			if (v2 > 0)
				minMoRate = Math.min(minMoRate, v2);
			maxMoRate = Math.max(maxMoRate, v1);
			maxMoRate = Math.max(maxMoRate, v2);
			
			diff.set(i, v1 - v2);
			maxDiff = Math.max(maxDiff, Math.abs(v1 - v2));
			ratio.set(i, v1/v2);
			pDiff.set(i, 100d*(v1 - v2)/v2);
			
			List<String> line = new ArrayList<>();
			line.add(i+"");
			Location loc = gridReg.getLocation(i);
			line.add((float)loc.getLatitude()+"");
			line.add((float)loc.getLongitude()+"");
			line.add((float)v1+"");
			line.add((float)v2+"");
			csv.addLine(line);
		}
		
		double logMaxMo = Math.ceil(Math.log10(maxMoRate));
		double logMinMo = Math.max(logMaxMo-10d, Math.floor(Math.log10(minMoRate)));
		
		CPT logMoCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMinMo, logMaxMo);
		logMoCPT.setNanColor(Color.WHITE);
		logMoCPT.setBelowMinColor(Color.WHITE);
		
		CPT linearMoCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0, maxMoRate);
		linearMoCPT.setNanColor(Color.WHITE);
		linearMoCPT.setBelowMinColor(Color.WHITE);
		
		maxDiff = Math.pow(10, Math.ceil(Math.min(16, Math.log10(maxDiff))));
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-maxDiff, maxDiff);
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		
		CPT logRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
//		logRatioCPT.setAboveMaxColor(Color.BLACK);
//		logRatioCPT.setBelowMinColor(Color.GREEN);
//		logRatioCPT.setAboveMaxColor(diffCPT.getMaxColor().brighter());
//		logRatioCPT.setBelowMinColor(diffCPT.getMinColor().brighter());
		logRatioCPT.setAboveMaxColor(Color.YELLOW);
		logRatioCPT.setBelowMinColor(Color.GREEN);
		
		GriddedGeoDataSet logNSHM23 = nshm23.copy();
		logNSHM23.log10();
		GriddedGeoDataSet logNSHM18 = nshm18.copy();
		logNSHM18.log10();
		
		GeographicMapMaker mapMaker23 = new RupSetMapMaker(nshm23sects, region);
		mapMaker23.setSectOutlineChar(null);
		PlotCurveCharacterstics sectChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(0, 0, 00, 127));
		mapMaker23.setSectTraceChar(sectChar);
		// re-index and remove dups
		HashSet<String> prevNames = new HashSet<>();
		for (int i=allSects18.size(); --i>=0;) {
			String name = allSects18.get(i).getSectionName();
			if (prevNames.contains(name))
				allSects18.remove(i);
			else
				prevNames.add(name);
		}
		for (int i=0; i<allSects18.size(); i++)
			allSects18.get(i).setSectionId(i);
		GeographicMapMaker mapMaker18 = new RupSetMapMaker(allSects18, region);
		mapMaker18.setSectOutlineChar(null);
		mapMaker18.setSectTraceChar(sectChar);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# NSHM23 vs NSHM18 Fault Moment Comparison"+(GEO_ONLY ? ", Geologic DM Only" : ""));
		lines.add("");
		
		File csvFile = new File(outputDir, "moment_rates.csv");
		csv.writeToFile(csvFile);
		lines.add("Download CSV: ["+csvFile.getName()+"]("+csvFile.getName()+")");
		lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("NSHM23", "NSHM18");
		table.initNewLine();
		mapMaker23.plotXYZData(zeroAsNan(nshm23), linearMoCPT, "Fault Moment Rate (N-m)");
		mapMaker23.plot(outputDir, "nshm23", " ", width);
		table.addColumn("![NSHM23](nshm23.png)");
		mapMaker18.plotXYZData(zeroAsNan(nshm18), linearMoCPT, "Fault Moment Rate (N-m)");
		mapMaker18.plot(outputDir, "nshm18", " ", width);
		table.addColumn("![NSHM18](nshm18.png)");
		table.finalizeLine();
		table.initNewLine();
		GriddedGeoDataSet.writeXYZFile(nshm23, new File(outputDir, "nshm23.xyz"));
		table.addColumn("[Download XYZ File](nshm23.xyz)");
		GriddedGeoDataSet.writeXYZFile(nshm18, new File(outputDir, "nshm18.xyz"));
		table.addColumn("[Download XYZ File](nshm18.xyz)");
		table.finalizeLine();
		table.initNewLine();
		mapMaker23.plotXYZData(logNSHM23, logMoCPT, "Log10(Fault Moment Rate) (N-m)");
		mapMaker23.plot(outputDir, "nshm23_log", " ", width);
		table.addColumn("![NSHM23](nshm23_log.png)");
		mapMaker18.plotXYZData(logNSHM18, logMoCPT, "Log10(Fault Moment Rate) (N-m)");
		mapMaker18.plot(outputDir, "nshm18_log", " ", width);
		table.addColumn("![NSHM18](nshm18_log.png)");
		table.finalizeLine();
		
		lines.addAll(table.build());
		lines.add("");
		
		
		// diff and ratios
		table = MarkdownUtils.tableBuilder();
		table.addLine("Difference", "Ratio");
		table.initNewLine();
		mapMaker23.plotXYZData(zeroAsNan(diff), diffCPT, "NSHM23 - NSHM18 (N-m)");
		mapMaker23.plot(outputDir, "diff", " ", width);
		table.addColumn("![Difference](diff.png)");
		
		// build log ratio, with saturation
		GriddedGeoDataSet logRatio = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<logRatio.size(); i++) {
			double val = ratio.get(i);
			double logVal = Math.log10(val);
			if (Double.isFinite(logVal)) {
				// neither one was zero
				// check to see if we should saturate it
				if ((float)logVal > logRatioCPT.getMaxValue())
					logVal = (double)logRatioCPT.getMaxValue();
				else if ((float)logVal < logRatioCPT.getMinValue())
					logVal = (double)logRatioCPT.getMinValue();
			}
			logRatio.set(i, logVal);
		}
		mapMaker23.plotXYZData(logRatio, logRatioCPT, "Log10(NSHM23 / NSHM18)");
		mapMaker23.plot(outputDir, "log_ratio", " ", width);
		table.addColumn("![Ratio](log_ratio.png)");
		table.finalizeLine();
		// pdiff
		mapMaker23.plotXYZData(pDiff, pDiffCPT, "NSHM23 vs NSHM18, % Difference");
		mapMaker23.plot(outputDir, "pdiff", " ", width);
		table.addLine(MarkdownUtils.boldCentered("% Difference"), "");
		table.addLine("![% Diff](pdiff.png)", "");
		lines.addAll(table.build());
		lines.add("");
		
		// plot slip rates
		CPT slipsCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 40d);
		CPT logSlipsCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 2);
		
		mapMaker23.clearXYZData();
		mapMaker18.clearXYZData();
		
		mapMaker23.setWriteGeoJSON(true);
		mapMaker18.setWriteGeoJSON(true);
		
		table = MarkdownUtils.tableBuilder();
		table.addLine("NSHM23 Slip Rates", "NSHM18 Slip Rates");
		table.initNewLine();
		mapMaker23.plotSectScalars(slipRates(nshm23sects), slipsCPT, "Slip Rate (mm/yr)");
		mapMaker23.plot(outputDir, "nshm23_slips", " ", width);
		table.addColumn("![NSHM23](nshm23_slips.png)");
		mapMaker18.plotSectScalars(slipRates(allSects18), slipsCPT, "Slip Rate (mm/yr)");
		mapMaker18.plot(outputDir, "nshm18_slips", " ", width);
		table.addColumn("![NSHM23](nshm18_slips.png)");
		table.finalizeLine();
		table.initNewLine();
		table.addColumn(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON", "nshm23_slips.geojson")
				+" "+"[Download GeoJSON](nshm23_slips.geojson)");
		table.addColumn(RupSetMapMaker.getGeoJSONViewerRelativeLink("View GeoJSON", "nshm18_slips.geojson")
				+" "+"[Download GeoJSON](nshm18_slips.geojson)");
		table.finalizeLine();
		table.initNewLine();
		mapMaker23.plotSectScalars(logSlipRates(nshm23sects), logSlipsCPT, "Log10 Slip Rate (mm/yr)");
		mapMaker23.plot(outputDir, "nshm23_slips_log", " ", width);
		table.addColumn("![NSHM23](nshm23_slips_log.png)");
		mapMaker18.plotSectScalars(logSlipRates(allSects18), logSlipsCPT, "Log10 Slip Rate (mm/yr)");
		mapMaker18.plot(outputDir, "nshm18_slips_log", " ", width);
		table.addColumn("![NSHM23](nshm18_slips_log.png)");
		table.finalizeLine();
		
		lines.addAll(table.build());
		lines.add("");
		
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
		
//		CPT moCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMinMo, logMaxMo);
//		moCPT.setNanColor(Color.WHITE);
	}
	
	private static double[] slipRates(List<? extends FaultSection> subSects) {
		double[] ret = new double[subSects.size()];
		for (int s=0; s<ret.length; s++)
			ret[s] = subSects.get(s).getOrigAveSlipRate();
		return ret;
	}
	
	private static double[] logSlipRates(List<? extends FaultSection> subSects) {
		double[] ret = new double[subSects.size()];
		for (int s=0; s<ret.length; s++)
			ret[s] = Math.log10(subSects.get(s).getOrigAveSlipRate());
		return ret;
	}
	
	private static GriddedGeoDataSet zeroAsNan(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		for (int i=0; i<xyz.size(); i++)
			if (xyz.get(i) == 0d)
				xyz.set(i, Double.NaN);
		return xyz;
	}
	
	private static void addFaultMoRates(GriddedGeoDataSet xyz, CubedGriddedRegion cgr,
			List<? extends FaultSection> sects, double scale) {
		FaultGridAssociations assoc;
		if (LINEAR_RAMP)
			assoc = new NSHM23_FaultCubeAssociations(sects, cgr, NSHM23_SingleRegionGridSourceProvider.DEFAULT_MAX_FAULT_NUCL_DIST);
		else
			assoc = FaultGridAssociations.getIntersectionAssociations(sects, xyz.getRegion());
		
		for (int s=0; s<sects.size(); s++) {
			FaultSection sect = sects.get(s);
			double area = sect.getArea(false); // already in m^2
			double slipRate = sect.getOrigAveSlipRate()*1e-3; // mm/yr -> m/yr
			double moRate = FaultMomentCalc.getMoment(area, slipRate);
			Map<Integer, Double> nodeFracts = assoc.getNodeFractions(s);
			for (int nodeIndex : nodeFracts.keySet()) {
				double fract = nodeFracts.get(nodeIndex)*scale;
				xyz.set(nodeIndex, xyz.get(nodeIndex) + fract*moRate);
			}
		}
	}
	
	public static GriddedGeoDataSet getMomentRatesNSHM18(GriddedRegion gridReg, CubedGriddedRegion cgr) throws IOException {
		return getMomentRatesNSHM18(gridReg, cgr, null);
	}
	
	public static GriddedGeoDataSet getMomentRatesNSHM18(GriddedRegion gridReg, CubedGriddedRegion cgr,
			List<FaultSection> sectsList) throws IOException {
		GriddedGeoDataSet nshm18 = new GriddedGeoDataSet(gridReg, false);

		// NSHM18 outside of CA
		NSHM18_DeformationModels dm18 = GEO_ONLY ? NSHM18_DeformationModels.GEOL : NSHM18_DeformationModels.BRANCH_AVERAGED;
		List<? extends FaultSection> nshm18sectsNoCA = dm18.build(NSHM18_FaultModels.NSHM18_WUS_NoCA);
		addFaultMoRates(nshm18, cgr, nshm18sectsNoCA, 1d);

		// UCERF3
		if (sectsList != null)
			sectsList.addAll(nshm18sectsNoCA);
		DeformationModels u3DM = GEO_ONLY ? DeformationModels.GEOLOGIC : DeformationModels.MEAN_UCERF3;
		for (FaultModels u3fm : FaultModels.values()) {
			double fmWeight = u3fm.getNodeWeight(null);
			if (fmWeight == 0d)
				continue;

			List<? extends FaultSection> fullDM = u3DM.build(u3fm);

			List<FaultSection> modDM = new ArrayList<>();
			for (FaultSection sect : fullDM) {
				if (sect.getParentSectionId() == 721 || sect.getParentSectionId() == 719)
					continue;
				sect.setSectionId(modDM.size());
				modDM.add(sect);
			}

			addFaultMoRates(nshm18, cgr, modDM, fmWeight);
			if (sectsList != null)
				sectsList.addAll(modDM);
		}
		
		return nshm18;
	}
	
	public static GriddedGeoDataSet getMomentRatesNSHM23(GriddedRegion gridReg, CubedGriddedRegion cgr) throws IOException {
		GriddedGeoDataSet nshm23 = new GriddedGeoDataSet(gridReg, false);
		
		NSHM23_DeformationModels dm23 = GEO_ONLY ? NSHM23_DeformationModels.GEOLOGIC : NSHM23_DeformationModels.AVERAGE;
		List<? extends FaultSection> nshm23sects = dm23.build(NSHM23_FaultModels.NSHM23_v2);
		
		addFaultMoRates(nshm23, cgr, nshm23sects, 1d);
		
		return nshm23;
	}

}
