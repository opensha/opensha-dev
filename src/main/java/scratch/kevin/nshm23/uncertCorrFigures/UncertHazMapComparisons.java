package scratch.kevin.nshm23.uncertCorrFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.nshm23.figures.MethodsAndIngredientsHazChangeFigures;

public class UncertHazMapComparisons {

	public static void main(String[] args) throws IOException {
		boolean redownload = false;
		String baseURL = "https://data.opensha.org/ftp/kmilner/markdown/batch_inversions/";
		File invsDir = new File("/data/kevin/nshm23/batch_inversions/");
		
		File outputDir = new File("/home/kevin/Documents/papers/2024_nshm23_uncert_correlation/figures/hazard_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		String refDir = "2024_02_02-nshm23_branches-WUS_FM_v3/";
		String refHazardDirName = "hazard_maps_full_gridded";
		
		String randHazardDirName = "hazard_maps_full_gridded_comp_correlated";
		String randSegBDir = "2023_11_16-nshm23_branches-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
		String randDMDir = "2023_11_17-nshm23_branches-dm_sampling-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
		String randDMSegBDir = "2023_11_20-nshm23_branches-dm_sampling-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
		
		GriddedRegion gridReg = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.1d, GriddedRegion.ANCHOR_0_0);
		
		String csvFileName = "pga_TWO_IN_50.csv";
		String imtLabel = "PGA (g), 2% in 50 years";
		String imtUnitlessLabel = "PGA, 2% in 50 years";
//		2023_11_16-nshm23_branches-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/hazard_maps_full_gridded_comp_correlated/resources/pga_TWO_IN_50.csv";
		
		MapCSV refMaps = loadCSV(invsDir, baseURL, refDir, refHazardDirName, csvFileName, gridReg, redownload);
		MapCSV randSegBMaps = loadCSV(invsDir, baseURL, randSegBDir, randHazardDirName, csvFileName, gridReg, redownload);
		MapCSV randDMMaps = loadCSV(invsDir, baseURL, randDMDir, randHazardDirName, csvFileName, gridReg, redownload);
		MapCSV randDMSegBMaps = loadCSV(invsDir, baseURL, randDMSegBDir, randHazardDirName, csvFileName, gridReg, redownload);
		
		boolean covToAverage = false;
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(NSHM23_RegionLoader.loadFullConterminousWUS());
		mapMaker.setFaultSections(NSHM23_FaultModels.WUS_FM_v3.getFaultSections());
		mapMaker.setSectOutlineChar(null);
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(0, 0, 0, 80)));
		
		CPT covCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0, 1d);
		covCPT.setNanColor(Color.LIGHT_GRAY);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT maskedPDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(pDiffCPT, 5d, 50d);
		maskedPDiffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.1d, 0.1d);
		diffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT covDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.3d, 0.3d);
		covDiffCPT.setNanColor(Color.LIGHT_GRAY);
		
		int refIndex = 0;
		MapCSV[] models = {
				refMaps,
				randSegBMaps,
				randDMMaps,
				randDMSegBMaps
		};
		String[] names = {
				"Fully Correlated",
				"Random b/Seg",
				"Random DM",
				"Random DM/b/Seg"
		};
		String[] prefixes = {
				"correlated",
				"rand_seg_b",
				"rand_dm",
				"rand_dm_seg_b"
		};
		
		for (int i=0; i<models.length; i++) {
			GriddedGeoDataSet cov = models[i].covMap;
			mapMaker.plotXYZData(cov, covCPT, names[i]+", COV, "+imtUnitlessLabel);
			mapMaker.plot(outputDir, "cov_"+prefixes[i], " ");
			
			if (i != refIndex) {
				if (covToAverage)
					// recalculate the COV using the average mean of each model
					cov = calcCovToAverageMean(models[i], refMaps);
				GriddedGeoDataSet covDiff = buildDiff(refMaps.covMap, cov);
				GriddedGeoDataSet covPDiff = buildPDiff(refMaps.covMap, cov);
				
				mapMaker.plotXYZData(covDiff, covDiffCPT, names[i]+", COV Difference, "+imtUnitlessLabel);
				mapMaker.plot(outputDir, "cov_diff_"+prefixes[i]+"_"+prefixes[refIndex], " ");
				mapMaker.plotXYZData(covPDiff, maskedPDiffCPT, names[i]+", COV % Change, "+imtUnitlessLabel);
				mapMaker.plot(outputDir, "cov_pDiff_"+prefixes[i]+"_"+prefixes[refIndex], " ");
				
				GriddedGeoDataSet diff = buildDiff(refMaps.meanMap, models[i].meanMap);
				GriddedGeoDataSet pDiff = buildPDiff(refMaps.meanMap, models[i].meanMap);
				
				mapMaker.plotXYZData(diff, diffCPT, names[i]+", Difference, "+imtLabel);
				mapMaker.plot(outputDir, "mean_diff_"+prefixes[i]+"_"+prefixes[refIndex], " ");
				mapMaker.plotXYZData(pDiff, maskedPDiffCPT, names[i]+", % Change, "+imtUnitlessLabel);
				mapMaker.plot(outputDir, "mean_pDiff_"+prefixes[i]+"_"+prefixes[refIndex], " ");
			}
		}
	}
	
	private static class MapCSV {
		
		public final GriddedGeoDataSet meanMap;
		public final GriddedGeoDataSet medianMap;
		public final GriddedGeoDataSet minMap;
		public final GriddedGeoDataSet maxMap;
		public final GriddedGeoDataSet sdMap;
		public final GriddedGeoDataSet covMap;
		
		public MapCSV(GriddedRegion gridReg, File csvFile) throws IOException {
			CSVFile<String> csv = CSVFile.readFile(csvFile, true);
			
			meanMap = new GriddedGeoDataSet(gridReg);
			medianMap = new GriddedGeoDataSet(gridReg);
			minMap = new GriddedGeoDataSet(gridReg);
			maxMap = new GriddedGeoDataSet(gridReg);
			sdMap = new GriddedGeoDataSet(gridReg);
			covMap = new GriddedGeoDataSet(gridReg);
			
			for (int row=1; row<csv.getNumRows(); row++) {
				int index = csv.getInt(row, 0);
				double lat = csv.getDouble(row, 1);
				double lon = csv.getDouble(row, 2);
				Preconditions.checkState(LocationUtils.areSimilar(gridReg.getLocation(index), new Location(lat, lon)));
				double mean = csv.getDouble(row, 3);
				double median = csv.getDouble(row, 4);
				double min = csv.getDouble(row, 5);
				double max = csv.getDouble(row, 6);
				double cov = csv.getDouble(row, 7);
				
				meanMap.set(index, mean);
				medianMap.set(index, median);
				minMap.set(index, min);
				maxMap.set(index, max);
				covMap.set(index, cov);
				sdMap.set(index, cov*mean);
			}
		}
	}
	
	private static GriddedGeoDataSet calcCovToAverageMean(MapCSV map, MapCSV other) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map.meanMap.getRegion());
		
		for (int i=0; i<ret.size(); i++) {
			double mean = 0.5*(map.meanMap.get(i) + other.meanMap.get(i));
			double cov = map.sdMap.get(i) / mean;
			ret.set(i, cov);
		}
		
		return ret;
	}
	
	private static MapCSV loadCSV(File baseDir, String baseURL, String dirName, String hazardDirName,
			String mapFileName, GriddedRegion gridReg, boolean redownload) throws IOException {
		String mapURL = mapURL(baseURL, dirName, hazardDirName, mapFileName);
		File mapFile = mapFile(baseDir, dirName, hazardDirName, mapFileName);
		return loadCSV(mapURL, mapFile, gridReg, redownload);
	}
	
	private static MapCSV loadCSV(String url, File file, GriddedRegion gridReg, boolean redownload) throws IOException {
		if (!file.exists() || redownload)
			// need to downoad
			FileUtils.downloadURL(url, file);
		return new MapCSV(gridReg, file);
	}
	
	private static String mapURL(String baseURL, String dirName, String hazardDirName, String mapFileName) {
		if (!baseURL.endsWith("/"))
			baseURL += "/";
		return baseURL+dirName+"/"+hazardDirName+"/resources/"+mapFileName;
	}
	
	private static File mapFile(File baseDir, String dirName, String hazardDirName, String mapFileName) {
		File invDir = new File(baseDir, dirName);
		Preconditions.checkState(invDir.exists(), "Doesn't exist: %s", invDir.getAbsolutePath());
		
		File hazardDir = new File(invDir, hazardDirName);
		Preconditions.checkState(hazardDir.exists() || hazardDir.mkdir());
		
		File resourcesDir = new File(hazardDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		return new File(resourcesDir, mapFileName);
	}
	
	private static GriddedGeoDataSet buildPDiff(GriddedGeoDataSet ref, GriddedGeoDataSet comp) {
		GriddedGeoDataSet diff = new GriddedGeoDataSet(ref.getRegion(), false);
		
		for (int i=0; i<ref.size(); i++) {
			double compVal = comp.get(i);
			double refVal = ref.get(i);
			double pDiff = 100d*(compVal-refVal)/refVal;
			diff.set(i, pDiff);
		}
		
		return diff;
	}
	
	private static GriddedGeoDataSet buildDiff(GriddedGeoDataSet ref, GriddedGeoDataSet comp) {
		GriddedGeoDataSet diff = new GriddedGeoDataSet(ref.getRegion(), false);
		
		for (int i=0; i<ref.size(); i++) {
			double compVal = comp.get(i);
			double refVal = ref.get(i);
			diff.set(i, compVal-refVal);
		}
		
		return diff;
	}

}
