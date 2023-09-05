package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.MathArrays;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.hazard.LogicTreeCurveAverager;
import org.opensha.sha.earthquake.faultSysSolution.hazard.LogicTreeHazardCompare;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

class CA_HazardDistributionFigures {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/u3_haz_distribution_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File modelDir = new File(invsDir, "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		LogicTree<?> modelTree = LogicTree.read(new File(modelDir, "logic_tree.json"));
		ZipFile modelZip = new ZipFile(new File(modelDir, "results_hazard_avg_gridded.zip"));
		FaultSystemRupSet rupSet23 = FaultSystemRupSet.load(new File(modelDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		File u3Dir = new File(invsDir, "2023_05_08-u3-full_dist-nshm23_gridded");
		LogicTree<?> u3Tree = LogicTree.read(new File(u3Dir, "logic_tree.json"));
		ZipFile u3Zip = new ZipFile(new File(u3Dir, "results_hazard_avg_gridded.zip"));
		
		File methodologyDir = new File(invsDir, "2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		LogicTree<?> methodologyTree = LogicTree.read(new File(methodologyDir, "logic_tree.json"));
		ZipFile methodologyZip = new ZipFile(new File(methodologyDir, "results_hazard_avg_nshm23_gridded.zip"));
		FaultSystemRupSet rupSetU3 = FaultSystemRupSet.load(new File(methodologyDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		
		String entryName = "map_pga_TWO_IN_50.txt";
		String hazLabel = "PGA, 2% in 50 yrs";
		String units = "(g)";
		
		ZipEntry regEntry = u3Zip.getEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME);
		System.out.println("Reading gridded region from zip file: "+regEntry.getName());
		BufferedReader bRead = new BufferedReader(new InputStreamReader(u3Zip.getInputStream(regEntry)));
		GriddedRegion gridReg = GriddedRegion.fromFeature(Feature.read(bRead));
		
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		System.out.println("Loading UCERF3");
		LightFixedXFunc[] u3NormCDFs = new LightFixedXFunc[gridReg.getNodeCount()];
		double[] u3StdDevs = new double[gridReg.getNodeCount()];
		GriddedGeoDataSet u3Mean = loadMaps(u3Zip, u3Tree, gridReg, entryName, u3NormCDFs, u3StdDevs, exec);
		
		System.out.println("Loading Ingredients");
		LightFixedXFunc[] methodologyNormCDFs = new LightFixedXFunc[gridReg.getNodeCount()];
		double[] methodologyStdDevs = new double[gridReg.getNodeCount()];
		GriddedGeoDataSet methodologyMean = loadMaps(methodologyZip, methodologyTree, gridReg, entryName, methodologyNormCDFs, methodologyStdDevs, exec);
		
		System.out.println("Loading Model");
		LightFixedXFunc[] modelNormCDFs = new LightFixedXFunc[gridReg.getNodeCount()];
		double[] modelStdDevs = new double[gridReg.getNodeCount()];
		GriddedGeoDataSet modelMean = loadMaps(modelZip, modelTree, gridReg, entryName, modelNormCDFs, modelStdDevs, exec);
		
		GeographicMapMaker mapMakerU3 = new RupSetMapMaker(rupSetU3, gridReg);
		GeographicMapMaker mapMaker23 = new RupSetMapMaker(rupSet23, gridReg);
		
		for (GeographicMapMaker mapMaker : new GeographicMapMaker[] {mapMakerU3, mapMaker23}) {
			mapMaker.setSectOutlineChar(null);
			mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(100, 100, 100, 127)));
		}
		
		CPT rainbow_0_1 = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 1d);
		rainbow_0_1.setNanColor(new Color(255, 255, 255, 0));
		CPT diffCPT_0p2 = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2, 0.2);
		diffCPT_0p2.setNanColor(new Color(255, 255, 255, 0));
		CPT diffCPT_0p5 = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.5, 0.5);
		diffCPT_0p5.setNanColor(new Color(255, 255, 255, 0));
		
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		Region[] caRegions = new Region[caOutlines.length];
		for (int i=0; i<caOutlines.length; i++) {
			LocationList outline = new LocationList();
			for (Point2D pt : caOutlines[i])
				outline.add(new Location(pt.getY(), pt.getX()));
			caRegions[i] = new Region(outline, BorderType.MERCATOR_LINEAR);
		}
		
		GriddedGeoDataSet u3COV = calcCOV(u3Mean, u3StdDevs);
		GriddedGeoDataSet modelCOV = calcCOV(modelMean, modelStdDevs);
		GriddedGeoDataSet methodologyCOV = calcCOV(methodologyMean, methodologyStdDevs);
		
		plotMap(outputDir, "cov_u3", mapMakerU3, caRegions, u3COV, rainbow_0_1, "UCERF3 COV");
		plotMap(outputDir, "cov_methodology", mapMakerU3, caRegions, methodologyCOV, rainbow_0_1, "NSHM23 Methodology COV");
		plotMap(outputDir, "cov_model", mapMaker23, caRegions, modelCOV, rainbow_0_1, "NSHM23 COV");
		
		plotMap(outputDir, "cov_diff_methodolgy_u3", mapMakerU3, caRegions, diff(methodologyCOV, u3COV), diffCPT_0p2,
				"NSHM23 Methodology COV - UCERF3 COV");
		plotMap(outputDir, "cov_diff_model_u3", mapMaker23, caRegions, diff(modelCOV, u3COV), diffCPT_0p2,
				"NSHM23 COV - UCERF3 COV");
		
		GriddedGeoDataSet u3Spread = logSpreadCOV(gridReg, u3NormCDFs);
		GriddedGeoDataSet modelSpread = logSpreadCOV(gridReg, modelNormCDFs);
		GriddedGeoDataSet methodologySpread = logSpreadCOV(gridReg, methodologyNormCDFs);
		
		plotMap(outputDir, "spread_u3", mapMakerU3, caRegions, u3Spread, rainbow_0_1, "UCERF3 Log10(max/min)");
		plotMap(outputDir, "spread_methodology", mapMakerU3, caRegions, methodologySpread, rainbow_0_1, "NSHM23 Methodology Log10(max/min)");
		plotMap(outputDir, "spread_model", mapMaker23, caRegions, modelSpread, rainbow_0_1, "NSHM23 Log10(max/min)");
		
		plotMap(outputDir, "spread_diff_methodolgy_u3", mapMakerU3, caRegions, diff(methodologySpread, u3Spread), diffCPT_0p5,
				"NSHM23 Methodology - UCERF3, Log10(max/min)");
		plotMap(outputDir, "spread_diff_model_u3", mapMaker23, caRegions, diff(modelSpread, u3Spread), diffCPT_0p5,
				"NSHM23 - UCERF3, Log10(max/min)");
		
		CPT percentileCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 100d);
		percentileCPT.setBelowMinColor(Color.LIGHT_GRAY);
		percentileCPT.setAboveMaxColor(Color.BLACK);
		percentileCPT.setNanColor(new Color(255, 255, 255, 0));
		
		GriddedGeoDataSet u3WithinMethodology = percentile(u3Mean, methodologyNormCDFs);
		GriddedGeoDataSet u3WithinModel = percentile(u3Mean, modelNormCDFs);
		
		plotMap(outputDir, "percentile_u3_within_methodology", mapMakerU3, caRegions, u3WithinMethodology,
				percentileCPT, "UCERF3 %-ile Within NSHM23 Methodology");
		plotMap(outputDir, "percentile_u3_within_model", mapMaker23, caRegions, u3WithinModel,
				percentileCPT, "UCERF3 %-ile Within NSHM23");
		
		exec.shutdown();
	}
	
	private static GriddedGeoDataSet loadMaps(ZipFile zip, LogicTree<?> tree, GriddedRegion reg, String mapEntryName,
			LightFixedXFunc[] normCDFs, double[] stdDevs, ExecutorService exec) throws IOException {
		List<GriddedGeoDataSet> maps = new ArrayList<>(tree.size());
		List<Double> weights = new ArrayList<>();
		double weightSum = 0d;
		for (LogicTreeBranch<?> branch : tree) {
			
			String entryName = branch.buildFileName()+"/"+mapEntryName;

			double weight = tree.getBranchWeight(branch);
			weightSum += weight;
			weights.add(weight);
			maps.add(readMap(zip, entryName, reg));
		}
		System.out.println("Loaded "+maps.size()+" maps");
		Preconditions.checkState(normCDFs.length == reg.getNodeCount());
		Preconditions.checkState(normCDFs.length == stdDevs.length);
		List<Future<CDF_Callable>> futures = new ArrayList<>(reg.getNodeCount());
		
		double[] weightsArray = MathArrays.normalizeArray(Doubles.toArray(weights), weights.size());
		for (int i=0; i<normCDFs.length; i++)
			futures.add(exec.submit(new CDF_Callable(maps, weights, weightsArray, i, weightSum)));
		
		GriddedGeoDataSet meanMap = readMap(zip, LogicTreeCurveAverager.MEAN_PREFIX+"_"+mapEntryName, reg);
		
		System.out.println("Waiting on "+futures.size()+" normCDF futures");
		for (int i=0; i<normCDFs.length; i++) {
			try {
				CDF_Callable call = futures.get(i).get();
				normCDFs[i] = call.normCDF;
				stdDevs[i] = call.stdDev;
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
			
		System.out.println("Done loading");
		
		return meanMap;
	}
	
	private static class CDF_Callable implements Callable<CDF_Callable> {
		
		// inputs
		private List<GriddedGeoDataSet> maps;
		private List<Double> weights;
		private double[] weightsArray;
		private int gridIndex;
		private double totWeight;
		
		// outputs
		private LightFixedXFunc normCDF;
		private double stdDev;

		public CDF_Callable(List<GriddedGeoDataSet> maps, List<Double> weights, double[] weightsArray, int gridIndex, double totWeight) {
			this.maps = maps;
			this.weights = weights;
			this.weightsArray = weightsArray;
			this.gridIndex = gridIndex;
			this.totWeight = totWeight;
		}

		@Override
		public CDF_Callable call() throws Exception {
			normCDF = LogicTreeHazardCompare.calcNormCDF(maps, weights, gridIndex, totWeight);
			Variance var = new Variance();
			double[] cellVals = new double[maps.size()];
			for (int j=0; j<cellVals.length; j++)
				cellVals[j] = maps.get(j).get(gridIndex);
			stdDev = Math.sqrt(var.evaluate(cellVals, weightsArray));
			
			return this;
		}
		
	}
	
	private static GriddedGeoDataSet readMap(ZipFile zip, String entryName, GriddedRegion reg) throws IOException {
		ZipEntry entry = zip.getEntry(entryName);
		Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(reg, false);
		String line = bRead.readLine();
		while (line != null) {
			line = line.trim();
			if (!line.startsWith("#")) {
				StringTokenizer tok = new StringTokenizer(line);
				double lon = Double.parseDouble(tok.nextToken());
				double lat = Double.parseDouble(tok.nextToken());
				double val = Double.parseDouble(tok.nextToken());
				Location loc = new Location(lat, lon);
				int index = reg.indexForLocation(loc);
				if (index > 0)
					xyz.set(index, val);
			}
			line = bRead.readLine();
		}
		return xyz;
	}
	
	private static GriddedGeoDataSet calcCOV(GriddedGeoDataSet meanMap, double[] stdDevs) {
		GriddedGeoDataSet cov = new GriddedGeoDataSet(meanMap.getRegion(), false);
		
		Preconditions.checkState(stdDevs.length == meanMap.size());
		
		for (int i=0; i<stdDevs.length; i++) {
			double val = stdDevs[i]/meanMap.get(i);
			if (val >= 1d)
				System.out.println("DEBUG >1: stdDev="+stdDevs[i]+", mean="+meanMap.get(i)+", cov="+val);
			cov.set(i, val);
		}
		
		return cov;
	}
	
	private static GriddedGeoDataSet logSpreadCOV(GriddedRegion reg, LightFixedXFunc[] normCDFs) {
		GriddedGeoDataSet cov = new GriddedGeoDataSet(reg, false);
		
		for (int i=0; i<normCDFs.length; i++)
			cov.set(i, Math.log10(normCDFs[i].getMaxX()/normCDFs[i].getMinX()));
		
		return cov;
	}
	
	private static GriddedGeoDataSet diff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion(), false);
		
		for (int i=0; i<ret.size(); i++)
			ret.set(i, map1.get(i) - map2.get(i));
		
		return ret;
	}
	
	private static GriddedGeoDataSet percentile(GriddedGeoDataSet meanMap, LightFixedXFunc[] ncdfs) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(meanMap.getRegion(), false);
		
		for (int i=0; i<ret.size(); i++) {
			double compVal = meanMap.get(i);
			
			double percentile;
			if (ncdfs[i].size() == 1) {
				if (compVal > 0d && (float)ncdfs[i].getX(0) == (float)compVal)
					percentile = 50d;
				else
					percentile = -1d;
			} else if (compVal < ncdfs[i].getMinX() || compVal > ncdfs[i].getMaxX()) {
				if ((float)compVal == (float)ncdfs[i].getMinX())
					percentile = 0d;
				else if ((float)compVal == (float)ncdfs[i].getMaxX())
					percentile = 100d;
				else if (compVal < ncdfs[i].getMinX())
					percentile = Double.NEGATIVE_INFINITY;
				else if (compVal > ncdfs[i].getMaxX())
					percentile = Double.POSITIVE_INFINITY;
				else
					percentile = Double.NaN;
			} else {
				percentile = 100d * ncdfs[i].getInterpolatedY(compVal);
			}
			ret.set(i, percentile);
		}
		
		return ret;
	}
	
	private static void plotMap(File outputDir, String prefix, GeographicMapMaker mapMaker, Region[] maskRegions,
			GriddedGeoDataSet map, CPT cpt, String label) throws IOException {
		System.out.println("Plotting "+label);
		map = MethodsAndIngredientsHazChangeFigures.mask(maskRegions, map);
		
		mapMaker.plotXYZData(map, cpt, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}

}
