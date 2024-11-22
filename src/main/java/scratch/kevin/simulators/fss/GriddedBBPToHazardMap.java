package scratch.kevin.simulators.fss;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.IMT;
import scratch.kevin.simCompare.SimulationHazardCurveCalc;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.LightweightBBP_CatalogSimZipLoader;

public class GriddedBBPToHazardMap {

	public static void main(String[] args) throws IOException {
		File bbpDir, outputDir;
		double gridSpacing;
		RSQSimCatalog catalog;
		if (args.length == 1 && args[0].equals("--hardcoded")) {
			System.out.println("HARDCODED");
			catalog = Catalogs.BRUCE_5895.instance();
			bbpDir = new File("/data/kevin/bbp/parallel/"
					+ "2024_11_08-rundir5895-all-m6.0-skipYears10000-minSubSects2-maxDist300-noHF-vmLA_BASIN_500-griddedSitesWUS-gridSpacing0.5");
			outputDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_11_12-rsqsim-wus-5895-bbp_0.5");
			gridSpacing = 0.5;
		} else {
			if (args.length != 4) {
				System.err.println("USAGE: <catalog-name> <bbp-dir> <output-dir> <grid-spacing>");
				System.exit(1);
			}
			catalog = Catalogs.valueOf(args[0]).instance();
			bbpDir = new File(args[1]);
			outputDir = new File(args[2]);
			gridSpacing = Double.parseDouble(args[3]);
		}
		Preconditions.checkState(outputDir.exists());
		
		File rotDFile = new File(bbpDir, "results_rotD.zip");
		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		GriddedRegion gridReg = new GriddedRegion(region, gridSpacing, GriddedRegion.ANCHOR_0_0);
		double[] periods = {2d, 3d, 5d};
		ReturnPeriods[] rps = SolHazardMapCalc.MAP_RPS;
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		
		Preconditions.checkState(gridReg.getNodeCount() == sites.size(),
				"Have %s sites but %s gridded region nodes, spacing=%s",
				sites.size(), gridReg.getNodeCount(), (float)gridSpacing);
		
		ZipFile rotDZip = new ZipFile(rotDFile);
		
		BiMap<BBP_Site, Site> gmpeSites = HashBiMap.create();
		for (BBP_Site site : sites)
			gmpeSites.put(site, site.buildGMPE_Site(vm));
		
		double catDurationYears = catalog.getDurationYears();
		String bbpDirName = bbpDir.getName();
		Preconditions.checkState(bbpDirName.contains("-skipYears"));
		String yearStr = bbpDirName.substring(bbpDirName.indexOf("-skipYears")+"-skipYears".length());
		if (yearStr.contains("-"))
			yearStr = yearStr.substring(0, yearStr.indexOf("-"));
		int skipYears = Integer.parseInt(yearStr);
		System.out.println("Detected BBP skipYears="+skipYears);
		catDurationYears -= skipYears;
		System.out.println("Duration: "+(float)catDurationYears);
		
		LightweightBBP_CatalogSimZipLoader bbpLoader = new LightweightBBP_CatalogSimZipLoader(
				rotDZip, sites, gmpeSites, catDurationYears) {

			@Override
			public double getMinimumCurvePlotRate(Site site) {
				return Double.NaN;
			}
		};
		
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		List<Future<?>> preloadFutures = new ArrayList<>(sites.size());
		for (BBP_Site site : sites) {
			preloadFutures.add(exec.submit(new Runnable() {
				
				@Override
				public void run() {
					System.out.println("Preloading for "+site.getName());
					try {
						bbpLoader.preloadAllRotD(gmpeSites.get(site));
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
				}
			}));
		}
		
		for (Future<?> future : preloadFutures) {
			try {
				future.get();
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		System.out.println("Calculating curves for "+sites.size()+" sites");
		
		List<Future<DiscretizedFunc[]>> curveFutures = new ArrayList<>(sites.size());
		
		SimulationHazardCurveCalc<Integer> calc = new SimulationHazardCurveCalc<>(bbpLoader);
		
		for (int i=0; i<sites.size(); i++) {
			Site site = gmpeSites.get(sites.get(i));
			curveFutures.add(exec.submit(new Callable<DiscretizedFunc[]>() {

				@Override
				public DiscretizedFunc[] call() throws Exception {
					DiscretizedFunc[] ret = new DiscretizedFunc[periods.length];
					for (int p=0; p<periods.length; p++) {
						IMT imt = IMT.forPeriod(periods[p]);
						ret[p] = calc.calc(site, imt, 1d);
					}
					return ret;
				}
			}));
		}
		
		List<DiscretizedFunc[]> curves = new ArrayList<>(curveFutures.size());
		for (Future<DiscretizedFunc[]> future : curveFutures) {
			try {
				curves.add(future.get());
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				System.exit(1);
			}
			System.out.print(".");
			if (curves.size() % 100 == 0 || curves.size() == curveFutures.size())
				System.out.println(" "+curves.size());
		}
		
		exec.shutdown();
		
		File resultsDir = new File(outputDir, "results");
		Preconditions.checkState(resultsDir.exists() || resultsDir.mkdir());
		File hazDir = new File(resultsDir, "hazard_"+(float)gridSpacing+"deg_grid_seis_EXCLUDE");
		Preconditions.checkState(hazDir.exists() || hazDir.mkdir());
		File hazFile = new File(outputDir, "results_hazard_EXCLUDE.zip");
		
		System.out.println("Writing results to:\n\tDir: "+hazDir.getAbsolutePath()+"\n\tFile: "+hazFile.getAbsolutePath());
		
		ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(hazFile)));
		for (int p=0; p<periods.length; p++) {
			DiscretizedFunc[] perCurves = new DiscretizedFunc[sites.size()];
			LocationList locs = new LocationList(sites.size());
			DiscretizedFunc refCurve = null;
			for (int i=0; i<sites.size(); i++) {
				locs.add(sites.get(i).getLoc());
				perCurves[i] = curves.get(i)[p];
				if (refCurve == null && perCurves[i] != null)
					refCurve = perCurves[i];
			}
			for (int i=0; i<sites.size(); i++) {
				if (perCurves[i] == null) {
					DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
					for (int j=0; j<refCurve.size(); j++)
						curve.set(refCurve.getX(j), 0d);
					perCurves[i] = curve;
				} else {
					Preconditions.checkState(perCurves[i].size() == refCurve.size());
				}
			}
			File curvesFile = new File(hazDir, SolHazardMapCalc.getCSV_FileName("curves", periods[p]));
			SolHazardMapCalc.writeCurvesCSV(curvesFile, perCurves, locs, true);
			
			for (ReturnPeriods rp : rps) {
				String mapFileName = MPJ_LogicTreeHazardCalc.mapPrefix(periods[p], rp)+".txt";
				
				double curveLevel = rp.oneYearProb;
				GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
				for (int i=0; i<xyz.size(); i++) {
					DiscretizedFunc curve = perCurves[i];
					double val;
					// curveLevel is a probability, return the IML at that probability
					if (curveLevel > curve.getMaxY())
						val = 0d;
					else if (curveLevel < curve.getMinY())
						// saturated
						val = curve.getMaxX();
					else
						val = curve.getFirstInterpolatedX_inLogXLogYDomain(curveLevel);
					xyz.set(i, val);
				}
				
				zout.putNextEntry(new ZipEntry(mapFileName));
				ArbDiscrGeoDataSet.writeXYZStream(xyz, zout);
				zout.flush();
				zout.closeEntry();
				ArbDiscrGeoDataSet.writeXYZFile(xyz, new File(hazDir, mapFileName));
			}
		}
		
		Feature feature = gridReg.toFeature();
		
		zout.putNextEntry(new ZipEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME));
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(zout));
		Feature.write(feature, writer);
		writer.flush();
		zout.closeEntry();
		
		zout.close();
		rotDZip.close();
	}

}
