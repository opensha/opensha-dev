package scratch.kevin.nshm23.wrapper;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.hazard.HazardCurveMetadata;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysHazardCalcSettings;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysHazardCalcSettings.CurveXValManager;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.calc.ReturnPeriod;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.imr.AttenRelSupplier;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.nshmp.shaded.model.NshmErf;
import mpi.MPI;

public class MPJ_WrapperHazardCalc extends MPJTaskCalculator {
	
	private File outputDir;
	private File outputFile;
	
	private static final double GRID_SPACING_DEFAULT = 0.1d;
	private double gridSpacing = GRID_SPACING_DEFAULT;
	
	private Map<TectonicRegionType, AttenRelSupplier> gmmSuppliers;
	private SourceFilterManager sourceFilters;
	private SourceFilterManager siteSkipSourceFilters;
	private boolean pointSourceOptimizations;
	private CurveXValManager xValManager;
	
//	private static final double[] PERIODS_DEFAULT = { 0d, 0.2d, 1d };
	public static final double[] PERIODS_DEFAULT = { 0d, 1d };
	private double[] periods = PERIODS_DEFAULT;
	
	private ReturnPeriod[] rps;
	
	private static final IncludeBackgroundOption GRID_SEIS_DEFAULT = IncludeBackgroundOption.INCLUDE;
	private IncludeBackgroundOption gridSeisOp = GRID_SEIS_DEFAULT;
	
	private GriddedSeismicitySettings gridSettings = GriddedSeismicitySettings.DEFAULT;
	
	private boolean noSubduction = false;
	private boolean noActive = false;
	private boolean noStable = false;
	
	private GridSourceProvider externalGridProv;
	
	private GriddedRegion gridRegion;

	private File nodesResultsDir;
	private File myResultsDir;
	
	private LocationList nodeLocs;
	private List<DiscretizedFunc[]> nodeCurves;
	
	private ExecutorService exec;
	
	private NshmErf erf;

	public MPJ_WrapperHazardCalc(CommandLine cmd) throws ZipException, IOException {
		super(cmd);
		
		this.shuffle = true;
		
		Path erfPath = Path.of(cmd.getOptionValue("input-dir"));
		
		outputDir = new File(cmd.getOptionValue("output-dir"));
		
		if (cmd.hasOption("faults-only"))
			gridSeisOp = IncludeBackgroundOption.EXCLUDE;
		if (cmd.hasOption("gridded-seis"))
			gridSeisOp = IncludeBackgroundOption.valueOf(cmd.getOptionValue("gridded-seis"));
		noSubduction = cmd.hasOption("no-subduction");
		noActive = cmd.hasOption("no-active");
		noStable = cmd.hasOption("no-stable");
		
		if (cmd.hasOption("grid-spacing"))
			gridSpacing = Double.parseDouble(cmd.getOptionValue("grid-spacing"));
		
		gmmSuppliers = FaultSysHazardCalcSettings.getGMMs(cmd);
		sourceFilters = FaultSysHazardCalcSettings.getSourceFilters(cmd);
		siteSkipSourceFilters = FaultSysHazardCalcSettings.getSiteSkipSourceFilters(sourceFilters, cmd);
		pointSourceOptimizations = FaultSysHazardCalcSettings.arePointSourceOptimizationsEnabled(cmd);
		xValManager = FaultSysHazardCalcSettings.getXValManager(cmd);
		gridSettings = FaultSysHazardCalcSettings.getGridSeisSettings(cmd);
		
		if (cmd.hasOption("periods")) {
			List<Double> periodsList = new ArrayList<>();
			String periodsStr = cmd.getOptionValue("periods");
			if (periodsStr.contains(",")) {
				String[] split = periodsStr.split(",");
				for (String str : split)
					periodsList.add(Double.parseDouble(str));
			} else {
				periodsList.add(Double.parseDouble(periodsStr));
			}
			periods = Doubles.toArray(periodsList);
		}
		
		File regFile = new File(cmd.getOptionValue("region"));
		Preconditions.checkState(regFile.exists(), "Supplied region file doesn't exist: %s", regFile.getAbsolutePath());
		Region region;
		if (regFile.getName().toLowerCase().endsWith(".zip")) {
			// it's a zip file, assume it's a prior hazard calc
			ZipFile zip = new ZipFile(regFile);
			ZipEntry regEntry = zip.getEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME);
			if (rank == 0) debug("Reading gridded region from zip file: "+regEntry.getName());
			BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(regEntry)));
			region = GriddedRegion.fromFeature(Feature.read(bRead));
			zip.close();
		} else {
			Feature feature = Feature.read(regFile);
			region = Region.fromFeature(feature);
		}
		if (region instanceof GriddedRegion) {
			gridRegion = (GriddedRegion)region;
			Preconditions.checkState(
					!cmd.hasOption("grid-spacing") || (float)gridSpacing == (float)gridRegion.getSpacing(),
					"Supplied a gridded region via the command line, cannont also specify grid spacing.");
			gridSpacing = gridRegion.getSpacing();
		} else {
			gridRegion = new GriddedRegion(region, gridSpacing, GriddedRegion.ANCHOR_0_0);
		}
		
		if (rank == 0) {
			MPJ_LogicTreeHazardCalc.waitOnDir(outputDir, 5, 1000);
			
			if (cmd.hasOption("output-file"))
				outputFile = new File(cmd.getOptionValue("output-file"));
			else
				outputFile = new File(outputDir.getParentFile(), "results_hazard.zip");
		}
		
		nodesResultsDir = new File(outputDir, "node_results");
		if (rank == 0) {
			if (nodesResultsDir.exists()) {
				// delete anything preexisting
				for (File file : nodesResultsDir.listFiles())
					Preconditions.checkState(FileUtils.deleteRecursive(file));
			} else {
				Preconditions.checkState(nodesResultsDir.mkdir());
			}
		}
		myResultsDir = new File(nodesResultsDir, "rank_"+rank);
		
		nodeLocs = new LocationList();
		nodeCurves = new ArrayList<>();
		
		if (cmd.hasOption("external-grid-prov")) {
			File gpFile = new File(cmd.getOptionValue("external-grid-prov"));
			Preconditions.checkState(gpFile.exists());
			ZipFile zip = new ZipFile(gpFile);
			
			if (FaultSystemSolution.isSolution(zip)) {
				externalGridProv = FaultSystemSolution.load(zip).requireModule(GridSourceProvider.class);
			} else {
				ModuleArchive<OpenSHA_Module> avgArchive = new ModuleArchive<>(zip);
				externalGridProv = avgArchive.requireModule(GridSourceProvider.class);
			}
			Preconditions.checkState(gridSeisOp != IncludeBackgroundOption.ONLY,
					"Background seismicity was set to ONLY, but is being overridden?");
			gridSeisOp = IncludeBackgroundOption.EXCLUDE;
			
			zip.close();
		}
		
		Set<TectonicRegionType> trts = EnumSet.noneOf(TectonicRegionType.class);
		if (!noActive) trts.add(TectonicRegionType.ACTIVE_SHALLOW);
		if (!noStable) trts.add(TectonicRegionType.STABLE_SHALLOW);
		if (!noSubduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}
		Preconditions.checkState(!trts.isEmpty(), "Must supply at least one TRT");
		if (gmmSuppliers.size() > 1)
			for (TectonicRegionType trt : trts)
				Preconditions.checkState(gmmSuppliers.containsKey(trt), "No GMM supplied for %s", trt);
		if (rank == 0) {
			debug("GMMs:");
			for (TectonicRegionType trt : gmmSuppliers.keySet())
				debug("\t"+trt+": "+gmmSuppliers.get(trt).getName());
			debug("Source filters: "+sourceFilters);
			debug("Site skip filters: "+siteSkipSourceFilters);
		}

		erf = new NshmErf(erfPath, trts, gridSeisOp);
		System.out.println("NSHM ERF size: " + erf.getNumSources());
		erf.getTimeSpan().setDuration(1.0);
		rps = ReturnPeriod.defaultsForCurveDuration(erf.getTimeSpan());
		erf.updateForecast();
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		
	}

	@Override
	protected int getNumTasks() {
		return gridRegion.getNodeCount();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		List<Future<CalcCallable>> futures = new ArrayList<>();
		
		debug("Calculating "+batch.length+" curves");
		for (int index : batch)
			futures.add(exec.submit(new CalcCallable(index)));
		
		for (Future<CalcCallable> future : futures) {
			CalcCallable call = future.get();
			
			nodeLocs.add(gridRegion.getLocation(call.index));
			nodeCurves.add(call.curves);
		}
		
		debug("DONE "+batch.length+" curves");
	}
	
	private ArrayDeque<Map<TectonicRegionType, ScalarIMR>> gmmsDeque;
	
	private synchronized Map<TectonicRegionType, ScalarIMR> checkOutGMMs() {
		if (gmmsDeque == null)
			gmmsDeque = new ArrayDeque<>();
		
		Map<TectonicRegionType, ScalarIMR> gmms;
		if (gmmsDeque.isEmpty()) {
			gmms = FaultSysHazardCalcSettings.getGmmInstances(gmmSuppliers);
		} else {
			gmms = gmmsDeque.pop();
		}
		
		return gmms;
	}
	
	private synchronized void checkInGMMs(Map<TectonicRegionType, ScalarIMR> gmms) {
		gmmsDeque.push(gmms);
	}
	
	private ArrayDeque<HazardCurveCalculator> calcDeque;
	
	private synchronized HazardCurveCalculator checkOutCurveCalc() {
		if (calcDeque == null)
			calcDeque = new ArrayDeque<>();
		
		HazardCurveCalculator calc;
		if (calcDeque.isEmpty()) {
			calc = new HazardCurveCalculator(sourceFilters);
			calc.setPointSourceOptimizationsEnabled(pointSourceOptimizations);
		} else {
			calc = calcDeque.pop();
		}
		
		return calc;
	}
	
	private synchronized void checkInCurveCalc(HazardCurveCalculator calc) {
		calcDeque.push(calc);
	}
	
	private class CalcCallable implements Callable<CalcCallable> {
		
		private int index;
		private DiscretizedFunc[] curves;
		
		public CalcCallable(int index) {
			super();
			this.index = index;
		}

		@Override
		public CalcCallable call() throws Exception {
			Location siteLoc = gridRegion.getLocation(index);
			Site site = new Site(siteLoc);
			if (canSkipSite(site)) {
				curves = new DiscretizedFunc[periods.length];
				for (int p=0; p<periods.length; p++) {
					curves[p] = xValManager.initLinearCurve(periods[p]);
				}
				return this;
			}
			
			Map<TectonicRegionType, ScalarIMR> gmms = checkOutGMMs();
			HazardCurveCalculator calc = checkOutCurveCalc();
			
			site.addParameterList(FaultSysHazardCalcSettings.getDefaultSiteParams(gmms));
			
			ExternalGridProvERF gridProvERF = null;
			if (externalGridProv != null)
				gridProvERF = new ExternalGridProvERF();
			
			curves = new DiscretizedFunc[periods.length];
			
			for (int p=0; p<periods.length; p++) {
				DiscretizedFunc logCurve = xValManager.initLogCurve(periods[p]);
				FaultSysHazardCalcSettings.setIMforPeriod(gmms, periods[p]);
				
				calc.getHazardCurve(logCurve, site, gmms, erf);
				
				if (externalGridProv != null) {
					// calc separately
					DiscretizedFunc gridCurve = xValManager.initLogCurve(periods[p]);
					calc.getHazardCurve(gridCurve, site, gmms, gridProvERF);
					
					// add it in
					for (int i=0; i<gridCurve.size(); i++)
						logCurve.set(i, 1d - (1d-gridCurve.getY(i)) * (1d-logCurve.getY(i)));
				}
				
				curves[p] = xValManager.remapToLinear(logCurve, periods[p]);
			}
			
			checkInGMMs(gmms);
			checkInCurveCalc(calc);
			
			return this;
		}
	}
	
	private boolean canSkipSite(Site site) {
		return SolHazardMapCalc.shouldSkipSite(site, siteSkipSourceFilters, erf, erf.getNumSources(), externalGridProv);
	}
	
	private class ExternalGridProvERF extends AbstractERF {

		@Override
		public int getNumSources() {
			return externalGridProv.getNumSources();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return externalGridProv.getSource(idx, 1d, null, gridSettings);
		}

		@Override
		public void updateForecast() {}

		@Override
		public String getName() {
			return "External Grid Source Provider";
		}
		
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		// write out all curves
		if (!nodeLocs.isEmpty()) {
			Preconditions.checkState(myResultsDir.exists() || myResultsDir.mkdir());
			for (int p=0; p<periods.length; p++) {
				DiscretizedFunc[] perCurves = new DiscretizedFunc[nodeLocs.size()];
				
				for (int s=0; s<nodeLocs.size(); s++)
					perCurves[s] = nodeCurves.get(s)[p];
				
				File curvesFile =
						new File(myResultsDir, SolHazardMapCalc.getCSV_FileName("curves", periods[p]) + ".gz");
				debug("writing node curves to "+curvesFile.getAbsolutePath());
				SolHazardMapCalc.writeCurvesCSV(curvesFile, perCurves, nodeLocs);
			}
		}
		
		if (!SINGLE_NODE_NO_MPJ)
			// wait for all to write out
			MPI.COMM_WORLD.Barrier();
		
		if (rank == 0) {
			// merge them
			
			List<String> zipFileNames = new ArrayList<>();
			HazardCurveMetadata curveMetadata = new HazardCurveMetadata(erf.getTimeSpan());
			File metadataFile = new File(outputDir, HazardCurveMetadata.FILE_NAME);
			curveMetadata.write(metadataFile);
			zipFileNames.add(metadataFile.getName());
			
			for (int p=0; p<periods.length; p++) {
				DiscretizedFunc[] perCurves = new DiscretizedFunc[gridRegion.getNodeCount()];
				
				String curvesFileName = SolHazardMapCalc.getCSV_FileName("curves", periods[p]) + ".gz";
				
				int curvesProcessed = 0;
				
				for (int oRank=0; oRank<size; oRank++) {
					DiscretizedFunc[] oCurves;
					LocationList oLocs;
					if (oRank == 0) {
						oCurves = new DiscretizedFunc[nodeLocs.size()];
						for (int s=0; s<nodeLocs.size(); s++)
							oCurves[s] = nodeCurves.get(s)[p];
						oLocs = nodeLocs;
					} else {
						File oResultsDir = new File(nodesResultsDir, "rank_"+oRank);
						if (oResultsDir.exists()) {
							File curvesFile = new File(oResultsDir, curvesFileName);
							Preconditions.checkState(curvesFile.exists());
							CSVFile<String> curvesCSV = CSVFile.readFile(curvesFile, true);
							oLocs = new LocationList();
							for (int row=1; row<curvesCSV.getNumRows(); row++) {
								double lat = curvesCSV.getDouble(row, 1);
								double lon = curvesCSV.getDouble(row, 2);
								oLocs.add(new Location(lat, lon));
							}
							oCurves = SolHazardMapCalc.loadCurvesCSV(curvesCSV, null);
						} else {
							continue;
						}
					}
					debug("Loaded "+oCurves.length+" curves from "+oRank);
					Preconditions.checkState(oCurves.length == oLocs.size());
					for (int i=0; i<oCurves.length; i++) {
						int nodeIndex = gridRegion.indexForLocation(oLocs.get(i));
						Preconditions.checkState(nodeIndex >= 0);
						Preconditions.checkState(perCurves[nodeIndex] == null, "Duplicate curve at index %s", nodeIndex);
						curvesProcessed++;
						perCurves[nodeIndex] = oCurves[i];
					}
				}
				File curvesFile = new File(outputDir, curvesFileName);
				Preconditions.checkState(curvesProcessed == gridRegion.getNodeCount(),
						"Processed %s curves, but gridded region has %s", curvesProcessed, gridRegion.getNodeCount());
				debug("Writing "+curvesProcessed+" curves to: "+curvesFile.getAbsolutePath());
				SolHazardMapCalc.writeCurvesCSV(curvesFile, perCurves, gridRegion.getNodeList());
				zipFileNames.add(curvesFileName);
				
				for (ReturnPeriod rp : rps) {
					GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridRegion, false);

					double curveLevel = rp.getProbability(1d);

					for (int i = 0; i < perCurves.length; i++) {
						DiscretizedFunc curve = perCurves[i];
						Preconditions.checkNotNull(curve, "Curve not calculated at index %s", i);
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

					String mapPrefix = MPJ_LogicTreeHazardCalc.mapPrefix(periods[p], rp);
					File mapFile = new File(outputDir, mapPrefix + ".txt");
					zipFileNames.add(mapFile.getName());

					GriddedGeoDataSet.writeXYZFile(xyz, mapFile);
				}
			}
			
			File gridFile = new File(outputDir, "gridded_region.geojson");
			Feature.write(gridRegion.toFeature(), gridFile);
			zipFileNames.add(gridFile.getName());
			
			debug("Writing results zip file: "+outputFile.getAbsolutePath());

			FileUtils.createZipFile(outputFile.getAbsolutePath(), outputDir.getAbsolutePath(), zipFileNames);
			
			debug("Done consolodating!");
		}
		
		exec.shutdown();
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		FaultSysHazardCalcSettings.addCommonOptions(ops, true);
		
		ops.addRequiredOption("id", "input-dir", true, "Path to inputs directory");
		ops.addRequiredOption("od", "output-dir", true, "Path to output directory");
		ops.addRequiredOption("r", "region", true, "Path to GeoJSON file containing a region for which we should "
				+ "compute hazard. Can be a gridded region or an outline. If a zip file is supplied, then it is assumed "
				+ "that the file is a prior hazard calculation zip file and the region will be reused from that prior "
				+ "calculation.");
		ops.addOption("of", "output-file", true, "Path to output zip file. Default will be based on the output directory");
		ops.addOption("sp", "grid-spacing", true, "Grid spacing in decimal degrees. Default: "+(float)GRID_SPACING_DEFAULT);
		ops.addOption("fo", "faults-only", false, "Flag to disable model gridded seismicity sources.");
		ops.addOption("gs", "gridded-seis", true, "Gridded seismicity option. One of "
				+FaultSysTools.enumOptions(IncludeBackgroundOption.class)+". Default: "+GRID_SEIS_DEFAULT.name());
		ops.addOption(null, "no-subduction", false, "Flag to disable subduction sources.");
		ops.addOption(null, "no-active", false, "Flag to disable active sources.");
		ops.addOption(null, "no-stable", false, "Flag to disable stable sources.");
		ops.addOption("egp", "external-grid-prov", true, "Path to external grid source provider to use for hazard "
				+ "calculations. Can be either a fault system solution, or a zip file containing just a grid source "
				+ "provider. Implies --faults-only.");
		
		return ops;
	}

	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_WrapperHazardCalc.class);
			
			MPJ_WrapperHazardCalc driver = new MPJ_WrapperHazardCalc(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
