package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.exceptions.ConstraintException;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.IntegerParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import mpi.MPI;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.erf.RSQSimSectBundledERF.RSQSimProbEqkRup;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.rotation.RSQSimRotatedRupVariabilityConfig;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

/**
 * Fake ERF for running CyberShake for RSQSim Rotated Rupture variability experiments
 * 
 * There is an individual source for each Site, Scenario, Distance, and RSQSim Event (sorted in that order).
 * 
 * @author kevin
 *
 */
public class RSQSimRotatedRuptureFakeERF extends AbstractERF {
	
	private RSQSimCatalog catalog;
	private Map<Scenario, RSQSimRotatedRupVariabilityConfig> configMap;
	
	private List<Site> sites;
	private List<Scenario> scenarios;
	private List<Float> distances;
	private Map<Scenario, List<Integer>> eventIDs;
	
	private List<RSQSimRotatedRuptureSource> sources = new ArrayList<>();
	
	private int totNumRuptures;

	public RSQSimRotatedRuptureFakeERF(RSQSimCatalog catalog, Map<Scenario, RSQSimRotatedRupVariabilityConfig> configMap) {
		this.catalog = catalog;
		this.configMap = configMap;
		
		// loop over enum for constant order
		scenarios = new ArrayList<>();
		eventIDs = new HashMap<>();
		for (Scenario scenario : Scenario.values()) {
			if (!configMap.containsKey(scenario))
				continue;
			scenarios.add(scenario);
			RSQSimRotatedRupVariabilityConfig config = configMap.get(scenario);
			if (this.sites == null) {
				// need to do it this way in order to get sites in rotation order
				sites = new ArrayList<>();
				Site prevSite = null;
				for (RotationSpec rotation : config.getRotations()) {
					if (prevSite == null || !prevSite.getName().equals(rotation.site.getName())) {
						sites.add(rotation.site);
						prevSite = rotation.site;
					}
				}
				System.out.println("Loaded "+sites.size()+" sites");
			} else {
				// verify that they're the same sites
				List<Site> mySites = config.getValues(Site.class, Quantity.SITE);
				Preconditions.checkState(sites.size() == mySites.size(), "Site list size mismatch in scenarios");
				for (Site site : sites)
					Preconditions.checkState(mySites.contains(site),
							"Config for scenario %s doesn't contain site %s", scenario, site);
			}
			
			if (this.distances == null) {
				this.distances = config.getValues(Float.class, Quantity.DISTANCE);
			} else {
				// verify that they're the same distances
				List<Float> myDists = config.getValues(Float.class, Quantity.DISTANCE);
				Preconditions.checkState(distances.size() == myDists.size(), "Distance list size mismatch in scenarios");
				for (Float distance : distances)
					Preconditions.checkState(myDists.contains(distance),
							"Config for scenario %s doesn't contain distance %s", scenario, distance);
			}
			
			eventIDs.put(scenario, config.getValues(Integer.class, Quantity.EVENT_ID));
		}
		
		totNumRuptures = 0;
		sources = new ArrayList<>();
		
		for (Site site : sites) {
			for (Scenario scenario : scenarios) {
				RSQSimRotatedRupVariabilityConfig config = configMap.get(scenario);
				for (Float distance : distances) {
					for (Integer eventID : eventIDs.get(scenario)) {
						List<RotationSpec> rotations = config.getRotationsForQuantities(Quantity.SITE, site,
								Quantity.DISTANCE, distance, Quantity.EVENT_ID, eventID);
						totNumRuptures += rotations.size();
						sources.add(new RSQSimRotatedRuptureSource(site, scenario, distance, eventID, rotations));
					}
				}
			}
		}
		
		initParams();
	}
	
	private void initParams() {
		int numScenarios = configMap.size();
		int numSourceAz = -1;
		int numSiteToSourceAz = -1;
		int numDistances = -1;
		
		for (RSQSimRotatedRupVariabilityConfig config : configMap.values()) {
			numSourceAz = config.getValues(Float.class, Quantity.SOURCE_AZIMUTH).size();
			numSiteToSourceAz = config.getValues(Float.class, Quantity.SITE_TO_SOURTH_AZIMUTH).size();
			numDistances = config.getValues(Float.class, Quantity.DISTANCE).size();
		}
		
		IntegerParameter scenariosParam = new IntegerParameter("Num Scenarios", numScenarios, numScenarios);
		scenariosParam.setValue(numScenarios);
		adjustableParams.addParameter(scenariosParam);
		
		IntegerParameter sourceAzParam = new IntegerParameter("Num Source Azimuths", numSourceAz, numSourceAz);
		sourceAzParam.setValue(numSourceAz);
		adjustableParams.addParameter(sourceAzParam);
		
		IntegerParameter siteToSourceAzParam = new IntegerParameter("Num Site-To-Source Azimuths", numSiteToSourceAz, numSiteToSourceAz);
		siteToSourceAzParam.setValue(numSiteToSourceAz);
		adjustableParams.addParameter(siteToSourceAzParam);
		
		IntegerParameter distsParam = new IntegerParameter("Num Distances", numDistances, numDistances);
		distsParam.setValue(numDistances);
		adjustableParams.addParameter(distsParam);
		
		DoubleParameter rupSurfResParam = new DoubleParameter(RSQSimSectBundledERF.RUP_SURF_RESOLUTION_PARAM_NAME, 0d, 100d);
		try {
			rupSurfResParam.setValue(Math.sqrt(catalog.getAveArea()));
		} catch (ConstraintException | ParameterException | IOException e) {
			e.printStackTrace();
		}
		adjustableParams.addParameter(rupSurfResParam);
		
		this.timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		this.timeSpan.setDuration(1d);
	}
	
	public Map<Scenario, RSQSimRotatedRupVariabilityConfig> getConfigMap() {
		return configMap;
	}
	
	public RSQSimCatalog getCatalog() {
		return catalog;
	}

	@Override
	public int getNumSources() {
		return sources.size();
	}

	@Override
	public RSQSimRotatedRuptureSource getSource(int index) {
		return sources.get(index);
	}
	
	private boolean loadRuptures = true;
	public void setLoadRuptures(boolean loadRuptures) {
		this.loadRuptures = loadRuptures;
	}
	
	public class RSQSimRotatedRuptureSource extends ProbEqkSource {
		
		private Site site;
		private Scenario scenario;
		private Float distance;
		private Integer eventID;
		private List<RotationSpec> rotations;

		public RSQSimRotatedRuptureSource(Site site, Scenario scenario, Float distance, Integer eventID, List<RotationSpec> rotations) {
			this.site = site;
			this.scenario = scenario;
			this.distance = distance;
			this.eventID = eventID;
			this.rotations = rotations;
			
			this.name = "Site: "+site.getName()+"; Scenario "+scenario.name()+"; Distance "+distance+"; EventID "+eventID;
		}

		@Override
		public LocationList getAllSourceLocs() {
			throw new UnsupportedOperationException("Not applicable");
		}

		@Override
		public RuptureSurface getSourceSurface() {
			throw new UnsupportedOperationException("Not applicable");
		}

		@Override
		public double getMinDistance(Site site) {
			if (LocationUtils.horzDistanceFast(this.site.getLocation(), site.getLocation()) < 0.1d)
				return distance;
			return Double.POSITIVE_INFINITY;
		}

		@Override
		public int getNumRuptures() {
			return rotations.size();
		}
		
		public RSQSimEvent getEvent(int nRupture) {
			checkLoadEvents();
			return configMap.get(scenario).getRotatedRupture(rotations.get(nRupture));
		}

		@Override
		public RSQSimProbEqkRup getRupture(int nRupture) {
			double prob = 1d/(double)totNumRuptures;
			RSQSimProbEqkRup rupture;
			Location hypo;
			if (loadRuptures) {
				RSQSimEvent rotated = getEvent(nRupture);
				hypo = RSQSimUtils.getHypocenter(rotated);
				rupture = new RSQSimProbEqkRup(rotated.getMagnitude(), Double.NaN, prob, hypo,
						eventID, rotated.getTime(), null, rotated.getAllElements());
			} else {
				hypo = null;
				rupture = new RSQSimProbEqkRup(Double.NaN, Double.NaN, prob, hypo,
						eventID, Double.NaN, null, null);
			}
			rupture.setRuptureSurface(new FakeModDistanceSurface(hypo, site, distance));
			return rupture;
		}
		
		public List<RotationSpec> getRotations() {
			return rotations;
		}
		
		public Site getSite() {
			return site;
		}
		
	}
	
	private synchronized void checkLoadEvents() {
		if (configMap.values().iterator().next().hasRuptures())
			return;
		// have to load them
		HashSet<Integer> eventIDs = new HashSet<>();
		for (RSQSimRotatedRupVariabilityConfig config : configMap.values())
			eventIDs.addAll(config.getValues(Integer.class, Quantity.EVENT_ID));
		
		System.out.println("Loading "+eventIDs.size()+" events...");
		List<RSQSimEvent> events;
		try {
			events = catalog.loader().byIDs(Ints.toArray(eventIDs));
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		System.out.println("Done loading "+events.size()+" events");
		Preconditions.checkState(events.size() == eventIDs.size(), 
				"Loaded %s events, expected %s", events.size(), eventIDs.size());
		for (RSQSimRotatedRupVariabilityConfig config : configMap.values())
			config.setRuptures(events);
	}
	
	/**
	 * Simple surface implementation hardcoded to return the given rRup distance for the given site,
	 * and positive infinity otherwise. Extends PointSurface for brevity.
	 * @author kevin
	 *
	 */
	private class FakeModDistanceSurface extends PointSurface {
		
		private Site site;
		private double distance;

		public FakeModDistanceSurface(Location hypocenter, Site site, double distance) {
			super(hypocenter);
			this.site = site;
			this.distance = distance;
		}

		@Override
		public double getDistanceRup(Location siteLoc) {
			if (LocationUtils.horzDistanceFast(siteLoc, site.getLocation()) < 0.1d)
				return distance;
			return Double.POSITIVE_INFINITY;
		}

		@Override
		public double getDistanceJB(Location siteLoc) {
			throw new UnsupportedOperationException("Not applicable");
		}

		@Override
		public double getDistanceSeis(Location siteLoc) {
			throw new UnsupportedOperationException("Not applicable");
		}

		@Override
		public double getDistanceX(Location siteLoc) {
			throw new UnsupportedOperationException("Not applicable");
		}
		
	}

	@Override
	public void updateForecast() {}

	@Override
	public String getName() {
		return "RSQSim Rotation ERF for "+catalog.getName()+" with "+scenarios.size()+" scenarios, "
				+distances.size()+" distances";
	}
	
	public void writeRupturePointFiles(File mainDir) throws IOException {
		RSQSimSectBundledERF.writeRupturePointFiles(this, mainDir);
	}
	
	public void writeRuptureSRFs(File mainDir, RSQSimCatalog catalog, double dt, SRFInterpolationMode interpMode,
			double momentPDiffThreshold) throws IOException {
		
		Map<Integer, Double> patchAreas = new HashMap<>();
		for (SimulatorElement elem : catalog.getElements())
			patchAreas.put(elem.getID(), elem.getArea());
		
		SummaryStatistics momentStats = new SummaryStatistics();
		for (int sourceID=0; sourceID<getNumSources(); sourceID++) {
			RSQSimRotatedRuptureSource source = getSource(sourceID);
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				RSQSimEvent event = source.getEvent(rupID);
				RSQSimSectBundledERF.writeRuptureSRF(mainDir, catalog, sourceID, rupID, event, dt, interpMode,
						momentPDiffThreshold, patchAreas, momentStats, false);
			}
		}
		
		System.out.println("Transition file event moment % differences:");
		System.out.println("\tMean: "+(float)momentStats.getMean());
		System.out.println("\tMin: "+(float)momentStats.getMin());
		System.out.println("\tMax: "+(float)momentStats.getMax());
	}
	
	private static List<RSQSimEvent> getScenarioEvents(RSQSimCatalog catalog, Scenario scenario, FilterMethod filter, int skipYears, int maxRuptures)
			throws IOException {
		List<RSQSimEvent> ruptures = scenario.getMatches(catalog, skipYears);
		if (maxRuptures > 0 && ruptures.size() > maxRuptures)
			ruptures = filter.filter(ruptures, maxRuptures, catalog, scenario.getMagnitude());
		return ruptures;
	}
	
	private static File getCSVFile(File outputDir, Scenario scenario) {
		return new File(outputDir, "rotation_config_"+scenario.getPrefix()+".csv");
	}
	
	public static Map<Scenario, RSQSimRotatedRupVariabilityConfig> loadRotationConfigs(RSQSimCatalog catalog, File dir, boolean loadRuptures)
			throws IOException {
		Map<Scenario, RSQSimRotatedRupVariabilityConfig> configsMap = new HashMap<>();
		List<Site> sites = null;
		HashSet<Integer> eventIDs = new HashSet<>();
		for (Scenario scenario : Scenario.values()) {
			File csvFile = getCSVFile(dir, scenario);
			if (csvFile.exists()) {
				System.out.println("Located CSV for "+scenario+": "+csvFile.getAbsolutePath());
				System.out.println("Loading CSV...");
				RSQSimRotatedRupVariabilityConfig config = RSQSimRotatedRupVariabilityConfig.loadCSV(catalog, csvFile, null, sites);
				if (sites == null)
					sites = config.getValues(Site.class, Quantity.SITE);
				System.out.println("Loaded "+config.getRotations().size()+" rotations");
				configsMap.put(scenario, config);
				eventIDs.addAll(config.getValues(Integer.class, Quantity.EVENT_ID));
			}
		}
		if (loadRuptures) {
			System.out.println("Loading "+eventIDs.size()+" events");
			List<RSQSimEvent> events = catalog.loader().byIDs(Ints.toArray(eventIDs));
			System.out.println("Loaded "+events.size()+" events");
			for (RSQSimRotatedRupVariabilityConfig config : configsMap.values())
				config.setRuptures(events);
		}
		return configsMap;
	}
	
	private static class MPJ_PointsAndSRFsWriter extends MPJTaskCalculator {
		
		private RSQSimRotatedRuptureFakeERF erf;
		private File csSourcesDir;
		private boolean writePoints;
		private boolean writeSRFs;
		private double dt;
		private SRFInterpolationMode interpMode;
		private double momentPDiffThreshold;
		private File scratchDir;
		
		private ExecutorService cleanExec;
		private List<Future<?>> cleanFutures;
		
		private Map<Integer, Double> patchAreas;

		public MPJ_PointsAndSRFsWriter(CommandLine cmd, RSQSimRotatedRuptureFakeERF erf, File csSourcesDir,
				boolean writePoints, boolean writeSRFs, double dt, SRFInterpolationMode interpMode, double momentPDiffThreshold)
						throws IOException {
			super(cmd);
			this.shuffle = false;
			
			this.erf = erf;
			this.csSourcesDir = csSourcesDir;
			this.writePoints = writePoints;
			this.writeSRFs = writeSRFs;
			this.dt = dt;
			this.interpMode = interpMode;
			this.momentPDiffThreshold = momentPDiffThreshold;
			
			if (rank == 0)
				Preconditions.checkState(csSourcesDir.exists() || csSourcesDir.mkdir());
			
			patchAreas = new HashMap<>();
			for (SimulatorElement elem : erf.catalog.getElements())
				patchAreas.put(elem.getID(), elem.getArea());
			
			scratchDir = Files.createTempDir();
			
			cleanExec = Executors.newFixedThreadPool(1);
			cleanFutures = new ArrayList<>();
		}

		@Override
		protected int getNumTasks() {
			return erf.getNumSources();
		}
		
		private void waitOnDir(File dir) {
			int retry = 0;
			while (!(dir.exists() || dir.mkdir())) {
				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				if (retry++ > 5)
					throw new IllegalStateException("Directory doesn't exist and couldn't be created after "
							+5+" retries: "+dir.getAbsolutePath());
			}
		}

		@Override
		protected void calculateBatch(int[] batch) throws Exception {
			for (int sourceID : batch) {
				File sourceDir = new File(scratchDir, sourceID+"");
				File zipFile = new File(csSourcesDir, sourceID+".zip");
				if (zipFile.exists()) {
					// try opening it
					debug("checking existing zip for "+sourceID);
					if (validateZipFile(zipFile)) {
						debug("zip is valid, skipping");
						continue;
					} else {
						debug("bad zip for "+sourceID+", will rebuild");
					}
				}
				if (writePoints || writeSRFs) {
					// build directories first
					waitOnDir(sourceDir);
					RSQSimRotatedRuptureSource source = erf.getSource(sourceID);
					for (int rupID=0; rupID<source.getNumRuptures(); rupID++)
						waitOnDir(new File(sourceDir, rupID+""));
				}
				if (writePoints) {
					debug("Writing point files for "+sourceID);
					RSQSimSectBundledERF.writeRupturePointFiles(erf, scratchDir, sourceID, false);
					debug("DONE writing point files for "+sourceID);
				}
				if (writeSRFs) {
					debug("Writing SRFs  for "+sourceID);
					RSQSimRotatedRuptureSource source = erf.getSource(sourceID);
					for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
						RSQSimEvent event = source.getEvent(rupID);
						RSQSimSectBundledERF.writeRuptureSRF(scratchDir, erf.catalog, sourceID, rupID, event, dt, interpMode,
								momentPDiffThreshold, patchAreas, null, false);
					}
					debug("DONE writing SRFs  for "+sourceID);
				}
				FileUtils.createZipFile(zipFile, sourceDir, false);
				cleanFutures.add(cleanExec.submit(new Runnable() {
					
					@Override
					public void run() {
						FileUtils.deleteRecursive(sourceDir);
					}
				}));
			}
		}
		
		private boolean validateZipFile(File zipFile) {
			try {
				new ZipFile(zipFile).close();
				// if we get here, we're good
				return true;
			} catch (Exception e) {
				return false;
			}
		}

		@Override
		protected void doFinalAssembly() throws Exception {
			for (Future<?> f : cleanFutures)
				f.get();
			cleanExec.shutdown();
		}
		
		protected static String[] initMPJ(String[] args) {
			return MPJTaskCalculator.initMPJ(args);
		}
		
		protected static Options createOptions() {
			return MPJTaskCalculator.createOptions();
		}
		
		protected static CommandLine parse(Options options, String args[], Class<?> clazz) {
			return MPJTaskCalculator.parse(options, args, clazz);
		}
		
		protected static void finalizeMPJ() {
			MPJTaskCalculator.finalizeMPJ();
		}
		
	}
	
	private static File copyCatalogDir(File catDir, File scratchDir) throws IOException {
		if (!scratchDir.exists()) {
			scratchDir.mkdir();
			MPJ_BBP_Utils.waitOnDir(scratchDir, 10, 2000);
		}
		File destDir = new File(scratchDir, catDir.getName());
		if (!destDir.exists()) {
			destDir.mkdir();
			MPJ_BBP_Utils.waitOnDir(destDir, 10, 2000);
		}
		List<File> filesToCopy = new ArrayList<>();
		for (File f : catDir.listFiles()) {
			String name = f.getName().toLowerCase();
			if (name.endsWith("list"))
				filesToCopy.add(f);
			else if (name.startsWith("trans.") && name.endsWith(".out"))
				filesToCopy.add(f);
			else if (name.startsWith("transv.") && name.endsWith(".out"))
				filesToCopy.add(f);
			else if (name.endsWith(".in"))
				filesToCopy.add(f);
			else if (name.endsWith(".flt"))
				filesToCopy.add(f);
		}
		if (filesToCopy.size() < 7) {
			String error = "Need at least 7 files: 4 list, trans, input, geom. Have "+filesToCopy.size()+":";
			for (File f : filesToCopy)
				error += "\n\t"+f.getAbsolutePath();
			throw new IllegalStateException(error);
		}
		
		System.out.println("copying "+filesToCopy.size()+" catalog files to "+destDir.getAbsolutePath());
		for (File f : filesToCopy) {
			File destFile = new File(destDir, f.getName());
			if (destFile.exists() && destFile.length() == f.length())
				// skip copy, already exists
				continue;
			Files.copy(f, destFile);
		}
		
		return destDir;
	}
	
	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog;
		
		boolean mpj = false;
		for (String arg : args)
			if (arg.startsWith("-") && arg.contains("dispatch"))
				mpj = true;
		
		CommandLine cmd = null;
		int rank = -1;
		if (mpj) {
			args = MPJ_PointsAndSRFsWriter.initMPJ(args);
			rank = MPI.COMM_WORLD.Rank();
			System.out.println("MPJ! "+rank);
			
			Options options = MPJ_PointsAndSRFsWriter.createOptions();
			
			cmd = MPJ_PointsAndSRFsWriter.parse(options, args, MPJ_PointsAndSRFsWriter.class);
			args = cmd.getArgs();
		}
		
		boolean writePoints, writeSRFs, buildConfigs;
		
		File origCatalogDir = null;
		
		if (args.length >= 1 && args.length < 5) {
			origCatalogDir = new File(args[0]);
			
			File catalogDir = origCatalogDir;
			
			if (mpj) {
				// copy catalog data over to shared scratch
				File sharedScratch = new File(USC_HPCC_ScriptWriter.SHARED_SCRATCH_DIR+"/kmilner/catalogs");
				if (rank == 0)
					copyCatalogDir(catalogDir, sharedScratch);
				MPI.COMM_WORLD.Barrier();
				catalogDir = new File(sharedScratch, catalogDir.getName());
			}
			
			catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
					null, null, null, FaultModels.FM3_1, DeformationModels.GEOLOGIC); // TODO
			
			if (args.length >= 2)
				writePoints = Boolean.parseBoolean(args[1]);
			else
				writePoints = true;
			if (args.length >= 3)
				writeSRFs = Boolean.parseBoolean(args[2]);
			else
				writeSRFs = true;
			if (args.length == 4)
				buildConfigs = Boolean.parseBoolean(args[3]);
			else
				buildConfigs = true;
		} else {
			System.out.println("Assuming hardcoded. Otherwise usage is:");
			System.out.println("<catalogDir> [writePoints writeSRFs buildConfigs]");
			File baseDir = new File("/data/kevin/simulators/catalogs");
			if (!baseDir.exists()) {
				System.err.println("hardcoded dir doesn't exist: "+baseDir.getAbsolutePath());
				System.exit(2);
			}
			catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//			catalog = Catalogs.BRUCE_2740.instance(baseDir);
			writePoints = false;
			writeSRFs = false;
			buildConfigs = true;
		}
		
		File outputDir = new File(origCatalogDir == null ? catalog.getCatalogDir() : origCatalogDir,
				"cybershake_rotation_inputs");
		System.out.println("Output dir: "+outputDir.getAbsolutePath());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int skipYears = 5000;
		
//		double dt = 0.05;
		double dt = 0.1;
		FilterMethod filter = BBP_PartBValidationConfig.FILTER_METHOD_DEFAULT;
		SRFInterpolationMode interpMode = SRFInterpolationMode.ADJ_VEL;
		double momentPDiffThreshold = 2;
		
		Map<Scenario, RSQSimRotatedRupVariabilityConfig> configsMap = new HashMap<>();
		boolean mpjBuild = mpj && buildConfigs;
		if (mpjBuild && rank > 0)
			buildConfigs = false;
		if (mpj) {
			System.out.println("rank "+rank+": build="+buildConfigs+", mpjBuild="+mpjBuild);
		}
		if (buildConfigs) {
			List<Site> sites = new ArrayList<>();
			
//			sites.add(new Site(new Location(34.0192, -118.286), "USC"));
//			sites.add(new Site(new Location(34.148426, -118.17119), "PAS"));
////			sites.add(new Site(new Location(34.55314, -118.72826), "s119"));
////			sites.add(new Site(new Location(34.37809, -118.34757), "s279"));
////			sites.add(new Site(new Location(34.15755, -117.87389), "s480"));
//			sites.add(new Site(new Location(33.93088, -118.17881), "STNI"));
//			sites.add(new Site(new Location(34.064986, -117.29201), "SBSM"));
//			sites.add(new Site(new Location(34.041824, -118.0653), "WNGC"));
////			sites.add(new Site(new Location(34.557, -118.125), "LAPD"));
//			sites.add(new Site(new Location(34.00909, -118.48939), "SMCA"));
//			sites.add(new Site(new Location(34.6145, -118.7235), "OSI"));
//			sites.add(new Site(new Location(34.44199, -118.58215), "PDE"));
//			sites.add(new Site(new Location(34.1717, -118.64971), "WSS"));
//			sites.add(new Site(new Location(33.86889, -118.33143), "LAF"));
//			sites.add(new Site(new Location(34.24505, -119.18086), "s022"));
			

			sites.add(new Site(new Location(34.0192, -118.286), "USC"));
			sites.add(new Site(new Location(34.064986, -117.29201), "SBSM"));
			sites.add(new Site(new Location(34.041824, -118.0653), "WNGC"));
			sites.add(new Site(new Location(33.93088, -118.17881), "STNI"));
			sites.add(new Site(new Location(34.00909, -118.48939), "SMCA"));
			sites.add(new Site(new Location(34.6145, -118.7235), "OSI"));
			sites.add(new Site(new Location(34.44199, -118.58215), "PDE"));
			sites.add(new Site(new Location(34.1717, -118.64971), "WSS"));
			sites.add(new Site(new Location(33.86889, -118.33143), "LAF"));
			sites.add(new Site(new Location(34.24505, -119.18086), "s022"));
			
			Scenario[] scenarios = { Scenario.M6p6_REVERSE, Scenario.M6p6_VERT_SS_SURFACE, Scenario.M7p2_VERT_SS_SURFACE };
			
			double[] distances = { 20d, 50d, 100d };
			
//			int numSourceAz = 18;
//			int numSiteToSourceAz = 36;
//			int maxRuptures = 100;
			
			int numSourceAz = 18;
			int numSiteToSourceAz = 18;
			int maxRuptures = 50;
			
			int numRotations = 0;
			
			long pointCount = 0;
			
			for (Scenario scenario : scenarios) {
				System.out.println("Loading ruptures for "+scenario);
				List<RSQSimEvent> ruptures = getScenarioEvents(catalog, scenario, filter, skipYears, maxRuptures);
				System.out.println("Loaded "+ruptures.size()+" ruptures");
				
				System.out.println("Building rotation config...");
				RSQSimRotatedRupVariabilityConfig config = new RSQSimRotatedRupVariabilityConfig(
						catalog, sites, ruptures, distances, numSourceAz, numSiteToSourceAz);
				System.out.println("Built "+config.getRotations().size()+" rotations");
				numRotations += config.getRotations().size();
				
				File csvFile = getCSVFile(outputDir, scenario);
				System.out.println("Writing config to: "+csvFile.getAbsolutePath());
				config.writeCSV(csvFile);
				configsMap.put(scenario, config);
				
				Map<Integer, RSQSimEvent> eventIDMap = new HashMap<>();
				for (RSQSimEvent event : ruptures)
					eventIDMap.put(event.getID(), event);
				for (RotationSpec rotation : config.getRotations())
					pointCount += eventIDMap.get(rotation.eventID).getNumElements();
			}
			
			int numRotationsPerSite = numRotations/sites.size();
			System.out.println("Created "+numRotations+", "+numRotationsPerSite+" per site");
			
			long pointsPerSite = pointCount/sites.size();
			System.out.println("Will have "+pointsPerSite+" points per site");
			
			if (mpjBuild)
				MPI.COMM_WORLD.Barrier();
		} else {
			if (mpjBuild)
				MPI.COMM_WORLD.Barrier();
			configsMap = loadRotationConfigs(catalog, outputDir, writePoints || writeSRFs);
			Preconditions.checkState(!configsMap.isEmpty(), "No configuration CSV files found in %s", outputDir.getAbsolutePath());
		}
		
		RSQSimRotatedRuptureFakeERF erf = new RSQSimRotatedRuptureFakeERF(catalog, configsMap);
		int numRups = 0;
		for (ProbEqkSource source : erf)
			numRups += source.getNumRuptures();
		System.out.println("ERF has "+erf.getNumSources()+" sources and "+numRups+" ruptures");
		
		File csSourcesDir = new File(outputDir, "cs_source_rup_files");
		
		if (mpj) {
			try {
				MPJ_PointsAndSRFsWriter mpjWrite = new MPJ_PointsAndSRFsWriter(cmd, erf, csSourcesDir, writePoints, writeSRFs,
						dt, interpMode, momentPDiffThreshold);
				mpjWrite.run();
				
				MPJ_PointsAndSRFsWriter.finalizeMPJ();
				
				System.exit(0);
			} catch (Throwable t) {
				MPJTaskCalculator.abortAndExit(t);
			}
		} else {
			if (writePoints) {
				System.out.println("Writing rupture points to: "+csSourcesDir.getAbsolutePath());
				Preconditions.checkState(csSourcesDir.exists() || csSourcesDir.mkdir());
				erf.writeRupturePointFiles(csSourcesDir);
			}
			if (writeSRFs) {
				System.out.println("Writing SRFs to: "+csSourcesDir.getAbsolutePath());
				Preconditions.checkState(csSourcesDir.exists() || csSourcesDir.mkdir());
				erf.writeRuptureSRFs(csSourcesDir, catalog, dt, interpMode, momentPDiffThreshold);
			}
		}
	}

}
