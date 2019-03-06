package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.IntegerParameter;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.erf.RSQSimSectBundledERF.RSQSimProbEqkRup;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

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
	private Map<Scenario, RotatedRupVariabilityConfig> configMap;
	
	private List<Site> sites;
	private List<Scenario> scenarios;
	private List<Float> distances;
	private Map<Scenario, List<Integer>> eventIDs;
	
	private List<RSQSimRotatedRuptureSource> sources = new ArrayList<>();

	public RSQSimRotatedRuptureFakeERF(RSQSimCatalog catalog, Map<Scenario, RotatedRupVariabilityConfig> configMap) {
		this.catalog = catalog;
		this.configMap = configMap;
		
		// loop over enum for constant order
		scenarios = new ArrayList<>();
		eventIDs = new HashMap<>();
		for (Scenario scenario : Scenario.values()) {
			if (!configMap.containsKey(scenario))
				continue;
			scenarios.add(scenario);
			RotatedRupVariabilityConfig config = configMap.get(scenario);
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
		
		sources = new ArrayList<>();
		
		for (Site site : sites) {
			for (Scenario scenario : scenarios) {
				RotatedRupVariabilityConfig config = configMap.get(scenario);
				for (Float distance : distances) {
					for (Integer eventID : eventIDs.get(scenario)) {
						List<RotationSpec> rotations = config.getRotationsForQuantities(Quantity.SITE, site,
								Quantity.DISTANCE, distance, Quantity.EVENT_ID, eventID);
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
		
		for (RotatedRupVariabilityConfig config : configMap.values()) {
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
			rupSurfResParam.setValue(catalog.getAveArea());
		} catch (ConstraintException | ParameterException | IOException e) {
			e.printStackTrace();
		}
		adjustableParams.addParameter(rupSurfResParam);
		
		this.timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		this.timeSpan.setDuration(1d);
	}
	
	public Map<Scenario, RotatedRupVariabilityConfig> getConfigMap() {
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
			if (LocationUtils.areSimilar(this.site.getLocation(), site.getLocation()))
				return distance;
			return Double.POSITIVE_INFINITY;
		}

		@Override
		public int getNumRuptures() {
			return rotations.size();
		}
		
		public RSQSimEvent getEvent(int nRupture) {
			return configMap.get(scenario).getRotatedRupture(rotations.get(nRupture));
		}

		@Override
		public RSQSimProbEqkRup getRupture(int nRupture) {
			RSQSimEvent rotated = getEvent(nRupture);
			Location hypo = RSQSimUtils.getHypocenter(rotated);
			RSQSimProbEqkRup rupture = new RSQSimProbEqkRup(rotated.getMagnitude(), Double.NaN, 1d, hypo,
					eventID, rotated.getTime(), null, rotated.getAllElements());
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
			if (LocationUtils.areSimilar(siteLoc, site.getLocation()))
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
						momentPDiffThreshold, patchAreas, momentStats);
			}
		}
		
		System.out.println("Transition file event moment % differences:");
		System.out.println("\tMean: "+(float)momentStats.getMean());
		System.out.println("\tMin: "+(float)momentStats.getMin());
		System.out.println("\tMax: "+(float)momentStats.getMax());
	}
	
	private static List<RSQSimEvent> getScenarioEvents(RSQSimCatalog catalog, Scenario scenario, int skipYears, int maxRuptures)
			throws IOException {
		List<RSQSimEvent> ruptures = scenario.getMatches(catalog, skipYears);
		if (maxRuptures > 0 && ruptures.size() > maxRuptures)
			ruptures = BBP_PartBValidationConfig.getBestMatches(scenario.getMagnitude(), ruptures, maxRuptures);
		return ruptures;
	}
	
	private static File getCSVFile(File outputDir, Scenario scenario) {
		return new File(outputDir, "rotation_config_"+scenario.getPrefix()+".csv");
	}
	
	public static Map<Scenario, RotatedRupVariabilityConfig> loadRotationConfigs(RSQSimCatalog catalog, File dir, boolean loadRuptures)
			throws IOException {
		Map<Scenario, RotatedRupVariabilityConfig> configsMap = new HashMap<>();
		List<Site> sites = null;
		for (Scenario scenario : Scenario.values()) {
			File csvFile = getCSVFile(dir, scenario);
			if (csvFile.exists()) {
				System.out.println("Located CSV for "+scenario+": "+csvFile.getAbsolutePath());
				List<RSQSimEvent> ruptures = null;
				if (loadRuptures) {
					System.out.println("Loading ruptures...");
					ruptures = getScenarioEvents(catalog, scenario, 0, -1);
					System.out.println("Loaded "+ruptures.size()+" ruptures");
				}
				
				System.out.println("Loading CSV...");
				RotatedRupVariabilityConfig config = RotatedRupVariabilityConfig.loadCSV(catalog, csvFile, ruptures, sites);
				if (sites == null)
					sites = config.getValues(Site.class, Quantity.SITE);
				System.out.println("Loaded "+config.getRotations().size()+" rotations");
				configsMap.put(scenario, config);
			}
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
		}

		@Override
		protected int getNumTasks() {
			return erf.getNumSources();
		}

		@Override
		protected void calculateBatch(int[] batch) throws Exception {
			for (int sourceID : batch) {
				if (writePoints) {
					debug("Writing point files for "+sourceID);
					RSQSimSectBundledERF.writeRupturePointFiles(erf, csSourcesDir, sourceID);
					debug("DONE writing point files for "+sourceID);
				}
				if (writeSRFs) {
					debug("Writing SRFs  for "+sourceID);
					RSQSimRotatedRuptureSource source = erf.getSource(sourceID);
					for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
						RSQSimEvent event = source.getEvent(rupID);
						RSQSimSectBundledERF.writeRuptureSRF(csSourcesDir, erf.catalog, sourceID, rupID, event, dt, interpMode,
								momentPDiffThreshold, patchAreas, null);
					}
					debug("DONE writing SRFs  for "+sourceID);
				}
			}
		}

		@Override
		protected void doFinalAssembly() throws Exception {}
		
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
	
	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog;
		
		boolean mpj = false;
		for (String arg : args)
			if (arg.startsWith("-") && arg.contains("dispatch"))
				mpj = true;
		
		CommandLine cmd = null;
		if (mpj) {
			System.out.println("MPJ!");
			args = MPJ_PointsAndSRFsWriter.initMPJ(args);
			
			Options options = MPJ_PointsAndSRFsWriter.createOptions();
			
			cmd = MPJ_PointsAndSRFsWriter.parse(options, args, MPJ_PointsAndSRFsWriter.class);
			args = cmd.getArgs();
		}
		
		boolean writePoints, writeSRFs, buildConfigs;
		
		if (args.length >= 1 && args.length < 5) {
			File catalogDir = new File(args[0]);
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
		
		File outputDir = new File(catalog.getCatalogDir(), "cybershake_rotation_inputs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int skipYears = 5000;
		
//		double dt = 0.05;
		double dt = 0.1;
		SRFInterpolationMode interpMode = SRFInterpolationMode.ADJ_VEL;
		double momentPDiffThreshold = 2;
		
		Map<Scenario, RotatedRupVariabilityConfig> configsMap = new HashMap<>();
		if (buildConfigs) {
			List<Site> sites = new ArrayList<>();
			
			sites.add(new Site(new Location(34.0192, -118.286), "USC"));
			sites.add(new Site(new Location(34.148426, -118.17119), "PAS"));
//			sites.add(new Site(new Location(34.55314, -118.72826), "s119"));
//			sites.add(new Site(new Location(34.37809, -118.34757), "s279"));
//			sites.add(new Site(new Location(34.15755, -117.87389), "s480"));
			sites.add(new Site(new Location(33.93088, -118.17881), "STNI"));
			sites.add(new Site(new Location(34.064986, -117.29201), "SBSM"));
			sites.add(new Site(new Location(34.041824, -118.0653), "WNGC"));
//			sites.add(new Site(new Location(34.557, -118.125), "LAPD"));
			sites.add(new Site(new Location(34.00909, -118.48939), "SMCA"));
			
			Scenario[] scenarios = { Scenario.M6p6_REVERSE, Scenario.M6p6_VERT_SS_SURFACE, Scenario.M7p2_VERT_SS_SURFACE };
			
			double[] distances = { 20d, 50d, 100d };
			
			int numSourceAz = 18;
			int numSiteToSourceAz = 36;
			int maxRuptures = 100;
			
			int numRotations = 0;
			
			long pointCount = 0;
			
			for (Scenario scenario : scenarios) {
				System.out.println("Loading ruptures for "+scenario);
				List<RSQSimEvent> ruptures = getScenarioEvents(catalog, scenario, skipYears, maxRuptures);
				System.out.println("Loaded "+ruptures.size()+" ruptures");
				
				System.out.println("Building rotation config...");
				RotatedRupVariabilityConfig config = new RotatedRupVariabilityConfig(
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
		} else {
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
