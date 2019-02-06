package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

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
				this.sites = config.getValues(Site.class, Quantity.SITE);
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
			if (site.equals(this.site))
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
			return new RSQSimProbEqkRup(rotated.getMagnitude(), Double.NaN, Double.NaN, RSQSimUtils.getHypocenter(rotated),
					eventID, rotated.getTime(), null, rotated.getAllElements());
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
	
	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog;
		
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
//			catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
			catalog = Catalogs.BRUCE_2740.instance(baseDir);
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
			sites.add(new Site(new Location(34.064986, -117.29201), "SBSM"));
			sites.add(new Site(new Location(34.041824, -118.0653), "WNGC"));
			sites.add(new Site(new Location(33.93088, -118.17881), "STNI"));
			sites.add(new Site(new Location(34.557, -118.125), "LAPD"));
			sites.add(new Site(new Location(34.55314, -118.72826), "s119"));
			sites.add(new Site(new Location(34.37809, -118.34757), "s279"));
			sites.add(new Site(new Location(34.15755, -117.87389), "s480"));
			sites.add(new Site(new Location(34.00909, -118.48939), "SMCA"));
			
			Scenario[] scenarios = Scenario.values();
			
			double[] distances = { 20d, 50d, 100d };
			
			int numSourceAz = 18;
			int numSiteToSourceAz = 18;
			int maxRuptures = 100;
			
			int numRotations = 0;
			
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
			}
			
			int numRotationsPerSite = numRotations/sites.size();
			System.out.println("Created "+numRotations+", "+numRotationsPerSite+" per site");
		} else {
			List<Site> sites = null;
			for (Scenario scenario : Scenario.values()) {
				File csvFile = getCSVFile(outputDir, scenario);
				if (csvFile.exists()) {
					System.out.println("Located CSV for "+scenario+": "+csvFile.getAbsolutePath());
					List<RSQSimEvent> ruptures = null;
					if (writePoints || writeSRFs) {
						System.out.println("Loading ruptures...");
						ruptures = getScenarioEvents(catalog, scenario, skipYears, -1);
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
			Preconditions.checkState(!configsMap.isEmpty(), "No configuration CSV files found in %s", outputDir.getAbsolutePath());
		}
		
		RSQSimRotatedRuptureFakeERF erf = new RSQSimRotatedRuptureFakeERF(catalog, configsMap);
		int numRups = 0;
		for (ProbEqkSource source : erf)
			numRups += source.getNumRuptures();
		System.out.println("ERF has "+erf.getNumSources()+" sources and "+numRups+" ruptures");
		
		File csSourcesDir = new File(outputDir, "cs_source_rup_files");
		
		if (writePoints) {
			System.out.println("Writing rupture points to: "+csSourcesDir.getAbsolutePath());
			Preconditions.checkState(csSourcesDir.exists() || csSourcesDir.mkdir());
			erf.writeRupturePointFiles(outputDir);
		}
		if (writeSRFs) {
			System.out.println("Writing SRFs to: "+csSourcesDir.getAbsolutePath());
			Preconditions.checkState(csSourcesDir.exists() || csSourcesDir.mkdir());
			erf.writeRuptureSRFs(outputDir, catalog, dt, interpMode, momentPDiffThreshold);
		}
	}

}
