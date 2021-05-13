package scratch.kevin.simulators.ruptures.rotation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.EventIDsRupIden;
import org.opensha.sha.simulators.iden.LinearRuptureIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import mpi.MPI;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.ruptures.AbstractMPJ_BBP_MultiRupSim;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.MPJ_BBP_PartBSim;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public class MPJ_BBP_RotatedRupVariabilityScenarioSim extends AbstractMPJ_BBP_MultiRupSim {

	private List<BBP_Site> sites;
	private Map<Site, BBP_Site> siteToBBPSites;
	
	private List<Config> configs;
	
	private MPJ_BBP_RotatedRupVariabilityScenarioSim(CommandLine cmd) throws IOException {
		super(cmd);
		this.shuffle = false; // do them in order for better caching
		
		int skipYears = -1;
		if (cmd.hasOption("skip-years"))
			skipYears = Integer.parseInt(cmd.getOptionValue("skip-years"));
		
		int maxRups = Integer.MAX_VALUE;
		if (cmd.hasOption("max-ruptures"))
			maxRups = Integer.parseInt(cmd.getOptionValue("max-ruptures"));
		
		File sitesFile = new File(cmd.getOptionValue("sites-file"));
		Preconditions.checkState(sitesFile.exists());
		sites = BBP_Site.readFile(sitesFile);
		siteToBBPSites = new HashMap<>();
		
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site bbpSite : sites) {
			Site site = bbpSite.buildGMPE_Site(null);
			siteToBBPSites.put(site, bbpSite);
			gmpeSites.add(site);
		}
		
		int numSourceAz = Integer.parseInt(cmd.getOptionValue("num-source-az"));
		int numSiteToSourceAz = Integer.parseInt(cmd.getOptionValue("num-site-source-az"));
		
		double[] distances = MPJ_BBP_PartBSim.getDistances(cmd);
		Scenario[] scenarios = MPJ_BBP_PartBSim.getScenarios(cmd);
		
		FilterMethod filter = MPJ_BBP_PartBSim.getFilterMethod(cmd);
		
		File[] csvFiles = new File[scenarios.length];
		RSQSimRotatedRupVariabilityConfig[] rotConfigs =
				new RSQSimRotatedRupVariabilityConfig[scenarios.length];
		List<List<RuptureIdentifier>> scenarioIdens = new ArrayList<>();
		List<RuptureIdentifier> loadIdens = new ArrayList<>();
		for (int s=0; s<scenarios.length; s++) {
			csvFiles[s] = new File(mainOutputDir, "rotation_config_"+scenarios[s].getPrefix()+".csv");
			List<RuptureIdentifier> idens = scenarios[s].getIdentifiers(catalog);
			scenarioIdens.add(idens);
			if (csvFiles[s].exists()) {
				debug("loading previous CSV: "+csvFiles[s].getName());
				rotConfigs[s] = RSQSimRotatedRupVariabilityConfig.loadCSV(catalog, csvFiles[s]);
				loadIdens.add(new EventIDsRupIden(
						Ints.toArray(rotConfigs[s].getValues(Integer.class, Quantity.EVENT_ID))));
			} else {
				for (RuptureIdentifier iden : idens)
					if (iden instanceof MagRangeRuptureIdentifier)
						loadIdens.add(iden);
			}
		}
		
		MPI.COMM_WORLD.Barrier();
		
		if (rank == 0)
			debug("Loading events for "+scenarios.length+" scenarios");
		List<RSQSimEvent> allEvents = catalog.loader().skipYears(skipYears).matches(new LogicalOrRupIden(loadIdens)).hasTransitions().load();
		debug("Loaded "+allEvents.size()+" potential matches");
		Map<Integer, RSQSimEvent> idEventMap = new HashMap<>();
		for (RSQSimEvent event : allEvents)
			idEventMap.put(event.getID(), event);
		
		configs = new ArrayList<>();
		
		for (int s=0; s<scenarios.length; s++) {
			Scenario scenario = scenarios[s];
			

			if (rank == 0)
				debug("Loading matches for scenario: "+scenario);
			
			List<RSQSimEvent> eventMatches;
			if (rotConfigs[s] == null) {
				eventMatches = new LogicalAndRupIden(scenarioIdens.get(s)).getMatches(allEvents);
				
				if (rank == 0)
					debug("Loaded "+eventMatches.size()+" matches for scenario: "+scenario);
				
				if (eventMatches.size() > maxRups) {
					eventMatches = filter.filter(eventMatches, maxRups, catalog, scenario.getMagnitude());
					debug("trimmed down to max of "+eventMatches.size()+" ruptures");
				}
				
				if (eventMatches.isEmpty()) {
					if (rank == 0)
						debug("skipping...");
					continue;
				}
				
				if (rank == 0)
					debug("Creating rotations for "+numSourceAz+" source azimuths, "+numSiteToSourceAz+" site-source azimuths, "
							+distances.length+" distances, and "+eventMatches.size()+" events");
				rotConfigs[s] = new RSQSimRotatedRupVariabilityConfig(
						catalog, gmpeSites, eventMatches, distances, numSourceAz, numSiteToSourceAz);
				if (rank == 0)
					rotConfigs[s].writeCSV(csvFiles[s]);
			} else {
				eventMatches = new ArrayList<>();
				for (int eventID : rotConfigs[s].getValues(Integer.class, Quantity.EVENT_ID)) {
					RSQSimEvent event = idEventMap.get(eventID);
					Preconditions.checkNotNull(event);
					eventMatches.add(event);
				}
				if (rank == 0)
					debug("Loaded "+eventMatches.size()+" matches for scenario: "+scenario);
				rotConfigs[s].setRuptures(allEvents);
			}
			
			File scenarioDir = new File(resultsDir, scenario.getPrefix());
			
			if (rank == 0)
				MPJ_BBP_Utils.waitOnDir(scenarioDir, 10, 2000);
			
			List<RotationSpec> rotations = rotConfigs[s].getRotations();
			if (rank == 0)
				debug("Created "+rotations.size()+" rotations");
			Table<Site, Float, File> siteDistDirTable = HashBasedTable.create();
			Table<File, Integer, File> eventDirTable = HashBasedTable.create();
			for (Site site : gmpeSites) {
				for (double distance : distances) {
					File siteDistDir = new File (scenarioDir, "site_"+site.getName()+"_dist_"+(float)distance);
					if (rank == 0)
						MPJ_BBP_Utils.waitOnDir(siteDistDir, 10, 2000);
					siteDistDirTable.put(site, (float)distance, siteDistDir);
					for (RSQSimEvent event : eventMatches) {
						File eventDir = new File (siteDistDir, "event_"+event.getID());
						eventDirTable.put(siteDistDir, event.getID(), eventDir);
						if (rank == 0)
							MPJ_BBP_Utils.waitOnDir(eventDir, 10, 2000);
					}
				}
			}
			for (RotationSpec rotation : rotConfigs[s].getRotations()) {
				File siteDistDir = siteDistDirTable.get(rotation.site, rotation.distance);
				File eventDir = eventDirTable.get(siteDistDir, rotation.eventID);
				Preconditions.checkNotNull(eventDir);
				configs.add(new Config(scenario, rotConfigs[s], rotation, eventDir));
			}
		}
	}
	
	private class Config {
		private final Scenario scenario;
		private final RSQSimRotatedRupVariabilityConfig config;
		private final RotationSpec rotation;
		private final File eventDir;
		
		public Config(Scenario scenario, RSQSimRotatedRupVariabilityConfig config, RotationSpec rotation, File eventDir) {
			this.scenario = scenario;
			this.config = config;
			this.rotation = rotation;
			this.eventDir = eventDir;
		}
	}
	
	@Override
	protected RSQSimEvent eventForIndex(int index) {
		Config config = configs.get(index);
		return config.config.getRotatedRupture(config.rotation);
	}

	@Override
	protected synchronized List<BBP_Site> sitesForIndex(int index) {
		Config config = configs.get(index);
		List<BBP_Site> sites = new ArrayList<>();
		sites.add(siteToBBPSites.get(config.rotation.site));
		return sites;
	}

	@Override
	protected File runDirForIndex(int index) {
		Config config = configs.get(index);
		return new File(config.eventDir, config.scenario.getPrefix()+"_"+config.rotation.getPrefix());
	}

	@Override
	protected int getNumTasks() {
		return configs.size();
	}
	
	static void addRotatedRupOptions(Options ops) {
		Option numSites = new Option("mr", "max-ruptures", true, "Maximum number of ruptures per scenario");
		numSites.setRequired(false);
		ops.addOption(numSites);
		
		Option skipYears = new Option("skip", "skip-years", true, "Skip the given number of years at the start");
		skipYears.setRequired(false);
		ops.addOption(skipYears);
		
		Option nSourceAz = new Option("nsrc", "num-source-az", true, "Num source centroid rotation azimuths");
		nSourceAz.setRequired(true);
		ops.addOption(nSourceAz);
		
		Option nSiteSourceAz = new Option("nsitesrc", "num-site-source-az", true, "Num site to source source rotation azimuths");
		nSiteSourceAz.setRequired(true);
		ops.addOption(nSiteSourceAz);
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		MPJ_BBP_Utils.addCommonOptions(ops, true, false, false, false);
		addCommonOptions(ops);
		
		MPJ_BBP_PartBSim.addPartB_ScenarioOptions(ops);
		
		addRotatedRupOptions(ops);
		
		return ops;
	}
	
	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_PartBSim.class);
			
			MPJ_BBP_RotatedRupVariabilityScenarioSim driver = new MPJ_BBP_RotatedRupVariabilityScenarioSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
