package scratch.kevin.simulators.ruptures.azimuthal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.ruptures.AbstractMPJ_BBP_MultiRupSim;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.MPJ_BBP_PartBSim;

public class MPJ_RSQSimAzimuthalSim extends AbstractMPJ_BBP_MultiRupSim {
	
	private List<Config> configs;

	private MPJ_RSQSimAzimuthalSim(CommandLine cmd) throws IOException {
		super(cmd);
		this.shuffle = false; // do them in order for better caching
		
		int skipYears = -1;
		if (cmd.hasOption("skip-years"))
			skipYears = Integer.parseInt(cmd.getOptionValue("skip-years"));
		
		int numRups = Integer.parseInt(cmd.getOptionValue("num-ruptures"));
		double buffer = Double.parseDouble(cmd.getOptionValue("buffer"));
		double spacing = Double.parseDouble(cmd.getOptionValue("spacing"));
		
		int sitesPerJob = cmd.hasOption("sites-per-job") ?
				Integer.parseInt(cmd.getOptionValue("sites-per-job")) : DEFAULT_MAX_SITES_PER_JOB;
		
		Scenario[] scenarios = MPJ_BBP_PartBSim.getScenarios(cmd);
		
		FilterMethod filter = MPJ_BBP_PartBSim.getFilterMethod(cmd);
		
		List<List<RuptureIdentifier>> scenarioIdens = new ArrayList<>();
		List<RuptureIdentifier> magIdens = new ArrayList<>();
		for (Scenario scenario : scenarios) {
			List<RuptureIdentifier> idens = scenario.getIdentifiers(catalog);
			scenarioIdens.add(idens);
			for (RuptureIdentifier iden : idens)
				if (iden instanceof MagRangeRuptureIdentifier)
					magIdens.add(iden);
		}
		
		if (rank == 0)
			debug("Loading events for "+scenarios.length+" scenarios");
		List<RSQSimEvent> allEvents = catalog.loader().skipYears(skipYears).matches(new LogicalOrRupIden(magIdens)).load();
		debug("Loaded "+allEvents.size()+" potential matches");
		
		configs = new ArrayList<>();
		
		for (int s=0; s<scenarios.length; s++) {
			Scenario scenario = scenarios[s];
			
			if (rank == 0)
				debug("Loading matches for scenario: "+scenario);
			
			List<RSQSimEvent> eventMatches = new LogicalAndRupIden(
					scenarioIdens.get(s)).getMatches(allEvents);
			
			if (rank == 0)
				debug("Loaded "+eventMatches.size()+" matches for scenario: "+scenario);
			
			if (eventMatches.size() > numRups) {
				eventMatches = filter.filter(eventMatches, numRups, catalog, scenario.getMagnitude());
				debug("trimmed down to max of "+eventMatches.size()+" ruptures");
			}
			
			if (eventMatches.isEmpty()) {
				if (rank == 0)
					debug("skipping...");
				continue;
			}
			
			File scenarioDir = new File(resultsDir, scenario.getPrefix());
			
			if (rank == 0)
				MPJ_BBP_Utils.waitOnDir(scenarioDir, 10, 2000);
			
			boolean positiveX = cmd.hasOption("positive-x");
			
			RSQSimAzimuthalSiteConfig config = new RSQSimAzimuthalSiteConfig(
					catalog, scenario, eventMatches, buffer, spacing, positiveX);
			
			int numSites = config.getGC2XYZ().size();
			
			for (RSQSimEvent event : eventMatches) {
				String prefix = scenario.getPrefix()+"_event_"+event.getID();
				if (numSites > sitesPerJob) {
					for (int start=0; start<numSites; start+=sitesPerJob) {
						int end = Integer.min(numSites, start + sitesPerJob);
						File runDir = new File(scenarioDir, prefix+"_s"+start+"_"+end);
						configs.add(new Config(config, event, runDir, start, end));
					}
				} else {
					File runDir = new File(scenarioDir, prefix);
					configs.add(new Config(config, event, runDir));
				}
			}
			
			config.writeGC2LocationCSV(new File(mainOutputDir, scenario.getPrefix()+"_locs.csv"));
			config.writeRupturesCSV(new File(mainOutputDir, scenario.getPrefix()+"_rups.csv"));
		}
	}
	
	private class Config {
		private final RSQSimAzimuthalSiteConfig config;
		private final RSQSimEvent event;
		private final File runDir;
		private final int startIndex, endIndex;
		
		public Config(RSQSimAzimuthalSiteConfig config, RSQSimEvent event, File runDir) {
			this(config, event, runDir, -1, -1);
		}
		
		public Config(RSQSimAzimuthalSiteConfig config, RSQSimEvent event, File runDir,
				int startIndex, int endIndex) {
			this.config = config;
			this.event = event;
			this.runDir = runDir;
			this.startIndex = startIndex;
			this.endIndex = endIndex;
		}
	}

	@Override
	protected RSQSimEvent eventForIndex(int index) {
		// needs to be oriented
		Config config = configs.get(index);
		return config.config.getOrientedEvent(config.event);
	}

	@Override
	protected List<BBP_Site> sitesForIndex(int index) {
		Config config = configs.get(index);
		List<BBP_Site> sites = config.config.getSites(config.event, vm);
		if (config.startIndex >= 0)
			return sites.subList(config.startIndex, config.endIndex);
		return sites;
	}

	@Override
	protected File runDirForIndex(int index) {
		return configs.get(index).runDir;
	}

	@Override
	protected int getNumTasks() {
		return configs.size();
	}
	
	static void addAzimuthalOptions(Options ops) {
		Option numSites = new Option("buf", "buffer", true, "Buffer around rupture surface in km");
		numSites.setRequired(true);
		ops.addOption(numSites);
		
		Option scenarios = new Option("scen", "scenarios", true, "BBP Part B Scenario names (comma separated). Default is all");
		scenarios.setRequired(true);
		ops.addOption(scenarios);
		
		Option numRups = new Option("nrup", "num-ruptures", true, "Num ruptures");
		numRups.setRequired(true);
		ops.addOption(numRups);
		
		Option spacing = new Option("sp", "spacing", true, "Site spacing in km");
		spacing.setRequired(true);
		ops.addOption(spacing);
		
		Option positiveX = new Option("px", "positive-x", false, "Only include positive RX sites");
		positiveX.setRequired(false);
		ops.addOption(positiveX);
		
		Option sitesPerJob = new Option("spj", "sites-per-job", true,
				"Maximum number of sites per job, default="+DEFAULT_MAX_SITES_PER_JOB);
		sitesPerJob.setRequired(false);
		ops.addOption(sitesPerJob);
	}
	
	static final int DEFAULT_MAX_SITES_PER_JOB = 50;
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		MPJ_BBP_Utils.addCommonOptions(ops, false, false, false, false);
		addCommonOptions(ops);
		
		Option filter = new Option("fil", "filter", true, "Filter method (for when there are more matches than allowed)");
		filter.setRequired(false);
		ops.addOption(filter);
		
		Option skipYears = new Option("skip", "skip-years", true, "Skip the given number of years at the start");
		skipYears.setRequired(false);
		ops.addOption(skipYears);
		
		addAzimuthalOptions(ops);
		
		return ops;
	}

	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_RSQSimAzimuthalSim.class);
			
			MPJ_RSQSimAzimuthalSim driver = new MPJ_RSQSimAzimuthalSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
