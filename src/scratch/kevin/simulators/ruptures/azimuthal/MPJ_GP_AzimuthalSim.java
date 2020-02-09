package scratch.kevin.simulators.ruptures.azimuthal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.sha.simulators.srf.SRF_PointData;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.ruptures.AbstractMPJ_BBP_Sim;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.MPJ_BBP_PartBSim;

public class MPJ_GP_AzimuthalSim extends AbstractMPJ_BBP_Sim {
	
	private List<Config> configs;

	private MPJ_GP_AzimuthalSim(CommandLine cmd) throws IOException {
		super(cmd);
		this.shuffle = false; // do them in order for better caching
		
		int numRups = Integer.parseInt(cmd.getOptionValue("num-ruptures"));
		double buffer = Double.parseDouble(cmd.getOptionValue("buffer"));
		double spacing = Double.parseDouble(cmd.getOptionValue("spacing"));
		
		int sitesPerJob = cmd.hasOption("sites-per-job") ?
				Integer.parseInt(cmd.getOptionValue("sites-per-job")) :
					MPJ_RSQSimAzimuthalSim.DEFAULT_MAX_SITES_PER_JOB;
		
		if (rank == 0)
			debug("buffer: "+buffer+", spacing: "+spacing);
		
		double patchArea = Double.parseDouble(cmd.getOptionValue("patch-area"));
		
		Scenario[] scenarios = MPJ_BBP_PartBSim.getScenarios(cmd);
		
		configs = new ArrayList<>();
		
		for (int s=0; s<scenarios.length; s++) {
			Scenario scenario = scenarios[s];
			
			File scenarioDir = new File(resultsDir, scenario.getPrefix());
			
			if (rank == 0)
				MPJ_BBP_Utils.waitOnDir(scenarioDir, 10, 2000);
			
			boolean positiveX = cmd.hasOption("positive-x");
			
			GPAzimuthalSiteConfig config = new GPAzimuthalSiteConfig(
					scenario, numRups, patchArea, buffer, spacing, positiveX);
			
			int numSites = config.getGC2XYZ().size();
			
			for (int i=0; i<numRups; i++) {
				String prefix = scenario.getPrefix()+"_event_"+i;
				if (numSites > sitesPerJob) {
					for (int start=0; start<numSites; start+=sitesPerJob) {
						int end = Integer.min(numSites, start + sitesPerJob);
						File runDir = new File(scenarioDir, prefix+"_s"+start+"_"+end);
						configs.add(new Config(config, i, runDir, start, end));
					}
				} else {
					File runDir = new File(scenarioDir, prefix);
					configs.add(new Config(config, i, runDir));
				}
			}
			
			config.writeGC2LocationCSV(new File(mainOutputDir, scenario.getPrefix()+"_locs.csv"));
			config.writeRupturesCSV(new File(mainOutputDir, scenario.getPrefix()+"_rups.csv"));
		}
	}
	
	private class Config {
		private final GPAzimuthalSiteConfig config;
		private final int eventID;
		private final File runDir;
		private final int startIndex, endIndex;
		
		public Config(GPAzimuthalSiteConfig config, int eventID, File runDir) {
			this(config, eventID, runDir, -1, -1);
		}
		
		public Config(GPAzimuthalSiteConfig config, int eventID, File runDir,
				int startIndex, int endIndex) {
			this.config = config;
			this.eventID = eventID;
			this.runDir = runDir;
			this.startIndex = startIndex;
			this.endIndex = endIndex;
		}
	}

	@Override
	protected List<BBP_Site> sitesForIndex(int index) {
		Config config = configs.get(index);
		List<BBP_Site> sites = config.config.getSites(config.eventID, vm);
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

	@Override
	protected List<SRF_PointData> getSRFPoints(int index) throws IOException {
		return null;
	}

	@Override
	protected BBP_SourceFile getBBPSource(int index) {
		Config config = configs.get(index);
		return config.config.getBBP_Source(config.eventID);
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		MPJ_BBP_Utils.addCommonOptions(ops, false, false, false, false);
		addCommonOptions(ops);
		
		Option patchArea = new Option("pa", "patch-area", true, "Patch area in square kilometers");
		patchArea.setRequired(true);
		ops.addOption(patchArea);
		
		MPJ_RSQSimAzimuthalSim.addAzimuthalOptions(ops);
		
		return ops;
	}

	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_GP_AzimuthalSim.class);
			
			MPJ_GP_AzimuthalSim driver = new MPJ_GP_AzimuthalSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
