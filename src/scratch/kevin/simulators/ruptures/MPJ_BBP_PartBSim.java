package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.Location;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.ruptures.BBP_CatalogPartBValidationConfig.Scenario;

public class MPJ_BBP_PartBSim extends AbstractMPJ_BBP_MultiRupSim {
	
	private List<File> runDirs;
	private List<RSQSimEvent> events;
	private List<List<BBP_Site>> sites;
	
	private static final int MAX_SITES_PER_SIM = 10;

	private MPJ_BBP_PartBSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		int skipYears = -1;
		if (cmd.hasOption("skip-years"))
			skipYears = Integer.parseInt(cmd.getOptionValue("skip-years"));
		
		int numSites = Integer.parseInt(cmd.getOptionValue("num-sites"));
		Preconditions.checkState(numSites > 0);
		
		boolean randomAz = cmd.hasOption("random-azimuth");
		
		runDirs = new ArrayList<>();
		events = new ArrayList<>();
		sites = new ArrayList<>();
		
		double[] distances = BBP_CatalogPartBValidationConfig.DISTANCES;
		
		for (Scenario scenario : Scenario.values()) {
			File scenarioDir = new File(resultsDir, scenario.getPrefix());
			
			if (rank == 0)
				MPJ_BBP_Utils.waitOnDir(scenarioDir, 10, 2000);
			
			if (rank == 0)
				debug("Loading matches for scenario: "+scenario);
			
			List<RSQSimEvent> eventMatches = scenario.getMatches(catalog, skipYears);
			
			if (rank == 0)
				debug("Loaded "+eventMatches.size()+" matches for scenario: "+scenario);
			
			if (rank == 0)
				debug("Creating site lists with "+numSites+" sites per event, "+distances.length+" distances");
			for (RSQSimEvent event : eventMatches) {
				File eventDir = new File (scenarioDir, "event_"+event.getID());
				if (rank == 0)
					MPJ_BBP_Utils.waitOnDir(eventDir, 10, 2000);
				for (double distance : distances) {
					Location[] siteLocs = BBP_CatalogPartBValidationConfig.selectSitesSites(
							numSites, distance, randomAz, catalog, event);
					List<List<BBP_Site>> siteBundles = new ArrayList<>();
					List<BBP_Site> curSites = null;
					for (int i=0; i<numSites; i++) {
						if (curSites == null) {
							curSites = new ArrayList<>();
						} else if (curSites.size() == MAX_SITES_PER_SIM) {
							siteBundles.add(curSites);
							curSites = new ArrayList<>();
						}
						BBP_Site site = new BBP_Site("s"+i, siteLocs[i], vm.getVs30(),
								RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ);
						curSites.add(site);
					}
					if (curSites != null)
						siteBundles.add(curSites);
					for (List<BBP_Site> siteBundle : siteBundles) {
						String dirName = "dist_"+(float)distance+"_"+siteBundle.get(0).getName()
								+"_"+siteBundle.get(siteBundle.size()-1).getName();
						File runDir = new File(eventDir, dirName);
						
						runDirs.add(runDir);
						events.add(event);
						sites.add(siteBundle);
					}
				}
			}
			
			if (rank == 0) {
			}
		}
	}

	@Override
	RSQSimEvent eventForIndex(int index) {
		return events.get(index);
	}

	@Override
	List<BBP_Site> sitesForIndex(int index) {
		// already written to sites file
		return null;
	}

	@Override
	File runDirForIndex(int index) {
		return runDirs.get(index);
	}

	@Override
	protected int getNumTasks() {
		return events.size();
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		MPJ_BBP_Utils.addCommonOptions(ops, false, false, false, false);
		addCommonOptions(ops);
		
		Option numSites = new Option("ns", "num-sites", true, "Minimum magnitude");
		numSites.setRequired(true);
		ops.addOption(numSites);
		
		Option skipYears = new Option("skip", "skip-years", true, "Skip the given number of years at the start");
		skipYears.setRequired(false);
		ops.addOption(skipYears);
		
		Option randomAz = new Option("rand", "random-azimuth", false, "Flag for random azimuth");
		randomAz.setRequired(false);
		ops.addOption(randomAz);
		
		return ops;
	}
	
	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_PartBSim.class);
			
			MPJ_BBP_PartBSim driver = new MPJ_BBP_PartBSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}
	
}
