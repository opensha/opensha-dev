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
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.ruptures.AbstractMPJ_BBP_MultiRupSim;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.MPJ_BBP_PartBSim;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityMagDistPageGen.RuptureType;

public class MPJ_BBP_RotatedRupVariabilityMagDistSim extends AbstractMPJ_BBP_MultiRupSim {

	private List<BBP_Site> sites;
	private Map<Site, BBP_Site> siteToBBPSites;
	
	private Map<RotatedRupVariabilityConfig, File> csvsToWrite;
	
	private List<Config> configs;
	
	private MPJ_BBP_RotatedRupVariabilityMagDistSim(CommandLine cmd) throws IOException {
		super(cmd);
		this.shuffle = false; // do them in order for better caching
		
		int skipYears = -1;
		if (cmd.hasOption("skip-years"))
			skipYears = Integer.parseInt(cmd.getOptionValue("skip-years"));
		
		int maxRups = Integer.MAX_VALUE;
		if (cmd.hasOption("max-ruptures"))
			maxRups = Integer.parseInt(cmd.getOptionValue("max-ruptures"));
		
		int minRups = 0;
		if (cmd.hasOption("min-ruptures"))
			minRups = Integer.parseInt(cmd.getOptionValue("min-ruptures"));
		
		// only 1 site, doesn't matter where
		BBP_Site singleSite = RSQSimBBP_Config.getStandardSites().get(0);
		// fix the Vs30 value, probably not important but why not
		singleSite = new BBP_Site(singleSite.getName(), singleSite.getLoc(), vm.getVs30(),
				singleSite.getLoPassFreq(), singleSite.getHiPassFreq());
		sites = new ArrayList<>();
		sites.add(singleSite);
		siteToBBPSites = new HashMap<>();
		
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site bbpSite : sites) {
			Site site = bbpSite.buildGMPE_Site();
			siteToBBPSites.put(site, bbpSite);
			gmpeSites.add(site);
		}
		
		int numSourceAz = Integer.parseInt(cmd.getOptionValue("num-source-az"));
		int numSiteToSourceAz = Integer.parseInt(cmd.getOptionValue("num-site-source-az"));
		
		double[] distances = loadRange(cmd, "min-distance", "num-distance", "delta-distance");
		double[] mags = loadRange(cmd, "min-mag", "num-mag", "delta-mag");
		double deltaMag = mags[1]-mags[0];
		RuptureType[] rupTypes = getRuptureTypes(cmd);
		
		configs = new ArrayList<>();
		
		List<RSQSimEvent> allEvents = catalog.loader().skipYears(skipYears).minMag(mags[0]-deltaMag).load();
		if (rank == 0)
			debug("Loaded a total of "+allEvents.size()+" events");
		
		if (rank == 0)
			csvsToWrite = new HashMap<>();
		
		for (int r=0; r<rupTypes.length; r++) {
			RuptureType rupType = rupTypes[r];
			
			for (int m=0; m<mags.length; m++) {
				double mag = mags[m];
				File configCSV = new File(mainOutputDir, "rotation_config_"+rupType.getMagPrefix(mag)+".csv");
				
				RotatedRupVariabilityConfig config;
				
				if (configCSV.exists()) {
					// load it in
					if (rank == 0)
						debug("Loading matches for rupture type: "+rupType+", M="+(float)mag+" from "+configCSV.getAbsolutePath());
					config = RotatedRupVariabilityConfig.loadCSV(catalog, configCSV, allEvents, gmpeSites);
					
					// validate it
					int numEvents = config.getValues(Integer.class, Quantity.EVENT_ID).size();
					int testNumSourceAz = config.getValues(Float.class, Quantity.SOURCE_AZIMUTH).size();
					int testNumSiteSourceAz = config.getValues(Float.class, Quantity.SITE_TO_SOURTH_AZIMUTH).size();
					int testNumDist = config.getValues(Float.class, Quantity.DISTANCE).size();
					int testNumSites = config.getValues(Site.class, Quantity.SITE).size();
					if (rank == 0)
						debug("Loaded "+config.getRotations().size()+" rotations for "+testNumSites+" sites, "+numEvents+" events, "
								+testNumSourceAz+" sourceAz, "+testNumSiteSourceAz+" siteSourceAz, "+testNumDist+" dists");
					Preconditions.checkState(numEvents <= maxRups && numEvents >= minRups, "Bad event count: %s", numEvents);
					Preconditions.checkState(testNumSites == 1, "Expected 1 site: %s", testNumSites);
					Preconditions.checkState(testNumSourceAz == numSourceAz, "Expected %s source az: %s",
							numSourceAz, testNumSourceAz);
					Preconditions.checkState(testNumSiteSourceAz == numSiteToSourceAz, "Expected %s site-source az: %s",
							numSiteToSourceAz, testNumSiteSourceAz);
					Preconditions.checkState(testNumDist == distances.length, "Expected %s dists: %s",
							distances.length, testNumDist);
				} else {
					if (rank == 0)
						debug("Loading matches for rupture type: "+rupType+", M="+(float)mag);
					
					List<RSQSimEvent> eventMatches = rupType.getMatches(catalog, allEvents, mag, deltaMag);
					
					if (rank == 0)
						debug("Loaded "+eventMatches.size()+" matches for scenario: "+rupType+", M="+(float)mag);
					
					if (eventMatches.size() > maxRups) {
						eventMatches = BBP_PartBValidationConfig.getClosestMagMatches(mag, eventMatches, maxRups);
						if (rank == 0)
							debug("trimmed down to max of "+eventMatches.size()+" ruptures");
					}
					
					if (eventMatches.isEmpty() || eventMatches.size() < minRups) {
						if (rank == 0)
							debug("skipping, not enough ruptures! have "+eventMatches.size());
						continue;
					}
					
					if (rank == 0)
						debug("Creating rotations for "+numSourceAz+" source azimuths, "+numSiteToSourceAz+" site-source azimuths, "
								+distances.length+" distances, and "+eventMatches.size()+" events");
					config = new RotatedRupVariabilityConfig(
							catalog, gmpeSites, eventMatches, distances, numSourceAz, numSiteToSourceAz);
					List<RotationSpec> rotations = config.getRotations();
					if (rank == 0)
						debug("Created "+rotations.size()+" rotations");
					if (rank == 0)
						csvsToWrite.put(config, new File(mainOutputDir, "rotation_config_"+rupType.getMagPrefix(mag)+".csv"));
				}
				
				File scenarioDir = new File(resultsDir, rupType.getMagPrefix(mag));
				
				if (rank == 0)
					MPJ_BBP_Utils.waitOnDir(scenarioDir, 10, 2000);
				
				List<Integer> eventIDs = config.getValues(Integer.class, Quantity.EVENT_ID);
				
				Table<Site, Float, File> siteDistDirTable = HashBasedTable.create();
				Table<File, Integer, File> eventDirTable = HashBasedTable.create();
				for (Site site : gmpeSites) {
					for (double distance : distances) {
						File siteDistDir = new File (scenarioDir, "site_"+site.getName()+"_dist_"+(float)distance);
						if (rank == 0)
							MPJ_BBP_Utils.waitOnDir(siteDistDir, 10, 2000);
						siteDistDirTable.put(site, (float)distance, siteDistDir);
						for (Integer eventID : eventIDs) {
							File eventDir = new File (siteDistDir, "event_"+eventID);
							eventDirTable.put(siteDistDir, eventID, eventDir);
							if (rank == 0)
								MPJ_BBP_Utils.waitOnDir(eventDir, 10, 2000);
						}
					}
				}
				for (RotationSpec rotation : config.getRotations()) {
					File siteDistDir = siteDistDirTable.get(rotation.site, rotation.distance);
					File eventDir = eventDirTable.get(siteDistDir, rotation.eventID);
					Preconditions.checkNotNull(eventDir);
					configs.add(new Config(rupType, mag, config, rotation, eventDir));
				}
			}
		}
	}
	
	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		if (rank == 0 && csvsToWrite != null) {
			debug("flushing CSVs to disk...");
			for (RotatedRupVariabilityConfig config : csvsToWrite.keySet())
				config.writeCSV(csvsToWrite.get(config));
			debug("done flushing CSVs to disk.");
			csvsToWrite = null;
		}
		super.calculateBatch(batch);
	}

	private static double[] loadRange(CommandLine cmd, String minArg, String numArg, String deltaArg) {
		double min = Double.parseDouble(cmd.getOptionValue(minArg));
		int num = Integer.parseInt(cmd.getOptionValue(numArg));
		Preconditions.checkArgument(num > 1, "Num must be >1");
		double delta = Double.parseDouble(cmd.getOptionValue(deltaArg));
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(min, num, delta);
		double[] ret = new double[func.size()];
		
		for (int i=0; i<ret.length; i++)
			ret[i] = func.getX(i);
		
		return ret;
	}
	
	static RuptureType[] getRuptureTypes(CommandLine cmd) {
		if (cmd.hasOption("rupture-type")) {
			String str = cmd.getOptionValue("rupture-type");
			String[] strs = str.split(",");
			RuptureType[] rupTypes = new RuptureType[strs.length];
			for (int i=0; i<strs.length; i++)
				rupTypes[i] = RuptureType.valueOf(strs[i]);
			return rupTypes;
		}
		return RuptureType.values();
	}
	
	private class Config {
		private final RuptureType rupType;
		private final double mag;
		private final RotatedRupVariabilityConfig config;
		private final RotationSpec rotation;
		private final File eventDir;
		
		public Config(RuptureType rupType, double mag, RotatedRupVariabilityConfig config, RotationSpec rotation, File eventDir) {
			this.rupType = rupType;
			this.mag = mag;
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
		return new File(config.eventDir, config.rupType.getMagPrefix(config.mag)+"_"+config.rotation.getPrefix());
	}

	@Override
	protected int getNumTasks() {
		return configs.size();
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		MPJ_BBP_Utils.addCommonOptions(ops, false, false, false, false);
		addCommonOptions(ops);
		
		MPJ_BBP_RotatedRupVariabilityScenarioSim.addRotatedRupOptions(ops);
		
		Option numSites = new Option("minr", "min-ruptures", true, "Minimum number of ruptures per scenario");
		numSites.setRequired(false);
		ops.addOption(numSites);
		
		Option minDistOp = new Option("md", "min-distance", true, "Minimum distance (km)");
		minDistOp.setRequired(true);
		ops.addOption(minDistOp);
		
		Option mumDistOp = new Option("nd", "num-distance", true, "Num distnaces");
		mumDistOp.setRequired(true);
		ops.addOption(mumDistOp);
		
		Option deltaDistOp = new Option("dd", "delta-distance", true, "Delta distnace (km)");
		deltaDistOp.setRequired(true);
		ops.addOption(deltaDistOp);
		
		Option minMagOp = new Option("mm", "min-mag", true, "Minimum magnitude");
		minMagOp.setRequired(true);
		ops.addOption(minMagOp);
		
		Option mumMagOp = new Option("nm", "num-mag", true, "Num magnitudes");
		mumMagOp.setRequired(true);
		ops.addOption(mumMagOp);
		
		Option deltaMagOp = new Option("dm", "delta-mag", true, "Delta magnitude");
		deltaMagOp.setRequired(true);
		ops.addOption(deltaMagOp);
		
		Option rupTypeOp = new Option("rt", "rupture-type", true, "Rupture type(s), comma separated");
		rupTypeOp.setRequired(true);
		ops.addOption(rupTypeOp);
		
		return ops;
	}
	
	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_PartBSim.class);
			
			MPJ_BBP_RotatedRupVariabilityMagDistSim driver = new MPJ_BBP_RotatedRupVariabilityMagDistSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
