package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.Region;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class MPJ_BBP_CatalogSim extends AbstractMPJ_BBP_MultiRupSim {
	
	private List<RSQSimEvent> events;
	
	public static final double CUTOFF_DIST = 200d;
	
	private List<BBP_Site> sites;
	private List<RegionIden> siteRegIdens;
	
	private int bundleSize;

	private MPJ_BBP_CatalogSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		if (resultsDir.exists()) {
			try {
				// for restarts
				bundleSize = detectBundleSize(resultsDir);
			} catch (FileNotFoundException e) {}
		}
		
		File sitesFile = new File(cmd.getOptionValue("sites-file"));
		Preconditions.checkState(sitesFile.exists());
		sites = BBP_Site.readFile(sitesFile);
		siteRegIdens = new ArrayList<>();
		for (BBP_Site site : sites)
			siteRegIdens.add(new RegionIden(new Region(site.getLoc(), CUTOFF_DIST)));
		
		// load the catalog
		Loader loader = catalog.loader().hasTransitions();
		if (cmd.hasOption("min-mag"))
			loader.minMag(Double.parseDouble(cmd.getOptionValue("min-mag")));
		if (cmd.hasOption("skip-years"))
			loader.skipYears(Integer.parseInt(cmd.getOptionValue("skip-years")));
		loader.matches(new LogicalOrRupIden(siteRegIdens));
		events = loader.load();
		
		int idSpan = events.get(events.size()-1).getID() - events.get(0).getID();
		double numEventsPerID = (double)events.size()/(double)(idSpan);
		if (bundleSize <= 0) {
			// new run, calculate
			double rangeFor1K = 1000d/numEventsPerID;
			double roundedLog = Math.round(Math.log10(rangeFor1K));
			bundleSize = (int)Math.pow(10, roundedLog);
			if (bundleSize < 1000)
				bundleSize = 1000;
		}
		if (rank == 0) {
			debug("loaded "+events.size()+" events");
			debug("initializing bundle dirs with bundleSize="+bundleSize
					+". Expected per bundle: "+(int)(bundleSize*numEventsPerID));
			// initialize parent directories
			for (SimulatorEvent e : events)
				getRunParentDir(e.getID());
			
			// wait a few seconds to make sure dir creation propagates through the NFS
			try {
				Thread.sleep(5000);
			} catch (InterruptedException e1) {}
		}
	}
	
	@Override
	protected RSQSimEvent eventForIndex(int index) {
		return events.get(index);
	}

	@Override
	protected List<BBP_Site> sitesForIndex(int index) {
		RSQSimEvent event = eventForIndex(index);
		List<BBP_Site> mySites = new ArrayList<>();
		for (int i=0; i<sites.size(); i++)
			if (siteRegIdens.get(i).isMatch(event))
				mySites.add(sites.get(i));
		return mySites;
	}

	@Override
	protected File runDirForIndex(int index) {
		int eventID = eventForIndex(index).getID();
		return getRunDir(eventID);
	}

	private File getRunDir(int eventID) {
		return getRunDir(resultsDir, eventID, true, bundleSize);
	}
	
	private File getRunParentDir(int eventID) {
		return getRunParentDir(resultsDir, eventID, true, bundleSize);
	}
	
	static int detectBundleSize(File resultsDir) throws FileNotFoundException {
		for (File dir : resultsDir.listFiles()) {
			if (dir.isDirectory() && dir.getName().startsWith("events_")) {
				String[] split = dir.getName().split("_");
				int end = Integer.parseInt(split[split.length-1]);
				int start = Integer.parseInt(split[split.length-2]);
				return end - start;
			}
		}
		throw new FileNotFoundException("Didn't find any event bundles");
	}
	
	static File getRunParentDir(File resultsDir, int eventID, boolean create, int bundleSize) {
		int eventBase = bundleSize*(int)(eventID / bundleSize);
		File parentDir = new File(resultsDir, "events_"+eventBase+"_"+(eventBase+bundleSize));
		if (create)
			MPJ_BBP_Utils.waitOnDir(parentDir, 10, 2000);
		return parentDir;
	}
	
	static File getRunDir(File resultsDir, int eventID, boolean create, int bundleSize) {
		File parentDir = getRunParentDir(resultsDir, eventID, create, bundleSize);
		File runDir = new File(parentDir, "event_"+eventID);
		if (create)
			MPJ_BBP_Utils.waitOnDir(runDir, 10, 2000);
		return runDir;
	}

	@Override
	protected int getNumTasks() {
		return events.size();
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		MPJ_BBP_Utils.addCommonOptions(ops, true, false, false, false);
		addCommonOptions(ops);
		
		Option mag = new Option("mag", "min-mag", true, "Minimum magnitude");
		mag.setRequired(true);
		ops.addOption(mag);
		
		Option skipYears = new Option("skip", "skip-years", true, "Skip the given number of years at the start");
		skipYears.setRequired(false);
		ops.addOption(skipYears);
		
		return ops;
	}
	
	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_CatalogSim.class);
			
			MPJ_BBP_CatalogSim driver = new MPJ_BBP_CatalogSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
