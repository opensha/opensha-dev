package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.SRF_PointData;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.MPJ_BBP_RupGenSim;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class MPJ_BBP_CatalogSim extends MPJTaskCalculator {
	
	private VelocityModel vm;
	private Method method;
	
	private double dt;
	private SRFInterpolationMode interp;
	
	private RSQSimCatalog catalog;
	private List<RSQSimEvent> events;
	
	public static final double CUTOFF_DIST = 200d;
	
	private File mainOutputDir;
	
	private List<BBP_Site> sites;
	private List<RegionIden> siteRegIdens;
	
	private File resultsDir;
	private boolean doHF = true;
	private boolean keepSRFs = false;
	private File bbpDataDir = null;
	private int bundleSize;
	
	private ExecutorService exec;

	public MPJ_BBP_CatalogSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		vm = VelocityModel.valueOf(cmd.getOptionValue("vm"));
		method = Method.valueOf(cmd.getOptionValue("method"));
		File catalogDir = new File(cmd.getOptionValue("catalog-dir"));
		Preconditions.checkState(catalogDir.exists());
		File sitesFile = new File(cmd.getOptionValue("sites-file"));
		Preconditions.checkState(sitesFile.exists());
		dt = Double.parseDouble(cmd.getOptionValue("time-step"));
		interp = SRFInterpolationMode.valueOf(cmd.getOptionValue("srf-interp"));
		mainOutputDir = new File(cmd.getOptionValue("output-dir"));
		resultsDir = new File(mainOutputDir, "results");
		if (rank == 0) {
			Preconditions.checkState((mainOutputDir.exists() && mainOutputDir.isDirectory()) || mainOutputDir.mkdir());
			Preconditions.checkState((resultsDir.exists() && resultsDir.isDirectory()) || resultsDir.mkdir());
		}
		if (resultsDir.exists()) {
			try {
				// for restarts
				bundleSize = detectBundleSize(resultsDir);
			} catch (FileNotFoundException e) {}
		}
		doHF = !cmd.hasOption("no-hf");
		
		if (cmd.hasOption("bbp-data-dir")) {
			bbpDataDir = new File(cmd.getOptionValue("bbp-data-dir"));
			if (rank == 0)
				debug("BBP data dir: "+bbpDataDir.getAbsolutePath());
		}
		
		sites = BBP_Site.readFile(sitesFile);
		siteRegIdens = new ArrayList<>();
		for (BBP_Site site : sites)
			siteRegIdens.add(new RegionIden(new Region(site.getLoc(), CUTOFF_DIST)));
		
		// load the catalog
		catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
				null, null, null, null, null);
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
		
		if (rank == 0)
			postBatchHook = new MasterZipHook();
		
		exec = Executors.newFixedThreadPool(getNumThreads());
	}
	
	private File getZipFile(int eventID) {
		File runDir = getRunDir(resultsDir, eventID, false, bundleSize);
		return new File(runDir.getParentFile(), runDir.getName()+".zip");
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
	
	private class MasterZipHook extends MPJ_BBP_Utils.MasterZipHook {

		public MasterZipHook() {
			super(new File(resultsDir.getParentFile(), resultsDir.getName()+".zip"),
					new File(resultsDir.getParentFile(), resultsDir.getName()+"_rotD.zip"));
		}

		@Override
		protected void debug(String message) {
			MPJ_BBP_CatalogSim.this.debug(message);
		}

		@Override
		protected void abortAndExit(int status) {
			MPJTaskCalculator.abortAndExit(status);
		}

		@Override
		protected File getSimZipFile(int index) {
			int eventID = events.get(index).getID();
			return getZipFile(eventID);
		}
		
	}

	@Override
	protected int getNumTasks() {
		return events.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		ArrayList<Future<?>> futures = new ArrayList<>();
		for (int index : batch) {
			futures.add(exec.submit(new Task(index)));
			Thread.sleep(100);
		}
		for (Future<?> f : futures)
			f.get();
	}
	
	private class Task implements Runnable {
		
		private int index;
		
		private Task(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			try {
				RSQSimEvent event = events.get(index);
				int eventID = event.getID();
				File zipFile = getZipFile(eventID);
				if (zipFile.exists()) {
					debug(eventID+" is already done, skipping");
					return;
				}
				File runDir = getRunDir(eventID);
				RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
				
				// write SRF
				debug("bulding/writing SRF for "+eventID);
				List<SRF_PointData> srfPoints = RSQSimSRFGenerator.buildSRF(func, event.getAllElements(), dt, interp);
				File srfFile = new File(runDir, runDir.getName()+".srf");
				SRF_PointData.writeSRF(srfFile, srfPoints, 1d);
				
				// write SRC
				debug("bulding/writing SRC for "+eventID);
				BBP_PlanarSurface bbpSurface = RSQSimBBP_Config.estimateBBP_PlanarSurface(event);
				BBP_SourceFile bbpSource = RSQSimBBP_Config.buildBBP_Source(event, func, bbpSurface, 12345);
				File srcFile = new File(runDir, runDir.getName()+".src");
				bbpSource.writeToFile(srcFile);
				
				// write sites file
				List<BBP_Site> mySites = new ArrayList<>();
				for (int i=0; i<sites.size(); i++)
					if (siteRegIdens.get(i).isMatch(event))
						mySites.add(sites.get(i));
				Preconditions.checkState(!mySites.isEmpty(), "Should be at least one site for each loaded rupture!");
				File sitesFile = new File(runDir, "sites.stl");
				BBP_Site.writeToFile(sitesFile, mySites);
				
				// run BBP
				debug("running BBP for "+eventID);
				BBP_Wrapper wrapper = new BBP_Wrapper(vm, method, srcFile, null, srfFile, sitesFile, runDir);
				wrapper.setDoHF(doHF);
				wrapper.setDataOnly(true);
				wrapper.setBBPDataDir(bbpDataDir);
				wrapper.run();
				
				if (!keepSRFs)
					srfFile.delete();
				
				// zip it
				FileUtils.createZipFile(zipFile, runDir, true);
				FileUtils.deleteRecursive(runDir);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		exec.shutdown();
		if (rank == 0) {
			debug("waiting for any post batch hook operations to finish");
			((MasterZipHook)postBatchHook).shutdown();
			debug("post batch hook done");
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJ_BBP_Utils.addCommonOptions(MPJTaskCalculator.createOptions(), true, false, false, false);
		
		Option dt = new Option("dt", "time-step", true, "SRF time step");
		dt.setRequired(true);
		ops.addOption(dt);
		
		Option interp = new Option("interp", "srf-interp", true, "SRF interpolation mode");
		interp.setRequired(true);
		ops.addOption(interp);
		
		Option catalogDir = new Option("cdir", "catalog-dir", true, "RSQSim catalog dir");
		catalogDir.setRequired(true);
		ops.addOption(catalogDir);
		
		Option mag = new Option("mag", "min-mag", true, "Minimum magnitude");
		mag.setRequired(true);
		ops.addOption(mag);
		
		Option skipYears = new Option("skip", "skip-years", true, "Skip the given number of years at the start");
		skipYears.setRequired(false);
		ops.addOption(skipYears);
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
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
