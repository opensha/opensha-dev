package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

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
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import mpi.MPI;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class MPJ_BBP_CatalogSim extends MPJTaskCalculator {
	
	private VelocityModel vm;
	private Method method;
	
	private double dt;
	private SRFInterpolationMode interp;
	
	private double timeScalarFactor = 1d;
	private boolean velocityScale = false;
	
	private RSQSimCatalog catalog;
	private List<RSQSimEvent> events;
	
	public static final double CUTOFF_DIST = 200d;
	
	private File mainOutputDir;
	
	private List<BBP_Site> sites;
	private List<RegionIden> siteRegIdens;
	
	private File resultsDir;
	private File resultsScratchDir;
	private boolean doHF = true;
	private boolean keepSRFs = false;
	
	private File bbpEnvFile = null;
	private File bbpDataDir = null;
	private File bbpGFDir = null;
	
	private int bundleSize;
	private int numRG = 0;
	
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
		
		File sharedScratch = null;
		if (cmd.hasOption("shared-scratch-dir")) {
			sharedScratch = new File(cmd.getOptionValue("shared-scratch-dir"));
			if (rank == 0)
				Preconditions.checkState(sharedScratch.exists() || sharedScratch.mkdir());
		}
		
		File nodeScratch = null;
		if (cmd.hasOption("node-scratch-dir")) {
			nodeScratch = new File(cmd.getOptionValue("node-scratch-dir"));
			Preconditions.checkState(nodeScratch.exists() || nodeScratch.mkdir());
		}
		
		if (nodeScratch != null) {
			// copy catalog data over to node scratch
			catalogDir = copyCatalogDir(catalogDir, nodeScratch);
			resultsScratchDir = new File(nodeScratch, "results_tmp");
			Preconditions.checkState(resultsScratchDir.exists() || resultsScratchDir.mkdir());
		} else if (sharedScratch != null) {
			// copy catalog data over to shared scratch
			if (rank == 0)
				copyCatalogDir(catalogDir, sharedScratch);
			resultsScratchDir = new File(sharedScratch, "results_tmp");
			if (rank == 0)
				MPJ_BBP_Utils.waitOnDir(resultsScratchDir, 10, 2000);
			MPI.COMM_WORLD.Barrier();
			catalogDir = new File(sharedScratch, catalogDir.getName());
		}
		
		if (resultsDir.exists()) {
			try {
				// for restarts
				bundleSize = detectBundleSize(resultsDir);
			} catch (FileNotFoundException e) {}
		}
		doHF = !cmd.hasOption("no-hf");
		
		if (cmd.hasOption("bbp-env")) {
			bbpEnvFile = new File(cmd.getOptionValue("bbp-env"));
			if (rank == 0) {
				debug("BBP env file: "+bbpEnvFile.getAbsolutePath());
				Preconditions.checkState(bbpEnvFile.exists(), "Env file doesn't exist: %s", bbpEnvFile.getAbsolutePath());
			}
		}
		
		if (cmd.hasOption("bbp-data-dir")) {
			bbpDataDir = new File(cmd.getOptionValue("bbp-data-dir"));
			if (!bbpDataDir.exists())
				bbpDataDir.mkdir();
			if (rank == 0)
				debug("BBP data dir: "+bbpDataDir.getAbsolutePath());
		}
		
		if (cmd.hasOption("bbp-gf-dir")) {
			bbpGFDir = new File(cmd.getOptionValue("bbp-gf-dir"));
			Preconditions.checkState(bbpGFDir.exists());
			if (rank == 0)
				debug("BBP GF dir: "+bbpGFDir.getAbsolutePath());
		}
		
		sites = BBP_Site.readFile(sitesFile);
		siteRegIdens = new ArrayList<>();
		for (BBP_Site site : sites)
			siteRegIdens.add(new RegionIden(new Region(site.getLoc(), CUTOFF_DIST)));
		
		// load the catalog
		catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
				null, null, null, FaultModels.FM3_1, DeformationModels.GEOLOGIC); // TODO
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
		
		if (cmd.hasOption("rup-gen-sims")) {
			numRG = Integer.parseInt(cmd.getOptionValue("rup-gen-sims"));
			if (rank == 0)
				debug("also doing "+numRG+" rupture generation simulations each");
		}
		
		if (cmd.hasOption("time-scalar"))
			timeScalarFactor = Double.parseDouble(cmd.getOptionValue("time-scalar"));
		if (cmd.hasOption("velocity-scale"))
			velocityScale = true;
		
		if (rank == 0)
			postBatchHook = new MasterZipHook();
		
		exec = Executors.newFixedThreadPool(getNumThreads());
	}
	
	private File copyCatalogDir(File catDir, File scratchDir) throws IOException {
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
		
		debug("copying "+filesToCopy.size()+" catalog files to "+destDir.getAbsolutePath());
		for (File f : filesToCopy) {
			File destFile = new File(destDir, f.getName());
			if (destFile.exists() && destFile.length() == f.length())
				// skip copy, already exists
				continue;
			Files.copy(f, destFile);
		}
		
		return destDir;
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
			super(numRG > 0 ? null : new File(resultsDir.getParentFile(), resultsDir.getName()+".zip"),
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
	
	private static void validateZip(File zipFile) throws Exception {
		ZipFile zip = new ZipFile(zipFile);
		Preconditions.checkState(zip.entries().hasMoreElements(), "No entries!");
		zip.close();
	}
	
	private class Task implements Runnable {
		
		private int index;
		
		private Task(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			RSQSimEvent event = null;
			try {
				event = events.get(index);
				int eventID = event.getID();
				File zipFile = getZipFile(eventID);
				if (zipFile.exists()) {
					try {
						validateZip(zipFile);
						debug(eventID+" is already done, skipping");
						return;
					} catch (Exception e) {
						e.printStackTrace();
						debug(eventID+" zip file exists, but won't open or is empty. Re-running. "+zipFile.getAbsolutePath());
					}
				}
				File runDir;
				if (resultsScratchDir == null) {
					runDir = getRunDir(eventID);
				} else {
					runDir = new File(resultsScratchDir, "event_"+eventID);
					MPJ_BBP_Utils.waitOnDir(runDir, 10, 2000);
				}
				RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
				
				if (timeScalarFactor != 1d)
					func = func.getTimeScaledFunc(timeScalarFactor, velocityScale);
				
				// write SRF
				debug("bulding/writing SRF for "+eventID);
				List<SRF_PointData> srfPoints = RSQSimSRFGenerator.buildSRF(func, event.getAllElements(), dt, interp);
				File srfFile = new File(runDir, runDir.getName()+".srf");
				SRF_PointData.writeSRF(srfFile, srfPoints, 1d);
				
				// write SRC
				debug("bulding/writing SRC for "+eventID);
				BBP_PlanarSurface bbpSurface;
				if (numRG > 0 && RSQSimBBP_Config.U3_SURFACES)
					bbpSurface = RSQSimBBP_Config.planarEquivalentU3Surface(catalog, event);
				else
					bbpSurface = RSQSimBBP_Config.estimateBBP_PlanarSurface(event);
				BBP_SourceFile bbpSource = RSQSimBBP_Config.buildBBP_Source(event, bbpSurface, 12345);
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
				wrapper.setDoFAS(false);
				wrapper.setDoRotD100(true);
				wrapper.setDoRotD50(false);
				wrapper.setBBPEnvFile(bbpEnvFile);
				wrapper.setBBPDataDir(bbpDataDir);
				wrapper.setBBPGFDir(bbpGFDir);
				wrapper.run();
				debug("DONE running BBP for "+eventID);
				
				if (!keepSRFs)
					srfFile.delete();
				
				if (numRG > 0) {
					Random r = new Random();
					for (int i=0; i<numRG; i++) {
						File rgDir = new File(runDir, "rup_gen_"+i);
						Preconditions.checkState(rgDir.exists() || rgDir.mkdir());
						debug("running BBP for "+eventID+", RG "+i+"/"+numRG);
						wrapper = new BBP_Wrapper(vm, method, srcFile, (long)r.nextInt(Short.MAX_VALUE), null, sitesFile, rgDir);
						wrapper.setDoHF(doHF);
						wrapper.setDataOnly(true);
						wrapper.setDoFAS(false);
						wrapper.setDoRotD100(true);
						wrapper.setDoRotD50(false);
						wrapper.setBBPEnvFile(bbpEnvFile);
						wrapper.setBBPDataDir(bbpDataDir);
						wrapper.setBBPGFDir(bbpGFDir);
						wrapper.run();
					}
				}
				
				// zip it
				FileUtils.createZipFile(zipFile, runDir, true);
				try {
					validateZip(zipFile);
				} catch (Exception e) {
					e.printStackTrace();
					debug(eventID+" zip file didn't validate, trying again");
					try {
						Thread.sleep(5000);
					} catch (InterruptedException e1) {}
					FileUtils.createZipFile(zipFile, runDir, true);
					try {
						validateZip(zipFile);
					} catch (Exception e1) {
						debug(eventID+" zip failed again after rety, bailing");
						throw ExceptionUtils.asRuntimeException(e1);
					}
				}
				FileUtils.deleteRecursive(runDir);
				debug("DONE cleanup for "+eventID);
			} catch (Exception e) {
				String message = "Exception for index "+index;
				if (event != null)
					message += ", event "+event.getID();
				message += ": "+e.getMessage();
				throw new RuntimeException(message, e);
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
		
		Option gp = new Option("rgs", "rup-gen-sims", true, "If supplied, do this amount of simulations with a rupture "
				+ "generator for each rupture. Only RotD50 will be archived.");
		gp.setRequired(false);
		ops.addOption(gp);
		
		Option nodeScratchDir = new Option("nsdir", "node-scratch-dir", true,
				"Node-local scratch directory to reduce I/O on network disks. Transition file will be copied here.");
		nodeScratchDir.setRequired(false);
		ops.addOption(nodeScratchDir);
		
		Option sharedScratchDir = new Option("ssdir", "shared-scratch-dir", true,
				"Shared scratch directory to reduce I/O on network disks. "
				+ "Transition file will be copied here if no node-local scratch. Unzipped results will be stored here.");
		sharedScratchDir.setRequired(false);
		ops.addOption(sharedScratchDir);
		
		Option timeScalar = new Option("ts", "time-scalar", true,
				"Time scalar for slip/time functions. 1.5 will increase rupture propagation velocities by 50%");
		timeScalar.setRequired(false);
		ops.addOption(timeScalar);
		
		Option velScale = new Option("vs", "velocity-scale", false,
				"Works with time scalar. If supplied, velocities will also be scaled (each slip event is faster and shorter). "
				+ "Otherwise slip events are sooner, but the same speed/duration.");
		velScale.setRequired(false);
		ops.addOption(velScale);
		
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
