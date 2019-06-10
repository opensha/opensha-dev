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

public abstract class AbstractMPJ_BBP_MultiRupSim extends MPJTaskCalculator {
	
	protected final VelocityModel vm;
	protected final Method method;
	
	private double dt;
	private SRFInterpolationMode interp;
	
	private double timeScalarFactor = 1d;
	private boolean velocityScale = false;
	
	protected final RSQSimCatalog catalog;
	
	protected final File mainOutputDir;
	
	protected final File resultsDir;
	private File resultsScratchDir;
	private boolean doHF = true;
	private boolean keepSRFs = false;
	
	private File bbpEnvFile = null;
	private File bbpDataDir = null;
	private File bbpGFDir = null;
	
	private int numRG = 0;
	
	private ExecutorService exec;

	public AbstractMPJ_BBP_MultiRupSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		vm = VelocityModel.valueOf(cmd.getOptionValue("vm"));
		method = Method.valueOf(cmd.getOptionValue("method"));
		File catalogDir = new File(cmd.getOptionValue("catalog-dir"));
		if (!catalogDir.exists()) {
			String catalogDirName = catalogDir.getName();
			if (rank == 0)
				debug("Catalog dir doesn't exist, searching for: "+catalogDirName);
			catalogDir = RSQSimCatalog.locateCatalog(catalogDirName);
			Preconditions.checkNotNull(catalogDir, "Couldn't locate catalog: "+catalogDirName);
			if (rank == 0)
				debug("Located catalog: "+catalogDir.getAbsolutePath());
		}
		Preconditions.checkState(catalogDir.exists());
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
		
		// load the catalog
		catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
				null, null, null, FaultModels.FM3_1, DeformationModels.GEOLOGIC); // TODO make option
		
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
	
	private File getZipFile(int index) {
		File runDir = getRunDir(index, false);
		return new File(runDir.getParentFile(), runDir.getName()+".zip");
	}
	
	private File getRunDir(int index, boolean create) {
		File runDir = runDirForIndex(index);
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
			AbstractMPJ_BBP_MultiRupSim.this.debug(message);
		}

		@Override
		protected void abortAndExit(int status) {
			MPJTaskCalculator.abortAndExit(status);
		}

		@Override
		protected File getSimZipFile(int index) {
			return getZipFile(index);
		}
		
	}
	
	protected abstract RSQSimEvent eventForIndex(int index);
	protected abstract List<BBP_Site> sitesForIndex(int index);
	protected abstract File runDirForIndex(int index);

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
				event = eventForIndex(index);
				int eventID = event.getID();
				File zipFile = getZipFile(index);
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
				File runDir = getRunDir(index, resultsScratchDir == null);
				if (resultsScratchDir != null) {
					// create a unique directory on scratch
					
					// need to include any parent directory names
					String dirName = runDir.getName();
					File parentFile = runDir.getParentFile();
					while (parentFile != null && !parentFile.equals(resultsDir) && !parentFile.getName().isEmpty()) {
						dirName = parentFile.getName()+"_"+dirName;
						parentFile = parentFile.getParentFile();
					}
					runDir = new File(resultsScratchDir, dirName);
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
				File sitesFile = new File(runDir, "sites.stl");
				List<BBP_Site> mySites = sitesForIndex(index);
				Preconditions.checkState(!mySites.isEmpty(), "Should be at least one site for each loaded rupture!");
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
	
	public static Options addCommonOptions(Options ops) {
		Option dt = new Option("dt", "time-step", true, "SRF time step");
		dt.setRequired(true);
		ops.addOption(dt);
		
		Option interp = new Option("interp", "srf-interp", true, "SRF interpolation mode");
		interp.setRequired(true);
		ops.addOption(interp);
		
		Option catalogDir = new Option("cdir", "catalog-dir", true, "RSQSim catalog dir");
		catalogDir.setRequired(true);
		ops.addOption(catalogDir);
		
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

}
