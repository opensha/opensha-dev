package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
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
import com.google.common.base.Stopwatch;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;

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

public abstract class AbstractMPJ_BBP_Sim extends MPJTaskCalculator {
	
	protected final VelocityModel vm;
	protected final Method method;
	
	protected final File mainOutputDir;
	
	protected final File resultsDir;
	protected File resultsScratchDir;
	protected File nodeScratch;
	protected File sharedScratch;
	private boolean doHF = true;
	private boolean keepSRFs = false;
	
	private File bbpEnvFile = null;
	private File bbpDataDir = null;
	private File bbpGFDir = null;
	
	protected ExecutorService exec;
	
	private HashSet<Integer> alreadyDones;
	
	protected int numRG = 0;

	public AbstractMPJ_BBP_Sim(CommandLine cmd) throws IOException {
		super(cmd);
		
		vm = VelocityModel.valueOf(cmd.getOptionValue("vm"));
		method = Method.valueOf(cmd.getOptionValue("method"));
		
		mainOutputDir = new File(cmd.getOptionValue("output-dir"));
		resultsDir = new File(mainOutputDir, "results");
		if (rank == 0) {
			Preconditions.checkState((mainOutputDir.exists() && mainOutputDir.isDirectory()) || mainOutputDir.mkdir());
			if (!resultsDir.exists())
				// don't need to search for already dones, first time
				alreadyDones = new HashSet<>();
			Preconditions.checkState((resultsDir.exists() && resultsDir.isDirectory()) || resultsDir.mkdir());
		}
		
		sharedScratch = null;
		if (cmd.hasOption("shared-scratch-dir")) {
			sharedScratch = new File(cmd.getOptionValue("shared-scratch-dir"));
			if (rank == 0)
				Preconditions.checkState(sharedScratch.exists() || sharedScratch.mkdir());
		}
		
		nodeScratch = null;
		if (cmd.hasOption("node-scratch-dir")) {
			nodeScratch = new File(cmd.getOptionValue("node-scratch-dir"));
			Preconditions.checkState(nodeScratch.exists() || nodeScratch.mkdir());
		}
		
		if (nodeScratch != null) {
			resultsScratchDir = new File(nodeScratch, "results_tmp");
			Preconditions.checkState(resultsScratchDir.exists() || resultsScratchDir.mkdir());
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
		
		if (cmd.hasOption("node-gf-dir")) {
			File localDir = new File(cmd.getOptionValue("node-gf-dir"));
			if (size >= 5) {
				// stagger them by 10s each
				try {
					Thread.sleep(10000l*rank);
				} catch (InterruptedException e) {}
			}
			rsyncGFs(localDir, bbpEnvFile, bbpGFDir, vm, rank, hostname);
			bbpGFDir = localDir;
		}

		if (cmd.hasOption("rup-gen-sims")) {
			numRG = Integer.parseInt(cmd.getOptionValue("rup-gen-sims"));
			if (rank == 0)
				debug("also doing "+numRG+" rupture generation simulations each");
		}
		
		if (rank == 0)
			postBatchHook = new MasterZipHook();
		
		exec = Executors.newFixedThreadPool(getNumThreads());
	}
	
	@Override
	protected synchronized Collection<Integer> getDoneIndexes() {
		if (alreadyDones == null && rank == 0) {
			debug("searching for simulations which are already complete...");
			int num = getNumTasks();
			alreadyDones = new HashSet<>();
			int validateFails = 0;
			Stopwatch totWatch = Stopwatch.createStarted();
			Stopwatch getZipWatch = Stopwatch.createUnstarted();
			Stopwatch validateWatch = Stopwatch.createUnstarted();
			for (int i=0; i<num; i++) {
				getZipWatch.start();
				File zipFile = getZipFile(i);
				getZipWatch.stop();
				if (zipFile.exists()) {
					try {
						validateWatch.start();
						validateZip(zipFile);
						validateWatch.stop();
						alreadyDones.add(i);
					} catch (Exception e) {
						validateFails++;
					}
				}
			}
			totWatch.stop();
			debug("Took "+totWatch.elapsed(TimeUnit.MILLISECONDS)/1000f+" s to search for old sims");
			debug("Took "+getZipWatch.elapsed(TimeUnit.MILLISECONDS)/1000f+" s to get zip files");
			debug("Took "+validateWatch.elapsed(TimeUnit.MILLISECONDS)/1000f+" s to validate old sim zip files");
			if (!alreadyDones.isEmpty()) {
				debug("found "+alreadyDones.size()+" simulations which are already complete (will skip dispatching)");
				if (validateFails > 0)
					debug("also found "+validateFails+" zip files which did not validate and will be ru-run");
				// send to post-batch hook for processing
				int[] doneBatch = Ints.toArray(alreadyDones);
				postBatchHook.batchProcessed(doneBatch, -1);
			} else {
				debug("didn't find any already complete simulations, will run everything");
			}
		}
		return alreadyDones;
	}

	public static void rsyncGFs(File localDir, File bbpEnvFile, File bbpGFDir, VelocityModel vm, int rank, String hostname)
			throws IOException {
		debug(rank, hostname, "rsyncing GFs to "+localDir.getAbsolutePath());
		MPJ_BBP_Utils.waitOnDir(localDir, 5, 100);
		File script = new File(localDir, "rsync.sh");
		FileWriter fw = new FileWriter(script);
		
		fw.write("#!/bin/bash\n");
		fw.write("\n");
		if (bbpEnvFile != null) {
			fw.write(". "+bbpEnvFile.getAbsolutePath()+"\n");
		} else {
			fw.write("if [ -f ~/.bash_profile ]; then\n");
			fw.write("    . ~/.bash_profile\n");
			fw.write("fi\n");
		}
		fw.write("\n");
		if (bbpGFDir != null) {
			fw.write("export BBP_GF_DIR="+bbpGFDir+"\n");
		} else {
			fw.write("if [ -z ${BBP_GF_DIR+x} ];then\n");
			fw.write("    echo \"BBP_GF_DIR is undefined\"\n");
			fw.write("    exit 2\n");
			fw.write("fi\n");
		}
		fw.write("VM_DIR=${BBP_GF_DIR}/"+vm.getDirName()+"\n");
		fw.write("if [ ! -e $VM_DIR ];then\n");
		fw.write("    echo \"VM_DIR is doesn't exist: $VM_DIR\"\n");
		fw.write("    exit 2\n");
		fw.write("fi\n");
		fw.write("\n");
		fw.write("VM_ZIP=\"${VM_DIR}.zip\"\n");
		fw.write("if [ -e $VM_ZIP ];then\n");
		fw.write("    echo \"Zip file exits, extracting: $VM_ZIP\"\n");
		fw.write("    cd "+localDir.getAbsolutePath()+"\n");
		fw.write("    unzip $VM_ZIP > /dev/null\n");
		fw.write("    RET=$?\n");
		fw.write("else\n");
		fw.write("    echo \"rsync-ing from $VM_DIR\"\n");
		fw.write("    rsync -a $VM_DIR "+localDir.getAbsolutePath()+" > /dev/null\n");
		fw.write("    RET=$?\n");
		fw.write("fi\n");
		fw.write("\n");
		fw.write("echo \"DONE. exitcode: $RET\"\n");
		fw.write("\n");
		fw.write("exit $RET\n");
		
		fw.close();
		
		int ret = BBP_Wrapper.runScript(script);
		
		Preconditions.checkState(ret == 0, "Couldn't copy GFs to node local storage");
		
		debug(rank, hostname, "done rsyncing GFs");
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
			AbstractMPJ_BBP_Sim.this.debug(message);
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
	
	protected abstract List<SRF_PointData> getSRFPoints(int index) throws IOException;
	protected abstract BBP_SourceFile getBBPSource(int index);
	
	private class Task implements Runnable {
		
		private int index;
		
		private Task(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			try {
				File zipFile = getZipFile(index);
				if (zipFile.exists()) {
					try {
						validateZip(zipFile);
						debug(index+" is already done, skipping");
						return;
					} catch (Exception e) {
						e.printStackTrace();
						debug(index+" zip file exists, but won't open or is empty. Re-running. "+zipFile.getAbsolutePath());
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
				
				List<SRF_PointData> srfPoints = getSRFPoints(index);
				File srfFile = new File(runDir, runDir.getName()+".srf");
				if (srfFile != null)
					SRF_PointData.writeSRF(srfFile, srfPoints, 1d);
				
				// write SRC
				BBP_SourceFile bbpSource = getBBPSource(index);
				File srcFile = new File(runDir, runDir.getName()+".src");
				bbpSource.writeToFile(srcFile);
				
				// write sites file
				File sitesFile = new File(runDir, "sites.stl");
				List<BBP_Site> mySites = sitesForIndex(index);
				Preconditions.checkState(!mySites.isEmpty(), "Should be at least one site for each loaded rupture!");
				BBP_Site.writeToFile(sitesFile, mySites);
				
				// run BBP
				debug("running BBP for index "+index);
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
				debug("DONE running BBP for "+index);
				
				if (!keepSRFs)
					srfFile.delete();
				
				if (numRG > 0) {
					Random r = new Random();
					for (int i=0; i<numRG; i++) {
						File rgDir = new File(runDir, "rup_gen_"+i);
						Preconditions.checkState(rgDir.exists() || rgDir.mkdir());
						debug("running BBP for "+index+", RG "+i+"/"+numRG);
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
					debug(index+" zip file didn't validate, trying again");
					try {
						Thread.sleep(5000);
					} catch (InterruptedException e1) {}
					FileUtils.createZipFile(zipFile, runDir, true);
					try {
						validateZip(zipFile);
					} catch (Exception e1) {
						debug(index+" zip failed again after rety, bailing");
						throw ExceptionUtils.asRuntimeException(e1);
					}
				}
				FileUtils.deleteRecursive(runDir);
				debug("DONE cleanup for "+index);
			} catch (Exception e) {
				String message = "Exception for index: "+e.getMessage();
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
		
		return ops;
	}

}
