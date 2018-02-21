package scratch.kevin.bbp;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class MPJ_BBP_RupGenSim extends MPJTaskCalculator {
	
	private VelocityModel vm;
	private Method method;
	
	private File srcFile;
	private File sitesFile;
	private File mainOutputDir;
	private int numSims;
	
	private File resultsDir;
	private boolean doHF = true;
	
	private List<BBP_Site> sites;
	private List<File> individualSiteFiles;
	
	private File bbpDataDir = null;
	
	private ExecutorService exec;
	
	private boolean keepSRFs = false;

	public MPJ_BBP_RupGenSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		vm = VelocityModel.valueOf(cmd.getOptionValue("vm"));
		method = Method.valueOf(cmd.getOptionValue("method"));
		srcFile = new File(cmd.getOptionValue("src-file"));
		Preconditions.checkState(srcFile.exists());
		sitesFile = new File(cmd.getOptionValue("sites-file"));
		Preconditions.checkState(sitesFile.exists());
		sites = BBP_Site.readFile(sitesFile);
		mainOutputDir = new File(cmd.getOptionValue("output-dir"));
		resultsDir = new File(mainOutputDir, "results");
		if (rank == 0) {
			Preconditions.checkState((mainOutputDir.exists() && mainOutputDir.isDirectory()) || mainOutputDir.mkdir());
			Preconditions.checkState((resultsDir.exists() && resultsDir.isDirectory()) || resultsDir.mkdir());
		}
		numSims = Integer.parseInt(cmd.getOptionValue("num-sims"));
		int numSeeds = numSims;
		doHF = !cmd.hasOption("no-hf");
		if (cmd.hasOption("bbp-data-dir")) {
			bbpDataDir = new File(cmd.getOptionValue("bbp-data-dir"));
			if (!bbpDataDir.exists())
				bbpDataDir.mkdir();
			if (rank == 0)
				debug("BBP data dir: "+bbpDataDir.getAbsolutePath());
		}
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		
		if (cmd.hasOption("split-sites")) {
			individualSiteFiles = new ArrayList<>();
			File subSiteDir = new File(mainOutputDir, "sites_individual");
			if (rank == 0)
				Preconditions.checkState((subSiteDir.exists() && subSiteDir.isDirectory()) || subSiteDir.mkdir());
			for (BBP_Site site : sites) {
				File subSiteFile = new File(subSiteDir, site.getName()+".stl");
				individualSiteFiles.add(subSiteFile);
				if (rank == 0) {
					ArrayList<BBP_Site> singleSite = new ArrayList<>();
					singleSite.add(site);
					BBP_Site.writeToFile(subSiteFile, singleSite);
				}
			}
			numSims *= individualSiteFiles.size();
		}
		
		// initialize results dirs, and create seeds
		if (rank == 0) {
			postBatchHook = new MasterZipHook();
			
			BBP_SourceFile sourceFile = BBP_SourceFile.readFile(srcFile);
			Random r = new Random(System.nanoTime());
			for (int i=0; i<numSeeds; i++) {
				int seed = r.nextInt(Short.MAX_VALUE);
				sourceFile.setSeed(seed);
				int[] indexes;
				if (individualSiteFiles != null) {
					indexes = new int[sites.size()];
					for (int s=0; s<sites.size(); s++)
						indexes[s] = i*sites.size() + s;
				} else {
					indexes = new int[] { i };
				}
				for (int index : indexes) {
					File runDir = getRunDir(index, true);
					sourceFile.writeToFile(new File(runDir, srcFile.getName()));
				}
			}
		}
	}
	
	private File getRunDir(int index, boolean create) {
		File runDir;
		if (individualSiteFiles != null) {
			int siteIndex = index % individualSiteFiles.size();
			index = index / individualSiteFiles.size();
			String siteName = sites.get(siteIndex).getName();
			runDir = new File(resultsDir, "run_"+index+"_"+siteName);
		} else {
			runDir = new File(resultsDir, "run_"+index);
		}
		if (create)
			MPJ_BBP_Utils.waitOnDir(runDir, 10, 2000);
		return runDir;
	}
	
	private File getZipFile(int index) {
		File runDir = getRunDir(index, false);
		return new File(runDir.getParentFile(), runDir.getName()+".zip");
	}

	@Override
	protected int getNumTasks() {
		return numSims;
	}
	
	private class MasterZipHook extends MPJ_BBP_Utils.MasterZipHook {

		public MasterZipHook() {
			super(new File(resultsDir.getParentFile(), resultsDir.getName()+".zip"),
					new File(resultsDir.getParentFile(), resultsDir.getName()+"_rotD.zip"));
		}

		@Override
		protected void debug(String message) {
			MPJ_BBP_RupGenSim.this.debug(message);
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

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		ArrayList<Future<?>> futures = new ArrayList<>();
		for (int index : batch) {
			Task task = new Task(index);
			futures.add(exec.submit(task));
			Thread.sleep(100);
		}
		for (Future<?> f : futures)
			f.get();
	}
	
	private class Task implements Runnable {
		
		private int index;

		public Task(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			try {
				File runDir = getRunDir(index, true);
				File zipFile = getZipFile(index);
				if (zipFile.exists()) {
					debug(index+" is already done, skipping");
					return;
				}
				File mySitesFile = sitesFile;
				String siteName = null;
				if (individualSiteFiles != null) {
					int siteIndex = index % individualSiteFiles.size();
					index = index / individualSiteFiles.size();
					mySitesFile = individualSiteFiles.get(siteIndex);
					siteName = sites.get(siteIndex).getName();
				}
				File srcFile = locateSrcFile(runDir);
				Preconditions.checkNotNull(srcFile, "No src file in %s", runDir.getAbsolutePath());
				BBP_Wrapper wrapper = new BBP_Wrapper(vm, method, srcFile, null, null, mySitesFile, runDir);
				wrapper.setDoHF(doHF);
				wrapper.setDoFAS(true);
				wrapper.setDoRotD100(true);
				wrapper.setDoRotD50(false);
				wrapper.setDataOnly(true);
				if (individualSiteFiles != null) {
					wrapper.setXMLFileName("inputs_"+siteName+".xml");
					wrapper.setScriptFileName("run_"+siteName+".sh");
				}
				
				// run BBP
				debug("running BBP for "+index);
				wrapper.setDoHF(doHF);
				wrapper.setBBPDataDir(bbpDataDir);
				wrapper.run();
				
				if (!keepSRFs) {
					for (File file : runDir.listFiles())
						if (file.getName().endsWith(".srf"))
							file.delete();
				}
				
				// zip it
				FileUtils.createZipFile(zipFile, runDir, true);
				FileUtils.deleteRecursive(runDir);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
	}
	
	private File locateSrcFile(File dir) {
		File srcFile = null;
		for (File file : dir.listFiles()) {
			if (file.getName().endsWith(".src")) {
				Preconditions.checkState(srcFile == null, "Multiple .src files exist in %s", dir.getAbsolutePath());
				srcFile = file;
			}
		}
		return srcFile;
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
		Options ops = MPJ_BBP_Utils.addCommonOptions(MPJTaskCalculator.createOptions(), true, true, false, false);
		
		Option numSims = new Option("n", "num-sims", true, "Number of simulations");
		numSims.setRequired(true);
		ops.addOption(numSims);
		
		Option splitSites = new Option("sp", "split-sites", false, "Flag to split sites into individual runs");
		splitSites.setRequired(false);
		ops.addOption(splitSites);
		
		return ops;
	}
	
	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_RupGenSim.class);
			
			MPJ_BBP_RupGenSim driver = new MPJ_BBP_RupGenSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
