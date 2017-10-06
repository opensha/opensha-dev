package scratch.kevin.bbp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
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

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class MPJ_BBP_Wrapper extends MPJTaskCalculator {
	
	private VelocityModel vm;
	private Method method;
	
	private File srcFile;
	private File sitesFile;
	private File mainOutputDir;
	private int numSims;
	
	private File resultsDir;
	private boolean doHF = true;
	
	private List<File> individualSiteFiles;
	
	private ExecutorService exec;

	public MPJ_BBP_Wrapper(CommandLine cmd) {
		super(cmd);
		
		vm = VelocityModel.valueOf(cmd.getOptionValue("vm"));
		method = Method.valueOf(cmd.getOptionValue("method"));
		srcFile = new File(cmd.getOptionValue("src-file"));
		Preconditions.checkState(srcFile.exists());
		sitesFile = new File(cmd.getOptionValue("sites-file"));
		Preconditions.checkState(sitesFile.exists());
		mainOutputDir = new File(cmd.getOptionValue("output-dir"));
		resultsDir = new File(mainOutputDir, "results");
		if (rank == 0) {
			Preconditions.checkState((mainOutputDir.exists() && mainOutputDir.isDirectory()) || mainOutputDir.mkdir());
			Preconditions.checkState((resultsDir.exists() && resultsDir.isDirectory()) || resultsDir.mkdir());
		}
		numSims = Integer.parseInt(cmd.getOptionValue("num-sims"));
		int numSeeds = numSims;
		doHF = !cmd.hasOption("no-hf");
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		
		if (cmd.hasOption("split-sites")) {
			individualSiteFiles = new ArrayList<>();
			File subSiteDir = new File(mainOutputDir, "sites_individual");
			if (rank == 0)
				Preconditions.checkState((subSiteDir.exists() && subSiteDir.isDirectory()) || subSiteDir.mkdir());
			try {
				for (String line : Files.readLines(sitesFile, Charset.defaultCharset())) {
					line = line.trim();
					if (line.startsWith("#") || line.isEmpty())
						continue;
					String[] split = line.split("\\s+");
					Preconditions.checkState(split.length > 2);
					String siteName = split[2];
					File subSiteFile = new File(subSiteDir, siteName+".stl");
					individualSiteFiles.add(subSiteFile);
					if (rank == 0) {
						debug("Writing individual site file: "+subSiteFile.getName());
						FileWriter fw = new FileWriter(subSiteFile);
						
						fw.write("#SLong    SLat     RSN   Vs30(m/s) LoPass_Freq(Hz) HiPass_Freq(Hz)\n");
						fw.write(line+"\n");
						
						fw.close();
					}
				}
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			numSims *= individualSiteFiles.size();
		}
		
		// initialize results dirs, and create seeds
		if (rank == 0) {
			Random r = new Random(System.nanoTime());
			for (int i=0; i<numSeeds; i++) {
				File runDir = getRunDir(i);
				long seed = r.nextInt(Short.MAX_VALUE);
				try {
					BBP_Wrapper.writeSrcFileNewSeed(srcFile, seed, runDir);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		}
	}
	
	private File getRunDir(int index) {
		File runDir = new File(resultsDir, "run_"+index);
		Preconditions.checkState(runDir.exists() || runDir.mkdir(),
				"Run dir could not be created: %s", runDir.getAbsolutePath());
		return runDir;
	}

	@Override
	protected int getNumTasks() {
		return numSims;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		ArrayList<Future<?>> futures = new ArrayList<>();
		for (int index : batch) {
			File sitesFile = this.sitesFile;
			String siteName = null;
			if (individualSiteFiles != null) {
				int siteIndex = index % individualSiteFiles.size();
				index = index / individualSiteFiles.size();
				sitesFile = individualSiteFiles.get(siteIndex);
				siteName = sitesFile.getName();
				siteName = siteName.substring(0, siteName.indexOf("."));
			}
			File subDir = getRunDir(index);
			if (isDone(subDir, siteName)) {
				debug(index+" is already done, skipping");
				continue;
			}
			File srcFile = locateSrcFile(subDir);
			Preconditions.checkNotNull(srcFile, "No src file in %s", subDir.getAbsolutePath());
			BBP_Wrapper wrapper = new BBP_Wrapper(vm, method, srcFile, null, null, sitesFile, subDir);
			wrapper.setDoHF(doHF);
			if (individualSiteFiles != null) {
				wrapper.setXMLFileName("inputs_"+siteName+".xml");
				wrapper.setScriptFileName("run_"+siteName+".sh");
			}
			futures.add(exec.submit(wrapper));
			Thread.sleep(5000);
		}
		for (Future<?> f : futures)
			f.get();
	}
	
	private boolean isDone(File dir, String doneCheckStr) {
		for (File file : dir.listFiles())
			if (file.getName().endsWith(".rd50"))
				return doneCheckStr == null || file.getName().contains(doneCheckStr);
		return false;
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
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option vmOp = new Option("vm", "vm", true, "Velocity model");
		vmOp.setRequired(true);
		ops.addOption(vmOp);
		
		Option methodOp = new Option("m", "method", true, "BBP method");
		methodOp.setRequired(true);
		ops.addOption(methodOp);
		
		Option srcFile = new Option("src", "src-file", true, "Source file");
		srcFile.setRequired(true);
		ops.addOption(srcFile);
		
		Option sitesFile = new Option("sites", "sites-file", true, "Sites file");
		sitesFile.setRequired(true);
		ops.addOption(sitesFile);
		
		Option outputDir = new Option("o", "output-dir", true, "Output dir");
		outputDir.setRequired(true);
		ops.addOption(outputDir);
		
		Option numSims = new Option("n", "num-sims", true, "Number of simulations");
		numSims.setRequired(true);
		ops.addOption(numSims);
		
		Option splitSites = new Option("sp", "split-sites", false, "Flag to split sites into individual runs");
		splitSites.setRequired(false);
		ops.addOption(splitSites);
		
		Option noHF = new Option("nhf", "no-hf", false, "Flag to disable high-frequency");
		noHF.setRequired(false);
		ops.addOption(noHF);
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_Wrapper.class);
			
			MPJ_BBP_Wrapper driver = new MPJ_BBP_Wrapper(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
