package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.dom4j.DocumentException;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;

public class MPJ_ETAS_CatalogEALCalculator extends MPJTaskCalculator {
	
	private List<File> resultsFiles;
	
	private Map<AttenRelRef, Double> imrWeightsMap;
	private FaultSystemSolution baSol;
	private CompoundFaultSystemSolution cfss;
	private FaultSystemSolution trueMeanSol;
	private Map<LogicTreeBranch, List<Integer>> branchMappings;
	private List<File> dataDirs;
	private boolean triggeredOnly;
	private boolean allSubDurations = false;
	
	private double[] durations = { 1d/365.25, 7d/365.25, 30/365.25, 1d, 10d };
	
	private double maxX = 200;

	public MPJ_ETAS_CatalogEALCalculator(CommandLine cmd) throws IOException, DocumentException {
		super(cmd);
		
		File inputFile = new File(cmd.getOptionValue("input-file"));
		Preconditions.checkArgument(inputFile.exists());
		
		File simsDir = new File(cmd.getOptionValue("sims-dir"));
		Preconditions.checkArgument(simsDir.exists());
		
		String resultsFileName = cmd.getOptionValue("results-file-name");
		
		resultsFiles = Lists.newArrayList();
		for (String line : Files.readLines(inputFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#"))
				continue;
			File dir = new File(simsDir, line);
			Preconditions.checkState(dir.exists(), "Input dir doesn't exist: %s", dir.getAbsolutePath());
			
			File resultsFile = new File(dir, resultsFileName);
			Preconditions.checkState(resultsFile.exists(), "Results file doesn't exist: %s", resultsFile.getAbsolutePath());
			resultsFiles.add(resultsFile);
		}
		
		File baSolFile = new File(cmd.getOptionValue("branch-avg-sol"));
		Preconditions.checkArgument(baSolFile.exists());
		baSol = FaultSystemIO.loadSol(baSolFile);
		
		File compoundSolFile = new File(cmd.getOptionValue("compound-sol"));
		Preconditions.checkArgument(compoundSolFile.exists());
		cfss = CompoundFaultSystemSolution.fromZipFile(compoundSolFile);
		
		File trueMeanSolFile = new File(cmd.getOptionValue("true-mean-sol"));
		Preconditions.checkArgument(trueMeanSolFile.exists());
		trueMeanSol = FaultSystemIO.loadSol(trueMeanSolFile);
		branchMappings = TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile);
		
		dataDirs = Lists.newArrayList();
		String dataDirsStr = cmd.getOptionValue("data-dirs");
		for (String dirStr : HazardCurvePlotter.commaSplit(dataDirsStr)) {
			File dir = new File(dirStr);
			Preconditions.checkState(dir.exists());
			dataDirs.add(dir);
		}
		
		Preconditions.checkState(!dataDirs.isEmpty());
		if (dataDirs.get(0).getName().contains("san-bernardino") || dataDirs.get(0).getName().contains("coachella"))
			maxX = 50;
		
		imrWeightsMap = Maps.newHashMap();
		imrWeightsMap.put(AttenRelRef.CB_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.CY_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.ASK_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.BSSA_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.IDRISS_2014, 0.12);
		
		triggeredOnly = cmd.hasOption("triggered-only");
		allSubDurations = cmd.hasOption("sub-durations");
	}
	
	private static String xAxisLabel = "$ (Billions)";
	
	// portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
	private static double thousandsToBillions = 1d/1e6;
	
	private static double inflationScalar = 1d/0.9d;
//	private static double deltaX = 1e6/inflationScalar;
	private static double deltaX = 1; // now delta is after scaling, so 1 here means 1 billion
	
	private static double xAxisScale = thousandsToBillions*inflationScalar;
	
	private static FaultModels fm = FaultModels.FM3_1;

	@Override
	protected int getNumTasks() {
		return resultsFiles.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int index : batch) {
			File resultsFile = resultsFiles.get(index);
			
			ETAS_CatalogEALCalculator.calculate(resultsFile, triggeredOnly, xAxisLabel, maxX, deltaX, xAxisScale,
					dataDirs, imrWeightsMap, fm, baSol, cfss, trueMeanSol, branchMappings, durations, allSubDurations);
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option inputFile = new Option("i", "input-file", true, "Input text file with list of directory names");
		inputFile.setRequired(true);
		ops.addOption(inputFile);
		
		Option simDir = new Option("sd", "sims-dir", true,
				"Directory containing all of the scenario directories listed in the input file");
		simDir.setRequired(true);
		ops.addOption(simDir);
		
		Option resultsFileName = new Option("rf", "results-file-name", true, "Results file name");
		resultsFileName.setRequired(true);
		ops.addOption(resultsFileName);
		
		Option baSolFile = new Option("ba", "branch-avg-sol", true, "Branch average solution file");
		baSolFile.setRequired(true);
		ops.addOption(baSolFile);
		
		Option compoundSolFile = new Option("cs", "compound-sol", true, "Compound solution file");
		compoundSolFile.setRequired(true);
		ops.addOption(compoundSolFile);
		
		Option trueMeanSolFile = new Option("tms", "true-mean-sol", true, "True mean solution file");
		trueMeanSolFile.setRequired(true);
		ops.addOption(trueMeanSolFile);
		
		Option dataDirs = new Option("d", "data-dirs", true, "Data directories (comma separated)");
		dataDirs.setRequired(true);
		ops.addOption(dataDirs);
		
		Option trigOnly = new Option("tr", "triggered-only", false, "If supplied, only triggered ruptures will be considered");
		trigOnly.setRequired(false);
		ops.addOption(trigOnly);
		
		Option allSubDurs = new Option("sub", "sub-durations", false, "If supplied, all sub durations will be calculated");
		allSubDurs.setRequired(false);
		ops.addOption(allSubDurs);
		
		return ops;
	}

	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_ETAS_CatalogEALCalculator.class);
			
			args = cmd.getArgs();
			
			MPJ_ETAS_CatalogEALCalculator driver = new MPJ_ETAS_CatalogEALCalculator(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
