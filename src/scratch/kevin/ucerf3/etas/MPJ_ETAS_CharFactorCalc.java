package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentMap;
import java.util.zip.ZipException;

import mpi.MPI;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class MPJ_ETAS_CharFactorCalc extends MPJTaskCalculator {
	
	private transient File outputDir;
	
	private transient FaultSystemSolutionFetcher cfss;
	
	// fault model: branch: factors
	private Map<FaultModels, List<double[]>> charFactorsMap;
	private Map<FaultModels, List<Double>> weightsMap;
	
	private transient ConcurrentMap<FaultModels, List<FaultSectionPrefData>> subSectsMap;
	
	private transient List<LogicTreeBranch> branches;
	
	private transient BranchWeightProvider weightProv;

	public MPJ_ETAS_CharFactorCalc(CommandLine cmd, File compoundFile, File outputDir) throws ZipException, IOException {
		super(cmd);
		
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55;
		weightProv = new APrioriBranchWeightProvider();
		
		cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
		if (cmd.hasOption("random-sample")) {
			int num = Integer.parseInt(cmd.getOptionValue("random-sample"));
			cfss = FaultSystemSolutionFetcher.getRandomSample(cfss, num);
		}
		
		this.outputDir = outputDir;
		
		charFactorsMap = Maps.newHashMap();
		weightsMap = Maps.newHashMap();
		
		branches = Lists.newArrayList(cfss.getBranches());
		subSectsMap = Maps.newConcurrentMap();
		
		if (cmd.hasOption("name-grep")) {
			List<String> greps = Lists.newArrayList(Splitter.on(",").split(cmd.getOptionValue("name-grep")));
			List<LogicTreeBranch> filtered = Lists.newArrayList();
			
			branchLoop:
			for (LogicTreeBranch branch : branches) {
				String fname = branch.buildFileName();
				for (String grep : greps) {
					if (!fname.contains(grep))
						continue branchLoop;
				}
				filtered.add(branch);
			}
			
			System.out.println("Filtered branches size: "+filtered.size()+"/"+branches.size());
			branches = filtered;
		}
		
		// sort to ensure order
		Collections.sort(branches);
	}

	@Override
	protected int getNumTasks() {
		return branches.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		List<CalcTask> tasks = Lists.newArrayList();
		
		for (int index : batch)
			tasks.add(new CalcTask(branches.get(index)));
		
		ThreadedTaskComputer comp = new ThreadedTaskComputer(tasks);
		
		comp.computeThreaded(getNumThreads());
		
		for (CalcTask task : tasks) {
			FaultModels fm = task.branch.getValue(FaultModels.class);
			
			if (!charFactorsMap.containsKey(fm)) {
				charFactorsMap.put(fm, new ArrayList<double[]>());
				weightsMap.put(fm, new ArrayList<Double>());
			}
			
			charFactorsMap.get(fm).add(task.result);
			weightsMap.get(fm).add(weightProv.getWeight(task.branch));
		}
	}
	
	private class CalcTask implements Task {
		
		private LogicTreeBranch branch;
		private double[] result;

		public CalcTask(LogicTreeBranch branch) {
			this.branch = branch;
		}

		@Override
		public void compute() {
			debug("Calculating: "+branch.buildFileName());
			// load solution
			InversionFaultSystemSolution sol = cfss.getSolution(branch);
			subSectsMap.putIfAbsent(branch.getValue(FaultModels.class), sol.getRupSet().getFaultSectionDataList());
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			erf.getTimeSpan().setDuration(1d);
			erf.updateForecast();
			
			List<? extends IncrementalMagFreqDist> subSeisMFDs = sol.getSubSeismoOnFaultMFD_List();
			SummedMagFreqDist[] supraSeisMFDs = FaultSysSolutionERF_Calc.calcNucleationMFDForAllSects(erf, 2.55, 8.95, 65);
			
			Preconditions.checkState(subSeisMFDs.size() == supraSeisMFDs.length);
			
			result = new double[supraSeisMFDs.length];
			
			for (int i=0; i<result.length; i++)
				result[i] = 1.0/ETAS_Utils.getScalingFactorToImposeGR_supraRates(supraSeisMFDs[i], subSeisMFDs.get(i), false);
			debug("Done calculating: "+branch.buildFileName());
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		Object[] factorsSendBuff = { charFactorsMap };
		Object[] factorsRecvBuff = null;
		if (rank == 0)
			factorsRecvBuff = new Object[size];
		MPI.COMM_WORLD.Gather(factorsSendBuff, 0, 1, MPI.OBJECT, factorsRecvBuff, 0, 1, MPI.OBJECT, 0);
		
		Object[] weightsSendBuff = { weightsMap };
		Object[] weightsRecvBuff = null;
		if (rank == 0)
			weightsRecvBuff = new Object[size];
		MPI.COMM_WORLD.Gather(weightsSendBuff, 0, 1, MPI.OBJECT, weightsRecvBuff, 0, 1, MPI.OBJECT, 0);
		
		if (rank == 0) {
			for (int i=1; i<size; i++) {
				Map<FaultModels, List<double[]>> oFactorsMap = (Map<FaultModels, List<double[]>>) factorsRecvBuff[i];
				Map<FaultModels, List<Double>> oWeightsMap = (Map<FaultModels, List<Double>>) weightsRecvBuff[i];
				
				for (FaultModels fm : oFactorsMap.keySet()) {
					if (!charFactorsMap.containsKey(fm)) {
						charFactorsMap.put(fm, new ArrayList<double[]>());
						weightsMap.put(fm, new ArrayList<Double>());
					}
					List<double[]> oFactors = oFactorsMap.get(fm);
					List<Double> oWeights = oWeightsMap.get(fm);
					Preconditions.checkState(oFactors.size() == oWeights.size());
					
					charFactorsMap.get(fm).addAll(oFactors);
					weightsMap.get(fm).addAll(oWeights);
				}
			}
			
			// now write results
			for (FaultModels fm : charFactorsMap.keySet()) {
				List<double[]> factors = charFactorsMap.get(fm);
				double[] weights = Doubles.toArray(weightsMap.get(fm));
				
				// look for a common set of branches
				LogicTreeBranch commonBranch = null;
				for (LogicTreeBranch branch : branches) {
					if (branch.getValue(FaultModels.class) != fm)
						continue;
					if (commonBranch == null) {
						commonBranch = (LogicTreeBranch) branch.clone();
					} else {
						for (int i=0; i<branch.size(); i++) {
							if (branch.getValue(i) != commonBranch.getValue(i))
								commonBranch.clearValue(i);
						}
					}
					if (!subSectsMap.containsKey(fm))
						// make sure we have the subsection list
						subSectsMap.put(fm, cfss.getSolution(branch).getRupSet().getFaultSectionDataList());
				}
				String prefix = null;
				for (int i=0; i<commonBranch.size(); i++) {
					LogicTreeBranchNode<?> val = commonBranch.getValue(i);
					if (val == null)
						continue;
					if (prefix == null)
						prefix = "";
					else
						prefix += "_";
					prefix += val.encodeChoiceString();
				}
				
				MatrixIO.doubleArraysListToFile(factors,
						new File(outputDir, prefix+"_char_factors.bin"));
				MatrixIO.doubleArrayToFile(weights,
						new File(outputDir, prefix+"_weights.bin"));
				
				writeCSV(new File(outputDir, prefix+"_char_factors.csv"), factors, weights, subSectsMap.get(fm));
			}
		}
	}
	
	private static void writeCSV(File outputFile, List<double[]> factors, double[] weights,
			List<FaultSectionPrefData> subSects) throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		
//		mean, stdDev, stdDevOfMean, median, and 0.05 and 0.95 fractiles
		double[] fracts = { 0.025, 0.16, 0.84, 0.975 };
		List<String> header = Lists.newArrayList("Index", "Name", "Mean", "Median");
		for (double fract : fracts)
			header.add((float)fract+"");
		header.add("Unweighted Mean");
		header.add("Unweighted Std. Dev.");
		header.add("Unweighted SDOM");
		csv.addLine(header);
		
		for (int s=0; s<subSects.size(); s++) {
			List<String> line = Lists.newArrayList(s+"", subSects.get(s).getName());
			
			double[] sectFracts = new double[weights.length];
			for (int i=0; i<sectFracts.length; i++)
				sectFracts[i] = factors.get(i)[s];
			
			ArbDiscrEmpiricalDistFunc func = new ArbDiscrEmpiricalDistFunc();
			for (int i=0; i<weights.length; i++)
				func.set(sectFracts[i], weights[i]);
			
			line.add(func.getMean()+"");
			line.add(func.getMedian()+"");
			
			for (double fract : fracts)
				line.add(func.getDiscreteFractile(fract)+"");
			
			double unWeightMean = StatUtils.mean(sectFracts);
			double unWeightStdDev = Math.sqrt(StatUtils.variance(sectFracts));
			double unWeightSDOM = unWeightStdDev/Math.sqrt(sectFracts.length);

			line.add(unWeightMean+"");
			line.add(unWeightStdDev+"");
			line.add(unWeightSDOM+"");
			
			csv.addLine(line);
		}
		
		csv.writeToFile(outputFile);
	}
	
	protected static Options createOptions() {
		Options options = MPJTaskCalculator.createOptions();
		
		Option randomSampleOption = new Option("rand", "random-sample", true,
				"If supplied, a random sample of the given size will be used.");
		randomSampleOption.setRequired(false);
		options.addOption(randomSampleOption);
		
		Option nameGrepsOption = new Option("ng", "name-grep", true,
				"If supplied, logic tree branches will be only be included out based on grepping for these comma separated " +
				"strings.");
		nameGrepsOption.setRequired(false);
		options.addOption(nameGrepsOption);
		
		return options;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_ETAS_Simulator.class);
			
			args = cmd.getArgs();
			
			if (args.length != 2) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(MPJ_ETAS_CharFactorCalc.class)
						+" [options] <compound-file> <output-dir>");
				abortAndExit(2);
			}
			
			File compoundFile = new File(args[0]);
			Preconditions.checkArgument(compoundFile.exists(),
					"Compound solution file doesn't exit: "+compoundFile.getAbsolutePath());
			
			File outputDir = new File(args[1]);
			Preconditions.checkArgument(outputDir.exists() || outputDir.mkdir(),
					"output directory doesn't exist: "+outputDir.getAbsolutePath());
			
			MPJ_ETAS_CharFactorCalc driver = new MPJ_ETAS_CharFactorCalc(cmd, compoundFile, outputDir);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
