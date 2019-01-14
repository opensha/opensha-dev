package scratch.kevin.ucerf3.eal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMM_Epistemic;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMMs;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_LogicTreeBranch;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_ProbModels;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_Vs30Model;

public class MPJ_UCERF3_EAL_Combiner extends MPJTaskCalculator {
	
	private FaultSystemSolution trueMeanSol;
	private Map<LogicTreeBranch, List<Integer>> mappings;
	private CompoundFaultSystemSolution cfss;
	
	private double erfProbsDuration;
	private Map<U3_EAL_ProbModels, ZipFile> probsZipFiles;
	
	private Map<U3_EAL_Vs30Model, File> vs30Dirs;
	
	private File resultsDir;
	private File resultsFile;
	
	private List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>();
	
	private Map<File, double[][]> rupLossesMap;
	private Map<File, DiscretizedFunc[]> griddedLossesMap;
	
	private ExecutorService exec;
	
	private boolean consolidateOnly;
	private CSVFile<String> resultsCSV;

	public MPJ_UCERF3_EAL_Combiner(CommandLine cmd, File outputDir) throws IOException, DocumentException {
		super(cmd);
		
		File trueMeanSolFile = new File(cmd.getOptionValue("true-mean-sol"));
		File compoundSolFile = new File(cmd.getOptionValue("compound-sol"));
		
		consolidateOnly = cmd.hasOption("consolidate-only");
		
		if (rank == 0)
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
					"Output dir doesn't exist and could not be created: %s", outputDir.getAbsolutePath());
		resultsDir = new File(outputDir, "results");
		if (consolidateOnly && rank == 0) {
			debug("Consolidating only");
			Preconditions.checkState(resultsDir.exists(), "Consolidate only but no results dir");
			size = 0;
			for (File file : resultsDir.listFiles())
				if (file.getName().startsWith("results_") && file.getName().endsWith(".csv"))
					size++;
			debug("new size="+size+" for consolidation");
		}
		
		if (rank == 0)
			debug("Loading true mean solution from: "+trueMeanSolFile.getAbsolutePath());
		trueMeanSol = FaultSystemIO.loadSol(trueMeanSolFile);
		// now load in the mappings
		if (rank == 0)
			debug("Loading true mean branch mappings");
		mappings = TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile);
		if (rank == 0)
			debug("Loading CFSS: "+compoundSolFile.getAbsolutePath());
		cfss = CompoundFaultSystemSolution.fromZipFile(compoundSolFile);
		
		if (rank == 0)
			Preconditions.checkState(resultsDir.exists() || resultsDir.mkdir(),
					"Results dir doesn't exist and could not be created: %s", resultsDir.getAbsolutePath());
		
		erfProbsDuration = Double.parseDouble(cmd.getOptionValue("erf-probs-duration"));
		
		File probsZipDir = new File(cmd.getOptionValue("erf-probs-dir"));
		Preconditions.checkState(probsZipDir.exists(), "Probs zip file doesn't exist: %s", probsZipDir.getAbsolutePath());
		
		probsZipFiles = new HashMap<>();
		for (U3_EAL_ProbModels probModel : U3_EAL_ProbModels.values()) {
			File probFile = new File(probsZipDir, "probs_"+(float)+erfProbsDuration+"yr_"+probModel.getShortName()+".zip");
			if (!probFile.exists() && (float)Math.round(erfProbsDuration) == (float)erfProbsDuration)
				probFile = new File(probsZipDir, "probs_"+(int)+erfProbsDuration+"yr_"+probModel.getShortName()+".zip");
			if (probFile.exists())
				probsZipFiles.put(probModel, new ZipFile(probFile));
		}
		Preconditions.checkState(!probsZipFiles.isEmpty(), "No prob zip files with duration=%s found in %s",
				(float)erfProbsDuration, probsZipDir.getAbsolutePath());
		
		vs30Dirs = new HashMap<>();
		if (cmd.hasOption("wills-dir")) {
			File willsDir = new File(cmd.getOptionValue("wills-dir"));
			Preconditions.checkState(willsDir.exists(), "Wills dir doesn't exist: %s", willsDir.getAbsolutePath());
			vs30Dirs.put(U3_EAL_Vs30Model.WILLS_2015, willsDir);
		}
		if (cmd.hasOption("wald-dir")) {
			File waldDir = new File(cmd.getOptionValue("wald-dir"));
			Preconditions.checkState(waldDir.exists(), "Wald dir doesn't exist: %s", waldDir.getAbsolutePath());
			vs30Dirs.put(U3_EAL_Vs30Model.WALD_ALLEN, waldDir);
		}
		Preconditions.checkArgument(!vs30Dirs.isEmpty(), "No Vs30 model directories supplied!");
		// TODO tracts
		
		if (rank == 0)
			debug("Building EAL branches");
		branches = new ArrayList<>();
		HashSet<U3_EAL_ProbModels> probModels = new HashSet<>();
		HashSet<U3_EAL_GMMs> gmms = new HashSet<>();
		HashSet<U3_EAL_GMM_Epistemic> gmmEpis = new HashSet<>();
		HashSet<U3_EAL_Vs30Model> vs30s = new HashSet<>();
		int numTI = mappings.keySet().size();
		for (U3_EAL_ProbModels probModel : probsZipFiles.keySet()) {
			for (U3_EAL_GMMs gmm : U3_EAL_GMMs.values()) {
				for (U3_EAL_GMM_Epistemic gmmEpi : U3_EAL_GMM_Epistemic.values()) {
					for (U3_EAL_Vs30Model vs30 : vs30Dirs.keySet()) {
						File vs30Dir = vs30Dirs.get(vs30);
						String prefix = gmm.getShortName();
						if (gmmEpi != U3_EAL_GMM_Epistemic.NONE)
							prefix += "_"+gmmEpi.getShortName();
						
						File fssFile = new File(vs30Dir, prefix+"_fss_index.bin");
						if (!fssFile.exists())
							fssFile = new File(vs30Dir, prefix+"_fss_index.bin.gz");
						if (!fssFile.exists())
							continue;
						
						File griddedFile = new File(vs30Dir, prefix+"_fss_gridded.bin");
						if (!griddedFile.exists())
							griddedFile = new File(vs30Dir, prefix+"_fss_gridded.bin.gz");
						if (!fssFile.exists())
							griddedFile = null;
						
						File tractDir = new File(vs30Dir, prefix+"_tract_results");
						if (!tractDir.exists())
							tractDir = null;
						
						probModels.add(probModel);
						gmms.add(gmm);
						gmmEpis.add(gmmEpi);
						vs30s.add(vs30);
						for (LogicTreeBranch tiBranch : mappings.keySet())
							branches.add(new U3_EAL_LogicTreeBranch(tiBranch, probModel, gmm, gmmEpi, vs30, fssFile, griddedFile, tractDir));
					}
				}
			}
		}
		if (rank == 0) {
			debug(numTI+" TI branches");
			debug("Prob Models: "+Joiner.on(", ").join(probModels));
			debug("GMMs: "+Joiner.on(", ").join(gmms));
			debug("GMM Epis: "+Joiner.on(", ").join(gmmEpis));
			debug("Vs30s: "+Joiner.on(", ").join(vs30s));
			int calcNum = numTI*probModels.size()*gmms.size()*gmmEpis.size()*vs30s.size();
			debug("Calculated number if fully specified: "+calcNum+" (fully specified ? "+(calcNum == branches.size())+")");
		}
		Collections.sort(branches);
		debug("Built "+branches.size()+" branches");
		Preconditions.checkState(!branches.isEmpty(), "No branches found!");
		
		rupLossesMap = new HashMap<>();
		griddedLossesMap = new HashMap<>();
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		
		resultsCSV = new CSVFile<>(true);
		U3_EAL_LogicTreeBranch branch0 = branches.get(0);
		List<String> header = Lists.newArrayList("Index", "Branch Weight", "Total EAL", "Fault EAL", "Gridded EAL");
		for (int i=0; i<branch0.size(); i++)
			header.add(branch0.getValue(i).getBranchLevelName());
		resultsCSV.addLine(header);
		resultsFile = new File(resultsDir, "results_"+rank+".csv");
	}

	@Override
	protected int getNumTasks() {
		if (consolidateOnly)
			return 1;
		return branches.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		if (consolidateOnly)
			return;
		List<Future<CalcTask>> futures = new ArrayList<>();
		
		for (int index : batch)
			futures.add(exec.submit(new CalcTask(index)));
		
		debug("Waiting on "+futures.size()+" futures");
		for (Future<CalcTask> future : futures) {
			CalcTask task = null;
			try {
				task = future.get();
			} catch (Exception e) {
				abortAndExit(e);
			}
			List<String> line = new ArrayList<>();
			line.add(task.index+"");
			line.add(task.branch.getAprioriBranchWt()+"");
			line.add(task.totalEAL+"");
			line.add(task.faultEAL+"");
			line.add(task.griddedEAL+"");
			for (LogicTreeBranchNode<?> node : task.branch)
				line.add(node.getShortName());
			resultsCSV.addLine(line);
		}
		debug("finished batch and flushing CSV");
		resultsCSV.writeToFile(resultsFile);
	}
	
	private synchronized double[][] getLossesFSS(File fssLossFile) throws IOException {
		double[][] expectedLosses = rupLossesMap.get(fssLossFile);
		if (expectedLosses != null)
			return expectedLosses;
		expectedLosses = MPJ_CondLossCalc.loadResults(fssLossFile);
		rupLossesMap.put(fssLossFile, expectedLosses);
		return expectedLosses;
	}
	
	private synchronized DiscretizedFunc[] getLossesGridded(File griddedLossFile) throws IOException {
		if (griddedLossFile == null)
			return null;
		DiscretizedFunc[] griddedFuncs = griddedLossesMap.get(griddedLossFile);
		if (griddedFuncs != null)
			return griddedFuncs;
		griddedFuncs = MPJ_CondLossCalc.loadGridSourcesFile(griddedLossFile,
				trueMeanSol.getGridSourceProvider().getGriddedRegion());
		griddedLossesMap.put(griddedLossFile, griddedFuncs);
		return griddedFuncs;
	}
	
	private class CalcTask implements Callable<CalcTask> {
		
		private int index;
		private U3_EAL_LogicTreeBranch branch;
		
		private double faultEAL;
		private double griddedEAL;
		private double totalEAL;

		public CalcTask(int index) {
			this.index = index;
			this.branch = branches.get(index);
		}

		@Override
		public CalcTask call() throws Exception {
			Map<LogicTreeBranch, List<Integer>> taskMappings = new HashMap<>();
			taskMappings.put(branch.getTIBranch(), mappings.get(branch.getTIBranch()));
			
			double[][] fssLosses = getLossesFSS(branch.getFSSIndexedBinFile());
			DiscretizedFunc[] griddedLosses = getLossesGridded(branch.getGriddedBinFile());
			ZipFile erfProbsZipFile = probsZipFiles.get(branch.getValue(U3_EAL_ProbModels.class));
			
			UCERF3_EAL_Combiner calc = new UCERF3_EAL_Combiner(cfss, taskMappings, trueMeanSol, fssLosses, griddedLosses,
					erfProbsZipFile, erfProbsDuration);
			
			faultEAL = calc.getFaultEALs()[0];
			griddedEAL = calc.getGriddedEALs()[0];
			totalEAL = calc.getTotalEALs()[0];
			
			return this;
		}
		
	}
	
	private class BranchNodeVals {
		private double weightTotal;
		private double meanLoss;
		private double meanFault;
		private double meanGridded;
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		exec.shutdown();
		rupLossesMap.clear();
		griddedLossesMap.clear();
		if (rank == 0) {
			debug("Consolidating CSVs");
			List<List<String>> allLines = new ArrayList<>();
			for (int i=0; i<branches.size(); i++)
				allLines.add(null);
			int loaded = 0;
			List<Map<LogicTreeBranchNode<?>, BranchNodeVals>> nodeVals = new ArrayList<>();
			for (int i=0; i<branches.get(0).size(); i++)
				nodeVals.add(new HashMap<>());
			double totalWeight = 0d;
			double totalMean = 0d;
			double faultMean = 0d;
			double griddedMean = 0d;
			for (int r=0; r<size; r++) {
				CSVFile<String> csv;
				if (r == 0 && !consolidateOnly) {
					csv = this.resultsCSV;
				} else {
					File resultsFile = new File(resultsDir, "results_"+r+".csv");
					if (!resultsFile.exists()) {
						debug("No results for rank "+r);
						continue;
					}
					csv = CSVFile.readFile(resultsFile, true);
				}
				for (int row=1; row<csv.getNumRows(); row++) {
					loaded++;
					int index = Integer.parseInt(csv.get(row, 0));
					Preconditions.checkState(allLines.get(index) == null, "Duplicate found for index "+index+" in rank "+r);
					allLines.set(index, csv.getLine(row));
					U3_EAL_LogicTreeBranch branch = branches.get(index);
					for (int i=0; i<branch.size(); i++) {
//						getBranchLevelName()
						String choice = branch.getValue(i).getShortName();
						String testChoice = csv.get(row, i+5);
						if (!choice.equals(testChoice)) {
							System.err.println("Branch mismatch for rank "+r+", index "+index);
							System.err.println("\tOriginal Branch: "+branch);
							List<String> line = csv.getLine(row);
							System.err.println("\tCSV Branch: "+Joiner.on(", ").join(line.subList(5, line.size())));
							System.err.flush();
							throw new IllegalStateException("Branch mismatch for rank "+r+", index "+index);
						}
					}
					double weight = branch.getAprioriBranchWt();
					double totEAL = Double.parseDouble(csv.get(row, 2));
					double faultEAL = Double.parseDouble(csv.get(row, 3));
					double griddedEAL = Double.parseDouble(csv.get(row, 4));
					totalWeight += weight;
					totalMean += totEAL*weight;
					faultMean += faultEAL*weight;
					griddedMean += griddedEAL*weight;
					for (int i=0; i<branch.size(); i++) {
						BranchNodeVals vals = nodeVals.get(i).get(branch.getValue(i));
						LogicTreeBranchNode<?> node = branch.getValue(i);
						if (vals == null) {
							vals = new BranchNodeVals();
							nodeVals.get(i).put(node, vals);
						}
						vals.weightTotal += weight;
						vals.meanLoss += totEAL*weight;
						vals.meanFault += faultEAL*weight;
						vals.meanGridded += griddedEAL*weight;
					}
				}
			}
			debug("Loaded "+loaded+" branches");
			debug("Total weight: "+totalWeight);
			Preconditions.checkState(loaded == branches.size(), "Did not load all branches. Expected %s, loaded %s", branches.size(), loaded);
			debug("Writing master CSV");
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine(resultsCSV.getLine(0)); // header
			csv.addAll(allLines);
			csv.writeToFile(new File("all_branch_results.csv"));
			debug("Writing branch levels CSV");
			csv = new CSVFile<>(true);
			csv.addLine("Branch Level", "Branch Choice", "Total Weight", "Weighted Total Mean EAL",
					"Weighted Fault Mean EAL", "Weighted Gridded Mean EAL");
			Comparator<LogicTreeBranchNode<?>> nodeComparator = new Comparator<LogicTreeBranchNode<?>>() {

				@Override
				public int compare(LogicTreeBranchNode<?> o1, LogicTreeBranchNode<?> o2) {
					return o1.getShortName().compareTo(o2.getShortName());
				}

			};
			
			for (int i=0; i<nodeVals.size(); i++) {
				Map<LogicTreeBranchNode<?>, BranchNodeVals> valsMap = nodeVals.get(i);
				List<LogicTreeBranchNode<?>> choices = new ArrayList<>(valsMap.keySet());
				choices.sort(nodeComparator);
				for (LogicTreeBranchNode<?> choice : choices) {
					List<String> line = new ArrayList<>();
					line.add(choice.getBranchLevelName());
					line.add(choice.getShortName());
					BranchNodeVals vals = valsMap.get(choice);
					// normalize
					vals.meanLoss /= vals.weightTotal;
					vals.meanFault /= vals.weightTotal;
					vals.meanGridded /= vals.weightTotal;
					line.add((float)(vals.weightTotal/totalWeight)+"");
					line.add(vals.meanLoss+"");
					line.add(vals.meanFault+"");
					line.add(vals.meanGridded+"");
					csv.addLine(line);
				}
			}
			// normalize totals
			totalMean /= totalWeight;
			faultMean /= totalWeight;
			griddedMean /= totalWeight;
			csv.addLine("COMPLETE MODEL", "MEAN", "1.0", totalMean+"", faultMean+"", griddedMean+"");
			csv.writeToFile(new File("branch_level_summary.csv"));
		}
	}
	
	public static Options createOptions() {
		Options options = MPJTaskCalculator.createOptions();
		
		Option willsDir = new Option("wills", "wills-dir", true, "Directory containing Wills 2015 results");
		willsDir.setRequired(false);
		options.addOption(willsDir);
		
		Option waldDir = new Option("wald", "wald-dir", true, "Directory containing Wald & Allen results");
		waldDir.setRequired(false);
		options.addOption(waldDir);
		
		Option erfProbsDir = new Option("probs", "erf-probs-dir", true, "ERF probabiltiies directory");
		erfProbsDir.setRequired(true);
		options.addOption(erfProbsDir);
		
		Option erfProbsDuration = new Option("dur", "erf-probs-duration", true, "ERF probabiltiies duration");
		erfProbsDuration.setRequired(true);
		options.addOption(erfProbsDuration);
		
		Option trueMeanSol = new Option("tms", "true-mean-sol", true, "True mean solution file (with mappings)");
		trueMeanSol.setRequired(true);
		options.addOption(trueMeanSol);
		
		Option compoundSol = new Option("cfss", "compound-sol", true, "Compound FSS File");
		compoundSol.setRequired(true);
		options.addOption(compoundSol);
		
		Option consolidateOnly = new Option("co", "consolidate-only", false, "Flag to consolidate only");
		consolidateOnly.setRequired(false);
		options.addOption(consolidateOnly);
		
		return options;
	}

	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_UCERF3_EAL_Combiner.class);
			
			args = cmd.getArgs();
			
			if (args.length != 1) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(MPJ_UCERF3_EAL_Combiner.class)
						+" <output-dir>");
				abortAndExit(2);
			}
			
			File outputDir = new File(args[0]);
			
			MPJ_UCERF3_EAL_Combiner driver = new MPJ_UCERF3_EAL_Combiner(cmd, outputDir);
			
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
