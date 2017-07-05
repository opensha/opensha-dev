package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import mpi.MPI;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.dists.RandomDistType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class MPJSynchLagRand extends MPJTaskCalculator {
	
	private List<double[][][]> gBarsList;
	
	private int numTrials;
	private List<? extends SimulatorEvent> events;
	private List<RuptureIdentifier> rupIdens;
	
	private RandomDistType dist = RandomDistType.ACTUAL;
	private double distSpacing = 10d;
	
	private int nDims;
	
	private static int[] lags = SynchParamCalculator.rangeInclusive(-30, 30);
	
	private File outputDir;
	private String setName;
	
	private EmpiricalMarkovChain origChain;

	public MPJSynchLagRand(CommandLine cmd, File outputDir) throws IOException {
		super(cmd);
		
		this.outputDir = outputDir;
		
		Preconditions.checkArgument(cmd.hasOption("trials"));
		numTrials = Integer.parseInt(cmd.getOptionValue("trials"));
		
		setName = cmd.getOptionValue("set");
		
		gBarsList = Lists.newArrayList();
		
		if (setName.equals("so_cal"))
			rupIdens = SynchIdens.getStandardSoCal();
		else if (setName.equals("nor_cal"))
			rupIdens = SynchIdens.getStandardNorCal();
		else
			throw new IllegalArgumentException("Unknown set: "+setName);
		
		events = new SimAnalysisCatLoader(true, rupIdens, false).getEvents();
		
		nDims = rupIdens.size();
	}

	@Override
	protected int getNumTasks() {
		return numTrials;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		checkBuildOrigChain();
		
		for (int item : batch) {
			double[][][] gBars = new double[nDims][nDims][lags.length];
			EmpiricalMarkovChain chain = SynchParamCalculator.createRandomizedChain(events, rupIdens, dist, distSpacing);
			
			List<SynchCalc> tasks = Lists.newArrayList();
			
			for (int m=0; m<nDims; m++)
				for (int n=m; n<nDims; n++)
					for (int l=0; l<lags.length; l++)
						tasks.add(new SynchCalc(chain, m, n, l, gBars));
			
			ThreadedTaskComputer comp = new ThreadedTaskComputer(tasks);
			comp.computeThreaded(getNumThreads());
			
			gBarsList.add(gBars);
		}
	}
	
	private void checkBuildOrigChain() {
		if (rank == 0 && origChain == null) {
			origChain = MarkovChainBuilder.build(distSpacing, events, rupIdens);
		}
	}
	
	static class SynchCalc implements Task {
		
		private EmpiricalMarkovChain chain;
		private int m, n, l;
		private double[][][] gBars;

		public SynchCalc(EmpiricalMarkovChain chain, int m, int n, int l, double[][][] gBars) {
			super();
			this.chain = chain;
			this.m = m;
			this.n = n;
			this.l = l;
			this.gBars = gBars;
		}

		@Override
		public void compute() {
			int lag = lags[l];
			double gBar = SynchParamCalculator.calcGBar(chain, m, n, lag);
			
			gBars[m][n][l] = gBar;
			gBars[n][m][l] = gBar;
		}
		
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	protected void doFinalAssembly() throws Exception {
		List[] sendbuf = new List[] {gBarsList};
		List[] recvbuf = null;
		
		if (rank == 0)
			recvbuf = new List[size];
		
		MPI.COMM_WORLD.Gather(sendbuf, 0, 1, MPI.OBJECT, recvbuf, 0, 1, MPI.OBJECT, 0);
		
		if (rank == 0) {
			for (int i=1; i<size; i++) {
				List<double[][][]> o = recvbuf[i];
				gBarsList.addAll(o);
			}
			
			Preconditions.checkState(numTrials == gBarsList.size());
			
			double[][][][] gBars = new double[nDims][nDims][numTrials][lags.length];
			
			for (int t=0; t<numTrials; t++) {
				double[][][] trialGBars = gBarsList.get(t);
				for (int m=0; m<nDims; m++)
					for (int n=0; n<nDims; n++)
						for (int l=0; l<lags.length; l++)
							gBars[m][n][t][l] = trialGBars[m][n][l];
			}
			
			File writeDir = new File(outputDir, SynchParamCalculator.getDirName());
			if (!writeDir.exists())
				writeDir.mkdir();
			
			writeDir = new File(writeDir, setName);
			if (!writeDir.exists())
				writeDir.mkdir();
			
			checkBuildOrigChain();
			SynchParamCalculator.doWriteSynchStdDevParams(writeDir, rupIdens, origChain, lags, numTrials, nDims, gBars);
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option trialsOption = new Option("trials", "num-trials", true,
				"Number of random trials");
		trialsOption.setRequired(true);
		ops.addOption(trialsOption);
		
		Option setOption = new Option("s", "set", true,
				"Set of faults to use (nor_cal or so_cal)");
		setOption.setRequired(true);
		ops.addOption(setOption);
		
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJSynchLagRand.class);
			
			args = cmd.getArgs();
			
			if (args.length != 1) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(MPJSynchLagRand.class)
						+" [options] <output-dir>");
				abortAndExit(2);
			}
			
			File outputDir = new File(args[0]);
			if (!outputDir.exists() && MPI.COMM_WORLD.Rank() == 0)
				Preconditions.checkState(outputDir.mkdir());
			
			MPJSynchLagRand calc = new MPJSynchLagRand(cmd, outputDir);
			calc.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
