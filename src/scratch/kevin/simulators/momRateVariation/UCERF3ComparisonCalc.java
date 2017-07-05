package scratch.kevin.simulators.momRateVariation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.EventRecord;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;

public class UCERF3ComparisonCalc {
	
	private static void runUCERF3Sim(FaultSystemSolution sol, int duration, File outputDir,
			MagDependentAperiodicityOptions cov) throws IOException, DocumentException {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
		erf.setParameter(MagDependentAperiodicityParam.NAME, cov);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		
		erf.getTimeSpan().setStartTime(2014);
		System.out.println("startYear: "+erf.getTimeSpan().getStartTimeYear());
		erf.eraseDatesOfLastEventAfterStartTime();
		erf.getTimeSpan().setDuration((double)duration);

		erf.updateForecast();
		
		ProbabilityModelsCalc calc = new ProbabilityModelsCalc(erf);
		
		calc.testER_NextXyrSimulation(outputDir, null, 1, false, null);
	}
	
	private static class CalcThread extends Thread {
		private FaultSystemSolution sol;
		private int duration;
		private File outputDir;
		private MagDependentAperiodicityOptions cov;
		
		public CalcThread(FaultSystemSolution sol, int duration, File outputDir,
				MagDependentAperiodicityOptions cov) {
			this.sol = sol;
			this.duration = duration;
			this.outputDir = outputDir;
			this.cov = cov;
		}

		@Override
		public void run() {
			try {
				runUCERF3Sim(sol, duration, outputDir, cov);
			} catch (Exception e) {
				e.printStackTrace();
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}
	
	private static void runUCERF3SimThreaded(File fssFile, int numThreads, final int duration, File outputDir,
			MagDependentAperiodicityOptions cov) throws IOException, DocumentException {
		List<Thread> threads = Lists.newArrayList();
		
		for (int i=0; i<numThreads; i++) {
			// need new one for each thread
			FaultSystemSolution sol = FaultSystemIO.loadSol(fssFile);
			
			File subDir = new File(outputDir, duration+"yr_run"+i);
			Preconditions.checkState(subDir.exists() || subDir.mkdir());;
			
			threads.add(new CalcThread(sol, duration, subDir, cov));
		}
		
		for (Thread thread : threads)
			thread.start();
		
		for (Thread thread : threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir, fssFile;
		int numThreads, duration;
		MagDependentAperiodicityOptions cov = MagDependentAperiodicityOptions.MID_VALUES;
		if (args.length == 0) {
			outputDir = new File("/home/kevin/Simulators/time_series/ucerf3_compare");
			fssFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
			
			numThreads = 8;
			duration = 10000;
		} else {
			if (args.length < 4 || args.length > 5) {
				System.err.println("Usage: "+ClassUtils.getClassNameWithoutPackage(UCERF3ComparisonCalc.class)
						+" <outputDir> <fssFile> <numThreads> <duration> [<cov>]");
				System.exit(2);
			}
			outputDir = new File(args[0]);
			fssFile = new File(args[1]);
			numThreads = Integer.parseInt(args[2]);
			duration = Integer.parseInt(args[3]);
			if (args.length == 5)
				cov = MagDependentAperiodicityOptions.valueOf(args[4]);
		}
		
//		int numThreads = 1;
//		int duration = 100;
		
		runUCERF3SimThreaded(fssFile, numThreads, duration, outputDir, cov);
	}

}
