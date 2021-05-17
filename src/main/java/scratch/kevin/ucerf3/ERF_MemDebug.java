package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.utils.FaultSystemIO;

public class ERF_MemDebug {
	
	public static void main(String[] args) throws IOException, DocumentException, InterruptedException {
		Thread.sleep(4000);
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/"
				+ "dev/scratch/UCERF3/data/scratch/InversionSolutions/"
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL.zip"));
		System.out.println("Sol only");
		printMemoryDebug();
		boolean timeDep = false;
		FaultSystemSolutionERF erf;
		erf = new FaultSystemSolutionERF(sol);
//		erf = new MeanUCERF3(sol);
//		((MeanUCERF3)erf).setCachingEnabled(false);
//		((MeanUCERF3)erf).setMeanParams(0.1d, true, 0d, MeanUCERF3.RAKE_BASIS_NONE);
		if (timeDep)
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
		erf.getTimeSpan().setDuration(50);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
//		erf.setUseFSSDateOfLastEvents(true);
		System.out.println("ERF Instantiation only");
		printMemoryDebug();
		erf.updateForecast();
		System.out.println("Done updating");
		System.out.println("ERF Updated");
		printMemoryDebug();
		while (true) {
			for (ProbEqkSource source : erf) {
				// do nothing
				source.getRupture(0).getMag();
			}
		}
	}
	
	private static void printMemoryDebug() {
		System.gc();
		try {
			Thread.sleep(100);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Runtime rt = Runtime.getRuntime();
		long totalMB = rt.totalMemory() / 1024 / 1024;
		long freeMB = rt.freeMemory() / 1024 / 1024;
		long usedMB = totalMB - freeMB;
		System.out.println("mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB);
	}

}
