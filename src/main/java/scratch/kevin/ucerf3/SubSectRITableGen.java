package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.zip.ZipException;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.analysis.CompoundFSSPlots;
import scratch.UCERF3.analysis.CompoundFSSPlots.SubSectRITable;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;

public class SubSectRITableGen {

	public static void main(String[] args) throws Exception {
		File dir = new File("/tmp");
		File compoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
		FaultModels fm = FaultModels.FM3_1;
//		DeformationModels dm = DeformationModels.GEOLOGIC;
		
		BranchWeightProvider weightProvider = new APrioriBranchWeightProvider();
		
//		FaultSystemSolutionFetcher fetch = FaultSystemSolutionFetcher.getSubset(cfss, fm, dm);
		FaultSystemSolutionFetcher fetch = FaultSystemSolutionFetcher.getSubset(cfss, fm);
		// for random sample
//		fetch = FaultSystemSolutionFetcher.getRandomSample(fetch, 8, FaultModels.FM3_1);
		System.out.println("Subset has "+fetch.getBranches().size()+" branches");
		
		SubSectRITable calc = new SubSectRITable(weightProvider);
		
		List<CompoundFSSPlots> plots = Lists.newArrayList();
		plots.add(calc);
		
		int threads = Runtime.getRuntime().availableProcessors();
		System.out.println("Calculating with "+threads+" threads");
		CompoundFSSPlots.batchPlot(plots, fetch, threads);
		System.out.println("Compiling/writing tables");
//		CompoundFSSPlots.batchWritePlots(plots, dir, fm.getShortName()+"_"+dm.getShortName(), true);
		CompoundFSSPlots.batchWritePlots(plots, dir, fm.getShortName(), true);
	}

}
