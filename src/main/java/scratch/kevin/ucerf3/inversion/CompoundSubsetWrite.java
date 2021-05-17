package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.zip.ZipException;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.analysis.CompoundFSSPlots;
import scratch.UCERF3.analysis.CompoundFSSPlots.BranchAvgFSSBuilder;
import scratch.UCERF3.analysis.CompoundFSSPlots.ParentSectMFDsPlot;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.UCERF3p2BranchWeightProvider;
import scratch.UCERF3.utils.DeformationModelFetcher;

public class CompoundSubsetWrite {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// YOU WILL NEED TO REVERT THE COULOMB FILES FOR THIS TO WORK!!!
		
		File origCompoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/" +
				"scratch/InversionSolutions/2013_01_14-stampede_3p2_production_runs_combined_COMPOUND_SOL.zip");
		File newCompoundFile = new File("/tmp/2013_01_14-stampede_3p2_production_runs_fm3p1_dm_scale_subset_COMPOUND_SOL.zip");
		File refCompoundFile = new File("/tmp/2013_04_28-refactored-test-runs_COMPOUND_SOL.zip");
		
		final CompoundFaultSystemSolution origCompound = CompoundFaultSystemSolution.fromZipFile(origCompoundFile);
		CompoundFaultSystemSolution refCompound = CompoundFaultSystemSolution.fromZipFile(refCompoundFile);
		
		System.out.println("Compound has "+origCompound.getBranches().size()+" branches");
		
		// we want ref branch
		LogicTreeBranch origRefBranch = LogicTreeBranch.fromValues(false, FaultModels.FM3_1, null,
				null, SlipAlongRuptureModels.UNIFORM, InversionModels.CHAR_CONSTRAINED,
				TotalMag5Rate.RATE_8p7, MaxMagOffFault.MAG_7p6, MomentRateFixes.NONE, SpatialSeisPDF.UCERF3);
		
		UCERF3p2BranchWeightProvider ucerf3p2WeightProvider = new UCERF3p2BranchWeightProvider();
		
		final List<LogicTreeBranch> subBranches = Lists.newArrayList();
		for (LogicTreeBranch branch : origCompound.getBranches()) {
			if (origRefBranch.matchesNonNulls(branch)) {
//			if (branch.matchesNonNulls(origRefBranch)) {
				Preconditions.checkState(ucerf3p2WeightProvider.getWeight(branch) > 0);
				System.out.println(branch);
				subBranches.add(branch);
			}
		}
		
		System.out.println("Found "+subBranches.size()+" sub branches");
		
		FaultSystemSolutionFetcher fetch = new FaultSystemSolutionFetcher() {
			
			@Override
			public Collection<LogicTreeBranch> getBranches() {
				return subBranches;
			}
			
			@Override
			protected InversionFaultSystemSolution fetchSolution(LogicTreeBranch branch) {
				return origCompound.getSolution(branch);
			}
		};
		
		DeformationModelFetcher.IMPERIAL_DDW_HACK = true;
		
//		CompoundFaultSystemSolution.toZipFile(newCompoundFile, fetch);
//		
//		MeanFSSBuilder builder = new MeanFSSBuilder(weightProv);
//		List<CompoundFSSPlots> plots = Lists.newArrayList();
//		plots.add(builder);
		
//		CompoundFSSPlots.batchPlot(plots, fetch, 1);
//		String prefix = "2013_01_14-stampede_3p2_production_runs_fm3p1_dm_scale_subset";
//		CompoundFSSPlots.batchWritePlots(plots, new File("/tmp"), prefix, true);
		
		// this is for generating comparison parent section MFD plots
		// first create for the new solution
		System.out.println("Calculating new MFD plots");
		APrioriBranchWeightProvider newWeightProvider = new APrioriBranchWeightProvider();
		ParentSectMFDsPlot refMFDs = new ParentSectMFDsPlot(newWeightProvider);
		List<CompoundFSSPlots> plots = Lists.newArrayList();
		plots.add(refMFDs);
		CompoundFSSPlots.batchPlot(plots, refCompound, 1);
		
		// now do for the old solution subset
		System.out.println("Calculating old MFD plots");
		ParentSectMFDsPlot oldMFDs = new ParentSectMFDsPlot(ucerf3p2WeightProvider);
		plots = Lists.newArrayList();
		plots.add(oldMFDs);
		CompoundFSSPlots.batchPlot(plots, fetch, 1);
		
		// now this is where it gets a little messy, add the mean curve in from the orig to the ref as a fractile
		System.out.println("Merging");
		refMFDs.addMeanFromExternalAsFractile(oldMFDs);
		
		System.out.println("Writing");
		File mfdOutDir = new File("/tmp/mfd_compare");
		CompoundFSSPlots.writeParentSectionMFDPlots(refMFDs, mfdOutDir);
	}

}
