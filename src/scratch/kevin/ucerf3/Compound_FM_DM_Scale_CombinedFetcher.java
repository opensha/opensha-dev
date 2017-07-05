package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;

public class Compound_FM_DM_Scale_CombinedFetcher extends FaultSystemSolutionFetcher {
	
	private Map<LogicTreeBranch, List<LogicTreeBranch>> branchesMap;
	private FaultSystemSolutionFetcher fetch;
	private BranchWeightProvider weightProv;
	
	public Compound_FM_DM_Scale_CombinedFetcher(FaultSystemSolutionFetcher fetch, BranchWeightProvider weightProv) {
		branchesMap = Maps.newHashMap();
		for (LogicTreeBranch branch : fetch.getBranches()) {
			LogicTreeBranch meanBranch = getMeanBranch(branch);
			List<LogicTreeBranch> subBranches = branchesMap.get(meanBranch);
			if (subBranches == null) {
				subBranches = Lists.newArrayList();
				branchesMap.put(meanBranch, subBranches);
			}
			subBranches.add(branch);
		}
		
		System.out.println("Found "+branchesMap.keySet().size()+" sub branches");
		
		this.fetch = fetch;
		this.weightProv = weightProv;
	}
	
	private static LogicTreeBranch getMeanBranch(LogicTreeBranch branch) {
		FaultModels fm = branch.getValue(FaultModels.class);
		DeformationModels dm = branch.getValue(DeformationModels.class);
		ScalingRelationships scale = branch.getValue(ScalingRelationships.class);
		
		return LogicTreeBranch.fromValues(fm, dm, scale, SlipAlongRuptureModels.MEAN_UCERF3,
				InversionModels.CHAR_CONSTRAINED, TotalMag5Rate.RATE_6p5,
				MaxMagOffFault.MAG_7p6, MomentRateFixes.NONE, SpatialSeisPDF.UCERF2);
	}

	@Override
	public Collection<LogicTreeBranch> getBranches() {
		return branchesMap.keySet();
	}

	@Override
	protected InversionFaultSystemSolution fetchSolution(LogicTreeBranch branch) {
		List<LogicTreeBranch> subBranches = branchesMap.get(branch);
		InversionFaultSystemSolution refSol = fetch.getSolution(subBranches.get(0));
		InversionFaultSystemRupSet rupSet = refSol.getRupSet();
		// change the branch
		rupSet = new InversionFaultSystemRupSet(rupSet, branch, rupSet.getLaughTestFilter(),
				rupSet.getAveSlipForAllRups(), rupSet.getCloseSectionsListList(),
				rupSet.getRupturesForClusters(), rupSet.getSectionsForClusters());
		double[] rates = new double[rupSet.getNumRuptures()];
		double totWeight = 0d;
		for (LogicTreeBranch subBranch : subBranches) {
			double weight = weightProv.getWeight(subBranch);
			double[] subRates = fetch.getRates(subBranch);
			for (int r=0; r<rates.length; r++)
				rates[r] += weight*subRates[r];
			totWeight += weight;
		}
		for (int r=0; r<rates.length; r++)
			rates[r] /= totWeight;
		InversionFaultSystemSolution meanSol = new InversionFaultSystemSolution(rupSet, rates);
		return meanSol;
	}

	public static void main(String[] args) throws ZipException, IOException {
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/tmp/comp_plots/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		Compound_FM_DM_Scale_CombinedFetcher subFetch =
				new Compound_FM_DM_Scale_CombinedFetcher(cfss, new APrioriBranchWeightProvider());
		
		subFetch.getSolution(subFetch.getBranches().iterator().next());
	}

}
