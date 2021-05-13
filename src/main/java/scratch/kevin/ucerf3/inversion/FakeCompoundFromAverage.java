package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.VariableLogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;

public class FakeCompoundFromAverage {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File dir = new File("/tmp");
		File inFile = new File(dir,
//				"FM3_1_ZENG_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
//				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
		File outFile = new File(dir, "3p3_convergence_compound.zip");
		
		final AverageFaultSystemSolution avgSol = FaultSystemIO.loadAvgInvSol(inFile);
		final LogicTreeBranch branch = LogicTreeBranch.fromFileName(inFile.getName());
		
		final int bundle = 1;
		
		final List<LogicTreeBranch> branches = Lists.newArrayList();
		
		if (bundle <= 1) {
			for (int s=0; s<avgSol.getNumSolutions(); s++)
				branches.add(new VariableLogicTreeBranch(branch, Lists.newArrayList("Var"+s)));
		} else {
			Preconditions.checkState(avgSol.getNumSolutions() % bundle == 0);
			int numBundles = avgSol.getNumSolutions() / bundle;
			for (int s=0; s<numBundles; s++)
				branches.add(new VariableLogicTreeBranch(branch, Lists.newArrayList("Var"+s)));
		}
		
		FaultSystemSolutionFetcher fetcher = new FaultSystemSolutionFetcher() {
			
			@Override
			public Collection<LogicTreeBranch> getBranches() {
				return branches;
			}
			
			@Override
			protected InversionFaultSystemSolution fetchSolution(LogicTreeBranch branch) {
				String var = ((VariableLogicTreeBranch)branch).getVariations().get(0);
				int s = Integer.parseInt(var.substring(3));
				if (bundle <= 1)
					return avgSol.getSolution(s);
				// this means we're bundling
				int bundleStart = s * bundle;
				List<double[]> ratesList = Lists.newArrayList();
				for (int i=bundleStart; i<bundleStart+bundle; i++)
					ratesList.add(avgSol.getRates(i));
				return new AverageFaultSystemSolution(avgSol.getRupSet(), ratesList);
			}
		};
		
		CompoundFaultSystemSolution.toZipFile(outFile, fetcher);
	}

}
