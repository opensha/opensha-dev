package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

public class BranchNonZeroCount {
	
	public static void main(String[] args) throws IOException {
		File resultsFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift3km-"
				+ "ThreshAvgIterRelGR-IncludeThruCreep/results.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		MinMaxAveTracker fractTrack = new MinMaxAveTracker();
		MinMaxAveTracker numTrack = new MinMaxAveTracker();
		
		for (LogicTreeBranch<?> branch : slt.getLogicTree()) {
			double[] rates = slt.loadRatesForBranch(branch);
			int numNonZero = 0;
			for (double rate : rates)
				if (rate > 0)
					numNonZero++;
			double fract = (double)numNonZero/(double)rates.length;
			fractTrack.addValue(fract);
			numTrack.addValue(numNonZero);
		}
		
		System.out.println(fractTrack);
		System.out.println(numTrack);
	}

}
