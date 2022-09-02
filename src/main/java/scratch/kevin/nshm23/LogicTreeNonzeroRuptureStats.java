package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.LogicTreeRateStatistics;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

public class LogicTreeNonzeroRuptureStats {
	
	public static void main(String[] args) throws IOException {
		File resultsFile, baFile;
		if (args.length == 0) {
			File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
					+ "2022_08_22-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
			resultsFile = new File(invDir, "results.zip");
			baFile = null;
		} else {
			if (args.length > 2) {
				System.err.println("USAGE: <results.zip> [<ba-file>]");
				System.exit(1);
			}
			resultsFile = new File(args[0]);
			if (args.length == 2)
				baFile = new File(args[1]);
			else
				baFile = null;
		}
		
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		LogicTreeRateStatistics stats = LogicTreeRateStatistics.forSolutionLogicTree(slt);
		System.out.println(stats.buildTable().toString());
		
		if (baFile != null) {
			FaultSystemSolution baSol = FaultSystemSolution.load(baFile);
			baSol.addModule(stats);
			baSol.write(baFile);
		}
	}

}
