package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;

public class BA_SolAvgGridProvAttach {

	public static void main(String[] args) throws IllegalStateException, IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_16-nshm23-all_ba-AVERAGE_NSHM23_Avg_AvgSupraB_NoRed_AverageFitPaleo_AvgSeg");
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File(dir, "solution.zip"));
		
		GridSourceProvider gridProv = NSHM23_InvConfigFactory.buildBranchAveragedGridSourceProv(
				sol, sol.getRupSet().requireModule(LogicTreeBranch.class));
		
		sol.setGridSourceProvider(gridProv);
		
		sol.write(new File(dir, "solution_gridded.zip"));
	}

}
