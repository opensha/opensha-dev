package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.NucleationRatePlot;

import com.google.common.base.Preconditions;

public class GriddedBranchMomentRatesFileWriter {

	public static void main(String[] args) throws IOException {
		File sltFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_gridded_branches_simplified.zip");
		SolutionLogicTree slt = SolutionLogicTree.load(sltFile);
		
		File outputDir = new File("/tmp/gridded_moment_rate_branches");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		for (LogicTreeBranch<?> branch : slt.getLogicTree()) {
			System.out.println("Branch: "+branch);
			GridSourceProvider gridProv =  slt.loadGridProvForBranch(branch);
			GriddedGeoDataSet xyz = NucleationRatePlot.calcGriddedNucleationMomentRates(gridProv);
			GriddedGeoDataSet.writeXYZFile(xyz, new File(outputDir, branch.buildFileName()+".xyz"));
		}
	}

}
