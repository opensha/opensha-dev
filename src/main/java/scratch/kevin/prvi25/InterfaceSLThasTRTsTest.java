package scratch.kevin.prvi25;

import java.io.IOException;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import scratch.kevin.prvi25.figures.PRVI_Paths;

public class InterfaceSLThasTRTsTest {

	public static void main(String[] args) throws IOException {
		SolutionLogicTree slt = SolutionLogicTree.load(PRVI_Paths.SUBDUCTION_SLT);
		
		for (LogicTreeBranch<?> branch : slt.getLogicTree().getBranches()) {
			FaultSystemSolution sol = slt.forBranch(branch);
			RupSetTectonicRegimes trts = sol.getRupSet().requireModule(RupSetTectonicRegimes.class);
			for (int i=0; i<sol.getRupSet().getNumRuptures(); i++)
				Preconditions.checkState(trts.get(i) == TectonicRegionType.SUBDUCTION_INTERFACE);
		}
	}

}
