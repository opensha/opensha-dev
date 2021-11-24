package scratch.kevin.nshm23;

import org.apache.commons.cli.CommandLine;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;

public interface InversionConfigurationFactory {
	
	public FaultSystemRupSet buildRuptureSet(LogicTreeBranch<?> branch);
	
	public InversionConfiguration buildInversionConfig(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch, CommandLine cmd);

}
