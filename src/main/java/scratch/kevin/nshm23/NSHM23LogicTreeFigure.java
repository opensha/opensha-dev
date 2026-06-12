package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;

public class NSHM23LogicTreeFigure {

	public static void main(String[] args) throws IOException {
		LogicTree<?> tree = LogicTree.buildExhaustive(NSHM23_LogicTreeBranch.levelsOnFault, true);
		LogicTreeFigureWriter treeFig = new LogicTreeFigureWriter(tree, false, true);
		treeFig.write(new File("/tmp/"), "nshm23_on_fault", true, false);
		
		tree = LogicTree.buildExhaustive(NSHM23_LogicTreeBranch.levelsOffFault, true);
		treeFig = new LogicTreeFigureWriter(tree, false, true);
		treeFig.write(new File("/tmp/"), "nshm23_off_fault", true, false);
		
		tree = LogicTree.buildExhaustive(NSHM23_LogicTreeBranch.levelsCombined, true);
		treeFig = new LogicTreeFigureWriter(tree, false, true);
		treeFig.write(new File("/tmp/"), "nshm23_combined", true, false);
	}

}
