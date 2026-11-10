package scratch.kevin.mfdInversion;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.PosteriorSectionBValueDistributions;

import com.google.common.base.Preconditions;

public class PosteriorModuleReprocess {

	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("USAGE: PosteriorModuleReprocess <results-dir> <logic_tree.json> <ba-sol-file.zip>");
		}
		File resultsDir = new File(args[0]);
		File logicTreeFile = new File(args[1]);
		File baSolFile = new File(args[2]);
		
		LogicTree<?> tree = LogicTree.read(logicTreeFile);
	
		AveragingAccumulator<PosteriorSectionBValueDistributions> accumulator = null;
		for (LogicTreeBranch<?> branch : tree) {
			System.out.println("Processing "+branch.buildFileName());
			double weight = tree.getBranchWeight(branch);
			File solDir = new File(resultsDir, branch.buildFileName());
			Preconditions.checkState(solDir.exists(), "Solution dir doesn't exist: %s", solDir.getAbsolutePath());
			File solFile = new File(solDir, "solution.zip");
			Preconditions.checkState(solFile.exists(), "Solution file doesn't exist: %s", solFile.getAbsolutePath());
			ModuleArchive<OpenSHA_Module> archive = new ModuleArchive<>(solFile);
			PosteriorSectionBValueDistributions posterior = archive.loadUnlistedModule(PosteriorSectionBValueDistributions.class, "ruptures/");
			Preconditions.checkNotNull(posterior);
			if (accumulator == null)
				accumulator = posterior.averagingAccumulator();
			accumulator.process(posterior, weight);
		}
		
		PosteriorSectionBValueDistributions baPosterior = accumulator.getAverage();
		FaultSystemSolution baSol = FaultSystemSolution.load(baSolFile);
		baSol.getRupSet().addModule(baPosterior);
		String prefix = baSolFile.getName();
		if (prefix.endsWith(".zip"))
			prefix = prefix.substring(0, prefix.indexOf(".zip"));
		baSol.write(new File(baSolFile.getParentFile(), prefix+"_posterior.zip"));
	}

}
