package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;

import com.google.common.base.Preconditions;

public class ResultsDirFaultModelModulesReprocess {

	public static void main(String[] args) throws IOException {
		if (args.length != 2 && args.length != 3) {
			System.err.println("USAGE: <results-dir> <logic-tree JSON> [<branch-average-output-file>]");
			System.exit(1);
		}
		File resultsDir = new File(args[0]);
		Preconditions.checkState(resultsDir.exists() && resultsDir.isDirectory());
		
		LogicTree<?> tree = LogicTree.read(new File(args[1]));
		BranchAverageSolutionCreator baCreator = args.length == 3 ? new BranchAverageSolutionCreator(tree.getWeightProvider()) : null;
		
		CompletableFuture<Void> writeFuture = null;
		CompletableFuture<Void> baFuture = null;
		
		int branchIndex = 0;
		for (LogicTreeBranch<?> branch : tree) {
			System.out.println("Branch "+(branchIndex++)+"/"+tree.size()+": "+branch);
			File solDir = new File(resultsDir, branch.buildFileName());
			Preconditions.checkState(solDir.exists(), "No dir for branch: %s", solDir.getAbsolutePath());
			File solFile = new File(solDir, "solution.zip");
			Preconditions.checkState(solFile.exists(), "No solution file for branch: %s", solFile.getAbsolutePath());
			System.out.println("Processing: "+solDir.getAbsolutePath());
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			
			RupSetFaultModel fm = branch.requireValue(RupSetFaultModel.class);
			fm.attachDefaultModules(sol.getRupSet());
			
			if (writeFuture != null)
				writeFuture.join();
			writeFuture = CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					try {
						sol.write(solFile);
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
				}
			});
			if (baFuture != null)
				baFuture.join();
			if (baCreator != null) {
				baFuture = CompletableFuture.runAsync(new Runnable() {
					
					@Override
					public void run() {
						baCreator.addSolution(sol, branch);
					}
				});
			}
		}
		if (writeFuture != null)
			writeFuture.join();
		if (baFuture != null)
			baFuture.join();
		System.out.println("Done reprocessing!");
		if (baCreator != null) {
			// build/write it
			System.out.println("Building branch-averaged solution...");
			FaultSystemSolution baSol = baCreator.build();
			System.out.println("Writing branch-averaged solution...");
			baSol.write(new File(args[2]));
			System.out.println("Done with BA!");
		}
	}

}
