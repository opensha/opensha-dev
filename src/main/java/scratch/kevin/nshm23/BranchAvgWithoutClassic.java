package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;

import com.google.common.base.Preconditions;

public class BranchAvgWithoutClassic {

	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("USAGE: <results-dir> <logic-tree.json> <output-file>");
			System.exit(1);
		}
		
		File resultsDir = new File(args[0]);
		
		LogicTree<?> tree = LogicTree.read(new File(args[1]));
		
		tree = tree.matchingNone(NSHM23_SegmentationModels.CLASSIC);
		
		System.out.println("Retained "+tree.size()+" logic tree branches");
		
		File outputFile = new File(args[2]);
		
		BranchAverageSolutionCreator baCreator = new BranchAverageSolutionCreator(tree.getWeightProvider());
		
		CompletableFuture<Void> processFuture = null;
		
		int doneCount = 0;
		
		for (LogicTreeBranch<?> branch : tree) {
			File dir = new File(resultsDir, branch.buildFileName());
			Preconditions.checkState(dir.exists());
			
			File solFile = new File(dir, "solution.zip");
			
			Preconditions.checkState(solFile.exists());
			
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			
			if (processFuture != null) {
				processFuture.join();
				doneCount++;
				System.out.println("DONE processing solution "+doneCount+"/"+tree.size());
			}
			
			final LogicTreeBranch<?> myBranch = branch;
			
			processFuture = CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					baCreator.addSolution(sol, myBranch);
				}
			});
		}
		
		processFuture.join();
		
		FaultSystemSolution baSol = baCreator.build();
		
		baSol.write(outputFile);
	}

}
