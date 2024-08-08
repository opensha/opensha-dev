package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

import com.google.common.base.Preconditions;

public class SolutionLTZipToResultsDir {

	public static void main(String[] args) throws IOException {
		if (args.length < 3 || args.length > 4) {
			System.err.println("USAGE: <slt file> <output results dir> <process?> [<rewrite?>]");
			System.exit(2);
		}
		try {
			ModuleContainer.VERBOSE_DEFAULT = false;
			File sltFile = new File(args[0]);
			Preconditions.checkState(sltFile.exists());
			File outputDir = new File(args[1]);
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
			
			boolean process = Boolean.parseBoolean(args[2]);
			System.out.println("Process? "+process);
			boolean rewrite = args.length < 4 ? false : Boolean.parseBoolean(args[3]);
			
			SolutionLogicTree slt = SolutionLogicTree.load(sltFile);
			
			ArrayDeque<CompletableFuture<Void>> writeFutures = new ArrayDeque<>();
			int maxConcurrentWrites = 10;
			
			for (int i=0; i<slt.getLogicTree().size(); i++) {
				LogicTreeBranch<?> branch = slt.getLogicTree().getBranch(i);
				System.out.println("Loading branch "+i+"/"+slt.getLogicTree().size()+": "+branch);
				String dirName = branch.buildFileName();
				File runDir = new File(outputDir, dirName);
				Preconditions.checkState(runDir.exists() || runDir.mkdir());
				
				File solFile = new File(runDir, "solution.zip");
				
				if (!rewrite && solFile.exists()) {
					System.out.println("Solution already exists, skipping because rewrite=false");
					continue;
				}
				
				FaultSystemSolution sol = slt.forBranch(branch, process);
				if (writeFutures.size() == maxConcurrentWrites)
					writeFutures.removeFirst().join();
				
				int branchIndex = i;
				writeFutures.add(CompletableFuture.runAsync(new Runnable() {
					
					@Override
					public void run() {
						try {
							sol.write(solFile);
							System.out.println("DONE writing branch "+branchIndex+"/"+slt.getLogicTree().size());
						} catch (IOException e) {
							e.printStackTrace();
							System.exit(1);
						}
					}
				}));
			}
			
			while (!writeFutures.isEmpty())
				writeFutures.removeFirst().join();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		System.out.println("DONE");
	}

}
