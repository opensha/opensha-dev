package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;

import com.google.common.base.Preconditions;

public class TargetMFDsBAReprocess {

	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("USAGE: TargetMFDsBAReprocess <results-dir> <logic_tree.json> <ba-sol-file.zip>");
		}
		File resultsDir = new File(args[0]);
		File logicTreeFile = new File(args[1]);
		File baSolFile = new File(args[2]);
		
		ExecutorService exec = Executors.newFixedThreadPool(5);
		
		LinkedList<CompletableFuture<WeightedModule>> futures = new LinkedList<>();
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		
		LogicTree<?> tree = LogicTree.read(logicTreeFile);
		
		// start loading the BA solution in the backround
		CompletableFuture<FaultSystemSolution> baFuture = CompletableFuture.supplyAsync(() ->{
			FaultSystemSolution baSol;
			try {
				baSol = FaultSystemSolution.load(baSolFile);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			baSol.getRupSet().loadAllAvailableModules();
			baSol.loadAllAvailableModules();
			return baSol;
		});
	
		for (int b=0; b<tree.size(); b++) {
			LogicTreeBranch<?> branch = tree.getBranch(b);
			final int branchIndex = b;
			double weight = tree.getBranchWeight(branch);
			File solDir = new File(resultsDir, branch.buildFileName());
			Preconditions.checkState(solDir.exists(), "Solution dir doesn't exist: %s", solDir.getAbsolutePath());
			File solFile = new File(solDir, "solution.zip");
			Preconditions.checkState(solFile.exists(), "Solution file doesn't exist: %s", solFile.getAbsolutePath());
			futures.add(CompletableFuture.supplyAsync(()->{
				System.out.println("Processing "+branchIndex+"/"+tree.size()+":\t"+branch.buildFileName());
				FaultSystemRupSet rupSet;
				try {
					rupSet = FaultSystemRupSet.load(solFile);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				return new WeightedModule(rupSet.requireModule(InversionTargetMFDs.class), weight);
			}, exec));
		}
		AveragingAccumulator<InversionTargetMFDs> accumulator = null;
		while (!futures.isEmpty()) {
			// doing it this way will allow old ones to be garbage collected after they're consumed
			CompletableFuture<WeightedModule> future = futures.removeFirst();
			WeightedModule module = future.join();
			if (accumulator == null)
				accumulator = module.module.averagingAccumulator();
			accumulator.process(module.module, module.weight);
		}
		
		exec.shutdown();
		
		InversionTargetMFDs baModule = accumulator.getAverage();
		FaultSystemSolution baSol = baFuture.join();
		baSol.getRupSet().addModule(baModule);
		String prefix = baSolFile.getName();
		if (prefix.endsWith(".zip"))
			prefix = prefix.substring(0, prefix.indexOf(".zip"));
		baSol.write(new File(baSolFile.getParentFile(), prefix+"_mfds.zip"));
	}
	
	private record WeightedModule(InversionTargetMFDs module, double weight) {}

}
