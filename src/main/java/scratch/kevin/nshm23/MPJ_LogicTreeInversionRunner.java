package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_LogicTreeInversionRunner extends MPJTaskCalculator {

	private File outputDir;
	
	private int runsPerBranch = 1;
	
	private InversionConfigurationFactory factory;

	private LogicTree<LogicTreeNode> tree;

	private CommandLine cmd;

	public MPJ_LogicTreeInversionRunner(CommandLine cmd) throws IOException {
		super(cmd);
		this.cmd = cmd;
		
		tree = LogicTree.read(new File(cmd.getOptionValue("logic-tree")));
		if (rank == 0)
			debug("Loaded "+tree.size()+" tree nodes");
		
		outputDir = new File(cmd.getOptionValue("output-dir"));
		
		if (rank == 0)
			waitOnDir(outputDir, 5, 1000);
		
		if (cmd.hasOption("runs-per-branch"))
			runsPerBranch = Integer.parseInt(cmd.getOptionValue("runs-per-branch"));
		
		try {
			@SuppressWarnings("unchecked")
			Class<? extends InversionConfigurationFactory> factoryClass = (Class<? extends InversionConfigurationFactory>)
					Class.forName(cmd.getOptionValue("inversion-factory"));
			factory = factoryClass.getDeclaredConstructor().newInstance();
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	private static void waitOnDir(File dir, int maxRetries, long sleepMillis) {
		int retry = 0;
		while (!(dir.exists() || dir.mkdir())) {
			try {
				Thread.sleep(sleepMillis);
			} catch (InterruptedException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (retry++ > maxRetries)
				throw new IllegalStateException("Directory doesn't exist and couldn't be created after "
						+maxRetries+" retries: "+dir.getAbsolutePath());
		}
	}

	@Override
	protected int getNumTasks() {
		return tree.size()*runsPerBranch;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int index : batch) {
			int branchIndex = index / runsPerBranch;
			int run = index % runsPerBranch;
			
			LogicTreeBranch<LogicTreeNode> branch = tree.getBranch(branchIndex);
			
			debug("index "+index+" is branch "+branchIndex+" run "+run+": "+branch);
			
			String dirName = branch.buildFileName();
			if (runsPerBranch > 1)
				dirName += "_run "+run;
			File runDir = new File(outputDir, dirName);
			Preconditions.checkState(runDir.exists() || runDir.mkdir());
			
			File solFile = new File(runDir, "solution.zip");
			
			if (solFile.exists()) {
				debug(solFile.getAbsolutePath()+" exists, testing loading...");
				try {
					FaultSystemSolution.load(solFile);
					debug("skipping "+index+" (already done)");
					continue;
				} catch (Exception e) {
					debug("Failed to load, re-inverting: "+e.getMessage());
				}
			}
			
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch);
			rupSet.addModule(branch);
			
			InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, cmd);
			
			debug("Running inversion for task "+index);
			FaultSystemSolution sol = Inversions.run(rupSet, config);
			
			sol.write(solFile);
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		// TODO Auto-generated method stub

	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		ops.addRequiredOption("lt", "logic-tree", true, "Path to logic tree JSON file");
		ops.addRequiredOption("od", "output-dir", true, "Path to output directory");
		ops.addOption("rpb", "runs-per-branch", true, "Runs per branch (default is 1)");
		ops.addRequiredOption("ifc", "inversion-factory", true, "Inversion configuration factory classname");
		
		for (Option op : InversionConfiguration.createSAOptions().getOptions())
			ops.addOption(op);
		
		return ops;
	}

	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_LogicTreeInversionRunner.class);
			
			MPJ_LogicTreeInversionRunner driver = new MPJ_LogicTreeInversionRunner(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
