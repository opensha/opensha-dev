package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.CompletableFuture;
import java.util.function.Consumer;
import java.util.function.Supplier;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.modules.RuptureSubSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.SolutionProcessor;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class CherawRemoveAndReprocess {
	
	private static final int CHERAW_ID = 2180;

	public static void main(String[] args) throws IOException {
		File sourceResultsDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/"
				+ "2022_12_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/results");
		File destResultsDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/results");
		
		Preconditions.checkState(sourceResultsDir.exists(), "Source doesn't exist: %s", sourceResultsDir.getAbsolutePath());
		Preconditions.checkState(destResultsDir.exists() || destResultsDir.mkdir(),
				"Destination doesn't exist and couldn't be created: %s", destResultsDir.getAbsolutePath());
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		
		factory.setCacheDir(new File("/project/scec_608/kmilner/nshm23/batch_inversions/cache/"));
		
		CompletableFuture<FaultSystemSolution> processFuture = null;
		CompletableFuture<Void> writeFuture = null;
		
		System.out.println("Listing directory...");
		File[] subDirs = sourceResultsDir.listFiles();
		System.out.println("Sorting...");
		Arrays.sort(subDirs, new FileNameComparator());
		System.out.println("Processing!");
		
		boolean mpj = false;
		for (String arg : args) {
			System.out.println("arg: "+arg);
			if (arg.trim().toLowerCase().contains("mpj"))
				mpj = true;
		}
		
		if (mpj) {
			System.out.println("Processing with MPJ");
			args = ParallelDriver.doInit(args);
			
			ParallelDriver driver = new ParallelDriver(factory, subDirs, destResultsDir);
			
			try {
				driver.run();
			} catch (IOException | InterruptedException e) {
				e.printStackTrace();
				System.exit(1);
			}
			
			ParallelDriver.doFinalize();
			
			System.exit(0);
		}
		
		for (File sourceDir : subDirs) {
			File sourceSolFile = new File(sourceDir, "solution.zip");
			if (!sourceSolFile.exists())
				continue;
			
			System.out.println("Processing "+sourceDir.getName());
			
			File destDir = new File(destResultsDir, sourceDir.getName());
			Preconditions.checkState(destDir.exists() || destDir.mkdir(),
					"Destination doesn't exist and couldn't be created: %s", destDir.getAbsolutePath());
			
			File destSolFile = new File(destDir, sourceSolFile.getName());
			
			FaultSystemSolution sourceSol = FaultSystemSolution.load(sourceSolFile);
			
			if (processFuture != null)
				processFuture.join();
			
			processFuture = CompletableFuture.supplyAsync(new Processor(sourceSol, factory));
			
			if (writeFuture != null)
				writeFuture.join();
			
			writeFuture = processFuture.thenAccept(new Writer(destSolFile));
		}
		
		if (processFuture != null)
			processFuture.join();
		
		if (writeFuture != null)
			writeFuture.join();
		
		System.out.println("DONE!");
	}
	
	private static class ParallelDriver extends MPJTaskCalculator {
		
		private NSHM23_InvConfigFactory factory;
		private File[] subDirs;
		private File destResultsDir;
		
		private static String[] doInit(String[] args) {
			return MPJTaskCalculator.initMPJ(args);
		}

		protected static void doFinalize() {
			finalizeMPJ();
		}

		public ParallelDriver(NSHM23_InvConfigFactory factory, File[] subDirs, File destResultsDir) {
			super(1, 1, 50, false);
			this.shuffle = true;
			this.factory = factory;
			this.subDirs = subDirs;
			this.destResultsDir = destResultsDir;
		}

		@Override
		protected int getNumTasks() {
			return subDirs.length;
		}

		@Override
		protected void calculateBatch(int[] batch) throws Exception {
			CompletableFuture<FaultSystemSolution> processFuture = null;
			CompletableFuture<Void> writeFuture = null;
			
			for (int index : batch) {
				File sourceDir = subDirs[index];
				
				File sourceSolFile = new File(sourceDir, "solution.zip");
				if (!sourceSolFile.exists())
					continue;
				
				debug("Processing "+sourceDir.getName());
				
				File destDir = new File(destResultsDir, sourceDir.getName());
				Preconditions.checkState(destDir.exists() || destDir.mkdir(),
						"Destination doesn't exist and couldn't be created: %s", destDir.getAbsolutePath());
				
				File destSolFile = new File(destDir, sourceSolFile.getName());
				
				if (destSolFile.exists())
					continue;
				
				FaultSystemSolution sourceSol = FaultSystemSolution.load(sourceSolFile);
				
				if (processFuture != null)
					processFuture.join();
				
				processFuture = CompletableFuture.supplyAsync(new Processor(sourceSol, factory));
				
				if (writeFuture != null)
					writeFuture.join();
				
				writeFuture = processFuture.thenAccept(new Writer(destSolFile));
			}
			
			if (processFuture != null)
				processFuture.join();
			
			if (writeFuture != null)
				writeFuture.join();
			
			debug("finished batch of size "+batch.length);
		}

		@Override
		protected void doFinalAssembly() throws Exception {}
		
	}
	
	private static class Processor implements Supplier<FaultSystemSolution> {
		
		private FaultSystemSolution sourceSol;
		private InversionConfigurationFactory factory;

		public Processor(FaultSystemSolution sourceSol, InversionConfigurationFactory factory) {
			super();
			this.sourceSol = sourceSol;
			this.factory = factory;
		}

		@Override
		public FaultSystemSolution get() {
			LogicTreeBranch<?> branch = null;
			try {
				// remove cheraw fault
				FaultSystemRupSet sourceRupSet = sourceSol.getRupSet();
				
				HashSet<Integer> retained = new HashSet<>();
				for (FaultSection sect : sourceRupSet.getFaultSectionDataList())
					if (sect.getParentSectionId() != CHERAW_ID)
						retained.add(sect.getSectionId());
				
				FaultSystemRupSet mappedRupSet = sourceRupSet.getForSectionSubSet(retained);
				RuptureSubSetMappings mappings = mappedRupSet.requireModule(RuptureSubSetMappings.class);
				double[] modRates = new double[mappedRupSet.getNumRuptures()];
				for (int r=0; r<modRates.length; r++)
					modRates[r] = sourceSol.getRateForRup(mappings.getOrigRupID(r));
				
				branch = sourceRupSet.requireModule(LogicTreeBranch.class);
				
				// load cached rupture set
				FaultSystemRupSet modRupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
				Preconditions.checkState(modRupSet.getNumSections() == mappedRupSet.getNumSections(),
						"Section count mismatch: %s != %s", modRupSet.getNumSections(), mappedRupSet.getNumSections());
				Preconditions.checkState(modRupSet.getNumRuptures() == mappedRupSet.getNumRuptures(),
						"Rupture count mismatch: %s != %s", modRupSet.getNumRuptures(), mappedRupSet.getNumRuptures());
				
				System.out.println("Processing rupture set");
				SolutionProcessor processor = factory.getSolutionLogicTreeProcessor();
				// build an inversion configuration so that any adjustments at that stage are processed
				factory.buildInversionConfig(modRupSet, branch, 8);
				
				FaultSystemSolution modSol = new FaultSystemSolution(modRupSet, modRates);
				// process the solution
				modSol = processor.processSolution(modSol, branch);
				return modSol;
			} catch (Throwable e) {
				System.err.println("Exception for "+branch+"!");
				e.printStackTrace();
				System.err.flush();
				System.exit(1);
				return null;
			}
		}
		
	}
	
	private static class Writer implements Consumer<FaultSystemSolution> {
		
		private File destFile;

		public Writer(File destFile) {
			super();
			this.destFile = destFile;
		}

		@Override
		public void accept(FaultSystemSolution sol) {
			try {
				// write out the reprocessed solution
				sol.write(destFile);
			} catch (Throwable e) {
				System.err.println("Exception writing "+destFile.getAbsolutePath()+"!");
				e.printStackTrace();
				System.err.flush();
				System.exit(1);
			}
		}
		
	}

}
