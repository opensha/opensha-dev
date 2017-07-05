package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.zip.ZipException;

import org.opensha.commons.util.ExceptionUtils;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.griddedSeismicity.GridSourceFileReader;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;

public class CompoundGriddedRecalc {

	public static void main(String[] args) throws ZipException, IOException {
		File compoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
		File writeDir = new File("/tmp/cfss_branches");
		if (!writeDir.exists())
			writeDir.mkdir();
		
		ArrayDeque<LogicTreeBranch> branchQueue = new ArrayDeque<LogicTreeBranch>();
		branchQueue.addAll(cfss.getBranches());
		
		int numThreads = 4;
		
		List<CalcThread> threads = Lists.newArrayList();
		
		for (int i=0; i<numThreads; i++)
			threads.add(new CalcThread(cfss, writeDir, branchQueue));
		
		for (CalcThread thread : threads)
			thread.start();
		
		for (CalcThread thread : threads)
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
	}
	
	private static synchronized LogicTreeBranch popBranch(ArrayDeque<LogicTreeBranch> branches) {
		if (branches.size() % 50 == 0)
			System.out.println(branches.size()+" braches left");
		return branches.pop();
	}
	
	private static boolean region_written = false;
	
	private static synchronized boolean isRegionWritten() {
		if (!region_written) {
			region_written = true;
			return false;
		}
		return true;
	}
	
	private static class CalcThread extends Thread {
		private ArrayDeque<LogicTreeBranch> branches;
		private File dir;
		private CompoundFaultSystemSolution cfss;
		
		public CalcThread(CompoundFaultSystemSolution cfss, File dir, ArrayDeque<LogicTreeBranch> branches) {
			this.branches = branches;
			this.dir = dir;
			this.cfss = cfss;
		}
		
		@Override
		public void run() {
			while (true) {
				LogicTreeBranch branch;
				try {
					branch = popBranch(branches);
				} catch (NoSuchElementException e) {
					// done
					return;
				}
				
				InversionFaultSystemSolution sol = cfss.getSolution(branch);
				GridSourceProvider gridSources = sol.getGridSourceProvider();
				
				File regXMLFile = null;
				if (!isRegionWritten())
					regXMLFile = new File(dir, CompoundFaultSystemSolution.getRemappedName("grid_sources_reg.xml", branch));
				File binFile = new File(dir, CompoundFaultSystemSolution.getRemappedName("grid_sources.bin", branch));
				try {
//					GridSourceFileReader.writeGriddedSeisFile(file, gridSources);
					GridSourceFileReader.writeGriddedSeisBinFile(binFile, regXMLFile, gridSources,
							AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF);
				} catch (IOException e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
			}
		}
	}

}
