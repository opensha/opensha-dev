package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.AnnealingProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

import com.google.common.base.Preconditions;

public class BranchWorseKeptCount {
	
	public static void main(String[] args) throws IOException {
		File resultsFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_06_10-nshm23_u3_hybrid_branches-full_sys_inv-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-"
				+ "Shift2km-ThreshAvgIterRelGR-IncludeThruCreep/results.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		ZipFile zip = new ZipFile(resultsFile);
		
		MinMaxAveTracker fractTrack = new MinMaxAveTracker();
		MinMaxAveTracker numTrack = new MinMaxAveTracker();
		
		for (LogicTreeBranch<?> branch : slt.getLogicTree()) {
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(AnnealingProgress.PROGRESS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += AnnealingProgress.PROGRESS_FILE_NAME;
			System.out.println("Loading "+entryName);
			ZipEntry entry = zip.getEntry(entryName);
			Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			AnnealingProgress progress = new AnnealingProgress(csv);
			
			int num = progress.size();
			long numIters = progress.getIterations(num-1);
			long worseKept = progress.getNumWorseKept(num-1);
			double fract = (double)worseKept/(double)numIters;
			fractTrack.addValue(fract);
			numTrack.addValue(worseKept);
		}
		zip.close();
		
		System.out.println(fractTrack);
		System.out.println(numTrack);
	}

}
