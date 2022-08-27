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

public class U3PerturbAndEnergyCalc {

	public static void main(String[] args) throws IOException {
		File resultsFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
//				+ "2022_03_25-u3_branches-orig_calc_params-FM3_1/results.zip");
//				+ "2022_03_24-u3_branches-FM3_1-2000ip/results.zip");
				+ "2022_03_25-u3_branches-orig_calc_params-new_avg-converged-FM3_1-2000ip/results.zip");
		
		ZipFile zip = new ZipFile(resultsFile);
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		MinMaxAveTracker pertubTrack = new MinMaxAveTracker();
		MinMaxAveTracker itersPerPerturbTrack = new MinMaxAveTracker();
		MinMaxAveTracker energyTrack = new MinMaxAveTracker();
		
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
			
			double totalE = progress.getEnergies(progress.size()-1)[0];
			long perturbs = progress.getNumPerturbations(progress.size()-1);
			long iters = progress.getIterations(progress.size()-1);
			pertubTrack.addValue(perturbs);
			itersPerPerturbTrack.addValue((double)iters/(double)perturbs);
			energyTrack.addValue(totalE);
		}
		
		System.out.println("Energy: "+energyTrack);
		System.out.println("Perturbations: "+pertubTrack);
		System.out.println("Iterations per perturbation: "+itersPerPerturbTrack);
		
		zip.close();
	}

}
