package scratch.kevin.nshm23.uncertCorrFigures;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;

import com.google.common.base.Preconditions;

public class SolutionSlipRateCOV_CSV_Writer {

	public static void main(String[] args) throws IOException {
		Preconditions.checkArgument(args.length == 2, "USAGE: <inversion-dir> <csv-file>");
		File invDir = new File(args[0]);
		Preconditions.checkState(invDir.exists() && invDir.isDirectory());
		File resultsDir = new File(invDir, "results");
		Preconditions.checkState(resultsDir.exists() && resultsDir.isDirectory(),
				"Results subdirectory doesn't exist: %s", resultsDir.getAbsolutePath());
		File logicTreeFile = new File(invDir, "logic_tree.json");
		Preconditions.checkState(logicTreeFile.exists(), "Logic tree file doesn't exist: %s", logicTreeFile.getAbsolutePath());
		
		LogicTree<?> tree = LogicTree.read(logicTreeFile);
		
		double[] branchWeights = new double[tree.size()];
		double[][] sectSolSlipRates = null;
		
		for (int i=0; i<tree.size(); i++) {
			LogicTreeBranch<?> branch = tree.getBranch(i);
			branchWeights[i] = tree.getBranchWeight(i);
			
			File branchDir = new File(resultsDir, branch.buildFileName());
			Preconditions.checkState(branchDir.exists(), "Branch dir doesn't exist: %s", branchDir.getAbsolutePath());
			
			File solFile = new File(branchDir, "solution.zip");
			Preconditions.checkState(solFile.exists(), "Solution file doesn't exist: %s", solFile.getAbsolutePath());
			
			ZipFile zip = new ZipFile(solFile);
			ZipEntry entry = zip.getEntry("solution/sol_slip_rates.csv");
			
			System.out.println("Reading for branch "+i+"/"+tree.size());
			
			CSVFile<String> csv = CSVFile.readStream(new BufferedInputStream(zip.getInputStream(entry)), true);
			
			int numSects = csv.getNumRows()-1;
			if (sectSolSlipRates == null)
				sectSolSlipRates = new double[numSects][tree.size()];
			else
				Preconditions.checkState(numSects == sectSolSlipRates.length);
			
			for (int row=1; row<csv.getNumRows(); row++) {
				int r = row-1;
				Preconditions.checkState(r == csv.getInt(row, 0),
						"Must be 0-based and in order. Expected %s at row %s", r, row);
				sectSolSlipRates[r][i] = csv.getDouble(row, 1);
			}
			zip.close();
		}
		
		System.out.println("Calculating weighted std. devs.");
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Section Index", "Mean Slip Rate (m/yr)", "Slip Rate Std. Dev. (m/yr)");
		
		for (int s=0; s<sectSolSlipRates.length; s++) {
			Variance var = new Variance();
			double sd = Math.sqrt(new Variance().evaluate(sectSolSlipRates[s], branchWeights));
			double mean = new Mean().evaluate(sectSolSlipRates[s], branchWeights);
			csv.addLine(s+"", mean+"", sd+"");
		}
		
		csv.writeToFile(new File(args[1]));
	}

}
