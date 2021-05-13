package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.util.FileUtils;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class CompoundAddToMean {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException {
		if (args.length != 4) {
			System.out.println("USAGE: CompoundAddToMean <mean-sol> <new-sol> <num orig> <out-file>");
			System.exit(2);
		}
		
		File meanFile = new File(args[0]);
		File newFile = new File(args[1]);
		int numOrig = Integer.parseInt(args[2]);
		Preconditions.checkArgument(numOrig > 0, "numOrig must be > 0");
		File outFile = new File(args[3]);
		
		CompoundFaultSystemSolution meanSol = new CompoundFaultSystemSolution(new ZipFile(meanFile));
		CompoundFaultSystemSolution newSol = new CompoundFaultSystemSolution(new ZipFile(newFile));
		
		File tempDir = FileUtils.createTempDir();
		FileUtils.unzipFile(meanFile, tempDir);
		
		for (LogicTreeBranch branch : meanSol.getBranches()) {
			double[] meanRates = meanSol.getRates(branch);
			double[] newRates = newSol.getRates(branch);
			
			double[] combinedRates = new double[meanRates.length];
			for (int i=0; i<combinedRates.length; i++)
				combinedRates[i] = (meanRates[i] * (double)numOrig + newRates[i]) / (double)(numOrig+1);
			
			String ratesFileName = CompoundFaultSystemSolution.getRemappedName("rates.bin", branch);
			File ratesFile = new File(tempDir, ratesFileName);
			Preconditions.checkState(ratesFile.exists());
			MatrixIO.doubleArrayToFile(combinedRates, ratesFile);
		}
		
		List<String> fileNames = Lists.newArrayList();
		for (File f : tempDir.listFiles()) {
			if (f.isDirectory())
				continue;
			if (f.getName().startsWith("."))
				continue;
			fileNames.add(f.getName());
		}
		FileUtils.createZipFile(outFile.getAbsolutePath(), tempDir.getAbsolutePath(), fileNames);
	}

}
