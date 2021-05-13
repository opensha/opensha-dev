package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.FaultSystemIO;

public class BatchSolutionWrite {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_30-fm2-a-priori-test/results/VarAPrioriZero_VarAPrioriWt1000_VarWaterlevel0");
		File odir = new File("/home/kevin/OpenSHA/UCERF3/eal/2012_05_03-eal-tests-apriori-1000");
		File rupSetFile = new File(dir, "rupSet.zip");
		
		InversionFaultSystemRupSet rupSet = FaultSystemIO.loadInvRupSet(rupSetFile);
		
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (!name.endsWith(".bin"))
				continue;
			if (!name.contains("run"))
				continue;
			
			File zipFile = new File(odir, name.replaceAll(".bin", ".zip"));
			if (zipFile.exists())
				continue;
			
			double[] rates = MatrixIO.doubleArrayFromFile(file);
			
			InversionFaultSystemSolution sol = new InversionFaultSystemSolution(rupSet, rates);
			FaultSystemIO.writeSol(sol, zipFile);
		}
	}

}
