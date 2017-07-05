package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.util.FileNameComparator;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.inversion.BatchPlotGen;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;

import java.util.Arrays;
import java.util.zip.ZipException;

public class MultiRunPlotBuilder {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws ZipException 
	 * @throws GMT_MapException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException, GMT_MapException {
		if (args.length != 4) {
			System.out.println("USAGE: <dir> <prefix> <rup set file> <num>");
			System.exit(2);
		}
		
		File mainDir = new File(args[0]);
		String prefix = args[1];
		File rupSetFile = new File(args[2]);
		int num = Integer.parseInt(args[3]);
		
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(rupSetFile);
		
		File[] files = mainDir.listFiles();
		Arrays.sort(files, new FileNameComparator());
		
		int count = 0;
		
		for (File dir : files) {
			if (!dir.isDirectory())
				continue;
			if (!dir.getName().startsWith(prefix))
				continue;
			if (dir.getName().contains("mean"))
				continue;
			
			if (count++ == num)
				break;
			
			File solFile = new File(dir, dir.getName()+"_sol.zip");
			if (!solFile.exists()) {
				File ratesFile = new File(dir, dir.getName()+".bin");
				
				double[] rates = MatrixIO.doubleArrayFromFile(ratesFile);
				FaultSystemSolution sol = FaultSystemSolution.buildSolAsApplicable(rupSet, rates);
				
				FaultSystemIO.writeSol(sol, solFile);
			}
			
			BatchPlotGen.handleDir(dir, null, 1);
		}
	}

}
