package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.dom4j.DocumentException;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.util.FileNameComparator;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class SolutionMomentCalc {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_03_07-moment-reduction-variations");
//		File dir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/");
		File dir = new File("/tmp/comp_plots");
		
		File[] files = dir.listFiles();
		Arrays.sort(files, new FileNameComparator());
		
		for (File file : files) {
			if (!file.isFile())
				continue;
			String name = file.getName();
			if (!name.endsWith("_sol.zip") && !name.endsWith("BRANCH_AVG_SOL.zip"))
					continue;
			// we have a solution
			FaultSystemSolution sol = FaultSystemIO.loadSol(file);
			FaultSystemRupSet rupSet = sol.getRupSet();
			
			// calculate the moment
			double totalSolutionMoment = 0;
			double totalSolutionRate = 0;
			double meanMag = 0;
			for (int rup=0; rup<rupSet.getNumRuptures(); rup++) { 
				totalSolutionMoment += sol.getRateForRup(rup)*MagUtils.magToMoment(rupSet.getMagForRup(rup));
				totalSolutionRate += sol.getRateForRup(rup);
				meanMag += rupSet.getMagForRup(rup);
			}
			meanMag /= (double)rupSet.getNumRuptures();
			System.out.println(name);
			System.out.println("Total moment rate of solution = "+totalSolutionMoment);
			System.out.println("Total rate of solution = "+totalSolutionRate);
			System.out.println("Mean mag of solution = "+meanMag);
			System.out.println();
		}
	}

}
