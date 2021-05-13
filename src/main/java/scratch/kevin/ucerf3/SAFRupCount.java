package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class SAFRupCount {

	public static void main(String[] args) throws IOException, DocumentException {
		int parentID = 301;
		FaultSystemSolution fss = FaultSystemIO.loadSol(new File(
				"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		System.out.println("Count: "+fss.getRupSet().getRupturesForParentSection(parentID).size());
	}

}
