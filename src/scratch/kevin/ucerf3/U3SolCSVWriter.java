package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class U3SolCSVWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip");
		File outputFile = new File("/tmp/fss_fm3p2.csv");
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
		CSVFile<String> csv = new CSVFile<>(false);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		csv.addLine("rupID", "mag", "annual rate", "aveRake", "area (m^2)", "length (m)", "numSectIDs", "sect1_ID", "sect2_ID", "...", "sectN_ID");
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			ArrayList<String> line = new ArrayList<>();
			line.add(r+"");
			line.add((float)rupSet.getMagForRup(r)+"");
			line.add((float)sol.getRateForRup(r)+"");
			line.add((float)rupSet.getAveRakeForRup(r)+"");
			line.add((float)rupSet.getAreaForRup(r)+"");
			line.add((float)rupSet.getLengthForRup(r)+"");
			List<Integer> sects = rupSet.getSectionsIndicesForRup(r);
			line.add(sects.size()+"");
			for (int s : sects)
				line.add(s+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(outputFile);
	}

}
