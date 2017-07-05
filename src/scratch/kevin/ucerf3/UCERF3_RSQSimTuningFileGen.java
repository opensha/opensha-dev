package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class UCERF3_RSQSimTuningFileGen {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		String solName = "fm3p1_branch_avg";
//		FaultSystemSolution sol = FaultSystemIO.loadSol(
//				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
//						+ "FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip"));
//		String solName = "fm3p1_geol";
		
		double minMag = 0d;
		
		double[] rates = sol.calcParticRateForAllSects(minMag, 10d);
		
		String magStr;
		String magFileStr;
		if (minMag > 0) {
			String m;
			if (minMag == (double)((int)minMag))
				m = (int)minMag+"";
			else
				m = (float)minMag+"";
			magStr = "M>="+m;
			magFileStr = "m"+m;
		} else {
			magStr = "Supra-Seismogenic";
			magFileStr = "supra_seis";
		}
		
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Subsection Index", "Parent Section ID", "Subsection Name",
				magStr+" Participation Rate", magStr+" Participation RI");
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		Preconditions.checkState(rates.length == rupSet.getNumSections());
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSectionPrefData sect = rupSet.getFaultSectionData(s);
			
			double ri = 1d/rates[s];
			
			csv.addLine(s+"", sect.getParentSectionId()+"", sect.getName(), rates[s]+"", ri+"");
		}
		
		String fName = solName+"_sub_sect_"+magFileStr+"_ris_for_tuning.csv";
		csv.writeToFile(new File("/tmp/"+fName));
	}

}
