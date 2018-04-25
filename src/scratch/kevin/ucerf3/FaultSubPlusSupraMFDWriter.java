package scratch.kevin.ucerf3;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.utils.FaultSystemIO;

public class FaultSubPlusSupraMFDWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		MeanUCERF3 erf = new MeanUCERF3();
		Presets[] presets = {Presets.FM3_1_MAG_VAR, Presets.FM3_2_MAG_VAR};
		String[] prefixes = {"fm3_1", "fm3_2"};
		File regDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions");
		File[] subSeismoSolFiles = {new File(regDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"),
				new File(regDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip")};
		File outputDir = new File("/tmp");
		
		int numMag = 90;
		double minMag = 0.05;
		double maxMag = 8.95;
		
		for (int p=0; p<presets.length; p++) {
			erf.setPreset(presets[p]);
			erf.updateForecast();
			FaultSystemSolution sol = erf.getSolution();
			FaultSystemRupSet rupSet = sol.getRupSet();
//			FaultSystemIO.writeSol(sol, new File("/tmp/"+prefixes[p]+"_sol.zip"));
			Map<String, IncrementalMagFreqDist> supraParentMFDs = new HashMap<>();
			Map<String, Integer> parentIDsMap = new HashMap<>();
			for (int s=0; s<rupSet.getNumSections(); s++) {
				FaultSectionPrefData fsd = rupSet.getFaultSectionData(s);
				parentIDsMap.put(fsd.getParentSectionName(), fsd.getParentSectionId());
			}
			for (String parentName : parentIDsMap.keySet())
				supraParentMFDs.put(parentName, sol.calcParticipationMFD_forParentSect(
						parentIDsMap.get(parentName), minMag, maxMag, numMag));
			// now load in sub siesmo MFDs from original sol
			FaultSystemSolution solWithSubSeis = FaultSystemIO.loadSol(subSeismoSolFiles[p]);
			FaultSystemRupSet rupSetWithSubSeis = solWithSubSeis.getRupSet();
			
			Map<String, SummedMagFreqDist> subParentMFDs = new HashMap<>();
			for (String parentName : supraParentMFDs.keySet())
				subParentMFDs.put(parentName, new SummedMagFreqDist(minMag, maxMag, numMag));
			
			for (int s=0; s<rupSetWithSubSeis.getNumSections(); s++) {
				String parentName = rupSetWithSubSeis.getFaultSectionData(s).getParentSectionName();
				IncrementalMagFreqDist subSeisMFD = solWithSubSeis.getSubSeismoOnFaultMFD_List().get(s);
//				if (s == 0)
//					System.out.println(subSeisMFD);
				subParentMFDs.get(parentName).addIncrementalMagFreqDist(subSeisMFD);
			}
			
			List<String> parentNames = new ArrayList<>(parentIDsMap.keySet());
			Collections.sort(parentNames);
			
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>();
			header.add("Parent Section Name");
			EvenlyDiscretizedFunc testMFD = supraParentMFDs.values().iterator().next().getCumRateDistWithOffset();
			for (Point2D pt : testMFD) {
				if (pt.getX() < 0.00001)
					header.add("0.0");
				else
					header.add((float)pt.getX()+"");
			}
			csv.addLine(header);
			for (String parentName : parentNames) {
				IncrementalMagFreqDist supra = supraParentMFDs.get(parentName);
				SummedMagFreqDist sub = subParentMFDs.get(parentName);
				SummedMagFreqDist total = new SummedMagFreqDist(minMag, maxMag, numMag);
				total.addIncrementalMagFreqDist(supra);
				total.addIncrementalMagFreqDist(sub);
				EvenlyDiscretizedFunc cumulative = total.getCumRateDistWithOffset();
				List<String> line = new ArrayList<>();
				line.add(parentName);
				for (Point2D pt : cumulative)
					line.add((float)pt.getY()+"");
				csv.addLine(line);
			}
			csv.writeToFile(new File(outputDir, prefixes[p]+"_sub_plus_supra_parent_cumulative_mfds.csv"));
		}
	}

}
