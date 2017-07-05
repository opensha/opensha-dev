package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.FaultSystemIO;

public class SAFRupFileWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		File inDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/");
		FaultModels[] fms = {FaultModels.FM3_1, FaultModels.FM3_2};
		
		boolean onlySAF = false;
		
		String onlyStr;
		if (onlySAF)
			onlyStr = "_only";
		else
			onlyStr = "";
		
		HashSet<Integer> ssafSects = new HashSet<Integer>(Lists.newArrayList(32,285,300,287,286,301,282,283,284,295));
		
		for (FaultModels fm : fms) {
			File solFile = new File(inDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_"
					+fm.name()+"_MEAN_BRANCH_AVG_SOL.zip");
			FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
			FaultSystemRupSet rupSet = sol.getRupSet();
			File outFile = new File("/tmp/"+fm.name()+"_ssaf"+onlyStr+"_rups.txt");
			
			HashSet<Integer> allParents = new HashSet<Integer>(fm.getNamedFaultsMapAlt().get("San Andreas"));
			
			HashSet<Integer> rups = new HashSet<Integer>();
			HashSet<Integer> skippedRups = new HashSet<Integer>();
			
			HashSet<Integer> tempRups = new HashSet<Integer>();
			for (int parentID : ssafSects) {
				List<Integer> parentRups = rupSet.getRupturesForParentSection(parentID);
				if (onlySAF) {
					// filter out ones that include anything off of SSAF
					List<Integer> filtered = Lists.newArrayList();
					for (int rup : parentRups) {
						if (rups.contains(rup) || skippedRups.contains(rup))
							continue;
						boolean skip = false;
						for (FaultSectionPrefData sect : rupSet.getFaultSectionDataForRupture(rup)) {
							if (!allParents.contains(sect.getParentSectionId())) {
								skip = true;
								break;
							}
						}
						if (skip)
							skippedRups.add(rup);
						else
							filtered.add(rup);
					}
					parentRups = filtered;
				}
				if (parentID == 301)
					tempRups.addAll(parentRups);
				rups.addAll(parentRups);
			}
			
			double rate7p8 = 0;
			for (int rupID : rupSet.getRupturesForParentSection(301))
				if (rups.contains(rupID) && rupSet.getMagForRup(rupID) >= 7.8)
					rate7p8 += sol.getRateForRup(rupID);
//			for (int rupID : tempRups) {
//				if (rupSet.getMagForRup(rupID) >= 7.8)
//					rate7p8 += sol.getRateForRup(rupID);
//			}
			System.out.println(fm.name()+" Mojave 7.8 rate: "+rate7p8+", 30yr Poisson: "+(1d-Math.exp(-rate7p8*30d)));
			
			List<Integer> rupsSorted = Lists.newArrayList(rups);
			Collections.sort(rupsSorted);
			
			FileWriter fw = new FileWriter(outFile);
			
			fw.write("#ID\tMag\tRate\tSubsection IDs");
			for (int rup : rupsSorted) {
				StringBuilder str = new StringBuilder();
				str.append(rup).append("\t").append(rupSet.getMagForRup(rup)).append("\t");
				str.append(sol.getRateForRup(rup));
				for (int sectIndex : rupSet.getSectionsIndicesForRup(rup))
					str.append("\t").append(sectIndex);
				str.append("\n");
				fw.write(str.toString());
			}
			
			fw.close();
		}
	}

}
