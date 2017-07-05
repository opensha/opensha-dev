package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;


public class ParentSectMultiFaultRupsTableGen {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File file = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/" +
				"2012_10_14-fm3-logic-tree-sample-x5_MEAN_BRANCH_AVG_SOL.zip");
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(file);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		int targetParent = 638;
		
		List<Integer> rups = sol.getRupSet().getRupturesForParentSection(targetParent);
		
		List<String> header = Lists.newArrayList("Rupture Index", "Parent Sections", "Rate", "Magnitude");
		csv.addLine(header);
		
		for (int r : rups) {
			HashSet<Integer> parentsForRup = new HashSet<Integer>();
			List<String> parentNames = Lists.newArrayList();
			
			for (FaultSectionPrefData sect : sol.getRupSet().getFaultSectionDataForRupture(r)) {
				Integer parentID = sect.getParentSectionId();
				if (!parentsForRup.contains(parentID)) {
					parentsForRup.add(parentID);
					parentNames.add(sect.getParentSectionName());
				}
			}
			
			String parents = "["+Joiner.on("],[").join(parentNames)+"]";
			
			List<String> line = Lists.newArrayList();
			
			line.add(r+"");
			line.add(parents);
			line.add(sol.getRateForRup(r)+"");
			line.add(sol.getRupSet().getMagForRup(r)+"");
			
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/hayward_rups.csv"));
	}

}
