package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.inversion.InversionInputGenerator;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoProbabilityModel;

public class RupJumpsTableGen {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
//		File solFile = new File("/tmp/FM3_1_NEOK_EllB_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip");
//		File solFile = new File("/tmp/FM3_1_NEOK_EllB_DsrTap_GRConst_M5Rate8.7_MMaxOff7.6_ApplyCC_SpatSeisU3_sol.zip");
		File solFile = new File("/tmp/FM3_1_GEOL_EllB_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_VarNoSubseismo_VarNone_VarAPrioriWt100_sol.zip");
		InversionFaultSystemSolution sol = FaultSystemIO.loadInvSol(solFile);
		InversionFaultSystemRupSet rupSet = sol.getRupSet();
		
		DeformationModelFetcher dmFetch = new DeformationModelFetcher(
				rupSet.getFaultModel(), rupSet.getFaultModel().getFilterBasis(),
				UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, InversionFaultSystemRupSetFactory.DEFAULT_ASEIS_VALUE);
		
		Map<IDPairing, Double> dists = dmFetch.getSubSectionDistanceMap(5d);
		
		String prefix = solFile.getName().replaceAll(".zip", "");
		PaleoProbabilityModel paleoProbModel = InversionInputGenerator.loadDefaultPaleoProbabilityModel();
		CommandLineInversionRunner.writeJumpPlot(sol, dists, solFile.getParentFile(), prefix, 1d, 7d, null);
		CommandLineInversionRunner.writeJumpPlot(sol, dists, solFile.getParentFile(), prefix, 1d, 0d, paleoProbModel);
		CommandLineInversionRunner.writeJumpPlot(sol, dists, solFile.getParentFile(), prefix, 1d, 0d, null);
		System.exit(0);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		csv.addLine("RupID", "Mag", "NumJumps", "NumJumps>1KM");
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			List<Integer> sects = rupSet.getSectionsIndicesForRup(r);
			
			int jumps = 0;
			int jumpsOver1 = 0;
			for (int i=1; i<sects.size(); i++) {
				int sect1 = sects.get(i-1);
				int sect2 = sects.get(i);
				
				int parent1 = rupSet.getFaultSectionData(sect1).getParentSectionId();
				int parent2 = rupSet.getFaultSectionData(sect2).getParentSectionId();
				
				if (parent1 != parent2) {
					jumps++;
					double dist = dists.get(new IDPairing(sect1, sect2));
					if (dist > 1)
						jumpsOver1++;
				}
			}
			
			csv.addLine(r+"", rupSet.getMagForRup(r)+"", jumps+"", jumpsOver1+"");
		}
		
		File csvFile = new File("/tmp/rup_jumps.csv");
		File txtFile = new File("/tmp/rup_jumps.txt");
		
		csv.writeToFile(csvFile);
		csv.writeToTabSeparatedFile(txtFile, 1);
	}

}
