package scratch.kevin.cybershake.etasCalcs;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;

import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.FindEquivUCERF2_Ruptures.FindEquivUCERF2_FM2pt1_Ruptures;

public class InvertedU2MappingFileCreator {

	public static void main(String[] args) throws IOException {
		LogicTreeBranch branch = LogicTreeBranch.fromFileName(
				"FM2_1_UC2ALL_AveU2_DsrUni_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU2");
		InversionFaultSystemRupSet origRupSet = InversionFaultSystemRupSetFactory.forBranch(branch);
		InversionFaultSystemRupSet rupSet = CommandLineInversionRunner.getUCERF2RupsOnly(origRupSet);
		
		File outputFile = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/u3_inverted_mappings.csv");
		
		FindEquivUCERF2_FM2pt1_Ruptures findUCERF2_Rups = new FindEquivUCERF2_FM2pt1_Ruptures(rupSet,
				UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR);
		
		CSVFile<String> mapping = new CSVFile<String>(true);
		
		mapping.addLine("UCERF2 Source Index", "UCERF2 Rupture Index", "Source Name", "Mag",
				"FSS Index", "FSS Start Sub", "FSS End Sub");
		ERF u2ERF = findUCERF2_Rups.getUCERF2_ERF();
		
		int r = 0;
		sourceLoop:
		for (int sourceIndex=0; sourceIndex<u2ERF.getNumSources(); sourceIndex++) {
			ProbEqkSource source = u2ERF.getSource(sourceIndex);
			for (int ruptureIndex=0; ruptureIndex<u2ERF.getNumRuptures(sourceIndex); ruptureIndex++) {
				ProbEqkRupture rup = source.getRupture(ruptureIndex);
				int fssIndex = findUCERF2_Rups.getEquivFaultSystemRupIndexForUCERF2_Rupture(r);
				if (fssIndex >= 0) {
					List<FaultSectionPrefData> subs = rupSet.getFaultSectionDataForRupture(fssIndex);
					mapping.addLine(sourceIndex+"", ruptureIndex+"", source.getName(), rup.getMag()+"",
							fssIndex+"", subs.get(0).getName(), subs.get(subs.size()-1).getName());
				}
				
				r++;
				if (r >= findUCERF2_Rups.getNumUCERF2_Ruptures())
					break sourceLoop;
			}
		}
		
		mapping.writeToFile(outputFile);
	}

}
