package scratch.kevin.cybershake.etasCalcs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.eq.MagUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.FindEquivUCERF2_Ruptures.FindEquivUCERF2_FM2pt1_Ruptures;
import scratch.UCERF3.utils.FindEquivUCERF2_Ruptures.FindEquivUCERF2_FM3_Ruptures;
import scratch.UCERF3.utils.FindEquivUCERF2_Ruptures.FindEquivUCERF2_Ruptures;

public class MappedU2SolCreator {

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas");
		
//		FaultModels fm = FaultModels.FM2_1;
//		DeformationModels dm = DeformationModels.UCERF2_ALL;
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		InversionFaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(
				LaughTestFilter.getDefault(), 0, fm, dm, ScalingRelationships.AVE_UCERF2,
				SlipAlongRuptureModels.TAPERED, InversionModels.CHAR_CONSTRAINED, SpatialSeisPDF.UCERF2);
		
		FindEquivUCERF2_Ruptures findUCERF2_Rups;
		if (fm == FaultModels.FM2_1) 
			findUCERF2_Rups = new FindEquivUCERF2_FM2pt1_Ruptures(rupSet,
					UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR);
		else
			findUCERF2_Rups = new FindEquivUCERF2_FM3_Ruptures(rupSet,
					UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, fm);
		
		ArrayList<double[]> ucerf2_magsAndRates = findUCERF2_Rups.getMagsAndRatesForRuptures();
		
		CSVFile<String> mapping = new CSVFile<String>(true);

		double[] mags = new double[ucerf2_magsAndRates.size()];
		double[] rates = new double[ucerf2_magsAndRates.size()];
		for (int i=0; i<ucerf2_magsAndRates.size(); i++) {
			double[] ucerf2_vals = ucerf2_magsAndRates.get(i);
			if (ucerf2_vals == null) {
				mags[i] = rupSet.getMagForRup(i);
				rates[i] = 0;
			} else {
				mags[i] = ucerf2_vals[0];
				rates[i] = ucerf2_vals[1];
			}
		}
		
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
		mapping.writeToFile(new File(outputDir, "mappings.csv"));
		
//		rupSet = new FaultSystemRupSet(rupSet.getFaultSectionDataList(), rupSet.getSlipRateForAllSections(),
//				rupSet.getSlipRateStdDevForAllSections(), rupSet.getAreaForAllSections(),
//				rupSet.getSectionIndicesForAllRups(), mags, rupSet.getAveRakeForAllRups(),
//				rupSet.getAreaForAllRups(), rupSet.getLengthForAllRups(), "Mapped UCERF2 for CyberShake");
//		
		InversionFaultSystemRupSet modRupSet = new InversionFaultSystemRupSet(rupSet, rupSet.getLogicTreeBranch(),
				rupSet.getLaughTestFilter(), new double[rupSet.getNumRuptures()], rupSet.getCloseSectionsListList(),
				rupSet.getRupturesForClusters(), rupSet.getSectionsForClusters());

		modRupSet.setMagForallRups(mags);

		InversionFaultSystemSolution sol = new InversionFaultSystemSolution(modRupSet, rates, null, null);
		
//		FaultSystemSolution sol = new FaultSystemSolution(rupSet, rates);
		
		// now load in UCERF3 mean sol just for grid sources
		FaultSystemSolution u3Sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/"
				+ "UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		sol.setGridSourceProvider(u3Sol.getGridSourceProvider());
		
		FaultSystemIO.writeSol(sol, new File(outputDir, "ucerf2_mapped_sol.zip"));
	}

}
