package scratch.kevin.mfdInversion;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import scratch.kevin.nshm23.figures.PaleoZTablesBuilder;
import scratch.kevin.nshm23.figures.SlipZTablesBuilder;
import scratch.kevin.nshm23.figures.SlipZTablesBuilder.SlipZRecord;

public class PaleoZComparisons {

	public static void main(String[] args) throws IOException {
		List<File> solFiles = new ArrayList<>();
		List<String> names = new ArrayList<>();
		
		Region[] regions = {
				NSHM23_RegionLoader.AnalysisRegions.CONUS_U3_RELM.load(),
				NSHM23_SingleStates.UT.loadRegion(),
		};
		String[] regionNames = {
				"CA",
				"UT"
		};
		File outputFile = new File("/tmp/posterior_b_z_scores.csv");
		
		File invDir = new File("/home/kevin/OpenSHA/fss_inversions");
		
		solFiles.add(new File(invDir, "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("NSHM23 Full BA");
		
		solFiles.add(new File(invDir, "2024_02_02-nshm23_branches-WUS_FM_v3/node_branch_averaged/PaleoUncert_EvenFitPaleo.zip"));
		names.add("NSHM23 Even-Fit Paleo");
		
		solFiles.add(new File(invDir, "2026_08_10-nshm23_branches-WUS_FM_v3-bPosterior10x/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("Posterior-b Full BA");
		
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		header.add("Model");
		for (int r=0; r<regions.length; r++) {
			header.add(regionNames[r]);
			header.add("Average slip z");
			header.add("Mo-weighted slip z");
			header.add("Paleo z");
			header.add("Paleo |z|");
		}
		csv.addLine(header);
		
		DecimalFormat zDF = new DecimalFormat("0.000");
		
		for (int s=0; s<solFiles.size(); s++) {
			FaultSystemSolution sol = FaultSystemSolution.load(solFiles.get(s));
			List<String> line = new ArrayList<>();
			line.add(names.get(s));
			
			for (int r=0; r<regions.length; r++) {
				line.add("");
				SlipZRecord slips = SlipZTablesBuilder.calcSlipZ(sol, NSHM23_FaultModels.WUS_FM_v3, NSHM23_DeformationModels.AVERAGE, true, regions[r]);
				line.add(zDF.format(slips.average()));
				line.add(zDF.format(slips.momentWeightedAverage()));
				double paleoZ = PaleoZTablesBuilder.avgZScore(sol, false, regions[r]);
				double paleoAbsZ = PaleoZTablesBuilder.avgZScore(sol, true, regions[r]);
				line.add(zDF.format(paleoZ));
				line.add(zDF.format(paleoAbsZ));
			}
			csv.addLine(line);
		}
		csv.writeToFile(outputFile);
	}

}
