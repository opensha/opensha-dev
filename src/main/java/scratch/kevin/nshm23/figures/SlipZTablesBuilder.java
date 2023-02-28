package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

public class SlipZTablesBuilder {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		boolean excludeClassic = true;
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		
		String prefix = "slip_misfits";
		File baSolFile;
		File nodeBADir;
		if (excludeClassic) {
			prefix += "_no_classic";
			baSolFile = new File(mainDir, "tmp_ba_no_classic.zip");
			nodeBADir = new File(mainDir, "node_branch_averaged_no_classic");
		} else {
			baSolFile = new File(mainDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
			nodeBADir = new File(mainDir, "node_branch_averaged");
		}
		
		NSHM23_DeformationModels[] dms = {
				NSHM23_DeformationModels.AVERAGE,
				NSHM23_DeformationModels.EVANS,
				NSHM23_DeformationModels.GEOLOGIC,
				NSHM23_DeformationModels.POLLITZ,
				NSHM23_DeformationModels.SHEN_BIRD,
				NSHM23_DeformationModels.ZENG
		};
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;

		double[] origUncertZs = new double[dms.length];
		double[] origUncertMoWeightedZs = new double[dms.length];
		double[] modUncertZs = new double[dms.length];
		double[] modUncertMoWeightedZs = new double[dms.length];
		
		for (int d=0; d<dms.length; d++) {
			NSHM23_DeformationModels dm = dms[d];
			FaultSystemSolution sol;
			if (dm == NSHM23_DeformationModels.AVERAGE) {
				sol = FaultSystemSolution.load(baSolFile);
			} else {
				File nodeFile = new File(nodeBADir, "DM_"+dm.getFilePrefix()+".zip");
				sol = FaultSystemSolution.load(nodeFile);
			}
			
			SectSlipRates sectSlips = sol.getRupSet().requireModule(SectSlipRates.class);
			SolutionSlipRates solSlips = sol.requireModule(SolutionSlipRates.class);
			
			for (boolean origUncert : new boolean[] {false,true}) {
				double[] slipSDs;
				if (origUncert) {
					NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = Double.NaN;
					List<? extends FaultSection> subSects = dm.build(fm);
					slipSDs = new double[subSects.size()];
					for (int s=0; s<slipSDs.length; s++)
						slipSDs[s] = subSects.get(s).getOrigSlipRateStdDev()*1e-3;
				} else {
					slipSDs = sectSlips.getSlipRateStdDevs();
				}
				
				double avgZ = 0d;
				double moWeightedAvgZ = 0d;
				double sumMo = 0d;
				
				MinMaxAveTracker zTrack = new MinMaxAveTracker();
				MinMaxAveTracker sdTrack = new MinMaxAveTracker();
				MinMaxAveTracker covTrack = new MinMaxAveTracker();
				System.out.println(dm.getShortName()+", orig="+origUncert);
				
				double maxZ = 0d;
				int maxIndex = -1;
				
				for (int s=0; s<slipSDs.length; s++) {
					FaultSection sect = sol.getRupSet().getFaultSectionData(s);
					double moRate = sect.calcMomentRate(false);
					
					double z = (solSlips.get(s) - sectSlips.getSlipRate(s))/slipSDs[s];
					z = Math.abs(z);
					
					zTrack.addValue(z);
					sdTrack.addValue(slipSDs[s]);
					covTrack.addValue(slipSDs[s]/sectSlips.getSlipRate(s));
					
					avgZ += z;
					moWeightedAvgZ += moRate*z;
					sumMo += moRate;
					
					if (z > maxZ) {
						maxIndex = s;
						maxZ = z;
					}
				}
				
				avgZ /= slipSDs.length;
				moWeightedAvgZ /= sumMo;
				
				System.out.println("\tabs z-scores: "+zTrack);
				System.out.println("\tsds: "+sdTrack);
				System.out.println("\tcovs: "+covTrack);
				System.out.println("\tmaxZ="+(float)maxZ+" for "+maxIndex+". "
						+sol.getRupSet().getFaultSectionData(maxIndex).getSectionName());
				System.out.println("\t\t("+(float)solSlips.get(maxIndex)+" - "
						+(float)sectSlips.getSlipRate(maxIndex)+") / "+(float)slipSDs[maxIndex]);
				
				if (origUncert ) {
					origUncertZs[d] = avgZ;
					origUncertMoWeightedZs[d] = moWeightedAvgZ;
				} else {
					modUncertZs[d] = avgZ;
					modUncertMoWeightedZs[d] = moWeightedAvgZ;
				}
			}
		}
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("");
		for (NSHM23_DeformationModels dm : dms)
			header.add(dm == null ? "Branch-Averaged" : dm.getShortName());
		csv.addLine(header);

		List<String> names = List.of("Original Uncertainties",
				"Original Uncertainties, Moment-Weighted",
				"10% Uncertainties",
				"10% Uncertainties, Moment-Weighted");
		List<double[]> arrays = List.of(
				origUncertZs, origUncertMoWeightedZs,
				modUncertZs, modUncertMoWeightedZs);
		
		for (int i=0; i<names.size(); i++) {
			List<String> line = new ArrayList<>();
			line.add(names.get(i));
			for (double z : arrays.get(i))
				line.add((float)z+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, prefix+".csv"));
	}

}
