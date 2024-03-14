package scratch.kevin.nshm23.uncertCorrFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.nshm23.figures.MethodsAndIngredientsHazChangeFigures;

public class BranchSlipRateCOVPlots {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/data/kevin/nshm23/batch_inversions/");
		
		File corrSolDir = new File(invsDir, "2024_02_02-nshm23_branches-WUS_FM_v3");
		File corrBAFile = new File(corrSolDir, "results_WUS_FM_v3_branch_averaged.zip");
		
		File randSolDir = new File(invsDir, "2023_11_17-nshm23_branches-dm_sampling-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File randBAFile = new File(randSolDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
		
		File outputDir = new File("/home/kevin/Documents/papers/2024_nshm23_uncert_correlation/figures/branch_slip_rate_covs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution corrBA = FaultSystemSolution.load(corrBAFile);
		FaultSystemSolution randBA = FaultSystemSolution.load(randBAFile);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(NSHM23_RegionLoader.loadFullConterminousWUS());
		mapMaker.setFaultSections(corrBA.getRupSet().getFaultSectionDataList());
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = Double.NaN;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = Double.NaN;
		List<? extends FaultSection> avgSects = NSHM23_DeformationModels.AVERAGE.build(fm);
		
		double[] targetSlips = new double[avgSects.size()];
		double[] targetSlipCOVs = new double[avgSects.size()];
		
		int numOver1 = 0;
		
		for (int s=0; s<avgSects.size(); s++) {
			FaultSection sect = avgSects.get(s);
			targetSlips[s] = sect.getReducedAveSlipRate();
			double sd = sect.getReducedSlipRateStdDev();
			targetSlipCOVs[s] = sd/targetSlips[s];
			if (targetSlipCOVs[s] >= 1d)
				numOver1++;
		}
		System.out.println(numOver1+"/"+avgSects.size()+" ("+(float)(100d*numOver1/avgSects.size())+" %) have target COV>1");
		
		double[] corrSolSlips = new double[avgSects.size()];
		double[] corrSolCOVs = new double[avgSects.size()];
		CSVFile<String> corrSlipCSV = CSVFile.readFile(new File(corrSolDir, "misc_plots/sol_slip_rate_sd.csv"), true);
		Preconditions.checkState(corrSlipCSV.getNumRows() == avgSects.size()+1);
		for (int s=0; s<avgSects.size(); s++) {
			corrSolSlips[s] = corrSlipCSV.getDouble(s+1, 1)*1e3;
			double sd = corrSlipCSV.getDouble(s+1, 2)*1e3;
			corrSolCOVs[s] = sd/corrSolSlips[s];
//			System.out.println(s+". slip="+(float)corrSolSlips[s]+"\tsd="+(float)sd+"\tcov="+(float)corrSolCOVs[s]);
		}
		
		double[] randSolSlips = new double[avgSects.size()];
		double[] randSolCOVs = new double[avgSects.size()];
		CSVFile<String> randSlipCSV = CSVFile.readFile(new File(randSolDir, "misc_plots/sol_slip_rate_sd.csv"), true);
		Preconditions.checkState(randSlipCSV.getNumRows() == avgSects.size()+1);
		for (int s=0; s<avgSects.size(); s++) {
			randSolSlips[s] = randSlipCSV.getDouble(s+1, 1)*1e3;
			double sd = randSlipCSV.getDouble(s+1, 2)*1e3;
			randSolCOVs[s] = sd/randSolSlips[s];
		}
		
		CPT covCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0, 2d);
		covCPT.setNanColor(Color.LIGHT_GRAY);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT maskedPDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(pDiffCPT, 5d, 50d);
		maskedPDiffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT covDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.3d, 0.3d);
		covDiffCPT.setNanColor(Color.LIGHT_GRAY);
		
		mapMaker.plotSectScalars(targetSlipCOVs, covCPT, "Branch Averaged Target Slip Rate COV");
		mapMaker.plot(outputDir, "target_cov", " ");
		mapMaker.plotSectScalars(corrSolCOVs, covCPT, "Original Model Solution Slip Rate COV");
		mapMaker.plot(outputDir, "correlated_sol_cov", " ");
		mapMaker.plotSectScalars(randSolCOVs, covCPT, "Randomly Sampled Model Solution Slip Rate COV");
		mapMaker.plot(outputDir, "rand_sol_cov", " ");

		mapMaker.plotSectScalars(diff(targetSlipCOVs, corrSolCOVs), covDiffCPT, "Original Solution - Target Slip Rate COV Difference");
		mapMaker.plot(outputDir, "correlated_sol_cov_diff", " ");
		mapMaker.plotSectScalars(pDiff(targetSlipCOVs, corrSolCOVs), maskedPDiffCPT, "Original Solution vs Target Slip Rate COV % Change");
		mapMaker.plot(outputDir, "correlated_sol_cov_pDiff", " ");
		
		mapMaker.plotSectScalars(diff(targetSlipCOVs, randSolCOVs), covDiffCPT, "Randomly Sampled Solution - Target Slip Rate COV Difference");
		mapMaker.plot(outputDir, "rand_sol_cov_diff", " ");
		mapMaker.plotSectScalars(pDiff(targetSlipCOVs, randSolCOVs), maskedPDiffCPT, "Randomly Sampled Solution vs Target Slip Rate COV % Change");
		mapMaker.plot(outputDir, "rand_sol_cov_pDiff", " ");
	}
	
	private static double[] diff(double[] ref, double[] comp) {
		double[] ret = new double[ref.length];
		for (int i=0; i<ret.length; i++)
			ret[i] = comp[i] - ref[i];
		return ret;
	}
	
	private static double[] pDiff(double[] ref, double[] comp) {
		double[] ret = new double[ref.length];
		for (int i=0; i<ret.length; i++)
			ret[i] = 100d*(comp[i] - ref[i])/ref[i];
		return ret;
	}

}
