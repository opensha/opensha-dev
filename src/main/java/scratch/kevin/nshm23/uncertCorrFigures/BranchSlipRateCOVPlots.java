package scratch.kevin.nshm23.uncertCorrFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.MathArrays;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
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
//		FaultSystemSolution randBA = FaultSystemSolution.load(randBAFile);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(NSHM23_RegionLoader.loadFullConterminousWUS());
		mapMaker.setFaultSections(corrBA.getRupSet().getFaultSectionDataList());
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = Double.NaN;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = Double.NaN;
		List<? extends FaultSection> avgSubSects = NSHM23_DeformationModels.AVERAGE.build(fm, null);
		List<List<? extends FaultSection>> subSectsList = new ArrayList<>();
		List<Double> dmWeights = new ArrayList<>();
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			double weight = dm.getNodeWeight(null);
			if (weight > 0d) {
				subSectsList.add(dm.build(fm, null));
				dmWeights.add(weight);
			}
		}
		
		int numSamples = 1000;
		
		int numSects = subSectsList.get(0).size();
		
		Random r = new Random(numSects*numSamples);

		double[] avgTargetSlipCOVs = new double[numSects];
		int numOver1 = 0;
		for (int s=0; s<numSects; s++) {
			FaultSection sect = avgSubSects.get(s);
			double targetSlip = sect.getReducedAveSlipRate();
			double sd = sect.getReducedSlipRateStdDev();
			avgTargetSlipCOVs[s] = sd/targetSlip;
			if (avgTargetSlipCOVs[s] >= 1d)
				numOver1++;
		}
		System.out.println(numOver1+"/"+numSects+" ("+(float)(100d*numOver1/numSects)+" %) have target average COV>1");
		
		double[] fullTargetSlipCOVs = new double[numSects];
		numOver1 = 0;
		for (int s=0; s<numSects; s++) {
			double[] samples = new double[numSamples*subSectsList.size()];
			double[] weights = new double[samples.length];
			int index = 0;
			for (int d=0; d<dmWeights.size(); d++) {
				List<? extends FaultSection> subSects = subSectsList.get(d);
				double dmWeight = dmWeights.get(d);
				FaultSection sect = subSects.get(s);
				double mean = sect.getOrigAveSlipRate();
				double sd = sect.getOrigSlipRateStdDev();
				for (int i=0; i<numSamples; i++) {
					double sample = mean + r.nextGaussian()*sd;
					Preconditions.checkState(Double.isFinite(sample), "Bad sample %s for mean=%s and sd=%s", sample, mean, sd);
					samples[index] = sample;
					weights[index++] = dmWeight;
				}
			}
			Preconditions.checkState(index == samples.length);
			weights = MathArrays.normalizeArray(weights, weights.length);
			double mean = new Mean().evaluate(samples, weights);
			double sd = Math.sqrt(new Variance().evaluate(samples, weights));
			fullTargetSlipCOVs[s] = sd/mean;
//			System.out.println("COV["+s+"] = "+(float)sd+" / "+(float)mean+" = "+(float)(targetSlipCOVs[s]));
			if (fullTargetSlipCOVs[s] >= 1d)
				numOver1++;
		}
		System.out.println(numOver1+"/"+numSects+" ("+(float)(100d*numOver1/numSects)+" %) have target full dist COV>1");
		
		double[] corrSolCOVs = new double[numSects];
		CSVFile<String> corrSlipCSV = CSVFile.readFile(new File(corrSolDir, "misc_plots/sol_slip_rate_sd.csv"), true);
		Preconditions.checkState(corrSlipCSV.getNumRows() == numSects+1);
		for (int s=0; s<numSects; s++) {
			double corrSolSlip = corrSlipCSV.getDouble(s+1, 1)*1e3;
			double sd = corrSlipCSV.getDouble(s+1, 2)*1e3;
			corrSolCOVs[s] = sd/corrSolSlip;
//			System.out.println(s+". slip="+(float)corrSolSlips[s]+"\tsd="+(float)sd+"\tcov="+(float)corrSolCOVs[s]);
		}
		
		double[] randSolCOVs = new double[numSects];
		CSVFile<String> randSlipCSV = CSVFile.readFile(new File(randSolDir, "misc_plots/sol_slip_rate_sd.csv"), true);
		Preconditions.checkState(randSlipCSV.getNumRows() == numSects+1);
		for (int s=0; s<numSects; s++) {
			double randSolSlip = randSlipCSV.getDouble(s+1, 1)*1e3;
			double sd = randSlipCSV.getDouble(s+1, 2)*1e3;
			randSolCOVs[s] = sd/randSolSlip;
		}
		
		System.out.println("Avg avg-target COV: "+StatUtils.mean(avgTargetSlipCOVs));
		System.out.println("Avg full-target COV: "+StatUtils.mean(fullTargetSlipCOVs));
		System.out.println("Avg corr model COV: "+StatUtils.mean(corrSolCOVs));
		System.out.println("Avg rand model COV: "+StatUtils.mean(randSolCOVs));
		
		System.out.println("Median avg-target COV: "+DataUtils.median(avgTargetSlipCOVs));
		System.out.println("Median full-target COV: "+DataUtils.median(fullTargetSlipCOVs));
		System.out.println("Median corr model COV: "+DataUtils.median(corrSolCOVs));
		System.out.println("Median rand model COV: "+DataUtils.median(randSolCOVs));
		
		CPT covCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0, 2d);
		covCPT.setNanColor(Color.LIGHT_GRAY);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT maskedPDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(pDiffCPT, 5d, 50d);
		maskedPDiffCPT.setNanColor(Color.LIGHT_GRAY);
		CPT covDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.3d, 0.3d);
		covDiffCPT.setNanColor(Color.LIGHT_GRAY);

		mapMaker.plotSectScalars(avgTargetSlipCOVs, covCPT, "Branch Averaged Target Slip Rate COV");
		mapMaker.plot(outputDir, "avg_target_cov", " ");
		mapMaker.plotSectScalars(fullTargetSlipCOVs, covCPT, "Target Slip Rate COV");
		mapMaker.plot(outputDir, "full_target_cov", " ");
		mapMaker.plotSectScalars(corrSolCOVs, covCPT, "NSHM23 Model Solution Slip Rate COV");
		mapMaker.plot(outputDir, "correlated_sol_cov", " ");
		mapMaker.plotSectScalars(randSolCOVs, covCPT, "Randomly Sampled Model Solution Slip Rate COV");
		mapMaker.plot(outputDir, "rand_sol_cov", " ");

//		mapMaker.plotSectScalars(diff(targetSlipCOVs, corrSolCOVs), covDiffCPT, "Original Solution - Target Slip Rate COV Difference");
//		mapMaker.plot(outputDir, "correlated_sol_cov_diff", " ");
//		mapMaker.plotSectScalars(pDiff(targetSlipCOVs, corrSolCOVs), maskedPDiffCPT, "Original Solution vs Target Slip Rate COV % Change");
//		mapMaker.plot(outputDir, "correlated_sol_cov_pDiff", " ");
//		
//		mapMaker.plotSectScalars(diff(targetSlipCOVs, randSolCOVs), covDiffCPT, "Randomly Sampled Solution - Target Slip Rate COV Difference");
//		mapMaker.plot(outputDir, "rand_sol_cov_diff", " ");
//		mapMaker.plotSectScalars(pDiff(targetSlipCOVs, randSolCOVs), maskedPDiffCPT, "Randomly Sampled Solution vs Target Slip Rate COV % Change");
//		mapMaker.plot(outputDir, "rand_sol_cov_pDiff", " ");
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
