package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.utils.U3FaultSystemIO;

public class BombayBeachSwarm2016Calcs {

	public static void main(String[] args) throws IOException, DocumentException {
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_09_26-2016_bombay_swarm-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont");
//		long ot = 1474920000000l;
//		File binFile = new File(dir, "results.bin");
		
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_09_29-2016_bombay_swarm-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont");
//				+ "2016_10_25-2016_bombay_swarm-24yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont-combined");
		File dir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/sims/"
				+ "2017_01_27-2016_bombay_swarm-10yr-u2mapped-full_td-gridSeisCorr-noSpont");
		long ot = 1474990200000l;
		File binFile = new File(dir, "results.bin");
//		File binFile = new File(dir, "results_300k.bin");
		
		
		long oneWeekOT = ot + 7*ProbabilityModelsCalc.MILLISEC_PER_DAY;
		
		List<? extends List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(binFile);
		
		File plotDir = new File(dir, "plots");
		Preconditions.checkState(plotDir.exists() || plotDir.mkdir());
		
		String name = "2016 Bombay Swarm";
		
		ETAS_MultiSimAnalysisTools.mfdMaxY = 1e2;
		ETAS_MultiSimAnalysisTools.mfdMinY = 1e-6;
		
		ETAS_MultiSimAnalysisTools.plotMagNum(catalogs, plotDir, name, "10_year", null, 0, null);
		
		ETAS_MultiSimAnalysisTools.plotMagNum(catalogs, plotDir, name, "one_week", null, 0, null, oneWeekOT);
		
		File fssFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution fss = U3FaultSystemIO.loadSol(fssFile);
		
		// use the old bombay beach scenario here, just used for finding nearby sections so it's fine to use
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fss);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		boolean rates = false;
		boolean subSects = false;
		ETAS_MultiSimAnalysisTools.plotScalesOfHazardChange(catalogs, null, TestScenario.BOMBAY_BEACH_M4pt8,
				ot, erf, plotDir, name, 10d, rates, subSects);
		
//		ETAS_MultiSimAnalysisTools.plotScalesOfHazardChange(catalogs, null, TestScenario.BOMBAY_BEACH_M4pt8, ot,
//				erf, outputDir, name, inputDuration, rates, subSects);
		
		/*
		
		double[] radii = { 200d, 150d, 100d, 50d, 25d };
		
		double duration = 7d/365.25;
		
		FaultSystemSolutionERF tdERF = MPJ_ETAS_Simulator.buildERF_millis(fss, false, duration, ot);
		tdERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		double startYear = 1970d + (double)ot/(double)ProbabilityModelsCalc.MILLISEC_PER_YEAR;
		tdERF.getParameter(HistoricOpenIntervalParam.NAME).setValue(startYear-1875d);
		tdERF.getTimeSpan().setStartTimeInMillis(ot+1);
		tdERF.getTimeSpan().setDuration(duration);
		tdERF.updateForecast();
		FaultSystemSolutionERF tiERF = MPJ_ETAS_Simulator.buildERF_millis(fss, true, duration, ot);
		tiERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		tiERF.getTimeSpan().setDuration(duration);
		tiERF.updateForecast();
		
		System.out.println("*** One Week Circular Region Nucleation Rates ***");
		
		for (double radius : radii) {
			Region region = new Region(new Location(33.298, -115.713), radius);
			
			IncrementalMagFreqDist tdNuclMFD = ERF_Calculator.getMagFreqDistInRegionFaster(tdERF, region, 4.05, 51, 0.1, true);
			IncrementalMagFreqDist tiNuclMFD = ERF_Calculator.getMagFreqDistInRegionFaster(tiERF, region, 4.05, 51, 0.1, true);
			// scale to 1 week
			tdNuclMFD.scale(duration);
			tiNuclMFD.scale(duration);
			
			EvenlyDiscretizedFunc tdNuclCmlMFD = tdNuclMFD.getCumRateDistWithOffset();
			EvenlyDiscretizedFunc tiNuclCmlMFD = tiNuclMFD.getCumRateDistWithOffset();
			
//			IncrementalMagFreqDist nuclMFD = fss.calcNucleationMFD_forRegion(region, 4.05, 9.05, 0.1, true);
//			// scale to 1 week
//			nuclMFD.scale(7d/365.25);
//			EvenlyDiscretizedFunc nuclCmlMFD = nuclMFD.getCumRateDistWithOffset();
			
			System.out.println("Radius: "+(float)radius+" km");
			System.out.println("\tUCERF3-TD:\tM>=4.3:\t"+(float)tdNuclCmlMFD.getY(4.3d)+"\tM>=7:\t"+(float)tdNuclCmlMFD.getY(7d));
			System.out.println("\tUCERF3-TI:\tM>=4.3:\t"+(float)tiNuclCmlMFD.getY(4.3d)+"\tM>=7:\t"+(float)tiNuclCmlMFD.getY(7d));
			
//			System.out.println("1 week Cumulative Nucleation MFD:");
//			System.out.println(nuclCmlMFD);
//			double m7 = nuclCmlMFD.getY(7d);
//			
//			System.out.println((float)radius+"km: "+(float)m7);
		}
		
		int parentID = 295;
		double tdCoachPartic = FaultSysSolutionERF_Calc.calcParticipationRateForParentSect(tdERF, parentID, 7d)*duration;
		double tdCoachNucl = FaultSysSolutionERF_Calc.calcNucleationRateForParentSect(tdERF, parentID, 7d)*duration;
		
		double tiCoachPartic = FaultSysSolutionERF_Calc.calcParticipationRateForParentSect(tiERF, parentID, 7d)*duration;
		double tiCoachNucl = FaultSysSolutionERF_Calc.calcNucleationRateForParentSect(tiERF, parentID, 7d)*duration;
		
		System.out.println("*** Coachella 1 Week M>=7 Rates ***");
		System.out.println("UCERF3-TD: partic: "+(float)tdCoachPartic+"\tnucl: "+(float)tdCoachNucl);
		System.out.println("UCERF3-TI: partic: "+(float)tiCoachPartic+"\tnucl: "+(float)tiCoachNucl);
		
		*/
	}

}
