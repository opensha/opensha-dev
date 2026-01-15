package scratch.kevin.tdProbModelPlayground;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.AperiodicityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.HistoricalOpenIntervals;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.RenewalModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.TimeDepFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.UCERF3_ProbabilityModel;
import org.opensha.sha.earthquake.param.AseismicityAreaReductionParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;

import com.google.common.base.Stopwatch;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class TD_Benchmark {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = TD_ERF_Example.fetchU3_BA();
		FaultSystemSolutionERF origERF = new FaultSystemSolutionERF(sol);
		TimeDepFaultSystemSolutionERF newERF = new TimeDepFaultSystemSolutionERF(sol);
		
		int spinUpTrials = 10;
		int trials = 30;
		
		double[] times = benchmark(origERF, true, spinUpTrials, trials);
		double origPoisson = DataUtils.median(times);
		times = benchmark(newERF, true, spinUpTrials, trials);
		double newPoisson = DataUtils.median(times);
		times = benchmark(origERF, false, spinUpTrials, trials);
		double origTD = DataUtils.median(times);
		times = benchmark(newERF, false, spinUpTrials, trials);
		double newTD = DataUtils.median(times);
		
		double tdSpeedup = origTD / newTD;

		double origTDExtra = origTD - origPoisson;
		double newTDExtra = newTD - newPoisson;

		// Guard against negative/zero “extra” times (can happen due to noise/caching)
		double tdExtraSpeedup = (newTDExtra > 0d && origTDExtra > 0d) ? (origTDExtra / newTDExtra) : Double.NaN;

		System.out.println();
		System.out.println("==== Median updateForecast() times (s) ====");
		System.out.println("Original ERF: Poisson = " + fmt(origPoisson) + ", TD = " + fmt(origTD));
		System.out.println("New ERF:      Poisson = " + fmt(newPoisson) + ", TD = " + fmt(newTD));
		System.out.println();

		System.out.println("==== Speedups (old/new) ====");
		System.out.println("TD total speedup: " + fmt(tdSpeedup) + "x");

		System.out.println();
		System.out.println("==== TD extra over Poisson (s) ====");
		System.out.println("Original TD extra: " + fmt(origTDExtra) + " (=" + fmt(origTD) + " - " + fmt(origPoisson) + ")");
		System.out.println("New TD extra:      " + fmt(newTDExtra) + " (=" + fmt(newTD) + " - " + fmt(newPoisson) + ")");

		System.out.println();
		if (Double.isNaN(tdExtraSpeedup)) {
			System.out.println("TD extra speedup: NaN (one of the extra times was <= 0; likely measurement noise/caching effects)");
		} else {
			System.out.println("TD extra speedup (old/new): " + fmt(tdExtraSpeedup) + "x");
		}
	}
	
	private static String fmt(double val) {
		if (Double.isNaN(val))
			return "NaN";
		if (Double.isInfinite(val))
			return val > 0 ? "Inf" : "-Inf";
		// 3 decimals is usually enough for seconds and ratios here
		return String.format("%.3f", val);
	}
	
	static void parameterizeERF(BaseFaultSystemSolutionERF erf, boolean poisson) {
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		if (erf.getAdjustableParameterList().containsParameter(UseRupMFDsParam.NAME))
			erf.setParameter(UseRupMFDsParam.NAME, false);
		
		if (erf instanceof FaultSystemSolutionERF) {
			if (poisson) {
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			} else {
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
				erf.setParameter(MagDependentAperiodicityParam.NAME, MagDependentAperiodicityOptions.MID_VALUES);
				erf.setParameter(HistoricOpenIntervalParam.NAME, 2014d-1875d);
				erf.setParameter(BPTAveragingTypeParam.NAME, BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
			}
		} else {
			TimeDepFaultSystemSolutionERF tdERF = (TimeDepFaultSystemSolutionERF)erf;
			if (poisson) {
				tdERF.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
			} else {
				tdERF.setProbabilityModelChoice(FSS_ProbabilityModels.UCERF3_BPT);
				UCERF3_ProbabilityModel probModel = (UCERF3_ProbabilityModel) tdERF.getProbabilityModel();
				probModel.setAperiodicityModelChoice(AperiodicityModels.UCERF3_MIDDLE);
				probModel.setHistOpenIntervalChoice(HistoricalOpenIntervals.UCERF3);
				probModel.setRenewalModelChoice(RenewalModels.BPT);
				probModel.setAveragingTypeChoice(BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
				probModel.setProbDistsDiscretization(9, 18000, false);
				probModel.setIntegrationNormCDFsDiscretization(5d, 501);
			}
		}
	}
	
	private static double[] benchmark(BaseFaultSystemSolutionERF erf, boolean poisson, int spinUpTrials, int trials) {
		parameterizeERF(erf, poisson);
		
		System.out.println("Spin up for "+ClassUtils.getClassNameWithoutPackage(erf.getClass())+", poisson="+poisson);
		runBenchmark(erf, spinUpTrials);
		
		System.out.println("Benchmarking "+ClassUtils.getClassNameWithoutPackage(erf.getClass())+", poisson="+poisson);
		return runBenchmark(erf, trials);
	}
	
	private static double[] runBenchmark(BaseFaultSystemSolutionERF erf, int trials) {
		double[] times = new double[trials];
		for (int t=0; t<trials; t++) {
			// make it update the fss sources
			erf.setParameter(AseismicityAreaReductionParam.NAME, false);
			erf.setParameter(AseismicityAreaReductionParam.NAME, true);
			
			Stopwatch watch = Stopwatch.createStarted();
			erf.updateForecast();
			watch.stop();
			times[t] = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
			System.out.println("\tfinished trial "+t+"/"+trials+" in "+(float)times[t]+" s");
		}
		return times;
	}

}
