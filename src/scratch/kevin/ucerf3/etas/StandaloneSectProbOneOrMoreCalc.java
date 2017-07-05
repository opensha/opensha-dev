package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.util.List;

import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.utils.FaultSystemIO;

public class StandaloneSectProbOneOrMoreCalc {

	public static void main(String[] args) {
		if (args.length != 5) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(StandaloneSectProbOneOrMoreCalc.class)
				+" <fss-file> <cat-bin-file> <start-year> <prob-duration> <output-dir>");
			System.exit(2);
		}
		
		File fssFile = new File(args[0]);
		File catalogsFile = new File(args[1]);
		int startYear = Integer.parseInt(args[2]);
		double durationForProb = Double.parseDouble(args[3]);
		File outputDir = new File(args[4]);
		
		try {
			FaultSystemSolution fss = FaultSystemIO.loadSol(fssFile);
			
			System.out.println("Creating ERF for comparisons");
//			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fss);
//			erf = new FaultSystemSolutionERF(fss);
//			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
//			erf.setParameter(MagDependentAperiodicityParam.NAME, MagDependentAperiodicityOptions.MID_VALUES);
//			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
//			erf.getTimeSpan().setDuration(1d);
//			erf.updateForecast();
			FaultSystemSolutionERF erf = MPJ_ETAS_Simulator.buildERF(fss, false, durationForProb, startYear);
			erf.updateForecast();
			
			Iterable<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.getBinaryCatalogsIterable(catalogsFile, 0d);
			
//			ETAS_MultiSimAnalysisTools.plotSectParticScatter(catalogs, duration, erf, outputDir);
			ETAS_MultiSimAnalysisTools.plotAndWriteSectProbOneOrMoreData(catalogs, durationForProb, erf, outputDir);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}

}
