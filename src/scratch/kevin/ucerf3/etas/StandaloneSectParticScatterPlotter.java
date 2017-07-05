package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.util.List;

import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.utils.FaultSystemIO;

public class StandaloneSectParticScatterPlotter {

	public static void main(String[] args) {
		if (args.length != 4) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(StandaloneSectParticScatterPlotter.class)
				+" <fss-file> <cat-bin-file> <duration> <output-dir>");
			System.exit(2);
		}
		
		File fssFile = new File(args[0]);
		File catalogsFile = new File(args[1]);
		double duration = Double.parseDouble(args[2]);
		File outputDir = new File(args[3]);
		
		try {
			FaultSystemSolution fss = FaultSystemIO.loadSol(fssFile);
			
			System.out.println("Creating ERF for comparisons");
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fss);
			erf = new FaultSystemSolutionERF(fss);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
			erf.getTimeSpan().setDuration(1d);
			erf.updateForecast();
			
			Iterable<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.getBinaryCatalogsIterable(catalogsFile, 0d);
			
			ETAS_MultiSimAnalysisTools.plotSectParticScatter(catalogs, duration, erf, outputDir);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}

}
