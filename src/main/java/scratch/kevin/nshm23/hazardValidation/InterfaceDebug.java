package scratch.kevin.nshm23.hazardValidation;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.opensha.nshmp.shaded.model.NshmErf;
import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.util.TectonicRegionType;

public class InterfaceDebug {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		File inputSolFile = new File("/data/kevin/fss_inversions/2026_10-nshm23-td_erf-solutions/nshm23-cascadia-all-hist_dole.zip");
//		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		NshmpHazardModel model = NshmpHazardModel.load(modelDir.toPath());
		
		NshmErf wrapper = new NshmErf(model, Set.of(TectonicRegionType.SUBDUCTION_INTERFACE), IncludeBackgroundOption.EXCLUDE);
		wrapper.getTimeSpan().setDuration(1d);
		wrapper.updateForecast();
		
		System.out.println("FSS has "+rupSet.getNumRuptures()+" ruptures with totalRate="
				+(float)sol.getTotalRateForAllFaultSystemRups()+" /yr");
		int numWrapperRups = 0;
		double totalWrapperRate = 0d;
		for (ProbEqkSource source : wrapper) {
			numWrapperRups += source.getNumRuptures();
			totalWrapperRate += source.computeTotalEquivMeanAnnualRate(1d);
		}
		System.out.println("Wrapper has "+wrapper.getNumSources()+" sources and "+numWrapperRups
				+" total ruptures with totalRate="+(float)totalWrapperRate+" /yr");
		
//		for (int rupIndex=0; rupIndex
	}

}
