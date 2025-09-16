package scratch.kevin.nshm23.timeDependence;

import java.io.File;
import java.io.IOException;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.MergedSolutionCreator;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoMappingAlgorithm;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class DOLE_DemosForNed {

	public static void main(String[] args) throws IOException {
		// load the fault system solution
		// this is NSHM23 ba, download from:
		// https://data.opensha.org/ftp/kmilner/markdown/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip"));
		
		// if you need to merge multiple solutions, you can use this util:
		// it doesn't yet do anything with gridded seismicity
//		sol = MergedSolutionCreator.merge(sol1, sol2);
		
		// load DOLE data
		DOLE_SubsectionMapper.mapDOLE(sol.getRupSet().getFaultSectionDataList(), PaleoMappingAlgorithm.FULL_PARENT, true); // boolean is for verbose
		
		// ERF
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		
		// prob model
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.setParameter(HistoricOpenIntervalParam.NAME, 2024d-1875d); // or whatever
		erf.getTimeSpan().setStartTime(2024);
		erf.getTimeSpan().setDuration(50d);
		
		erf.updateForecast();
	}

}
