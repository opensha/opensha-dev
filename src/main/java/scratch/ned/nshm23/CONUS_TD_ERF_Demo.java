package scratch.ned.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.MergedSolutionCreator;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.DOLE_MappingAlgorithm;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.ned.nshm23.CEUS_FSS_creator.FaultModelEnum;

public class CONUS_TD_ERF_Demo {

	public static void main(String[] args) throws IOException {
		
		
		// load the fault system solution
		// this is NSHM23 ba, download from:
		// https://data.opensha.org/ftp/kmilner/markdown/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip
		FaultSystemSolution conusSol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip"));
		// load DOLE data
		DOLE_SubsectionMapper.mapDOLE(conusSol.getRupSet().getFaultSectionDataList(), DOLE_MappingAlgorithm.FULL_PARENT, true); // boolean is for verbose
		
		
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.b.4/";
		ArrayList<FaultSystemSolution>  ceusSolList = CEUS_FSS_creator.getFaultSystemSolutionList(nshmModelDirPath,CEUS_FSS_creator.FaultModelEnum.PREFERRED);

		
		FaultSystemSolution cascadiaSol = scratch.ned.nshm23.Cascadia_FSS_creator.getFaultSystemSolution(nshmModelDirPath, Cascadia_FSS_creator.FaultModelEnum.MIDDLE);
		
		// if you need to merge multiple solutions, you can use this util:
		// it doesn't yet do anything with gridded seismicity
		FaultSystemSolution sol = MergedSolutionCreator.merge(conusSol, ceusSolList.get(0), ceusSolList.get(1), cascadiaSol);
//		FaultSystemSolution sol = MergedSolutionCreator.merge(conusSol);
		
		System.out.println("DONE with sol");

		// ERF
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		
		// prob model
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.setParameter(HistoricOpenIntervalParam.NAME, 2024d-1875d); // or whatever
		erf.getTimeSpan().setStartTime(2024);
		erf.getTimeSpan().setDuration(50d);
		
		erf.updateForecast();
		
		System.out.println("DONE");
	}

}
