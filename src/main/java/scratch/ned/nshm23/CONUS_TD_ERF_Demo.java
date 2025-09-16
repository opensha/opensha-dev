package scratch.ned.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.GregorianCalendar;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.MergedSolutionCreator;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoMappingAlgorithm;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen.DataToInclude;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.ned.nshm23.AK_FSS_creator.DeformationModelEnum;
import scratch.ned.nshm23.CEUS_FSS_creator.FaultModelEnum;

/**
 * NOTES:
 * 
 * Parent section names are null for AK fault system solution faults (info not readily available, but IDs are good)
 */
public class CONUS_TD_ERF_Demo {
	
	private static FaultSystemSolution getFull_FSS(String full_FSS_fileName) {
		
		
		FaultSystemSolution sol=null;
		try {
			//WUS
			FaultSystemSolution wusSol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip"));
			String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.1.2/";
			// CEUS
			ArrayList<FaultSystemSolution>  ceusSolList = CEUS_FSS_creator.getFaultSystemSolutionList(nshmModelDirPath,CEUS_FSS_creator.FaultModelEnum.PREFERRED);
			// AK
			String akModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-3.0.1/";
			ArrayList<FaultSystemSolution> ak_fssList = AK_FSS_creator.getFaultSystemSolutionList(akModelDirPath, DeformationModelEnum.GEO);
			FaultSystemSolution solAK = MergedSolutionCreator.merge(ak_fssList);
			
//			// DEBUG null parent names
//			int n=0;
//			for(FaultSection sect:solAK.getRupSet().getFaultSectionDataList())
//				if(sect.getParentSectionName() == null) {
//					System.out.println(sect.getParentSectionName()+"\t"+sect.getSectionName());
//					n+=1;
//				}
//			System.out.println(n+" of "+solAK.getRupSet().getFaultSectionDataList().size());
//			System.exit(0);

			
			// full/merged solution
			sol = MergedSolutionCreator.merge(wusSol, ceusSolList.get(0), ceusSolList.get(1), solAK);
		} catch (IOException e) {
			e.printStackTrace();
		}		
		
		try {
			sol.write(new File("/Users/field/nshm-haz_data/full_FSS_test.zip"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return sol;

	}


	public static void main(String[] args) throws IOException {
		
		long currentTimeEpoch = System.currentTimeMillis();
		String dateString = new java.text.SimpleDateFormat("MM_dd_yyyy").format(new java.util.Date (currentTimeEpoch)); // Epoch in seconds, remove '*1000' for milliseconds
//		System.out.println(dateString);
//		System.exit(-1);
		
		File tdMainDir = new File("/Users/field/markdown/nshm23_time_dependence_"+dateString);

		String full_FSS_fileName = "/Users/field/nshm-haz_data/full_FSS_test.zip";
		
		FaultSystemSolution sol;
		File file = new File("/Users/field/nshm-haz_data/full_FSS_test.zip");
		if(file.exists())
			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
		else
			sol = getFull_FSS(full_FSS_fileName);
		
		if(!tdMainDir.exists()) tdMainDir.mkdir();
		

		
		TimeDependentReportPageGen.generatePage(new File(tdMainDir, "allDOLE_fullParent"), sol, PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.ALL_DATA);

		// recreating solution to avoid propagating previous DOLE settings
//		sol = FaultSystemSolution.load(new File(full_FSS_fileName));
		TimeDependentReportPageGen.generatePage(new File(tdMainDir, "allDOLE_neighbors"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.ALL_DATA);

//		sol = FaultSystemSolution.load(new File(full_FSS_fileName));
		TimeDependentReportPageGen.generatePage(new File(tdMainDir, "forDebugging_onlyPaleoDOLE_nearestSubsect"), sol, PaleoMappingAlgorithm.CLOSEST_SECT, DataToInclude.PALEO_ONLY);
		
//		sol = FaultSystemSolution.load(new File(full_FSS_fileName));
		TimeDependentReportPageGen.generatePage(new File(tdMainDir, "onlyHistoricRupDOLE"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.HIST_RUPS_ONLY);

//		sol = FaultSystemSolution.load(new File(full_FSS_fileName));
		TimeDependentReportPageGen.generatePage(new File(tdMainDir, "onlyPaleoDOLE_fullParent"), sol, PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.PALEO_ONLY);

//		sol = FaultSystemSolution.load(new File(full_FSS_fileName));
		TimeDependentReportPageGen.generatePage(new File(tdMainDir, "onlyPaleoDOLE_neighbors"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.PALEO_ONLY);

		
//		TimeDependentReportPageGen.generatePage(new File ("pageGenTestRightHere"), sol, PaleoMappingAlgorithm.CLOSEST_SECT);

//		System.out.println("DONE with sol");
//		String mappingSummary = DOLE_SubsectionMapper.mapDOLE(sol.getRupSet().getFaultSectionDataList(), PaleoMappingAlgorithm.CLOSEST_SECT, false); // boolean is for verbose
//		System.out.println(mappingSummary);
//		System.out.println("DONE");

		
		
		// load the WUS fault system solution
		// this is NSHM23 ba, download from:
		// https://data.opensha.org/ftp/kmilner/markdown/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip
//		FaultSystemSolution wusSol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip"));		
//		System.out.println("DONE with wusSol");
		
		// I NEED TO REVIEW AND CHECK SUBSECTIONING INDICES IN NON-WUS FSSs BELOW BEFORE CONTINUING ********************

//		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.1.2/";
//		ArrayList<FaultSystemSolution>  ceusSolList = CEUS_FSS_creator.getFaultSystemSolutionList(nshmModelDirPath,CEUS_FSS_creator.FaultModelEnum.PREFERRED);
//		System.out.println("DONE with ceusSolList; writing parent fault sections:");
//		for(FaultSystemSolution fss:ceusSolList) {
//			for(FaultSection fs:fss.getRupSet().getFaultSectionDataList()) {
//				System.out.println("\t"+fs.getParentSectionId()+"\t"+fs.getParentSectionName());
//			}
//		}
//		System.exit(0);

//		FaultSystemSolution cascadiaSol = scratch.ned.nshm23.Cascadia_FSS_creator.getFaultSystemSolution(nshmModelDirPath, Cascadia_FSS_creator.FaultModelEnum.MIDDLE);
		
//		// Alaska test
//		String akModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-3.0.1/";
//		ArrayList<FaultSystemSolution> ak_fssList = AK_FSS_creator.getFaultSystemSolutionList(akModelDirPath, DeformationModelEnum.GEO);
//		System.out.println("ak_fssList.size() = "+ak_fssList.size());
//		FaultSystemSolution mrgedAK_Sol = MergedSolutionCreator.merge(ak_fssList);
//		FaultSystemSolution solAK = MergedSolutionCreator.merge(ak_fssList);



		// *************************************************
		
//		FaultSystemSolution sol;
		// this doesn't yet do anything with gridded seismicity
//		sol = wusSol;
//		sol = MergedSolutionCreator.merge(ceusSolList.get(0), ceusSolList.get(1));
//		sol = MergedSolutionCreator.merge(wusSol, ceusSolList.get(0), ceusSolList.get(1));
//		sol = MergedSolutionCreator.merge(wusSol, ceusSolList.get(0), ceusSolList.get(1), solAK);
		

//		// This shows no duplicate parent IDs for Alaska
//		ArrayList<Integer> solParID_List = new ArrayList<Integer>();
//		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
//			if(!solParID_List.contains(sect.getParentSectionId())) {
//					solParID_List.add(sect.getParentSectionId());
//					System.out.println("Sol Parent: "+sect.getParentSectionId());
//			}
//		}
//		ArrayList<Integer> solAK_ParID_List = new ArrayList<Integer>();
//		for (FaultSection sect : solAK.getRupSet().getFaultSectionDataList()) {
//			if(!solAK_ParID_List.contains(sect.getParentSectionId())) {
//				solAK_ParID_List.add(sect.getParentSectionId());
//				System.out.println("AK Parent: "+sect.getParentSectionId());
//			}
//
//		}
//		for(Integer id:solParID_List) {
//			if(solAK_ParID_List.contains(id))
//				System.out.println("DUPLICATE: "+id+"\t");
//		}

				
				
//		System.out.println("DONE with sol");
//		System.exit(0);
		//

		// load DOLE data
//		DOLE_SubsectionMapper.mapDOLE(conusSol.getRupSet().getFaultSectionDataList(), DOLE_MappingAlgorithm.FULL_PARENT, true); // boolean is for verbose
//		DOLE_SubsectionMapper.mapDOLE(sol.getRupSet().getFaultSectionDataList(), PaleoMappingAlgorithm.CLOSEST_SECT, true); // boolean is for verbose


//		// ERF
//		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
//		
//		// prob model
//		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
//		erf.setParameter(HistoricOpenIntervalParam.NAME, 2024d-1875d); // or whatever
//		erf.getTimeSpan().setStartTime(2024);
//		erf.getTimeSpan().setDuration(50d);
//		
//		erf.updateForecast();
		
		
	}

}
