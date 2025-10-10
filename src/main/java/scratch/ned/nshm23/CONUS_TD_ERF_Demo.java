package scratch.ned.nshm23;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.util.MergedSolutionCreator;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.HistoricalRupture;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoDOLE_Data;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoMappingAlgorithm;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen.DataToInclude;
import org.opensha.sha.faultSurface.FaultSection;


import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.ned.nshm23.AK_FSS_creator.DeformationModelEnum;
import scratch.ned.nshm23.CEUS_FSS_creator.FaultModelEnum;

/**
 * NOTES:
 * 
 * Parent section names are null for AK fault system solution faults (info not readily available, but IDs are good)
 */
public class CONUS_TD_ERF_Demo {
	
	
	
	/**
	 * This confirms that the implied sum of rates for the ERF rups for each FSS source 
	 * are consistent with the FSS rates.
	 */
	private static void testPoisERF_ProbsVsFSS_rates() {
		
		// test threshold of numerical accuracy for prob --> rate conversions
		
		double t=1; // duration
		System.out.println("Numerical rate-->prob-->rate test");

		for(int i=1; i<20; i++) {
			double rate = Math.pow(10,-i);
			double prob = 1-Math.exp(-rate*t); // T=1
//			if(rate <= 1e-8)
//				prob=rate;
			double rate2 = -Math.log(1 - prob)/t;
			System.out.println((float)rate+"\t"+(float)(rate2/rate)+"\t"+(float)prob);
		}
		
		String fssFileName = "results_WUS_FM_v3_branch_averaged_gridded_simplified.zip";
		double forecastDurationYears = 50;
		
		// get solution
		FaultSystemSolution sol = null;
		try {
			sol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/"+fssFileName));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(forecastDurationYears);
		erf.updateForecast();
		
		// this shows that the mfd for each source (more than one rupture) exactly equals the FSS rates
		// I had to temporarily make getFaultSysRupMFD() public to get this to work
		double minRatio=Double.MAX_VALUE;
		double maxRatio=Double.MIN_VALUE;
		for(int s=0;s<erf.getNumFaultSystemSources();s++) {
			int fsrIndex = erf.getFltSysRupIndexForSource(s);
			RupMFDsModule rupMFDs = ((BaseFaultSystemSolutionERF)erf).getSolution().getModule(RupMFDsModule.class);
			DiscretizedFunc mfd = rupMFDs.getRuptureMFD(fsrIndex);
			if(mfd == null)
				continue;
			double rate=mfd.calcSumOfY_Vals();
			double rate2 = sol.getRateForRup(fsrIndex);
			double ratio = rate/rate2;
			if(maxRatio<ratio) maxRatio=ratio;
			if(minRatio>ratio) minRatio=ratio;
		}
		System.out.println("DONE verifying rup mfd vs FSS rate; minRatio="+(float)minRatio+", maxRatio="+(float)maxRatio);

		double numLowProbRups=0;
		double numRups=0;
		for(int s=0;s<erf.getNumFaultSystemSources();s++) {
			int fsrIndex = erf.getFltSysRupIndexForSource(s);
			double rate=0;
			for(int r=0;r<erf.getSource(s).getNumRuptures();r++) {
				rate+=erf.getSource(s).getRupture(r).getMeanAnnualRate(forecastDurationYears);
				numRups+=1;
				if(erf.getSource(s).getRupture(r).getProbability()<1e-15)
					numLowProbRups+=1;
			}
			double rate2 = sol.getRateForRup(fsrIndex);
			double ratio = rate/rate2;
			if(ratio > 1.0001 || ratio < 0.9999) {
				if(rate2>1e-13)
					System.out.println("Diff\t"+s+"\t"+rate+"\t"+rate2+"\t"+ratio);
			}
		}
		double perc = 100d*(double)numLowProbRups/(double)numRups;
		System.out.println("DONE verifying ERF vs FSS rate; numerical problems for FSS rates less the 1e-13 for T=50");
		System.out.println(numLowProbRups+" Ruptures have prob<1e-15 (out of "+numRups+ ", which is "+(float)perc+"%)");
	}

	public static FaultSystemSolution getCEUS_FSS(String ceus_FSS_fileName) {

		// try reading from file
		if(ceus_FSS_fileName != null) {
			File file = new File(ceus_FSS_fileName);
			if(file.exists()) {
				try {
					return FaultSystemSolution.load(file);
				} catch (IOException e) {
					e.printStackTrace();
				};
			}
		}	

		// create FSS
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.1.2/";
		ArrayList<FaultSystemSolution>  ceusSolList = CEUS_FSS_creator.getFaultSystemSolutionList(nshmModelDirPath,CEUS_FSS_creator.FaultModelEnum.PREFERRED);
		FaultSystemSolution sol = MergedSolutionCreator.merge(ceusSolList.get(0), ceusSolList.get(1));
		if(ceus_FSS_fileName != null) {
			try {
				sol.write(new File(ceus_FSS_fileName));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return sol;
	}
	
	
	public static FaultSystemSolution getAK_FSS(String alaska_FSS_fileName) {

		// try reading from file
		if(alaska_FSS_fileName != null) {
			File file = new File(alaska_FSS_fileName);
			if(file.exists()) {
				try {
					return FaultSystemSolution.load(file);
				} catch (IOException e) {
					e.printStackTrace();
				};
			}
		}	

		// create FSS
		String akModelDirPath = "/Users/field/nshm-haz_data/nshm-alaska-3.0.1/";
		ArrayList<FaultSystemSolution> ak_fssList = AK_FSS_creator.getFaultSystemSolutionList(akModelDirPath, DeformationModelEnum.GEO);
		FaultSystemSolution sol = MergedSolutionCreator.merge(ak_fssList);
		if(alaska_FSS_fileName != null) {
			try {
				sol.write(new File(alaska_FSS_fileName));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return sol;
	}


	
	public static FaultSystemSolution getFull_FSS(String full_FSS_fileName) {
		
		FaultSystemSolution sol=null;
		
		if(full_FSS_fileName != null) {
			File file = new File(full_FSS_fileName);
			if(file.exists()) {
				try {
					return FaultSystemSolution.load(file);
				} catch (IOException e) {
					e.printStackTrace();
				};
			}
		}	
		
		// Make from scratch
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

//			System.out.println("Testing rup areas for each region:");
//			testRupAreas(wusSol);
//			for(FaultSystemSolution fss:ceusSolList) testRupAreas(fss);
//			for(FaultSystemSolution fss:ak_fssList) testRupAreas(fss);
//
//			System.out.println("Testing rup areas for merged sol:");
//			testRupAreas(sol);
			
			System.out.println("wusSol\t"+wusSol.getRupSet().getMagForRup(0)+"\t"+wusSol.getRupSet().getAreaForRup(0));
			System.out.println("ceusSolList\t"+ceusSolList.get(0).getRupSet().getMagForRup(0)+"\t"+ceusSolList.get(0).getRupSet().getAreaForRup(0));
			System.out.println("ak_fssList\t"+ak_fssList.get(0).getRupSet().getMagForRup(0)+"\t"+ak_fssList.get(0).getRupSet().getAreaForRup(0));

		} catch (IOException e) {
			e.printStackTrace();
		}		
		
		if(full_FSS_fileName != null) {
			try {
				sol.write(new File(full_FSS_fileName));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return sol;
	}
	
	

	public static void testRupAreas(FaultSystemSolution sol) {
		int numZeroAreaRups=0;
		int numZeroLengthRups=0;
		for(int r=0;r<sol.getRupSet().getNumRuptures();r++) {
			if(sol.getRupSet().getAreaForRup(r) == 0d)
				numZeroAreaRups+=1;
			if(sol.getRupSet().getLengthForRup(r) == 0d)
				numZeroLengthRups+=1;
		}
		System.out.println("numZeroAreaRups="+numZeroAreaRups+" for "+sol.getName()+"\t"+
				sol.getInfoString()+"\ttotNumRups="+sol.getRupSet().getNumRuptures());
		System.out.println("numZeroLengthRups="+numZeroLengthRups+" for "+sol.getName()+"\t"+
				sol.getInfoString()+"\ttotNumRups="+sol.getRupSet().getNumRuptures());
	}
	

	public static void main(String[] args) throws IOException {
	
		getFull_FSS(null);
		
//		testRupAreas(getFull_FSS("/Users/field/nshm-haz_data/full_FSS_test.zip"));		

		
//		testPoisERF_ProbsVsFSS_rates();
//		System.exit(0);
		

		
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
