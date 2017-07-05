package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator;
import scratch.UCERF3.erf.ETAS.FaultSystemSolutionERF_ETAS;
import scratch.UCERF3.erf.ETAS.ETAS_Params.ETAS_ParameterList;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.RELM_RegionUtils;

public class CacheSpeedTester {

	public static void main(String[] args) throws IOException, DocumentException {
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55;
		double size = 0d;
		File mainResultsDir = null;
		int forceRunNum = -1;
		if (args.length == 0) {
			size = 2d;
			mainResultsDir = new File("/home/kevin/OpenSHA/UCERF3/etas/cache_tests");
		} else if (args.length == 2 || args.length == 3) {
			mainResultsDir = new File(args[0]);
			size = Double.parseDouble(args[1]);
			if (args.length == 3)
				forceRunNum = Integer.parseInt(args[2]);
		} else {
			System.err.println("USAGE: <mian-dir> <size gb> <run #>");
			System.exit(2);
		}
		File inputDir = mainResultsDir.getParentFile();
		System.setProperty("etas.cache.size.gb", size+"");
//		File solFile = new File(args[0]);
//		FaultSystemSolution fss = FaultSystemIO.loadSol(solFile);
//		FaultSystemSolutionERF_ETAS erf = MPJ_ETAS_Simulator.buildERF(fss, false, 1d);
//		FaultSystemSolutionERF_ETAS erf = ETAS_Simulator.getU3_ETAS_ERF(fss);
		FaultSystemSolutionERF_ETAS erf = ETAS_Simulator.getU3_ETAS_ERF(2014d,1d);
		erf.updateForecast();
		
		int numRuns = 5;
		if (numRuns <= forceRunNum)
			numRuns = forceRunNum+1;
		
		File resultsDir = new File(mainResultsDir, "cache_"+(float)size+"gb");
		if (!resultsDir.exists())
			resultsDir.mkdir();
		GriddedRegion reg = RELM_RegionUtils.getGriddedRegionInstance();
		
		long randSeed = 1408453138855l;
		
		boolean includeIndirectTriggering = true;
		boolean includeSpontEvents = true;
		
		double gridSeisDiscr = 0.1;
		
		ObsEqkRupList histQkList = new ObsEqkRupList();
		
		ETAS_Simulator.D = false;
		
//		ETAS_EqkRupture mainshockRup = null;
		ETAS_EqkRupture mainshockRup = ETAS_Simulator.buildScenarioRup(TestScenario.MOJAVE_M7, erf);
//		ETAS_EqkRupture mainshockRup = new ETAS_EqkRupture();
//		long ot = Math.round((2014.0-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR); // occurs at 2014
//		mainshockRup.setOriginTime(ot);
//		
//		// Mojave M 7.05 rupture
//		int fssScenarioRupID = 30473;
//		mainshockRup.setAveRake(fss.getRupSet().getAveRakeForRup(fssScenarioRupID));
//		mainshockRup.setMag(fss.getRupSet().getMagForRup(fssScenarioRupID));
//		mainshockRup.setRuptureSurface(fss.getRupSet().getSurfaceForRupupture(fssScenarioRupID, 1d, false));
//		mainshockRup.setID(0);
//		erf.setFltSystemSourceOccurranceTimeForFSSIndex(fssScenarioRupID, ot);
//		
//		erf.updateForecast();
		
		File fractionSrcAtPointListFile = new File(inputDir, "fractionSrcAtPointList.bin");
		File srcAtPointListFile = new File(inputDir, "srcAtPointList.bin");
		File isCubeInsideFaultPolygonFile = new File(inputDir, "isCubeInsideFaultPolygon.bin");
		Preconditions.checkState(fractionSrcAtPointListFile.exists(),
				"cache file not found: "+fractionSrcAtPointListFile.getAbsolutePath());
		Preconditions.checkState(srcAtPointListFile.exists(),
				"cache file not found: "+srcAtPointListFile.getAbsolutePath());
		Preconditions.checkState(isCubeInsideFaultPolygonFile.exists(),
				"cache file not found: "+isCubeInsideFaultPolygonFile.getAbsolutePath());
		List<float[]> fractionSrcAtPointList = MatrixIO.floatArraysListFromFile(fractionSrcAtPointListFile);
		List<int[]> srcAtPointList = MatrixIO.intArraysListFromFile(srcAtPointListFile);
		int[] isCubeInsideFaultPolygon = MatrixIO.intArrayFromFile(isCubeInsideFaultPolygonFile);
		
		for (int run=0; run<numRuns; run++) {
			if (forceRunNum >= 0 && run != forceRunNum)
				continue;
			File subdir = new File(resultsDir, "run_"+run);
			if (!subdir.exists())
				subdir.mkdir();
			System.gc();
			ETAS_Simulator.testETAS_Simulation(subdir, erf, reg, mainshockRup, histQkList, includeSpontEvents,
					includeIndirectTriggering, gridSeisDiscr, null, randSeed, fractionSrcAtPointList, srcAtPointList,
					isCubeInsideFaultPolygon, new ETAS_ParameterList());
		}
		
//		ETAS_Simulator.testETAS_Simulation(resultsDir, erf, reg, mainshockRup, histQkList, includeSpontEvents,
//				includeIndirectTriggering, includeEqkRates, gridSeisDiscr, null, randSeed, new ETAS_ParameterList());
	}

}
