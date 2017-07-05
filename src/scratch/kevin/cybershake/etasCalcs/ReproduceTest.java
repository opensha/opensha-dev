package scratch.kevin.cybershake.etasCalcs;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator;
import scratch.UCERF3.erf.ETAS.FaultSystemSolutionERF_ETAS;
import scratch.UCERF3.erf.ETAS.ETAS_Params.ETAS_ParameterList;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;
import scratch.kevin.ucerf3.etas.MPJ_ETAS_Simulator;

public class ReproduceTest {

	public static void main(String[] args) throws IOException, DocumentException {
		File csDir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/");
		FaultSystemSolution fss = FaultSystemIO.loadSol(new File(csDir, "ucerf2_mapped_sol.zip"));
		File resultsDir = new File("/tmp");
		FaultSystemSolutionERF_ETAS erf = MPJ_ETAS_Simulator.buildERF(fss, false, 1d);
		GriddedRegion reg = new CaliforniaRegions.RELM_GRIDDED();
		
		long randSeed = 1408453138855l;
		
		boolean includeEqkRates = true;
		boolean includeIndirectTriggering = true;
		boolean includeSpontEvents = false;
		
		double gridSeisDiscr = 0.1;
		
		ObsEqkRupList histQkList = new ObsEqkRupList();
		
		List<float[]> fractionSrcAtPointList = MatrixIO.floatArraysListFromFile(
				new File(csDir, "fractionSrcAtPointList.bin"));
		List<int[]> srcAtPointList = MatrixIO.intArraysListFromFile(new File(csDir, "srcAtPointList.bin"));
		int[] isCubeInsideFaultPolygon = MatrixIO.intArrayFromFile(new File(csDir, "isCubeInsideFaultPolygon.bin"));
		
		ETAS_Simulator.D = false;
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55;
		
		long ot = Math.round((2014.0-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR); // occurs at 2014
		ETAS_EqkRupture mainshockRup = new ETAS_EqkRupture();
		mainshockRup.setOriginTime(ot);
		
		// Mojave M 7.05 rupture
		int fssScenarioRupID = 30473;
		mainshockRup.setAveRake(fss.getRupSet().getAveRakeForRup(fssScenarioRupID));
		mainshockRup.setMag(fss.getRupSet().getMagForRup(fssScenarioRupID));
		mainshockRup.setRuptureSurface(fss.getRupSet().getSurfaceForRupupture(fssScenarioRupID, 1d, false));
		mainshockRup.setID(0);
		erf.setFltSystemSourceOccurranceTimeForFSSIndex(fssScenarioRupID, ot);
		
		erf.updateForecast();
		
		ETAS_Simulator.testETAS_Simulation(resultsDir, erf, reg, mainshockRup, histQkList, includeSpontEvents,
				includeIndirectTriggering, gridSeisDiscr, null, randSeed,
				fractionSrcAtPointList, srcAtPointList, isCubeInsideFaultPolygon, new ETAS_ParameterList());
	}

}
