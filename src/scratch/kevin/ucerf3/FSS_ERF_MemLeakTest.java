package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.zip.ZipException;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;

public class FSS_ERF_MemLeakTest {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ZipException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws ZipException, IOException, InterruptedException {
		System.out.println("Start the profiler!");
		Thread.sleep(15000);
		
		File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/" +
				"InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		CompoundFaultSystemSolution fetch = CompoundFaultSystemSolution.fromZipFile(solFile);
		
		List<LogicTreeBranch> branches = Lists.newArrayList(fetch.getBranches());
		Collections.shuffle(branches);
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		GriddedRegion griddedRegion = new GriddedRegion(region, 1d, null);
		Random r = new Random();
		
		Site site = new Site(new Location(34, -119));
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		ScalarIMR imr = AttenRelRef.CB_2008.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(PGA_Param.NAME);
		
		site.addParameterList(imr.getSiteParams());
		
		int cnt = 0;
		for (LogicTreeBranch branch : branches) {
			FaultSystemSolution sol = fetch.getSolution(branch);
			System.gc();
			System.out.println("*** Loaded solution "+(cnt++));
			
//			UCERF3_FaultSysSol_ERF erf = new UCERF3_FaultSysSol_ERF(new InversionFaultSystemSolution(fetch.getSolution(branch)));
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			erf.updateForecast();
//			System.out.println("Calculating MFDs");
//			ERF_Calculator.getMagFreqDistInRegion(erf, region, 5d, 41, 0.1d, r.nextBoolean());
//			System.out.println("Calculating Participation Rates");
//			ERF_Calculator.getParticipationRatesInRegion(erf, griddedRegion, 6.7d, 10d);
			System.out.println("DONE");
			
			System.out.println("*** Calculating Hazard Curve ***");
			DiscretizedFunc curve = IMT_Info.getUSGS_PGA_Function();
			
			ArbitrarilyDiscretizedFunc logHazFunction = new ArbitrarilyDiscretizedFunc();
			for (int i = 0; i < curve.size(); ++i)
				logHazFunction.set(Math.log(curve.getX(i)), 1);
			calc.getHazardCurve(logHazFunction, site, imr, erf);
			System.out.println("*** DONE ***");
			
//			if (cnt == 3)
//				Thread.sleep(1000000000);
		}
	}

}
