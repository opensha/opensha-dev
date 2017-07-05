package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.TimeUnit;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrapper;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.base.Stopwatch;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class NGAW2_WrapperSpeedTest {

	public static void main(String[] args) throws IOException, DocumentException {
		
		int trials = 10;
		
		Stopwatch directWatch = Stopwatch.createUnstarted();
		Stopwatch paramWatch = Stopwatch.createUnstarted();
		
		ASK_2014 directGMPE = new ASK_2014();
		ASK_2014 paramGMPE = new ASK_2014();
		AttenuationRelationship directWrapper = new NGAW2_Wrapper("TestDirect", directGMPE);
		AttenuationRelationship paramWrapper = new NGAW2_WrapperFullParam("TestParam", paramGMPE, true);
		
		directWrapper.setParamDefaults();
		paramWrapper.setParamDefaults();
		
		directWrapper.setIntensityMeasure(PGA_Param.NAME);
		paramWrapper.setIntensityMeasure(PGA_Param.NAME);
		
		DiscretizedFunc xVals = new IMT_Info().getUSGS_PGA_Function();
		
		File fssFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_"
				+ "MEAN_BRANCH_AVG_SOL.zip");
		
		FaultSystemSolution directFSS = FaultSystemIO.loadSol(fssFile);
		FaultSystemSolution paramFSS = FaultSystemIO.loadSol(fssFile);
		
		FaultSystemSolutionERF directERF = new FaultSystemSolutionERF(directFSS);
		directERF.updateForecast();
		FaultSystemSolutionERF paramERF = new FaultSystemSolutionERF(paramFSS);
		paramERF.updateForecast();
		
		Site site = new Site();
		Vs30_Param vs30 = new Vs30_Param();
		vs30.setValue(760d);
		site.addParameter(vs30);
		Vs30_TypeParam vs30Type = new Vs30_TypeParam();
		vs30Type.setValue(Vs30_TypeParam.VS30_TYPE_INFERRED);
		site.addParameter(vs30Type);
		site.addParameter(new DepthTo1pt0kmPerSecParam());
		site.addParameter(new DepthTo2pt5kmPerSecParam());
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		for (int i=0; i<trials; i++) {
			Location loc = new Location(35d+Math.random(), -118+Math.random());
			site.setLocation(loc);
			
			DiscretizedFunc directCurve, paramCurve;
			if (Math.random() < 0.5) {
				// do direct first
				directCurve = calc(directWrapper, directERF, site, calc, xVals, directWatch);
				paramCurve = calc(paramWrapper, paramERF, site, calc, xVals, paramWatch);
			} else {
				// do param first
				paramCurve = calc(paramWrapper, paramERF, site, calc, xVals, paramWatch);
				directCurve = calc(directWrapper, directERF, site, calc, xVals, directWatch);
			}
			
			// now make sure they're equal
			double maxDiscrep = 0d;
			for (int j=0; j<xVals.size(); j++) {
				double discrep = Math.abs(paramCurve.getY(j) - directCurve.getY(j));
				maxDiscrep = Math.max(maxDiscrep, discrep);
			}
			System.out.println("Trial "+i+" max discrep: "+maxDiscrep);
			if (maxDiscrep > 0.05) {
//				System.out.println(directCurve);
//				System.out.println(paramCurve);
				
				// now check individual param vals
				System.out.println("Direct:\t"+directGMPE.toString());
				System.out.println("Param:\t"+paramGMPE.toString());
			}
		}
		
		long directSecs = directWatch.elapsed(TimeUnit.SECONDS);
		long paramSecs = paramWatch.elapsed(TimeUnit.SECONDS);
		
		System.out.println("Direct time: "+directSecs);
		System.out.println("Param time: "+paramSecs);
	}
	
	private static DiscretizedFunc calc(AttenuationRelationship gmpe, FaultSystemSolutionERF erf,
			Site site, HazardCurveCalculator calc, DiscretizedFunc xVals, Stopwatch watch) {
		DiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xVals);
		watch.start();
		calc.getHazardCurve(curve, site, gmpe, erf);
		watch.stop();
		return HazardCurveSetCalculator.unLogFunction(xVals, curve);
	}

}
