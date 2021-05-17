package scratch.kevin;

import static org.opensha.commons.geo.GeoTools.TO_RAD;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.TimeUnit;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CB_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.imr.attenRelImpl.ngaw2.IMT;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.base.Stopwatch;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class IMRParamSpeedTest {

	public static void main(String[] args) throws IOException, DocumentException {
		File baSolFile = new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(FaultSystemIO.loadSol(baSolFile));
		erf.updateForecast();
		
		boolean direct = true;
		boolean cb08 = true;
		
		ScalarIMR imr;
		if (cb08) {
			if (direct)
				imr = new TestCB_2008_UsesPrimitives_AttenRel(null);
			else
				imr = AttenRelRef.CB_2008.instance(null);
		} else {
			imr = AttenRelRef.CB_2014.instance(null);
		}
		CB_2014 cb = new CB_2014();
		imr.setIntensityMeasure(PGA_Param.NAME);
		imr.setParamDefaults();
		
		Site site = new Site(new Location(34, -120));
		site.addParameter(imr.getParameter(Vs30_Param.NAME));
		if (!cb08)
			site.addParameter(imr.getParameter(Vs30_TypeParam.NAME));
		if (!cb08)
			site.addParameter(imr.getParameter(DepthTo1pt0kmPerSecParam.NAME));
		site.addParameter(imr.getParameter(DepthTo2pt5kmPerSecParam.NAME));
		
		System.out.println("Starting test: direct="+direct+" cb08="+cb08);
		Stopwatch watch = Stopwatch.createStarted();
		
		double totMean = 0;
		double totStdDev = 0;
		
		for (ProbEqkSource src : erf) {
			for (ProbEqkRupture rup : src) {
				double mean, stdDev;
				if (!direct || cb08) {
					imr.setAll(rup, site, imr.getIntensityMeasure());
					mean = imr.getMean();
					stdDev = imr.getStdDev();
				} else {
					RuptureSurface surf = rup.getRuptureSurface();
					Location siteLoc = site.getLocation();
					double rJB = surf.getDistanceJB(siteLoc);
					double rRup = surf.getDistanceRup(siteLoc);
					double rX = surf.getDistanceX(siteLoc);
					double dip = surf.getAveDip();
					double width = surf.getAveWidth();
					double zTop = surf.getAveRupTopDepth();
					double zHyp = surf.getAveRupTopDepth() +
							Math.sin(surf.getAveDip() * TO_RAD) * surf.getAveWidth() /
							2.0;
					double vs30 = (Double)imr.getParameter(Vs30_Param.NAME).getValue();
					Double z2p5 = (Double)imr.getParameter(DepthTo2pt5kmPerSecParam.NAME).getValue();
					if (z2p5 == null)
						z2p5 = Double.NaN;
					double rake = rup.getAveRake();
					// assert in range [-180 180]
					FaultUtils.assertValidRake(rake);
					FaultStyle style = getFaultStyle(rake);
					ScalarGroundMotion gm = cb.calc(IMT.PGA, rup.getMag(), rJB, rRup, rX, dip, width, zTop, zHyp, vs30, z2p5, style);
					mean = gm.mean();
					stdDev = gm.stdDev();
				}
				totMean += mean;
				totStdDev += stdDev;
			}
		}
		watch.stop();
		System.out.println(watch.elapsed(TimeUnit.SECONDS)+"s");
		System.out.println("Tot mean: "+totMean);
		System.out.println("Tot stdDev: "+totStdDev);
	}
	
	static FaultStyle getFaultStyle(double rake) {
		if (rake >= 135 || rake <= -135)
			// right lateral
			return FaultStyle.STRIKE_SLIP;
		else if (rake >= -45 && rake <= 45)
			// left lateral
			return FaultStyle.STRIKE_SLIP;
		else if (rake >= 45 && rake <= 135)
			return FaultStyle.REVERSE;
		else
			return FaultStyle.NORMAL;
	}

}
