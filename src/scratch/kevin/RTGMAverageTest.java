package scratch.kevin;

import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.calc.mcer.AbstractMCErProbabilisticCalc;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.collect.Lists;

public class RTGMAverageTest {
	
	public static void main(String[] args) {
		Location loc = new Location(35, -118);
		Site site = new Site(loc);
		
		Vs30_Param vs30 = new Vs30_Param();
		vs30.setValue(760);
		site.addParameter(vs30);
		Vs30_TypeParam vs30TypeParam = new Vs30_TypeParam();
		vs30TypeParam.setValue(Vs30_TypeParam.VS30_TYPE_INFERRED);
		site.addParameter(vs30TypeParam);
		DepthTo1pt0kmPerSecParam z10param = new DepthTo1pt0kmPerSecParam();
		z10param.setValue(null);
		site.addParameter(z10param);
		DepthTo2pt5kmPerSecParam z25param = new DepthTo2pt5kmPerSecParam(null, 0d, 100000d, true);
		z25param.setValue(null);
		site.addParameter(z25param);
		
		DiscretizedFunc xVals = IMT_Info.getUSGS_SA_Function();
		
		ERF erf = new MeanUCERF2();
		erf.updateForecast();
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		List<ScalarIMR> gmpes = Lists.newArrayList();
		gmpes.add(AttenRelRef.ASK_2014.instance(null));
		gmpes.add(AttenRelRef.BSSA_2014.instance(null));
		gmpes.add(AttenRelRef.CB_2014.instance(null));
		gmpes.add(AttenRelRef.CY_2014.instance(null));
		
		for (int j=0; j<2; j++) {
			double weightEach = 1d/gmpes.size();
			DiscretizedFunc meanCurve = null;
			double[] rtgmVals = new double[gmpes.size()];

			for (int n=0; n<gmpes.size(); n++) {
				ScalarIMR gmpe = gmpes.get(n);
				gmpe.setParamDefaults();
				gmpe.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 3d);

				DiscretizedFunc hazFunction = HazardCurveSetCalculator.getLogFunction(xVals);
				calc.getHazardCurve(hazFunction, site, gmpe, erf);
				hazFunction = HazardCurveSetCalculator.unLogFunction(xVals, hazFunction);

				if (meanCurve == null) {
					meanCurve = new ArbitrarilyDiscretizedFunc();
					for (int i=0; i<xVals.size(); i++)
						meanCurve.set(xVals.getX(i), 0d);
				}

				for (int i=0; i<meanCurve.size(); i++)
					meanCurve.set(i, meanCurve.getY(i) + weightEach*hazFunction.getY(i));

//				System.out.println(hazFunction);
				rtgmVals[n] = AbstractMCErProbabilisticCalc.calcRTGM(hazFunction);
			}

			double rtgmFromMean = AbstractMCErProbabilisticCalc.calcRTGM(meanCurve);
			double meanRTGM = StatUtils.mean(rtgmVals);

			System.out.println("RTGM From Mean Curve: "+rtgmFromMean);
			System.out.println("Mean RTGM: "+meanRTGM);
		}
	}

}
