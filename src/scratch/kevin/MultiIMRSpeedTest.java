package scratch.kevin;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.attenRelImpl.AS_2008_AttenRel;
import org.opensha.sha.imr.attenRelImpl.BA_2008_AttenRel;
import org.opensha.sha.imr.attenRelImpl.CB_2008_AttenRel;
import org.opensha.sha.imr.attenRelImpl.CY_2008_AttenRel;
import org.opensha.sha.imr.attenRelImpl.NGA_2008_Averaged_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

public class MultiIMRSpeedTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		NGA_2008_Averaged_AttenRel multi = new NGA_2008_Averaged_AttenRel(null);
		multi.setParamDefaults();
		
		AS_2008_AttenRel as = new AS_2008_AttenRel(null);
		as.setParamDefaults();
		CB_2008_AttenRel cb = new CB_2008_AttenRel(null);
		cb.setParamDefaults();
		BA_2008_AttenRel ba = new BA_2008_AttenRel(null);
		ba.setParamDefaults();
		CY_2008_AttenRel cy = new CY_2008_AttenRel(null);
		cy.setParamDefaults();
		
		HazardCurveCalculator hc = new HazardCurveCalculator();
		
		MeanUCERF2 ucerf = new MeanUCERF2();
		ucerf.updateForecast();
		
		String imt = PGA_Param.NAME;
		multi.setIntensityMeasure(imt);
		as.setIntensityMeasure(imt);
		ba.setIntensityMeasure(imt);
		cb.setIntensityMeasure(imt);
		cy.setIntensityMeasure(imt);
		
		IMT_Info imtInfo = new IMT_Info();
		
		Site site = new Site(new Location(34, -118));
		site.addParameter(cb.getParameter(Vs30_Param.NAME));
		site.addParameter(cb.getParameter(DepthTo2pt5kmPerSecParam.NAME));
		site.addParameter(as.getParameter(DepthTo1pt0kmPerSecParam.NAME));
		site.addParameter(as.getParameter(Vs30_TypeParam.NAME));
		
		long startMulti = System.currentTimeMillis();
		hc.getHazardCurve(imtInfo.getDefaultHazardCurve(imt), site, multi, ucerf);
		long multiMillis = System.currentTimeMillis() - startMulti;
		double multiSecs = (double)multiMillis / 1000d;
		System.out.println("Multi: " + multiSecs);
		
		long startIndv = System.currentTimeMillis();
		hc.getHazardCurve(imtInfo.getDefaultHazardCurve(imt), site, as, ucerf);
		hc.getHazardCurve(imtInfo.getDefaultHazardCurve(imt), site, ba, ucerf);
		hc.getHazardCurve(imtInfo.getDefaultHazardCurve(imt), site, cb, ucerf);
		hc.getHazardCurve(imtInfo.getDefaultHazardCurve(imt), site, cy, ucerf);
		long indvMillis = System.currentTimeMillis() - startIndv;
		
		double indvSecs = (double)indvMillis / 1000d;
		
		System.out.println("Indv: " + indvSecs);
	}

}
