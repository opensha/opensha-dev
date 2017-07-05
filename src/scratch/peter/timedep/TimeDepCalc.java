package scratch.peter.timedep;

import java.util.concurrent.Callable;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.imr.NSHMP08_WUS;
import org.opensha.nshmp2.util.Period;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import scratch.peter.nshmp.NSHMP_UCERF2_ERF;

/**
 * Standalone calculator class for NSHMP_ListERFs. Assumes Poissonian.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class TimeDepCalc implements Callable<HazardResult> {

	private ERF erf;
	private SourceIMR imrRef;
	private Site site;
	private Period period;
	private boolean epiUncert;
	private HazardCurveCalculator calc;
	
	private TimeDepCalc() {}
	
	public static TimeDepCalc create(ERF erf, SourceIMR imrRef, Site site,
			Period period) {
		TimeDepCalc hc = new TimeDepCalc();
		hc.erf = erf;
		hc.site = site;
		hc.period = period;
		hc.epiUncert = false;
		hc.imrRef = imrRef;
		return hc;
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public HazardResult call() {
		initSite(site); // ensure required site parameters are set
		ScalarIMR imr = imrRef.instance(period);
		imr.getParameter(NSHMP08_WUS.IMR_UNCERT_PARAM_NAME).setValue(epiUncert);
		imr.setSite(site);
		
		calc = new HazardCurveCalculator(); // init calculator
		calc.setPtSrcDistCorrType(PtSrcDistCorr.Type.NSHMP08);
		calc.setMaxSourceDistance(erf instanceof NSHMP_UCERF2_ERF ? 200.0 : 300.0);
		DiscretizedFunc curve = period.getFunction(); // utility function
		curve = calc.getHazardCurve(curve, site, imr, erf);

		return new HazardResult(period, site.getLocation(), curve, null);
	}
	
	
	
	public static void initSite(Site s) {
		
		// CY AS
		DepthTo1pt0kmPerSecParam d10p = new DepthTo1pt0kmPerSecParam(null,
			0, 1000, true);
		d10p.setValueAsDefault();
		s.addParameter(d10p);
		// CB
		DepthTo2pt5kmPerSecParam d25p = new DepthTo2pt5kmPerSecParam(null,
			0, 1000, true);
		d25p.setValueAsDefault();
		s.addParameter(d25p);
		// all
		Vs30_Param vs30p = new Vs30_Param(760);
		vs30p.setValueAsDefault();
		s.addParameter(vs30p);
		// AS CY
		Vs30_TypeParam vs30tp = new Vs30_TypeParam();
		vs30tp.setValueAsDefault();
		s.addParameter(vs30tp);
		
	}

}
