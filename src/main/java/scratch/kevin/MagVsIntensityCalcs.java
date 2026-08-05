package scratch.kevin;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;
import org.opensha.sha.imr.param.IntensityMeasureParams.MMI_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;

public class MagVsIntensityCalcs {

	public static void main(String[] args) {
		
		double mileToKm = 1.60934; // 10 miles
		
		double[] distsMiles = { 10d, 20d, 50d, 100d };
//		ScalarIMR gmm = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.get();
		ScalarIMR gmm = AttenRelRef.ASK_2014.get();
		
		String imt = PGA_Param.NAME;
//		String imt = PGV_Param.NAME;
//		String imt = MMI_Param.NAME;
		
		Location siteLoc = new Location(0d, 0d);
		
		PointSurface m5surf = new PointSurface(siteLoc);
		m5surf.setDepth(5d);
		
		Site site = new Site(siteLoc);
		site.addParameterList(gmm.getSiteParams());
		
		EqkRupture m5Rup = new EqkRupture(5d, 0d, m5surf, null);
		
		double mmi5 = calcGM(gmm, site, m5Rup, imt);
		System.out.println("M5 underneath "+imt+": "+(float)mmi5);
		
		for (double distMiles : distsMiles) {
			double distKM = distMiles * mileToKm;
			Location rupLoc = LocationUtils.location(siteLoc, 0d, distKM);
			PointSurface surf = new PointSurface(rupLoc);
			surf.setDepth(1d);
			
			double closestMag = Double.NaN;
			double closestMMI = Double.NaN;
			double closestMMIdiff = Double.POSITIVE_INFINITY;
			EvenlyDiscretizedFunc mags = new EvenlyDiscretizedFunc(5d, 51, 0.1);
			for (int m=0; m<mags.size(); m++) {
				double mag = mags.getX(m);
				EqkRupture rup = new EqkRupture(mag, 0d, surf, null);
				double mmi = calcGM(gmm, site, rup, imt);
				double diff = Math.abs(mmi - mmi5);
				if (diff < closestMMIdiff) {
					closestMMIdiff = diff;
					closestMMI = mmi;
					closestMag = mag;
				}
			}
			
			System.out.println("Mag for matching "+imt+" at "+(int)distMiles+" miles: M"+(float)closestMag+": "+(float)closestMMI);
		}
	}
	
	private static double calcGM(ScalarIMR gmm, Site site, EqkRupture rup, String imt) {
		gmm.setSite(site);
		gmm.setEqkRupture(rup);
		if (imt.equals(MMI_Param.NAME)) {
			gmm.setIntensityMeasure(PGA_Param.NAME);
			double pga = gmm.getMean();
			gmm.setIntensityMeasure(PGV_Param.NAME);
			double pgv = gmm.getMean();
			return Wald_MMI_Calc.getMMI(Math.exp(pga), Math.exp(pgv));
		} else {
			gmm.setIntensityMeasure(imt);
			return Math.exp(gmm.getMean());
		}
	}

}
