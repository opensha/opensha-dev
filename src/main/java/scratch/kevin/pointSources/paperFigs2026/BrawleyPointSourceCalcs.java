package scratch.kevin.pointSources.paperFigs2026;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

public class BrawleyPointSourceCalcs {

	public static void main(String[] args) {
		Location loc = new Location(32.986, -115.516, 15.4);
		double rake = 0d;
		PointSurface surf = new PointSurface(loc);
		surf.setDepth(loc.depth);
		surf.setAveDip(90);
		EqkRupture rup = new EqkRupture(4.3, rake, surf, loc);
		Site site = new Site(new Location(32.991515021954726, -115.51342661196108));  // brawley airport
		
		double sitePGA = 0.25;
		double logPGA = Math.log(sitePGA);
		
		ScalarIMR[] gmms = {
				AttenRelRef.ASK_2014.get(),
				AttenRelRef.BSSA_2014.get(),
				AttenRelRef.CB_2014.get(),
				AttenRelRef.CY_2014.get(),
				AttenRelRef.NGAWest_2014_AVG_NOIDRISS.get()
		};
		
		for (ScalarIMR gmm : gmms) {
			gmm.setParamDefaults();
			for (Parameter<?> param : gmm.getSiteParams())
				if (!site.containsParameter(param))
					site.addParameter(param);
			gmm.setIntensityMeasure(PGA_Param.NAME);
			gmm.setSite(site);
			gmm.setEqkRupture(rup);
			
			double logMean = gmm.getMean();
			double sd = gmm.getStdDev();
			double z = (logPGA - logMean)/sd;
			double mean = Math.exp(logMean);
			System.out.println(gmm.getShortName());
			System.out.println("\t"+gmm.getAllParamMetadata());
			System.out.println("\tMean:\tlog="+(float)logMean+", linear="+(float)mean);
			System.out.println("\tSigma:\t"+(float)sd);
			System.out.println("\tz:\t"+(float)z);
			System.out.println("\tEpsilon:\t"+(float)gmm.getEpsilon(logPGA));
		}
	}

}
