package scratch.kevin.simulators.ruptures;

import java.io.IOException;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import scratch.kevin.bbp.BBP_Site;

public class RSQSimBBP_HazardCurveCalc {
	
	private BBP_CatalogSimZipLoader zip;
	private double durationYears;
	private DiscretizedFunc xVals;
	
	private static DiscretizedFunc getDefaultHazardCurve(int xValMult) {
		ArbitrarilyDiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
		if (xValMult > 0) {
			ArbitrarilyDiscretizedFunc newXValues = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<xValues.size()-1; i++) {
				double x0 = Math.log(xValues.getX(i));
				double x1 = Math.log(xValues.getX(i+1));
				double dx = (x1 - x0)/xValMult;
				for (int j=0; j<xValMult; j++)
					newXValues.set(Math.exp(x0 + j*dx), 1d);
			}
			newXValues.set(xValues.getMaxX(), 1d);
			xValues = newXValues;
		}
		return xValues;
	}
	
	public RSQSimBBP_HazardCurveCalc(BBP_CatalogSimZipLoader zip, double durationYears) {
		this(zip, durationYears, getDefaultHazardCurve(4));
	}
	
	public RSQSimBBP_HazardCurveCalc(BBP_CatalogSimZipLoader zip, double durationYears, DiscretizedFunc xVals) {
		this.zip = zip;
		this.durationYears = durationYears;
		this.xVals = xVals;
	}
	
	public DiscretizedFunc calc(BBP_Site site, double period, double curveDuration) throws IOException {
		double rateEach = 1d/durationYears;
		
		// annual rate curve
		DiscretizedFunc curve = xVals.deepClone();
		for (int i=0; i<curve.size(); i++)
			curve.set(i, 0d);
		for (int eventID : zip.getEventIDs(site)) {
			DiscretizedFunc spectra = zip.readRotD50(site, eventID);
			double rd50 = spectra.getInterpolatedY(period);
			for (int i=0; i<curve.size(); i++)
				if (curve.getX(i) <= rd50)
					curve.set(i, curve.getY(i)+rateEach);
		}
		
		// now probabilities
		for (int i=0; i<curve.size(); i++) {
			double rate = curve.getY(i);
			double prob = 1d - Math.exp(-rate*curveDuration);
			curve.set(i, prob);
		}
		
		return curve;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
