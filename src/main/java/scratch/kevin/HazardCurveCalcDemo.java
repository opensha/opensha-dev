package scratch.kevin;

import java.awt.geom.Point2D;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ERF_Ref;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

public class HazardCurveCalcDemo {

	public static void main(String[] args) {
		// choose an ERF, easy way is from the ERF_Ref enum (which contains all of them used in our apps):
		ERF erf = (ERF)ERF_Ref.MEAN_UCERF3.instance(); // casting here makes sure that it's a regular ERF and not an epistemic list erf
		
		System.out.println("ERF: "+erf.getName());
		// let's view the adjustable parameters
		System.out.println("ERF Params");
		for (Parameter<?> param : erf.getAdjustableParameterList())
			System.out.println("\t"+param.getName()+": "+param.getValue());
		// same for the time span
		TimeSpan timeSpan = erf.getTimeSpan();
		System.out.println("Duration: "+timeSpan.getDuration());
		
		// if you want to change any parameters or the time span, do them now
		
		// finally, call the updateForecast() method which builds the forecast for the current set of parameters
		System.out.println("Updating forecast...");
		erf.updateForecast();
		System.out.println("DONE updating forecast");
		
		// now choose a GMM, easiest way is from the AttenRelRef enum
		ScalarIMR gmm = AttenRelRef.ASK_2014.instance(null);
		
		System.out.println("GMM: "+gmm.getName());
		
		// set paramter defaults
		gmm.setParamDefaults();
		
		// set the intensity measure type
		gmm.setIntensityMeasure(SA_Param.NAME);
		
		// if SA, set the period
		SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), 1d);
		
		// build a site
		Site site = new Site(new Location(34, -118));
		System.out.println("Site location: "+site.getLocation());
		
		// need to add site parameters to the site that will be used by the GMM
		System.out.println("Site parameters:");
		for (Parameter<?> param : gmm.getSiteParams()) {
			// if you're be doing calculation with multiple sites, best to clone it so that you can set a parameter
			// for one site without changing it in another
			param = (Parameter<?>)param.clone();
			// could set it's value here
			System.out.println("\t"+param.getName()+": "+param.getValue());
			// add that parameter to the site
			site.addParameter(param);
		}
		
		// create a hazard curve calculator
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		// get default x-values for hazard calculation
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmm.getIntensityMeasure());
		
		// create the same but in ln spacing
		DiscretizedFunc lnXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			lnXVals.set(Math.log(pt.getX()), 1d);
		
		// calculate curve, will fill in the y-values from the above
		calc.getHazardCurve(lnXVals, site, gmm, erf);
		
		// now combine those y values with the original linear x values
		DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			curve.set(xVals.getX(i), lnXVals.getY(i));
		
		System.out.println("Hazard curve:\n"+curve);
	}

}
