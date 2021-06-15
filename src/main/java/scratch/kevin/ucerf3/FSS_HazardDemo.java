package scratch.kevin.ucerf3;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.MaxDistanceParam;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class FSS_HazardDemo {

	public static void main(String[] args) throws IOException, DocumentException {
		File fssFile = new File("/home/kevin/workspace/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		System.out.println("Loading solution");
		FaultSystemSolution fss = FaultSystemIO.loadSol(fssFile);
		
		System.out.println("Building ERF");
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fss);
		erf.getTimeSpan().setDuration(50d); // 50 years
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.updateForecast();
		System.out.println("ERF has "+erf.getNumSources()+" sources");
		
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		// for PGA (units: g)
//		gmpe.setIntensityMeasure(PGA_Param.NAME);
		// for 1s SA (units: g)
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 1d); 
		
		Site site = new Site(new Location(34, -118));
		// simplest: use GMPE defaults
//		site.addParameterList(gmpe.getSiteParams());
		// custom Vs30
		System.out.println("Site parameters:");
		for (Parameter<?> siteParam : gmpe.getSiteParams()) {
			// need to clone it if you're going to have multiple sites in memory at once
			siteParam = (Parameter<?>)siteParam.clone();
			if (siteParam.getName().equals(Vs30_Param.NAME))
				// set Vs30 to 600 m/s
				((Parameter<Double>)siteParam).setValue(new Double(600d));
			site.addParameter(siteParam);
			System.out.println(siteParam.getName()+": "+siteParam.getValue());
		}
		
		// x-values for the hazard curve
		DiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		System.out.println("Default x-values:\n"+xValues);
		
		// need natural log x-values for curve calculation
		DiscretizedFunc logHazCurve = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xValues)
			logHazCurve.set(Math.log(pt.getX()), 1d); // y values don't matter yet
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		// default maximum source-site distance to include is 200 km, you can change it here if you want
//		calc.getAdjustableParams().getParameter(Double.class, MaxDistanceParam.NAME).setValue(300d);
		
		// Calculate the curve
		System.out.println("Calculating hazard curve");
		// this actually stores the y-values directly in logHazCurve
		calc.getHazardCurve(logHazCurve, site, gmpe, erf);
		
		// can convert back to linear if you want
		DiscretizedFunc linearHazCurve = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logHazCurve)
			linearHazCurve.set(Math.exp(pt.getX()), pt.getY());
		
		System.out.println("Final hazard curve:");
		System.out.println(linearHazCurve);
		
		// this was a 50 year curve, if you want to pull out the 2% value you can do this
		double imlAt2percent = linearHazCurve.getFirstInterpolatedX_inLogXLogYDomain(0.02);
		System.out.println("2% in "+(float)erf.getTimeSpan().getDuration()+" yr hazard: "+imlAt2percent);
	}

}
