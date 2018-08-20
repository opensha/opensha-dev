/**
 * 
 */
package scratch.aftershockStatistics;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import static org.opensha.commons.geo.GeoTools.TO_DEG;
import static org.opensha.commons.geo.GeoTools.TO_RAD;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

/**
 * @author field
 *
 */
public class AftershockStatsCalc {
	
	public final static double MILLISEC_PER_YEAR = 1000*60*60*24*365.25;
	public final static long MILLISEC_PER_DAY = 1000*60*60*24;




	/**
	 * Compute the function value (x^z - y^z)/z
	 * This function is designed to avoid numerical cancellation when z is small,
	 * and deliver a correct result when z is zero.
	 * Requires: x > 0 and y > 0.
	 *
	 * Implementation: The function uses the identity
	 *  (x^z - y^z)/z == (2/z)*(sqrt(x*y)^z)*sinh(z*log(sqrt(x/y)))
	 * It also uses the fact that if abs(v) < 1.0e-8 then sinh(v) equals v to numerical precision.
	 */
	public static double pow_diff_div (double x, double y, double z) {
		if (!( x > 0.0 && y > 0.0 )) {
			throw new RuntimeException("AftershockStatsCalc.pow_diff_div: Parameters x and y are not positive");
		}
		double sqrt_x = Math.sqrt(x);
		double sqrt_y = Math.sqrt(y);
		double u = Math.log(sqrt_x / sqrt_y);
		double v = z * u;
		if (Math.abs(v) < 1.0e-8) {
			return 2.0 * Math.pow(sqrt_x * sqrt_y, z) * u;
		}
		return 2.0 * Math.pow(sqrt_x * sqrt_y, z) * Math.sinh(v) / z;
	}



	
	/**
	 * This computes the log-likelihood for the given modified Omori parameters according to 
	 * equation (6) of Ogata (1983; J. Phys. Earth, 31,115-124).
	 * @param k = Omori k-parameter (amplitude).
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tMinDays - the start time of of the catalog in days from main shock
	 * @param tMaxDays - the end time of of the catalog in days from main shock
	 * @param relativeEventTimes - in order of occurrence relative to main shock,
	 *        this is an array whose i-th element is the time the i-th aftershock
	 *        occurred, measured in days since the mainshock.
	 * @return
	 * The return value is log(L(k,p)), where L(k,p) is the likelihood of the parameter
	 * values k and p, taking c to be fixed.
	 *  log(L(k,p)) = SUM(log(lambda(t_i))) - INTEGRAL(lambda(t)*dt)
	 * where
	 *  t_i = Time of the i-th aftershock, in days since the mainshock.
	 *  lamda(t) = k*((t + c)^(-p))     [Omori formula]
	 *  SUM runs over the all the aftershocks.
	 *  INTEGRAL runs over times t = tMinDays to t = tMaxDays.
	 */
	public static double getLogLikelihoodForOmoriParams(double k, double p, double c, double tMinDays, double tMaxDays, double[] relativeEventTimes) {
		double funcA=Double.NaN;
		int n=relativeEventTimes.length;

//		if(p == 1)
//			funcA = Math.log(tMaxDays+c) - Math.log(tMinDays+c);
//		else
//			funcA = (Math.pow(tMaxDays+c,1-p) - Math.pow(tMinDays+c,1-p)) / (1-p);
		funcA = pow_diff_div (tMaxDays+c, tMinDays+c, 1.0 - p);

		double sumLn_t = 0;
		for(double t : relativeEventTimes)
			sumLn_t += Math.log(t+c);
//double tempLL = n*Math.log(k) - p*sumLn_t - k*funcA;
//System.out.println(n*Math.log(k)+"\t"+(-p*sumLn_t)+"\t"+(-k*funcA)+"\t"+tempLL);
		return n*Math.log(k) - p*sumLn_t - k*funcA;
	}
	
	

	
	/**
	 * This computes the maximum likelihood k values for constrained 
	 * values of p and c as given.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tMinDays - the start time of of the catalog in days from main shock
	 * @param tMaxDays - the end time of of the catalog in days from main shock
	 * @param numEvents - the number of events in the catalog
	 * @return
	 * The return value is the maximum-likelihood value of k, for given Omori parameters
	 * p and c which are taken to be fixed, and for a given number of aftershocks.
	 */
	public static double getMaxLikelihood_k(double p, double c, double tMinDays, double tMaxDays, int numEvents) {
		double funcA=Double.NaN;

//		if(p == 1)
//			funcA = Math.log(tMaxDays+c) - Math.log(tMinDays+c);
//		else
//			funcA = (Math.pow(tMaxDays+c,1-p) - Math.pow(tMinDays+c,1-p)) / (1-p);
		funcA = pow_diff_div (tMaxDays+c, tMinDays+c, 1.0 - p);
		
// System.out.println("getMaxLikelihood_k: \t"+p+"\t"+c+"\t"+tMinDays+"\t"+tMaxDays+"\t"+funcA+"\t"+numEvents+"\t"+(numEvents/funcA));
		return (double)numEvents/funcA;
	}
	

	
	
	/**
	 * This converts the productivity value from "a" to "k"
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @return k = Omori k-parameter (amplitude).
	 *  k = 10^(a + b*(magMain - magMin))
	 */
	public static double convertProductivityTo_k(double a, double b, double magMain, double magMin) {
		return Math.pow(10.0, a+b*(magMain-magMin));
	}



	
	/**
	 * This converts the productivity value from "k" to "a"
	 * @param k = Omori k-parameter (amplitude).
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @return a = Reasenberg-Jones productivity parameter.
	 *  k = 10^(a + b*(magMain - magMin))
	 */
	public static double convertProductivityTo_a(double k, double b, double magMain, double magMin) {
		return Math.log10(k) - b*(magMain-magMin);
	}




	/**
	 * This returns the Reasenberg Jones (1989, 1994) expected number of primary aftershocks 
	 * between time tMinDays and tMaxDays (days after the mainshock) for  the given arguments.
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tMinDays = Beginning of forecast time window (since origin time), in days.
	 * @param tMaxDays = End of forecast time window (since origin time), in days.
	 * @return
	 * According to R&J, the rate of aftershocks of magnitude >= magMin is
	 *  lambda(t) = k * (t + c)^(-p)
	 * where
	 *  k = 10^(a + b*(magMain - magMin))
	 * The value returned by this function is the integral of lambda(t) from t=tMinDays to t=tMaxDays.
	 */
	public static double getExpectedNumEvents(double a, double b, double magMain, double magMin, double p, double c, double tMinDays, double tMaxDays) {
		double k = convertProductivityTo_k(a, b, magMain, magMin);

//		if(p!=1) {
//			double oneMinusP= 1-p;
//			return (k/oneMinusP)*(Math.pow(c+tMaxDays,oneMinusP) - Math.pow(c+tMinDays,oneMinusP));
//		}
//		else {
//			return k*(Math.log(c+tMaxDays) - Math.log(c+tMinDays));
//		}

		return k * pow_diff_div (c+tMaxDays, c+tMinDays, 1.0 - p);
	}




	/**
	 * This returns the Reasenberg Jones (1989, 1994) rate of expected number of primary aftershocks 
	 * at time tDays (days after the mainshock) for the given arguments.
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tDays = Time (since origin time), in days.
	 * @return
	 * According to R&J, the rate of aftershocks of magnitude >= magMin is
	 *  lambda(t) = k * (t + c)^(-p)
	 * where
	 *  k = 10^(a + b*(magMain - magMin))
	 * The value returned by this function is lambda(t) at t=tDays.
	 */
	public static double getExpectedEventsRate(double a, double b, double magMain, double magMin, double p, double c, double tDays) {
		double k = convertProductivityTo_k(a, b, magMain, magMin);
		return k * Math.pow(c+tDays, -p);
	}

	
	
	
	/**
	 * This returns the expected number of primary aftershocks as a function of time
	 * 
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tMin = Start of time range, in days after the mainshock.
	 * @param tMax = End of time range, in days after the mainshock.
	 * @param tDelta = Spacing between time values in the returned function.
	 * @return
	 * Returns a discrete function with:
	 *  x = Time (in days after mainshock).
	 *  y = Expected number of aftershocks within the corresponding time interval.
	 * The range [tMin, tMax] is partitioned into equal-sized intervals, with the
	 * width of each interval equal to approximately tDelta.  The interval width is
	 * adjusted so that a whole number of intervals fit within the range [tMin, tMax].
	 * For each such interval, the x value is the midpoint of the interval (and so
	 * the first and last x values are tMin+x_delta/2 and tMax-x_delta/2,
	 * where x_delta is the adjusted interval width).  The y value is the expected
	 * number of aftershocks within that interval according to the R&J formula.
	 */
	public static  EvenlyDiscretizedFunc getExpectedNumWithTimeFunc(double a, double b, double magMain, double magMin, double p, double c,  double tMin, double tMax, double tDelta) {
		
		// Get the function spacing

		int num = Math.max((int)Math.round((tMax-tMin)/tDelta), 1);
		double x_delta = (tMax-tMin)/num;
		double x_min;
		double x_max;
		if (num == 1) {
			x_min = (tMin + tMax)*0.5;
			x_max = x_min;	// x_max and x_min must be precisely equal for num==1, otherwise EvenlyDiscretizedFunc throws an exception
		} else {
			x_min = tMin+x_delta/2;
			x_max = tMax-x_delta/2;
		}

		// Construct the function
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(x_min, x_max, num);
		for(int i=0;i<func.size();i++) {
			double binTmin = func.getX(i) - x_delta/2;
			double binTmax = func.getX(i) + x_delta/2;
			double yVal = getExpectedNumEvents(a, b, magMain, magMin, p, c, binTmin, binTmax);
			func.set(i,yVal);
		}

		// Insert info

		func.setName("Expected Number of Primary Aftershocks for "+x_delta+"-day intervals");
		func.setInfo("for a="+a+", b="+b+", p="+p+", c="+c+", magMain="+magMain+", magMin="+magMin);
		return func;
	}
	
	
	/**
	 * This returns the expected cumulative number of primary aftershocks as a function of time
	 * 
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset).
	 * @param tMin = Start of time range, in days after the mainshock.
	 * @param tMax = End of time range, in days after the mainshock.
	 * @param tDelta = Spacing between time values in the returned function.
	 * @return
	 * Returns a discrete function with:
	 *  x = Time (in days after mainshock).
	 *  y = Cumulative expected number of aftershocks for the corresponding time interval.
	 * The range [tMin, tMax] is partitioned into equal-sized intervals, with the
	 * width of each interval equal to approximately tDelta.  The interval width is
	 * adjusted so that a whole number of intervals fit within the range [tMin, tMax].
	 * For each such interval, the x value is the midpoint of the interval (and so
	 * the first and last x values are tMin+x_delta/2 and tMax-x_delta/2,
	 * where x_delta is the adjusted interval width).  The y value is the expected
	 * number of aftershocks within that interval plus all preceding intervals according
	 * to the R&J formula.
	 */
	public static  EvenlyDiscretizedFunc getExpectedCumulativeNumWithTimeFunc(double a, double b, double magMain, double magMin, double p, double c,  double tMin, double tMax, double tDelta) {
		
		// Get the function spacing

		int num = Math.max((int)Math.round((tMax-tMin)/tDelta), 1);
		double x_delta = (tMax-tMin)/num;
		double x_min;
		double x_max;
		if (num == 1) {
			x_min = (tMin + tMax)*0.5;
			x_max = x_min;	// x_max and x_min must be precisely equal for num==1, otherwise EvenlyDiscretizedFunc throws an exception
		} else {
			x_min = tMin+x_delta/2;
			x_max = tMax-x_delta/2;
		}

		// Construct the function
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(x_min, x_max, num);
		for(int i=0;i<func.size();i++) {
			double binTmin = func.getX(0) - x_delta/2;		// this line is the only change from getExpectedNumWithTimeFunc
			double binTmax = func.getX(i) + x_delta/2;
			double yVal = getExpectedNumEvents(a, b, magMain, magMin, p, c, binTmin, binTmax);
			func.set(i,yVal);
		}

		return func;
	}



	
	/**
	 * This returns the poisson probability given the expected number of events
	 * @param expectedNum
	 * @return
	 */
	public static double getPoissonProbability(double expectedNum) {
		return 1.0-Math.exp(-expectedNum);
	}



	
	/**
	 * This returns the maximum-likelihood b-value defined by Aki (1965, Bull. Earthq. Res. Inst., 43, 237-239)
	 * @param rups - obs eqk rupture list
	 * @param magComplete - the magnitude above which no events have gone undetected
	 * @param magPrecision - the degree to which magnitude have been rounded
	 * @return
	 */
	public static double getMaxLikelihood_b_value(List<ObsEqkRupture> rups, double magComplete,
			double magPrecision) {
		double magMean = 0d;
		int num = 0;
		for (ObsEqkRupture rup : rups) {
			if (rup.getMag() >= magComplete) {
				num++;
				magMean += rup.getMag();
			}
		}
		Preconditions.checkState(num > 0, "No ruptures above mc="+magComplete);
		magMean /= (double)num;
		return getMaxLikelihood_b_value(magMean, magComplete, magPrecision);
	}



	
	/**
	 * This does not check for negative values
	 * @param mainShock
	 * @param aftershockList
	 * @return
	 * Given a mainshock and a list of aftershocks, this function returns an array
	 * whose i-th element contains the elapsed time from the mainshock to the i-th
	 * aftershock, in days.
	 */
	public static double[] getDaysSinceMainShockArray(ObsEqkRupture mainShock, List<ObsEqkRupture> aftershockList) {
		double[] relativeEventTimesDays = new double[aftershockList.size()];
		for(int i=0; i<aftershockList.size();i++) {
			long epochDiff = aftershockList.get(i).getOriginTime()-mainShock.getOriginTime();
			relativeEventTimesDays[i] = (double)(epochDiff) / (double)MILLISEC_PER_DAY;
		}
		return relativeEventTimesDays;
	}



	
	/**
	 * This returns the maximum-likelihood b-value defined by Aki (1965, Bull. Earthq. Res. Inst., 43, 237-239)
	 * @param magMean - mean magnitude above magComplete
	 * @param magComplete - the magnitude above which no events have gone undetected
	 * @param magPrecision - the degree to which magnitude have been rounded
	 * @return
	 */
	public static double getMaxLikelihood_b_value(double magMean, double magComplete, double magPrecision) {
		return Math.log10(Math.E) /(magMean - (magComplete-0.5*magPrecision));
	}




	/**
	 * This returns the Page et al. (2016) time-dependent magnitude of completeness 
	 * at time tDays (days after the mainshock) for the given arguments.
	 * @param magMain = Magnitude of mainshock.
	 * @param magCat = Magnitude of completeness when there has not been a mainshock.
	 * @param capG = The "G" parameter in the time-dependent magnitude of completeness model. 
	 *               As a special case, if capG == 10.0 then the return value is always magCat.
	 * @param capH = The "H" parameter in the time-dependent magnitude of completeness model.
	 * @param tDays = Time (since origin time), in days.
	 * @return
	 * According to Page et al. the magnitude of completeness is
	 *  magMin(t) = Max(magMain/2 - G - H*log10(t), magCat)
	 * where t is in days.
	 * This formula reflects the fact that in the aftermath of a mainshock, the catalog's
	 * magnitude of completeness temporarily increases, and then gradually returns to normal.
	 * The value returned by this function is magMin(t) at t=tDays.
	 * Note: The magnitude of completeness is the minimum magnitude at which the catalog contains
	 * essentially all earthquakes (which Page et al. defines as 95% of earthquakes).
	 */
	public static double getPageMagCompleteness(double magMain, double magCat, double capG, double capH, double tDays) {
		if (capG > 9.999) {
			return magCat;
		}
		return Math.max (magCat, 0.5 * magMain - capG - capH * Math.log10 (Math.max (tDays, Double.MIN_NORMAL)));
	}




	/**
	 * This returns the Page et al. (2016) time of completeness, in days since the mainshock.
	 * @param magMain = Magnitude of mainshock.
	 * @param magCat = Magnitude of completeness when there has not been a mainshock.
	 * @param capG = The "G" parameter in the time-dependent magnitude of completeness model. 
	 *               As a special case, if capG == 10.0 then the return value is always magCat.
	 * @param capH = The "H" parameter in the time-dependent magnitude of completeness model.
	 * @return
	 * This function returns the time (in days since the mainshock) at which the
	 * magnitude of completness equals the catalog's normal magnitude of completeness.
	 * According to Page et al. the magnitude of completeness is
	 *  magMin(t) = Max(magMain/2 - G - H*log10(t), magCat)
	 * where t is in days.
	 * This function returns the value of t at which magMin(t) == magCat, which is
	 *  t = 10^((magMain/2 - G - magCat)/H)
	 * However this must be evaluated carefully to avoid overflow or divide-by-zero.
	 */
	public static double getPageTimeOfCompleteness(double magMain, double magCat, double capG, double capH) {
		if (capG > 9.999) {
			return 0.0;
		}

		if (!( capH >= 0.0 )) {
			throw new RuntimeException("AftershockStatsCalc.getPageTimeOfCompleteness: H parameter is negative");
		}

		double x = 0.5 * magMain - capG - magCat;

		if (x <= -8.0 * capH) {
			return 0.0;			// less than about 1 millisecond, just return 0
		}
		if (x >= 12.0 * capH) {
			return 1.0e12;		// more than about 3 billion years
		}

		return Math.pow(10.0, x/capH);
	}




	/**
	 * This returns the Reasenberg Jones (1989, 1994) rate of expected number of primary aftershocks 
	 * at time tDays (days after the mainshock) for the given arguments,
	 * with the Page et al. (2016) time-dependent magnitude of completeness.
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magCat = Magnitude of completeness when there has not been a mainshock.
	 * @param capG = The "G" parameter in the time-dependent magnitude of completeness model. 
	 *               As a special case, if capG == 10.0 then the return value is always magCat.
	 * @param capH = The "H" parameter in the time-dependent magnitude of completeness model.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tDays = Time (since origin time), in days.
	 * @return
	 * According to R&J, the rate of aftershocks of magnitude >= magMin is
	 *  lambda(t) = k * (t + c)^(-p)
	 * where
	 *  k = 10^(a + b*(magMain - magMin))
	 * According to Page et al. the magnitude of completeness is
	 *  magMin(t) = Max(magMain/2 - G - H*log10(t), magCat)
	 * In these formulas, t is measured in days.
	 * The value returned by this function is lambda(t) at t=tDays.
	 */
	public static double getPageExpectedEventsRate(double a, double b, double magMain, double magCat, double capG, double capH, double p, double c, double tDays) {
		double magMin = getPageMagCompleteness(magMain, magCat, capG, capH, tDays);
		return getExpectedEventsRate(a, b, magMain, magMin, p, c, tDays);
	}




	/**
	 * This returns the Reasenberg Jones (1989, 1994) expected number of primary aftershocks 
	 * between time tMinDays and tMaxDays (days after the mainshock) for the given arguments,
	 * with the Page et al. (2016) time-dependent magnitude of completeness.
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magCat = Magnitude of completeness when there has not been a mainshock.
	 * @param capG = The "G" parameter in the time-dependent magnitude of completeness model. 
	 *               As a special case, if capG == 10.0 then the return value is always magCat.
	 * @param capH = The "H" parameter in the time-dependent magnitude of completeness model.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tMinDays = Beginning of forecast time window (since origin time), in days.
	 * @param tMaxDays = End of forecast time window (since origin time), in days.
	 * @return
	 * According to R&J, the rate of aftershocks of magnitude >= magMin is
	 *  lambda(t) = k * (t + c)^(-p)
	 * where
	 *  k = 10^(a + b*(magMain - magMin))
	 * According to Page et al. the magnitude of completeness is
	 *  magMin(t) = Max(magMain/2 - G - H*log10(t), magCat)
	 * In these formulas, t is measured in days.
	 * The value returned by this function is the integral of lambda(t) from t=tMinDays to t=tMaxDays.
	 *
	 * Implementation note: This function uses numerical integration for times t < tPage,
	 * and an analytic formula for times t > tPage, where tPage is the time when the
	 * magnitude of completeness first becomes equal to magCat.  Aside from the gain in
	 * efficiency, it is necessary to break the domain of integration into two parts because
	 * the rate function is non-differentiable at t = tPage.
	 *
	 * Note: It is possible to force the use of numeric integration for an RJ distribution
	 * with constant magMin by choosing parameters so that:
	 *  magMain/2 - G == magMin
	 *  H == 0
	 *  magCat < magMin
	 * This is useful for testing the numeric integration code, by comparing to the analytic formula.
	 */
	public static double getPageExpectedNumEvents(double a, double b, double magMain, double magCat, double capG, double capH, double p, double c, double tMinDays, double tMaxDays) {
		
		// Transition time, when magnitude of completeness first becomes equal to magCat

		double tPage = getPageTimeOfCompleteness(magMain, magCat, capG, capH);
		
		// Integral value

		double s = 0.0;

		// Numeric integration for times before tPage

		if (tPage > tMinDays) {

			// Set up functional object
		
			funcExpectedEventsRate func = new funcExpectedEventsRate (a, b, magMain, magCat, capG, capH, p, c);

			// Upper limit of integration

			double tUpper = Math.min(tMaxDays, tPage);

			// Force at least 30 points to be sampled, but don't force steps smaller than about 1 second

			double max_h = Math.max(1.0e-5, (tUpper - tMinDays) / 30.0);

			// Error tolerances

			double abs_tol = 0.0;
			double rel_tol = 1.0e-7;

			// Do the integration

			s += adapQuadSimpson (func, tMinDays, tUpper, abs_tol, rel_tol, max_h);
		}

		// Analytic formula for times after tPage

		if (tPage < tMaxDays) {
			s += getExpectedNumEvents(a, b, magMain, magCat, p, c, Math.max(tMinDays, tPage), tMaxDays);
		}

		return s;
	}




	/**
	 * This produces a randomly-generated simulated aftershock sequence, following the
	 * Reasenberg-Jones (1989, 1994) statistics formula, with the Page et al. (2016)
	 * time-dependent magnitude of completeness.
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magCat = Magnitude of completeness when there has not been a mainshock.
	 * @param capG = The "G" parameter in the time-dependent magnitude of completeness model. 
	 *               As a special case, if capG == 10.0 then the magnitude of completeness is always magCat.
	 * @param capH = The "H" parameter in the time-dependent magnitude of completeness model.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * @param tMinDays = Beginning of time span (since origin time), in days.
	 * @param tMaxDays = End of time span (since origin time), in days.
	 * @return
	 * Returns a randomly-generated sequence of aftershocks.
	 * Each entry in the returned list contains a time (in milliseconds since the mainshock)
	 * and magnitude; no other information is placed in the list.
	 *
	 * Implementation notes:
	 * The given time span is paritioned into intervals where there is about 3 expected aftershock
	 * (according to the R&J formula, with a lower bound for time-dependent magnitude of completeness).
	 * In each interval, a Poisson distribution is used to randomly select the number of aftershocks.
	 * Then, a rejection sampling technique is used to randomly select a time within the interval
	 * for each aftershock, according to the R&J density function.
	 * Then, an exponential distribution is used to randomly select a magnitude for each aftershock.
	 * Then, any aftershocks below the time-dependent magnitude of completeness are discarded.
	 *
	 * The documentation for UniformRealDistribution is here:
	 * http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/distribution/UniformRealDistribution.html
	 * The algorithm used for generating the Poisson distribution is the first algorithm given here:
	 * https://en.wikipedia.org/wiki/Poisson_distribution
	 * The algorithm used for generating the exponential distribution is given here:
	 * https://en.wikipedia.org/wiki/Exponential_distribution
	 */
	public static ObsEqkRupList simAftershockSequence(double a, double b, double magMain, double magCat, double capG, double capH, double p, double c, double tMinDays, double tMaxDays) {
		if (!( b > 0.0 )) {
			throw new RuntimeException("AftershockStatsCalc.simAftershockSequence: b parameter is negative or zero");
		}
		if (!( capH >= 0.0 )) {
			throw new RuntimeException("AftershockStatsCalc.simAftershockSequence: H parameter is negative");
		}
		if (!( p > 0.0 )) {
			throw new RuntimeException("AftershockStatsCalc.simAftershockSequence: p parameter is negative or zero");
		}
		if (!( c >= 0.0 )) {
			throw new RuntimeException("AftershockStatsCalc.simAftershockSequence: c parameter is negative");
		}
		if (!( tMinDays + c > 0.0 && tMinDays < tMaxDays )) {
			throw new RuntimeException("AftershockStatsCalc.simAftershockSequence: invalid time span");
		}

		// Maximum number of aftershocks per interval that we allow

		final int max_aftershocks = 100;

		// Array to hold the time for each aftershock within an interval

		double[] t_aftershock = new double[max_aftershocks];

		// Random number generator, produces random numbers between 0.0 (inclusive) and 1.0 (exclusive)

		UniformRealDistribution rangen = new UniformRealDistribution();

		// The list of aftershocks

		ObsEqkRupList aftershock_list = new ObsEqkRupList();

		// The minimum magnitude we need to consider is the magnitude of completeness at the end of the time span,
		// which is a lower bound for magnitude of completeness throughout the time span

		double magMin = getPageMagCompleteness(magMain, magCat, capG, capH, tMaxDays);

		// The start time of the current interval, in days

		double t_now = tMinDays;

		// Flag used to control interval generation

		boolean f_continue = true;			// true if there are more intervals to do

		// Loop until all intervals are done

		while (f_continue) {

			// The current aftershock rate is an upper bound for aftershock rate in the interval

			double rate_now = getExpectedEventsRate(a, b, magMain, magMin, p, c, t_now);

			// Get the approximate time interval in which 3 aftershocks are expected
			// (this is an upper bound because the rate is decreasing)

			double t_delta = 3.0 / rate_now;

			// If it extends past the end of the time span, clip to end of time span and make it the last interval

			if (t_now + t_delta >= tMaxDays) {
				t_delta = tMaxDays - t_now;
				f_continue = false;
			}

			// If it extends almost to the end of the time span, go halfway to the end

			else if (t_now + t_delta*1.5 >= tMaxDays) {
				t_delta = (tMaxDays - t_now) * 0.5;
			}

			// The minimum magnitude we need to consider is the magnitude of completeness at the end of the interval,
			// which is a lower bound for magnitude of completeness throughout the interval

			double magMinInt = getPageMagCompleteness(magMain, magCat, capG, capH, t_now + t_delta);

			// Get the expected number of aftershocks in the interval from t_now to t_now + t_delta

			double expected_aftershocks = getExpectedNumEvents(a, b, magMain, magMinInt, p, c, t_now, t_now + t_delta);

			// If the expected number of aftershocks is less than 0.75, double the interval until it's larger

			if (f_continue) {
				while (expected_aftershocks < 0.75 && t_now + t_delta*4.0 < tMaxDays) {
					double new_t_delta = t_delta * 2.0;
					double new_magMinInt = getPageMagCompleteness(magMain, magCat, capG, capH, t_now + new_t_delta);
					double new_expected_aftershocks = getExpectedNumEvents(a, b, magMain, new_magMinInt, p, c, t_now, t_now + new_t_delta);

					if (new_expected_aftershocks > 3.0) {
						break;
					}

					t_delta = new_t_delta;
					magMinInt = new_magMinInt;
					expected_aftershocks = new_expected_aftershocks;
				}
			}

			// Apply the Poisson distribution to select the actual number of aftershocks

			int actual_aftershocks = -1;
			double pd_l = Math.exp (-expected_aftershocks);
			double pd_p = 1.0;

			do {
				++actual_aftershocks;
				pd_p *= rangen.sample();
			} while (pd_p > pd_l && actual_aftershocks < max_aftershocks - 1);

			// If there are aftershocks in this interval ...

			if (actual_aftershocks > 0) {

				// The current aftershock rate is an upper bound for aftershock rate in the interval

				double rate_ub = getExpectedEventsRate(a, b, magMain, magMinInt, p, c, t_now);

				// Loop over aftershocks within this interval

				for (int i = 0; i < actual_aftershocks; ++i) {

					// Use rejection sampling technique to select a time for this aftershock.
					// This works by sampling points uniformly in the rectangle t_now <= t <= t_now + t_delta
					// and 0 <= h <= rate_ub, then rejecting those that lie above the R&J probability density.

					double t;
					double h;
					double r;

					do {
						t = t_now + t_delta * rangen.sample();
						h = rate_ub * rangen.sample();
						r = getExpectedEventsRate(a, b, magMain, magMinInt, p, c, t);
					} while (h > r);

					t_aftershock[i] = t;
				}

				// Sort the aftershocks into temporal order

				Arrays.sort (t_aftershock, 0, actual_aftershocks);

				// Loop over aftershocks within this interval

				for (int i = 0; i < actual_aftershocks; ++i) {

					// Time of this aftershock

					double t = t_aftershock[i];

					// Use the exponential distribution to get the magnitude

					double u = 1.0 - rangen.sample();
					double mag = magMinInt - Math.log10(Math.max(u, Double.MIN_NORMAL)) / b;

					// If the magnitude is at least the magnitude of completeness ...

					if (mag >= getPageMagCompleteness(magMain, magCat, capG, capH, t)) {
					
						// Add the aftershock to the list

						double timeMillis = t*((double)MILLISEC_PER_DAY);
						int eventId = aftershock_list.size() + 1;
						aftershock_list.add(new ObsEqkRupture(Integer.toString(eventId), (long)timeMillis, null, mag));
					}
				}
			}

			// Advance to next time interval

			t_now += t_delta;
		}

		// Return the resulting list of aftershocks

		return aftershock_list;
	}
	
	
	public static double[] readAndysFile() {
		
		try {
			BufferedReader buffRead = new BufferedReader(new InputStreamReader(
					AftershockStatsCalc.class.getResourceAsStream("AndysSimulationData.txt")));
			ArrayList<Double> eventTimeList = new ArrayList<Double>();
			String line = buffRead.readLine();
			while (line != null) {
				StringTokenizer tok = new StringTokenizer(line);
				while(tok.hasMoreTokens()) {
					eventTimeList.add(Double.parseDouble(tok.nextToken()));
				}
				line = buffRead.readLine();
			}
			double[] eventTimeArray = new double[eventTimeList.size()];
			for(int i=0;i<eventTimeList.size();i++) {
				eventTimeArray[i] = eventTimeList.get(i);
//				System.out.println(eventTimeArray[i]);
			}
			buffRead.close();
			return eventTimeArray;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private static ObsEqkRupList readJeannesFile() {
		
		try {
			BufferedReader buffRead = new BufferedReader(new InputStreamReader(
					AftershockStatsCalc.class.getResourceAsStream("JeannesSimulationData.txt")));
			ObsEqkRupList aftershockList = new ObsEqkRupList();
			String line = buffRead.readLine();
			int eventId = 0;
			while (line != null) {
				StringTokenizer tok = new StringTokenizer(line);
				double mag = Double.parseDouble(tok.nextToken());
				double timeMillis = Double.parseDouble(tok.nextToken())*(double)MILLISEC_PER_DAY;
				if(eventId != 0) // skip main shock
					aftershockList.add(new ObsEqkRupture(Integer.toString(eventId), (long)timeMillis, null, mag));
				eventId+=1;
				line = buffRead.readLine();
			}
			buffRead.close();
			return aftershockList;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	
	
	public static void plot2D_PDF(EvenlyDiscrXYZ_DataSet pdf2D, String title,
			String xAxisLabel, String yAxisLabel, String zAxisLabel) {
		CPT cpt=null;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(pdf2D.getMinZ(), pdf2D.getMaxZ());
		} catch (IOException e) {
			e.printStackTrace();
		}
		XYZPlotSpec logLikeSpec = new XYZPlotSpec(pdf2D, cpt, title, xAxisLabel, yAxisLabel, zAxisLabel);
		XYZPlotWindow window_logLikeSpec = new XYZPlotWindow(logLikeSpec, new Range(pdf2D.getMinX(),pdf2D.getMaxX()), new Range(pdf2D.getMinY(),pdf2D.getMaxY()));
	}
	
	
	
	public static void testJeanneCalc() {
		
//		Here's a synthetic dataset, with magnitude and time following a M7.5 mainshock.  
//		The synthetics were created with the parameters:  a=-1.67, b_in=0.91, c_in=0.05 days, 
//		and p_in=1.08.  The completeness is described by Mcat=2.5, G=4.5, and H=0.75.  
//		I solved just for a and p, fixing everything else to the correct value, and found 
//		a=-1.69, and p=1.05.  This will force you to search in the vicinity of p=1.


		ObsEqkRupture mainShock = new ObsEqkRupture("0", 0l,null, 7.5);
		ObsEqkRupList aftershockList = readJeannesFile();
		System.out.println("Num aShocks = "+aftershockList.size());
		double magCat = 2.5;
		double capG=4.5;
		double capH=0.75;
		double b=0.91;
		double dataStartTimeDays=0;
		double dataEndTimeDays=30;
		
		double min_a = -2.0;
		double max_a = -1.0;
		int num_a = 101;

		double min_p = 0.9; 
		double max_p = 1.2; 
		int num_p = 31;
//		double min_p = 1.0; 
//		double max_p = 1.0; 
//		int num_p = 1;
		
		double min_c=0.05;
		double max_c=0.05;
		int num_c=1;

		
//		ReasenbergJonesAftershockModel solution = new ReasenbergJonesAftershockModel(mainShock, aftershockList, magCat, capG, capH, b, dataStartTimeDays, dataEndTimeDays,
//				min_a, max_a, num_a, min_p, max_p, num_p, min_c, max_c, num_c);
//		
//		plot2D_PDF(solution.get2D_PDF_for_a_and_p(), "PDF for a vs p", "a", "p", "density");
//		
//		GraphWindow graph = new GraphWindow(solution.getPDF_a(), "a-value PDF"); 
//		graph.setX_AxisLabel("a-axis");
//		graph.setY_AxisLabel("DensityF");
////		graph.setX_AxisRange(-4, 3);
////		graph.setY_AxisRange(1e-3, graph.getY_AxisRange().getUpperBound());
//		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLACK));
//		graph.setPlotChars(plotChars);
//		graph.setPlotLabelFontSize(18);
//		graph.setAxisLabelFontSize(16);
//		graph.setTickLabelFontSize(14);
//
//		
//		GraphWindow graph2 = new GraphWindow(solution.getPDF_p(), "p-value PDF"); 
//		graph2.setX_AxisLabel("p-axis");
//		graph2.setY_AxisLabel("DensityF");
//		graph2.setPlotChars(plotChars);
//		graph2.setPlotLabelFontSize(18);
//		graph2.setAxisLabelFontSize(16);
//		graph2.setTickLabelFontSize(14);


	}

	
	
	
	
	
	public static double getMmaxC(IncrementalMagFreqDist mfd) {
		List<Double> magsAtMax = Lists.newArrayList();
		double max = 0d;
		
		for (Point2D pt : mfd) {
			if (pt.getY() == max)
				magsAtMax.add(pt.getX());
			else if (pt.getY() > max) {
				// start over
				magsAtMax = Lists.newArrayList(pt.getX());
				max = pt.getY();
			}
		}
		
		double mmaxc = 0;
		for (double mag : magsAtMax)
			mmaxc += mag;
		mmaxc /= (double)magsAtMax.size();
		System.out.println("Mmaxc="+(float)mmaxc+" from MFD mode(s): "+Joiner.on(",").join(magsAtMax));
		return mmaxc;
	}





	/**
	 * Functional object to calculate R&J aftershock rate with time-dependent magnitude of completeness.
	 * This object stores all its parameters.
	 * Then, a call to the value(t) function returns the aftershock rate at time t (days since the mainshock).
	 * @param a = Reasenberg-Jones productivity parameter.
	 * @param b = Gutenberg-Richter b-parameter.
	 * @param magMain = Magnitude of mainshock.
	 * @param magCat = Magnitude of completeness when there has not been a mainshock.
	 * @param capG = The "G" parameter in the time-dependent magnitude of completeness model. 
	 *               As a special case, if capG == 10.0 then the magnitude of completeness is always magCat.
	 * @param capH = The "H" parameter in the time-dependent magnitude of completeness model.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 *
	 * The documentation for UnivariateFunction is here:
	 * http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/analysis/UnivariateFunction.html
	 */
	public static class funcExpectedEventsRate implements UnivariateFunction {

		// Parameters for R&J and magnitude of completeness

		public double a;
		public double b;
		public double magMain;
		public double magCat;
		public double capG;
		public double capH;
		public double p;
		public double c;

		// Constructor saves all the parameters.

		public funcExpectedEventsRate (double a, double b, double magMain, double magCat, double capG, double capH, double p, double c) {
			this.a = a;
			this.b = b;
			this.magMain = magMain;
			this.magCat = magCat;
			this.capG = capG;
			this.capH = capH;
			this.p = p;
			this.c = c;
		}

		// Get the aftershock rate at time tDays (in days since the mainshock).

		@Override
		public double value(double tDays) {
			return AftershockStatsCalc.getPageExpectedEventsRate(a, b, magMain, magCat, capG, capH, p, c, tDays);
		}
	}



	
	/**
	 * Recursive subroutine for adaptive quadrature using Simpson's rule.
	 * This implements a 5-point Simpson's rule at equally spaced points x0,x1,x2,x3,x4.
	 * It estimates error by comparing to the 3-point rule for x0,x2,x4.
	 * If error exceeds tolerance, it calls itself recursively for the two subintervals.
	 * @param func = Function being integrated.
	 * @param x0 = Abscissa 0.
	 * @param x4 = Abscissa 4.
	 * @param f0 = Function value for x0.
	 * @param f2 = Function value for x2.
	 * @param f4 = Function value for x4.
	 * @param abs_tol = Absolute error tolerance (must be >= 0.0).
	 * @param rel_tol = Relative error tolerance (must be >= 0.0).
	 * @param min_level = Minimum level of recursion before accepting results (must be <= max_level).
	 * @param max_level = Maximum allowed level of recursion.
	 * @return
	 * Returns the integral of the function from x0 to x4.
	 */
	private static double adapQuadSimpsonRecur (UnivariateFunction func,
			double x0, double x4, double f0, double f2, double f4,
			double abs_tol, double rel_tol, int min_level, int max_level) {

		// Get the function values for x1 and x3

		double x2 = (x0 + x4) * 0.5;
		double x1 = (x0 + x2) * 0.5;
		double x3 = (x2 + x4) * 0.5;

		double f1 = func.value(x1);
		double f3 = func.value(x3);

		// If satisfies minimum level ...

		if (min_level <= 0) {

			// Calculate the 3-point and 5-point Simpson's rules

			double h = x4 - x0;

			double s3pt = (f0 + 4.0*f2 + f4) * h/6.0;
			double s5pt = (f0 + 4.0*f1 + 2.0*f2 + 4.0*f3 + f4) * h/12.0;

			double sdelta = s5pt - s3pt;

			// If we're at maximum level, return the result

			if (max_level <= 0) {
				return s5pt + (sdelta / 15.0);
			}

			// Estimate the magnitude of the integral

			double smag = (Math.abs(f0) + 4.0*Math.abs(f1) + 2.0*Math.abs(f2) + 4.0*Math.abs(f3) + Math.abs(f4)) * h/12.0;

			// If error criterion is satisfied, return the result

			if (Math.abs(sdelta) <= 15.0 * (Math.max(abs_tol, rel_tol*smag) + Double.MIN_NORMAL)) {
				return s5pt + (sdelta / 15.0);
			}
		}

		// Recursive invocation

		return adapQuadSimpsonRecur (func, x0, x2, f0, f1, f2, abs_tol * 0.5, rel_tol, min_level - 1, max_level - 1)
		     + adapQuadSimpsonRecur (func, x2, x4, f2, f3, f4, abs_tol * 0.5, rel_tol, min_level - 1, max_level - 1);
	}



	
	/**
	 * Adaptive quadrature using Simpson's rule.
	 * Calculate the integral of the function from a to b.
	 * @param func = Function being integrated.
	 * @param a = Lower limit of range of integration.
	 * @param b = Upper limit of range of integration.
	 * @param abs_tol = Absolute error tolerance.
	 *   Attempt to make the absolute error less that abs_tol.
	 * @param rel_tol = Relative error tolerance.
	 *   Attempt to make the relative error less than rel_tol.
	 *   If rel_tol and abs_tol are both nonzero then the more lenient one is used.
	 * @param max_h = Maximum h for subintervals (h = length of subinterval) .
	 *   This can be used to ensure that sufficient samples are taken.
	 *   It is not an error for max_h to be larger than b - a.
	 *   As a special case, max_h == 0 means there is no maximum.
	 * @return
	 * Returns the integral of the function from a to b.
	 */
	private static double adapQuadSimpson (UnivariateFunction func,
			double a, double b, double abs_tol, double rel_tol, double max_h) {

		// Check for zero length interval

		double h = b - a;

		if (Math.abs(h) <= Double.MIN_NORMAL) {
			return 0.0;
		}

		// Get the initial function values

		double x2 = (a + b) * 0.5;

		double f0 = func.value(a);
		double f2 = func.value(x2);
		double f4 = func.value(b);

		// Choose a maximum level to retain about 20 significant bits in h (2.10e6 ~ 2^21, 8.59e9 ~ 2^33)

		double min_h = Math.max(Double.MIN_NORMAL * 2.10e6, Math.max(Math.abs(a), Math.abs(b)) / 8.59e9);

		int max_level = 0;
		double trial_h = Math.abs(b - a) * 0.5;
		while (max_level < 28 && trial_h > min_h) {
			++max_level;
			trial_h *= 0.5;
		}

		// Choose a minimum level with sufficient samples

		int min_level = 0;
		if (Math.abs(max_h) > Double.MIN_NORMAL) {
			trial_h = Math.abs(b - a) * 0.25;
			while (min_level < max_level && trial_h > Math.abs(max_h)) {
				++min_level;
				trial_h *= 0.5;
			}
		}

		// Invoke recursive adaptive routine

		return adapQuadSimpsonRecur (func, a, b, f0, f2, f4, Math.abs(abs_tol), Math.abs(rel_tol), min_level, max_level);
	}


	

	// [DEPRECATED]
    public static Location getCentroid(ObsEqkRupture mainshock, List<ObsEqkRupture> aftershocks) {
		// now works across prime meridian
		List<Location> locs = Lists.newArrayList(mainshock.getHypocenterLocation());
		for (ObsEqkRupture aftershock : aftershocks)
			locs.add(aftershock.getHypocenterLocation());
		List<Double> lats = Lists.newArrayList();
		List<Double> lons = Lists.newArrayList();
		for (Location loc : locs) {
			lats.add(loc.getLatitude());
			lons.add(loc.getLongitude());
		}
		double lat = FaultUtils.getAngleAverage(lats);
		while (lat > 90)
			lat -= 360;
		double lon = FaultUtils.getAngleAverage(lons);
//		System.out.println("Mainshock loc: "+mainshock.getHypocenterLocation());
//		System.out.println("Orig centroid lon: "+lon);
		while (lon > 180)
			lon -= 360;
		// now make sure longitude is in the same domain as the input event
		if (Math.abs(lon - mainshock.getHypocenterLocation().getLongitude()) > 270)
			lon += 360;
		Location centroid = new Location(lat, lon);
		double dist = LocationUtils.horzDistanceFast(mainshock.getHypocenterLocation(), centroid);
		System.out.println("Centroid: "+(float)lat+", "+(float)lon+" ("+(float)dist+" km from epicenter)");
		return centroid;
	}


	

	// Calculate the centroid of an aftershock sequence, using spherical geometry.
	// Note: This routine implicitly assumes that the mainshock itself is not in the list of aftershocks.

    public static Location getSphCentroid(ObsEqkRupture mainshock, List<ObsEqkRupture> aftershocks) {
		
		// Convert spherical to rectangular coordinates, and sum the unit vectors

		double x = 0.0;
		double y = 0.0;
		double z = 0.0;

		double lat;
		double lon;

		lat = mainshock.getHypocenterLocation().getLatRad();
		lon = mainshock.getHypocenterLocation().getLonRad();

		x += (Math.cos(lat) * Math.cos(lon));
		y += (Math.cos(lat) * Math.sin(lon));
		z += Math.sin(lat);

		for (ObsEqkRupture aftershock : aftershocks) {
			lat = aftershock.getHypocenterLocation().getLatRad();
			lon = aftershock.getHypocenterLocation().getLonRad();

			x += (Math.cos(lat) * Math.cos(lon));
			y += (Math.cos(lat) * Math.sin(lon));
			z += Math.sin(lat);
		}

		// If the resulting vector is very small, just return the mainshock location

		if (x*x + y*y + z*z < 1.0e-4) {
			lat = mainshock.getHypocenterLocation().getLatitude();
			lon = mainshock.getHypocenterLocation().getLongitude();
			if (lon > 180.0) {
				lon -= 360.0;
			}
		}

		// Otherwise, convert rectangular to spherical coordinates

		else {
			lat = Math.atan2 (z, Math.hypot(x, y)) * TO_DEG;
			lon = Math.atan2 (y, x) * TO_DEG;
		}

		// Make sure the angles are in range, since they were converted from radians

		if (lat > 90.0) {lat = 90.0;}
		if (lat < -90.0) {lat = -90.0;}
		if (lon > 180.0) {lon = 180.0;}
		if (lon < -180.0) {lon = -180.0;}

		// Centroid

		Location centroid = new Location (lat, lon);
		double dist = LocationUtils.horzDistance (mainshock.getHypocenterLocation(), centroid);
		System.out.println (String.format ("Centroid: %.5f, %.5f (%.3f km from epicenter)", lat, lon, dist));
		return centroid;
	}




	/**
	 * @param args
	 */
	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("AftershockStatsCalc : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1
		// This is the pre-existing test for the class.

		if (args[0].equalsIgnoreCase ("test1")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("AftershockStatsCalc : Invalid 'test1' subcommand");
				return;
			}

			// Run the test
		
//			double mean = 0.78;
			double mean = 11.5079;
		
			PoissonDistribution poissDist = new PoissonDistribution(mean);
			System.out.println(poissDist.inverseCumulativeProbability(0.025));
			System.out.println(poissDist.inverseCumulativeProbability(0.975));
		
			int minInt=0;
			int maxInt=poissDist.inverseCumulativeProbability(0.999);
		
			HistogramFunction hist = new HistogramFunction((double)minInt,(double)maxInt,maxInt+1);
			HistogramFunction histCum = new HistogramFunction((double)minInt,(double)maxInt,maxInt+1);
			for(int i=0;i<hist.size();i++) {
				hist.set(i, poissDist.probability(i));
				histCum.set(i, poissDist.cumulativeProbability(i));
			}
		
			int lowBound = (int)Math.round(histCum.getClosestXtoY(0.025));
			if(histCum.getY(lowBound)<0.025)
				lowBound += 1;
		
			int highBound = (int)Math.round(histCum.getClosestXtoY(0.975));
			if(histCum.getY(highBound)<0.975)
				highBound += 1;
		
			System.out.println(lowBound);
			System.out.println(highBound);
		
			ArrayList<HistogramFunction> funcList = new ArrayList<HistogramFunction>();
			funcList.add(hist);
			funcList.add(histCum);
		
			ArrayList<PlotCurveCharacterstics> plotCharList = new ArrayList<PlotCurveCharacterstics>();
			plotCharList.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLACK));
			plotCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
			GraphWindow graph = new GraphWindow(funcList, "Poisson Distributions",plotCharList); 
			graph.setX_AxisLabel("Num");
			graph.setY_AxisLabel("Probability");
		
//			System.out.println(getExpectedNumEvents(0d, 1d, 1d, 1d, 1.025, 0.05, 0d, 1e6));
		
//			double val = Math.log(Double.MAX_VALUE);
//			System.out.println(val);
//			System.out.println(Math.exp(700.0));
		
//			testJeanneCalc();
//			testAndyCalc();

			return;
		}




		// Subcommand : Test #2
		// Command format:
		//  test2
		// Test the implementation of getExpectedNumEvents using pow_diff_div
		// by comparing to the old way.

		if (args[0].equalsIgnoreCase ("test2")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("AftershockStatsCalc : Invalid 'test2' subcommand");
				return;
			}

			// Construct a list of p values

			ArrayList<Double> pList = new ArrayList<Double>();

			pList.add(0.9);
			pList.add(0.99);
			pList.add(0.999);
			pList.add(0.9999);
			pList.add(0.99999);
			pList.add(0.999999);
			pList.add(0.9999999);
			pList.add(1.0);
			pList.add(1.0000001);
			pList.add(1.000001);
			pList.add(1.00001);
			pList.add(1.0001);
			pList.add(1.001);
			pList.add(1.01);
			pList.add(1.1);

			// Other parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double magMain = 7.5;
			double magMin = 2.5;
			double tMinDays = 10.0;
			double tMaxDays = 11.0;

			// Calculate results for each p

			for (Double pd : pList) {
				double p = pd;

				// New value

				double newVal = getExpectedNumEvents(a, b, magMain, magMin, p, c, tMinDays, tMaxDays);

				// Old value

				double oldVal;
				double k = convertProductivityTo_k(a, b, magMain, magMin);

				if(p!=1) {
					double oneMinusP= 1-p;
					oldVal = (k/oneMinusP)*(Math.pow(c+tMaxDays,oneMinusP) - Math.pow(c+tMinDays,oneMinusP));
				}
				else {
					oldVal = k*(Math.log(c+tMaxDays) - Math.log(c+tMinDays));
				}

				// Display results

				System.out.println (String.format ("%.7f  %.15e  %.15e", p, newVal, oldVal));
			}

			return;
		}




		// Subcommand : Test #3
		// Command format:
		//  test3
		// Test the implementation of simAftershockSequence by producing a simulated
		// aftershock with Jeanne's parameters for 0 to 30 days.
		// The resulting list of (magnitude, time in days) is written to standard output.
		// There should be about 3000 aftershocks.

		if (args[0].equalsIgnoreCase ("test3")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("AftershockStatsCalc : Invalid 'test3' subcommand");
				return;
			}

			// Parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
//			double capG = 4.5;
			double capH = 0.75;
			double tMinDays = 0.0;
			double tMaxDays = 30.0;

			// The completeness parameters don't match Jeanne's simulation data, here is my best guess
			// (note this value yields a time of completeness of exactly 1 day)

			double capG = 1.25;

			// Run the simulation

			ObsEqkRupList aftershock_list = simAftershockSequence(a, b, magMain, magCat, capG, capH, p, c, tMinDays, tMaxDays);

			// Output the number of aftershocks
					
			System.out.println (String.format ("Generated %d aftershocks", aftershock_list.size()));

			// Output the results, but stop at 10,000

			for (int i = 0; i < aftershock_list.size(); ++i) {
				if (i == 10000) {
					System.out.println (String.format ("... plus %d more aftershocks", aftershock_list.size() - i));
					break;
				}
				ObsEqkRupture rupture = aftershock_list.get(i);
				double mag = rupture.getMag();
				double tDays = ((double)(rupture.getOriginTime())) / ((double)MILLISEC_PER_DAY);
				System.out.println (String.format ("%.2f  %.8f", mag, tDays));
			}

			return;
		}




		// Subcommand : Test #4
		// Command format:
		//  test4
		// Test the numeric integration code.

		if (args[0].equalsIgnoreCase ("test4")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("AftershockStatsCalc : Invalid 'test4' subcommand");
				return;
			}

			// Parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
			double capG = 1.25;
			double capH = 0.75;
			double tMinDays = 0.0;
			double tMaxDays = 30.0;

			// R&J, direct calculation

			double rj_direct = getExpectedNumEvents(a, b, magMain, magCat, p, c, tMinDays, tMaxDays);

			// R&J, analytic

			double rj_analytic = getPageExpectedNumEvents(a, b, magMain, magCat, 10.0, 0.0, p, c, tMinDays, tMaxDays);

			// R&J, numeric

			double rj_numeric = getPageExpectedNumEvents(a, b, magMain, 1.00, 1.25, 0.0, p, c, tMinDays, tMaxDays);

			// R&J, page completeness

			double rj_page = getPageExpectedNumEvents(a, b, magMain, magCat, capG, capH, p, c, tMinDays, tMaxDays);

			// Output the results
					
			System.out.println (String.format ("rj_direct   = %.15e", rj_direct));
			System.out.println (String.format ("rj_analytic = %.15e", rj_analytic));
			System.out.println (String.format ("rj_numeric  = %.15e", rj_numeric));
			System.out.println (String.format ("rj_page     = %.15e", rj_page));

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("AftershockStatsCalc : Unrecognized subcommand : " + args[0]);
		return;

	}

}
