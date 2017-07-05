/**
 * 
 */
package scratch.aftershockStatisticsETAS;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;



/**
 * @author field
 * @author van der Elst
 *
 */
public class ETAS_StatsCalc {
	
	public final static double MILLISEC_PER_YEAR = 1000*60*60*24*365.25;
	public final static long MILLISEC_PER_DAY = 1000*60*60*24;

	
	/**
	 * This computes the log-likelihood for the given modified Omori parameters according to 
	 * equation (6) of Ogata (1983; J. Phys. Earth, 31,115-124).
	 * @param k
	 * @param p
	 * @param c
	 * @param tMinDays - the start time of of the catalog in days from main shock
	 * @param tMaxDays - the end time of of the catalog in days from main shock
	 * @param relativeEventTimes - in order of occurrence relative to main shock
	 * @return
	 */
	public static double getLogLikelihoodForOmoriParams(double k, double p, double c, double tMinDays, double tMaxDays, double[] relativeEventTimes) {
		double funcA=Double.NaN;
		int n=relativeEventTimes.length;
		if(p == 1)
			funcA = Math.log(tMaxDays+c) - Math.log(tMinDays+c);
		else
			funcA = (Math.pow(tMaxDays+c,1-p) - Math.pow(tMinDays+c,1-p)) / (1-p);
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
	 * @param p
	 * @param c
	 * @param tMinDays - the start time of of the catalog in days from main shock
	 * @param tMaxDays - the end time of of the catalog in days from main shock
	 * @param numEvents - the number of events in the catalog
	 * @return
	 */
	public static double getMaxLikelihood_k(double p, double c, double tMinDays, double tMaxDays, int numEvents) {
		double funcA=Double.NaN;
		if(p == 1)
			funcA = Math.log(tMaxDays+c) - Math.log(tMinDays+c);
		else
			funcA = (Math.pow(tMaxDays+c,1-p) - Math.pow(tMinDays+c,1-p)) / (1-p);
		
// System.out.println("getMaxLikelihood_k: \t"+p+"\t"+c+"\t"+tMinDays+"\t"+tMaxDays+"\t"+funcA+"\t"+numEvents+"\t"+(numEvents/funcA));
		return (double)numEvents/funcA;
	}
	
	/**
	 * This computes the log-likelihood for the given ETAS parameters according to 
	 * Ogata (1992). This method includes an extra k-parameter in addition to the a-parameter.
	 * 		k = 10^{a + alpha(Mmain - Mref};
	 * This allows for the mainshock productivity to be fit independent of the productivity 
	 * of the aftershocks, which may be desired, given alleatory variability in mainshock
	 * productivity. Implementation is a bit awkward in that productivity is computed outside the
	 * likelihood calculator. Not the way I would do it.
	 * @param k
	 * @param a
	 * @param p
	 * @param c
	 * @param alpha
	 * @param Mref - reference magnitude (usually the minimum, or completeness, magnitude)
	 * @param tMinDays - the start time of of the catalog in days from main shock
	 * @param tMaxDays - the end time of of the catalog in days from main shock
	 * @param relativeEventTimes - in order of occurrence relative to main shock
	 * @param aftershockMagnitudes - magnitudes of mainshock and observed aftershocks
	 * @return
	 * 
	 */

	
	/**
	 * This converts the productivity value from "a" to "k"
	 * @param a
	 * @param b
	 * @param magMain
	 * @param magMin
	 * @return
	 */
	public static double convertProductivityTo_k(double a, double b, double magMain, double magMin) {
		return Math.pow(10.0, a+b*(magMain-magMin));
	}
	
	/**
	 * This converts the productivity value from "k" to "a"
	 * @param k
	 * @param b
	 * @param magMain
	 * @param magMin
	 * @return
	 */
	public static double convertProductivityTo_a(double k, double b, double magMain, double magMin) {
		return Math.log10(k) - b*(magMain-magMin);
	}

	/**
	 * This returns the Reasenberg Jones (1989, 1994) expected number of primary aftershocks 
	 * between time tMin and tMax (days) for 
	 * the given arguments.
	 * @param a
	 * @param b
	 * @param magMain
	 * @param magMin
	 * @param p 
	 * @param c - days
	 * @param tMin - beginning of forecast time window (since origin time), in days
	 * @param tMax - end of forecast time window (since origin time), in days
	 * @return
	 */
	public static double getExpectedNumEvents(double a, double b, double magMain, double magMin, double p, double c, double tMinDays, double tMaxDays) {
		double k = convertProductivityTo_k(a, b, magMain, magMin);
		if(p!=1) {
			double oneMinusP= 1-p;
			return (k/oneMinusP)*(Math.pow(c+tMaxDays,oneMinusP) - Math.pow(c+tMinDays,oneMinusP));
		}
		else {
			return k*(Math.log(c+tMaxDays) - Math.log(c+tMinDays));
		}
	}
	
	
	
	/**
	 * This returns the expected number of primary aftershocks as a function of time
	 * 
	 * @param a
	 * @param b
	 * @param magMain
	 * @param magMin
	 * @param p
	 * @param c
	 * @param tMin
	 * @param tMax
	 * @param tDelta
	 * @return
	 */
	public static  EvenlyDiscretizedFunc getExpectedNumWithTimeFunc(double a, double b, double magMain, double magMin, double p, double c,  double tMin, double tMax, double tDelta) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(tMin+tDelta/2, tMax-tDelta/2, (int)Math.round((tMax-tMin)/tDelta));
		for(int i=0;i<func.size();i++) {
			double binTmin = func.getX(i) - tDelta/2;
			double binTmax = func.getX(i) + tDelta/2;
			double yVal = getExpectedNumEvents(a, b, magMain, magMin, p, c, binTmin, binTmax);
			func.set(i,yVal);
		}
		func.setName("Expected Number of Primary Aftershocks for "+tDelta+"-day intervals");
		func.setInfo("for a="+a+", b="+b+", p="+p+", c="+c+", magMain="+magMain+", magMin="+magMin);
		return func;
	}
	
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @param magMain
	 * @param magMin
	 * @param p
	 * @param c
	 * @param tMin
	 * @param tMax
	 * @param tDelta
	 * @return
	 */
	public static  EvenlyDiscretizedFunc getExpectedCumulativeNumWithTimeFunc(double a, double b, double magMain, double magMin, double p, double c,  double tMin, double tMax, double tDelta) {
		EvenlyDiscretizedFunc cumFunc = new EvenlyDiscretizedFunc(tMin+tDelta/2, tMax-tDelta/2, (int)Math.round((tMax-tMin)/tDelta));
		EvenlyDiscretizedFunc numWithTimeFunc = getExpectedNumWithTimeFunc(a, b, magMain, magMin, p, c,  tMin, tMax, tDelta);
		double sum = 0;
		for(int i=0;i<cumFunc.size();i++) {
			sum += numWithTimeFunc.getY(i);
			cumFunc.set(i,sum);
		}
		return cumFunc;
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
	public static double getMaxLikelihood_b_value(ObsEqkRupList rups, double magComplete,
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
	 */
	public static double[] getDaysSinceMainShockArray(ObsEqkRupture mainShock, ObsEqkRupList aftershockList) {
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
	
	
	public static double[] readAndysFile() {
		
		try {
			BufferedReader buffRead = new BufferedReader(new InputStreamReader(
					ETAS_StatsCalc.class.getResourceAsStream("AndysSimulationData.txt")));
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
					ETAS_StatsCalc.class.getResourceAsStream("JeannesSimulationData.txt")));
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

	
	
	public static double[] linspace(double min, double max, int npts){
		double[] vec = new double[npts];
		double dx = (max-min)/(npts-1);
		for(int i = 0; i < npts-1; i++){
			vec[i] = min + i*dx;
		}
		vec[npts-1] = max;
		return vec;
	}
	
	public static double[] logspace(double min, double max, int npts){
		double[] vec = new double[npts];
		double log_min = Math.log10(min);
		double log_max = Math.log10(max);
		double dx = (log_max-log_min)/npts;
		for(int i = 0; i < npts-1; i++){
			vec[i] = Math.pow(10, log_min + i*dx);
			
		}
		vec[npts-1] = max;
		
		return vec;
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
	
// MC(t) NOT IMPLEMENTED
//	   /**
//	    *  adaptive quadrature integration.  This is put here because it is not clear how general/robust this code is
//	    *  (e.g., with respect to singularities in the function).
//	    *  
//	    *  TODO In fact, this should be tested over the viable range of rate-decay functions that could be given
//	    * @param func
//	    * @param startTime
//	    * @param endTime
//	    * @return
//	    */
//    public static double adaptiveQuadratureIntegration(ETAS_AftershockModel_SequenceSpecific func, double startTime, double endTime) {
//    	final double EPSILON = 1e-6;
//    	double a=startTime;
//    	double b=endTime;
//        double h = b - a;
//        double c = (a + b) / 2.0;
//        double d = (a + c) / 2.0;
//        double e = (b + c) / 2.0;
//        double Q1 = h/6  * (func.getRateAboveMagCompleteAtTime(a) + 4*func.getRateAboveMagCompleteAtTime(c) + func.getRateAboveMagCompleteAtTime(b));
//        double Q2 = h/12 * (func.getRateAboveMagCompleteAtTime(a) + 4*func.getRateAboveMagCompleteAtTime(d) + 2*func.getRateAboveMagCompleteAtTime(c) + 4*func.getRateAboveMagCompleteAtTime(e) + func.getRateAboveMagCompleteAtTime(b));
//        if (Math.abs(Q2 - Q1) <= EPSILON)
//            return Q2 + (Q2 - Q1) / 15;
//        else
//            return adaptiveQuadratureIntegration(func,a, c) + adaptiveQuadratureIntegration(func,c, b);
//    }

    public static Location getCentroid(ObsEqkRupture mainshock, ObsEqkRupList aftershocks) {
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
		while (lon > 180)
			lon -= 360;
		Location centroid = new Location(lat, lon);
		double dist = LocationUtils.horzDistanceFast(mainshock.getHypocenterLocation(), centroid);
		System.out.println("Centroid: "+(float)lat+", "+(float)lon+" ("+(float)dist+" km from epicenter)");
		return centroid;
	}

    /** returns radius in km for magnitude and stressDrop in MPa 
     * 
     */
    public static double magnitude2radius(double magnitude, double stressDrop){
    	double r = Math.pow( 7.0/16.0, 1.0/3.0 ) * Math.pow(10.0, 0.5*magnitude - 1.0/3.0*Math.log10(stressDrop) - 2.0);
    	
		return r; 
	}
    
    /** returns magnitude for radius in km and stressDrop in MPa
     * 
     */
    public static double radius2magnitude(double radius, double stressDrop){
		//return (Math.log10(Math.pow(radius*1000,3) * stressDrop * 1e6 * 16/7) - 9)/1.5;
		double mag = 2.0 * (Math.log10(radius * Math.pow(16.0/7.0, 1.0/3.0)) + 1.0/3.0*Math.log10(stressDrop) + 2.0) ;
		return mag ;
	}
    
    /** fitMainshockLineSource fits a mainshock line source to the aftershocks in aftershockFitList, consisting of npts points
     * 
     */
    public static ObsEqkRupList fitMainshockLineSource(ObsEqkRupture mainshock, ObsEqkRupList aftershockFitList, double stressDrop){
    	double q = 0.68; //quantile of the AS zone to try to capture
        
    	//mainshock coords (lat0 and lon0 not actaully used...)
    	double mag0, lat0, lon0;
        mag0 = mainshock.getMag();
        lat0 = mainshock.getHypocenterLocation().getLatitude();
        lon0 = mainshock.getHypocenterLocation().getLongitude();
        
        // aftershock coords
        double[] mag = new double[aftershockFitList.size()];
        double[] lat = new double[mag.length];
        double[] lon = new double[mag.length];
        double[] w = new double[mag.length];
        for(int i = 0; i < mag.length; i++){
        	ObsEqkRupture as = aftershockFitList.get(i);
        	mag[i] = as.getMag();
        	lat[i] = as.getHypocenterLocation().getLatitude();
        	lon[i] = as.getHypocenterLocation().getLongitude();
        }

        // compute weighted centroid location
        double lonwsum = 0;
        double latwsum = 0;
        double wsum = 0;
    	for(int i = 0; i < mag.length; i++){
    		w[i] = Math.pow(10, 0.5*mag[i]);	//weight by lengthscale
    		lonwsum += lon[i]*w[i];
    		latwsum += lat[i]*w[i];
    		wsum += w[i];
    	}
        double lonc = lonwsum / wsum;
        double latc = latwsum / wsum;

        System.out.println("Centroid location: " + lonc + " " + latc);
        
        // compute x, y locs relative to centroid
        double[] dx = new double[mag.length];
        double[] dy = new double[mag.length];
        double[] r = new double[mag.length];
        double[] th = new double[mag.length];
        for(int i = 0; i < mag.length; i++){
        	dy[i] = (lat[i] - latc) * 111.111;
        	dx[i] = (lon[i] - lonc) * Math.cos(Math.toRadians(latc)) * 111.111;
        	th[i] = Math.atan2(dy[i], dx[i]);
        	r[i] = Math.sqrt(dx[i]*dx[i] + dy[i]*dy[i]);
        }
        	
        // compute an initial guess for orientation based on L1 norm solution
        double omega0, sinsum = 0, cossum = 0;
        for(int i = 0; i< mag.length; i++){
        	 sinsum += Math.sin(2*dy[i]/dx[i]);
        	 cossum += Math.cos(2*dy[i]/dx[i]);
        }
        omega0 = 0.5*Math.atan2(sinsum, cossum);	// initial guess
        
        // find the best-fitting line with a grid search
        int searchDepth = 3; //number of times to refine grid
        double deltaOm = Math.PI/2;	//search range
        double dOm = deltaOm/10;	//search increment
        double[] SSE = new double[21];
        
        double[] phi = new double[th.length];
        double[] p = new double[th.length];
        double[] d = new double[th.length];
    	
        System.out.println("Seeking the optimal fit line");
        for (int n = 0; n < searchDepth; n++){
        	System.out.println(n + " out of " + searchDepth + " grid refinements...");
        	double SSEmin = Double.POSITIVE_INFINITY;
        	int iBest = 1;
        	for(int i = 0; i < SSE.length; i++){
            	double om = omega0 - deltaOm + dOm*i;
                
                // compute Sum of Squared distances from test line
                for(int j = 0; j < th.length; j++){
            		phi[j] = th[j] - om;
            	  	p[j] = r[j] * Math.sin(phi[j]);
//            	  	d[j] = r[j] * Math.cos(phi[j]);
            	  	SSE[i] += p[j]*p[j];
            	}
                
            	if(SSE[i] < SSEmin){
            		SSEmin = SSE[i];
            		iBest = i;
            	}
        	}
        	omega0 = omega0 - deltaOm + dOm*iBest;
        	deltaOm = dOm;
        	dOm = dOm/10;
        }
         
        System.out.println("Best fit mainshock strike: " + omega0*180/Math.PI + " north of east");
        
        // compute coordinates relative to best-fit line
        for(int i = 0; i < th.length; i++){
    		phi[i] = th[i] - omega0;
    	  	p[i] = r[i] * Math.sin(phi[i]);
    	  	d[i] = r[i] * Math.cos(phi[i]);
    	}

        //coordinate transformation matrix
        double[][] coeff = new double[2][2];
        coeff[0][0] = Math.cos(omega0);
        coeff[0][1] = -Math.sin(omega0);
        coeff[1][0] = Math.sin(omega0);
        coeff[1][1] = Math.cos(omega0);

        //find quantiles of x and y.
        java.util.Arrays.sort(d);
        double xmin = d[(int) Math.floor( (0.5 - q/2)*d.length )];
        double xmax = d[(int) Math.floor( (0.5 + q/2)*d.length )];
        java.util.Arrays.sort(p);
        double ymin = p[(int) Math.floor( (0.5 - q/2)*p.length )];
        double ymax = p[(int) Math.floor( (0.5 + q/2)*p.length )];
       
        System.out.println("Quantiles " + xmin + " " + xmax + " " + ymin + " " + ymax);
        
        //compute aspect ratio of aftershocks
        double Dx = xmax - xmin;
        double Dy = ymax - ymin;
        double c = Dx/Dy;
        double a = Math.sqrt(c) * magnitude2radius(mag0, stressDrop);
        double b = a/c;
        
        System.out.println("Dimensions of ellipse: " + a + " " + " " + b + "; EQ mag: " + mag0 + "; DS:" + stressDrop + "; EQ radius: " + magnitude2radius(mag0, stressDrop) );
        
        //build set of equivalent sources to represent line source
        double rMSequivalent, magEquivalent;
        int nEquivalent;
        
        if (a < b)
        	rMSequivalent = a;
        else
        	rMSequivalent = b;
        
        // get the ideal number of equivalent sources based on the true magEquivalent
        magEquivalent = radius2magnitude(rMSequivalent, stressDrop);
    	if (magEquivalent > mag0)
    		System.out.println("somethings wrong with the line-fitting code");
    	nEquivalent = (int) Math.round(Math.pow(10, 1.0*(mag0 - magEquivalent)));
    	
    	System.out.println("Replacing Mainshocks with " + nEquivalent + " M=" + magEquivalent + " sources");
        // now find the equivalent magnitude for an integer number of sources
    	if(nEquivalent < 1)
    		nEquivalent = 1;
    	magEquivalent = mag0 - Math.log10(nEquivalent);
    
    	// set up equivalent sources
    	double dEquivalent, pEquivalent, dxEquivalent, dyEquivalent, latEquivalent, lonEquivalent;
    	ObsEqkRupList mainshockEquivalent = new ObsEqkRupList();
    	
    	for(int n = 0; n < nEquivalent; n++){
    		if (nEquivalent == 1)
    			dEquivalent = 0;
    		else
    			dEquivalent = -b + 2*b/(nEquivalent - 1) * n;

    		// rotate p and d into dx and dy, convert to lat lon
    		dxEquivalent = dEquivalent * Math.cos(omega0);
    		dyEquivalent = dEquivalent * Math.sin(omega0);

    		latEquivalent = latc + dyEquivalent/111.111;
    		lonEquivalent = lonc + dxEquivalent/111.111 / Math.cos(Math.toRadians(latc));

    		Location loc = new Location(latEquivalent,lonEquivalent);
    		ObsEqkRupture newSource = new ObsEqkRupture(); 
    		newSource.setMag(magEquivalent);
    		newSource.setHypocenterLocation(loc);
    		newSource.setOriginTime(mainshock.getOriginTime());

    		// add Equivalent source to list
    		mainshockEquivalent.add(newSource);
    	}
    	
    	return mainshockEquivalent; 
    }
    
    
    
    		
    
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
//		double mean = 0.78;
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

		
		
//		System.out.println(getExpectedNumEvents(0d, 1d, 1d, 1d, 1.025, 0.05, 0d, 1e6));
		
//		double val = Math.log(Double.MAX_VALUE);
//		System.out.println(val);
//		System.out.println(Math.exp(700.0));
		
//		testJeanneCalc();
//		testAndyCalc();
	}

}
