package scratch.aftershockStatisticsETAS;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.sql.Time;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;


/**
 * This represents an ETAS (Ogata 1992, 2006) aftershock model.
 * 
 * This is a modification of the RJ_AftershockModel class that implements the forecast through an ETAS model 
 * Subclasses should normalize Likelihood values so they sum to 1.0 over the range of parameter-values specified;
 * this can be accomplished using the convertArrayToLikelihood(double maxLogLikeVal) method.
 * 
 * @author field
 * @author van der Elst
 *
 */
public abstract class ETAS_AftershockModel {
	
	Boolean D=false;	// debug flag
	
	double base_a = Double.NaN;  
	double mean_a, sigma_a;
	double mean_p, sigma_p;
	double mean_c, sigma_logc;
	
	double dataStartTimeDays, dataEndTimeDays;
	double testTotalLikelihood;
	double forecastMinDays;
	double forecastMaxDays;
	
	ObsEqkRupList aftershockList;
	ObsEqkRupture mainShock;
	
	double b, magMain, magComplete; 
	double[] magAftershocks;
	double[] relativeTimeAftershocks;
	double alpha, refMag, mu;
	double min_a, max_a, delta_a=0, min_p, max_p, delta_p=0, min_c, max_c, delta_c=0;	//grid search not used for forecast
	int num_a = 101, num_p = 101, num_c = 101;
	double[][][] likelihood;
	double[][][] amsMatrix;
	int max_a_index=-1;
	int max_p_index=-1;
	int max_c_index=-1;
	
	double[] a_vec, p_vec, c_vec;
	
	
//	double[] a_vec;
//	double[] p_vec;
//	double[] c_vec;
//	double[][][] likelihood;
//	
	ArbDiscrEmpiricalDistFunc numMag5_DistributionFunc = null;
	ArbDiscrEmpiricalDistFunc num_DistributionFunc = null;
	double tMinDaysCurrent=-1, tMaxDaysCurrent=-1;

	double maxMag;
	int maxGenerations;
	int nSims;
	ETAScatalog simulatedCatalog;	//results of the stochastic simulations
	
	
	
	/**
	 * This converts the likelihood from log-likelihood to likelihood values, making sure NaNs do not occur 
	 * due to high logLikelihoods, and re-normalizes the likelihood so as values sum to 1.0.
	 * 
	 * @param maxLogLikeVal - the maximum values in the input likelihood
	 * @return
	 */
	protected double convertLogLikelihoodArrayToLikelihood(double maxLogLikeVal) {
		double total=0;
		double corr = 0;
		if(maxLogLikeVal>100.0)	// values above ~700 cause NaNs from Math.exp(), so we subtract some number from all values to avoid such numerical problems
			corr = maxLogLikeVal-100;
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			for(int pIndex=0;pIndex<num_p;pIndex++) {
				for(int cIndex=0;cIndex<num_c;cIndex++) {
					double like = Math.exp(likelihood[aIndex][pIndex][cIndex]-corr);
					likelihood[aIndex][pIndex][cIndex] = like;
					total += likelihood[aIndex][pIndex][cIndex];
				}
			}
		}

		// now re-normalize all likelihood values so they sum to 1.0
		double testTotal=0;
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			for(int pIndex=0;pIndex<num_p;pIndex++) {
				for(int cIndex=0;cIndex<num_c;cIndex++) {
					likelihood[aIndex][pIndex][cIndex] /= total;
					testTotal += likelihood[aIndex][pIndex][cIndex];
				}
			}
		}
		return testTotal;
	}
	
	
	/**
	 * This converts the likelihood from log-likelihood to likelihood values, making sure NaNs do not occur 
	 * due to high logLikelihoods, but does NOT re-normalize the likelihood so as values sum to 1.0.
	 * 
	 * @param maxLogLikeVal - the maximum values in the input likelihood
	 * @return
	 */
	protected double convertLogLikelihoodArrayToLikelihood_nonNormalized(double[][][] likelihood, double maxLogLikeVal) {
		double total=0;
		double corr = 0;
		if(maxLogLikeVal>100.0)	// values above ~700 cause NaNs from Math.exp(), so we subtract some number from all values to avoid such numerical problems
			corr = maxLogLikeVal-100;
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			for(int pIndex=0;pIndex<num_p;pIndex++) {
				for(int cIndex=0;cIndex<num_c;cIndex++) {
					double like = Math.exp(likelihood[aIndex][pIndex][cIndex]-corr);
					likelihood[aIndex][pIndex][cIndex] = like;
					total += likelihood[aIndex][pIndex][cIndex];
				}
			}
		}

		return total;
	}
	
	public double getMaxLikelihood_a() { return get_a(max_a_index);}
	
	public double getMaxLikelihood_p() { return get_p(max_p_index);}
	
	public double getMaxLikelihood_c() { return get_c(max_c_index);}
		
	protected double get_a(int aIndex) { return a_vec[aIndex];}
	
	protected double get_p(int pIndex) { return p_vec[pIndex];}
	
	protected double get_c(int cIndex) { return c_vec[cIndex];}
	
	public double[][][] getLikelihood() { return likelihood; };
	
	public void setMagComplete(double magComplete) {
		double prevMc = this.magComplete;
		this.magComplete = magComplete;
		this.mean_a += Math.log10( (this.maxMag - prevMc)/(this.maxMag - magComplete) );
	}
	
	/**
	 * This computes the distribution of the number of M≥5.0 events given all a, p, and c values, as well as the associated
	 * weight for each set of values. 
	 * @param tMinDays
	 * @param tMaxDays
	 */
	public ArbDiscrEmpiricalDistFunc computeNumMag5_DistributionFunc(double tMinDays, double tMaxDays) {
		//special case with forecastMag = 5.0; to be removed.
		if(tMinDaysCurrent == tMinDays && tMaxDaysCurrent == tMaxDays && numMag5_DistributionFunc != null) // already computed
			return numMag5_DistributionFunc;
		
		tMinDaysCurrent = tMinDays;
		tMaxDaysCurrent = tMaxDays;
		numMag5_DistributionFunc = new ArbDiscrEmpiricalDistFunc();
		
		numMag5_DistributionFunc = computeNum_DistributionFunc(tMinDays, tMaxDays, 5.0);
		
		return numMag5_DistributionFunc;
	}
	
	
	/**
	 * This computes the distribution of the number of M≥forecastMag events for the suite of simulated ETAS catalogs
	 * @param tMinDays
	 * @param tMaxDays
	 * @param forecastMag
	 */
	public ArbDiscrEmpiricalDistFunc computeNum_DistributionFunc(double tMinDays, double tMaxDays, double forecastMag) {
		
		num_DistributionFunc = new ArbDiscrEmpiricalDistFunc();
		
		int[] numM = new int[simulatedCatalog.nSims];
		
		List<double[]> eqCat = new ArrayList<double[]>();
		
		Point2D pt = new Point2D.Double();
		
		//cycle through the simulated catalogs
		for(int i = 0; i < simulatedCatalog.nSims; i++){
			eqCat = simulatedCatalog.getETAScatalog(i); 	//double[] eqCat = {relativeTime, magnitude, generationNumber}
			numM[i] = 0;
			//count all events in time window and magnitude range in this catalog
			for(double[] eq : eqCat){
				if(eq[0] >= tMinDays && eq[0] < tMaxDays && eq[1] >= forecastMag)
					numM[i] ++;
			}
			pt.setLocation(numM[i], 1d/simulatedCatalog.nSims);
			num_DistributionFunc.set(pt);	//increment the distribution
		}
		
		
		if(D) {
			System.out.println("N≥5 mean = "+num_DistributionFunc.getMean());
			System.out.println("N≥5 mode = "+num_DistributionFunc.getApparentMode());
			System.out.println("N≥5 median = "+num_DistributionFunc.getDiscreteFractile(0.5)); //discrete median
			System.out.println("N≥5 2.5 Percentile = "+num_DistributionFunc.getDiscreteFractile(0.025));
			System.out.println("N≥5 97.5 Percentile = "+num_DistributionFunc.getDiscreteFractile(0.975));
		}
		return num_DistributionFunc;
	}
	
	
	/**
	 * This gives the number of aftershocks associated with the maximum likelihood a/p/c parameters (which represents the mode
	 * in the number of events space) above the given minimum magnitude and over the specified site span.  A GR distribution with
	 * no upper bound is assumed.
	 * @param magMin
	 * @param tMinDays
	 * @param tMaxDays
	 * @return
	 */
	public double getModalNumEvents(double magMin, double tMinDays, double tMaxDays) {
//		num_DistributionFunc = computeNum_DistributionFunc( tMinDays,  tMaxDays,  magMin);
		return getFractileNumEvents(magMin, tMinDays, tMaxDays, 0.5);
	}
	
	public double getFractileNumEvents(double magMin, double tMinDays, double tMaxDays, double fractile) {
		num_DistributionFunc = computeNum_DistributionFunc( tMinDays,  tMaxDays,  magMin);
		return num_DistributionFunc.getDiscreteFractile(fractile);
	}
	
	
	
	/**
	 * This gives the number of aftershocks as a function of time for the maximum likelihood a/p/c parameters 
	 * (which represents the mode in the number of events space) above the given minimum magnitude and over 
	 * the specified time span.  A GR distribution with no upper bound is assumed.
	 * @param magMin
	 * @param tMinDays - left edge of first time interval
	 * @param tMaxDays - right edge of last time interval
	 * @param tDelta
	 * @return
	 */
	public EvenlyDiscretizedFunc getModalCumNumEventsWithTime(double magMin, double tMinDays, double tMaxDays, double tDelta) {
		
		
				
//		EvenlyDiscretizedFunc cumFunc = new EvenlyDiscretizedFunc(tMinDays+tDelta/2, tMaxDays-tDelta/2, (int)Math.round((tMaxDays-tMinDays)/tDelta));
		EvenlyDiscretizedFunc cumFunc = new EvenlyDiscretizedFunc(tMinDays, tMaxDays, (int)Math.round((tMaxDays-tMinDays)/tDelta )+1);
		double count = 0, obsCount = 0;
		
		int i = 0;
		for(double t = tMinDays; t <= tMaxDays; t += tDelta){
			if(t <= this.dataEndTimeDays){
				//compute expected number for comparison with known quakes (plot fit)
				obsCount = getExpectedNumEvents(magMin, tMinDays, t);
				
				cumFunc.set(i++, obsCount);
			}
			else{
				//report simulated number of events
				count = getModalNumEvents( magMin,  tMinDays,  t);
				cumFunc.set(i++, count + obsCount);
			}
			
			
		}
		return cumFunc;
	}
		
public ArbitrarilyDiscretizedFunc getModalCumNumEventsWithLogTime(double magMin, double tMinDays, double tMaxDays, int numPts) {
	
	return getFractileCumNumEventsWithLogTime(magMin,tMinDays,tMaxDays,numPts,0.5);
//		ArbitrarilyDiscretizedFunc cumFunc = new ArbitrarilyDiscretizedFunc();
//		
//		double count = 0, obsCount = 0, expCount = 0;
//		
//		double[] tvec = ETAS_StatsCalc.logspace(tMinDays, tMaxDays, numPts);
//		for(double t : tvec){
//			if(t <= this.dataEndTimeDays){
//				//compute expected number for comparison with known quakes (plot fit)
//				expCount = getExpectedNumEvents(magMin, tMinDays, t);
//				
//				cumFunc.set(t, expCount);
//			}
//			
//			else{
//				//get number of events observed so far
//				ObsEqkRupList subList = aftershockList.getRupsAboveMag(magMin);
//				obsCount= subList.size();
//				System.out.println(obsCount);
//				subList = subList.getRupsBefore((long)(mainShock.getOriginTime() + tMaxDays*ETAS_StatsCalc.MILLISEC_PER_DAY));
//				obsCount = subList.size();
//				
//				System.out.println(obsCount);
//				
//				//report simulated number of events
//				count = getModalNumEvents( magMin,  tMinDays,  t);
//				cumFunc.set(t, count + obsCount);
//			}
//			
//			
//		}
//		return cumFunc;
	}
		

public ArbitrarilyDiscretizedFunc getFractileCumNumEventsWithLogTime(double magMin, double tMinDays, double tMaxDays, int numPts, double fractile) {
	
	ArbitrarilyDiscretizedFunc cumFunc = new ArbitrarilyDiscretizedFunc();
	
	double count = 0, baseCount = 0, obsCount = 0, expCount = 0;
	
	double[] tvec = ETAS_StatsCalc.logspace(tMinDays, tMaxDays, numPts);
	
	// get number of events observed prior to forecastWindow
	ObsEqkRupList subListMc = aftershockList.getRupsAboveMag(magMin);
	ObsEqkRupList subList = subListMc.getRupsBefore((long)(mainShock.getOriginTime() + tMinDays*ETAS_StatsCalc.MILLISEC_PER_DAY));
	baseCount = subList.size();
	
	
	subList = subListMc.getRupsBefore((long)(mainShock.getOriginTime() + forecastMinDays*ETAS_StatsCalc.MILLISEC_PER_DAY));
	obsCount = subList.size();
	
	for(double t : tvec){
		if(t <= forecastMinDays){
			//compute expected number for comparison with known quakes (plot fit)
			expCount = getExpectedNumEvents(magMin, tMinDays, t);
			
			
			cumFunc.set(t, baseCount + expCount);
		}
		else{
		
			
			//report simulated number of events (plot forecast)
			count = getFractileNumEvents( magMin,  tMinDays,  t, fractile);
			cumFunc.set(t, obsCount + count);
		}
		
		
	}
	return cumFunc;
}
	
	public double getExpectedNumEvents(double magMin, double tMinDays, double tMaxDays){

		double[] relativeTimeAftershocks = new double[this.relativeTimeAftershocks.length];
		double[] magAftershocks = new double[this.magAftershocks.length];
		
		int Nas = 0;
		for(int i = 0; i < magAftershocks.length; i++){
			if(this.magAftershocks[i] >= magMin){
				magAftershocks[Nas] = this.magAftershocks[i];
				relativeTimeAftershocks[Nas] = this.relativeTimeAftershocks[i];
				Nas++;
			}
		}
		relativeTimeAftershocks = Arrays.copyOf(relativeTimeAftershocks, Nas);
		magAftershocks = Arrays.copyOf(magAftershocks, Nas);
		
		
		
		// now compute numbers
		double a = getMaxLikelihood_a();
		double p = getMaxLikelihood_p();
		double c = getMaxLikelihood_c();
		double k = Math.pow(10, a + alpha*(mainShock.getMag() - refMag) + b*(refMag - magComplete));
		
		System.out.println("a, k_ms, p, c, alpha, refMag, magComplete, magMain, Nas:");
		System.out.println(a +" "+ k +" "+ p +" "+ c + " "+ alpha +" "+ refMag +" "+ magComplete +" "+ magMain +" "+ Nas);
		double timeIntegral = Double.NaN;
//		int Nas=relativeTimeAftershocks.length;	//the number of aftershocks, not counting the mainshock
				
		//compute total number at end of window due to mainshock
		double Ntot;
		if (p == 1)
			timeIntegral = Math.log(tMaxDays + c) - Math.log(tMinDays + c);
		else
			timeIntegral = (Math.pow(tMaxDays + c, 1-p) - Math.pow(tMinDays + c, 1-p)) / (1-p);
		Ntot = k*timeIntegral;		//mainshock contribution

		double[] productivity = new double[Nas];
		
		for(int i=0; i<Nas; i++){
			//compute productivity for this aftershock
			productivity[i] = Math.pow(10, a + alpha*(magAftershocks[i] - refMag) + b*(refMag - magComplete));	//productivity of this aftershock
			
			//compute number at end of window due to this aftershock
			if(relativeTimeAftershocks[i] <= tMaxDays){
				if(p == 1)
					timeIntegral = Math.log(tMaxDays - relativeTimeAftershocks[i] + c) - Math.log(c); 
				else
					timeIntegral = (Math.pow(tMaxDays - relativeTimeAftershocks[i] + c, 1-p) - Math.pow(c, 1-p)) / (1-p);

				Ntot += productivity[i]*timeIntegral;	//aftershock Contributions
			}
			//System.out.format(" %d", (int) Ntot); 
			
			
		}
		
		return Ntot;
	}
	
	/**
	 * This returns the cumulative MFD associated with the maximum likelihood a/p/c parameters (which represents the mode
	 * in the number of events space) and over the specified site span.   A GR distribution with no upper bound is assumed.
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays
	 * @param tMaxDays
	 * @return
	 */
	public EvenlyDiscretizedFunc getModalCumNumMFD(double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc mfd = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
		
		//get the number of forecast events with mag >= Mc 
		System.out.println("magComplete = " + magComplete);	//debug
		double baseNum = getModalNumEvents(mfd.getInterpolatedY(magComplete), tMinDays, tMaxDays);
		
		for(int i=0;i<mfd.size();i++) {
			double scaledNum = baseNum*Math.pow(10, -b*(mfd.getX(i) - magComplete));	 
			mfd.set(i, scaledNum);
//			mfd.set(i, getModalNumEvents(mfd.getX(i), tMinDays, tMaxDays));
		}
		mfd.setName("Modal Num Events");
		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
//		double totExpNum = getExpectedNumEvents(minMag, tMinDays, tMaxDays);
//		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(b, totExpNum, minMag, maxMag, numMag);
//		mfd.setName("Expected Num Incr. MFD");
//		mfd.setInfo("Total Expected Num = "+totExpNum);
		return mfd;
	}
	
	
	
	/**
	 * This returns the mean cumulative MFD given each a/p/c parameter set, and the associated likelihood, 
	 * for the specified site span.   A GR distribution with no upper bound is assumed.
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays
	 * @param tMaxDays
	 * @return
	 */
	public EvenlyDiscretizedFunc getMeanCumNumMFD(double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc mfd = getModalCumNumMFD(minMag, maxMag, numMag, tMinDays, tMaxDays);
		double m5val = mfd.getInterpolatedY(5.0);
		computeNumMag5_DistributionFunc(tMinDays, tMaxDays);
		mfd.scale(numMag5_DistributionFunc.getMean()/m5val);
		mfd.setName("Mean Num Events");
		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		return mfd;
	}
	
	/**
	 * This returns the fractile MFD, associated with the fractile number forecast. The algorithm is approximate.
	 * We start with the event count for a given fractile, considering only the total count above the reference magnitude 
	 * (i.e. the minimum magnitude used in the simulation). Since magnitudes above this value are distributed GR, we don't
	 * have to use the simulation (with its poor sampling at large magnitudes) -- we can use the GR distribution directly.
	 * Using the fractile event count, scaled by GR statistics, gives the range due to epistemic variability alone.  We also 
	 * want to capture aleatory variability -- the probability of getting N events above magnitude M, given N0 events above 
	 * magnitude M0. This is treated as a Bernoulli trial. For each magnitude bin, we find the number N that meets the given
	 * fractile, using a normal approximation as an initial guess for that number N. This is only approximate, because we should
	 * be computing the fractile N given the entire distribution of N0, and not just the fractile value of N0.  
	 *
	 * @param fractile - the fractile (percentile/100) for the distribution
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays
	 * @param tMaxDays
	 * @return
	 */
	public EvenlyDiscretizedFunc getCumNumMFD_Fractile(double fractile, double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc mfd = getModalCumNumMFD(minMag, maxMag, numMag, tMinDays, tMaxDays);
		
		computeNumMag5_DistributionFunc(tMinDays, tMaxDays);
		
		
//		ArbDiscrEmpiricalDistFunc numSpan = numMag5_DistributionFunc.deepClone(); 
		
				
		
//		double baseNum = numMag5_DistributionFunc.getInterpolatedFractile(fractile);
		
//		double baseNum = getModalNumEvents(mfd.getInterpolatedY(refMag), tMinDays, tMaxDays);
		double Mref=0, probLarger=0, m5=0, m5_scaled=0;
//		double mMax, mMax_scaled;

	
		
		ArbitrarilyDiscretizedFunc num_CumulativeDistFunc = num_DistributionFunc.getCumDist();
		
		if(num_CumulativeDistFunc.getY(1) > fractile)
			m5 = 0;
		else
			m5 = num_CumulativeDistFunc.getFirstInterpolatedX(fractile);
		
		m5_scaled = m5 * Math.pow(10, -b*(minMag - magComplete));
		
		for(int i = 0; i < mfd.size(); i++){
			Mref = mfd.getX(i);
			mfd.set(i, 0d);
			probLarger = Math.pow(10, -b*( Mref - minMag ));

			BinomialDistribution binoDist = new BinomialDistribution((int) Math.ceil(m5_scaled), probLarger);
			double cumProb;
			int nTest;
			if(fractile < 0.5)
				nTest = (int) Math.round(m5_scaled*probLarger - 2*Math.sqrt(m5_scaled*probLarger*(1-probLarger)));
			else
				nTest = (int) Math.round(m5_scaled*probLarger + 2*Math.sqrt(m5_scaled*probLarger*(1-probLarger)));			
			if(nTest<0) nTest=0;
			
			cumProb = binoDist.cumulativeProbability(nTest);
			
//			System.out.println(i +" "+ Mref +" "+ probLarger +" "+ m5 +" "+ m5_scaled +" "+ (nTest) +" "+ cumProb);
			
			if(cumProb > fractile){
				while(cumProb > fractile){
					cumProb = binoDist.cumulativeProbability(--nTest);
//					System.out.println(i +" "+ Mref +" "+ probLarger +" "+ m5 +" "+ m5_scaled*probLarger +" "+ nTest +" "+ cumProb);
				}
				mfd.set(i, nTest+1);
			} else {
				while(cumProb < fractile){
					cumProb =  binoDist.cumulativeProbability(++nTest);;
//					System.out.println(i +" "+ Mref +" "+ probLarger +" "+ m5 +" "+ m5_scaled*probLarger +" "+ nTest +" "+ cumProb);
				}
				mfd.set(i, nTest-1);
			}
			
			
			;
		}
		
			
	
		mfd.setName(fractile+" Fractile for Num Events");
		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		return mfd;
		
//		double m5val;
//		if (minMag > 5 || maxMag < 5)
//			m5val = getModalCumNumMFD(5d, 5d, 1, tMinDays, tMaxDays).getY(0);
//		else
//			m5val = mfd.getInterpolatedY(5.0);
//		computeNumMag5_DistributionFunc(tMinDays, tMaxDays);
//		mfd.scale(numMag5_DistributionFunc.getInterpolatedFractile(fractile)/m5val);
//		mfd.setName(fractile+" Fractile for Num Events");
//		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
//		return mfd;
		
	}

	
	
	
	/**
	 * This returns the fractile MFDs implied by each a/p/c parameter set and their associated likelihoods, 
	 * for the specified site span.   A GR distribution with no upper bound is assumed.  Both epistemic uncertainty
	 * and aleatory variability are considered.
	 * @param fractileArray - the fractile (percentile/100) for the distribution
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays
	 * @param tMaxDays
	 * @return EvenlyDiscretizedFunc[]
	 */
	public EvenlyDiscretizedFunc[] getCumNumMFD_FractileWithAleatoryVariability(double[] fractileArray, double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
			
		
		EvenlyDiscretizedFunc[] mfdArray = new EvenlyDiscretizedFunc[fractileArray.length];
		
		
		
		for(int i=0;i<fractileArray.length;i++) {
//			mfdArray[i] =  getCumNumMFD_Fractile(fractileArray[i], minMag, maxMag, numMag, tMinDays, tMaxDays);
			mfdArray[i] = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
			mfdArray[i].setName(fractileArray[i]+" Fractile for Num Events, including aleatory variability");
			mfdArray[i].setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		}
		
		
		for(int i=0;i<numMag;i++) {
			double mag = mfdArray[0].getX(i);	// any MFD will do, as they all have the same x-axis values
			double[] valsArray = getCumNumFractileWithAleatory(fractileArray, mag, tMinDays, tMaxDays);
			for(int j=0;j<fractileArray.length;j++) {
				mfdArray[j].set(i,valsArray[j]);
//				System.out.println("\tworking on "+mag);
			}

		}
		return mfdArray;
	}

	
	
	public double[] getCumNumFractileWithAleatory(double[] fractileArray, double mag, double tMinDays, double tMaxDays) {
		computeNum_DistributionFunc(tMinDays, tMaxDays, mag);
//		
		double[] fractValArray = new double[fractileArray.length];
		
		for(int i = 0; i < fractileArray.length; i++){
			if(mag < magComplete){
				double fractValComplete = getSubSampleInterpolatedFractile(num_DistributionFunc, fractileArray[i]);
				
				if(fractValComplete <= 4 && Math.abs(fractileArray[i] - 0.5) < 1e-9){
					double probOne = 1 - num_DistributionFunc.getY(0);
					fractValComplete = -Math.log(1-probOne);
					System.out.println(fractValComplete);
				}
				fractValArray[i] = (int) Math.round(Math.pow(10, -b*(mag - magComplete)) * fractValComplete);	
				
			}
			else
				fractValArray[i] = num_DistributionFunc.getDiscreteFractile(fractileArray[i]);
			
			
		}
		
		return fractValArray;
	}

	public double getCumNumMeanWithAleatory(double mag, double tMinDays, double tMaxDays) {
		computeNum_DistributionFunc(tMinDays, tMaxDays, mag);
//		
		double mean;
		
			if(mag < magComplete){
				mean =  num_DistributionFunc.getMean() * Math.pow(10, -b*(mag - magComplete));	
			}
			else
				mean =  num_DistributionFunc.getMean();
			
		return mean;
	}

	
	public double getSubSampleInterpolatedFractile(ArbDiscrEmpiricalDistFunc num_DistributionFunc, double fractile){
		double xval = 0;
		double yval = 0, ycum = 0;
		

		int i = 0;
		while(ycum < fractile){
			xval = num_DistributionFunc.getX(i);
			yval = num_DistributionFunc.getY(i);
			ycum += yval;
			i++;
		}
		if(i == 1 && Math.abs(ycum - 1) < 1e-6) //it's all zeros
			return 0;
		else{
			double xpre = xval;
			double xpost = xval+1;
			double ypre = ycum-yval;
			double ypost = ycum;

			return(xpre + (xpost-xpre)* (fractile - ypre)/(ypost - ypre));
		}

		
	}
	

	
	/**
	 * This returns the PDF of a, which is a marginal distribution if either c or p 
	 * are unconstrained (either num_p or num_c not equal to 1). Null is returned if
	 * a is constrained (num_a=1). (how important is this? I want to show pdf with one value)
	 * @return
	 */
	public HistogramFunction getPDF_a() {
//		if(num_a == 1) {
		if(num_a == 0) {
			// changed this to still return a 1-bar histogram.
			return null;
		}
		else {
			System.out.println(min_a +" "+ num_a +" "+ delta_a);
			HistogramFunction hist = new HistogramFunction(min_a, max_a, num_a);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
//						hist.add(get_a(aIndex), likelihood[aIndex][pIndex][cIndex]);
//						hist.add(a_vec[aIndex], likelihood[aIndex][pIndex][cIndex]);
						hist.add(aIndex,  likelihood[aIndex][pIndex][cIndex]);
						
					}
				}
			}
			String name = "PDF of a-value";
			if(num_a == 1)
				name += " (constrained)"; 
			else if(num_p !=1 || num_c != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) {
//				System.out.println("PDF of a-value:  "+hist);
				System.out.println("PDF of a-value: totalTest = "+hist.calcSumOfY_Vals());
			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}
	
	/**
	 * This returns the PDF of p, which is a marginal distribution if either a or c 
	 * are unconstrained (either num_a or num_c not equal to 1). Null is returned if
	 * p is constrained (num_p=1).
	 * @return
	 */
	public HistogramFunction getPDF_p() {
		if(num_p == 0) {
			// changed this to still return a 1-bar histogram
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_p, max_p, num_p);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
//						hist.add(get_p(pIndex), likelihood[aIndex][pIndex][cIndex]);
						hist.add(pIndex, likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
			String name = "PDF of p-value";
			if(num_p == 1)
				name += " (constrained)";
			else if(num_a !=1 || num_c != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) {
				System.out.println("PDF of p-value: totalTest = "+hist.calcSumOfY_Vals() + " Max = " + hist.getMaxY());
			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}
	
	/**
	 * This returns the PDF of c, which is a marginal distribution if either a or p 
	 * are unconstrained (either num_a or num_p not equal to 1). Null is returned if
	 * c is constrained (num_c=1).
	 * @return
	 */
	public HistogramFunction getPDF_c() {
		if(num_c == 0) {
			// still return a 1-bar histogram
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_c, max_c, num_c);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
//						hist.add(get_c(cIndex), likelihood[aIndex][pIndex][cIndex]);
						hist.add(cIndex, likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
			String name = "PDF of c-value";
			if(num_c == 1)
				name += " (constrained)";
			else if(num_a !=1 || num_p != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) {
				System.out.println("PDF of c-value: totalTest = "+hist.calcSumOfY_Vals());
			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}

	
	/**
	 * This returns the PDF of logc, which is a marginal distribution if either a or p 
	 * are unconstrained (either num_a or num_p not equal to 1). Null is returned if
	 * c is constrained (num_c=1).
	 * @return
	 */
	public HistogramFunction getPDF_logc() {
		if(num_c == 0) {
			
			return null;
		}
		else {
			double min_logc = Math.log10(min_c);
			double max_logc = Math.log10(max_c);
//			double delta_logc = (max_logc - min_logc)/num_c;

			HistogramFunction hist = new HistogramFunction(min_logc, max_logc, num_c);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
//						hist.add(get_c(cIndex), likelihood[aIndex][pIndex][cIndex]);
						hist.add(cIndex, likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
			String name = "PDF of logc-value";
			if(num_c == 1)
				name += " (constrained)";
			else if(num_a !=1 || num_p != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) {
				System.out.println("PDF of logc-value: totalTest = "+hist.calcSumOfY_Vals());
			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}

	/**
	 * This returns a 2D PDF for a and p, which is a marginal distribution if c 
	 * is unconstrained (num_c not equal to 1). Null is returned if either
	 * a or p are constrained (num_a=1 or num_p=1).
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_a_and_p() {
		if(num_a == 1 || num_p == 1) {
			return null;
		}
		else {
			double delta_a = (max_a - min_a)/((double)num_a - 1);
			double delta_p = (max_p - min_p)/((double)num_p - 1);
			
			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_a, num_p, min_a, min_p, delta_a, delta_p);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						double prevVal = hist2D.get(aIndex,pIndex);
						hist2D.set(aIndex,pIndex, prevVal+likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
//			String name = "2D PDF of a vs p";
//			if(num_c != 1)
//				name += " (marginal)";
			if(D) {
				System.out.println("2D PDF of a vs p: totalTest = "+hist2D.getSumZ());
			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}

	
	
//	/**
//	 * This returns a 2D PDF for a and c, which is a marginal distribution if p 
//	 * is unconstrained (num_p not equal to 1). Null is returned if either
//	 * a or c are constrained (num_a=1 or num_c=1).
//	 * @return
//	 */
//	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_a_and_c() {
//		if(num_a == 1 || num_c == 1) {
//			return null;
//		}
//		else {
//			double delta_c = (max_c - min_c)/((double)num_c - 1);
//			double delta_a = (max_a - min_a)/((double)num_a - 1);
//			
//			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_a, num_c, min_a, min_c, delta_a, delta_c);
//			for(int aIndex=0;aIndex<num_a;aIndex++) {
//				for(int pIndex=0;pIndex<num_p;pIndex++) {
//					for(int cIndex=0;cIndex<num_c;cIndex++) {
//						double prevVal = hist2D.get(aIndex,cIndex);
//						hist2D.set(aIndex,cIndex, prevVal+likelihood[aIndex][pIndex][cIndex]);
//					}
//				}
//			}
////			String name = "2D PDF of a vs c";
////			if(num_p != 1)
////				name += " (marginal)";
//			if(D) {
//				System.out.println("2D PDF of a vs c: totalTest = "+hist2D.getSumZ());
//			}
//			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
//			return hist2D;
//		}
//	}

	/**
	 * This returns a 2D PDF for a and c, which is a marginal distribution if p 
	 * is unconstrained (num_p not equal to 1). Null is returned if either
	 * a or c are constrained (num_a=1 or num_c=1).
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_a_and_logc() {
		if(num_a == 1 || num_c == 1) {
			return null;
		}
		else {
			double min_logc = Math.log10(min_c);
			double max_logc = Math.log10(max_c);
			double delta_logc = (max_logc - min_logc)/num_c;
			double delta_a = (max_a - min_a)/((double)num_a - 1);
			
			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_a, num_c, min_a, min_logc, delta_a, delta_logc);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						double prevVal = hist2D.get(aIndex,cIndex);
						hist2D.set(aIndex,cIndex, prevVal+likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
//			String name = "2D PDF of a vs c";
//			if(num_p != 1)
//				name += " (marginal)";
			if(D) {
				System.out.println("2D PDF of a vs logc: totalTest = "+hist2D.getSumZ());
			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}

	/**
	 * This returns a 2D PDF for c and p, which is a marginal distribution if a 
	 * is unconstrained (num_a not equal to 1). Null is returned if either
	 * c or p are constrained (num_c=1 or num_p=1).
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_logc_and_p() {
		if(num_c == 1 || num_p == 1) {
			return null;
		}
		else {
			double min_logc = Math.log10(min_c);
			double max_logc = Math.log10(max_c);
			double delta_logc = (max_logc - min_logc)/((double)num_c - 1);
			double delta_p = (max_p - min_p)/((double)num_p - 1);
	
			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_c, num_p, min_logc, min_p, delta_logc, delta_p);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						double prevVal = hist2D.get(cIndex,pIndex);
						hist2D.set(cIndex,pIndex, prevVal+likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
//			String name = "2D PDF of c vs p";
//			if(num_a != 1)
//				name += " (marginal)";
			if(D) {
				System.out.println("2D PDF of logc vs p: totalTest = "+hist2D.getSumZ());
			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}
	
	/**
	 * This sets the likelihood and maximum likelihood values from the a-values in the given discretized function,
	 *  and holding the other parameters fixed at the values given.  The likelihood is normalized so values sum
	 *  to 1.0.
	 * @param aValueFunc
	 * @param b
	 * @param p
	 * @param c
	 */
	protected void setArrayAndMaxLikelyValuesFrom_aValueFunc(EvenlyDiscretizedFunc aValueFunc, double b, double p, double c) {
		EvenlyDiscretizedFunc pValueFunc = new EvenlyDiscretizedFunc(p, p, 1);
		EvenlyDiscretizedFunc cValueFunc = new EvenlyDiscretizedFunc(c, c, 1);
		
		setArrayAndMaxLikelyValuesFrom_aValueFunc(aValueFunc, b, pValueFunc, cValueFunc);
	}
	
	protected void setArrayAndMaxLikelyValuesFrom_aValueFunc(EvenlyDiscretizedFunc aValueFunc, EvenlyDiscretizedFunc pValueFunc, EvenlyDiscretizedFunc cValueFunc) {
		setArrayAndMaxLikelyValuesFrom_aValueFunc(aValueFunc, b, pValueFunc, cValueFunc);
	}
	
	protected void setArrayAndMaxLikelyValuesFrom_aValueFunc(EvenlyDiscretizedFunc aValueFunc, double b, EvenlyDiscretizedFunc pValueFunc, EvenlyDiscretizedFunc cValueFunc) {
			
		this.delta_a = aValueFunc.getDelta();
		this.min_a = aValueFunc.getMinX();
		this.max_a = aValueFunc.getMaxX();
		this.num_a = aValueFunc.size();

		this.num_p=1;
		this.min_p = mean_p;
		this.max_p = mean_p;
		this.max_p_index=0;		
		
		this.num_c=1;
		this.min_c = mean_c;
		this.max_c = mean_c;
		this.max_c_index=0;
		
		likelihood = new double[num_a][num_p][num_c];
		double maxWt= Double.NEGATIVE_INFINITY;
		double totWt=0;
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			double wt = aValueFunc.getY(aIndex);
			likelihood[aIndex][0][0] = wt;
			totWt+=wt;
			if(wt>maxWt) {
				max_a_index=aIndex;
				maxWt=wt;
			}
		}
		// now normalize so that it sums to 1.0
		double totTest=0;
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			likelihood[aIndex][0][0] /= totWt;
			totTest += likelihood[aIndex][0][0];
		}

		if(D) {
			System.out.println("a-values range:\t"+min_a+"\t"+max_a+"\t"+num_a+"\t"+delta_a);
			System.out.println("totTest="+(float)totTest);
			System.out.println("getMaxLikelihood_a()="+getMaxLikelihood_a());
			System.out.println("getMaxLikelihood_p()="+getMaxLikelihood_p());
			System.out.println("getMaxLikelihood_c()="+getMaxLikelihood_c());
		}
		
	}
	
	public double getMainShockMag() {return magMain;}
	
	public double[] getAftershockMags() {return magAftershocks;}
	
	/**
	 * Returns an array of aftershock magnitudes from an ObsEqkRupList. Useful for computing ETAS forecasts.
	 * 
	 * @param rupList
	 */
	public double[] getAftershockMags(ObsEqkRupList rupList){
		
//		System.out.println(rupList);
		int Nas = rupList.size();
		
		double[] magnitudes = new double[Nas];
		ObsEqkRupture rup;
		
		for(int i=0 ; i<Nas ; i++){
			rup = rupList.get(i);
			magnitudes[i] = rup.getMag();
		}
		return magnitudes;
	}
	
	public double get_b() {return b;}


	// getters/setters commented out until needed:
//	public double getMin_a() { return min_a;}
//	public double getMin_p() { return min_p;}
//	public double getMin_c() { return min_c;}
//	public double getMax_a() { return max_a;}
//	public double getMax_p() { return max_p;}
//	public double getMax_c() { return max_c;}
//	public double getDelta_a() { return delta_a;}
//	public double getDelta_p() { return delta_p;}
//	public double getDelta_c() { return delta_c;}
//	public int getNum_a() { return num_a;}
//	public int getNum_p() { return num_p;}
//	public int getNum_c() { return num_c;}

	/**
	 * Returns a 2d grid of earthquake rates based on epicenters. 
	 * 
	 */
	public GriddedGeoDataSet getRateModel2D(double spacing, double stressDrop, double mainshockFitDuration){
		
		// TODO: assign seismogenicDepth as variable somewhere
		double seismogenicDepth = 10;	//in km
		
		//load ETAS solution
		double a = this.getMaxLikelihood_a();
		double p = this.getMaxLikelihood_p();
		double c = this.getMaxLikelihood_c();
		double mref = this.magComplete;
			
		ObsEqkRupList aftershockPlotList = new ObsEqkRupList();
		ObsEqkRupList aftershockFitList = new ObsEqkRupList();
		ObsEqkRupList equivalentMainshockSources = new ObsEqkRupList();
		
		aftershockPlotList = this.aftershockList.getRupsBefore((long) (this.mainShock.getOriginTime() + this.forecastMinDays* 24*60*60*1000));
		aftershockFitList = this.aftershockList.getRupsBefore((long) (this.mainShock.getOriginTime() + mainshockFitDuration*24*60*60*1000));
		
		// set up grid centered on mainshock
		Location centerLocation = ETAS_StatsCalc.getCentroid(this.mainShock, aftershockPlotList);
		double lat0 = centerLocation.getLatitude();
		double lon0 = centerLocation.getLongitude();
		double geomFactor = Math.cos(Math.toRadians(lat0));
		double mainshockRadius = ETAS_StatsCalc.magnitude2radius(this.mainShock.getMag(), stressDrop);	//in km

		// make grid extend to 10 times source radius or at least 1 degree 
		double latmin, latmax, lonmin, lonmax;
		if (mainshockRadius*10 < 111.111){
			latmin = lat0 - 1;
			latmax = lat0 + 1;
			lonmin = lon0 - 1/geomFactor;
			lonmax = lon0 + 1/geomFactor;
		}else{
			latmin = lat0 - mainshockRadius*10/111.111;
			latmax = lat0 + mainshockRadius*10/111.111;
			lonmin = lon0 - mainshockRadius*10/geomFactor/111.111;
			lonmax = lon0 + mainshockRadius*10/geomFactor/111.111;
		}
		GriddedRegion griddedRegion = new GriddedRegion(new Location(latmin, lonmin),
				new Location(latmax, lonmax), spacing, null);
		GriddedGeoDataSet gridData = new GriddedGeoDataSet(griddedRegion, false);
		
		System.out.println(latmin + " " + latmax + " " + lonmin + " " + lonmax + " " + griddedRegion.getNodeCount());
		
		// fit finite mainshock source to early aftershocks
		System.out.println("Fitting " + aftershockFitList.size() + " early aftershocks, out of " + aftershockList.size() + " total aftershocks.");
		if (aftershockFitList.size() >= 3){ 
			equivalentMainshockSources = ETAS_StatsCalc.fitMainshockLineSource(mainShock, aftershockFitList, stressDrop);
		} else {
			equivalentMainshockSources.add(mainShock);
		}
		for (ObsEqkRupture rup:equivalentMainshockSources)
			System.out.println(rup.getMag() + " " + rup.getHypocenterLocation().getLatitude() + " " + rup.getHypocenterLocation().getLongitude());
		
		// compute rates at each point in the rate map for the mainshock equivalent sources
		double prevVal, x, y, newVal, x0, y0, mag0, t0;
		System.out.println("computing MS rate integral from day " + forecastMinDays + " to " + forecastMaxDays);
		for (ObsEqkRupture rup : equivalentMainshockSources){
			System.out.println(rup);
			x0 = rup.getHypocenterLocation().getLongitude();
			y0 = rup.getHypocenterLocation().getLatitude();
			mag0 = rup.getMag();
			t0 = (rup.getOriginTime() - mainShock.getOriginTime()) / ETAS_StatsCalc.MILLISEC_PER_DAY;
			
			for (int i=0; i<gridData.size(); i++) {
				Location gridLoc = gridData.getLocation(i);
				x = gridLoc.getLongitude();
				y = gridLoc.getLatitude();
				
				newVal = rateXY(x,y,t0,mag0,x0,y0, stressDrop, forecastMinDays, forecastMaxDays, seismogenicDepth);
				
				gridData.set(i, gridData.get(i) + newVal);
			}
		}
		
		// compute rates at each point in the rate map for the aftershock sources
		System.out.println("computing AS rate integral for time " + forecastMinDays + " to " + forecastMaxDays);
		for (ObsEqkRupture rup : aftershockList){
			System.out.println(rup);
			x0 = rup.getHypocenterLocation().getLongitude();
			y0 = rup.getHypocenterLocation().getLatitude();
			mag0 = rup.getMag();
			t0 = (rup.getOriginTime() - mainShock.getOriginTime()) / ETAS_StatsCalc.MILLISEC_PER_DAY;
			
			for (int i=0; i<gridData.size(); i++) {
				Location gridLoc = gridData.getLocation(i);
				x = gridLoc.getLongitude();
				y = gridLoc.getLatitude();
				
				newVal = rateXY(x,y,t0,mag0,x0,y0, stressDrop, forecastMinDays, forecastMaxDays, seismogenicDepth);
				
				gridData.set(i, gridData.get(i) + newVal);
			}
		}

		return gridData;
	}
	
	/**	compute rate at one point for one source
     * 
     */
    private double rateXY(double lon, double lat, double t0, double mag0, double lon0, double lat0, double stressDrop, double ts, double te, double H){
    	// quick convert to x,y
    	double dx = (lon-lon0)*Math.cos(Math.toRadians(lat0))*111.111;
    	double dy = (lat-lat0)*111.111;
    	
    	// constants
    	double r = Math.sqrt(dx*dx + dy*dy);
    	double d = ETAS_StatsCalc.magnitude2radius(mag0, stressDrop);
    	
    	// compute productivity for weighting this event
    	double productivity = Math.pow(10d, this.getMaxLikelihood_a() + 1d*(mag0 - this.refMag));
    	// compute rate, integrated over seismogenic depth H
    	double spatialDecay = H / (d*d + r*r) / Math.pow(H*H/4d + r*r + d*d, 1d/2d) * d/(2*Math.PI);
    	
    	double p = this.getMaxLikelihood_p();
    	double c = this.getMaxLikelihood_c();
    	double timeIntegral = 1d/(1d - p) * ( Math.pow(te - t0 + c, 1d-p) - Math.pow(ts - t0 + c, 1d-p) );  
//        
        double rate = productivity *  timeIntegral * spatialDecay;
        
        return rate;
    }
		
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}


	
}
