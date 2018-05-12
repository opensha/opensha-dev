package scratch.aftershockStatisticsETAS;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariateOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.DoubleArray;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.gui.infoTools.CalcProgressBar;

/**
 * This represents an ETAS (Ogata 1992, 2006) aftershock model.
 * 
 * This class implements the forecast through an ETAS model 
 * Subclasses should (sometimes) normalize Likelihood values so they sum to 1.0 over the range of parameter-values specified;
 * this can be accomplished using the convertArrayToLikelihood(double maxLogLikeVal) method.
 * 
 * @author field
 * @author van der Elst
 *
 */
public abstract class ETAS_AftershockModel {

	private Boolean D=false;	// debug flag

	protected ArbDiscrEmpiricalDistFunc num_DistributionFunc = null;
	protected CalcProgressBar progress;
	
	// generic/prior model parameters
	protected TectonicRegime regime;
	protected double mean_ams, sigma_ams;
	protected double mean_a, sigma_a;
	protected double mean_p, sigma_p;
	protected double mean_c, sigma_logc;
	protected double b, bSigma;
	protected double alpha, refMag, mu;
	protected double ac; // this is the productivity parameter used for time-dependentMc calculations
	
	// forecast parameters
	protected double dataStartTimeDays, dataEndTimeDays;
	protected double forecastMinDays, forecastMaxDays;

	// sequence data
	protected ObsEqkRupList aftershockList;
	protected ObsEqkRupture mainShock;
	protected double magMain, magComplete; 
	protected double[] magAftershocks;
	protected double[] relativeTimeAftershocks;
	
	// likelihood search parameters
	protected double min_ams, max_ams, delta_ams=0, min_a, max_a, delta_a=0, min_p, max_p, delta_p=0, min_c, max_c, delta_c=0;	//grid search not used for forecast
	protected int num_ams = 101, num_a = 1, num_p = 1, num_c = 1;
	protected double[] ams_vec, a_vec, p_vec, c_vec;
	protected double[][][][] likelihood;
	protected double[][][][] epiLikelihood;
	protected int max_ams_index = -1;
	protected int max_a_index=-1;
	protected int max_p_index=-1;
	protected int max_c_index=-1;

	// simulation parameters
	protected double maxMag;
	protected int maxGenerations;
	protected int nSims;
	protected ETAScatalog simulatedCatalog;	//results of the stochastic simulations
	protected Boolean timeDependentMc = false;

	/**
	 * This converts the likelihood from log-likelihood to likelihood values, making sure NaNs do not occur 
	 * due to high logLikelihoods, and re-normalizes the likelihood so as values sum to 1.0.
	 * 
	 * @param maxLogLikeVal - the maximum values in the input likelihood
	 * @return
	 */
	protected double convertLogLikelihoodArrayToLikelihood(double maxLogLikeVal) {
		double total=0;
		//		double corr = 0;
		//		if(maxLogLikeVal>100.0)	// values above ~700 cause NaNs from Math.exp(), so we subtract some number from all values to avoid such numerical problems
		//			corr = maxLogLikeVal-100;

		double corr = maxLogLikeVal;
		double loglike, like;
		for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						loglike = likelihood[amsIndex][aIndex][pIndex][cIndex] - corr;
						if(loglike < -20)
							like = 0;
						else
							like = Math.exp(loglike);

						//						if(Double.isNaN(like))
						//							System.out.println("NaN encountered in likelihood for index [" + amsIndex + " " + aIndex + " " + pIndex + " " + cIndex + "]" );

						likelihood[amsIndex][aIndex][pIndex][cIndex] = like;
						total += likelihood[amsIndex][aIndex][pIndex][cIndex];
					}
				}
			}
		}

		// now re-normalize all likelihood values so they sum to 1.0
		double testTotal=0;
		for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						likelihood[amsIndex][aIndex][pIndex][cIndex] /= total;
						testTotal += likelihood[amsIndex][aIndex][pIndex][cIndex];
					}
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
	protected double convertLogLikelihoodArrayToLikelihood_nonNormalized(double[][][][] likelihood, double maxLogLikeVal) {
		double total = 0;
		//		double corr;
		double like, loglike;
		//		if(maxLogLikeVal>100.0)	// values above ~700 cause NaNs from Math.exp(), so we subtract some number from all values to avoid such numerical problems
		//			corr = maxLogLikeVal-100;
		double corr = maxLogLikeVal;
		for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						loglike = likelihood[amsIndex][aIndex][pIndex][cIndex]-corr;
						if(loglike < -20)
							like = 0;
						else
							like = Math.exp(loglike);

						likelihood[amsIndex][aIndex][pIndex][cIndex] = like;
						total += likelihood[amsIndex][aIndex][pIndex][cIndex];
					}
				}
			}
		}

		return total;
	}

	public double getMaxLikelihood_ams() { return get_ams(max_ams_index);}

	public double getMaxLikelihood_a() { return get_a(max_a_index);}

	public double getMaxLikelihood_p() { return get_p(max_p_index);}

	public double getMaxLikelihood_c() { return get_c(max_c_index);}

	protected double get_ams(int amsIndex) { return ams_vec[amsIndex];}

	protected double get_a(int aIndex) { return a_vec[aIndex];}

	protected double get_p(int pIndex) { return p_vec[pIndex];}

	protected double get_c(int cIndex) { return c_vec[cIndex];}

	public double getMainShockMag() {return magMain;}

	public double get_b() {return b;}

	public void set_b(double b){ this.b = b; }
	
	public double get_bSigma() {return bSigma;}
	
	public void set_bSigma(double bSigma){ this.bSigma = bSigma; }

	public void set_alpha(double alpha){ this.alpha = alpha; }

	public double[][][][] getLikelihood() { return likelihood; };

	public double getForecastMinDays(){
		return forecastMinDays;
	};

	public double getForecastMaxDays(){
		return forecastMaxDays;
	};

	public void setRegime(TectonicRegime regime){
		this.regime = regime;
	}

	public TectonicRegime getRegime(){
		return regime;
	}

	public int get_nSims(){
		return nSims;
	}

	public void setMagComplete(double magComplete) {
		double prevMc = this.magComplete;
		this.magComplete = magComplete;
		this.mean_a += Math.log10( (this.maxMag - prevMc)/(this.maxMag - magComplete) );
		this.mean_ams += Math.log10( (this.maxMag - prevMc)/(this.maxMag - magComplete) );
	}
	
	//  additional getters/setters commented out until needed:
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
	 * This computes the distribution of the number of M≥forecastMag events for the suite of simulated ETAS catalogs
	 * @param tMinDays
	 * @param tMaxDays
	 * @param forecastMag
	 */
	public ArbDiscrEmpiricalDistFunc computeNum_DistributionFunc(double tMinDays, double tMaxDays, double forecastMag) {

		num_DistributionFunc = new ArbDiscrEmpiricalDistFunc();

		int[] numM = new int[simulatedCatalog.nSims];

		List<float[]> eqCat = new ArrayList<float[]>();

		Point2D pt = new Point2D.Double();

		//cycle through the simulated catalogs
		for(int i = 0; i < simulatedCatalog.nSims; i++){
			eqCat = simulatedCatalog.getETAScatalog(i); 	//double[] eqCat = {relativeTime, magnitude, generationNumber}
			numM[i] = 0;
			//count all events in time window and magnitude range in this catalog
			for(float[] eq : eqCat){
				if(eq[0] > tMinDays && eq[0] <= tMaxDays && eq[1] >= forecastMag)
					numM[i] ++;
			}
			pt.setLocation(numM[i], 1d/simulatedCatalog.nSims);
			num_DistributionFunc.set(pt);	//increment the distribution
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
		if (num_DistributionFunc.size() > 0)
			return num_DistributionFunc.getDiscreteFractile(fractile);
		else
			return 0;
	}

	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
		if(D){
			System.out.println("Computing Forecast with " + nSims + " simulations (Data window: " + dataMinDays +" "+ dataMaxDays + ", Forecast window: "+ forecastMinDays +" "+ forecastMaxDays + ")");
			System.out.println("Model Params: "+ getMaxLikelihood_ams() +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete + " " + refMag);
		}

		ETAScatalog simulatedCatalog;
		try{
			simulatedCatalog = new ETAScatalog(ams_vec, a_vec, p_vec, c_vec, epiLikelihood, alpha, b, refMag, 
					mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
		} catch(Exception e) {
			e.printStackTrace();
			System.err.println("The Java Virtual Machine may have run out of memory.\n"
			+" Increase Mc by one unit and try calculating the forecast again.");
			simulatedCatalog = null;
		}

		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.simulatedCatalog = simulatedCatalog;
		this.nSims = nSims;
	}

	public ArbitrarilyDiscretizedFunc getModalCumNumEventsWithLogTime(double magMin, double tMinDays, double tMaxDays, int numPts) {
		return getFractileCumNumEventsWithLogTime(magMin,tMinDays,tMaxDays,numPts,0.5);
	}


	// get forecast number from simulation
	public ArbitrarilyDiscretizedFunc getFractileCumNumEventsWithLogTime(double magMin, double tMinDays, double tMaxDays, int numPts, double fractile) {

		ArbitrarilyDiscretizedFunc cumFunc = new ArbitrarilyDiscretizedFunc();

		double count = 0, obsCount = 0;

		// set up vector of times
		double[] tvec;
		if (tMinDays > 1e-3){
			tvec = ETAS_StatsCalc.logspace(tMinDays, tMaxDays, numPts);
		} else if (tMaxDays > 1e-3){
			tvec = ETAS_StatsCalc.logspace(1e-3, tMaxDays, numPts);
		} else {
			tvec = ETAS_StatsCalc.linspace(tMinDays, tMaxDays, numPts);
		}

		// get number of events observed prior to forecastWindow
		ObsEqkRupList subList = aftershockList.getRupsAboveMag(magMin)
				.getRupsBefore((long)(mainShock.getOriginTime() + forecastMinDays*ETAS_StatsCalc.MILLISEC_PER_DAY));
		obsCount = subList.size();

		for(double t : tvec){
			count = getFractileNumEvents(magMin,  tMinDays,  t, fractile);
			cumFunc.set(t, obsCount + count);
		}
		return cumFunc;
	}

	
	
	
	// get expected number (straight from the ETAS model)
	public ArbitrarilyDiscretizedFunc getExpectedNumEventsWithLogTime(double magMin, double tMinDays, double tMaxDays, int numPts) {

		ArbitrarilyDiscretizedFunc cumFunc = new ArbitrarilyDiscretizedFunc();

		double expCount = 0;

		// get number of events observed prior to forecastWindow
		ObsEqkRupList subList = aftershockList.getRupsAboveMag(magMin)
				.getRupsBefore((long)(mainShock.getOriginTime() + tMinDays*ETAS_StatsCalc.MILLISEC_PER_DAY));
		int baseCount = subList.size();
		
		// get number of aftershocks in forecastWindow
		subList = aftershockList.getRupsAboveMag(magMin)
				.getRupsBetween((long)(tMinDays*ETAS_StatsCalc.MILLISEC_PER_DAY), (long)(tMaxDays*ETAS_StatsCalc.MILLISEC_PER_DAY));
		
		double[] tvec, tas, tvecLog, tvecLin;
		
		tas = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, subList);
		int totalCount = tas.length;
		
		// use a minimum of numPts pts, but start with the times of the actual aftershocks
		if (totalCount < numPts) {
			// pad it with log and lin-spaced points
			if(D) System.out.println("padding time vector");
			int remainingCount = numPts - totalCount;
			int logCount = remainingCount/2;
			int linCount = remainingCount - logCount;
			
			if (tMinDays > 1e-3){
				tvecLog = ETAS_StatsCalc.logspace(tMinDays, tMaxDays, logCount);
			} else if (tMaxDays > 1e-3){
				tvecLog = ETAS_StatsCalc.logspace(1e-3, tMaxDays, logCount);
			} else {
				tvecLog = new double[0];
				linCount = remainingCount;
			}
			
			tvecLin = ETAS_StatsCalc.linspace(tMinDays, tMaxDays, linCount);
			
			tvec = ArrayUtils.addAll(ArrayUtils.addAll(tas, tvecLin), tvecLog);
			Arrays.sort(tvec);
		} else {
			//randomly thin the aftershock list
			if(D) System.out.println("thinning time vector");
			Collections.shuffle(Arrays.asList(tas));
			tvec = ArrayUtils.subarray(tas, 0, numPts);
		}

//		
//		if (tMinDays > 1e-3){
//			tvec = ETAS_StatsCalc.logspace(tMinDays, tMaxDays, numPts);
//		} else if (tMaxDays > 1e-3){
//			tvec = ETAS_StatsCalc.logspace(1e-3, tMaxDays, numPts);
//		} else {
//			tvec = ETAS_StatsCalc.linspace(tMinDays, tMaxDays, numPts);
//		}
		
		for(double t : tvec){
			//compute expected number for comparison with known quakes (plot fit)
			expCount = getExpectedNumEvents(magMin, tMinDays, t);
			cumFunc.set(t, baseCount + expCount);
		}
		return cumFunc;
	}

	// get expected number (straight from the ETAS model)
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
		double ams = getMaxLikelihood_ams();
		double a = getMaxLikelihood_a();
		double p = getMaxLikelihood_p();
		double c = getMaxLikelihood_c();
//		double k = Math.pow(10, ams + alpha*(mainShock.getMag() - refMag) + b*(refMag - magComplete));
		double k = Math.pow(10, ams + alpha*(mainShock.getMag() - magComplete)); //update 4/2/18
		double kc, cms;
		
		double timeIntegral;

		//compute total number at end of window due to mainshock
		double Ntot;
		if (timeDependentMc) {
			kc = Math.pow(10, ac-ams);
			cms = c*Math.pow(k*kc,1d/p);
			if (p == 1) {
				timeIntegral = Math.log(tMaxDays + cms) - Math.log(tMinDays + cms);
			} else {
				timeIntegral = (Math.pow(tMaxDays + cms, 1d-p) - Math.pow(tMinDays + cms, 1d-p)) / (1d-p);
			}
		} else {
			if (p == 1) 
				timeIntegral = Math.log(tMaxDays + c) - Math.log(tMinDays + c);
			else
				timeIntegral = (Math.pow(tMaxDays + c, 1d-p) - Math.pow(tMinDays + c, 1d-p)) / (1d-p);
		}
		Ntot = k*timeIntegral;		//mainshock contribution
		
		
		double[] productivity = new double[Nas];

		kc = Math.pow(10, ac-a);
		for(int i=0; i<Nas; i++){
			//compute productivity for this aftershock
//			productivity[i] = Math.pow(10, a + alpha*(magAftershocks[i] - refMag) + b*(refMag - magComplete));
			productivity[i] = Math.pow(10, a + alpha*(magAftershocks[i] - magComplete));	//update 4/2/18

			//compute number at end of window due to this aftershock
			if(relativeTimeAftershocks[i] <= tMaxDays){
				if (timeDependentMc) {
					cms = c*Math.pow(kc*productivity[i],1d/p);
					if(p == 1) {
						timeIntegral = Math.log(tMaxDays - relativeTimeAftershocks[i] + cms) - Math.log(cms); 
					} else {
						timeIntegral = (Math.pow(tMaxDays - relativeTimeAftershocks[i] + cms, 1d-p) - Math.pow(cms, 1d-p)) / (1d-p);
					}
					Ntot += productivity[i]*timeIntegral;	//aftershock Contributions
				} else {
					if(p == 1)
						timeIntegral = Math.log(tMaxDays - relativeTimeAftershocks[i] + c) - Math.log(c); 
					else
						timeIntegral = (Math.pow(tMaxDays - relativeTimeAftershocks[i] + c, 1d-p) - Math.pow(c, 1d-p)) / (1d-p);

					Ntot += productivity[i]*timeIntegral;	//aftershock Contributions
				}

			}
		}

		return Ntot;
	}

	public EvenlyDiscretizedFunc[] getCumNumMFD_FractileWithAleatoryVariability(double[] fractileArray, double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc[] mfdArray = new EvenlyDiscretizedFunc[fractileArray.length];

		for(int i=0;i<fractileArray.length;i++) {
			mfdArray[i] = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
			mfdArray[i].setName(fractileArray[i]+" Fractile for Num Events, including aleatory variability");
			mfdArray[i].setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		}

		// get fractile rates at magComplete
		//		double mag0 = mfdArray[0].getX(0);	// any MFD will do, as they all have the same x-axis values
		double mag0 = magComplete;
		double[] valsArray0 = getCumNumFractileWithAleatory(fractileArray, mag0, tMinDays, tMaxDays);

		double mag, val;
		for(int i=0;i<numMag;i++) {
			mag = mfdArray[0].getX(i);	// any MFD will do, as they all have the same x-axis values
			for(int j=0;j<fractileArray.length;j++) {
				val = valsArray0[j]*Math.pow(10, -b*(mag - mag0));
				mfdArray[j].set(i,val);
			}
		}
		return mfdArray;
	}

	public double[] getCumNumFractileWithAleatory(double[] fractileArray, double mag, double tMinDays, double tMaxDays) {
		computeNum_DistributionFunc(tMinDays, tMaxDays, mag);
		double[] fractValArray = new double[fractileArray.length];
		 
		for(int i = 0; i < fractileArray.length; i++){
			double fractValComplete = num_DistributionFunc.getDiscreteFractile(fractileArray[i]);
			// fractValComplete is an integer (the number of earthquakes for which f% of the simulations have fewer)
			// this can be zero, which is no good for scaling the rate to magnitudes lower than mag (or magComplete)
			// It is therefore necessary to produce some kind of fractional count. We do this by computing the
			// Poisson rate that gives the observed probability of zero events in the simulations. This is approximate, but
			// do you have a better idea?
			if(fractValComplete <= 4){					
				double probOne = 1 - num_DistributionFunc.getY(0); // (1 - probability of zero events)
				double lambda = -Math.log(1-probOne);
				fractValComplete = poissQuantile(lambda, fractileArray[i]);
			}
			// the preceding gives the number of magComplete used in the simulation. We want the number of mag. Scale:
			if(mag < magComplete)
				fractValArray[i] = Math.pow(10, -b*(mag - magComplete)) * fractValComplete;	
			else
				fractValArray[i] = fractValComplete;
		}
		return fractValArray;
	}

	/** Computes the "Poissonian" quantile using a continuous gamma distribution to get non-integer quantiles.
	 *  This is approximate. Ideas welcomed.
	 *  
	 *  @author Nicholas van der Elst
	 **/
	private double poissQuantile(double lambda, double fractile){
		GammaIncInverse gammaIncInverse = new GammaIncInverse();
		gammaIncInverse.setParams(lambda, fractile);
		UnivariateOptimizer fminbnd = new BrentOptimizer(1e-6,1e-9);

		UnivariatePointValuePair optimum = fminbnd.optimize(
				new UnivariateObjectiveFunction(gammaIncInverse), new MaxEval(100), GoalType.MINIMIZE,
				new SearchInterval(0,20));
		return optimum.getPoint();
	}

	private class GammaIncInverse implements UnivariateFunction{
		private double lambda, fractile;

		public double value(double x) {
			return Math.abs(fractile - Gamma.regularizedGammaQ(x, lambda));
		}

		public void setParams(double lambda, double fractile){
			this.lambda = lambda;
			this.fractile = fractile;
		}

	}

	/** This one gives probability as a function of magnitude for the specified time range
	 * 
	 **/
	public EvenlyDiscretizedFunc getMagnitudePDFwithAleatoryVariability(double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc magnitudePDF = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
		magnitudePDF.setName("probabiltiy of exceeding magnitude M");
		magnitudePDF.setInfo("Cumulative distribution (greater than or equal to each magnitude)");

		for(int i=0;i<numMag;i++) {
			double mag = magnitudePDF.getX(i);	// any MFD will do, as they all have the same x-axis values
			double value = getProbabilityWithAleatory(mag, tMinDays, tMaxDays);
			magnitudePDF.set(i,value);
		}
		return magnitudePDF;
	}

	public double getProbabilityWithAleatory(double mag, double tMinDays, double tMaxDays) {
		computeNum_DistributionFunc(tMinDays, tMaxDays, mag);

		double probOne = 1 - num_DistributionFunc.getY(0);

		// the above probability is the fraction of simulations with events above max(magComplete, mag),
		// so if mag<magComplete, we need to scale up the probability. We do this with a Poisson rate assumption.
		if(mag < magComplete){
//			double r1 = -Math.log(1d - probOne);
//			double r2 = r1*Math.pow(10, -b*(mag - magComplete));
//			probOne = 1d - Math.exp(-r2);
			probOne = 1 - Math.pow(1-probOne, Math.pow(10, -b*(mag-magComplete)));
		}

		return probOne;
	}

	/**
	 * This returns the PDF of ams, which is a marginal distribution if either c or p 
	 * are unconstrained (either num_p or num_c not equal to 1). 
	 * @return
	 */
	public HistogramFunction getPDF_ams() {
		if(num_ams == 0) {
			// changed this to still return a 1-bar histogram.
			return null;
		} else {
			if(D) System.out.println("ams: " + getMaxLikelihood_ams() + " " + min_ams + " " + max_ams + " " + num_ams);
			HistogramFunction hist = new HistogramFunction(min_ams, max_ams, num_ams);

			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							hist.add(amsIndex,  likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}

			String name = "PDF of ams-value";
			if(num_ams == 1)
				name += " (constrained)";
			else if(num_p !=1 || num_c != 1 || num_a != 1)
				name += " (marginal)";
			
			hist.setName(name);
			if(D) System.out.println("PDF of ams-value: totalTest = "+hist.calcSumOfY_Vals());

			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}

	/**
	 * This returns the PDF of a, which is a marginal distribution if either c or p 
	 * are unconstrained (either num_p or num_c not equal to 1). Null is returned if
	 * a is constrained (num_a=1). (how important is this? I want to show pdf with one value)
	 * @return
	 */
	public HistogramFunction getPDF_a() {
		if(num_a == 0) {
			return null;
		}
		else {
			if (D) System.out.println(min_a +" "+ max_a +" "+ num_a);
			
			HistogramFunction hist = new HistogramFunction(min_a, max_a, num_a);
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							hist.add(aIndex,  likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}
			String name = "PDF of a-value";
			if(num_a == 1)
				name += " (constrained)"; 
			else if(num_p !=1 || num_c != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) System.out.println("PDF of a-value: totalTest = "+hist.calcSumOfY_Vals());

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
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_p, max_p, num_p);
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							hist.add(pIndex, likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}
			String name = "PDF of p-value";
			if(num_p == 1)
				name += " (constrained)";
			else if(num_a !=1 || num_c != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) System.out.println("PDF of p-value: totalTest = "+hist.calcSumOfY_Vals() + " Max = " + hist.getMaxY());

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
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_c, max_c, num_c);
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							hist.add(cIndex, likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}
			String name = "PDF of c-value";
			if(num_c == 1)
				name += " (constrained)";
			else if(num_a !=1 || num_p != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) System.out.println("PDF of c-value: totalTest = "+hist.calcSumOfY_Vals());

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

			HistogramFunction hist = new HistogramFunction(min_logc, max_logc, num_c);
			
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							hist.add(cIndex, likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}
			String name = "PDF of logc-value";
			if(num_c == 1)
				name += " (constrained)";
			else if(num_a !=1 || num_p != 1)
				name += " (marginal)";
			hist.setName(name);
			if(D) System.out.println("PDF of logc-value: totalTest = "+hist.calcSumOfY_Vals());

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
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							double prevVal = hist2D.get(aIndex,pIndex);
							hist2D.set(aIndex,pIndex, prevVal+likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
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


	/**
	 * This returns a 2D PDF for a and c, which is a marginal distribution if p 
	 * is unconstrained (num_p not equal to 1). Null is returned if either
	 * a or c are constrained (num_a=1 or num_c=1).
	 * this one isn't ever used because c is log-normal. Use get2D_PDF_for_a_and_logc() instead.
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_a_and_c() {
		if(num_a == 1 || num_c == 1) {
			return null;
		}
		else {
			double delta_c = (max_c - min_c)/((double)num_c - 1);
			double delta_a = (max_a - min_a)/((double)num_a - 1);

			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_a, num_c, min_a, min_c, delta_a, delta_c);
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {

					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							double prevVal = hist2D.get(aIndex,cIndex);
							hist2D.set(aIndex,cIndex, prevVal+likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}

			if(D) {
				System.out.println("2D PDF of a vs c: totalTest = "+hist2D.getSumZ());
			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}

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
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							double prevVal = hist2D.get(aIndex,cIndex);
							hist2D.set(aIndex,cIndex, prevVal+likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}
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
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					for(int pIndex=0;pIndex<num_p;pIndex++) {
						for(int cIndex=0;cIndex<num_c;cIndex++) {
							double prevVal = hist2D.get(cIndex,pIndex);
							hist2D.set(cIndex,pIndex, prevVal+likelihood[amsIndex][aIndex][pIndex][cIndex]);
						}
					}
				}
			}
			if(D) {
				System.out.println("2D PDF of logc vs p: totalTest = "+hist2D.getSumZ());
			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}

	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}



}
