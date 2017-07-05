package scratch.aftershockStatisticsETAS;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.primitives.Doubles;


/**
 * This computes an ETAS aftershock model from aftershock data and with an assumed time-dependent 
 * magnitude completeness model described in the constructor.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests in order to cover all cases.
 *
 * @author van der Elst
 *
 */
public class ETAS_AftershockModel_SequenceSpecific extends ETAS_AftershockModel implements UnivariateFunction {
	
	Boolean D=false;	// debug flag
//	double alpha = 1;
//	double mu = 1;
	double rmax;
	double ams, a, k, p, c;	
//	
	
	
	/**
	 * PLACEHOLDER: NOT IMPLEMENTED.
	 * Use this constructor to apply a time-independent magnitude of completeness.  This is faster because it
	 * uses an analytical solution for the integral. 
	 * 
	 */
	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
				double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
				double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
				double maxMag, int maxGenerations, int nSims) {
		
		this(mainshock, aftershocks, Double.POSITIVE_INFINITY, aVec, pVec, cVec, alpha, b, refMag,
				dataStartTimeDays, dataEndTimeDays, forecastMinDays, forecastMaxDays, 
				maxMag, maxGenerations, nSims);
	}

	
	
	/**
	 * This solves for the ETAS parameters from the given mainShock, aftershockList,
	 * and other specifications as described below, and for a time-dependent magnitude of completeness
	 * model defined by a maximum rate. 
	 * Likelihood values are normalized so they sum to 1.0 over the range of parameter-values specified.
	 * @param mainShock
	 * @param aftershockList - events with mag below magComplete will be filtered out
	 * @param rmax - maximum rate (default 200 events/day) assume events in excess of this rate are missed.
	 * @param b - assumed b value
	 * @param min_a \
	 * @param max_a  | - range of a-values for grid search (set min=max and num=1 to constraint to single value)
	 * @param num_a /
	 * @param min_p \
	 * @param max_p  | - range of p-values for grid search (set min=max and num=1 to constraint to single value)
	 * @param num_p /
	 * @param min_c \
	 * @param max_c  | - range of c-values for grid search (set min=max and num=1 to constraint to single value)
	 * @param num_c /
	 */
	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainShock, ObsEqkRupList aftershockList,	double rmax,
			 								double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			 								double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			 								double maxMag, int maxGenerations, int nSims) {
		
		// check range values
		if(num_a == 1 && min_a != max_a) {
			throw new RuntimeException("Problem: num_a == 1 && min_a != max_a");
		}
		if(num_p == 1 && min_p != max_p) {
			throw new RuntimeException("Problem: num_p == 1 && min_p != max_p");
		}
		if(num_c == 1 && min_c != max_c) {
			throw new RuntimeException("Problem: num_c == 1 && min_c != max_c");
		}
		if(min_a > max_a) {
			throw new RuntimeException("Problem: min_a > max_a");
		}
		if(min_p > max_p) {
			throw new RuntimeException("Problem: min_p > max_p");
		}
		if(min_c > max_c) {
			throw new RuntimeException("Problem: min_c > max_c");
		}

		this.min_a = aVec[0];
		this.max_a = aVec[aVec.length-1];
		this.num_a = aVec.length;
		this.min_p = pVec[0];
		this.max_p = pVec[pVec.length-1];
		this.num_p = pVec.length;
		this.min_c = cVec[0];
		this.max_c = cVec[cVec.length-1];
		this.num_c = cVec.length;
		
		this.a_vec = aVec;
		this.p_vec = pVec;
		this.c_vec = cVec;
		
		this.b = b;
		this.alpha = alpha;
		this.refMag = refMag;
		this.magComplete = refMag;
		this.aftershockList=aftershockList;
//		this.aftershockList = new ObsEqkRupList();
		this.mainShock=mainShock;
		this.magMain = mainShock.getMag();
		this.dataStartTimeDays=dataStartTimeDays;
		this.dataEndTimeDays=dataEndTimeDays;
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.rmax = rmax;
		this.nSims = nSims;
		this.maxMag = maxMag;
		this.maxGenerations = maxGenerations;
		

//		if(num_a>1) // otherwise defaults to zero
//			delta_a = (max_a-min_a)/((double)num_a - 1.);
//		if(num_p>1)
//			delta_p = (max_p-min_p)/((double)num_p - 1.);
//		if(num_c>1)
//			delta_c = (max_c-min_c)/((double)num_c - 1.);

		
		
		if(D) {
			System.out.println("a-values range:\t"+min_a+"\t"+max_a+"\t"+num_a);
			System.out.println("p-values range:\t"+min_p+"\t"+max_p+"\t"+num_p+"\t");
			System.out.println("log c-values range:\t"+min_c+"\t"+max_c+"\t"+num_c+"\t");
			System.out.println("rmax:\t"+rmax);
			System.out.println("magComplete:\t"+magComplete);
		}
		
//		if(Double.isNaN(capG))
		
		// Find the max like parameters by grid search
		System.out.println("finding maximum likelihood parameters...");
		getLikelihoodMatrixGrid();
			
	
//		else
			//computeSequenceSpecificParams();
		
		if (D) {
			System.out.println("testTotalLikelihood="+testTotalLikelihood);
			System.out.println("getMaxLikelihood_a()="+getMaxLikelihood_a());
			System.out.println("getMaxLikelihood_p()="+getMaxLikelihood_p());
			System.out.println("getMaxLikelihood_c()="+getMaxLikelihood_c());
		}
		
		
		// get aftershock times and mags and store as simple doubles[]
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = getAftershockMags(aftershockList.getRupsAboveMag(magComplete));

		List<double[]> sortedEQlist = new ArrayList<double[]>();

		System.out.println(relativeEventTimes.length + " "+ dataEndTimeDays +" "+ magComplete);

		for(int i = 0; i < relativeEventTimes.length; i++){
			double[] temp = new double[]{relativeEventTimes[i], magAftershocks[i]};
			if(temp[0] <= dataEndTimeDays)
				sortedEQlist.add(temp);
		}

		//sort double[] of times and magnitudes
		Collections.sort(sortedEQlist, new java.util.Comparator<double[]>() {
			public int compare(double[] a, double[] b) {
				return Double.compare(a[0], b[0]);
			}
		});

		for(int i = 0; i < relativeEventTimes.length; i++){
			double[] temp = sortedEQlist.get(i);
			relativeEventTimes[i] = temp[0];
			magAftershocks[i] = temp[1];
			//					sortedEQlist.add(temp);
		}

		this.magAftershocks = magAftershocks;
		this.relativeTimeAftershocks = relativeEventTimes;
		
//		// DEBUG
//		for(int i = 0; i < num_a; i++){
//			for(int j = 0; j < num_p; j++){
//				for(int k = 0; k < num_c; k++){
//					System.out.format("%4.2f ", this.likeArray[i][j][k]);
//				}
//				System.out.format("\n");
//			}
//			System.out.format("\n\n");
//		}
//		
		
//		
		// generate forecast object
		computeNewForecast(dataStartTimeDays, dataEndTimeDays, forecastMinDays, forecastMaxDays, nSims);
		
	}		

	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
		System.out.println("Data/Forecast duration: " + dataMinDays +" "+ dataMaxDays +" "+ forecastMinDays +" "+ forecastMaxDays +" "+ nSims);
		System.out.println("Params: "+ getProductivityMag() +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete);
		
		ETAScatalog simulatedCatalog = new ETAScatalog(a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
				mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
		
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.simulatedCatalog = simulatedCatalog;
		this.nSims = nSims;
	}

	
	/**
	 * Get likelihood matrix with no time dependent Mc. Checks for supercriticality and gives a warning if too many supercritical parameter sets are found;
	 */
	private void getLikelihoodMatrixGrid() {
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = getAftershockMags(aftershockList.getRupsAboveMag(magComplete));
				
		List<double[]> sortedEQlist = new ArrayList<double[]>();
		
		for(int i = 0; i < relativeEventTimes.length; i++){
			double[] temp = new double[]{relativeEventTimes[i], magAftershocks[i]};
			sortedEQlist.add(temp);
		}
		
		//sort times and magnitudes
		Collections.sort(sortedEQlist, new java.util.Comparator<double[]>() {
		    public int compare(double[] a, double[] b) {
		        return Double.compare(a[0], b[0]);
		    }
		});

		for(int i = 0; i < relativeEventTimes.length; i++){
			double[] temp = sortedEQlist.get(i);
			relativeEventTimes[i] = temp[0];
			magAftershocks[i] = temp[1];
			sortedEQlist.add(temp);
		}

		// instantiate two likelihood matrices -- one to hold subcritical likelihoods, one to hold supercritical likelihoods (for bookkeeping).
		likelihood = new double[num_a][num_p][num_c];
		amsMatrix = new double[num_a][num_p][num_c];
		double[][][] superCriticalLikelihood = new double[num_a][num_p][num_c];
		double[][][] subCriticalLikelihood = new double[num_a][num_p][num_c];
		
		double maxVal= Double.NEGATIVE_INFINITY;
		long startTime = System.currentTimeMillis();
		double warnTime = 3;
		boolean longRunFlag = false;
		
		double logLike;
		boolean subCritFlag;
		
//		System.out.println(mainShock.getMag());
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			double a = a_vec[aIndex];
			if(!longRunFlag){
				double toc = (System.currentTimeMillis() - startTime) / 1000;
				if(toc > warnTime){
					longRunFlag = true;
					double timeEstimate = (double)toc/(aIndex+1) * num_a;
					System.out.format("This might take a while. Probably about %d seconds...\n", (int) timeEstimate);
				}
			}
			
			for(int pIndex=0;pIndex<num_p;pIndex++) {
				double p = p_vec[pIndex];
				
				for(int cIndex=0;cIndex<num_c;cIndex++) {
					double c = c_vec[cIndex];

					//check for supercritical parameters over the forecast time window
					double n;	//branching ratio
					if( p == 1)
						n = b * Math.log(10) * (maxMag - magComplete) * Math.pow(10, a) * ( Math.log(forecastMaxDays + c) - Math.log(c) );
					else
						n = b * Math.log(10) * (maxMag - magComplete) * Math.pow(10, a)/(1-p) * ( Math.pow(forecastMaxDays + c, 1-p) - Math.pow(c, 1-p) );
					subCritFlag = ( n<1 );
//					System.out.println(n);	//debug

					// optimize the mainshock productivity magnitude for this a p c, but do it behind the scenes
					double sigma_ams = 0.4; //set this as a parameter somewhere
					int num_ams = 101;	//set this parameter somewhere
					double mean_ams = a;
					double min_ams = mean_ams - 3d*sigma_ams;
					double max_ams = mean_ams + 3d*sigma_ams;
										
					double[] ams_vec = ETAS_StatsCalc.linspace(min_ams, max_ams, num_ams); 
					
					double[] logLike_vec = new double[num_ams];
					double[] amsLogLike_vec = new double[num_ams];
					
					// compute likelihoods for range of ms productivities 
					for(int amsIndex = 0 ; amsIndex < num_ams ; amsIndex++){
						double ams = ams_vec[amsIndex];
						double k = Math.pow(10, ams + alpha*(mainShock.getMag()-magComplete) );
						
						//ETAS log-likelihood component
						logLike_vec[amsIndex] = getLogLikelihoodForETASParams( k,  a,  p,  c,  alpha,  magComplete,
								dataStartTimeDays,  dataEndTimeDays, magAftershocks, relativeEventTimes); 
						
						//MS mag log-likelihood component
						amsLogLike_vec[amsIndex] = -1d/2d * Math.log(2d*Math.PI*sigma_ams*sigma_ams) - (ams - mean_ams)*(ams - mean_ams ) / (2d*sigma_ams*sigma_ams); 
//								1/Math.sqrt(2*Math.PI)/sigma_ams * Math.exp( -(ams - mean_ams) * (ams - mean_ams ) / (2*sigma_ams*sigma_ams) );
						
					}
					
					//find the max likelihood index
					double maxLike = Double.NEGATIVE_INFINITY;
					int maxLLIndex = 0;
					for (int amsIndex = 0; amsIndex < num_ams; amsIndex++){
						logLike = logLike_vec[amsIndex] + amsLogLike_vec[amsIndex];
						if(logLike > maxLike){
							maxLike = logLike;
							maxLLIndex = amsIndex;
						}
					}
					logLike = maxLike;
					amsMatrix[aIndex][pIndex][cIndex] = ams_vec[maxLLIndex];

//					DEBUG
//					System.out.println(aIndex + " " + ams_vec[maxLLIndex] + " " + sigma_ams +  " "  + logLike_vec[maxLLIndex] + " " + amsLogLike_vec[maxLLIndex]);
					
					
					// fill out the likelihood matrices with the joint likelihood
					if(Doubles.isFinite(logLike)){
						if (subCritFlag){
							likelihood[aIndex][pIndex][cIndex] = logLike;
							subCriticalLikelihood[aIndex][pIndex][cIndex] = logLike;
							superCriticalLikelihood[aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;

							if(maxVal<logLike ) {
								maxVal=logLike;
								max_a_index=aIndex;
								max_p_index=pIndex;
								max_c_index=cIndex;
							}
						} else {
							likelihood[aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
							subCriticalLikelihood[aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
							superCriticalLikelihood[aIndex][pIndex][cIndex] = logLike;
						}
					}else{
						likelihood[aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
						subCriticalLikelihood[aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
						superCriticalLikelihood[aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
					}

				}
			}
		}
		

		// convert array from log-likelihood to likelihood
		testTotalLikelihood = convertLogLikelihoodArrayToLikelihood(maxVal);
		System.out.println("Total likelihood  = " + testTotalLikelihood); //debug
		
		//measure the proportion of supercritical combinations
		double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
		double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
		double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
		
		double toc = (System.currentTimeMillis() - startTime) / 1000;
		System.out.format("Grid search took %d seconds.\n", (int)toc);
		System.out.format("%3.2f percent of the solution space is subcritical.\n", fractionSubCritical*100);
		System.out.format("Max likelihood mainshock productivity manigude is %2.2f\n" , getProductivityMag());
	}

	/*
	 * To make this run fast, give it sorted vectors of magnitudes and times. Actually, it will not work otherwise.
	 * 
	 * todo: check for sorted vectors and sort if not. 
	 */
	public double getLogLikelihoodForETASParams(double k, double a, double p, double c, double alpha, double Mc,
			double tMinDays, double tMaxDays, double[] aftershockMagnitudes, double[] relativeEventTimes) {
		
		double LL;
		double timeIntegral = Double.NaN;
		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
				
		//compute total number at end of fit window due to mainshock (NEEDS TO BE FIXED FOR RMAX)
		// compute completeness time for this mainshock
		double tcompleteMS = Math.pow(k/rmax,1d/p) - c;
		if(tcompleteMS < 0) tcompleteMS = 0;
		
		double Ntot;
		if (p == 1)
			timeIntegral = Math.log(tMaxDays + c) - Math.log(tMinDays + c);
		else
			timeIntegral = (Math.pow(tMaxDays + c, 1-p) - Math.pow(tMinDays + c, 1-p)) / (1-p);
		Ntot = k*timeIntegral;		//mainshock contribution
//		System.out.println("mainshock: " + Ntot + " " + Nas +" "+ tMinDays + " " + tMaxDays + " " + Mref); //debug
	
		double[] productivity = new double[Nas];
		double[] lambda = new double[Nas];
		double[] tcomplete = new double[Nas];
		
		for(int i=0; i<Nas; i++){
			// compute productivity for this aftershock
			productivity[i] = Math.pow(10, a + alpha*(aftershockMagnitudes[i] - Mc));	//productivity of this aftershock

			// compute completeness time for this aftershock
			tcomplete[i] = Math.pow(productivity[i]/rmax,1d/p) - c;
			if(tcomplete[i] < 0) tcomplete[i]=0;
			
			// compute intensity at this moment due to previous earthquakes
			lambda[i] = k/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
			if(lambda[i] > rmax) lambda[i] = rmax;
			
			for(int j = 0; j < i; j++){//from the aftershocks
				if(relativeEventTimes[j] < relativeEventTimes[i]){
					double lambdasub = productivity[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p);
					if(lambdasub > rmax) lambdasub = rmax;
					lambda[i] += lambdasub;
				}
			}
	
			//compute number at end of window due to this aftershock
			if(relativeEventTimes[i] < tMaxDays){
				double Nbase, Ntotsub, tIntegralStart;
				
				if(tcomplete[i] > 0){	//if there's an incomplete period
					if(relativeEventTimes[i] + tcomplete[i] < tMaxDays){	//if the incomplete period ends before catalog end
						//first compute the integral over the maxed out period
						Nbase = rmax * tcomplete[i];
						tIntegralStart = tcomplete[i] + relativeEventTimes[i]; //set integral start time to the completeness time
					}else{ //if there's an incomplete period that lasts the whole duration of the catalog
						// compute the integral over the maxed out period
						Nbase = rmax * (tMaxDays - relativeEventTimes[i]);
						tIntegralStart = tMaxDays; //set integral start time to catalog end time
					}
//					System.out.format("rmax exceeded for %f days, %f events\n", tcomplete[i], Nbase); 
				}else{//otherwise, Nbase = 0 and tIntegralStart is the same as the earthquake origin time
					Nbase = 0;
					tIntegralStart = relativeEventTimes[i];
				}
			
				//compute the remaining integral (may be zero)	
				if(p == 1)
					timeIntegral = Math.log(tMaxDays - relativeEventTimes[i] + c) - Math.log(tIntegralStart - relativeEventTimes[i] + c); 
				else
					timeIntegral = (Math.pow(tMaxDays - relativeEventTimes[i] + c, 1-p) - Math.pow(tIntegralStart - relativeEventTimes[i] + c, 1-p)) / (1-p);
				
				Ntotsub = productivity[i]*timeIntegral;
				
				Ntot += Nbase + Ntotsub; 
			}
			//System.out.format(" %d", (int) Ntot); 
		}
			
		//compute likelihood (sum of the log-lambdas at times of all aftershocks minus the total number expected)
		LL = -Ntot;
		for(int i = 0; i<Nas; i++)
			LL += Math.log(lambda[i]);
		
		//debug mode
//		System.out.println("for k="+k+", a="+a+", p="+p+", c="+c+", Nsim="+ Ntot +", Nas="+Nas+", LL="+ LL);
		
		return LL;
	}
	
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	public double getProductivityMag(){
		return magMain + amsMatrix[max_a_index][max_p_index][max_c_index] - a_vec[max_a_index];
	}

@Override
public double value(double arg0) {
	// TODO Auto-generated method stub
	return 0;
}

}
