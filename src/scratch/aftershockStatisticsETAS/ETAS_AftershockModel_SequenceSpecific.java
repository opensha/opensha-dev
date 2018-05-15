package scratch.aftershockStatisticsETAS;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.gui.infoTools.CalcProgressBar;

import com.google.common.base.Stopwatch;
import com.google.common.primitives.Doubles;

/**
 * This computes an ETAS aftershock model from aftershock data and with an assumed time-dependent 
 * magnitude completeness model described in the constructor.
 * 
 */
public class ETAS_AftershockModel_SequenceSpecific extends ETAS_AftershockModel {
	
	Boolean D = false;	// debug flag
	private boolean fitMSProductivity;
	private volatile boolean stopRequested;
	private volatile boolean pauseRequested;
	private ETAS_AftershockModel_Generic priorModel;

	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double[] amsVec, double amsSigma, double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			double minMag, double maxMag, int maxGenerations, int nSims, boolean fitMSProductivity) {

		this(mainshock, aftershocks, amsVec, amsSigma, aVec, pVec, cVec, alpha, b, refMag,
				dataStartTimeDays, dataEndTimeDays, forecastMinDays, forecastMaxDays, 
				minMag, maxMag, maxGenerations, nSims, fitMSProductivity, false, null, null);
	}

	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double[] amsVec, double amsSigma, double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			double minMag, double maxMag, int maxGenerations, int nSims, boolean fitMSProductivity, boolean timeDependentMc, ETAS_AftershockModel_Generic priorModel){

		this(mainshock, aftershocks,
				amsVec, amsSigma, aVec, pVec, cVec, alpha, b, refMag,
				dataStartTimeDays, dataEndTimeDays, forecastMinDays, forecastMaxDays, 
				minMag, maxMag, maxGenerations, nSims, fitMSProductivity, timeDependentMc,
				priorModel, null);	
	}

	/**
	 * This solves for the ETAS parameters from the given mainShock, aftershockList,
	 * and other specifications as described below
	 * @param mainShock
	 * @param aftershockList - events with mag below magComplete will be filtered out
	 * @param b - assumed b value
	 * @param amsVec: values of ams to use in grid search 
	 * @param aVec: values of ams to use in grid search
	**/
	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainShock, ObsEqkRupList aftershockList,
			 								double[] amsVec, double amsSigma, double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			 								double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			 								double magComplete, double maxMag, int maxGenerations, int nSims, boolean fitMSProductivity, boolean timeDependentMc,
			 								ETAS_AftershockModel_Generic priorModel,
			 								CalcProgressBar progress) {
		
		if(fitMSProductivity){
			this.min_ams = amsVec[0];
			this.max_ams = amsVec[amsVec.length-1];
			this.num_ams = amsVec.length;
			this.ams_vec = amsVec;
			this.sigma_ams = amsSigma;
		}else{
			this.min_ams = amsVec[0];
			this.max_ams = amsVec[amsVec.length-1];
			this.num_ams = 1;
			this.ams_vec = ETAS_StatsCalc.linspace(min_ams, max_ams, num_ams);
			this.sigma_ams = 0d;
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
		this.magComplete = magComplete;
		this.aftershockList=aftershockList;
		this.mainShock=mainShock;
		this.magMain = mainShock.getMag();
		this.dataStartTimeDays=dataStartTimeDays;
		this.dataEndTimeDays=dataEndTimeDays;
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.nSims = nSims;
		this.maxMag = maxMag;
		this.maxGenerations = maxGenerations;
		this.fitMSProductivity = fitMSProductivity;
		this.timeDependentMc = timeDependentMc;
		this.priorModel = priorModel;
		this.progress = progress;
		
		if(D) {
			System.out.println("ams-values range:\t"+min_ams+"\t"+max_ams+"\t"+num_ams);
			System.out.println("a-values range:\t"+min_a+"\t"+max_a+"\t"+num_a);
			System.out.println("p-values range:\t"+min_p+"\t"+max_p+"\t"+num_p+"\t");
			System.out.println("c-values range:\t"+min_c+"\t"+max_c+"\t"+num_c+"\t");
			System.out.println("Mc(t):\t" + timeDependentMc + " Mc = " + magComplete);
		}
		
		// Find the max like parameters by grid search
		if(D) System.out.println("finding maximum likelihood parameters...");
		if (timeDependentMc && (dataEndTimeDays > dataStartTimeDays) && priorModel != null) {
			if(D) System.out.println("Time-dependent Mc is active, the fit will be done twice...");
			this.ac = priorModel.ac;
			getLikelihoodMatrixGridFastMc();
			this.ac = getMaxLikelihood_ams();
			getLikelihoodMatrixGridFastMc();
		} else { 
			getLikelihoodMatrixGridFast();
		}
		
		if (D) {
			System.out.println("getMaxLikelihood_ams()="+getMaxLikelihood_ams());
			System.out.println("getMaxLikelihood_a()="+getMaxLikelihood_a());
			System.out.println("getMaxLikelihood_p()="+getMaxLikelihood_p());
			System.out.println("getMaxLikelihood_c()="+getMaxLikelihood_c());
		}
		
		
		// get aftershock times and mags and store as simple doubles[]
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = ETAS_StatsCalc.getAftershockMags(aftershockList.getRupsAboveMag(magComplete));

		List<double[]> sortedEQlist = new ArrayList<double[]>();
		
		if(D)	System.out.println("numAS: " + relativeEventTimes.length + " "+ dataEndTimeDays +" "+ magComplete);

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
		
		// generate forecast object
		computeNewForecast(dataStartTimeDays, dataEndTimeDays, forecastMinDays, forecastMaxDays, nSims);
		
	}		
	
	/**
	 * Get likelihood matrix with or without time dependent Mc. Checks for supercriticality 
	 * and gives a warning if too many supercritical parameter sets are found;
	 */
	private void getLikelihoodMatrixGridFastMc() {
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = ETAS_StatsCalc.getAftershockMags(aftershockList.getRupsAboveMag(magComplete));
				
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
		likelihood = new double[num_ams][num_a][num_p][num_c];
		double[][][][] superCriticalLikelihood = new double[num_ams][num_a][num_p][num_c];
		double[][][][] subCriticalLikelihood = new double[num_ams][num_a][num_p][num_c];
		double[][][][] priorLikelihood = new double[num_ams][num_a][num_p][num_c];
		
		if (priorModel != null && priorModel.likelihood != null)
			priorLikelihood = priorModel.get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, false);
		else{
//			timeDependentMc = false; // can't do time-dependent Mc in a timely fashion without fixing incompleteness to typical values
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++)
				for(int aIndex=0;aIndex<num_a;aIndex++)
					for(int pIndex=0;pIndex<num_p;pIndex++)
						for(int cIndex=0;cIndex<num_c;cIndex++){
							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = 1;
						}
		}
		
		double maxVal= Double.NEGATIVE_INFINITY;
		double logLike;
		boolean subCritFlag;
		
		double ams, a, k, kms, p, c, c0;
		double kc = Math.pow(10, ac);
		int amsIndex, aIndex, pIndex, cIndex;
		double timeIntegralMS, NtotAS, NtotMS, Ntot;
		
		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
		double productivityMS;
		double[] productivityAS = new double[Nas];
		double[] lambda = new double[Nas];
		double[] timeDecayMS = new double[Nas];
		double[] timeDecayAS = new double[Nas];
		double timeIntegral;
		
		
		//do productivities
		productivityMS = Math.pow(10, alpha*(mainShock.getMag()-magComplete) );
		for(int i=0; i<Nas; i++) //compute productivity for this aftershock
			productivityAS[i] = Math.pow(10, alpha*(magAftershocks[i] - magComplete));	//productivity of this aftershock


		double tStartIntegration;
	
		// set up timer/time estimator
		double toc, timeEstimate, n;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		String initialMessageString = "Estimating sequence-specific model. ";
		
		for(pIndex=0;pIndex<num_p;pIndex++) {
			p = p_vec[pIndex];

			for(cIndex=0;cIndex<num_c;cIndex++) {
				c0 = c_vec[cIndex];

				tStartIntegration = 0;
				
				//compute total number at end of fit window for mainshock (unscaled by a)
				c = c0*Math.pow(kc*productivityMS,1/p);
				if (p == 1){
					timeIntegralMS = Math.log(dataEndTimeDays + c) - Math.log(tStartIntegration + c);
				} else {
					timeIntegralMS = (Math.pow(dataEndTimeDays + c, 1-p) - Math.pow(tStartIntegration + c, 1-p)) / (1-p);
				}
				NtotMS = productivityMS*timeIntegralMS;
				
				//compute instantaneous intensities and total number for aftershocks (unscaled by a)
				NtotAS = 0;
				for(int i=0; i<Nas; i++){
					//compute intensity at this moment due to mainshock (unscaled by ams)
					c = c0*Math.pow(kc*productivityMS,1/p);
					timeDecayMS[i] = productivityMS/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
					
					//compute intensity at this moment due to previous aftershocks (unscaled by a)
					timeDecayAS[i] = 0;
					for(int j = 0; j < i; j++){
						if(relativeEventTimes[j] < relativeEventTimes[i])
							c = c0*Math.pow(kc*productivityAS[j],1/p);
							timeDecayAS[i] += productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p);	//from the aftershocks
					}
				
					tStartIntegration = relativeEventTimes[i];
					
					//compute total number at end of fit window due to this aftershock (unscaled by a)
					c = c0*Math.pow(kc*productivityAS[i],1/p);
					if(relativeEventTimes[i] < dataEndTimeDays){
						if(p == 1){
							timeIntegral = Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(tStartIntegration - relativeEventTimes[i] + c); 
						} else {
							timeIntegral = (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1-p) - Math.pow(tStartIntegration - relativeEventTimes[i] + c, 1-p)) / (1-p);
						}
						
						NtotAS += productivityAS[i]*timeIntegral;	//aftershock Contributions
					}
				}
				
				// loop over productivities
				for(amsIndex=0;amsIndex<num_ams;amsIndex++) {
					//			double ams;	//assigned in next loop
					ams = ams_vec[amsIndex];
					
					for(aIndex=0;aIndex<num_a;aIndex++) {
						a = a_vec[aIndex];

						logLike = 0;
						
						// now put in the productivity terms and compute likelihood
						kms = Math.pow(10,ams);
						k = Math.pow(10, a);
						
						Ntot = kms*NtotMS + k*NtotAS;
						
						logLike += -Ntot;
						
						for(int i=0; i<Nas; i++){
							//compute intensity at this moment due to previous earthquakes
							lambda[i] = kms*timeDecayMS[i]; //from the mainshock
							lambda[i] += k*timeDecayAS[i];	//from the aftershocks
							logLike += Math.log(lambda[i]);
						}
						
						//add prior regularization
						logLike += Math.log(priorLikelihood[amsIndex][aIndex][pIndex][cIndex]);
						
						//check for supercritical parameters over the forecast time window
						c = c0*Math.pow(kc,1/p);
						n = ETAS_StatsCalc.calculateBranchingRatio(a, p, c, alpha, b, forecastMaxDays, magComplete, maxMag);
						
						subCritFlag = ( n < 1 );

						
						// fill out the likelihood matrices with the joint likelihood
						if(Doubles.isFinite(logLike)){
							if (subCritFlag){
								likelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
								subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
								superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;

								if(maxVal<logLike ) {
									maxVal=logLike;
									max_ams_index=amsIndex;
									max_a_index=aIndex;
									max_p_index=pIndex;
									max_c_index=cIndex;
								}
							} else {
								likelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
								subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
								superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
							}
						}else{
							likelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
							subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
							superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
						}

						// run the timer to see how long this is going to take
						toc = watch.elapsed(TimeUnit.SECONDS);
						if(toc > warnTime){
							//check for stopRequested?
							if (stopRequested){
								System.out.println("Parameter estimation terminated prematurely.");
								return;
							}
								

							warnTime += 10;

							long count = (pIndex)*(num_c*num_ams*num_a) + (cIndex)*(num_ams*num_a) + (amsIndex)*(num_a) + aIndex;
							long total = (num_p*num_c*num_ams*num_a);
							timeEstimate = toc * total/count;
							System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
							initialMessageString = "...";
							if (progress != null){
								progress.updateProgress(count, total, String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100), (int) ((timeEstimate - toc))));
//								progress.setProgressMessage(String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc))));
								progress.repaint();
								try {
									Thread.sleep(100);
								} catch (InterruptedException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}
						}
					}
				}
			}
		}

		// if not fitting MS productivity, make sure the ams_vector reflects this constraint
		if(!fitMSProductivity){
			this.min_ams = getMaxLikelihood_a();
			this.max_ams = getMaxLikelihood_a();
			this.num_ams = 1;
			this.ams_vec = ETAS_StatsCalc.linspace(min_ams, max_ams, num_ams);
		};

		// convert array from log-likelihood to likelihood
		double testTotalLikelihood = convertLogLikelihoodArrayToLikelihood(maxVal); // this converts likelihood from log to linear
		toc = watch.elapsed(TimeUnit.SECONDS);
		
		if(D) System.out.format("Grid search took %d seconds.\n", (int)toc);
		if(D) System.out.println("Total likelihood  = " + testTotalLikelihood); //debug
		
		//measure the proportion of supercritical combinations
		double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
		double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
		double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
		
		toc = watch.elapsed(TimeUnit.SECONDS);
		
		if(D) System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
		if(D) System.out.format("Mainshock productivity magnitude: %2.2f\n" , getProductivityMag());
		watch.stop();
		this.epiLikelihood = likelihood;
	}

	/**
	 * Get likelihood matrix with or without time dependent Mc. Checks for supercriticality 
	 * and gives a warning if too many supercritical parameter sets are found;
	 */
	private void getLikelihoodMatrixGridFast() {
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = ETAS_StatsCalc.getAftershockMags(aftershockList.getRupsAboveMag(magComplete));
				
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
		likelihood = new double[num_ams][num_a][num_p][num_c];
		double[][][][] superCriticalLikelihood = new double[num_ams][num_a][num_p][num_c];
		double[][][][] subCriticalLikelihood = new double[num_ams][num_a][num_p][num_c];
		double[][][][] priorLikelihood = new double[num_ams][num_a][num_p][num_c];
		
		if (priorModel != null && priorModel.likelihood != null)
			priorLikelihood = priorModel.get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, false);
		else{
//			timeDependentMc = false; // can't do time-dependent Mc in a timely fashion without fixing incompleteness to typical values
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++)
				for(int aIndex=0;aIndex<num_a;aIndex++)
					for(int pIndex=0;pIndex<num_p;pIndex++)
						for(int cIndex=0;cIndex<num_c;cIndex++){
							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = 1;
						}
		}
		
		double maxVal= Double.NEGATIVE_INFINITY;
		double logLike;
		boolean subCritFlag;
		
		double ams, a, k, kms, p, c;
		int amsIndex, aIndex, pIndex, cIndex;
		double timeIntegralMS, NtotAS, NtotMS, Ntot;
		
		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
		double productivityMS;
		double[] productivityAS = new double[Nas];
		double[] lambda = new double[Nas];
		double[] timeDecayMS = new double[Nas];
		double[] timeDecayAS = new double[Nas];
		double timeIntegral;
		
		
		//do productivities
		productivityMS = Math.pow(10, alpha*(mainShock.getMag()-magComplete) );
		for(int i=0; i<Nas; i++) //compute productivity for this aftershock
			productivityAS[i] = Math.pow(10, alpha*(magAftershocks[i] - magComplete));	//productivity of this aftershock


		double tStartIntegration;
	
		// set up timer/time estimator
		double toc, timeEstimate, n;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		String initialMessageString = "Estimating sequence-specific model. ";
		
		for(pIndex=0;pIndex<num_p;pIndex++) {
			p = p_vec[pIndex];

			for(cIndex=0;cIndex<num_c;cIndex++) {
				c = c_vec[cIndex];

				tStartIntegration = 0;
				
				
				//compute total number at end of fit window for mainshock (unscaled by a)
				if (p == 1){
					timeIntegralMS = Math.log(dataEndTimeDays + c) - Math.log(tStartIntegration + c);
				} else {
					timeIntegralMS = (Math.pow(dataEndTimeDays + c, 1-p) - Math.pow(tStartIntegration + c, 1-p)) / (1-p);
				}
				NtotMS = productivityMS*timeIntegralMS;
				
				//compute instantaneous intensities and total number for aftershocks (unscaled by a)
				NtotAS = 0;
				for(int i=0; i<Nas; i++){
					//compute intensity at this moment due to mainshock (unscaled by ams)
					timeDecayMS[i] = productivityMS/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
					
					//compute intensity at this moment due to previous aftershocks (unscaled by a)
					timeDecayAS[i] = 0;
					for(int j = 0; j < i; j++){
						if(relativeEventTimes[j] < relativeEventTimes[i])
							timeDecayAS[i] += productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p);	//from the aftershocks
					}
				
					tStartIntegration = relativeEventTimes[i];
					
					//compute total number at end of fit window due to this aftershock (unscaled by a)
					if(relativeEventTimes[i] < dataEndTimeDays){
						if(p == 1){
							timeIntegral = Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(tStartIntegration - relativeEventTimes[i] + c); 
						} else {
							timeIntegral = (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1-p) - Math.pow(tStartIntegration - relativeEventTimes[i] + c, 1-p)) / (1-p);
						}
						
						NtotAS += productivityAS[i]*timeIntegral;	//aftershock Contributions
					}
				}
				
				// loop over productivities
				for(amsIndex=0;amsIndex<num_ams;amsIndex++) {
					//			double ams;	//assigned in next loop
					ams = ams_vec[amsIndex];
					
					for(aIndex=0;aIndex<num_a;aIndex++) {
						a = a_vec[aIndex];

						logLike = 0;
						
						// now put in the productivity terms and compute likelihood
						kms = Math.pow(10,ams);
						k = Math.pow(10, a);
						
						Ntot = kms*NtotMS + k*NtotAS;
						
						logLike += -Ntot;
						
						for(int i=0; i<Nas; i++){
							//compute intensity at this moment due to previous earthquakes
							lambda[i] = kms*timeDecayMS[i]; //from the mainshock
							lambda[i] += k*timeDecayAS[i];	//from the aftershocks
							logLike += Math.log(lambda[i]);
						}
						
						//add prior regularization
						logLike += Math.log(priorLikelihood[amsIndex][aIndex][pIndex][cIndex]);
						
						//check for supercritical parameters over the forecast time window
						n = ETAS_StatsCalc.calculateBranchingRatio(a, p, c, alpha, b, forecastMaxDays, magComplete, maxMag);
						
						subCritFlag = ( n < 1 );

						
						// fill out the likelihood matrices with the joint likelihood
						if(Doubles.isFinite(logLike)){
							if (subCritFlag){
								likelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
								subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
								superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;

								if(maxVal<logLike ) {
									maxVal=logLike;
									max_ams_index=amsIndex;
									max_a_index=aIndex;
									max_p_index=pIndex;
									max_c_index=cIndex;
								}
							} else {
								likelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
								subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
								superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
							}
						}else{
							likelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
							subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
							superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
						}

						// run the timer to see how long this is going to take
						toc = watch.elapsed(TimeUnit.SECONDS);
						if(toc > warnTime){
							//check for stopRequested?
							if (stopRequested){
								System.out.println("Parameter estimation terminated prematurely.");
								return;
							}
								

							warnTime += 10;

							long count = (pIndex)*(num_c*num_ams*num_a) + (cIndex)*(num_ams*num_a) + (amsIndex)*(num_a) + aIndex;
							long total = (num_p*num_c*num_ams*num_a);
							timeEstimate = toc * total/count;
							System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
							initialMessageString = "...";
							if (progress != null){
								progress.updateProgress(count, total, String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100), (int) ((timeEstimate - toc))));
//								progress.setProgressMessage(String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc))));
								progress.repaint();
								try {
									Thread.sleep(100);
								} catch (InterruptedException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}
						}
					}
				}
			}
		}

		// if not fitting MS productivity, make sure the ams_vector reflects this constraint
		if(!fitMSProductivity){
			this.min_ams = getMaxLikelihood_a();
			this.max_ams = getMaxLikelihood_a();
			this.num_ams = 1;
			this.ams_vec = ETAS_StatsCalc.linspace(min_ams, max_ams, num_ams);
		};

		// convert array from log-likelihood to likelihood
		double testTotalLikelihood = convertLogLikelihoodArrayToLikelihood(maxVal); // this converts likelihood from log to linear
		toc = watch.elapsed(TimeUnit.SECONDS);
		
		if(D) System.out.format("Grid search took %d seconds.\n", (int)toc);
		if(D) System.out.println("Total likelihood  = " + testTotalLikelihood); //debug
		
		//measure the proportion of supercritical combinations
		double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
		double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
		double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
		
		toc = watch.elapsed(TimeUnit.SECONDS);
		
		if(D) System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
		if(D) System.out.format("Mainshock productivity magnitude: %2.2f\n" , getProductivityMag());
		watch.stop();
		this.epiLikelihood = likelihood;
	}
	
	
//	/**
//	 * Get likelihood matrix with no shortcuts (computing p and c terms inside a and ams loops) SLOOOOOOW!
//	 */
//	private void getLikelihoodMatrixGridSlow() {
//		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
//		double[] magAftershocks = ETAS_StatsCalc.getAftershockMags(aftershockList.getRupsAboveMag(magComplete));
//
//		List<double[]> sortedEQlist = new ArrayList<double[]>();
//
//		for(int i = 0; i < relativeEventTimes.length; i++){
//			double[] temp = new double[]{relativeEventTimes[i], magAftershocks[i]};
//			sortedEQlist.add(temp);
//		}
//
//		//sort times and magnitudes
//		Collections.sort(sortedEQlist, new java.util.Comparator<double[]>() {
//			public int compare(double[] a, double[] b) {
//				return Double.compare(a[0], b[0]);
//			}
//		});
//
//		for(int i = 0; i < relativeEventTimes.length; i++){
//			double[] temp = sortedEQlist.get(i);
//			relativeEventTimes[i] = temp[0];
//			magAftershocks[i] = temp[1];
//			sortedEQlist.add(temp);
//		}
//
//		// instantiate two likelihood matrices -- one to hold subcritical likelihoods, one to hold supercritical likelihoods (for bookkeeping).
//		likelihood = new double[num_ams][num_a][num_p][num_c];
//		double[][][][] superCriticalLikelihood = new double[num_ams][num_a][num_p][num_c];
//		double[][][][] subCriticalLikelihood = new double[num_ams][num_a][num_p][num_c];
//		double[][][][] priorLikelihood = new double[num_ams][num_a][num_p][num_c];
//
//		if (priorModel != null && priorModel.likelihood != null){
//			priorLikelihood = priorModel.get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, false);
//		}else{
//			for(int amsIndex=0;amsIndex<num_ams;amsIndex++)
//				for(int aIndex=0;aIndex<num_a;aIndex++)
//					for(int pIndex=0;pIndex<num_p;pIndex++)
//						for(int cIndex=0;cIndex<num_c;cIndex++)
//							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = 1;
//		}
//
//		double maxVal= Double.NEGATIVE_INFINITY;
//		double logLike;
//		boolean subCritFlag;
//
//		double ams, a, k, kms, p, c;
//		//	double kDefault = Math.pow(10, mean_a);
//		int amsIndex, aIndex, pIndex, cIndex;
//		double timeIntegral, Ntot;
//
//		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
//		double productivityMS;
//		double[] productivityAS = new double[Nas];
//		double[] lambda = new double[Nas];
//		double rateAS;
//
//
//		//do productivities
//		productivityMS = Math.pow(10, alpha*(mainShock.getMag()-magComplete) );
//		for(int i=0; i<Nas; i++) //compute productivity for this aftershock
//			productivityAS[i] = Math.pow(10, alpha*(magAftershocks[i] - magComplete));	//productivity of this aftershock
//
//		// set up timer/time estimator
//		double toc, timeEstimate, n;
//		Stopwatch watch = Stopwatch.createStarted();
//		int warnTime = 3;
//		String initialMessageString = "Computing aftershock model parameters. ";
//		// simpson integrator
//		//	SimpsonIntegrator integral = new SimpsonIntegrator();
//		//	RateIntegrand integrand = new RateIntegrand();
//
//		double timeCompleteMS = 0;
//
//		if(D) System.out.println("Mc(t) = " + timeDependentMc +" Data Start: " + dataStartTimeDays + " Data End: " + dataEndTimeDays + " Mc = " + magComplete);
//		
//		double tStartIntegration;
//		// loop over productivities
//		for(amsIndex=0;amsIndex<num_ams;amsIndex++) {
//			ams = ams_vec[amsIndex];
//
//			for(aIndex=0;aIndex<num_a;aIndex++) {
//				a = a_vec[aIndex];
//
//				// productivity scaling
//				kms = Math.pow(10,ams);
//				k = Math.pow(10, a);
//
//				for(pIndex=0;pIndex<num_p;pIndex++) {
//					p = p_vec[pIndex];
//
//					for(cIndex=0;cIndex<num_c;cIndex++) {
//						c = c_vec[cIndex];
//
//						logLike = 0;
//
//						tStartIntegration = timeCompleteMS;
//
//						//compute total number at end of fit window for mainshock 
//						if (p == 1)
//							timeIntegral = kms*productivityMS * (Math.log(dataEndTimeDays + c) - Math.log(tStartIntegration + c));
//						else
//							timeIntegral = kms*productivityMS * (Math.pow(dataEndTimeDays + c, 1d-p) - Math.pow(tStartIntegration + c, 1d-p)) / (1d-p);
//						Ntot = timeIntegral;
//
//						//compute instantaneous intensities and total numbers for aftershocks 
//						for(int i=0; i<Nas; i++){
//							// compute total number for this aftershock 
//
//							// start either at event time or complete time
//							tStartIntegration = Math.max(timeCompleteMS, relativeEventTimes[i]);
//
//							if (p == 1)
//								timeIntegral = k*productivityAS[i] * (Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(tStartIntegration - relativeEventTimes[i] + c));
//							else
//								timeIntegral = k*productivityAS[i] * (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1d-p) - Math.pow(tStartIntegration  - relativeEventTimes[i] + c, 1d-p)) / (1d-p);
//
//							Ntot += timeIntegral;
//
//							//compute intensity at this moment due to mainshock (scaled by ams)
//
//							if (relativeEventTimes[i] > timeCompleteMS){
//								rateAS = kms*productivityMS/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
//								lambda[i] = rateAS;
//
//								//compute intensity at this moment due to previous aftershocks (unscaled by a)
//								for(int j = 0; j < i; j++){
//									if(relativeEventTimes[j] < relativeEventTimes[i]){
//										rateAS = k*productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p); //from the aftershocks
//										lambda[i] += rateAS;
//									}
//								}
//							}
//						}
//
//						logLike -= Ntot;
//						for(int i=0; i<Nas; i++){
//							//compute intensity at this moment due to previous earthquakes
//							if (relativeEventTimes[i] > timeCompleteMS)
//								logLike += Math.log(lambda[i]);
//						}
//
//						//add prior regularization
//						logLike += Math.log(priorLikelihood[amsIndex][aIndex][pIndex][cIndex]);
//
//						//check for supercritical parameters over the forecast time window
//						n = ETAS_StatsCalc.calculateBranchingRatio(a, p, c, alpha, b, forecastMaxDays, magComplete, maxMag);
//
//						subCritFlag = ( n<1 );
//
//						// fill out the likelihood matrices with the joint likelihood
//						if(Doubles.isFinite(logLike)){
//							if (subCritFlag){
//								likelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
//								subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
//								superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
//
//								if(maxVal<logLike ) {
//									maxVal=logLike;
//									max_ams_index=amsIndex;
//									max_a_index=aIndex;
//									max_p_index=pIndex;
//									max_c_index=cIndex;
//								}
//							} else {
//								likelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
//								subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
//								superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = logLike;
//							}
//						}else{
//							likelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
//							subCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
//							superCriticalLikelihood[amsIndex][aIndex][pIndex][cIndex] = Double.NEGATIVE_INFINITY;
//						}
//
//						// run the timer to see how long this is going to take
//						toc = watch.elapsed(TimeUnit.SECONDS);
//						if(toc > warnTime){
//							warnTime += 10;
//
//							long count = ((amsIndex)*(num_c*num_p*num_a) + (aIndex)*(num_p*num_c) + (pIndex)*(num_c) + cIndex);
//							long total = (num_p*num_c*num_ams*num_a);
//							timeEstimate = toc * total/count;
//							System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
//							initialMessageString = "...";
//							
//							if (progress != null){
//								//									progress.updateProgress(count, total, String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc)/1000)));
//								progress.setProgressMessage(String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc))));
//								progress.pack();
//							}
//						}
//
//					}
//				}
//			}
//		}
//
//
//		// if not fitting MS productivity, make sure the ams_vector reflects this constraint
//		if(!fitMSProductivity){
//			this.min_ams = getMaxLikelihood_a();
//			this.max_ams = getMaxLikelihood_a();
//			this.num_ams = 1;
//			this.ams_vec = ETAS_StatsCalc.linspace(min_ams, max_ams, num_ams);
//		};
//
//		//measure the proportion of supercritical combinations
//		double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
//		double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
//		double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
//
//		toc = watch.elapsed(TimeUnit.SECONDS);
//
//		if(D) System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
//
//
//		// convert array from log-likelihood to likelihood
//		double testTotalLikelihood = convertLogLikelihoodArrayToLikelihood(maxVal); // this converts likelihood from log to linear
//		toc = watch.elapsed(TimeUnit.SECONDS);
//
//		if(D) System.out.format("Grid search took %d seconds.\n", (int)toc);
//		if(D) System.out.println("Total likelihood  = " + testTotalLikelihood); //debug
//
//		//	//measure the proportion of supercritical combinations
//		//	double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
//		//	double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
//		//	double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
//		//	
//		//	toc = watch.elapsed(timeUnit.SECONDS);
//		//	
//		//	System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
//		if(D) System.out.format("Mainshock productivity magnitude: %2.2f\n" , getProductivityMag());
//
//		watch.stop();
//		this.epiLikelihood = likelihood;
//	}

	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	public double getProductivityMag(){
		if(fitMSProductivity)
			return magMain + ams_vec[max_ams_index] - a_vec[max_a_index];
		else
			return magMain;
	}

}
