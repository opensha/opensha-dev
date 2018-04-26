package scratch.aftershockStatisticsETAS;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.analysis.UnivariateFunction;
//import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
//import org.opensha.commons.data.function.HistogramFunction;
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
	
	Boolean D = true;	// debug flag
	private double rmax;
	private boolean fitMSProductivity;
	private boolean timeDependentMc;
	private volatile boolean stopRequested;
	private volatile boolean pauseRequested;
	private ETAS_AftershockModel_Generic priorModel;
//	private final CalcProgressBar subProgress;

	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double rmax, double[] amsVec, double amsSigma, double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			double minMag, double maxMag, int maxGenerations, int nSims, boolean fitMSProductivity) {

		this(mainshock, aftershocks, Double.NaN, amsVec, amsSigma, aVec, pVec, cVec, alpha, b, refMag,
				dataStartTimeDays, dataEndTimeDays, forecastMinDays, forecastMaxDays, 
				minMag, maxMag, maxGenerations, nSims, fitMSProductivity, false, null, null);
	}

	public ETAS_AftershockModel_SequenceSpecific(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double rmax, double[] amsVec, double amsSigma, double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			double minMag, double maxMag, int maxGenerations, int nSims, boolean fitMSProductivity, boolean timeDependentMc, ETAS_AftershockModel_Generic priorModel){

		this(mainshock, aftershocks,
				Double.NaN,
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
			 								double rmax,
			 								double[] amsVec, double amsSigma, double[] aVec, double[] pVec, double[] cVec, double alpha, double b, double refMag, 	
			 								double dataStartTimeDays, double dataEndTimeDays, double forecastMinDays, double forecastMaxDays, 
			 								double magComplete, double maxMag, int maxGenerations, int nSims, boolean fitMSProductivity, boolean timeDependentMc,
			 								ETAS_AftershockModel_Generic priorModel,
			 								CalcProgressBar progress) {
		
//		// check range values
//		if(num_ams == 1 && min_ams != max_ams) {
//			throw new RuntimeException("Problem: num_ams == 1 && min_ams != max_ams");
//		}
//		if(num_a == 1 && min_a != max_a) {
//			throw new RuntimeException("Problem: num_a == 1 && min_a != max_a");
//		}
//		if(num_p == 1 && min_p != max_p) {
//			throw new RuntimeException("Problem: num_p == 1 && min_p != max_p");
//		}
//		if(num_c == 1 && min_c != max_c) {
//			throw new RuntimeException("Problem: num_c == 1 && min_c != max_c");
//		}
//		if(min_a > max_a) {
//			throw new RuntimeException("Problem: min_a > max_a");
//		}
//		if(min_p > max_p) {
//			throw new RuntimeException("Problem: min_p > max_p");
//		}
//		if(min_c > max_c) {
//			throw new RuntimeException("Problem: min_c > max_c");
//		}

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
//		this.ams_vec = amsVec;
		this.p_vec = pVec;
		this.c_vec = cVec;
		
		this.b = b;
		this.alpha = alpha;
		this.refMag = refMag;
		this.magComplete = magComplete;
		this.rmax = rmax;
//		this.magComplete = magCat;
		this.aftershockList=aftershockList;
		this.mainShock=mainShock;
		this.magMain = mainShock.getMag();
		this.dataStartTimeDays=dataStartTimeDays;
		this.dataEndTimeDays=dataEndTimeDays;
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
//		this.capG=capG;
//		this.capH=capH;
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
			System.out.println("Mc(t):\t" + timeDependentMc + "\trmax:\t" + rmax + " Mc = " + magComplete);
		}
		
//		if(Double.isNaN(capG))
		
		// Find the max like parameters by grid search
		if(D) System.out.println("finding maximum likelihood parameters...");
		
//		// set up the new progress window
//		subProgress = new CalcProgressBar(progress.getOwner(), "Progress", "Calculating maximum likelihood parameters...", false);
////		SwingUtilities.invokeLater(new Runnable() {
////
////			@Override
////			public void run() {
//				subProgress.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
////				subProgress.setModalityType(ModalityType.APPLICATION_MODAL);
//				subProgress.updateProgress(0,1);
//				subProgress.pack();
//				WindowAdapter wl = new WindowAdapter() {
//					@Override
//					public void windowClosing(WindowEvent e){
//						System.out.println("Calculation interrupted");
//						stopRequested = true;
//					}
//				};
//				subProgress.addWindowListener(wl);
//				subProgress.setVisible(true);
////			};
////		});

//		pauseRequested = true;
//		SwingUtilities.invokeLater(new Runnable() {
//
//				@Override
//				public void run() {
//					System.out.println("Doing it");
					
					if (timeDependentMc && (dataEndTimeDays > dataStartTimeDays))
						getLikelihoodMatrixGridFast();
					else 
						getLikelihoodMatrixGridFast();
					
//					subProgress.setVisible(false);
//					pauseRequested = false;
//				};
//			});
		
//		while (pauseRequested) {
//			try {
//				Thread.sleep(1000);
//				System.out.println("zzz...");
//			} catch (InterruptedException e1) {
//				// TODO Auto-generated catch block
//				System.err.println("Calculation was interrupted");
//			}
//		}
		
		
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
		
		if (D)	System.out.println("numAS: " + relativeEventTimes.length + " "+ dataEndTimeDays +" "+ magComplete);

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
//
//	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
//		if(D){
//			System.out.println("Computing Forecast with " + nSims + " simulations (Data window: " + dataMinDays +" "+ dataMaxDays + ", Forecast window: "+ forecastMinDays +" "+ forecastMaxDays + ")");
//			System.out.println("SeqSpec Params: "+ getMaxLikelihood_ams() +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete);
//		}
//		
//		ETAScatalog simulatedCatalog;
//		try{
//			simulatedCatalog = new ETAScatalog(ams_vec, a_vec, p_vec, c_vec, epiLikelihood, alpha, b, refMag, 
//				mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
//		} catch(InterruptedException e) {
//			simulatedCatalog = null;
//		}
//		
//		this.forecastMinDays = forecastMinDays;
//		this.forecastMaxDays = forecastMaxDays;
//		this.simulatedCatalog = simulatedCatalog;
//		this.nSims = nSims;
//	}
//
//	
//		
	/**
	 * Get likelihood matrix with no time dependent Mc. Checks for supercriticality 
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
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++)
				for(int aIndex=0;aIndex<num_a;aIndex++)
					for(int pIndex=0;pIndex<num_p;pIndex++)
						for(int cIndex=0;cIndex<num_c;cIndex++){
							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = 1;
							//			double amsLike = Math.exp(-Math.pow(ams_vec[amsIndex]- mean_a,2)/2/Math.pow(sigma_a, 2));
							//			priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = amsLike;
						}
		}
		
		double maxVal= Double.NEGATIVE_INFINITY;
		double logLike;
		boolean subCritFlag;
		
		double ams, a, k, kms, p, c;
		int amsIndex, aIndex, pIndex, cIndex;
		double timeIntegralMS, NtotAS, NtotMS, Ntot;
		double timeIntegralMSc, NtotASc, NtotMSc, Ntotc;
		
		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
		double productivityMS;
		double[] productivityAS = new double[Nas];
		double[] lambda = new double[Nas];
		double[] timeDecayMS = new double[Nas];
		double[] timeDecayAS = new double[Nas];
		double timeIntegral, timeIntegralc;
		
		
		//do productivities
		productivityMS = Math.pow(10, alpha*(mainShock.getMag()-magComplete) );
		for(int i=0; i<Nas; i++) //compute productivity for this aftershock
			productivityAS[i] = Math.pow(10, alpha*(magAftershocks[i] - magComplete));	//productivity of this aftershock


		double timeCompleteMS, tStartIntegration;
		double[] weights = new double[relativeEventTimes.length];
		double weightMS = 1;
		if (timeDependentMc && priorModel != null) {
			timeCompleteMS = Math.pow(Math.pow(10d, priorModel.mean_ams)*productivityMS/(double)rmax, 1d/priorModel.mean_p) - priorModel.mean_c;
			if (timeCompleteMS < 0) timeCompleteMS = 0;
			if (timeCompleteMS > dataEndTimeDays) timeCompleteMS = dataEndTimeDays;
			if(D) System.out.println("timeComplete = " + timeCompleteMS);
			
			double Mc;
			for (int i = 0; i < relativeEventTimes.length; i++){
				Mc = magComplete - 1/priorModel.b*Math.log10(rmax/Math.pow(10, priorModel.mean_ams)/productivityMS) - priorModel.mean_p/priorModel.b*Math.log10(relativeEventTimes[i]+ priorModel.mean_c);
//				if (Mc > magComplete) {
				if (relativeEventTimes[i] < timeCompleteMS){
					weights[i] = Math.pow(10,-priorModel.b*(Mc - magComplete));
//					weights[i] = 0;
				}
				else
					weights[i] = 1;
			}
			Mc = magComplete - 1/priorModel.b*Math.log10(rmax/Math.pow(10, priorModel.mean_ams)/productivityMS) - priorModel.mean_p/priorModel.b*Math.log10(priorModel.mean_c);
			weightMS = Math.pow(10,-priorModel.b*(Mc - magComplete));
//			weightMS = 0;
		} else { 
			timeCompleteMS = 0;
			for (int i = 0; i < relativeEventTimes.length; i++){
				weights[i] = 1;
			}
			weightMS = 1;
		}
	
		// set up timer/time estimator
		double toc, timeEstimate, n;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;

		
		for(pIndex=0;pIndex<num_p;pIndex++) {
			p = p_vec[pIndex];

			for(cIndex=0;cIndex<num_c;cIndex++) {
				c = c_vec[cIndex];

				tStartIntegration = timeCompleteMS;
				
				
				//compute total number at end of fit window for mainshock (unscaled by a)
				if (p == 1){
					timeIntegralMSc = Math.log(tStartIntegration + c) - Math.log(c);
					timeIntegralMS = Math.log(dataEndTimeDays + c) - Math.log(tStartIntegration + c);
				} else {
					timeIntegralMSc = (Math.pow(tStartIntegration + c, 1-p) - Math.pow(c, 1-p)) / (1-p);
					timeIntegralMS = (Math.pow(dataEndTimeDays + c, 1-p) - Math.pow(tStartIntegration + c, 1-p)) / (1-p);
				}
				NtotMS = productivityMS*timeIntegralMS;
				NtotMSc = productivityMS*timeIntegralMSc*weightMS;
				
				//compute instantaneous intensities and total number for aftershocks (unscaled by a)
				NtotAS = 0;
				NtotASc = 0;
				for(int i=0; i<Nas; i++){
					//compute intensity at this moment due to mainshock (unscaled by ams)
					timeDecayMS[i] = productivityMS/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
					if (relativeEventTimes[i] < timeCompleteMS)
						timeDecayMS[i] *= weightMS;
					
					//compute intensity at this moment due to previous aftershocks (unscaled by a)
					timeDecayAS[i] = 0;
					for(int j = 0; j < i; j++){
						if (relativeEventTimes[i] < timeCompleteMS){
							if(relativeEventTimes[j] < relativeEventTimes[i]){
								//							timeDecayAS[i] += productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p);	//from the aftershocks
								timeDecayAS[i] += weights[j]*productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p);	//from the aftershocks
							}
						} else {
							if(relativeEventTimes[j] < relativeEventTimes[i])
								timeDecayAS[i] += productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p);	//from the aftershocks
						}

					}
				
					tStartIntegration = Math.max(timeCompleteMS, relativeEventTimes[i]);
					
					//compute total number at end of fit window due to this aftershock (unscaled by a)
					if(relativeEventTimes[i] < dataEndTimeDays){
						if(p == 1){
							timeIntegralc = Math.log(tStartIntegration - relativeEventTimes[i] + c) - Math.log(c);
							timeIntegral = Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(tStartIntegration - relativeEventTimes[i] + c); 
						} else {
							timeIntegralc = (Math.pow(tStartIntegration - relativeEventTimes[i] + c, 1-p) - Math.pow(c, 1-p)) / (1-p);
							timeIntegral = (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1-p) - Math.pow(tStartIntegration - relativeEventTimes[i] + c, 1-p)) / (1-p);
						}
						
						NtotAS += productivityAS[i]*timeIntegral;	//aftershock Contributions
						NtotASc += productivityAS[i]*timeIntegralc*weights[i];	//aftershock Contributions
					}
				}
				
				// loop over productivities
				for(amsIndex=0;amsIndex<num_ams;amsIndex++) {
					//			double ams;	//assigned in next loop
					ams = ams_vec[amsIndex];
					
					for(aIndex=0;aIndex<num_a;aIndex++) {
						a = a_vec[aIndex];

						logLike = 0;
						logLike = 0;
						
						// now put in the productivity terms and compute likelihood
						kms = Math.pow(10,ams);
						k = Math.pow(10, a);
						
						Ntot = kms*NtotMS + k*NtotAS;
						Ntotc = kms*NtotMSc + k*NtotASc;
						
						logLike += -Ntot;
						logLike += -Ntotc;
						
						for(int i=0; i<Nas; i++){
//							if (relativeEventTimes[i] > timeCompleteMS){ 
								//compute intensity at this moment due to previous earthquakes
								lambda[i] = kms*timeDecayMS[i]; //from the mainshock
								lambda[i] += k*timeDecayAS[i];	//from the aftershocks
								logLike += Math.log(lambda[i]);
//								logLikec += Math.log(lambda[i])*weights[i];
//							}
						}
						
//						logLike = logLike + logLikec;
						
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
								

							warnTime += 3;

							long count = (pIndex)*(num_c*num_ams*num_a) + (cIndex)*(num_ams*num_a) + (amsIndex)*(num_a) + aIndex;
							long total = (num_p*num_c*num_ams*num_a);
							timeEstimate = toc * total/count;
							System.out.format("This might take a while. Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
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
	 * Get likelihood matrix with time dependent Mc. This is *much* slower than without Mc(t)
	 */
//	private void getLikelihoodMatrixGrid(CancelableProgressBar progress) {
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
//		    public int compare(double[] a, double[] b) {
//		        return Double.compare(a[0], b[0]);
//		    }
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
//		if (priorModel != null)
//			priorLikelihood = priorModel.get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, false);
//		else{
//			for(int amsIndex=0;amsIndex<num_ams;amsIndex++)
//				for(int aIndex=0;aIndex<num_a;aIndex++)
//					for(int pIndex=0;pIndex<num_p;pIndex++)
//						for(int cIndex=0;cIndex<num_c;cIndex++){
//							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = 1;
////							double amsLike = Math.exp(-Math.pow(ams_vec[amsIndex]- mean_a,2)/2/Math.pow(sigma_a, 2));
////							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = amsLike;
//						}
//		}
//		
//		double maxVal= Double.NEGATIVE_INFINITY;
//		double logLike;
//		boolean subCritFlag;
//		
//		double ams, a, k, kms, p, c;
////		double kDefault = Math.pow(10, mean_a);
//		int amsIndex, aIndex, pIndex, cIndex;
//		double timeIntegralAS, timeIntegralMS, NtotAS, NtotMS, Ntot;
//		
//		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
//		double productivityMS;
//		double[] productivityAS = new double[Nas];
//		double[] lambda = new double[Nas];
//		double[] timeDecayMS = new double[Nas];
//		double[] timeDecayAS = new double[Nas];
//		double rateAS;
//		double rateMS;
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
//		
//		// simpson integrator
//		SimpsonIntegrator integral = new SimpsonIntegrator();
//		RateIntegrand integrand = new RateIntegrand();
//		
//		if(D) System.out.println("Mc(t) = " + timeDependentMc +" Data Start: " + dataStartTimeDays + " Data End: " + dataEndTimeDays + " Mc = " + magComplete);
//		
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
//						if (timeDependentMc) {
//							// time-dependent integrator
//							double timeCompleteMS = Math.pow(kms*productivityMS/rmax, 1d/p) - c;
//							double tStartIntegration;
//							double Nbase;
//
//							if (timeCompleteMS > 0){
//
//								integrand.setParams(kms*productivityMS, p, c);
//								double t2 = Math.min(timeCompleteMS, dataEndTimeDays);
//								//						
//								Nbase = integral.integrate((int) 1e6, integrand, 0, t2);
//								tStartIntegration = Math.min(timeCompleteMS, dataEndTimeDays);
//
//							} else {
//								tStartIntegration = 0d;
//								Nbase = 0d;
//							}
//
//							//compute total number at end of fit window for mainshock 
//							if (p == 1)
//								timeIntegralMS = Nbase + kms*productivityMS * Math.log(dataEndTimeDays + c) - Math.log(tStartIntegration + c);
//							else
//								timeIntegralMS = Nbase + kms*productivityMS * (Math.pow(dataEndTimeDays + c, 1d-p) - Math.pow(tStartIntegration + c, 1d-p)) / (1d-p);
//
//							NtotMS = timeIntegralMS;
//						} else {
//							//compute total number at end of fit window for mainshock 
//							if (p == 1)
//								timeIntegralMS = Math.log(dataEndTimeDays + c) - Math.log(dataStartTimeDays + c);
//							else
//								timeIntegralMS = (Math.pow(dataEndTimeDays + c, 1d-p) - Math.pow(dataStartTimeDays + c, 1d-p)) / (1d-p);
//
//							NtotMS = kms*productivityMS*timeIntegralMS;
//						}
//
//						//compute instantaneous intensities and total number for aftershocks 
//						NtotAS = 0;
//						for(int i=0; i<Nas; i++){
//							//compute intensity at this moment due to mainshock (scaled by ams)
//							rateMS = kms*productivityMS/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
//
//							if (timeDependentMc){
//								// scale by detection rate
//								timeDecayMS[i] = rateMS * Math.exp(-rateMS/rmax * Math.pow(10, -b*(magAftershocks[i] - refMag)));
//							} else {
//								timeDecayMS[i] = rateMS; //from the mainshock
//							}
//
//							//compute intensity at this moment due to previous aftershocks (unscaled by a)
//							timeDecayAS[i] = 0;
//							for(int j = 0; j < i; j++){
//								if(relativeEventTimes[j] < relativeEventTimes[i]){
//									rateAS = k*productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p); //from the aftershocks
//
//									if (timeDependentMc){
//										timeDecayAS[i] += rateAS * Math.exp(-rateAS/rmax * Math.pow(10, -b*(magAftershocks[i] - refMag)));
//									} else {
//										timeDecayAS[i] += rateAS;
//									}
//
//								}
//							}
//
//							//compute total number at end of fit window due to this aftershock (unscaled by a)
//							if(relativeEventTimes[i] < dataEndTimeDays){
//								if (timeDependentMc) {
//									/*
//									 * insert time-dependent integrator
//									 */
//									double timeCompleteAS = Math.pow(k*productivityAS[i]/rmax, 1d/p) - c;
//									double tStartIntegration;
//									double Nbase;
//
//									if (timeCompleteAS > 0){
//										// simpson integrator
//										integrand.setParams(k*productivityAS[i], p, c);
//
//										double t2 = Math.min(timeCompleteAS, dataEndTimeDays - relativeEventTimes[i]);
//										Nbase = integral.integrate((int) 1e6, integrand, 0, t2);
//
//										tStartIntegration = Math.min(relativeEventTimes[i] + timeCompleteAS, dataEndTimeDays);
//									}else{
//										tStartIntegration = relativeEventTimes[i];
//										Nbase = 0d;
//									}
//
//									//compute total number at end of fit window for mainshock (unscaled by a)
//									if (p == 1)
//										timeIntegralAS = Nbase + k*productivityAS[i] * Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(tStartIntegration - relativeEventTimes[i] + c);
//									else
//										timeIntegralAS = Nbase + k*productivityAS[i] * (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1d-p) - Math.pow(tStartIntegration - relativeEventTimes[i] + c, 1d-p)) / (1d-p);
//
//									NtotAS += timeIntegralAS;	//aftershock Contributions
//
//								} else {
//									if(p == 1)
//										timeIntegralAS = Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(c); 
//									else
//										timeIntegralAS = (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1d-p) - Math.pow(c, 1d-p)) / (1d-p);
//								
//									NtotAS += k*productivityAS[i]*timeIntegralAS;	//aftershock Contributions
//								}
//
//								
//							}
//						}
//
//
//						Ntot = NtotMS + NtotAS;
//
//						logLike += -Ntot;
//						for(int i=0; i<Nas; i++){
//							//compute intensity at this moment due to previous earthquakes
//							lambda[i] = timeDecayMS[i]; //from the mainshock
//							lambda[i] += timeDecayAS[i];	//from the aftershocks
//							logLike += Math.log(lambda[i]);
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
//							warnTime += 3;
//
//							long count = ((amsIndex)*(num_c*num_p*num_a) + (aIndex)*(num_p*num_c) + (pIndex)*(num_c) + cIndex);
//							long total = (num_p*num_c*num_ams*num_a);
//							timeEstimate = toc * total/count;
//							System.out.format("This might take a while. Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
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
//		System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
//		
//		
//		// convert array from log-likelihood to likelihood
//		testTotalLikelihood = convertLogLikelihoodArrayToLikelihood(maxVal); // this converts likelihood from log to linear
//		toc = watch.elapsed(TimeUnit.SECONDS);
//		
//		if(D) System.out.format("Grid search took %d seconds.\n", (int)toc);
//		if(D) System.out.println("Total likelihood  = " + testTotalLikelihood); //debug
//		
////		//measure the proportion of supercritical combinations
////		double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
////		double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
////		double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
////		
////		toc = watch.elapsed(timeUnit.SECONDS);
////		
////		System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
//		System.out.format("Mainshock productivity magnitude: %2.2f\n" , getProductivityMag());
//		
//		watch.stop();
//		this.epiLikelihood = likelihood;
//	}
	
	/**
	 * Get likelihood matrix with time dependent Mc. This is *much* slower than without Mc(t)
	 */
	private void getLikelihoodMatrixGridMc() {
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

		if (priorModel != null && priorModel.likelihood != null){
			priorLikelihood = priorModel.get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, false);
		}else{
			for(int amsIndex=0;amsIndex<num_ams;amsIndex++)
				for(int aIndex=0;aIndex<num_a;aIndex++)
					for(int pIndex=0;pIndex<num_p;pIndex++)
						for(int cIndex=0;cIndex<num_c;cIndex++)
							priorLikelihood[amsIndex][aIndex][pIndex][cIndex] = 1;
		}

		double maxVal= Double.NEGATIVE_INFINITY;
		double logLike;
		boolean subCritFlag;

		double ams, a, k, kms, p, c;
		//	double kDefault = Math.pow(10, mean_a);
		int amsIndex, aIndex, pIndex, cIndex;
		double timeIntegral, Ntot;

		int Nas = relativeEventTimes.length;	//the number of aftershocks, not counting the mainshock
		double productivityMS;
		double[] productivityAS = new double[Nas];
		double[] lambda = new double[Nas];
		double rateAS;


		//do productivities
		productivityMS = Math.pow(10, alpha*(mainShock.getMag()-magComplete) );
		for(int i=0; i<Nas; i++) //compute productivity for this aftershock
			productivityAS[i] = Math.pow(10, alpha*(magAftershocks[i] - magComplete));	//productivity of this aftershock

		// set up timer/time estimator
		double toc, timeEstimate, n;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;

		// simpson integrator
		//	SimpsonIntegrator integral = new SimpsonIntegrator();
		//	RateIntegrand integrand = new RateIntegrand();

		double timeCompleteMS;
		if (timeDependentMc && priorModel != null) {
			// time-dependent integrator
			timeCompleteMS = Math.pow(Math.pow(10d, priorModel.mean_ams)*productivityMS/(double)rmax, 1d/priorModel.mean_p) - priorModel.mean_c;
			
			
			if (timeCompleteMS < 0) timeCompleteMS = 0;
			if (timeCompleteMS > dataEndTimeDays) timeCompleteMS = dataEndTimeDays;
		} else
			timeCompleteMS = 0;

		if(D) System.out.println("Mc(t) = " + timeDependentMc +" Data Start: " + dataStartTimeDays + " Data End: " + dataEndTimeDays + " Mc = " + magComplete);
		if(D) System.out.println("rmax = " + rmax + " tc = " + timeCompleteMS);
		
		double tStartIntegration;
		// loop over productivities
		for(amsIndex=0;amsIndex<num_ams;amsIndex++) {
			ams = ams_vec[amsIndex];

			for(aIndex=0;aIndex<num_a;aIndex++) {
				a = a_vec[aIndex];

				// productivity scaling
				kms = Math.pow(10,ams);
				k = Math.pow(10, a);

				for(pIndex=0;pIndex<num_p;pIndex++) {
					p = p_vec[pIndex];

					for(cIndex=0;cIndex<num_c;cIndex++) {
						c = c_vec[cIndex];

						logLike = 0;

						tStartIntegration = timeCompleteMS;

						//compute total number at end of fit window for mainshock 
						if (p == 1)
							timeIntegral = kms*productivityMS * (Math.log(dataEndTimeDays + c) - Math.log(tStartIntegration + c));
						else
							timeIntegral = kms*productivityMS * (Math.pow(dataEndTimeDays + c, 1d-p) - Math.pow(tStartIntegration + c, 1d-p)) / (1d-p);
						Ntot = timeIntegral;

						//compute instantaneous intensities and total numbers for aftershocks 
						for(int i=0; i<Nas; i++){
							// compute total number for this aftershock 

							// start either at event time or complete time
							tStartIntegration = Math.max(timeCompleteMS, relativeEventTimes[i]);

							if (p == 1)
								timeIntegral = k*productivityAS[i] * (Math.log(dataEndTimeDays - relativeEventTimes[i] + c) - Math.log(tStartIntegration - relativeEventTimes[i] + c));
							else
								timeIntegral = k*productivityAS[i] * (Math.pow(dataEndTimeDays - relativeEventTimes[i] + c, 1d-p) - Math.pow(tStartIntegration  - relativeEventTimes[i] + c, 1d-p)) / (1d-p);

							Ntot += timeIntegral;

							//compute intensity at this moment due to mainshock (scaled by ams)

							if (relativeEventTimes[i] > timeCompleteMS){
								rateAS = kms*productivityMS/Math.pow(relativeEventTimes[i] + c, p); //from the mainshock
								lambda[i] = rateAS;

								//compute intensity at this moment due to previous aftershocks (unscaled by a)
								for(int j = 0; j < i; j++){
									if(relativeEventTimes[j] < relativeEventTimes[i]){
										rateAS = k*productivityAS[j]/Math.pow(relativeEventTimes[i] - relativeEventTimes[j] + c, p); //from the aftershocks
										lambda[i] += rateAS;
									}
								}
							}
						}

						logLike -= Ntot;
						for(int i=0; i<Nas; i++){
							//compute intensity at this moment due to previous earthquakes
							if (relativeEventTimes[i] > timeCompleteMS)
								logLike += Math.log(lambda[i]);
						}

						//add prior regularization
						logLike += Math.log(priorLikelihood[amsIndex][aIndex][pIndex][cIndex]);

						//check for supercritical parameters over the forecast time window
						n = ETAS_StatsCalc.calculateBranchingRatio(a, p, c, alpha, b, forecastMaxDays, magComplete, maxMag);

						subCritFlag = ( n<1 );

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
							warnTime += 3;

							long count = ((amsIndex)*(num_c*num_p*num_a) + (aIndex)*(num_p*num_c) + (pIndex)*(num_c) + cIndex);
							long total = (num_p*num_c*num_ams*num_a);
							timeEstimate = toc * total/count;
							System.out.format("This might take a while. Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
							if (progress != null){
								//									progress.updateProgress(count, total, String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc)/1000)));
								progress.setProgressMessage(String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc))));
								progress.pack();
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

		//measure the proportion of supercritical combinations
		double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
		double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
		double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);

		toc = watch.elapsed(TimeUnit.SECONDS);

		System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);


		// convert array from log-likelihood to likelihood
		double testTotalLikelihood = convertLogLikelihoodArrayToLikelihood(maxVal); // this converts likelihood from log to linear
		toc = watch.elapsed(TimeUnit.SECONDS);

		if(D) System.out.format("Grid search took %d seconds.\n", (int)toc);
		if(D) System.out.println("Total likelihood  = " + testTotalLikelihood); //debug

		//	//measure the proportion of supercritical combinations
		//	double totalSubCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(subCriticalLikelihood, maxVal);
		//	double totalSuperCriticalLikelihood = convertLogLikelihoodArrayToLikelihood_nonNormalized(superCriticalLikelihood, maxVal);
		//	double fractionSubCritical = totalSubCriticalLikelihood/(totalSubCriticalLikelihood + totalSuperCriticalLikelihood);
		//	
		//	toc = watch.elapsed(timeUnit.SECONDS);
		//	
		//	System.out.format("Sequence Specific Model: %3.2f%% subcritical.\n", fractionSubCritical*100);
		System.out.format("Mainshock productivity magnitude: %2.2f\n" , getProductivityMag());

		watch.stop();
		this.epiLikelihood = likelihood;
	}

	// for time-dependent Mc calculations
	private class RateIntegrand implements UnivariateFunction{
		private double prod, c, p;
		
		public double value(double t) {
			double r = prod / Math.pow(t + c, p);
			double y = rmax * (1d - Math.exp(-r/rmax));
			return y;
		}
		
		public void setParams(double prod, double p, double c){
			this.prod = prod;
			this.c = c;
			this.p = p;
		}
	}
	
	
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
