package scratch.aftershockStatisticsETAS;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;

import com.google.common.base.Stopwatch;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;


public class OgataMagFreqDist{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private boolean D = false; //debug
	
	private double bPrior_mean;
	private double bPrior_sigma;
	private double b_Ogata;
	private double mu_Ogata;
	private double sigma_Ogata;
	private double[] magnitudes;
	private double[][][] logLikelihood;
	private int npts = 101;
	double[] b_span;
	double[] mu_span;
	double[] sig_span;

	public OgataMagFreqDist(ObsEqkRupList earthquakes, double bPrior_mean, double bPrior_sigma){
		this.bPrior_mean = bPrior_mean;
		this.bPrior_sigma = bPrior_sigma;

		this.magnitudes = ETAS_StatsCalc.getAftershockMags(earthquakes);
		
		calculateParametersFast();
	}

	public double calculateBWithFixedMc(double mc){
		// this one caluclates b assuming Mc = whatever
		double logDistributionLike;
		double logLikeIncrement;
		
		// do the aki method
		double thisLike;
		double maxLike = Double.NEGATIVE_INFINITY;
		int b_index = 0; 
		
		double logBayesLike;
		double B3 = Math.log(1d/Math.sqrt(2*Math.PI)/bPrior_sigma);

		for (int j = 0; j < b_span.length; j++){
			double B = b_span[j]*Math.log(10);
			double logB = Math.log(B); 

			logDistributionLike = 0;
			for(double mag : magnitudes){
				if (mag >= mc){
					logLikeIncrement = logB - B*(mag - mc);
					logDistributionLike += logLikeIncrement;
				}
			}
			logBayesLike = B3 - (b_span[j] - bPrior_mean)*(b_span[j] - bPrior_mean)/2/bPrior_sigma/bPrior_sigma; //Gaussian prior on b-value
			
			thisLike = logDistributionLike + logBayesLike;

			if (thisLike > maxLike){
				maxLike = thisLike;
				b_index = j;
			}
		}
		
		this.b_Ogata = b_span[b_index];
		this.mu_Ogata = mc;
		this.sigma_Ogata = 0;
		
		return b_span[b_index];
	}


	public double calculateMcWithFixedB(double b){
		return calculateMcWithFixedB(magnitudes, b);
	}
	
	public double calculateMcWithFixedB(double[] magnitudes, double b){
		// this one caluclates Mc assuming b = whatever
		
		double logDistributionLike;
		double logLikeIncrement;
		double logBayesLike;
		double B3 = Math.log(1d/Math.sqrt(2*Math.PI)/bPrior_sigma);

		double[] mc_span = ETAS_StatsCalc.linspace(Math.min(getMin(magnitudes),1.0), getMax(magnitudes), npts);
		
		
		int mc_index = 0;
		for (; mc_index < mc_span.length; mc_index++){
			double mc = mc_span[mc_index];
			
//			b = calculateBWithFixedMc(mc);
			
			double[] b_span = new double[]{b - bPrior_sigma, b, b + bPrior_sigma};
//			if(D) System.out.println("b: " + b_span[0] + " " + b_span[1] + " " + b_span[2]);
			
			double[] thisLike = new double[b_span.length];
//			double thisLike;
//			double N;
//			double prevLike;
//			double bprev = b;
//			int Nprev = 0;
//			
			for (int j = 0; j < b_span.length; j++){
				double B = b_span[j]*Math.log(10);
				double logB = Math.log(B); 
				
				logDistributionLike = 0;

//				N = 0;
				for(double mag : magnitudes){
					if (mag >= mc){
						logLikeIncrement = logB - B*(mag - mc);
						logDistributionLike += logLikeIncrement;
//						N++;
					}
				}
				
				logBayesLike = B3 - (b_span[j] - bPrior_mean)*(b_span[j] - bPrior_mean)/2/bPrior_sigma/bPrior_sigma; //Gaussian prior on b-value
				
				thisLike[j] = logDistributionLike + logBayesLike;
				
			}
//			if(D) System.out.println(mc + ": " + thisLike[0] + " "+ thisLike[1] + " " + thisLike[2]);
			
			if (thisLike[1] > thisLike[0] && thisLike[1] > thisLike[2])
				break;
			else{
//				bprev = b;
//				Nprev = N;
//				prevLike = thisLike;
			}
			
		}
		b_Ogata = b;
		mu_Ogata = mc_span[mc_index]+0.5;
		sigma_Ogata = 0.0; 
		return mc_span[mc_index]+0.5; //sets to highest index for which b not consistent with 1.
	}
		
	private void calculateParametersFast(){
			
		logLikelihood = new double[npts][npts][npts];

		//if the Bayesian prior on b is super narrow, call it fixed and don't waste time looping.
		if (bPrior_sigma > 0.01){
			b_span = ETAS_StatsCalc.linspace(Math.max(bPrior_mean - 3*bPrior_sigma, 0.5), Math.min(bPrior_mean + 5*bPrior_sigma,2), npts);
		} else {
			bPrior_sigma = 0.01;
			b_span = new double[]{bPrior_mean};
		}

		sig_span = ETAS_StatsCalc.linspace(0.01, 1, npts);
		mu_span = ETAS_StatsCalc.linspace(Math.min(getMin(magnitudes),1.0), getMax(magnitudes), npts);
	
		double thisLike;
		double maxLike = Double.NEGATIVE_INFINITY;

		int bIndex = -1;
		int muIndex = -1;
		int sigIndex = -1;

		/////
		double b;
		double mu;
		double sig;
		
		double logLikeIncrement;
		double logBayesLike;
		double logDetectionLike;
		double logDistributionLike;

		double B;
		double B2;
		double logB; 

		double B3 = Math.log(1d/Math.sqrt(2*Math.PI)/bPrior_sigma);

		// set up timer/time estimator
		double toc, timeEstimate;
		int warnTime = 3; //ms
		Stopwatch watch = Stopwatch.createStarted();
		String initialMessageString = "Calculating MFD model. ";
		
		// 3D grid search for maximum likelihood b, mu, sig
		for(int j = 0; j < mu_span.length; j++){
			mu = mu_span[j];
				
			for(int k = 0; k < sig_span.length; k++){
				sig = sig_span[k];
				
				NormalDistribution detectFunction = new NormalDistribution(mu, sig);

				logDetectionLike = 0;
				for(double mag : magnitudes){
						logLikeIncrement = Math.log(detectFunction.cumulativeProbability(mag)); //GR with detection likelihood
						logDetectionLike += logLikeIncrement;
				}

				for(int i = 0; i < b_span.length; i++){
					b = b_span[i];
					B = b*Math.log(10);
					B2 = B*B*sig*sig/2;
					logB = Math.log(B); 

					logDistributionLike = 0;
					for(double mag : magnitudes){
							logLikeIncrement = logB - B*(mag - mu) - B2;
							logDistributionLike += logLikeIncrement;
					}
					
					logBayesLike = B3 - (b - bPrior_mean)*(b - bPrior_mean)/2/bPrior_sigma/bPrior_sigma; //Gaussian prior on b-value
					thisLike = logDetectionLike + logDistributionLike + logBayesLike;
					
					if(thisLike > maxLike){
						maxLike = thisLike;
						bIndex = i;
						muIndex = j;
						sigIndex = k;
					}

					
					logLikelihood[i][j][k] = thisLike;
//					if(D) System.out.println(i + " " + j + " " + k + " " + " " + logLikelihood[i][j][k] + " " + b + " " + mu + " " + sig );
					
					
					// run the timer to see how long this is going to take
					toc = watch.elapsed(TimeUnit.SECONDS);
					if(toc > warnTime){
						// 							longRunFlag = true;	//mark as done. don't do it again
						warnTime += 10;

						timeEstimate = toc * (npts*npts*b_span.length)/((j)*(npts*b_span.length) + (k)*(b_span.length) + i);
						System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
						initialMessageString = "...";
					}
				}
			}
		}
		watch.stop();

		if(D) System.out.println("max like: " +maxLike);
		///scale it all up?
		for(int i = 0; i<npts; i++){
			for(int j = 0; j<npts; j++){
				for(int k = 0; k<npts; k++){
					logLikelihood[i][j][k] += -maxLike;
				}
			}
		}
		
		b_Ogata = b_span[bIndex];
		mu_Ogata = mu_span[muIndex];
		sigma_Ogata = sig_span[sigIndex];
		
	}

	private double getMin(double[] mags){
		double minMag = Double.POSITIVE_INFINITY;
		for (int i = 0; i < mags.length; i++){
			if (mags[i] < minMag)
				minMag = mags[i];
		}
		return minMag;
	}
	
	private double getMax(double[] mags){
		double maxMag = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < mags.length; i++){
			if (mags[i] > maxMag)
				maxMag = mags[i];
		}
		return maxMag;
	}
	
	double computeLogLikelihood(double b, double mu, double sig){

		double mag;
		double logLikeIncrement;
		double bayesLikeIncrement;
		double logLikelihood = 0;

		double B = b*Math.log(10);
		double B2 = B*B*sig*sig/2;
		double logB = Math.log(B); 

		double B3 = Math.log(1d/Math.sqrt(2*Math.PI)/bPrior_sigma);

		NormalDistribution detectFunction = new NormalDistribution(mu, sig);

		for(int i = 0; i < magnitudes.length; i++){
			mag = magnitudes[i];
			logLikeIncrement = logB - B*(mag - mu) - B2 + Math.log(detectFunction.cumulativeProbability(mag)); //GR with detection likelihood

			logLikelihood += logLikeIncrement;
		}
		bayesLikeIncrement = B3 - (b - bPrior_mean)*(b - bPrior_mean)/2/bPrior_sigma/bPrior_sigma; //Gaussian prior on b-value

		return logLikelihood + bayesLikeIncrement;
	};

	/*
	 * Returns magComplete as 2-sigma above 50% detection magnitude
	 */
	public double get_Mc(){
		return get_Mc(2d);
	}

	/*
	 * Returns magComplete as nstds-sigma above 50% detection magnitude
	 */
	public double get_Mc(double nstds){
		return get_Mc(mu_Ogata, sigma_Ogata, nstds);
	}

	public double get_Mc(double mu, double sig, double nstds){
		return mu + nstds*sig;
	}

	public double get_bValue(){
		return b_Ogata;
	}

	@Override
	public String toString() {
		return "b-Value: " + get_bValue() + " mu: " + mu_Ogata + " sigma: " + sigma_Ogata + " Mc: " + get_Mc();
	}

	public double getDensity(double mag){
		double B = b_Ogata*Math.log(10);
		double B2 = B*B*sigma_Ogata*sigma_Ogata/2;
		double logB = Math.log(B); 
		
		if (sigma_Ogata>0){
			NormalDistribution detectFunction = new NormalDistribution(mu_Ogata, sigma_Ogata);
			return Math.exp(logB - B*(mag - mu_Ogata) - B2) * detectFunction.cumulativeProbability(mag);
		} else {
			if (mag >= mu_Ogata)
				return Math.exp(logB - B*(mag - mu_Ogata) - B2);
			else
				return 0d;
		}
	}

	private class RateIntegrand implements UnivariateFunction{
		public double value(double mag) {
			double y = getDensity(mag);
			return y;
		}
	}

	public double getRate(double mag, double magIncrement){
		return getDensity(mag) * magnitudes.length * magIncrement;
	}

	public EvenlyDiscretizedFunc getMND(double minMag, double maxMag, double magIncrement){
		EvenlyDiscretizedFunc mnd = new EvenlyDiscretizedFunc(minMag, maxMag, (int) (1 + (maxMag - minMag)/magIncrement));

		for (int i = 0; i < mnd.size(); i++) {
			double val = getRate(mnd.getX(i), magIncrement);
			mnd.set(i, val);
		}

		return mnd;
	}

	public EvenlyDiscretizedFunc getCumulativeMND(double minMag, double maxMag, double magIncrement){
		EvenlyDiscretizedFunc cumMND = new EvenlyDiscretizedFunc(minMag, maxMag, (int) (1.5 + (maxMag - minMag)/magIncrement));
		
		if (sigma_Ogata == 0){
			double mc = get_Mc();
			int count = 0;
			for (double mag : magnitudes){
				if (mag >= mc) count++;
			}
			
			for (int i = 0; i < cumMND.size(); i++) {
				double val = Math.pow(10,-b_Ogata*(cumMND.getX(i) - magIncrement/2 - get_Mc()));
				cumMND.set(i, val*count);
			}
		} else {
			SimpsonIntegrator cumSum = new SimpsonIntegrator();
			UnivariateFunction integrand = new RateIntegrand();

			for (int i = 0; i < cumMND.size(); i++) {
				double val = cumSum.integrate((int) 1e6, integrand, cumMND.getX(i) - magIncrement/2, maxMag);
				cumMND.set(i, val*magnitudes.length);
			} 
		}
		return cumMND;
	}

//	public EvenlyDiscretizedFunc getMcCumulativeDensity(){
//		double Mc;
//		double like = 0;
//		
//		HistogramFunction mcfuncEven = new HistogramFunction(mu_span[0], mu_span[mu_span.length-1], mu_span.length);
//	
//		double likeSum = 0;
//		if (logLikelihood != null){
//			for (int imu = 0; imu <mu_span.length; imu++){
//				for (int isig = 0; isig <sig_span.length; isig++){
//					Mc = mu_span[imu] + 2d*sig_span[isig];
//					
//					if(Mc <= mu_span[mu_span.length-1]){
//						like = 0;
//						for (int ib = 0; ib<b_span.length; ib++){
//							if(logLikelihood[ib][imu][isig] > -50){
//								like += Math.exp(logLikelihood[ib][imu][isig]);
//							}
//						}
//						likeSum += like;
//						mcfuncEven.add(Mc, like);
//
//					}
//				}
//			}
//			mcfuncEven.scale(1.0/likeSum);
//			if(D) System.out.println("likelihood sum = " + likeSum);
//			
//			// make cumulative;
//			double runningTot = 0;
//			for (int i = 0; i < mcfuncEven.size(); i++){
//				runningTot += mcfuncEven.getY(i);
//				mcfuncEven.set(i, runningTot);
//			}
//			if(D) System.out.println(runningTot);
//			return mcfuncEven;
//		} else {
//			return null;
//		}
//	}
	
	public double goodnessOfFit(){
		double gof;
		
		Arrays.sort(magnitudes);
		
		// get numbers
		int index = 0;
		double[] b = new double[mu_span.length];
		int muIndex = 0;
		for (double mu:mu_span){
			while (magnitudes[index++] < mu);
			index--;

			b[muIndex] = magnitudes.length-index;
			muIndex++;
		}
		
		if(D) System.out.println(this);
		EvenlyDiscretizedFunc Nexp = getCumulativeMND(mu_span[0], mu_span[mu_span.length-1], mu_span[1]-mu_span[0]);
		if(D) System.out.println(Nexp);
		
		double gofNum = 0;
		double gofDen = 1;
		for (int i = 0; i < b.length; i++){
			if (mu_span[i] >= get_Mc()){
				
				if(D) System.out.println(mu_span[i] + ": " + b[i] + " " + Nexp.getY(i));				
				
				gofNum += Math.abs(b[i] - Nexp.getY(i));
				gofDen += b[i];
			}
		}
		
		gof = 1 - gofNum/gofDen;
		if(D) System.out.println(gof);
		return gof;
		
	}
	public double verify(double mc){
//		EvenlyDiscretizedFunc cumDensityMc = getMcCumulativeDensity();
//		if (mc < cumDensityMc.getMinX()){
//			System.out.println("Warning: setting Mc to a lower value than supported by the data \n"
//					+ " will cause aftershock activity to be underestimated.");
//		}
//			try {
//				if (mc < cumDensityMc.getMinX()){
//					System.err.println("Warning: Mc is outside model 95% confidence bound: " + String.format("%3.2f", cumDensityMc.getFirstInterpolatedX(0.025)));
//					mc = cumDensityMc.getFirstInterpolatedX(0.025);
//				} else if (mc > cumDensityMc.getMaxX()) {
//					// do nothing
//				} else if (cumDensityMc.getInterpolatedY(mc) < 0.025){
//					System.err.println("Warning: Mc is outside model 95% confidence bound: " + String.format("%3.2f", cumDensityMc.getFirstInterpolatedX(0.025)));
//					mc = cumDensityMc.getFirstInterpolatedX(0.025);
//				}
//			} catch (Exception e){
//				System.err.println("Warning: could not validate choice of mc parameter.");
//			}
		return mc;
	}
}
