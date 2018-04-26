package scratch.aftershockStatisticsETAS;

import java.util.concurrent.TimeUnit;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;


/**
 * This represents an ETAS aftershock model where the solution is the Bayesian combination of the generic and sequence-specific models.
 * The generic model is treated as if the covariance between a, p, and c were zero. This means that the generic model uncertainties given 
 * in the generic model data file should be based on the marginal distributions of the parameters. 
 * 
 * a-value discretization is the same as for the sequence-specific model
 * 
 * Note also that the Gaussian distribution is renormalized so that values sum to 1.0 over the final range of
 * a-values represented.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests that cover all cases.
 *
 * @author field, vdElst
 *
 */
public class ETAS_AftershockModel_Bayesian extends ETAS_AftershockModel{
	
	Boolean D=true;	// debug flag
	
	/**
	 * This instantiates a Bayesian combination of the two given ETAS models
	 * @param g_model
	 * @param ss_model
	 */
	public ETAS_AftershockModel_Bayesian(ETAS_AftershockModel_SequenceSpecific ss_model, ETAS_AftershockModel_Generic g_model) {
		
		this.magMain = g_model.getMainShockMag();
		this.b = ss_model.get_b();	//using the seq-specific b-value, which is already Bayesian.
		
		// check similarity of these with the second model.
		Preconditions.checkArgument(areModelsEquivalent(g_model, ss_model),
				"Models are not equivalent so Bayesian combination impossible.");
				
		sigma_ams = g_model.sigma_ams;
		sigma_a = g_model.sigma_a;
		sigma_logc = g_model.sigma_logc;
		sigma_p = g_model.sigma_p;
		mean_ams = g_model.mean_ams;
		mean_a = g_model.mean_a;
		mean_p = g_model.mean_p;
		mean_c = g_model.mean_c;
		double mean_logc = Math.log10(mean_c);
		ams_vec = ss_model.ams_vec;
		a_vec = ss_model.a_vec;
		c_vec = ss_model.c_vec;
		p_vec = ss_model.p_vec;
				
//		double g_likelihood[][][][] = g_model.get_likelihoodMatrix(ams_vec, a_vec, p_vec, c_vec);
		double g_likelihood[][][][] = g_model.get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, true);
		double ss_likelihood[][][][] = ss_model.likelihood;
		double like1, like2;
		double[][][][] likelihood = new double[ams_vec.length][a_vec.length][p_vec.length][c_vec.length];
		double cumSum = 0;

		//debug/timer stuff

		double toc, timeEstimate;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		
//		double amsDiff, aDiff, pDiff, cDiff;
//		double ams, a, logc, p;

		if(D) System.out.println("Computing Bayesian update based on sequence specific fit");
//		for(int h = 0; h < ams_vec.length; h++ ){
//			// the sequence-specific solution already has the ams prior accounted for through the regularization in the seq-specific fit
//			amsDiff = 0;
//			
//			for(int i = 0; i < a_vec.length; i++ ){
//				a = ss_model.a_vec[i];
//				aDiff = (-(a - mean_a)*(a - mean_a)/2/sigma_a/sigma_a);
//			
//				for(int j = 0; j < p_vec.length; j++ ){
//					p = ss_model.p_vec[j];
//					pDiff = (-(p - mean_p)*(p - mean_p)/2/sigma_p/sigma_p);
//				
//					for(int k = 0; k < c_vec.length; k++){
//						logc = Math.log10(c_vec[k]);
//						cDiff = (-(logc - mean_logc)*(logc - mean_logc)/2/sigma_logc/sigma_logc);
//						
//						like1 = Math.exp(amsDiff + aDiff + pDiff + cDiff);
////						like2 = ss_likelihood[h][i][j][k];
//						like2 = 1;
//						likelihood[h][i][j][k] = like1*like2;
//
//						cumSum += likelihood[h][i][j][k];
//						
//						// run the timer to see how long this is going to take
//						if(D){
//							toc = watch.elapsed(TimeUnit.SECONDS);
//							if(toc > warnTime){
//	 							warnTime += 3;
//	 							
//								timeEstimate = toc * (p_vec.length*c_vec.length*ams_vec.length*a_vec.length)/((h)*(c_vec.length*ams_vec.length*a_vec.length) + (i)*(ams_vec.length*a_vec.length) + (j)*(a_vec.length) + k);
//								System.out.format("This might take a while. Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
//							}
//						}
//					}
//				}
//			}
//		}

		for(int h = 0; h < ams_vec.length; h++ ){
			for(int i = 0; i < a_vec.length; i++ ){
				for(int j = 0; j < p_vec.length; j++ ){
					for(int k = 0; k < c_vec.length; k++){
						like1 = g_likelihood[h][i][j][k];
						like2 = ss_likelihood[h][i][j][k];
//						like1 = 1;
						
						likelihood[h][i][j][k] = like1*like2;

						cumSum += likelihood[h][i][j][k];
						
						// run the timer to see how long this is going to take
						if(D){
							toc = watch.elapsed(TimeUnit.SECONDS);
							if(toc > warnTime){
	 							warnTime += 3;
	 							
								timeEstimate = toc * (p_vec.length*c_vec.length*ams_vec.length*a_vec.length)/((h)*(c_vec.length*ams_vec.length*a_vec.length) + (i)*(ams_vec.length*a_vec.length) + (j)*(a_vec.length) + k);
								System.out.format("This might take a while. Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
							}
						}
					}
				}
			}
		}

		if(D) {
			toc = watch.elapsed(TimeUnit.SECONDS);
			System.out.format("It took %d seconds to compute the Bayesian Likelihood matrix\n", (int) toc);
		}
		watch.stop();
		
		//normalize likelihood and find maximum
		double maxVal = Double.NEGATIVE_INFINITY;
		for(int h = 0; h < ams_vec.length; h++ ){
			for(int i = 0; i < a_vec.length; i++ ){
				for(int j = 0; j < p_vec.length; j++ ){
					for(int k = 0; k < c_vec.length; k++){
						likelihood[h][i][j][k] /= cumSum;

						if(maxVal<likelihood[h][i][j][k]){
							maxVal=likelihood[h][i][j][k];
							this.max_ams_index=h;
							this.max_a_index=i;
							this.max_p_index=j;
							this.max_c_index=k;
						}

					}
				}
			}
		}
		
		this.likelihood = likelihood;
		this.mainShock = ss_model.mainShock;
		this.magMain = ss_model.magMain;
		this.aftershockList = ss_model.aftershockList;
		this.dataEndTimeDays = ss_model.dataEndTimeDays;
		this.dataStartTimeDays = ss_model.dataStartTimeDays;
		this.forecastMinDays = ss_model.forecastMinDays;
		this.forecastMaxDays = ss_model.forecastMaxDays;
		this.maxMag = ss_model.maxMag;
		this.magComplete = ss_model.magComplete;
		this.maxGenerations = ss_model.maxGenerations;
		this.nSims = ss_model.nSims;
		this.ams_vec = ss_model.ams_vec;
		this.a_vec = ss_model.a_vec;
		this.p_vec = ss_model.p_vec;
		this.c_vec = ss_model.c_vec;
		this.alpha = ss_model.alpha;
		this.refMag = ss_model.refMag;
		this.min_ams = ams_vec[0];
		this.max_ams = ams_vec[ams_vec.length-1];
		this.num_ams = ams_vec.length;
		this.min_a = a_vec[0];
		this.max_a = a_vec[a_vec.length-1];
		this.num_a = a_vec.length;
		this.min_p = p_vec[0];
		this.max_p = p_vec[p_vec.length-1];
		this.num_p = p_vec.length;
		this.min_c = c_vec[0];
		this.max_c = c_vec[c_vec.length-1];
		this.num_c = c_vec.length;
		this.magAftershocks = ss_model.magAftershocks;
		this.relativeTimeAftershocks = ss_model.relativeTimeAftershocks;
		
		computeNewForecast(dataStartTimeDays,dataEndTimeDays, forecastMinDays,  forecastMaxDays, nSims);
		
	}		

	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
		if(D) System.out.println("Computing Forecast with " + nSims + " simulations (Data window: " + dataMinDays +" "+ dataMaxDays + ", Forecast window: "+ forecastMinDays +" "+ forecastMaxDays + ")");
		if(D) System.out.println("Bayesian params: "+ getMaxLikelihood_ams() +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete);
		
		ETAScatalog simulatedCatalog;
		try{
			simulatedCatalog = new ETAScatalog(ams_vec, a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
				mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
		}catch(InterruptedException e){
			System.out.println("Parameters are supercritical, failed to compute forecast.");
			simulatedCatalog = null;
		}
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.simulatedCatalog = simulatedCatalog;
		this.nSims = nSims;
	}

	
	private boolean areModelsEquivalent(ETAS_AftershockModel g_model, ETAS_AftershockModel ss_model) {
		return true; //might want to revisit this sometime
	}
	
}
