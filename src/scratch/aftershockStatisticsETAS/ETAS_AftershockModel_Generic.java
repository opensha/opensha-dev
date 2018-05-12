package scratch.aftershockStatisticsETAS;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

//import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.gui.infoTools.CalcProgressBar;


/**
 * This represents an ETAS aftershock model where a-values are assumed to 
 * be Gaussian distributed and p and c are held fixed at the values given.
 * 
 * a-value discretization is hard-coded as 0.01 (the resolution of value given for various regions by
 * Page et al. (2016).
 * 
 * Note also that the Gaussian distribution is renormalized so that values sum to 1.0 over the range of
 * a-values represented.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests that cover all cases.
 *
 * @author field
 * @author van der Elst
 *
 */
public class ETAS_AftershockModel_Generic extends ETAS_AftershockModel {
	
	private static Boolean D=false;	// debug flag
	
	private double[][] covariance;
	private double[][] covInverse;
	private double covDeterminant;
	private double[][] priorCovariance;
	private double[][] priorCovInverse;
	private double priorCovDeterminant;
		
//	public double[][][][] priorLikelihood; //is this used?
	
	public ETAS_AftershockModel_Generic(GenericETAS_Parameters genericParams){
		this.mean_a = genericParams.get_a();
		this.mean_ams = genericParams.get_ams();
		this.mean_c = genericParams.get_c();
		this.mean_p = genericParams.get_p();
		this.alpha = genericParams.get_alpha();
		this.b = genericParams.get_b();
		this.ac = genericParams.get_a();
	};
	
	public ETAS_AftershockModel_Generic(
			ObsEqkRupture mainShock, ObsEqkRupList aftershockList, GenericETAS_Parameters genericETAS_Parameters,
			double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, double Mc, 
			double maxSimMag, int maxNumGenerations, int nSims, boolean fitMSProductivity, boolean timeDependentMc) {

		this(mainShock,  aftershockList,  genericETAS_Parameters,
				dataMinDays,  dataMaxDays,  forecastMinDays,  forecastMaxDays,  Mc, 
				maxSimMag,  maxNumGenerations,  nSims,  fitMSProductivity, timeDependentMc, null); 
	}

	
	/**
	 * This instantiates a generic ETAS model from a GenericETAS_Parameters object, where aValueMin and aValueMax
	 * are set as -4.5 and -0.5, respectively.
	 * @param mainShock
	 * @param aftershockList
	 * @param genericETAS_Parameters
	 */
	public ETAS_AftershockModel_Generic(
			ObsEqkRupture mainShock, ObsEqkRupList aftershockList, GenericETAS_Parameters genericETAS_Parameters,
			double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, double Mc, 
			double maxSimMag, int maxNumGenerations, int nSims, boolean fitMSProductivity, boolean timeDependentMc, CalcProgressBar progress) {
		
		this(mainShock, aftershockList,
				genericETAS_Parameters.get_a(), genericETAS_Parameters.get_aSigma(), 
				genericETAS_Parameters.get_p(), genericETAS_Parameters.get_pSigma(),
				genericETAS_Parameters.get_c(), genericETAS_Parameters.get_logcSigma(), 
				genericETAS_Parameters.get_covariance(), genericETAS_Parameters.get_priorCovariance(),
				genericETAS_Parameters.get_alpha(), genericETAS_Parameters.get_b(), genericETAS_Parameters.get_refMag(),
				dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, Mc,
				maxSimMag, maxNumGenerations, nSims, fitMSProductivity, timeDependentMc, progress); 
	}
	
	
	/**
	 * This instantiates a generic ETAS model for the values given
	 * @param magMain - main shock magnitude
	 * @param mean_a - mean a-value for the Gaussian distribution
	 * @param sigma_a - a-value standard deviation for the Gaussian distribution
	 * @param b - b-value
	 * @param p - p-value
	 * @param c - c-value
	 */
	public ETAS_AftershockModel_Generic(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double mean_a, double sigma_a, 
			double mean_p, double sigma_p,
			double mean_c, double sigma_logc,
			double[][] covariance, double[][] priorCovariance,
			double alpha, double b, double refMag,
			double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, 
			double Mc, double maxMag,
			int maxGenerations, int nSims,
			boolean fitMSProductivity,
			boolean timeDependentMc,
			CalcProgressBar progress) {

		//		this.simulatedCatalog = simulatedCatalog;
		
		if(D) System.out.println("Initializing Generic model...");

		double bSigma = 0.1;						//add this to generic param set

		this.progress = progress;
		this.mainShock=mainshock;
		this.aftershockList=aftershocks;
		this.magMain = mainshock.getMag();
		this.magAftershocks = ETAS_StatsCalc.getAftershockMags(aftershocks);
		this.b=b;
		this.bSigma = bSigma;
		this.timeDependentMc = timeDependentMc;

		//correct the a-value if Mc is not the same as refMag
		mean_a += Math.log10((maxMag - refMag)/(maxMag - Mc));

		this.ac = mean_a;
		this.mean_ams=mean_a;
		this.sigma_ams=sigma_a;
		this.mean_a=mean_a;
		this.sigma_a=sigma_a; //use this for prior and for generating mainshock variability for forecast
		this.mean_p = mean_p;
		this.sigma_p = sigma_p;
		this.mean_c = mean_c;
		this.sigma_logc = sigma_logc;
		this.covariance = covariance; //use this for generating forecasts
		this.priorCovariance = priorCovariance; //use this for generating forecasts
		this.alpha = alpha;
		this.refMag = refMag;
		this.magComplete = Mc;
		this.nSims = nSims;
		this.maxMag = maxMag;
		this.maxGenerations = maxGenerations;
		this.dataStartTimeDays = dataMinDays;
		this.dataEndTimeDays = dataMaxDays;
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
				
		// set up solution vectors
		double[] ams_vec;
		if(fitMSProductivity && num_ams > 1){
			min_ams = mean_ams - 3*sigma_ams;
			max_ams = mean_ams + 3*sigma_ams;
			ams_vec = ETAS_StatsCalc.linspace(min_ams, max_ams, num_ams);
		}else{
			min_ams = mean_ams;
			max_ams = mean_ams;
			num_ams = 1;
			ams_vec = new double[]{mean_a};
		}
//		System.out.println("ams_vec: " + min_ams +" "+ max_ams +" "+ num_ams);
		
		double[] a_vec;
		if(sigma_a == 0 || num_a == 1){
			min_a = mean_a;
			max_a = mean_a;
			num_a = 1;
			a_vec = new double[]{mean_a};
		}else{
			min_a = mean_a - 3*sigma_a;
			max_a = mean_a + 3*sigma_a;
			a_vec = ETAS_StatsCalc.linspace(min_a, max_a, num_a);
		}
//		System.out.println("a_vec: " + min_a +" "+ max_a +" "+ num_a);
		
		double[] p_vec;
		if(sigma_p == 0 || num_p == 1){
			min_p = mean_p;
			max_p = mean_p;
			num_p = 1;
			p_vec = new double[]{mean_p};
		}else{
			min_p = mean_p - 3*sigma_p;
			max_p = mean_p + 3*sigma_p;
			p_vec = ETAS_StatsCalc.linspace(min_p, max_p, num_p);
		}
//		System.out.println("p_vec: " + min_p +" "+ max_p +" "+ num_p);

		double[] c_vec;
		if(sigma_logc == 0 || num_c == 1){
			c_vec = new double[]{mean_c};
			min_c = mean_c;
			max_c = mean_c;
			num_c = 1;
		}else{
			min_c = Math.pow(10, Math.log10(mean_c) - 3*sigma_logc);
			max_c = Math.pow(10, Math.log10(mean_c) + 3*sigma_logc);
			c_vec = ETAS_StatsCalc.logspace(min_c, max_c, num_c);
		}
//		System.out.println("c_vec: " + min_c +" "+ max_c +" "+ num_c);

		this.ams_vec = ams_vec;
		this.a_vec = a_vec;
		this.p_vec = p_vec;
		this.c_vec = c_vec;
		
		this.max_ams_index = Math.round((ams_vec.length-1)/2);
		this.max_a_index = Math.round((a_vec.length-1)/2);
		this.max_p_index = Math.round((p_vec.length-1)/2);
		this.max_c_index = Math.round((c_vec.length-1)/2);
		
		if(sigma_a<=0){
			throw new RuntimeException("Problem: sigma_a must be greater than 0");
		}
		
		//set up some matrices for easy calculations later
		//invert covariance matrix
		RealMatrix COV = MatrixUtils.createRealMatrix(covariance);
		RealMatrix iCOV = MatrixUtils.blockInverse(COV, 1);
		double[][] icovariance = new double[3][3];
		for(int i = 0; i < 3 ; i++){
			for(int j = 0; j < 3 ; j++){
				icovariance[i][j] = iCOV.getEntry(i, j);
			}
		}
		this.covInverse = icovariance;
		
		//find determinant of covariance matrix
		double determinant = new CholeskyDecomposition(COV).getDeterminant();
//		System.out.println("Determinant: " + determinant); //debug
		this.covDeterminant = determinant;
		
		//set up some matrices for easy calculations later
		//invert covariance matrix
		RealMatrix pCOV = MatrixUtils.createRealMatrix(priorCovariance);
		RealMatrix ipCOV = MatrixUtils.blockInverse(pCOV, 1);
		double[][] ipCovariance = new double[3][3];
		for(int i = 0; i < 3 ; i++){
			for(int j = 0; j < 3 ; j++){
				ipCovariance[i][j] = ipCOV.getEntry(i, j);
			}
		}
		this.priorCovInverse = ipCovariance;

		//find determinant of covariance matrix
		double pDeterminant = new CholeskyDecomposition(pCOV).getDeterminant();
		//				System.out.println("Determinant: " + determinant); //debug
		this.priorCovDeterminant = pDeterminant;

		//get likelihoods (epistemic and prior("likelihood))
		epiLikelihood = get_likelihoodMatrix(ams_vec, a_vec, p_vec, c_vec);		//the tight covariance
		likelihood = get_priorLikelihoodMatrix(ams_vec, a_vec, p_vec, c_vec, true);		//the prior covariance
		
		// get aftershock times and mags and store as simple doubles[]
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = ETAS_StatsCalc.getAftershockMags(aftershockList.getRupsAboveMag(magComplete));

		List<double[]> sortedEQlist = new ArrayList<double[]>();

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

		if(D) System.out.println(relativeEventTimes.length + " events found " + sortedEQlist.size() + " events in list");
		
		if (relativeEventTimes.length > 0) {
			for(int i = 0; i < relativeEventTimes.length; i++){
				double[] temp = sortedEQlist.get(i);
				relativeEventTimes[i] = temp[0];
				magAftershocks[i] = temp[1];
				//					sortedEQlist.add(temp);
			}
		}
		
		this.magAftershocks = magAftershocks;
		this.relativeTimeAftershocks = relativeEventTimes;
		
		// initialize a forecast
		computeNewForecast(dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, nSims);
		
	}		

//	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
//		if(D) System.out.println("Computing Forecast with " + nSims + " simulations (Data window: " + dataMinDays +" "+ dataMaxDays + ", Forecast window: "+ forecastMinDays +" "+ forecastMaxDays + ")");
//		if(D) System.out.println("Generic Params: "+ mean_ams +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete);
//		
//		ETAScatalog simulatedCatalog;
//		try{
//			simulatedCatalog = new ETAScatalog(ams_vec, a_vec, p_vec, c_vec, epiLikelihood, alpha, b, refMag, 
//				mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
//		}catch(InterruptedException e){
//			simulatedCatalog = null;
//		}
//		this.forecastMinDays = forecastMinDays;
//		this.forecastMaxDays = forecastMaxDays;
//		this.simulatedCatalog = simulatedCatalog;
//		this.nSims = nSims;
//	}
//

	/**
	 * Returns likelihood matrix for vectors [a,p,c] assuming a MULTIVARIATEGaussian distribution on each parameter (specified by covariance matrix.)
	 */
	public double[][][][] get_likelihoodMatrix(double[] ams_vec, double[] a_vec, double[] p_vec, double[] c_vec){
		double[][][][] likelihood = new double[ams_vec.length][a_vec.length][p_vec.length][c_vec.length];
		double cumSum = 0;

		if(D) System.out.println("generating Generic likelihood matrix...");
		if(D) System.out.println("Mc-adjusted params: " + mean_ams + " (" +sigma_ams +") " + mean_a +" ("+sigma_a+") "+ mean_p+" ("+sigma_p+") "+ mean_c + " ("+sigma_logc+") ");
		long tic = System.currentTimeMillis();
		long toc;
		
		for(int h = 0; h < ams_vec.length ; h++ ){
			for(int i = 0; i < a_vec.length ; i++ ){
				for(int j = 0; j < p_vec.length ; j++ ){
					for(int k = 0; k < c_vec.length ; k++ ){
						likelihood[h][i][j][k] = get_likelihood(ams_vec[h], a_vec[i], p_vec[j], c_vec[k]);
						cumSum += likelihood[h][i][j][k];
					}
				}
			}
		}
		
		for(int h = 0; h < ams_vec.length ; h++ ){
			for(int i = 0; i < a_vec.length ; i++ ){
				for(int j = 0; j < p_vec.length ; j++ ){
					for(int k = 0; k < c_vec.length ; k++ ){
						likelihood[h][i][j][k] /= cumSum;
					}
				}
			}
		}
		
		toc = System.currentTimeMillis();
		if (D) System.out.println("It took " + (toc-tic)/1000 + " seconds to compute generic likelihood matrix.");
		
		return likelihood;
	}
		
	/**
	 * Returns likelihood for given [a,p,c] assuming a 3D Gaussian distribution)
	 */
	public double get_likelihood(double ams, double a, double p, double c){
		double like;
		
		double logc = Math.log10(c);
		double mean_logc = Math.log10(mean_c);
		
		double[] X = {a, p, logc};
		double[] MU = {mean_a, mean_p, mean_logc};
		
		double D = 0;
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				D += (X[i] - MU[i]) * (X[j] - MU[j]) * covInverse[j][i];
			}
		}
		like = Math.exp(-0.5*D)/Math.sqrt(Math.pow(2d*Math.PI,3) * covDeterminant ) ;
		
		double amsLike = 1/Math.sqrt(2*Math.PI)/sigma_ams * Math.exp(-(ams - mean_a)*(ams - mean_a)/2/sigma_ams/sigma_ams);
		
		like *= amsLike;
		return like;
		
	}


	/**
	 * Returns likelihood matrix for vectors [a,p,c] assuming a Gaussian distribution on each parameter (specified by mean_a, sigma_a, mean_p, sigma_p, etc.)
	 */
	public double[][][][] get_priorLikelihoodMatrix(double[] ams_vec, double[] a_vec, double[] p_vec, double[] c_vec, boolean normalized){
		double[][][][] likelihood = new double[ams_vec.length][a_vec.length][p_vec.length][c_vec.length];
		double cumSum = 0;

		if(D )System.out.println("generating Generic prior...");
		long tic = System.currentTimeMillis();
		long toc;

		for(int h = 0; h < ams_vec.length ; h++ ){
			for(int i = 0; i < a_vec.length ; i++ ){
				for(int j = 0; j < p_vec.length ; j++ ){
					for(int k = 0; k < c_vec.length ; k++ ){
						likelihood[h][i][j][k] = get_priorLikelihood(ams_vec[h], a_vec[i], p_vec[j], c_vec[k]);
						cumSum += likelihood[h][i][j][k];
					}
				}
			}
		}

		if (normalized){
			for(int h = 0; h < ams_vec.length ; h++ ){
				for(int i = 0; i < a_vec.length ; i++ ){
					for(int j = 0; j < p_vec.length ; j++ ){
						for(int k = 0; k < c_vec.length ; k++ ){
							likelihood[h][i][j][k] /= cumSum;
						}
					}
				}
			}
		}

		toc = System.currentTimeMillis();
		if(D) System.out.println("It took " + (toc-tic)/1000 + " seconds to compute the Generic prior.");

		return likelihood;
	}
	
	/**
	 * Returns likelihood for given [a,p,c] assuming a Gaussian distribution on each parameter (specified by mean_a, sigma_a, mean_p, sigma_p, etc.)
	 */
	public double get_priorLikelihood(double ams, double a, double p, double c){
//		double like;
//
//		double logc = Math.log10(c);
//		double mean_logc = Math.log10(mean_c);
//
//		double amsLike = 1/Math.sqrt(2*Math.PI)/sigma_ams * Math.exp(-(ams - mean_a)*(ams - mean_a)/2/sigma_ams/sigma_ams);
//		double aLike = 1/Math.sqrt(2*Math.PI)/sigma_a * Math.exp(-(a - mean_a)*(a - mean_a)/2/sigma_a/sigma_a);
//		double pLike = 1/Math.sqrt(2*Math.PI)/sigma_p * Math.exp(-(p - mean_p)*(p - mean_p)/2/sigma_p/sigma_p);
//		double cLike = 1/Math.sqrt(2*Math.PI)/sigma_logc * Math.exp(-(logc - mean_logc)*(logc - mean_logc)/2/sigma_logc/sigma_logc);
//		
//		like = amsLike*aLike*pLike*cLike;
//		
//		return like;
		double like;
		
		double logc = Math.log10(c);
		double mean_logc = Math.log10(mean_c);
		
		double[] X = {a, p, logc};
		double[] MU = {mean_a, mean_p, mean_logc};
		
		double D = 0;
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				D += (X[i] - MU[i]) * (X[j] - MU[j]) * priorCovInverse[j][i];
			}
		}
		like = Math.exp(-0.5*D)/Math.sqrt(Math.pow(2d*Math.PI,3) * priorCovDeterminant ) ;
//		double aLike = ???
		double amsLike = 1/Math.sqrt(2*Math.PI)/sigma_ams * Math.exp(-(ams - mean_a)*(ams - mean_a)/2/sigma_ams/sigma_ams);
		
		like *= amsLike;
		return like;
		
	}

	public HistogramFunction getPriorPDF_ams() {
		// this overrides the grid search vector to return a full prior distribution regardless of the constraints being used
		
		double min_ams = mean_a - 3*sigma_ams;
		double max_ams = mean_a + 3* sigma_ams;
		int num_ams = 101;
		double[] likelihood = new double[num_ams];
		double likelihoodTotal = 0;
		
		HistogramFunction hist = new HistogramFunction(min_ams, max_ams, num_ams);
		String name = "ams-value prior distribution";
		
		for(int amsIndex=0;amsIndex<num_ams;amsIndex++) {
			likelihood[amsIndex] = 1/Math.sqrt(2*Math.PI)/sigma_ams * Math.exp(- (mean_a - hist.getX(amsIndex))*(mean_a - hist.getX(amsIndex))/2/sigma_ams/sigma_ams);   
			likelihoodTotal += likelihood[amsIndex];
		}
		
		for (int amsIndex=0;amsIndex<num_ams;amsIndex++) {
			likelihood[amsIndex] /= likelihoodTotal;
			hist.add(amsIndex,  likelihood[amsIndex]);
		}
		
		hist.setName(name);
		if(D) System.out.println("PDF of ams-value: totalTest = "+hist.calcSumOfY_Vals());
		
		hist.scale(1d/hist.getDelta());
		return hist;
	}
	
	public HistogramFunction getPriorPDF_a() {
		// this overrides the grid search vector to return a full prior distribution regardless of the constraints being used
		
		double min_a = mean_a - 3*sigma_a;
		double max_a = mean_a + 3* sigma_a;
		int num_a = 101;
		double[] likelihood = new double[num_a];
		double likelihoodTotal = 0;
		
		HistogramFunction hist = new HistogramFunction(min_a, max_a, num_a);
		String name = "a-value prior distribution";
		
		for(int aIndex=0;aIndex<num_a;aIndex++) {
			likelihood[aIndex] = 1/Math.sqrt(2*Math.PI)/sigma_a * Math.exp(- (mean_a - hist.getX(aIndex))*(mean_a - hist.getX(aIndex))/2/sigma_a/sigma_a);   
			likelihoodTotal += likelihood[aIndex];
		}
		
		for (int aIndex=0;aIndex<num_a;aIndex++) {
			likelihood[aIndex] /= likelihoodTotal;
			hist.add(aIndex,  likelihood[aIndex]);
		}
		
		hist.setName(name);
		if(D) System.out.println("PDF of a-value: totalTest = "+hist.calcSumOfY_Vals());
		
		hist.scale(1d/hist.getDelta());
		return hist;
	}
	
	public HistogramFunction getPriorPDF_p() {
		// this overrides the grid search vector to return a full prior distribution regardless of the constraints being used
		
		double min_p = mean_p - 3*sigma_p;
		double max_p = mean_p + 3* sigma_p;
		int num_p = 101;
		double[] likelihood = new double[num_p];
		double likelihoodTotal = 0;
		
//		System.out.println("p: " + getMaxLikelihood_p() + " " + min_p + " " + max_p + " " + num_p);
		HistogramFunction hist = new HistogramFunction(min_p, max_p, num_p);
		String name = "p-value prior distribution";
		
		for(int pIndex=0;pIndex<num_p;pIndex++) {
			likelihood[pIndex] = 1/Math.sqrt(2*Math.PI)/sigma_p * Math.exp(- (mean_p - hist.getX(pIndex))*(mean_p - hist.getX(pIndex))/2/sigma_p/sigma_p);   
			likelihoodTotal += likelihood[pIndex];
		}
		
		for (int pIndex=0;pIndex<num_p;pIndex++) {
			likelihood[pIndex] /= likelihoodTotal;
			hist.add(pIndex,  likelihood[pIndex]);
		}
		
		hist.setName(name);
		if(D) System.out.println("PDF of p-value: totalTest = "+hist.calcSumOfY_Vals());
		
		hist.scale(1d/hist.getDelta());
		return hist;
	}
	
	public HistogramFunction getPriorPDF_logc() {
		// this overrides the grid search vector to return a full prior distribution regardless of the constraints being used
		
		double mean_logc = Math.log10(mean_c);
		double min_logc = mean_logc - 3*sigma_logc;
		double max_logc = mean_logc + 3*sigma_logc;
		int num_c = 101;
		double[] likelihood = new double[num_c];
		double likelihoodTotal = 0;
		
		if(D) System.out.println("c: " + getMaxLikelihood_c() + " " + sigma_logc + " " + min_logc + " " + max_logc + " " + num_c);
		HistogramFunction hist = new HistogramFunction(min_logc, max_logc, num_c);
		String name = "logc-value prior distribution";
		
		for(int cIndex=0;cIndex<num_c;cIndex++) {
			likelihood[cIndex] = 1/Math.sqrt(2*Math.PI)/sigma_logc * Math.exp(- (mean_logc - hist.getX(cIndex))*(mean_logc - hist.getX(cIndex))/2/sigma_logc/sigma_logc);   
			likelihoodTotal += likelihood[cIndex];
		}
		
		for (int cIndex=0;cIndex<num_c;cIndex++) {
			likelihood[cIndex] /= likelihoodTotal;
			hist.add(cIndex,  likelihood[cIndex]);
		}
		
		hist.setName(name);
		if(D) System.out.println("PDF of logc-value: totalTest = "+hist.calcSumOfY_Vals());
		
		hist.scale(1d/hist.getDelta());
		return hist;
	}
	
}
