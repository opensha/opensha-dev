package scratch.aftershockStatisticsETAS;

public class GenericETAS_Parameters implements java.io.Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 578249832338308677L;
	
	/**
	 * 
	 */
	
	// see how we have all these parameters defined in this parameter set? Change genericETASparams to just ETASParams!
	double amsValue_mean;
	double amsValue_sigma;	// magnitude independent sigma
	double aValue_mean;
	double aValue_sigma;	// magnitude independent sigma
	double pValue;
	double pValue_sigma;
	double cValue;
	double logcValue_sigma;
	double[][] covariance;
	double[][] priorCovariance;
	double alpha;
	double bValue;
	double bValue_sigma;
	double refMag;
	double maxMag;
	
	
	/**
	 * This class is a container for the Generic ETAS parameters defined by van der Elst and Page (in prep).
	 * In development -- The distribution of a-values is assumed
	 * to be Gaussian, with a mean of aValue_mean and a standard deviation of aValue_sigma.  
	 * For now, alpha (the magnitude scaling exponent) is set to 1.0. Model is limited to k,c, and p. 
	 * 
	 * @param aValue_mean
	 * @param aValue_sigma
	 * @param pValue
	 * @param pValue_sigma
	 * @param cValue
	 * @param logcValue_sigma
	 * @param refMag
	 * @param alpha
	 * @param bValue
	 */
	public GenericETAS_Parameters() {
		//if called without arguments, initialize with global average parameters
		this(-2.423, 0.395,	0.966, 0.2, Math.pow(10,-2.565), 0.7, 7.13E-06,	7.67E-06, 8.99E-04,	2.63E-06, 3.75E-05,	4.50E-05, 2300, 1, 1, 4.5);

	}
	
	public GenericETAS_Parameters(double aValue_mean, double aValue_sigma,  double pValue, double pValue_sigma, double cValue, double logcValue_sigma,
			double covaa, double covpp, double covcc, double covap, double covac, double covcp, double nsamples,
			double alpha, double bValue, double refMag) {
	
		this(aValue_mean, aValue_sigma, pValue, pValue_sigma, cValue, logcValue_sigma,
				covaa, covpp, covcc, covap, covac, covcp, nsamples, alpha, bValue, 0.1, refMag, 9.5);
	}
		
	public GenericETAS_Parameters(double aValue_mean, double aValue_sigma,  double pValue, double pValue_sigma, double cValue, double logcValue_sigma,
			double covaa, double covpp, double covcc, double covap, double covac, double covcp, double n,
			double alpha, double bValue, double bValue_sigma, double refMag, double maxMag) {

		double[][] covariance = new double[][]{{covaa, covap, covac},{covap, covpp, covcp},{covac, covcp, covcc}};  
		//			double[][] priorCovariance = new double[][]{{covaa*n, covap*n, covac*n},{covap*n, covpp*n, covcp*n},{covac*n, covcp*n, covcc*n}};
		double[][] priorCovariance = new double[][]{{aValue_sigma*aValue_sigma, covap*n, covac*n},{covap*n, covpp*n, covcp*n},{covac*n, covcp*n, covcc*n}};

		this.amsValue_mean = aValue_mean;
		this.amsValue_sigma = aValue_sigma;	// magnitude independent sigma
		this.aValue_mean = aValue_mean;
		this.aValue_sigma = aValue_sigma;	// magnitude independent sigma
		this.pValue = pValue;
		this.pValue_sigma = pValue_sigma;
		this.cValue = cValue;
		this.logcValue_sigma = logcValue_sigma;
		this.alpha = alpha;
		this.bValue = bValue;
		this.bValue_sigma = bValue_sigma;
		this.covariance = covariance;
		this.priorCovariance = priorCovariance;
		this.refMag = refMag;
		this.maxMag = maxMag;
	}



	/**
	 * This returns the mean a-value (aValue_mean).
	 * @return
	 */
	public double get_ams() {return amsValue_mean;}

	/**
	 * This returns the magnitude-independent a-value standard deviation used in the Bayesian prior (aValue_sigma).
	 * @return
	 */
	public double get_amsSigma(){return amsValue_sigma;}
	
	public double get_a(){return aValue_mean;}
	
	public double get_aSigma(){return aValue_sigma;}

	public double get_p() {return pValue;}

	public double get_pSigma() {return pValue_sigma;}

	public double get_c() {return cValue;}
	
	public double get_logc() {return Math.log10(cValue);}
	
	public double get_logcSigma() {return logcValue_sigma;}

	/**
	 * This returns the covariance matrix used for epistemic uncertainty in the generic forecasts.
	 * @return
	 */
	public double[][] get_covariance() {return covariance;}

	public double[][] get_priorCovariance() {return priorCovariance;}

	public double get_alpha() {return alpha;}

	public double get_b() {return bValue;}
	
	public double get_bSigma() {return bValue_sigma;}

	public double get_refMag() {return refMag;}

	public double get_maxMag() {return maxMag;}

	@Override
	public String toString() {
		return "ETAS_Params[a=" + get_a() + ", aSigma="+get_aSigma() + ","
				+ " p="+get_p() + ", pSigma="+get_pSigma() + ","
				+ " c="+get_c() + ", logcSigma="+get_logcSigma() + ","
				+ " alpha="+get_alpha() + " b="+get_b() + ", refMag="+get_refMag() + "]";
	}

}
