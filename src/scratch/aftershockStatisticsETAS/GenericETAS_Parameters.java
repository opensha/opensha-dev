package scratch.aftershockStatisticsETAS;

public class GenericETAS_Parameters implements java.io.Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 578249832338308677L;
	
	/**
	 * 
	 */
	
	double aValue_mean;
	double aValue_sigma;	// magnitude independent sigma
	double bValue;
	double pValue;
	double cValue;
	double alpha;
	double refMag;
	double[][] covariance;
	int numSequences;
	
	/**
	 * This class is a container for the Generic ETAS parameters defined by van der Elst and Page (in prep).
	 * In development -- The distribution of a-values is assumed
	 * to be Gaussian, with a mean of aValue_mean and a standard deviation of aValue_sigma.  
	 * For now, alpha (the magnitude scaling exponent) is set to 1.0. Model is limited to k,c, and p. 
	 * 
	 * @param aValue_mean
	 * @param aValue_sigma
	 * @param bValue
	 * @param alpha
	 * @param pValue
	 * @param cValue
	 * @param refMag
	 */
	public GenericETAS_Parameters() {
		//if called without arguments, initialize with global average parameters
		this(-2.423,	0.395,	-2.565,	0.966,	7.13E-06,	7.67E-06,	8.99E-04,	2.63E-06,	3.75E-05,	4.50E-05, 2099, 1, 1, 4.5);

		//		this(-2.43, 1, 0.96, 0.13, Math.pow(10,  -2.565), 0.75, 1, 2.63e-6, 3.75e-5, 4.50e-5, 2099, 1, 1, 4.5);
		//		this.aValue_mean = -2.43;
		//		this.aValue_sigma = 0.4;	// magnitude independent sigma
		//		this.bValue = 1;
		//		this.pValue = 0.96;
		//		this.pValue_sigma = 0.13;
		//		this.cValue = Math.pow(10, -2.64);
		//		this.log_cValue_sigma = 0.75;
		//		this.alpha = 1;
		//		this.refMag = 4.5;
		//		this.covariance = covariance;
	}
	
	public GenericETAS_Parameters(double aValue_mean, double aValue_sigma,  double pValue, double cValue, 
			double covaa, double covpp, double covcc, double covap, double covac, double covcp,
			int numSequences, double alpha, double bValue, double refMag) {
	
			double[][] covariance = new double[][]{{covaa, covap, covac},{covap, covpp, covcp},{covac, covcp, covpp}};  
			
			
			this.aValue_mean = aValue_mean;
			this.aValue_sigma = aValue_sigma;	// magnitude independent sigma
			this.bValue = bValue;
			this.pValue = pValue;
			this.cValue = cValue;
			this.alpha = alpha;
			this.covariance = covariance;
			this.refMag = refMag;
	}
	


	/**
	 * This returns the mean a-value (aValue_mean).
	 * @return
	 */
	public double get_aValueMean() {return aValue_mean;}
	
	/**
	 * This returns the magnitude-independent a-value standard deviation (aValue_sigma).
	 * @return
	 */
	public double get_aValueSigma() {return aValue_sigma;}
	
	/**
	 * This returns the b-value.
	 * @return
	 */
	public double get_bValue() {return bValue;}

	/**
	 * This returns the p-value.
	 * @return
	 */
	public double get_pValue() {return pValue;}

	public double get_pValueSigma() {
		double pValue_sigma = Math.sqrt(covariance[1][1] * (double) numSequences); 
		return pValue_sigma;
	}

	/**
	 * This returns the b-value.
	 * @return
	 */
	public double get_cValue() {return cValue;}

	public double get_logcValueSigma() {
		double cValue_sigma = Math.sqrt(covariance[2][2] * (double) numSequences); 
		return cValue_sigma;
	}

	
	public double get_alpha() {return alpha;}
	
	public double get_refMag() {return refMag;}

	@Override
	public String toString() {
		return "ETAS_Params[a="+get_aValueMean()+", aSigma="+get_aValueSigma()+","
				+ " b="+get_bValue()+", p="+get_pValue()+", pSigma="+get_pValueSigma()+","
				+ " c="+get_cValue()+", logcSigma="+get_logcValueSigma()+","
				+ " alpha="+get_alpha()+", refMag="+get_refMag()+"]";
	}
	
}
