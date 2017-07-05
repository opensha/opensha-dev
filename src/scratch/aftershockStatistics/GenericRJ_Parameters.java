package scratch.aftershockStatistics;

public class GenericRJ_Parameters {
	
	double aValue_mean;
	double aValue_sigma;	// magnitude independent sigma
	double aValue_sigma0;	// param for mag-depednen
	double aValue_sigma1;
	double bValue;
	double pValue;
	double cValue;
	
	/**
	 * This class is a container for the Generic Reasenberg-Jones parameters defined by 
	 * Page et al. (2016, FILL IN REF AFTER PUBLICATION).  The distribution of a-values is assumed
	 * to be Gaussian, with a mean of aValue_mean and a standard deviation of aValue_sigma.  
	 * A magnitude-dependent sigma is an option, and can be computed as:
	 * 
	 * 		sigma(M) = (sigma0^2 + sigma1^2/10^M)^0.5 if M≥6
	 * or
	 * 		sigma(M) = (sigma0^2 + sigma1^2/10^6)^0.5 if M<6
	 * 
	 * @param aValue_mean
	 * @param aValue_sigma
	 * @param aValue_sigma0
	 * @param aValue_sigma1
	 * @param bValue
	 * @param pValue
	 * @param cValue
	 */
	public GenericRJ_Parameters(double aValue_mean, double aValue_sigma, double aValue_sigma0, double aValue_sigma1, 
			double bValue, double pValue, double cValue) {
		this.aValue_mean = aValue_mean;
		this.aValue_sigma = aValue_sigma;	// magnitude independent sigma
		this.aValue_sigma0 = aValue_sigma0;	// param for mag-depednen
		this.aValue_sigma1 = aValue_sigma1;
		this.bValue = bValue;
		this.pValue = pValue;
		this.cValue = cValue;
	}

	public GenericRJ_Parameters(){}
	
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
	 * This returns the magnitude-dependent a-value standard deviation.
	 * @return
	 */
	public double get_aValueSigma(double magnitude) {
		if(magnitude>=6.0)
			return Math.sqrt(aValue_sigma0*aValue_sigma0 + aValue_sigma1*aValue_sigma1/Math.pow(10, magnitude));
		else
			return Math.sqrt(aValue_sigma0*aValue_sigma0 + aValue_sigma1*aValue_sigma1/Math.pow(10, 6.0));
	}
	
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

	/**
	 * This returns the b-value.
	 * @return
	 */
	public double get_cValue() {return cValue;}

	@Override
	public String toString() {
		return "RJ_Params[a="+get_aValueMean()+", aSigma="+get_aValueSigma()+", aSigma0="+aValue_sigma0
			+", aSigma1="+aValue_sigma1+", b="+get_bValue()+", p="+get_pValue()+", c="+get_cValue()+"]";
	}
	
}
