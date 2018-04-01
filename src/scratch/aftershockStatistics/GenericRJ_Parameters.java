package scratch.aftershockStatistics;

public class GenericRJ_Parameters {
	
	private double aValue_mean;
	private double aValue_sigma;	// magnitude independent sigma
	private double aValue_sigma0;	// param for mag-depednen
	private double aValue_sigma1;
	private double bValue;
	private double pValue;
	private double cValue;
	private double aValue_min;		// minimum a-value to consider
	private double aValue_max;		// maximum a-value to consider
	private double aValue_delta;	// spacing between a-values to consider
	
	/**
	 * This class is a container for the Generic Reasenberg-Jones parameters defined by 
	 * Page et al. (2016, BSSA).  The distribution of a-values is assumed
	 * to be Gaussian, with a mean of aValue_mean and a standard deviation of aValue_sigma.  
	 * A magnitude-dependent sigma is an option, and can be computed as:
	 * 
	 * 		sigma(M) = (sigma0^2 + sigma1^2/10^M)^0.5 if M >= 6
	 * or
	 * 		sigma(M) = (sigma0^2 + sigma1^2/10^6)^0.5 if M < 6
	 * 
	 * @param aValue_mean
	 * @param aValue_sigma
	 * @param aValue_sigma0
	 * @param aValue_sigma1
	 * @param bValue
	 * @param pValue
	 * @param cValue
	 * @param aValue_min
	 * @param aValue_max
	 * @param aValue_delta
	 */
	public GenericRJ_Parameters(double aValue_mean, double aValue_sigma, double aValue_sigma0, double aValue_sigma1, 
			double bValue, double pValue, double cValue, double aValue_min, double aValue_max, double aValue_delta) {
		this.aValue_mean = aValue_mean;
		this.aValue_sigma = aValue_sigma;	// magnitude independent sigma
		this.aValue_sigma0 = aValue_sigma0;	// param for mag-depednen
		this.aValue_sigma1 = aValue_sigma1;
		this.bValue = bValue;
		this.pValue = pValue;
		this.cValue = cValue;
		this.aValue_min = aValue_min;
		this.aValue_max = aValue_max;
		this.aValue_delta = aValue_delta;
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
	
	/**
	 * This returns the minimum a-value to consider (aValue_min).
	 * @return
	 */
	public double get_aValue_min() {return aValue_min;}
	
	/**
	 * This returns the maximum a-value to consider (aValue_max).
	 * @return
	 */
	public double get_aValue_max() {return aValue_max;}
	
	/**
	 * This returns the spacing between a-values to consider (aValue_delta).
	 * @return
	 */
	public double get_aValue_delta() {return aValue_delta;}

	@Override
	public String toString() {
		return "RJ_Params[a=" + get_aValueMean()
				+ ", aSigma=" + get_aValueSigma()
				+ ", aSigma0=" + aValue_sigma0
				+ ", aSigma1=" + aValue_sigma1
				+ ", b=" + get_bValue()
				+ ", p=" + get_pValue()
				+ ", c=" + get_cValue()
				+ ", aMin=" + get_aValue_min()
				+ ", aMax=" + get_aValue_max()
				+ ", aDelta=" + get_aValue_delta()
				+ "]";
	}




	//----- Marshaling -----

	// Marshal version number.

	public static final long MARSHAL_VER = 1001L;

	// Marshal object.

	public void marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_VER);

		// Contents

		writer.marshalDouble (aValue_mean  );
		writer.marshalDouble (aValue_sigma );
		writer.marshalDouble (aValue_sigma0);
		writer.marshalDouble (aValue_sigma1);
		writer.marshalDouble (bValue       );
		writer.marshalDouble (pValue       );
		writer.marshalDouble (cValue       );
		writer.marshalDouble (aValue_min   );
		writer.marshalDouble (aValue_max   );
		writer.marshalDouble (aValue_delta );
	
		return;
	}

	// Unmarshal object.

	public static GenericRJ_Parameters unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_VER, MARSHAL_VER);

		return new GenericRJ_Parameters (ver, reader);
	}

	private GenericRJ_Parameters (long ver, MarshalReader reader) {

		// Contents

		aValue_mean   = reader.unmarshalDouble();
		aValue_sigma  = reader.unmarshalDouble();
		aValue_sigma0 = reader.unmarshalDouble();
		aValue_sigma1 = reader.unmarshalDouble();
		bValue        = reader.unmarshalDouble();
		pValue        = reader.unmarshalDouble();
		cValue        = reader.unmarshalDouble();
		aValue_min    = reader.unmarshalDouble();
		aValue_max    = reader.unmarshalDouble();
		aValue_delta  = reader.unmarshalDouble();
	}
	
}
