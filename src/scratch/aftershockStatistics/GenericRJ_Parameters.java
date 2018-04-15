package scratch.aftershockStatistics;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

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

	private static final int MARSHAL_VER_1 = 1001;

	private static final String M_VERSION_NAME = "GenericRJ_Parameters";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 1000;
	protected static final int MARSHAL_GENERIC_RJ = 1001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_GENERIC_RJ;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalDouble ("aValue_mean"  , aValue_mean  );
		writer.marshalDouble ("aValue_sigma" , aValue_sigma );
		writer.marshalDouble ("aValue_sigma0", aValue_sigma0);
		writer.marshalDouble ("aValue_sigma1", aValue_sigma1);
		writer.marshalDouble ("bValue"       , bValue       );
		writer.marshalDouble ("pValue"       , pValue       );
		writer.marshalDouble ("cValue"       , cValue       );
		writer.marshalDouble ("aValue_min"   , aValue_min   );
		writer.marshalDouble ("aValue_max"   , aValue_max   );
		writer.marshalDouble ("aValue_delta" , aValue_delta );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		aValue_mean   = reader.unmarshalDouble ("aValue_mean"  );
		aValue_sigma  = reader.unmarshalDouble ("aValue_sigma" );
		aValue_sigma0 = reader.unmarshalDouble ("aValue_sigma0");
		aValue_sigma1 = reader.unmarshalDouble ("aValue_sigma1");
		bValue        = reader.unmarshalDouble ("bValue"       );
		pValue        = reader.unmarshalDouble ("pValue"       );
		cValue        = reader.unmarshalDouble ("cValue"       );
		aValue_min    = reader.unmarshalDouble ("aValue_min"   );
		aValue_max    = reader.unmarshalDouble ("aValue_max"   );
		aValue_delta  = reader.unmarshalDouble ("aValue_delta" );

		return;
	}

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public GenericRJ_Parameters unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, GenericRJ_Parameters obj) {

		writer.marshalMapBegin (name);

		if (obj == null) {
			writer.marshalInt (M_TYPE_NAME, MARSHAL_NULL);
		} else {
			writer.marshalInt (M_TYPE_NAME, obj.get_marshal_type());
			obj.do_marshal (writer);
		}

		writer.marshalMapEnd ();

		return;
	}

	// Unmarshal object, polymorphic.

	public static GenericRJ_Parameters unmarshal_poly (MarshalReader reader, String name) {
		GenericRJ_Parameters result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("GenericRJ_Parameters.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_GENERIC_RJ:
			result = new GenericRJ_Parameters();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}
	
}
