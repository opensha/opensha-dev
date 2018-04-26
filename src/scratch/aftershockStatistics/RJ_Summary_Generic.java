package scratch.aftershockStatistics;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

public class RJ_Summary_Generic extends RJ_Summary {

	// Summary values, see RJ_AftershockModel_Generic for description.

	private double mean_a              ;
	private double sigma_a             ;

	
	/**
	 * Default constructor.
	 */
	public RJ_Summary_Generic() {
		super();
	}

	
	/**
	 * Construct from R&J model.
	 */
	public RJ_Summary_Generic(RJ_AftershockModel_Generic model) {
		super (model);

		this.mean_a               = model.mean_a              ;
		this.sigma_a              = model.sigma_a             ;
	}

	
	/**
	 * Getters.
	 */
	public double get_mean_a              () {return mean_a              ;}
	public double get_sigma_a             () {return sigma_a             ;}



	@Override
	public String toString() {
		return super.toString() + "\n" +
			"mean_a               = " + get_mean_a              () + "\n" +
			"sigma_a              = " + get_sigma_a             () ;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 5001;

	private static final String M_VERSION_NAME = "RJ_Summary_Generic";

	// Get the type code.

	@Override
	protected int get_marshal_type () {
		return MARSHAL_RJGEN;
	}

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalDouble ("mean_a"              , mean_a              );
		writer.marshalDouble ("sigma_a"             , sigma_a             );
	
		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Superclass

		super.do_umarshal (reader);

		// Contents

		mean_a               = reader.unmarshalDouble ("mean_a"              );
		sigma_a              = reader.unmarshalDouble ("sigma_a"             );

		return;
	}

	// Marshal object.

	@Override
	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	@Override
	public RJ_Summary_Generic unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, RJ_Summary_Generic obj) {

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

	public static RJ_Summary_Generic unmarshal_poly (MarshalReader reader, String name) {
		RJ_Summary_Generic result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("RJ_Summary_Generic.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_RJGEN:
			result = new RJ_Summary_Generic();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}
	
}
