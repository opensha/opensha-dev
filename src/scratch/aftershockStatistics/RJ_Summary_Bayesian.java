package scratch.aftershockStatistics;

public class RJ_Summary_Bayesian extends RJ_Summary {

	// Summary values, see RJ_AftershockModel_Bayesian for description.


	
	/**
	 * Default constructor.
	 */
	public RJ_Summary_Bayesian() {
		super();
	}

	
	/**
	 * Construct from R&J model.
	 */
	public RJ_Summary_Bayesian(RJ_AftershockModel_Bayesian model) {
		super (model);

	}

	
	/**
	 * Getters.
	 */



	@Override
	public String toString() {
		return super.toString();
	}




	//----- Marshaling -----

	// Marshal version number.

	public static final long MARSHAL_RJBAY_NULL = 6000L;
	public static final long MARSHAL_RJBAY_VER = 6001L;

	// Marshal object.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_RJBAY_VER);

		// Superclass

		super.do_marshal (writer);

		// Contents

	
		return;
	}

	// Marshal object.

	public static void marshal (MarshalWriter writer, RJ_Summary_Bayesian obj) {

		if (obj == null) {
			writer.marshalLong (MARSHAL_RJBAY_NULL);
		} else {
			obj.do_marshal (writer);
		}

		return;
	}

	// Unmarshal object.

	public static RJ_Summary_Bayesian unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_RJBAY_NULL, MARSHAL_RJBAY_VER);

		if (ver == MARSHAL_RJBAY_NULL) {
			return null;
		}

		return new RJ_Summary_Bayesian (ver, reader);
	}

	protected RJ_Summary_Bayesian (MarshalReader reader) {
		this (reader.unmarshalLong (MARSHAL_RJBAY_VER, MARSHAL_RJBAY_VER), reader);
	}

	private RJ_Summary_Bayesian (long ver, MarshalReader reader) {
		super (reader);

		// Contents

	}
	
}
