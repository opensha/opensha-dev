package scratch.aftershockStatistics;

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

	public static final long MARSHAL_RJGEN_NULL = 5000L;
	public static final long MARSHAL_RJGEN_VER = 5001L;

	// Marshal object.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_RJGEN_VER);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalDouble (mean_a              );
		writer.marshalDouble (sigma_a             );
	
		return;
	}

	// Marshal object.

	public static void marshal (MarshalWriter writer, RJ_Summary_Generic obj) {

		if (obj == null) {
			writer.marshalLong (MARSHAL_RJGEN_NULL);
		} else {
			obj.do_marshal (writer);
		}

		return;
	}

	// Unmarshal object.

	public static RJ_Summary_Generic unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_RJGEN_NULL, MARSHAL_RJGEN_VER);

		if (ver == MARSHAL_RJGEN_NULL) {
			return null;
		}

		return new RJ_Summary_Generic (ver, reader);
	}

	protected RJ_Summary_Generic (MarshalReader reader) {
		this (reader.unmarshalLong (MARSHAL_RJGEN_VER, MARSHAL_RJGEN_VER), reader);
	}

	private RJ_Summary_Generic (long ver, MarshalReader reader) {
		super (reader);

		// Contents

		mean_a               = reader.unmarshalDouble ();
		sigma_a              = reader.unmarshalDouble ();
	}
	
}
