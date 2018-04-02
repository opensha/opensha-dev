package scratch.aftershockStatistics;

public class SeqSpecRJ_Parameters {

	// Parameter values, see RJ_AftershockModel_SequenceSpecific for description.

	private double b;
	private double min_a;
	private double max_a;
	private int    num_a;
	private double min_p;
	private double max_p;
	private int    num_p;
	private double min_c;
	private double max_c;
	private int    num_c;

	
	/**
	 * This class is a container for the Sequence-Specific Reasenberg-Jones parameters.
	 */
	public SeqSpecRJ_Parameters(
			double b,
			double min_a,
			double max_a,
			int    num_a,
			double min_p,
			double max_p,
			int    num_p,
			double min_c,
			double max_c,
			int    num_c) {

		this.b     = b;
		this.min_a = min_a;
		this.max_a = max_a;
		this.num_a = num_a;
		this.min_p = min_p;
		this.max_p = max_p;
		this.num_p = num_p;
		this.min_c = min_c;
		this.max_c = max_c;
		this.num_c = num_c;
	}

	
	/**
	 * Default constructor.
	 */
	public SeqSpecRJ_Parameters(){}

	
	/**
	 * Construct from a set of generic parameters.
	 */
	public SeqSpecRJ_Parameters(GenericRJ_Parameters param) {

		this.b     = param.get_bValue();

		double delta_a = param.get_aValue_delta();
		this.min_a = ((double)Math.round(param.get_aValue_min()/delta_a))*delta_a;	// round to nearest multiple of delta_a
		this.max_a = ((double)Math.round(param.get_aValue_max()/delta_a))*delta_a;	// round to nearest multiple of delta_a
		this.num_a = (int)Math.round((this.max_a-this.min_a)/delta_a) + 1;

		this.min_p = param.get_pValue();
		this.max_p = param.get_pValue();
		this.num_p = 1;

		this.min_c = param.get_cValue();
		this.max_c = param.get_cValue();
		this.num_c = 1;
	}

	
	/**
	 * Getters.
	 */
	public double get_b    () {return b    ;}
	public double get_min_a() {return min_a;}
	public double get_max_a() {return max_a;}
	public int    get_num_a() {return num_a;}
	public double get_min_p() {return min_p;}
	public double get_max_p() {return max_p;}
	public int    get_num_p() {return num_p;}
	public double get_min_c() {return min_c;}
	public double get_max_c() {return max_c;}
	public int    get_num_c() {return num_c;}


	@Override
	public String toString() {
		return "SQ_Params["
				+       "b=" + get_b    ()
				+ ", min_a=" + get_min_a()
				+ ", max_a=" + get_max_a()
				+ ", num_a=" + get_num_a()
				+ ", min_p=" + get_min_p()
				+ ", max_p=" + get_max_p()
				+ ", num_p=" + get_num_p()
				+ ", min_c=" + get_min_c()
				+ ", max_c=" + get_max_c()
				+ ", num_c=" + get_num_c()
				+ "]";
	}




	//----- Marshaling -----

	// Marshal version number.

	public static final long MARSHAL_NULL = 3000L;
	public static final long MARSHAL_VER = 3001L;

	// Marshal object.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_VER);

		// Contents

		writer.marshalDouble (b    );
		writer.marshalDouble (min_a);
		writer.marshalDouble (max_a);
		writer.marshalInt    (num_a);
		writer.marshalDouble (min_p);
		writer.marshalDouble (max_p);
		writer.marshalInt    (num_p);
		writer.marshalDouble (min_c);
		writer.marshalDouble (max_c);
		writer.marshalInt    (num_c);
	
		return;
	}

	// Marshal object.

	public static void marshal (MarshalWriter writer, SeqSpecRJ_Parameters obj) {

		if (obj == null) {
			writer.marshalLong (MARSHAL_NULL);
		} else {
			obj.do_marshal (writer);
		}

		return;
	}

	// Unmarshal object.

	public static SeqSpecRJ_Parameters unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_NULL, MARSHAL_VER);

		if (ver == MARSHAL_NULL) {
			return null;
		}

		return new SeqSpecRJ_Parameters (ver, reader);
	}

	private SeqSpecRJ_Parameters (long ver, MarshalReader reader) {

		// Contents

		b     = reader.unmarshalDouble();
		min_a = reader.unmarshalDouble();
		max_a = reader.unmarshalDouble();
		num_a = reader.unmarshalInt(1);
		min_p = reader.unmarshalDouble();
		max_p = reader.unmarshalDouble();
		num_p = reader.unmarshalInt(1);
		min_c = reader.unmarshalDouble();
		max_c = reader.unmarshalDouble();
		num_c = reader.unmarshalInt(1);
	}
	
}
