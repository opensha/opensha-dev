package scratch.aftershockStatistics;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

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

	public static final int MARSHAL_VER_1 = 3001;

	private static final String M_VERSION_NAME = "SeqSpecRJ_Parameters";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 3000;
	protected static final int MARSHAL_SEQ_SPEC_RJ = 3001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_SEQ_SPEC_RJ;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalDouble ("b"    , b    );
		writer.marshalDouble ("min_a", min_a);
		writer.marshalDouble ("max_a", max_a);
		writer.marshalInt    ("num_a", num_a);
		writer.marshalDouble ("min_p", min_p);
		writer.marshalDouble ("max_p", max_p);
		writer.marshalInt    ("num_p", num_p);
		writer.marshalDouble ("min_c", min_c);
		writer.marshalDouble ("max_c", max_c);
		writer.marshalInt    ("num_c", num_c);
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		b     = reader.unmarshalDouble ("b"    );
		min_a = reader.unmarshalDouble ("min_a");
		max_a = reader.unmarshalDouble ("max_a");
		num_a = reader.unmarshalInt    ("num_a", 1);
		min_p = reader.unmarshalDouble ("min_p");
		max_p = reader.unmarshalDouble ("max_p");
		num_p = reader.unmarshalInt    ("num_p", 1);
		min_c = reader.unmarshalDouble ("min_c");
		max_c = reader.unmarshalDouble ("max_c");
		num_c = reader.unmarshalInt    ("num_c", 1);

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

	public SeqSpecRJ_Parameters unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, SeqSpecRJ_Parameters obj) {

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

	public static SeqSpecRJ_Parameters unmarshal_poly (MarshalReader reader, String name) {
		SeqSpecRJ_Parameters result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("SeqSpecRJ_Parameters.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_SEQ_SPEC_RJ:
			result = new SeqSpecRJ_Parameters();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}
	
}
