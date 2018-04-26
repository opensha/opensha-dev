package scratch.aftershockStatistics;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

public class RJ_Summary {

	// Summary values, see RJ_AftershockModel for description.

	private double b                   ;
	private double magMain             ;
	private double min_a               ;
	private double max_a               ;
	private double delta_a             ;
	private int    num_a               ;
	private double min_p               ;
	private double max_p               ;
	private double delta_p             ;
	private int    num_p               ;
	private double min_c               ;
	private double max_c               ;
	private double delta_c             ;
	private int    num_c               ;
	private int    apc_total_size      ;
	private int    apc_support_size    ;
	private double apc_support_total   ;
	private double apc_max_tail_element;
	private int    a_support_lo        ;
	private int    a_support_hi        ;
	private int    p_support_lo        ;
	private int    p_support_hi        ;
	private int    c_support_lo        ;
	private int    c_support_hi        ;
	private double stat_a_mean         ;
	private double stat_a_sdev         ;
	private double stat_a_like         ;
	private double stat_p_mean         ;
	private double stat_p_sdev         ;
	private double stat_p_like         ;
	private double stat_c_mean         ;
	private double stat_c_sdev         ;
	private double stat_c_like         ;

	
	/**
	 * Default constructor.
	 */
	public RJ_Summary(){}

	
	/**
	 * Construct from R&J model.
	 */
	public RJ_Summary(RJ_AftershockModel model) {

		this.b                    = model.b                   ;
		this.magMain              = model.magMain             ;
		this.min_a                = model.min_a               ;
		this.max_a                = model.max_a               ;
		this.delta_a              = model.delta_a             ;
		this.num_a                = model.num_a               ;
		this.min_p                = model.min_p               ;
		this.max_p                = model.max_p               ;
		this.delta_p              = model.delta_p             ;
		this.num_p                = model.num_p               ;
		this.min_c                = model.min_c               ;
		this.max_c                = model.max_c               ;
		this.delta_c              = model.delta_c             ;
		this.num_c                = model.num_c               ;
		this.apc_total_size       = model.apc_total_size      ;
		this.apc_support_size     = model.apc_support_size    ;
		this.apc_support_total    = model.apc_support_total   ;
		this.apc_max_tail_element = model.apc_max_tail_element;
		this.a_support_lo         = model.a_support_lo        ;
		this.a_support_hi         = model.a_support_hi        ;
		this.p_support_lo         = model.p_support_lo        ;
		this.p_support_hi         = model.p_support_hi        ;
		this.c_support_lo         = model.c_support_lo        ;
		this.c_support_hi         = model.c_support_hi        ;
		this.stat_a_mean          = model.stat_a_mean         ;
		this.stat_a_sdev          = model.stat_a_sdev         ;
		this.stat_a_like          = model.stat_a_like         ;
		this.stat_p_mean          = model.stat_p_mean         ;
		this.stat_p_sdev          = model.stat_p_sdev         ;
		this.stat_p_like          = model.stat_p_like         ;
		this.stat_c_mean          = model.stat_c_mean         ;
		this.stat_c_sdev          = model.stat_c_sdev         ;
		this.stat_c_like          = model.stat_c_like         ;
	}

	
	/**
	 * Getters.
	 */
	public double get_b                   () {return b                   ;}
	public double get_magMain             () {return magMain             ;}
	public double get_min_a               () {return min_a               ;}
	public double get_max_a               () {return max_a               ;}
	public double get_delta_a             () {return delta_a             ;}
	public int    get_num_a               () {return num_a               ;}
	public double get_min_p               () {return min_p               ;}
	public double get_max_p               () {return max_p               ;}
	public double get_delta_p             () {return delta_p             ;}
	public int    get_num_p               () {return num_p               ;}
	public double get_min_c               () {return min_c               ;}
	public double get_max_c               () {return max_c               ;}
	public double get_delta_c             () {return delta_c             ;}
	public int    get_num_c               () {return num_c               ;}
	public int    get_apc_total_size      () {return apc_total_size      ;}
	public int    get_apc_support_size    () {return apc_support_size    ;}
	public double get_apc_support_total   () {return apc_support_total   ;}
	public double get_apc_max_tail_element() {return apc_max_tail_element;}
	public int    get_a_support_lo        () {return a_support_lo        ;}
	public int    get_a_support_hi        () {return a_support_hi        ;}
	public int    get_p_support_lo        () {return p_support_lo        ;}
	public int    get_p_support_hi        () {return p_support_hi        ;}
	public int    get_c_support_lo        () {return c_support_lo        ;}
	public int    get_c_support_hi        () {return c_support_hi        ;}
	public double get_stat_a_mean         () {return stat_a_mean         ;}
	public double get_stat_a_sdev         () {return stat_a_sdev         ;}
	public double get_stat_a_like         () {return stat_a_like         ;}
	public double get_stat_p_mean         () {return stat_p_mean         ;}
	public double get_stat_p_sdev         () {return stat_p_sdev         ;}
	public double get_stat_p_like         () {return stat_p_like         ;}
	public double get_stat_c_mean         () {return stat_c_mean         ;}
	public double get_stat_c_sdev         () {return stat_c_sdev         ;}
	public double get_stat_c_like         () {return stat_c_like         ;}



	@Override
	public String toString() {
		return
			"b                    = " + get_b                   () + "\n" +
			"magMain              = " + get_magMain             () + "\n" +
			"min_a                = " + get_min_a               () + "\n" +
			"max_a                = " + get_max_a               () + "\n" +
			"delta_a              = " + get_delta_a             () + "\n" +
			"num_a                = " + get_num_a               () + "\n" +
			"min_p                = " + get_min_p               () + "\n" +
			"max_p                = " + get_max_p               () + "\n" +
			"delta_p              = " + get_delta_p             () + "\n" +
			"num_p                = " + get_num_p               () + "\n" +
			"min_c                = " + get_min_c               () + "\n" +
			"max_c                = " + get_max_c               () + "\n" +
			"delta_c              = " + get_delta_c             () + "\n" +
			"num_c                = " + get_num_c               () + "\n" +
			"apc_total_size       = " + get_apc_total_size      () + "\n" +
			"apc_support_size     = " + get_apc_support_size    () + "\n" +
			"apc_support_total    = " + get_apc_support_total   () + "\n" +
			"apc_max_tail_element = " + get_apc_max_tail_element() + "\n" +
			"a_support_lo         = " + get_a_support_lo        () + "\n" +
			"a_support_hi         = " + get_a_support_hi        () + "\n" +
			"p_support_lo         = " + get_p_support_lo        () + "\n" +
			"p_support_hi         = " + get_p_support_hi        () + "\n" +
			"c_support_lo         = " + get_c_support_lo        () + "\n" +
			"c_support_hi         = " + get_c_support_hi        () + "\n" +
			"stat_a_mean          = " + get_stat_a_mean         () + "\n" +
			"stat_a_sdev          = " + get_stat_a_sdev         () + "\n" +
			"stat_a_like          = " + get_stat_a_like         () + "\n" +
			"stat_p_mean          = " + get_stat_p_mean         () + "\n" +
			"stat_p_sdev          = " + get_stat_p_sdev         () + "\n" +
			"stat_p_like          = " + get_stat_p_like         () + "\n" +
			"stat_c_mean          = " + get_stat_c_mean         () + "\n" +
			"stat_c_sdev          = " + get_stat_c_sdev         () + "\n" +
			"stat_c_like          = " + get_stat_c_like         () ;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 4001;

	private static final String M_VERSION_NAME = "RJ_Summary";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 4000;
	protected static final int MARSHAL_RJBASE = 4001;
	protected static final int MARSHAL_RJGEN = 5001;
	protected static final int MARSHAL_RJBAY = 6001;
	protected static final int MARSHAL_RJSEQ = 7001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_RJBASE;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalDouble ("b"                   , b                   );
		writer.marshalDouble ("magMain"             , magMain             );
		writer.marshalDouble ("min_a"               , min_a               );
		writer.marshalDouble ("max_a"               , max_a               );
		writer.marshalDouble ("delta_a"             , delta_a             );
		writer.marshalInt    ("num_a"               , num_a               );
		writer.marshalDouble ("min_p"               , min_p               );
		writer.marshalDouble ("max_p"               , max_p               );
		writer.marshalDouble ("delta_p"             , delta_p             );
		writer.marshalInt    ("num_p"               , num_p               );
		writer.marshalDouble ("min_c"               , min_c               );
		writer.marshalDouble ("max_c"               , max_c               );
		writer.marshalDouble ("delta_c"             , delta_c             );
		writer.marshalInt    ("num_c"               , num_c               );
		writer.marshalInt    ("apc_total_size"      , apc_total_size      );
		writer.marshalInt    ("apc_support_size"    , apc_support_size    );
		writer.marshalDouble ("apc_support_total"   , apc_support_total   );
		writer.marshalDouble ("apc_max_tail_element", apc_max_tail_element);
		writer.marshalInt    ("a_support_lo"        , a_support_lo        );
		writer.marshalInt    ("a_support_hi"        , a_support_hi        );
		writer.marshalInt    ("p_support_lo"        , p_support_lo        );
		writer.marshalInt    ("p_support_hi"        , p_support_hi        );
		writer.marshalInt    ("c_support_lo"        , c_support_lo        );
		writer.marshalInt    ("c_support_hi"        , c_support_hi        );
		writer.marshalDouble ("stat_a_mean"         , stat_a_mean         );
		writer.marshalDouble ("stat_a_sdev"         , stat_a_sdev         );
		writer.marshalDouble ("stat_a_like"         , stat_a_like         );
		writer.marshalDouble ("stat_p_mean"         , stat_p_mean         );
		writer.marshalDouble ("stat_p_sdev"         , stat_p_sdev         );
		writer.marshalDouble ("stat_p_like"         , stat_p_like         );
		writer.marshalDouble ("stat_c_mean"         , stat_c_mean         );
		writer.marshalDouble ("stat_c_sdev"         , stat_c_sdev         );
		writer.marshalDouble ("stat_c_like"         , stat_c_like         );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		b                    = reader.unmarshalDouble ("b"                   );
		magMain              = reader.unmarshalDouble ("magMain"             );
		min_a                = reader.unmarshalDouble ("min_a"               );
		max_a                = reader.unmarshalDouble ("max_a"               );
		delta_a              = reader.unmarshalDouble ("delta_a"             );
		num_a                = reader.unmarshalInt    ("num_a"               , 1);
		min_p                = reader.unmarshalDouble ("min_p"               );
		max_p                = reader.unmarshalDouble ("max_p"               );
		delta_p              = reader.unmarshalDouble ("delta_p"             );
		num_p                = reader.unmarshalInt    ("num_p"               , 1);
		min_c                = reader.unmarshalDouble ("min_c"               );
		max_c                = reader.unmarshalDouble ("max_c"               );
		delta_c              = reader.unmarshalDouble ("delta_c"             );
		num_c                = reader.unmarshalInt    ("num_c"               , 1);
		apc_total_size       = reader.unmarshalInt    ("apc_total_size"      , 1);
		apc_support_size     = reader.unmarshalInt    ("apc_support_size"    , 0);
		apc_support_total    = reader.unmarshalDouble ("apc_support_total"   );
		apc_max_tail_element = reader.unmarshalDouble ("apc_max_tail_element");
		a_support_lo         = reader.unmarshalInt    ("a_support_lo"        , 0);
		a_support_hi         = reader.unmarshalInt    ("a_support_hi"        , 0);
		p_support_lo         = reader.unmarshalInt    ("p_support_lo"        , 0);
		p_support_hi         = reader.unmarshalInt    ("p_support_hi"        , 0);
		c_support_lo         = reader.unmarshalInt    ("c_support_lo"        , 0);
		c_support_hi         = reader.unmarshalInt    ("c_support_hi"        , 0);
		stat_a_mean          = reader.unmarshalDouble ("stat_a_mean"         );
		stat_a_sdev          = reader.unmarshalDouble ("stat_a_sdev"         );
		stat_a_like          = reader.unmarshalDouble ("stat_a_like"         );
		stat_p_mean          = reader.unmarshalDouble ("stat_p_mean"         );
		stat_p_sdev          = reader.unmarshalDouble ("stat_p_sdev"         );
		stat_p_like          = reader.unmarshalDouble ("stat_p_like"         );
		stat_c_mean          = reader.unmarshalDouble ("stat_c_mean"         );
		stat_c_sdev          = reader.unmarshalDouble ("stat_c_sdev"         );
		stat_c_like          = reader.unmarshalDouble ("stat_c_like"         );

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

	public RJ_Summary unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, RJ_Summary obj) {

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

	public static RJ_Summary unmarshal_poly (MarshalReader reader, String name) {
		RJ_Summary result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("RJ_Summary.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_RJBASE:
			result = new RJ_Summary();
			result.do_umarshal (reader);
			break;

		case MARSHAL_RJGEN:
			result = new RJ_Summary_Generic();
			result.do_umarshal (reader);
			break;

		case MARSHAL_RJBAY:
			result = new RJ_Summary_Bayesian();
			result.do_umarshal (reader);
			break;

		case MARSHAL_RJSEQ:
			result = new RJ_Summary_SequenceSpecific();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}
	
}
