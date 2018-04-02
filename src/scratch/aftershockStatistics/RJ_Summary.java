package scratch.aftershockStatistics;

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

	public static final long MARSHAL_RJBASE_NULL = 4000L;
	public static final long MARSHAL_RJBASE_VER = 4001L;

	// Marshal object.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_RJBASE_VER);

		// Contents

		writer.marshalDouble (b                   );
		writer.marshalDouble (magMain             );
		writer.marshalDouble (min_a               );
		writer.marshalDouble (max_a               );
		writer.marshalDouble (delta_a             );
		writer.marshalInt    (num_a               );
		writer.marshalDouble (min_p               );
		writer.marshalDouble (max_p               );
		writer.marshalDouble (delta_p             );
		writer.marshalInt    (num_p               );
		writer.marshalDouble (min_c               );
		writer.marshalDouble (max_c               );
		writer.marshalDouble (delta_c             );
		writer.marshalInt    (num_c               );
		writer.marshalInt    (apc_total_size      );
		writer.marshalInt    (apc_support_size    );
		writer.marshalDouble (apc_support_total   );
		writer.marshalDouble (apc_max_tail_element);
		writer.marshalInt    (a_support_lo        );
		writer.marshalInt    (a_support_hi        );
		writer.marshalInt    (p_support_lo        );
		writer.marshalInt    (p_support_hi        );
		writer.marshalInt    (c_support_lo        );
		writer.marshalInt    (c_support_hi        );
		writer.marshalDouble (stat_a_mean         );
		writer.marshalDouble (stat_a_sdev         );
		writer.marshalDouble (stat_a_like         );
		writer.marshalDouble (stat_p_mean         );
		writer.marshalDouble (stat_p_sdev         );
		writer.marshalDouble (stat_p_like         );
		writer.marshalDouble (stat_c_mean         );
		writer.marshalDouble (stat_c_sdev         );
		writer.marshalDouble (stat_c_like         );
	
		return;
	}

	// Marshal object.

	public static void marshal (MarshalWriter writer, RJ_Summary obj) {

		if (obj == null) {
			writer.marshalLong (MARSHAL_RJBASE_NULL);
		} else {
			obj.do_marshal (writer);
		}

		return;
	}

	// Unmarshal object.

	public static RJ_Summary unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_RJBASE_NULL, MARSHAL_RJBASE_VER);

		if (ver == MARSHAL_RJBASE_NULL) {
			return null;
		}

		return new RJ_Summary (ver, reader);
	}

	protected RJ_Summary (MarshalReader reader) {
		this (reader.unmarshalLong (MARSHAL_RJBASE_VER, MARSHAL_RJBASE_VER), reader);
	}

	private RJ_Summary (long ver, MarshalReader reader) {

		// Contents

		b                    = reader.unmarshalDouble ();
		magMain              = reader.unmarshalDouble ();
		min_a                = reader.unmarshalDouble ();
		max_a                = reader.unmarshalDouble ();
		delta_a              = reader.unmarshalDouble ();
		num_a                = reader.unmarshalInt    (1);
		min_p                = reader.unmarshalDouble ();
		max_p                = reader.unmarshalDouble ();
		delta_p              = reader.unmarshalDouble ();
		num_p                = reader.unmarshalInt    (1);
		min_c                = reader.unmarshalDouble ();
		max_c                = reader.unmarshalDouble ();
		delta_c              = reader.unmarshalDouble ();
		num_c                = reader.unmarshalInt    (1);
		apc_total_size       = reader.unmarshalInt    (1);
		apc_support_size     = reader.unmarshalInt    (0);
		apc_support_total    = reader.unmarshalDouble ();
		apc_max_tail_element = reader.unmarshalDouble ();
		a_support_lo         = reader.unmarshalInt    (0);
		a_support_hi         = reader.unmarshalInt    (0);
		p_support_lo         = reader.unmarshalInt    (0);
		p_support_hi         = reader.unmarshalInt    (0);
		c_support_lo         = reader.unmarshalInt    (0);
		c_support_hi         = reader.unmarshalInt    (0);
		stat_a_mean          = reader.unmarshalDouble ();
		stat_a_sdev          = reader.unmarshalDouble ();
		stat_a_like          = reader.unmarshalDouble ();
		stat_p_mean          = reader.unmarshalDouble ();
		stat_p_sdev          = reader.unmarshalDouble ();
		stat_p_like          = reader.unmarshalDouble ();
		stat_c_mean          = reader.unmarshalDouble ();
		stat_c_sdev          = reader.unmarshalDouble ();
		stat_c_like          = reader.unmarshalDouble ();
	}
	
}
