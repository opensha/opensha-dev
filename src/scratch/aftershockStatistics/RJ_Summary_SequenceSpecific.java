package scratch.aftershockStatistics;

public class RJ_Summary_SequenceSpecific extends RJ_Summary {

	// Summary values, see RJ_AftershockModel_SequenceSpecific for description.

	private double capG                ;
	private double capH                ;
	private double magCat              ;
	private double dataStartTimeDays   ;
	private double dataEndTimeDays     ;
	private int    numAftershocks      ;

	
	/**
	 * Default constructor.
	 */
	public RJ_Summary_SequenceSpecific() {
		super();
	}

	
	/**
	 * Construct from R&J model.
	 */
	public RJ_Summary_SequenceSpecific(RJ_AftershockModel_SequenceSpecific model) {
		super (model);

		this.capG                 = model.capG                 ;
		this.capH                 = model.capH                 ;
		this.magCat               = model.magCat               ;
		this.dataStartTimeDays    = model.dataStartTimeDays    ;
		this.dataEndTimeDays      = model.dataEndTimeDays      ;
		this.numAftershocks       = model.numAftershocks       ;
	}

	
	/**
	 * Getters.
	 */
	public double get_capG                () {return capG                ;}
	public double get_capH                () {return capH                ;}
	public double get_magCat              () {return magCat              ;}
	public double get_dataStartTimeDays   () {return dataStartTimeDays   ;}
	public double get_dataEndTimeDays     () {return dataEndTimeDays     ;}
	public int    get_numAftershocks      () {return numAftershocks      ;}



	@Override
	public String toString() {
		return super.toString() + "\n" +
			"capG                 = " + get_capG                () + "\n" +
			"capH                 = " + get_capH                () + "\n" +
			"magCat               = " + get_magCat              () + "\n" +
			"dataStartTimeDays    = " + get_dataStartTimeDays   () + "\n" +
			"dataEndTimeDays      = " + get_dataEndTimeDays     () + "\n" +
			"numAftershocks       = " + get_numAftershocks      () ;
	}




	//----- Marshaling -----

	// Marshal version number.

	public static final long MARSHAL_RJSEQ_NULL = 7000L;
	public static final long MARSHAL_RJSEQ_VER = 7001L;

	// Marshal object.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_RJSEQ_VER);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalDouble (capG                );
		writer.marshalDouble (capH                );
		writer.marshalDouble (magCat              );
		writer.marshalDouble (dataStartTimeDays   );
		writer.marshalDouble (dataEndTimeDays     );
		writer.marshalInt    (numAftershocks      );
	
		return;
	}

	// Marshal object.

	public static void marshal (MarshalWriter writer, RJ_Summary_SequenceSpecific obj) {

		if (obj == null) {
			writer.marshalLong (MARSHAL_RJSEQ_NULL);
		} else {
			obj.do_marshal (writer);
		}

		return;
	}

	// Unmarshal object.

	public static RJ_Summary_SequenceSpecific unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_RJSEQ_NULL, MARSHAL_RJSEQ_VER);

		if (ver == MARSHAL_RJSEQ_NULL) {
			return null;
		}

		return new RJ_Summary_SequenceSpecific (ver, reader);
	}

	protected RJ_Summary_SequenceSpecific (MarshalReader reader) {
		this (reader.unmarshalLong (MARSHAL_RJSEQ_VER, MARSHAL_RJSEQ_VER), reader);
	}

	private RJ_Summary_SequenceSpecific (long ver, MarshalReader reader) {
		super (reader);

		// Contents

		capG                 = reader.unmarshalDouble ();
		capH                 = reader.unmarshalDouble ();
		magCat               = reader.unmarshalDouble ();
		dataStartTimeDays    = reader.unmarshalDouble ();
		dataEndTimeDays      = reader.unmarshalDouble ();
		numAftershocks       = reader.unmarshalInt    (0);
	}
	
}
