package scratch.aftershockStatistics;

public class MagCompPage_Parameters {
	
	private double magCat;
	private double capG;
	private double capH;
	
	/**
	 * This class is a container for the magnitude of completeness parameters defined by 
	 * Page et al. (2016, BSSA).
	 * According to Page et al. the magnitude of completeness is
	 *  magMin(t) = Max(magMain/2 - G - H*log10(t), magCat)
	 * where t is measured in days.
	 * As a special case, if capG == 10.0 then the magnitude of completeness is always magCat
	 * (in which case it is recommended that capH == 0.0).
	 * 
	 * @param magCat
	 * @param capG
	 * @param capH
	 */
	public MagCompPage_Parameters(double magCat, double capG, double capH) {
		this.magCat = magCat;
		this.capG = capG;
		this.capH = capH;
	}

	public MagCompPage_Parameters(){}
	
	/**
	 * This returns the catalog magnitude of completeness (magCat).
	 * @return
	 */
	public double get_magCat() {return magCat;}
	
	/**
	 * This returns the G parameter (capG).
	 * @return
	 */
	public double get_capG() {return capG;}
	
	/**
	 * This returns the H parameter (capH).
	 * @return
	 */
	public double get_capH() {return capH;}

	@Override
	public String toString() {
		return "Page_Params[magCat="+get_magCat()+", capG="+get_capG()+", capH="+get_capH()+"]";
	}




	//----- Marshaling -----

	// Marshal version number.

	public static final long MARSHAL_VER = 2001L;

	// Marshal object.

	public void marshal (MarshalWriter writer) {

		// Version

		writer.marshalLong (MARSHAL_VER);

		// Contents

		writer.marshalDouble (magCat);
		writer.marshalDouble (capG  );
		writer.marshalDouble (capH  );
	
		return;
	}

	// Unmarshal object.

	public static MagCompPage_Parameters unmarshal (MarshalReader reader) {
	
		// Version

		long ver = reader.unmarshalLong (MARSHAL_VER, MARSHAL_VER);

		return new MagCompPage_Parameters (ver, reader);
	}

	private MagCompPage_Parameters (long ver, MarshalReader reader) {

		// Contents

		magCat = reader.unmarshalDouble();
		capG   = reader.unmarshalDouble();
		capH   = reader.unmarshalDouble();
	}
	
}
