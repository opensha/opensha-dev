package scratch.aftershockStatistics;

public class MagCompPage_Parameters {
	
	double magCat;
	double capG;
	double capH;
	
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
	
}
