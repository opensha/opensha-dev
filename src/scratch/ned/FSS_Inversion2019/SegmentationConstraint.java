/**
 * 
 */
package scratch.ned.FSS_Inversion2019;

/**
 * 
 * This class specifies the rate at which two fault sections co-rupture.  e.g., a strict segmentation constraint
 * exists between two neighboring sections that have a co-rupture rate of zero.
 * @author Field
 *
 */
public class SegmentationConstraint implements java.io.Serializable {
	private static final long serialVersionUID = 1L;
	private String faultName; // fault name
	private int sect1_Index=-1; // section index
	private int sect2_Index=-1; // section index
	private double meanCoRuptureRate = Double.NaN; // mean segment rate
	private double stdDevOfMean = Double.NaN; // Standard deviation of the mean
	private double lower95Conf = Double.NaN; // Lower 95% confidence
	private double upper95Conf = Double.NaN; // Upper 95% confidence
	
	/**
	 * Constructor
	 * @param faultName
	 * @param sect1_Index
	 * @param sect2_Index
	 * @param meanCoRuptureRate
	 * @param stdDevtoMean
	 * @param lower95Conf
	 * @param upper95Conf
	 */
	public SegmentationConstraint(String name, int sect1_Index, int sect2_Index, double meanCoRuptureRate, 
			double stdDevtoMean, double lower95Conf, double upper95Conf) {
		this.faultName = name;
		this.sect1_Index = sect1_Index;
		this.sect2_Index = sect2_Index;
		this.meanCoRuptureRate = meanCoRuptureRate;
		this.stdDevOfMean = stdDevtoMean;
		this.lower95Conf = lower95Conf;
		this.upper95Conf = upper95Conf;

	}
	
	/**
	 * Constructor
	 * @param faultName
	 * @param sect1_Index
	 * @param sect2_Index
	 * @param meanCoRuptureRate
	 * @param stdDevtoMean
	 */
	public SegmentationConstraint(String name, int sect1_Index, int sect2_Index, double meanCoRuptureRate, 
			double stdDevtoMean) {
		this.faultName = name;
		this.sect1_Index = sect1_Index;
		this.sect2_Index = sect2_Index;
		this.meanCoRuptureRate = meanCoRuptureRate;
		this.stdDevOfMean = stdDevtoMean;

	}

	
	/**
	 * Get the fault name
	 * @return
	 */
	public String getFaultName() {
		return faultName;
	}
	
	/**
	 * Get the 1st section index
	 * @return
	 */
	public int getSect1_Index() {
		return sect1_Index;
	}
	
	/**
	 * Get the 2nd section index
	 * @return
	 */
	public int getSect2_Index() {
		return sect2_Index;
	}
	
	/**
	 * Get mean section rate
	 * @return
	 */
	public double getMeanCoRuptureRate() {
		return meanCoRuptureRate;
	}
	
	/**
	 * Get StdDev to mean for the rate
	 * @return
	 */
	public double getStdDevOfMean() {
		return this.stdDevOfMean;
	}
	
	public double getLower95Conf() {
		return lower95Conf;
	}

	public double getUpper95Conf() {
		return upper95Conf;
	}

}
