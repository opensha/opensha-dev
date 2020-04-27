/**
 * 
 */
package scratch.ned.FSS_Inversion2019;

/**
 * 
 * This class specifies the rate or slip-rate at which two fault sections co-rupture.  
 * e.g., a strict segmentation constraint exists between two neighboring sections that 
 * have a co-rupture rate of zero.
 * @author Field
 *
 */
public class SegmentationConstraint implements java.io.Serializable {
	private static final long serialVersionUID = 1L;
	private String faultName; // fault name
	private boolean isSlipRateConstraint = false;
	private int sect1_Index=-1; // section index
	private int sect2_Index=-1; // section index
	private double meanJointRate = Double.NaN; // mean segment rate
	private double stdDevOfMean = Double.NaN; // Standard deviation of the mean
	private double lower95Conf = Double.NaN; // Lower 95% confidence
	private double upper95Conf = Double.NaN; // Upper 95% confidence
	
	/**
	 * Constructor - this assumes it's an event-rate constraint (not a slip-rate constraint)
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
		this.meanJointRate = meanCoRuptureRate;
		this.stdDevOfMean = stdDevtoMean;
		this.lower95Conf = lower95Conf;
		this.upper95Conf = upper95Conf;

	}
	
	/**
	 * Constructor - this assumes it's an event-rate constraint (not a slip-rate constraint)
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
		this.meanJointRate = meanCoRuptureRate;
		this.stdDevOfMean = stdDevtoMean;

	}

	
	/**
	 * Constructor - general case
	 * @param faultName
	 * @param sect1_Index
	 * @param sect2_Index
	 * @param meanCoRuptureRate
	 * @param isSlipRateConstraint - indicates whether this is a slip- are event-rate constraint
	 * @param stdDevtoMean
	 * @param lower95Conf
	 * @param upper95Conf
	 */
	public SegmentationConstraint(String name, int sect1_Index, int sect2_Index, double meanJointRate, 
			double stdDevtoMean, double lower95Conf, double upper95Conf, boolean isSlipRateConstraint) {
		this.faultName = name;
		this.sect1_Index = sect1_Index;
		this.sect2_Index = sect2_Index;
		this.meanJointRate = meanJointRate;
		this.stdDevOfMean = stdDevtoMean;
		this.lower95Conf = lower95Conf;
		this.upper95Conf = upper95Conf;
		this.isSlipRateConstraint = isSlipRateConstraint;
	}
	
	/**
	 * Constructor - general case
	 * @param faultName
	 * @param sect1_Index
	 * @param sect2_Index
	 * @param meanCoRuptureRate
	 * @param isSlipRateConstraint - indicates whether this is a slip- are event-rate constraint
	 * @param stdDevtoMean
	 */
	public SegmentationConstraint(String name, int sect1_Index, int sect2_Index, double meanJointRate, 
			double stdDevtoMean, boolean isSlipRateConstraint) {
		this.faultName = name;
		this.sect1_Index = sect1_Index;
		this.sect2_Index = sect2_Index;
		this.meanJointRate = meanJointRate;
		this.stdDevOfMean = stdDevtoMean;
		this.isSlipRateConstraint = isSlipRateConstraint;
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
	 * Get mean joint rate (co-rupture or slip-rate at boundary)
	 * @return
	 */
	public double getMeanJointRate() {
		return meanJointRate;
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
	
	/**
	 * @return true if it's an event-rate constraint and false if it's a slip-rate constraint
	 */
	public boolean isSlipRateConstraint() {
		return isSlipRateConstraint;
	}

}
