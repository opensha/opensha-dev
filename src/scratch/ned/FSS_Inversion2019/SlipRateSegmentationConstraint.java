/**
 * 
 */
package scratch.ned.FSS_Inversion2019;

/**
 * 
 * This class constitutes a segmentation constraint by specifying a slip rate (and stdev) reduction for
 * a fault section.
 * 
 * @author Field
 *
 */
public class SlipRateSegmentationConstraint implements java.io.Serializable {
	private static final long serialVersionUID = 1L;
	private String name; // fault name
	private int sectIndex=-1; // section index
	private double slipRateReductionFactor = Double.NaN; // slip rate reduction factor
	private double slipRateStdevReductionFactor = Double.NaN; // slip rate stdev reduction factor
	
	/*
	 * Constructor
	 */
	public SlipRateSegmentationConstraint(String name, int sectIndex, double slipRateReductionFactor, 
			double slipRateStdevReductionFactor) {
		this.name = name;
		this.sectIndex = sectIndex;
		this.slipRateReductionFactor = slipRateReductionFactor;
		this.slipRateStdevReductionFactor = slipRateStdevReductionFactor;
	}
	
	
	/**
	 * Get the fault name
	 * @return
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * Get the section index
	 * @return
	 */
	public int getSectIndex() {
		return sectIndex;
	}
	
	
	/**
	 * Get section slipRateReductionFactor
	 * @return
	 */
	public double getSlipRateReductionFactor() {
		return slipRateReductionFactor;
	}
	
	/**
	 * Get section slipRateStdevReductionFactor
	 * @return
	 */
	public double getSlipRateStdevReductionFactor() {
		return slipRateStdevReductionFactor;
	}
}
