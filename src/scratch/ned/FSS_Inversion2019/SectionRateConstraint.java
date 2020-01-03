/**
 * 
 */
package scratch.ned.FSS_Inversion2019;

import java.util.ArrayList;


/**
 * This class represents a fault section rate constraint
 * @TODO add ability to compute confidence bounds from stdom (and vice versa), and 
 * use log-normal assumption to avoid negative rates?
 * @author Field
 *
 */
public class SectionRateConstraint implements java.io.Serializable {
	private static final long serialVersionUID = 1L;
	private String faultName; // fault name
	private int sectIndex=-1; // section index
	private double meanSectRate = Double.NaN; // mean segment rate
	private double stdDevOfMean = Double.NaN; // Standard deviation of the mean
	private double lower95Conf = Double.NaN; // Lower 95% confidence
	private double upper95Conf = Double.NaN; // Upper 95% confidence
	
	/**
	 * Constructor
	 * @param faultName
	 * @param sectIndex
	 * @param meanRate
	 * @param stdDevtoMean
	 * @param lower95Conf
	 * @param upper95Conf
	 */
	public SectionRateConstraint(String faultName, int sectIndex, double meanRate, double stdDevtoMean, double lower95Conf, double upper95Conf) {
		this.faultName = faultName;
		this.sectIndex = sectIndex;
		this.meanSectRate = meanRate;
		this.stdDevOfMean = stdDevtoMean;
		this.lower95Conf = lower95Conf;
		this.upper95Conf = upper95Conf;

	}
	
	/**
	 * Constructor
	 * @param faultName
	 * @param sectIndex
	 * @param meanRate
	 * @param stdDevtoMean
	 */
	public SectionRateConstraint(String faultName, int sectIndex, double meanRate, double stdDevtoMean) {
		this.faultName = faultName;
		this.sectIndex = sectIndex;
		this.meanSectRate = meanRate;
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
	 * Get the section index
	 * @return
	 */
	public int getSectIndex() {
		return sectIndex;
	}
	
	/**
	 * Get mean section rate
	 * @return
	 */
	public double getMeanRate() {
		return meanSectRate;
	}
	
	/**
	 * Get StdDev to mean for the rate
	 * @return
	 */
	public double getStdDevOfMean() {
		return this.stdDevOfMean;
	}
	
	/**
	 * This returns the weight average of the provided constraints (this does not check that faultNa,e and sectIndex
	 * are the save in the given constraints)
	 * @param sectRateConstraintList
	 * @return
	 */
	  public static SectionRateConstraint getWeightAveConstraint(ArrayList<SectionRateConstraint> sectRateConstraintList) {
		  double total = 0;
		  double sigmaTotal = 0;
		  String faultName=null;
		  int sectIndex = -1;
		  for(int i=0; i<sectRateConstraintList.size(); ++i) {
			  SectionRateConstraint sectRateConstraint = sectRateConstraintList.get(i);
			  faultName = sectRateConstraint.getFaultName();
			  sectIndex = sectRateConstraint.getSectIndex();
			  double sigmaSq = 1.0/(sectRateConstraint.getStdDevOfMean()*sectRateConstraint.getStdDevOfMean());
			  sigmaTotal+=sigmaSq;
			  total+=sigmaSq*sectRateConstraint.getMeanRate();
		  }
		  SectionRateConstraint finalSegRateConstraint = new SectionRateConstraint(faultName, sectIndex, total/sigmaTotal, Math.sqrt(1.0/sigmaTotal), Double.NaN, Double.NaN);
		  return finalSegRateConstraint;
	  }

	public double getLower95Conf() {
		return lower95Conf;
	}

	public double getUpper95Conf() {
		return upper95Conf;
	}

}
