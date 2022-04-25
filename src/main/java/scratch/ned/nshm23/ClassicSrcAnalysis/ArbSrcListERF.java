package scratch.ned.nshm23.ClassicSrcAnalysis;

import java.util.ArrayList;

import org.opensha.commons.data.TimeSpan;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;


/**
 * <p>Title: ArbSrcListERF</p>
 * <p>Description: 
 * </UL><p>
 * @author Ned Field
 * Date : March, 2022
 * @version 1.0
 */

public class ArbSrcListERF extends AbstractERF{

  //for Debug purposes
  private static String  C = "ArbSrcListERF";
  private boolean D = false;

  //name for this classs
  public final static String  NAME = "ArbSrcListERF";

  // this is the source (only 1 for this ERF)
  private ArrayList<ProbEqkSource> sourceList;


  /**
   * Constructor for this source
   */
  public ArbSrcListERF(ArrayList<ProbEqkSource> sourceList, double duration) {
	  timeSpan = new TimeSpan(TimeSpan.NONE,TimeSpan.YEARS);
	  timeSpan.setDuration(duration);
	  this.sourceList = sourceList;
  }



   /**
    * update the source based on the paramters (only if a parameter value has changed)
    */
   public void updateForecast(){
     String S = C + "updateForecast::";

     if(parameterChangeFlag) {
    	 throw new RuntimeException("this class has no adjustable prameters; problem with timeSpan?");
 //      parameterChangeFlag = false;
     }

   }


   /**
    * Return the earhthquake source at index i.   Note that this returns a
    * pointer to the source held internally, so that if any parameters
    * are changed, and this method is called again, the source obtained
    * by any previous call to this method will no longer be valid.
    *
    * @param iSource : index of the desired source (only "0" allowed here).
    *
    * @return Returns the ProbEqkSource at index i
    *
    */
   public ProbEqkSource getSource(int iSource) {

    return sourceList.get(iSource);
   }


   /**
    * Returns the number of earthquake sources (always "1" here)
    *
    * @return integer value specifying the number of earthquake sources
    */
   public int getNumSources(){
     return sourceList.size();
   }


    /**
     *  This returns a list of sources (contains only one here)
     *
     * @return ArrayList of Prob Earthquake sources
     */
    public ArrayList  getSourceList(){
      return sourceList;
    }


  /**
   * Return the name for this class
   *
   * @return : return the name for this class
   */
   public String getName(){
     return NAME;
   }
   
//   /**
//    * This overides the parent method to ignore whatever is passed in
//    * (because timeSpan is always null in this class)
//    * @param timeSpan : TimeSpan object
//    */
//   public void setTimeSpan(TimeSpan time) {
//     // do nothing
//   }


}
