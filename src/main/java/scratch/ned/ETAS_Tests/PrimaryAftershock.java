package scratch.ned.ETAS_Tests;

import org.opensha.commons.exceptions.InvalidRangeException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.faultSurface.AbstractEvenlyGriddedSurface;


/**
 *
 * <b>Title:</b> EqkRupture<br>
 * <b>Description:</b> <br>
 *
 * @author Ned Field
 * @version 1.0
 */

public class PrimaryAftershock extends EqkRupture {
	
	private int id, parentID, srcIndex, rupIndex, generation;
	private double originTime, distToParent=Double.NaN;
	
	public PrimaryAftershock() {};
	
	public PrimaryAftershock(ProbEqkRupture probRup) {
		this.setMag(probRup.getMag());
		this.setRuptureSurface(probRup.getRuptureSurface());
		this.setAveRake(probRup.getAveRake());
	}
	
	/**
	 * The ID of the parent that spawned this primary aftershock
	 * @return
	 */
	public int getParentID() {
		return parentID;
	}
	
	/**
	 * Sets the ID of the parent that spawned this primary aftershock
	 * @return
	 */
	public void setParentID(int parentID) {
		this.parentID = parentID;
	}
	
	/**
	 * This returns the distance to the parent (NaN if never set)
	 * @return
	 */
	public double getDistanceToParent() {
		return distToParent;
	}
	
	/**
	 * This sets the distance to parent
	 * @param distToParent
	 */
	public void setDistanceToParent(double distToParent) {
		this.distToParent = distToParent;
	}
	
	/**
	 * The ID of this event
	 * @return
	 */
	public int getID() {
		return id;
	}

	
	/**
	 * Sets the ID of this event
	 * @return
	 */
	public void setID(int id) {
		this.id = id;
	}

	/**
	 * The ID of this event
	 * @return
	 */
	public int getGeneration() {
		return generation;
	}

	
	/**
	 * Sets the ID of this event
	 * @return
	 */
	public void setGeneration(int generation) {
		this.generation = generation;
	}


	
	/**
	 * The origin time of this event (units are unspecified)
	 * @return
	 */
	public double getOriginTime(){
		return originTime;
	}
	
	/**
	 * Sets the origin time of this event (units are unspecified)
	 * @return
	 */
	public void setOriginTime(double originTime){
		this.originTime = originTime;
	}

	public void setERF_SourceIndex(int srcIndex){
		this.srcIndex = srcIndex;
	}

	public int getERF_SourceIndex(){
		return srcIndex;
	}
	
	public void setERF_RupIndex(int rupIndex){
		this.rupIndex = rupIndex;
	}

	public int getERF_RupIndex(){
		return rupIndex;
	}
}
