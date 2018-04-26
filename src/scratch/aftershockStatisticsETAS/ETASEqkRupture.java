package scratch.aftershockStatisticsETAS;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;

/**
 * <p>Title: ObsEqkRupture </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author rewritten by Ned Field and Nicholas vdE
 * @version 1.0
 */
public class ETASEqkRupture extends ObsEqkRupture implements java.io.Serializable{

	
	private static final long serialVersionUID = 1L;
	
	protected double spatialKernelDistance;
	protected double stressDrop;
	protected FaultTrace faultTrace;
	protected String name;
	protected String description;
	
	public ETASEqkRupture(){}

	public ETASEqkRupture(ObsEqkRupture mainshock, double stressDrop){
		
		this.setHypocenterLocation(mainshock.getHypocenterLocation());
		this.setMag(mainshock.getMag());
		this.setEventId(mainshock.getEventId());
		this.setOriginTime(mainshock.getOriginTime());
		this.setStressDrop(stressDrop);
		this.setSpatialKernelDistance(ETAS_StatsCalc.magnitude2radius(mainshock.getMag(), stressDrop));
		faultTrace = new FaultTrace("Mainshock");
		faultTrace.add(getHypocenterLocation());
		this.setFaultTrace(faultTrace);
	}
	
	
	/**
	 * This constructor sets the rupture surface as a point source at the given hypocenter
	 * @param eventId
	 * @param originTimeInMillis
	 * @param hypoLoc
	 * @param mag
	 */
	public ETASEqkRupture(String eventId, long originTimeInMillis, 
			Location hypoLoc, double mag, double stressDrop) {
		super(eventId, originTimeInMillis, hypoLoc, mag);
		
		setSpatialKernelDistance(ETAS_StatsCalc.magnitude2radius(mag, stressDrop));
	}

	
	public void setFaultTrace(FaultTrace faultTrace){
		this.faultTrace = faultTrace;

	}
	
	public FaultTrace getFaultTrace(){
		return faultTrace;

	}

	public void setStressDrop(double stressDrop){
		this.stressDrop = stressDrop;
	}
	
	public double getStressDrop(){
		return stressDrop;
	}
	
	public double getSpatialKernelDistance(){
		return spatialKernelDistance;
	}
	
	public void setSpatialKernelDistance(double mag, double stressDrop){
		setSpatialKernelDistance(ETAS_StatsCalc.magnitude2radius(mag, stressDrop));
	}

	public void setSpatialKernelDistance(double distance){
		this.spatialKernelDistance = distance;
	}

	/**
	 * Returns an info String for the ETAS Observed EqkRupture
	 * @return String
	 */
	public String getInfo(){
		String obsEqkInfo = super.getInfo();
		obsEqkInfo += "EventId ="+eventId+"\n";
		obsEqkInfo += "OriginTimeInMillis ="+originTimeInMillis+"\n";
		obsEqkInfo += "SpatialKernelDistance ="+spatialKernelDistance+"\n";
		return obsEqkInfo;
	}


	/**
	 * Clones the eqk rupture and returns the new cloned object
	 * @return
	 */
	public Object clone() {
		ETASEqkRupture eqkEventClone = new ETASEqkRupture();
		eqkEventClone.setEventId(eventId);
		eqkEventClone.setMag(mag);
		eqkEventClone.setRuptureSurface(getRuptureSurface());
		eqkEventClone.setHypocenterLocation(hypocenterLocation);
		eqkEventClone.setOriginTime(originTimeInMillis);
		eqkEventClone.setAveRake(aveRake);
		eqkEventClone.setSpatialKernelDistance(spatialKernelDistance);
		return eqkEventClone;
	}

	public void setName(String string) {
		this.name = string;
	}

	public String getName() {
		return name;
	}

	public void setDescription(String string) {
		this.description = string;
	}

	public String getDescription() {
		return description;
	}

}
