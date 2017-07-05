package scratch.aftershockStatistics.cmu;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import scratch.aftershockStatistics.AftershockStatsCalc;
import scratch.aftershockStatistics.ComcatAccessor;


public class RegionParams {
	/** define various region types and center types*/
    public static String CIRCULAR="Circular";
    public static String WCCIRCULAR="WC Circular 1994";
    public static String RECTANGULAR="Rectangular";
    public static String EPICENTER="Epicenter";
    public static String CUSTOM ="Custom Location";
    public static String CENTROID="Centroid";
    private String regionType;
    private String centerType;
    private double radius;
    private Location circleCenter;
    private Location minLocation;
    private Location maxLocation;

    public RegionParams(){}

    public RegionParams(ComcatMessage message)
    {
        this.regionType=message.regionType;
        this.centerType=message.centerType;
        this.radius=message.radius;
        this.circleCenter=null;
        this.minLocation=new Location(message.minLocation.lat,message.minLocation.lon,message.minLocation.depth);
        this.maxLocation=new Location(message.maxLocation.lat,message.maxLocation.lon,message.maxLocation.depth);
    }




    public void setRegionType(String regionType)
    {
        this.regionType=regionType;
    }

    public String getRegionType()
    {
        return regionType;
    }

    public void setCenterType(String centerType)
    {

        this.centerType=centerType;
    }

    public String getCenterType()
    {
        return this.centerType;
    }

    public void setRadius(double radius)
    {
        this.radius=radius;
    }
    public double getRadius()
    {
        return this.radius;
    }

    public void setCicleCenter(Location circleCenter)
    {
        this.circleCenter=circleCenter;
    }

    public Location getCicleCenter()
    {
        return this.circleCenter;
    }
    public void setminLoc(Location minLoc)
    {
        this.minLocation=minLoc;
    }

    public Location getminLoc()
    {
        return this.minLocation;
    }

    public void setmaxLoc(Location maxLoc)
    {
        this.maxLocation=maxLoc;
    }

    public Location getmaxLoc()
    {
        return this.maxLocation;
    }

    /**
     * returns a region circular or rectangular depending on region type
     * @param void
     * @return void
     * */
    public Region buildRegion()
    {
        if(this.regionType.equals("Rectangular"))
        {
            return new Region(this.minLocation,this.maxLocation);
        }
        else
        {
            return new Region(this.circleCenter,this.radius);
        }
    }
    
    /**
     * returns a circular region based on centroid calculation
     * @param mainshock, midays, maxDays
     * @return Region to examine for forecast
     * */

    public Region buildRegion(ObsEqkRupture mainShock, double dataMinDays, double dataMaxDays)
    {
        double minDepth = this.minLocation.getDepth();
        double maxDepth = this.maxLocation.getDepth();
        Region region = new Region(this.circleCenter,this.radius);
        ComcatAccessor accessor = new ComcatAccessor();
        ObsEqkRupList afs = accessor.fetchAftershocks(mainShock, dataMinDays, dataMaxDays, minDepth, maxDepth, region);
        if (!afs.isEmpty())
        {
            Location centroid = AftershockStatsCalc.getCentroid(mainShock, afs);
            this.circleCenter=centroid;
            if(this.regionType.equals(WCCIRCULAR))
            {
                WC1994_MagLengthRelationship wcMagLen = new WC1994_MagLengthRelationship();
                this.setRadius(wcMagLen.getMedianLength(mainShock.getMag()));
            }
            region = new Region(centroid, this.radius);
        }
        return region;

    }
    /**
     * tells whether the region under examination is circular or no
     * @param void
     * @return void
     * */

    public boolean isCircular()
    {
        if(this.regionType.equals(CIRCULAR)|| this.regionType.equals(WCCIRCULAR))
            return true;
        return false;
    }
}