package scratch.aftershockStatistics.cmu;

public class LocationWrapper {

    public double lat;
    public double lon;
    public double depth;

    public LocationWrapper(){}

    public LocationWrapper(double lat, double lon, double depth)
    {
        this.lat=lat;
        this.lon=lon;
        this.depth=depth;
    }

    public void setLatitude(double lat)
    {
        this.lat=lat;
    }

    public void setLongitude(double lon)
    {
        this.lon=lon;
    }


    public void setDepth(double depth)
    {
        this.depth=depth;
    }


}