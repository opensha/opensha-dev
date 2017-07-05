package scratch.aftershockStatistics.cmu;


/**
 * Created by clark on 11/11/2016.
 */

/**
 * parse parameters received in the message from activemq queue
 * */

public class ComcatMessage {
    public String eventId;
    public String regionType;
    public String centerType;
    public double radius;
    public LocationWrapper circleCenter;
    public LocationWrapper minLocation;
    public LocationWrapper maxLocation;
    public boolean persistent;
    public boolean emulate;
    public double dataMaxDays;
    public double dataMinDays;
    public double minA;
    public double maxA;
    public double minP;
    public double maxP;
    public double minC;
    public double maxC;
    public double b;
    public double h;
    public double g;
    public double magCat;

    public ComcatMessage(){}
    public ComcatMessage(String eventId, double dataMaxDays,double dataMinDays, String regionType,String centerType, LocationWrapper circleCenter,double radius, LocationWrapper minLoc, LocationWrapper maxLoc, boolean isPersistent, boolean emulate, double minA, double maxA, double minP, double maxP, double minC, double maxC, double b, double g, double h, double magCat ){
        this.eventId = eventId;
        this.regionType=regionType;
        this.centerType=centerType;
        this.circleCenter=circleCenter;
        this.radius=radius;
        this.minLocation=minLoc;
        this.maxLocation=maxLoc;
        this.persistent=isPersistent;
        this.emulate=emulate;
        this.dataMaxDays=dataMaxDays;
        this.dataMinDays=dataMinDays;
        this.minA=minA;
        this.maxA=maxA;
        this.minP=minP;
        this.maxP=maxP;
        this.minC=minC;
        this.maxC=maxC;
        this.b = b;
        this.h = h;
        this.g = g;
        this.magCat = magCat;
    }


}