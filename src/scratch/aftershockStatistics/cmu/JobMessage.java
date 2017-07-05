package scratch.aftershockStatistics.cmu;

import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import com.fasterxml.jackson.databind.ser.std.ToStringSerializer;
import org.bson.types.ObjectId;

/**
 * Created by clark on 11/11/2016.
 */
public class JobMessage {
    @JsonSerialize(using = ToStringSerializer.class)
    public ObjectId id;
    public double dataMaxDays;
    public double dataMinDays;
    public double minA;
    public double maxA;
    public double minP;
    public double maxP;
    public double minC;
    public double maxC;
    public double b;
    public double g;
    public double h;
    public double magCat;

    public JobMessage(){}

    public JobMessage(ObjectId id,double dataMaxDays,double dataMinDays,double minA,double maxA,double minP,double maxP,double minC,double maxC, double b, double g, double h, double magCat){
        this.id = id;
        this.dataMaxDays=dataMaxDays;
        this.dataMinDays=dataMinDays;
        this.minA=minA;
        this.maxA=maxA;
        this.minC=minC;
        this.maxC=maxC;
        this.minP=minP;
        this.maxP=maxP;
        this.b=b;
        this.g=g;
        this.h=h;
        this.magCat=magCat;
    }
}