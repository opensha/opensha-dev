package scratch.aftershockStatistics.cmu;

import java.util.Collections;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupOrigTimeComparator;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.util.EqkRuptureMagComparator;

import org.bson.types.ObjectId;

import org.mongodb.morphia.annotations.Entity;
import org.mongodb.morphia.annotations.*;




@Entity("Aftershock")
/**
 * <p>Title: ObsEqkSequence </p>
 *
 * <p>Description: Modified implementation of ObsEqkRuptureList</p>
 *
 * <p>Copyright: Copyright (c) 2016</p>
 *
 * <p>Company: CMU</p>
 *
 * @author rewritten by Daniel Chen
 * @version 1.0
 */
public class ObsEqkSequence implements java.io.Serializable{

    private String mainShockId;
    private String foreShockId;
    private RegionParams mainShockRegion;
    private RegionParams foreShockRegion;
    private boolean isActive;
    private ObsEqkRupList eqkList;
    private long sampleTime;

    @Id
    private ObjectId id;

    public ObsEqkSequence(){}

    public ObsEqkSequence(String mainShockId){
        this.mainShockId = mainShockId;
        this.foreShockId = null;
        this.isActive=false;
        this.eqkList=new ObsEqkRupList();
    }


    public ObjectId getId() {
        return id;
    }

    public void setId(ObjectId id) {
        this.id = id;
    }

    public void setMainShockRegion(RegionParams region)
    {
        this.mainShockRegion=region;
    }
    public void setForeShockRegion(RegionParams region)
    {
        this.foreShockRegion=region;
    }
    public RegionParams getMainShockRegion()
    {
        return this.mainShockRegion;
    }

    public RegionParams getForeShockRegion()
    {
        return this.foreShockRegion;
    }

    public ObsEqkRupture getMainShock() {
        for(ObsEqkRupture eqk : this.eqkList){
            if (eqk.getEventId().equals(mainShockId)){
                return eqk;
            }
        }
        return null;
    }
    /**
     * Sets the mainshock ruptures.
     * @param mainShockId
     * @return void
     */
    public void setMainShockId(String mainShockId) {
        this.mainShockId = mainShockId;
    }

    public ObsEqkRupture getForeShock() {
        if(foreShockId!=null){
            for(ObsEqkRupture eqk : this.eqkList){
                if (eqk.getEventId().equals(foreShockId)){
                    return eqk;
                }
            }
        }
        return null;
    }

    public void setForeShockId(String foreShockId) {
        this.foreShockId = foreShockId;
    }

    /**
     * Sets the list of aftershock ruptures.
     * @param list
     * @return void
     */
    public void setEqkList(ObsEqkRupList list) {
        this.eqkList = list;
    }

    public ObsEqkRupList getEqkList()
    {
        return this.eqkList;
    }

    public long getSampleTime()
    {
        return this.sampleTime;
    }

    public void setSampleTime(long t)
    {
        this.sampleTime=t;
    }

    public boolean getisActive()
    {
        return this.isActive;
    }

    public void setisActive(boolean active)
    {
        this.isActive=active;
    }


    /**
     * Returns the list of aftershock ruptures.
     * @return the list of aftershock ruptures.
     */
    public ObsEqkRupList getAfterShocks() {
        for(ObsEqkRupture eqk : this.eqkList){
            if (eqk.getEventId().equals(mainShockId)){
                return (ObsEqkRupList) this.getRupsAfter(eqk.getOriginTime());
            }
        }
        return null;
    }

    /**
     * Returns the list of aftershock ruptures.
     * @return the list of aftershock ruptures.
     */
    public ObsEqkRupList getForeShocks() {
        ObsEqkRupList list=null;
        Long foreShockOriginTime = 0L;
        Long mainShockOriginTime = 0L;
        if(foreShockId!=null){
            for(ObsEqkRupture eqk : this.eqkList){
                if (eqk.getEventId().equals(foreShockId)){
                    list = new ObsEqkRupList();
                    list.add(eqk);
                    foreShockOriginTime=eqk.getOriginTime();

                }
                else if(eqk.getEventId().equals(mainShockId))
                {
                    mainShockOriginTime=eqk.getOriginTime();
                    break;
                }
            }
            list.addAll(this.getRupsBetween(foreShockOriginTime, mainShockOriginTime));
            return list;
        }
        return null;
    }


    /**
     * Returns the list of the Observed events above/at the given magnitude.
     * @param mag double Magnitude
     * @return the subset of total observed events as ObsEqkRupList list
     * above a given magnitude
     */
    public ObsEqkRupList getRupsAboveMag(double mag) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            if (eqkRup.getMag() >= mag)
                obsEventList.add(eqkRup);
        }
        return obsEventList;
    }

    /**
     * Returns the list of the Observed events below the given magnitude.
     * @param mag double Magnitude
     * @return the subset of total observed events as ObsEqkRupList list
     * below a given magnitude
     */
    public ObsEqkRupList getRupsBelowMag(double mag) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            if (eqkRup.getMag() < mag)
                obsEventList.add(eqkRup);

        }
        return obsEventList;
    }

    /**
     * Returns the list of the Observed events between 2 given magnitudes.
     * It includes lower magnitude in the range but excludes the upper magnitude.
     * @param mag1 double lower magnitude
     * @param mag2 double upper magnitude
     * @return the subset of total observed events as ObsEqkRupList list
     * between 2 given magnitudes.
     */
    public ObsEqkRupList getRupsBetweenMag(double mag1, double mag2) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            double eventMag = eqkRup.getMag();
            if (eventMag >= mag1 && eventMag < mag2)
                obsEventList.add(eqkRup);
        }
        return obsEventList;

    }

    /**
     * Returns the list of Observed events before a given time in milliseconds (epoch)
     * @param timeInMillis - what returned by GregorianCalendar.getTimeInMillis()
     * @return the subset of total observed events as ObsEqkRupList list
     * before a given time period
     */
    public ObsEqkRupList getRupsBefore(long timeInMillis) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            long eventTime = eqkRup.getOriginTime();
            if (eventTime < timeInMillis)
                obsEventList.add(eqkRup);
        }
        return obsEventList;
    }



    /**
     * Returns the list of Observed events after a given time in milliseconds (epoch)
     * @param timeInMillis - what returned by GregorianCalendar.getTimeInMillis()
     * @return the subset of total observed events as ObsEqkRupList list
     * after a given time period
     */
    public ObsEqkRupList getRupsAfter(long timeInMillis) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            long eventTime = eqkRup.getOriginTime();
            if (eventTime > timeInMillis)
                obsEventList.add(eqkRup);
        }
        return obsEventList;

    }

    /**
     * Returns the list of the Observed events between 2 given time periods.
     * @param timeInMillis1 GregorianCalendar Time Period
     * @param timeInMillis2 GregorianCalendar Time Period
     * @return the subset of total observed events as ObsEqkRupList list
     * between 2 given time periods.
     */
    public ObsEqkRupList getRupsBetween(long timeInMillis1,
                                        long timeInMillis2) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            long eventTime = eqkRup.getOriginTime();
            if (eventTime > timeInMillis1 && eventTime < timeInMillis2)
                obsEventList.add(eqkRup);
        }
        return obsEventList;

    }

    /**
     * Returns the list of the Observed events inside a given geographic region
     * @param region Region
     * @return the subset of total observed events as ObsEqkRupList list
     * inside a given region.
     */
    public ObsEqkRupList getRupsInside(Region region) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            Location loc = eqkRup.getHypocenterLocation();
            if(region.contains(loc))
                obsEventList.add(eqkRup);
        }
        return obsEventList;

    }

    /**
     * Returns the list of the Observed events outside a given geographic region
     * @param region Region
     * @return the subset of total observed events as ObsEqkRupList list
     * outside a given region.
     */
    public ObsEqkRupList getRupsOutside(Region region) {
        ObsEqkRupList obsEventList = new ObsEqkRupList();
        for (ObsEqkRupture eqkRup : this.eqkList) {
            Location loc = eqkRup.getHypocenterLocation();
            if (!region.contains(loc))
                obsEventList.add(eqkRup);
        }
        return obsEventList;
    }


    /**
     * Sorts the Observed Eqk Rupture Event list based on the magitude.
     *
     */
    public void sortByMag(){
        Collections.sort(this.eqkList, new EqkRuptureMagComparator());
    }

    /**
     * Sorts the Observed Eqk Rupture Event list based on the Origin time.
     *
     */
    public void sortByOriginTime() {
        Collections.sort(this.eqkList, new ObsEqkRupOrigTimeComparator());
    }

}