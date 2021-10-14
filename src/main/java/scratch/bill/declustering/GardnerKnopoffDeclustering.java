package scratch.bill.declustering;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import java.util.ArrayList;

/**
 * Implements Gardner Knopoff Gardner-Knopoff Declustering  (Gardner and Knopoff, 1974, Bulletin of the
 * Seismological Society of America, vol. 64, pages 1363-1367) using a Matlab algorithm provided by Andreas LLenos
 * in Oct. 2021.
 *
 * @author wsavran
 *
 */
public class GardnerKnopoffDeclustering extends CatalogDeclustering {

    /**
     * Instantiates the class, but does not run the declustering method. This can be considered a 'lazy' implementation
     * because the decluster() method must be called before accessing any of the class members.
     *
     * @param inputCatalog Catalog to be declustered
     */
    public GardnerKnopoffDeclustering(ObsEqkRupList inputCatalog) {
        // initialize arrays
        this.inputCatalog = inputCatalog;
        this.ruptureIsAftershockList = new boolean[this.inputCatalog.size()];
        this.numAftershocksForRupture = new int[this.inputCatalog.size()];
        this.numAftershocksForMainshock = new ArrayList<>();
        this.aftershockCatalog = new ObsEqkRupList();
        this.mainshockCatalog = new ObsEqkRupList();
        this.isDeclustered = false;
    }

    @Override
    public ObsEqkRupList getMainshockCatalog() {
        if (!this.isDeclustered) {
            this.decluster();
        }
        return this.mainshockCatalog;
    }

    /**
     * Implements the algorithm provicded by Andrea Llenos in October 2021 from her Matlab script gkdecluster.m. This
     * implementation does not make the distinction between aftershocks and foreshocks, but simply classifies everything
     * as either a mainshock or aftershock.
     *
     * Todo: this function needs to be validated against calculations made by Andrea using her Matlab script.
     */
    @Override
    public void decluster() {
        // ned removes events below 2.5 and above 8.0. any reason why we do this in this function in particular?
        int numEvents = this.inputCatalog.size();
        int numMainshocks = 0;
        for (int i=0; i<numEvents-1; i++) {
            ObsEqkRupture parent = this.inputCatalog.get(i);
            double distanceWindow = distanceWindowInKm(parent.getMag());
            double timeWindow = timeWindowInMillis(parent.getMag());
            // search for children based on magnitude and space windows
            for (int j=i+1; j<numEvents; j++) {
                ObsEqkRupture child = this.inputCatalog.get(j);
                // catalog must be ordered in time
                if (child.getOriginTime() < parent.getOriginTime()) {
                    throw new RuntimeException("Error: catalog must be in chronological order");
                }
                // apply time window
                if (child.getOriginTime() > timeWindow) {
                    break;
                }
                // apply distance window
                double dist = distDegToKm(parent.getHypocenterLocation(), child.getHypocenterLocation());
                if (dist > distanceWindow) {
                    continue;
                }
                if (parent.getMag() >= child.getMag()) {
                    ruptureIsAftershockList[j] = true;
                    numAftershocksForRupture[j] += 1;
                } else {
                    ruptureIsAftershockList[i] = true;
                    numAftershocksForRupture[i] += 1;
                }
            }
        }
        // classify events as either mainshocks or aftershocks
        for (int i=0; i<numEvents; i++) {
            boolean isAftershock = ruptureIsAftershockList[i];
            ObsEqkRupture rupture = this.inputCatalog.get(i);
            if (isAftershock) {
                this.aftershockCatalog.add(rupture);
            } else {
                this.mainshockCatalog.add(rupture);
                this.numAftershocksForMainshock.add(this.numAftershocksForRupture[i]);
            }
        }
        this.isDeclustered = true;
    }

    @Override
    public ObsEqkRupList getAftershockCatalog() {
        if (!this.isDeclustered) {
            this.decluster();
        }
        return this.aftershockCatalog;
    }

    @Override
    public ArrayList<Integer> getNumAftershocksForMainshock() {
        if (!this.isDeclustered) {
            this.decluster();
        }
        return this.numAftershocksForMainshock;
    }

    /**
     * Implements the great circle distance calculation from Andrea LLenos provided in her Matlab script distdeg2km. This
     * function contains her implementation of the Matlab function distance.m
     * 
     * @param loc1 Location of earthquake 1
     * @param loc2 Location of earthquake 2
     * @return     Distance in kilometers between loc1 and loc2
     */
    private double distDegToKm(Location loc1, Location loc2) {
        // define useful constants
        double radius = 6378.0;
        double pi = 4.0*Math.atan(1.0);
        double deg2rad = pi / 180.0;
        // convert to radians
        double lat1 = loc1.getLatitude()*deg2rad;
        double lat2 = loc2.getLatitude()*deg2rad;
        double lon1 = loc1.getLongitude()*deg2rad;
        double lon2 = loc2.getLongitude()*deg2rad;
        // compute to distance
        double latSinTerm = Math.pow(Math.sin((lat2-lat1)/2.0), 2.0);
        double cosTerm = Math.cos(lat1)*Math.cos(lat2);
        double lonSinTerm = Math.pow(Math.sin(lon2-lon1)/2.0, 2.0);
        double a = latSinTerm + cosTerm + lonSinTerm;
        return radius*2*Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
    }

    /**
     * Distance window function provided by Andrea Llenos matlab function getWindow.m
     *
     * @param mag Magnitude of potential parent event
     * @return    Distance window in Km
     */
    private double distanceWindowInKm(double mag) {
        double exponent = 0.1238*mag+0.983;
        return Math.pow(10.0, exponent);
    }

    /**
     * Time window function provided by Andrea Lleons matlab function getWindow.m
     *
     * @param mag Magtnidue of potential parent event
     * @return    Time window in milliseconds
     */
    private long timeWindowInMillis(double mag) {
        double exponent;
        long daysToMillis = 24*60*60*1000;
        if (mag >= 6.5) {
            exponent = 0.06166*mag+2.5117;
        } else {
            exponent = 0.556*mag-0.6027;
        }
        return (long) Math.pow(10.0, exponent)*daysToMillis;
    }
}
