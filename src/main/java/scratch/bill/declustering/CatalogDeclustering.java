package scratch.bill.declustering;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;

import java.util.ArrayList;

public abstract class CatalogDeclustering {

    boolean[] ruptureIsAftershockList;
    int[] numAftershocksForRupture;
    protected ArrayList<Integer> numAftershocksForMainshock;
    protected ObsEqkRupList inputCatalog;
    protected ObsEqkRupList mainshockCatalog, aftershockCatalog;
    boolean isDeclustered = false;

    /**
     * Getter for mainshock catalog
     *
     * @return Declustered catalog only containing mainshocks
     */
    abstract public ObsEqkRupList getMainshockCatalog();

    /**
     * Specific declustering algorithm should be implemented here. This function should assign the member variables
     * defined in this abstract class.
     */
    abstract public void decluster();

    /**
     * Getter for aftershock catalog
     *
     * @return Catalog containing only aftershocks
     */
    abstract public ObsEqkRupList getAftershockCatalog();

    /**
     * Getter for list containing the number of aftershocks for each mainshock
     *
     * @return Aftershock count for each mainshock
     */
    abstract public ArrayList<Integer> getNumAftershocksForMainshock();
}

