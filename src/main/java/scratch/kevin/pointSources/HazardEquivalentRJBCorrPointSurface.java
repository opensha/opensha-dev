package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.util.Random;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.faultSurface.FiniteApproxPointSurface;

import com.google.common.base.Preconditions;

public class HazardEquivalentRJBCorrPointSurface extends FiniteApproxPointSurface {
	
	private EvenlyDiscrXYZ_DataSet rJBtable;
	private double mag;
	private Location loc;

	public HazardEquivalentRJBCorrPointSurface(EvenlyDiscrXYZ_DataSet rJBtable, Location loc, double dip, double zTop,
			double zBot, double length, double mag, boolean footwall) {
		super(loc, dip, zTop, zBot, footwall, length);
		this.rJBtable = rJBtable;
		this.loc = loc;
		this.mag = mag;
	}
	
	@Override
	public double getDistanceJB(Location siteLoc) {
		return TableBackedDistCorrPointSurface.getDist(rJBtable, LocationUtils.horzDistanceFast(loc, siteLoc), mag);
	}

}
