package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.util.Random;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;

import com.google.common.base.Preconditions;

public class TableBackedDistCorrPointSurface extends FiniteApproxPointSurface {
	
	public static class DistanceCorrTables {
		public final EvenlyDiscrXYZ_DataSet rJB;
		public final EvenlyDiscrXYZ_DataSet rRup;
		public final EvenlyDiscrXYZ_DataSet rX;
		public final EvenlyDiscrXYZ_DataSet rXPositive;
		private DistanceCorrTables(EvenlyDiscrXYZ_DataSet rJB, EvenlyDiscrXYZ_DataSet rRup, EvenlyDiscrXYZ_DataSet rX,
				EvenlyDiscrXYZ_DataSet rXPositive) {
			super();
			this.rJB = rJB;
			this.rRup = rRup;
			this.rX = rX;
			this.rXPositive = rXPositive;
		}
	}
	
	static DistanceCorrTables loadCorrTables(File dir, String prefix) throws IOException {
		EvenlyDiscrXYZ_DataSet rJB = loadXYZ_CSV(new File(dir, prefix+"_rJB.csv"));
		EvenlyDiscrXYZ_DataSet rRup = loadXYZ_CSV(new File(dir, prefix+"_rRup.csv"));
		EvenlyDiscrXYZ_DataSet rX = loadXYZ_CSV(new File(dir, prefix+"_rX.csv"));
		EvenlyDiscrXYZ_DataSet rXPositive = loadXYZ_CSV(new File(dir, prefix+"_rX_positive.csv"));
		return new DistanceCorrTables(rJB, rRup, rX, rXPositive);
	}
	
	static EvenlyDiscrXYZ_DataSet loadXYZ_CSV(File csvFile) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		double minMag = csv.getDouble(0, 1);
		double magSpacing = csv.getDouble(0, 2) - minMag;
		int numMag = csv.getNumCols()-1;
		
		double minDist = csv.getDouble(1, 0);
		double distSpacing = csv.getDouble(2, 0) - minDist;
		int numDist = csv.getNumRows()-1;
		
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(numDist, numMag, minDist, minMag, distSpacing, magSpacing);
		for (int d=0; d<numDist; d++) {
			Preconditions.checkState((float)xyz.getX(d) == (float)csv.getDouble(d+1, 0));
			for (int m=0; m<numMag; m++) {
				Preconditions.checkState(d > 0 || (float)xyz.getY(m) == (float)csv.getDouble(0, m+1));
				xyz.set(d, m, csv.getDouble(d+1, m+1));
			}
		}
		return xyz;
	}
	
	private DistanceCorrTables corrTables;
	private double mag;
	private Location loc;
	
	private Random rand;

	public TableBackedDistCorrPointSurface(DistanceCorrTables corrTables, Location loc, double dip, double zTop,
			double zBot, double length, double mag, Random rand) {
		super(loc, dip, zTop, zBot, true, length);
		this.corrTables = corrTables;
		this.loc = loc;
		this.mag = mag;
		this.rand = rand;
	}
	
	static double getDist(EvenlyDiscrXYZ_DataSet distMagFunc, double dist, double mag) {
		// buffer this check by 10km as we can have supersampled locations that are a little further
		Preconditions.checkState((float)dist >= (float)distMagFunc.getMinX()
				&& (float)dist <= (float)(distMagFunc.getMaxX()+10), "Bad dist=%s with range=[%s, %s]",
				(float)dist, (float)distMagFunc.getMinX(), (float)distMagFunc.getMaxX());
		if (dist < distMagFunc.getMinX())
			dist = distMagFunc.getMinX();
		else if (dist > distMagFunc.getMaxX())
			dist = distMagFunc.getMaxX();
		Preconditions.checkState((float)mag >= (float)distMagFunc.getMinY()
				&& (float)mag <= (float)distMagFunc.getMaxY(), "Bad mag=%s with range=[%s, %s]",
						(float)mag, (float)distMagFunc.getMinY(), (float)distMagFunc.getMaxY());
		if (mag < distMagFunc.getMinY())
			mag = distMagFunc.getMinY();
		else if (mag > distMagFunc.getMaxY())
			mag = distMagFunc.getMaxY();
		return distMagFunc.bilinearInterpolation(dist, mag);
	} 

	@Override
	public double getDistanceX(Location siteLoc) {
		double dist = LocationUtils.horzDistanceFast(this.loc, siteLoc);
		double rX = getDist(corrTables.rX, dist, mag);
		boolean positive = rand.nextDouble() < getDist(corrTables.rXPositive, dist, mag);
		return positive ? rX : -rX;
	}
	
	@Override
	public double getDistanceJB(Location siteLoc) {
		return getDist(corrTables.rJB, LocationUtils.horzDistanceFast(loc, siteLoc), mag);
	}

	@Override
	public double getDistanceRup(Location siteLoc) {
		return getDist(corrTables.rRup, LocationUtils.horzDistanceFast(loc, siteLoc), mag);
	}

	@Override
	public double getDistanceSeis(Location loc) {
		throw new IllegalStateException();
	}

}
