package scratch.kevin.ucerf3.eal.spatialCorr;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

public class RandomFieldLoader {
	
	protected static boolean D = false;
	
	private EvenlyDiscrXYZ_DataSet field;
	private double centerX, centerY;
	
	private double wrapMinX, wrapMaxX, wrapMinY, wrapMaxY;
	
	private long numWrappedCalcs;
	private long totalNumCalcs;

	private RandomFieldLoader(EvenlyDiscrXYZ_DataSet field) {
		this.field = field;
		
		centerX = 0.5*(field.getMinX() + field.getMaxX());
		centerY = 0.5*(field.getMinY() + field.getMaxY());
		
		double deltaX = field.getGridSpacingX();
		double deltaY = field.getGridSpacingY();
		
		wrapMinX = field.getMinX() - 0.5*deltaX;
		wrapMaxX = field.getMaxX() + 0.5*deltaX;
		wrapMinY = field.getMinY() - 0.5*deltaY;
		wrapMaxY = field.getMaxY() + 0.5*deltaY;
	}
	
	public double getValue(Location siteLoc, Location ruptureCentroid) {
		LocationVector vector = LocationUtils.vector(ruptureCentroid, siteLoc);
		double dist = vector.getHorzDistance();
		double az = vector.getAzimuthRad();
		double x = centerX + dist*Math.sin(az);
		double y = centerY + dist*Math.cos(az);
		
		if (D) System.out.println("GET "+totalNumCalcs);
		if (D) System.out.println("\tx="+(float)x+"\ty="+(float)y);
		
//		field.getX(xIndex)
		boolean wrapped = false;
		while (x > wrapMaxX) {
			wrapped = true;
			x -= (wrapMaxX-wrapMinX);
		}
		while (x < wrapMinX) {
			wrapped = true;
			x += (wrapMaxX-wrapMinX);
		}
		while (y > wrapMaxY) {
			wrapped = true;
			y -= (wrapMaxY-wrapMinY);
		}
		while (y < wrapMinX) {
			wrapped = true;
			y += (wrapMaxY-wrapMinY);
		}
		
		if (wrapped == true)
			numWrappedCalcs++;
		totalNumCalcs++;
		
		int xInd = field.getXIndex(x);
		int yInd = field.getYIndex(y);
		
		if (D && wrapped)
			System.out.println("\t\twrapped! new x="+(float)x+"\ty="+(float)y);
		
		if (D) System.out.println("\txInd="+xInd+"\tyInd="+yInd);
		
		double ret = field.get(xInd, yInd);
		
		if (D) System.out.println("\tRET: "+ret);
		
		return ret;
	}

	public long getNumWrappedCalcs() {
		return numWrappedCalcs;
	}

	public long getTotalNumCalcs() {
		return totalNumCalcs;
	}
	
	public static RandomFieldLoader load(File csvFile, double gridSpacing) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		int minXIndex = csv.getInt(1, 0);
		int minYIndex = csv.getInt(1, 1);
		int maxXIndex = csv.getInt(csv.getNumRows()-1, 0);
		int maxYIndex = csv.getInt(csv.getNumRows()-1, 1);
		
		int nx = 1 + maxXIndex - minXIndex;
		int ny = 1 + maxYIndex - minYIndex;
		
		double minX = 0d;
		double minY = 0d;
		
		EvenlyDiscrXYZ_DataSet field = new EvenlyDiscrXYZ_DataSet(nx, ny, minX, minY, gridSpacing);
		for (int row=1; row<csv.getNumRows(); row++) {
			int x = csv.getInt(row, 0)-minXIndex;
			int y = csv.getInt(row, 1)-minYIndex;
			field.set(x, y, csv.getDouble(row, 2));
		}
		
		return new RandomFieldLoader(field);
	}

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/eal/random_fields/sa10_1km_800x800");
		
		RandomFieldLoader loader = load(new File(dir, "800x800SA10_001.csv"), 1d);
		System.out.println(loader.field.getNumX()+" "+loader.field.getNumY());
		System.out.println(loader.field.getMaxX()+" "+loader.field.getMaxY());
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (int i=0; i<loader.field.size(); i++)
			track.addValue(loader.field.get(i));
		System.out.println("Field stats: "+track);
		
		Location centroid = new Location(34, -118);
		
		D = true;

		loader.getValue(new Location(34, -118), centroid);
		loader.getValue(new Location(33.99, -118), centroid);
		loader.getValue(new Location(0, -0), centroid);
	}

}
