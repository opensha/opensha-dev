package scratch.kevin.prvi25;

import java.awt.geom.Point2D;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import scratch.kevin.prvi25.figures.MapSourceTypeDisagg;

public class QuickHazardLandMask {

	public static void main(String[] args) {
//		List<Region> landRegions = new ArrayList<>();
//		for (XY_DataSet polBound : PoliticalBoundariesData.loadDefaultOutlines(PRVI25_RegionLoader.loadPRVI_MapExtents())) {
//			LocationList list = new LocationList();
//			for (Point2D pt : polBound)
//				list.add(new Location(pt.getY(), pt.getX()));
//			Region reg = new Region(list, BorderType.MERCATOR_LINEAR);
//			landRegions.add(reg);
//		}
//		
//		CSVFile<String> csv = CSVFile.readFile(new File("/tmp/1.0s_TWO_IN_50.csv"), false);
//		
//		MinMaxAveTracker covTrack = new MinMaxAveTracker();
//		for (int row=1; row<csv.getNumRows(); row++) {
//			Location loc = new Location(csv.getDouble(row, 1), csv.getDouble(row, 2));
//			boolean inside = false;
//			for (Region reg : caRegions) {
//				if (reg.contains(loc)) {
//					inside = true;
//					break;
//				}
//			}
//			if (inside) {
//				double cov = csv.getDouble(row, 7);
//				covTrack.addValue(cov);
//			}
//		}
//		System.out.println("CA COV: "+covTrack);
	}

}
