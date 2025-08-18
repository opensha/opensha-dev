package scratch.kevin.prvi25;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
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

public class CA_COV_Compare {

	public static void main(String[] args) throws IOException {
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		List<Region> caRegions = new ArrayList<>();
		for (int i=0; i<caOutlines.length; i++) {
			LocationList outline = new LocationList();
			for (Point2D pt : caOutlines[i])
				outline.add(new Location(pt.getY(), pt.getX()));
			caRegions.add(new Region(outline, BorderType.MERCATOR_LINEAR));
		}
		
		CSVFile<String> csv = CSVFile.readFile(new File("/tmp/1.0s_TWO_IN_50.csv"), false);
		
		MinMaxAveTracker covTrack = new MinMaxAveTracker();
		for (int row=1; row<csv.getNumRows(); row++) {
			Location loc = new Location(csv.getDouble(row, 1), csv.getDouble(row, 2));
			boolean inside = false;
			for (Region reg : caRegions) {
				if (reg.contains(loc)) {
					inside = true;
					break;
				}
			}
			if (inside) {
				double cov = csv.getDouble(row, 7);
				covTrack.addValue(cov);
			}
		}
		System.out.println("CA COV: "+covTrack);
	}

}
