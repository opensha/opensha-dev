package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import scratch.kevin.simulators.ruptures.rotation.RuptureRotationUtils;

public class CatalogCSVWriter {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_4983_STITCHED.instance();
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		csv.addLine("Event ID", "Occurrence Time (s)", "Magnitude", "Moment (N-m)", "Area (m^2)",
				"Number of Participating Elements", "Average Slip (m)", "Average Element Slip Rate (m/yr)",
				"Hypocenter Latitude", "Hypocenter Longitude", "Hypocenter Depth (km)",
				"Centroid Latitude", "Centroid Longitude", "Centroid Depth (km)",
				"Upper Depth (km)", "Lower Depth (km)");
		
		CSVFile<String> csvM6p5 = new CSVFile<>(true);
		csvM6p5.addLine(csv.getLine(0));
		
		for (RSQSimEvent e : catalog.loader().skipYears(65000).iterable()) {
			List<String> line = new ArrayList<>();
			line.add(e.getID()+"");
			line.add(e.getTime()+"");
			line.add((float)e.getMagnitude()+"");
			double moment = 0d;
			for (EventRecord rec : e)
				moment += rec.getMoment();
			line.add((float)moment+"");
			line.add((float)e.getArea()+"");
			double avgSlip = 0d;
			double avgSlipRate = 0d;
			double sumArea = 0d;
			double[] slips = e.getAllElementSlips();
			List<SimulatorElement> elems = e.getAllElements();
			MinMaxAveTracker depthTrack = new MinMaxAveTracker();
			for (int i=0; i<slips.length; i++) {
				SimulatorElement elem = elems.get(i);
				double elemArea = elem.getArea();
				avgSlip += slips[i]*elemArea;
				avgSlipRate += elem.getSlipRate()*elemArea;
				sumArea += elemArea;
				for (Location loc : elem.getVertices())
					depthTrack.addValue(loc.getDepth());
			}
			avgSlip /= sumArea;
			avgSlipRate /= sumArea;
			line.add(elems.size()+"");
			line.add((float)avgSlip+"");
			line.add((float)avgSlipRate+"");
			Location hypo = RSQSimUtils.getHypocenter(e);
			line.add((float)hypo.getLatitude()+"");
			line.add((float)hypo.getLongitude()+"");
			line.add((float)hypo.getDepth()+"");
			Location centroid = RuptureRotationUtils.calcRuptureCentroid(e);
			line.add((float)centroid.getLatitude()+"");
			line.add((float)centroid.getLongitude()+"");
			line.add((float)centroid.getDepth()+"");
			line.add((float)depthTrack.getMin()+"");
			line.add((float)depthTrack.getMax()+"");
			csv.addLine(line);
			if (e.getMagnitude() >= 6.5d)
				csvM6p5.addLine(line);
		}
		
		csvM6p5.writeToFile(new File(catalog.getCatalogDir(), catalog.getCatalogDir().getName()+"_catalog_m6.5.csv"));
		csv.writeToFile(new File(catalog.getCatalogDir(), catalog.getCatalogDir().getName()+"_catalog.csv"));
	}
}
