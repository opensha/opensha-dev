package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FaultUtils.AngleAverager;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import scratch.kevin.simulators.ruptures.rotation.RuptureRotationUtils;

public class CatalogCSVWriter {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_4983_STITCHED.instance();
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		boolean addLengths = true;
		boolean addRakes = true;
		boolean addDips = true;
		
		List<String> header = new ArrayList<>(List.of("Event ID", "Occurrence Time (s)", "Magnitude", "Moment (N-m)",
				"Area (m^2)", "Number of Participating Elements", "Average Slip (m)", "Average Element Slip Rate (m/yr)",
				"Hypocenter Latitude", "Hypocenter Longitude", "Hypocenter Depth (km)", "Centroid Latitude",
				"Centroid Longitude", "Centroid Depth (km)", "Upper Depth (km)", "Lower Depth (km)"));

		if (addRakes) {
			header.add("Average Rake");
			header.add("Average Absolute Rake Deviation");
			header.add("Max Absolute Rake Deviation");
		}
		if (addDips) {
			header.add("Average Dip");
			header.add("Average Absolute Dip Deviation");
			header.add("Max Absolute Dip Deviation");
		}
		if (addLengths) {
			for (SlipAlongSectAlgorithm type : SlipAlongSectAlgorithm.values())
				header.add(type.toString()+" (km)");
		}
		csv.addLine(header);
		
		CSVFile<String> csvM6p5 = new CSVFile<>(true);
		csvM6p5.addLine(csv.getLine(0));
		
		RSQSimSubSectionMapper mapper = null;
		if (addLengths) {
			mapper = catalog.getSubSectMapper();
			mapper.trackSlipOnSections();
		}
		
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
			
			if (addRakes) {
				AngleAverager avg = new AngleAverager();
				
				for (int i=0; i<slips.length; i++)
					// weighted by slip
					avg.add(elems.get(i).getFocalMechanism().getRake(), slips[i]);
				
				double aveRake = FaultUtils.getInRakeRange(avg.getAverage());
				double aveDeviation = 0d;
				double weightSum = 0d;
				double maxDeviation = 0d;
				for (int i=0; i<slips.length; i++) {
					double rake = elems.get(i).getFocalMechanism().getRake();
					double diff = angleDiff(aveRake, rake);
					maxDeviation = Math.max(diff, maxDeviation);
					aveDeviation += diff*slips[i];
					weightSum += slips[i];
				}
				aveDeviation /= weightSum;
				line.add((float)aveRake+"");
				line.add((float)aveDeviation+"");
				line.add((float)maxDeviation+"");
			}
			
			if (addDips) {
				double aveDip = 0d;
				double weightSum = 0d;
				
				for (int i=0; i<slips.length; i++) {
					aveDip += elems.get(i).getFocalMechanism().getDip()*slips[i];
					weightSum += slips[i];
				}
				aveDip /= weightSum;
				
				double aveDeviation = 0d;
				double maxDeviation = 0d;
				for (int i=0; i<slips.length; i++) {
					double diff =  Math.abs(elems.get(i).getFocalMechanism().getDip() - aveDip);
					aveDeviation += diff * slips[i];
					maxDeviation = Math.max(diff, maxDeviation);
				}
				aveDeviation /= weightSum;
				line.add((float)aveDip+"");
				line.add((float)aveDeviation+"");
				line.add((float)maxDeviation+"");
			}
			
			if (addLengths) {
				List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings(e);
				for (SlipAlongSectAlgorithm alg : SlipAlongSectAlgorithm.values()) {
					double length = 0d;
					for (List<SubSectionMapping> bundle : mappings)
						for (SubSectionMapping mapping : bundle)
							length += mapping.getLengthForSlip(alg);
					line.add((float)length+"");
				}
			}
			
			csv.addLine(line);
			if (e.getMagnitude() >= 6.5d)
				csvM6p5.addLine(line);
		}
		
		String prefix = catalog.getCatalogDir().getName();
		if (addRakes)
			prefix += "_rakes";
		if (addDips)
			prefix += "_dips";
		if (addLengths)
			prefix += "_lengths";
		
		csvM6p5.writeToFile(new File(catalog.getCatalogDir(), prefix+"_m6.5.csv"));
		csv.writeToFile(new File(catalog.getCatalogDir(), prefix+".csv"));
	}
	
	private static double angleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		return Math.abs(angleDiff);
	}
}
