package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoProbabilityModel;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class PaleoSiteCSVWriter {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
		
		File outputDir = new File("/tmp/rsqsim_paleo");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<String> siteNames = new ArrayList<>();
		List<Location> siteLocs = new ArrayList<>();
		
		siteNames.add("HL");
		siteLocs.add(new Location(33.6153, -116.7091));
		
		siteNames.add("ML");
		siteLocs.add(new Location(33.9, -117.1));
		
		siteNames.add("WW");
		siteLocs.add(new Location(34.3697, -117.6680));
		
		siteNames.add("PC");
		siteLocs.add(new Location(34.4556, -117.8870));
		
		siteNames.add("FP");
		siteLocs.add(new Location(34.8122, -118.9034));
		
		List<SimulatorElement> elems = catalog.getElements();
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		
		List<SimulatorElement> mappings = new ArrayList<>();
		List<CSVFile<String>> siteCSVs = new ArrayList<>();
		
		Map<Integer, CSVFile<String>> elemCSVs = new HashMap<>();
		for (int i=0; i<siteLocs.size(); i++) {
			Location siteLoc = siteLocs.get(i);
			// find closest
			SimulatorElement closest = null;
			double minDist = Double.POSITIVE_INFINITY;
			for (SimulatorElement elem : elems) {
				double dist = Double.POSITIVE_INFINITY;
				for (Location loc : elem.getVertices())
					dist = Math.min(dist, LocationUtils.linearDistanceFast(loc, siteLoc));
				if (dist < minDist) {
					minDist = dist;
					closest = elem;
				}
			}
			System.out.println("Closest for "+siteNames.get(i)+" at "+siteLoc+": element "+closest.getID()+", "+(float)minDist+" km away");
			System.out.println("\tMapped U3 section: "+mapper.getMappedSection(closest).getSectionName());
			System.out.println("\tElement vertexes:");
			for (Location loc : closest.getVertices())
				System.out.println("\t\t"+loc);
			
			mappings.add(closest);
			
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Event ID", "Occurrence Time (s)", "Magnitude", "Average Slip (m)", "Slip On Matching Element (m)", "U3 Paleo Obs Prob");
			siteCSVs.add(csv);
			
			elemCSVs.put(closest.getID(), csv);
		}
		
		UCERF3_PaleoProbabilityModel paleoProbModel = UCERF3_PaleoProbabilityModel.load();
		
		int numRescuedProbs = 0;
		int numTotProbs = 0;
		for (RSQSimEvent e : catalog.loader().skipYears(65000).iterable()) {
			double aveSlip = 0d;
			int[] ids = e.getAllElementIDs();
			double[] slips = e.getAllElementSlips();
			double sumArea = 0d;
			for (int i=0; i<ids.length; i++) {
				double area = elems.get(i).getArea();
				aveSlip += area*slips[i];
				sumArea += area;
			}
			aveSlip /= sumArea;
			
			for (int i=0; i<ids.length; i++) {
				CSVFile<String> csv = elemCSVs.get(ids[i]);
				if (csv != null) {
					// ruptures this site
					List<String> line = new ArrayList<>();
					line.add(e.getID()+"");
					line.add(e.getTime()+"");
					line.add(e.getMagnitude()+"");
					line.add(aveSlip+"");
					line.add(slips[i]+"");
					RSQSimSubSectEqkRupture mapped = catalog.getMappedSubSectRupture(e);
					int mappedSectID = mapper.getMappedSection(elems.get(ids[i])).getSectionId();
					double obsProb;
					try {
						obsProb = paleoProbModel.getProbPaleoVisible(e.getMagnitude(), mapped.getSubSections(), mappedSectID);
					} catch (IllegalStateException e1) {
						List<List<SubSectionMapping>> allSects  = mapper.getAllSubSectionMappings(e);
						List<FaultSection> rupSects = new ArrayList<>();
						for (List<SubSectionMapping> bundle : allSects)
							for (SubSectionMapping mapping : bundle)
								rupSects.add(mapping.getSubSect());
						numRescuedProbs++;
						obsProb = paleoProbModel.getProbPaleoVisible(e.getMagnitude(), rupSects, mappedSectID);
					}
					numTotProbs++;
					line.add(obsProb+"");
					csv.addLine(line);
				}
			}
		}
		System.out.println("Rescued "+numRescuedProbs+"/"+numTotProbs+" paleo probs");
		
		for (int i=0; i<siteCSVs.size(); i++)
			siteCSVs.get(i).writeToFile(new File(outputDir, siteNames.get(i)+".csv"));
	}

}
