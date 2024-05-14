package scratch.kevin.nshm23.prvi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class InitialDefModConvert {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/fault_models/initial");
		File inputFile = new File(dir, "NSHM2025_GeoDefModel_PRVI.geojson");
		File outputFile = new File(dir, "NSHM2025_GeoDefModel_PRVI_mod.geojson");
		List<Feature> features = new ArrayList<>(FeatureCollection.read(inputFile).features);
		
		features.sort(new Comparator<Feature>() {

			@Override
			public int compare(Feature o1, Feature o2) {
				String name1 = o1.properties.getString("FaultID");
				String name2 = o2.properties.getString("FaultID");
				Preconditions.checkState(name1.startsWith("FS") || name1.startsWith("FP"));
				Preconditions.checkState(name2.startsWith("FS") || name2.startsWith("FP"));
				if (name1.startsWith("FS") && name2.startsWith("FS")) {
					int id1 = Integer.parseInt(name1.substring(2));
					int id2 = Integer.parseInt(name2.substring(2));
					return Integer.compare(id1, id2);
				} else if (name1.startsWith("FP") && name2.startsWith("FP")) {
					int id1 = Integer.parseInt(name1.substring(2));
					int id2 = Integer.parseInt(name2.substring(2));
					return Integer.compare(id1, id2);
				} else if (name1.startsWith("FS")) {
					return -1;
				} else {
					return 1;
				}
			}
		});
		
		HashSet<String> prevNames = new HashSet<>();
		
		for (int i=0; i<features.size(); i++) {
			Feature feature = features.get(i);
			feature.properties.remove("id");
			
			int id = i+1;
			feature = Feature.getWithID(feature, id);
			
			FeatureProperties props = feature.properties;
			
			String name = props.getString("Fault");
			props.remove("Fault");
			props.set(GeoJSONFaultSection.FAULT_NAME, name);
			
			Preconditions.checkState(!prevNames.contains(name), "Duplicate name: %s", name);
			prevNames.add(name);
			
			props.set("FaultID", id);
			
			double upDepth = props.getDouble("UpDep", Double.NaN);
			props.remove("UpDep");
			props.set(GeoJSONFaultSection.UPPER_DEPTH, upDepth);
			
			double lowDepth = props.getDouble("LowerDep", Double.NaN);
			props.remove("LowerDep");
			props.set(GeoJSONFaultSection.LOW_DEPTH, lowDepth);
			
			double dip = props.getDouble("FaultDip", Double.NaN);
			props.remove("FaultDip");
			props.set(GeoJSONFaultSection.DIP, dip);
			
			props.remove("DipDir");
			
			double slip = props.getDouble("PrefSlpRat", Double.NaN);
			props.remove("PrefSlpRat");
			props.set(GeoJSONFaultSection.SLIP_RATE, slip);
			
			double lowSlip = props.getDouble("MinSlpRat", Double.NaN);
			props.remove("MinSlpRat");
			props.set("LowRate", lowSlip);
			
			double highSlip = props.getDouble("MaxSlpRat", Double.NaN);
			props.remove("MaxSlpRat");
			props.set("HighRate", highSlip);
			
			features.set(i, feature);
		}
		
		FeatureCollection.write(new FeatureCollection(features), outputFile);
	}

}
