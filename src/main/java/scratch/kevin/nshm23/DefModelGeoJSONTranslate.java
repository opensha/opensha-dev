package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;

public class DefModelGeoJSONTranslate {

	public static void main(String[] args) throws IOException {
		Map<String, String> propRemaps = new HashMap<>();
		propRemaps.put("PrefRate", "SlipRate");
		propRemaps.put("PrefRate_OnPlane", "SlipRate");
		
		propRemaps.put("LowRate_OnPlane", "LowRate");
		propRemaps.put("HighRate_OnPlane", "HighRate");
		
		propRemaps.put("Stddev", "SlipRateStdDev");
		propRemaps.put("StdDev", "SlipRateStdDev");
		propRemaps.put("Stddev", "SlipRateStdDev");
		propRemaps.put("StdMinus", "SlipRateStdDev");
		propRemaps.put("Stdminus", "SlipRateStdDev");
		propRemaps.put("StdDev_OnPlane", "SlipRateStdDev");
		
		Map<String, Class<?>> forcePropTypes = new HashMap<>();
		forcePropTypes.put("Treat", String.class);
		
//		File inputFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/"
//				+ "def_models/geologic/v1p3/NSHM23_GeolDefMod_v1p3.geojson.orig");
//		File outputFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/"
//				+ "def_models/geologic/v1p3/NSHM23_GeolDefMod_v1p3.geojson");
		File inputFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/"
				+ "def_models/geologic/v1p4/NSHM23_GeolDefMod_v1p4.geojson.orig");
		File outputFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/"
				+ "def_models/geologic/v1p4/NSHM23_GeolDefMod_v1p4.geojson");
		
		FeatureCollection features = FeatureCollection.read(inputFile);
		
		for (Feature feature : features) {
			FeatureProperties props = feature.properties;
			
			System.out.println(props.get("FaultName"));
			for (String name : propRemaps.keySet()) {
				if (props.containsKey(name)) {
					String remap = propRemaps.get(name);
					Object val = props.get(name);
					props.remove(name);
					props.put(remap, val);
					System.out.println("\t"+name+"->"+remap+": "+val);
				}
			}
			for (String name : forcePropTypes.keySet()) {
				if (props.containsKey(name)) {
					Object val = props.get(name);
					Class<?> targetClass = forcePropTypes.get(name);
					if (!targetClass.isAssignableFrom(val.getClass())) {
						System.out.println("\tRemoving '"+name+"', wrong type ("+val.getClass()+")");
						props.remove(name);
					}
				}
			}
		}
		
		FeatureCollection.write(features, outputFile);
	}

}
