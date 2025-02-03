package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.GeoJSON_Type;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.GeometryCollection;
import org.opensha.commons.geo.json.Geometry.LineString;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.geo.json.Geometry.MultiPolygon;
import org.opensha.commons.geo.json.Geometry.Polygon;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class CrustalDefModConvert {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/fault_models/crustal");
//		File dir = new File("C:\\Users\\kmilner\\git\\opensha\\opensha\\src\\main\\resources\\data\\erf\\prvi25\\fault_models\\initial");
//		File inputFile = new File(dir, "NSHM2025_GeoDefModel_PRVI_v1-1.geojson");
//		File inputPolys = new File(dir, "NSHM2025_FaultPolygons_PRVI_v1.geojson");
//		File outputFile = new File(dir, "NSHM2025_GeoDefModel_PRVI_v1-1_mod.geojson");
		File inputFile = new File(dir, "NSHM2025_GeoDefModel_PRVI_v1-1_ProjRates.geojson");
		File inputPolys = new File(dir, "NSHM2025_FaultPolygons_PRVI_v1.geojson");
		File outputFile = new File(dir, "NSHM2025_GeoDefModel_PRVI_v1-2_ProjRates_mod.geojson");
		List<Feature> features = new ArrayList<>(FeatureCollection.read(inputFile).features);
		List<Feature> polyFeatures = new ArrayList<>(FeatureCollection.read(inputPolys).features);
		
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
			
			String origID = props.getString("FaultID");
			Preconditions.checkNotNull(origID);
			
			String name = props.getString("FaultName");
			props.remove("FaultName");
			props.set(GeoJSONFaultSection.FAULT_NAME, name);
			
			System.out.println("Processing "+id+". "+origID+": "+name);
			
			Preconditions.checkState(!prevNames.contains(name), "Duplicate name: %s", name);
			prevNames.add(name);
			
			props.set("FaultID", id);
			
			double upDepth = props.getDouble("UpDep", Double.NaN);
			props.remove("UpDep");
			props.set(GeoJSONFaultSection.UPPER_DEPTH, upDepth);

			double lowDepth = props.getDouble("LowerDepth", Double.NaN);
			props.remove("LowerDepth");
//			double lowDepth = props.getDouble("LowerDep", Double.NaN);
//			props.remove("LowerDep");
			props.set(GeoJSONFaultSection.LOW_DEPTH, lowDepth);
			
			double dip = props.getDouble("FaultDip", Double.NaN);
			FaultUtils.assertValidDepth(dip);
			props.remove("FaultDip");
			props.set(GeoJSONFaultSection.DIP, dip);
			
			props.remove("DipDir");
			
			double slip = props.getDouble("PrefRate_Projected", Double.NaN);
			props.remove("PrefRate_Projected");
//			double slip = props.getDouble("PrefSlpRat", Double.NaN);
//			props.remove("PrefSlpRat");
			props.set(GeoJSONFaultSection.SLIP_RATE, slip);

			double lowSlip = props.getDouble("MinRate_Projected", Double.NaN);
			props.remove("MinRate_Projected");
//			double lowSlip = props.getDouble("MinSlpRat", Double.NaN);
//			props.remove("MinSlpRat");
			props.set("LowRate", lowSlip);

			double highSlip = props.getDouble("MaxRate_Projected", Double.NaN);
			props.remove("MaxRate_Projected");
//			double highSlip = props.getDouble("MaxSlpRat", Double.NaN);
//			props.remove("MaxSlpRat");
			props.set("HighRate", highSlip);
			
			double rake = props.getDouble("Rake", Double.NaN);
			FaultUtils.assertValidRake(rake);
			props.set("Rake", rake);
			
			LineString trace;
			if (feature.geometry.type == GeoJSON_Type.MultiLineString) {
				MultiLineString multiLine = (MultiLineString)feature.geometry;
				Preconditions.checkState(multiLine.lines.size() == 1,
						"Expected 1 line, have %1", multiLine.lines.size());
				trace = new LineString(multiLine.lines.get(0));
			} else {
				Preconditions.checkState(feature.geometry.type == GeoJSON_Type.LineString);
				trace = (LineString)feature.geometry;
			}
			
			Geometry geom = trace;
			
			// see if we have a polygon
			Feature polyFeature = null;
			for (Feature test : polyFeatures)
				if (origID.equalsIgnoreCase(test.properties.getString("PolyID")))
					polyFeature = test;
			if (polyFeature != null) {
				FeatureProperties polyProps = polyFeature.properties;
				String polyName = polyProps.getString("FZ_Name");
				System.out.println("\tmatched polygon for "+name+": "+polyName);
				Preconditions.checkState(name.startsWith(polyName), "'%s' != '%s'", name, polyName);
				double polyDip = polyProps.getDouble("FaultDip", Double.NaN);
				Preconditions.checkState((float)dip == (float)polyDip, "Poly dip mismatch: %s != %s",
						(float)dip, (float)polyDip);
				double polyRake = polyProps.getDouble("Rake", Double.NaN);
				Preconditions.checkState((float)rake == (float)polyRake, "Poly rake mismatch: %s != %s",
						(float)rake, (float)polyRake);
				double polyUpper = polyProps.getDouble("UpDep", Double.NaN);
				Preconditions.checkState((float)upDepth == (float)polyUpper, "Poly up depth mismatch: %s != %s",
						(float)upDepth, (float)polyUpper);
				double polyLower = polyProps.getDouble("LowerDep", Double.NaN);
				Preconditions.checkState((float)lowDepth == (float)polyLower, "Poly low depth mismatch: %s != %s",
						(float)lowDepth, (float)polyLower);
				
				Preconditions.checkState(polyFeature.geometry.type == GeoJSON_Type.Polygon
						|| polyFeature.geometry.type == GeoJSON_Type.MultiPolygon);
				Polygon poly;
				if (polyFeature.geometry.type == GeoJSON_Type.MultiPolygon) {
					MultiPolygon multiPoly = (MultiPolygon)polyFeature.geometry;
					Preconditions.checkState(multiPoly.polygons.size() == 1,
							"Expected 1 poly, have %1", multiPoly.polygons.size());
					poly = multiPoly.polygons.get(0);
				} else {
					poly = (Polygon)polyFeature.geometry;
				}
				
				Region reg = Region.fromFeature(polyFeature);
				GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
				System.out.println("\tGridded polygon has "+gridReg.getNodeCount()+" nodes if 0.1 degrees");
				gridReg = new GriddedRegion(reg, 0.05, GriddedRegion.ANCHOR_0_0);
				System.out.println("\tGridded polygon has "+gridReg.getNodeCount()+" nodes if 0.05 degrees");
				
				props.set(GeoJSONFaultSection.PROXY, true);
				
				geom = new GeometryCollection(List.of(trace, poly));
			}
			
			features.set(i, new Feature(id, geom, props));
		}
		
		FeatureCollection.write(new FeatureCollection(features), outputFile);
	}

}
