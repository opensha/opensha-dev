package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.modules.PolygonFaultGridAssociations;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class ModPolyTests {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/tmp/u3_poly");
		
		FaultModels fm = FaultModels.FM3_2;
//		File refDir = new File(baseDir, "orig_buffers_"+fm.encodeChoiceString());
		File compDir = new File(baseDir, "mod_buffers_orig_loc_"+fm.encodeChoiceString());
//		File compDir = new File(baseDir, "mod_buffers_new_loc_"+fm.encodeChoiceString());
		
		File cachedFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/ucerf3/seismicityGrids/"
				+fm.name().toLowerCase()+"_fault_polygon_grid_node_associations.zip");
		PolygonFaultGridAssociations polys = new ModuleArchive<>(cachedFile).requireModule(PolygonFaultGridAssociations.class);
		
		int numSects = 0;
		int numSizeMismatch = 0;
		int numNotExactlyEqual = 0;
		int numNotApproxEqual = 0;
		
		double approxDist = 0.1;
		
		List<Feature> features = new ArrayList<>();
		
		for (File file : compDir.listFiles()) {
			if (!file.getName().endsWith(".txt") || !file.getName().startsWith("poly_"))
				continue;
			int id = Integer.parseInt(file.getName().substring(5, file.getName().indexOf(".txt")));
			System.out.println("Comparing "+file.getName()+" = "+id);
			List<String> compLines = Files.readLines(file, Charset.defaultCharset());
			LocationList locs = new LocationList();
			for (int i=0; i<compLines.size(); i++) {
				double[] compPt = load(compLines.get(i));
				locs.add(new Location(compPt[1], compPt[0]));
			}
			Region reg = throughJSON(new Region(locs, null));
			locs = reg.getBorder();
			Region origReg = throughJSON(polys.getPoly(id));
			LocationList origLocs = origReg.getBorder();
			boolean write = false;
			System.out.println("\tShape equals: "+reg.getShape().equals(origReg.getShape()));
			if (origLocs.size() != locs.size()) {
				write = true;
				System.err.println("\tSIZE MISMATCH! "+origLocs.size()+" != "+locs.size());
				numSizeMismatch++;
//				for (int i=0; i<Integer.max(locs.size(), origLocs.size()); i++) {
//					System.err.print(i+"\t");
//					if (i < locs.size())
//						System.err.println(locs.get(i)+"\tNEW");
//					else
//						System.err.println("NO NEW");
//					if (i < origLocs.size())
//						System.err.println("\t"+origLocs.get(i)+"\tORIG");
//					else
//						System.err.println("\tNO ORIG");
//				}
			} else {
				boolean exactMatch = true;
				boolean appoxMatch = true;
				for (int i=0; i<locs.size(); i++) {
					Location loc1 = locs.get(i);
					Location loc2 = origLocs.get(i);
					if (loc1.getLatitude() != loc2.getLatitude() || loc1.getLongitude() != loc2.getLongitude()) {
						exactMatch = false;
						double dist = LocationUtils.horzDistance(loc1, loc2);
						if (dist > approxDist)
							appoxMatch = false;
					}
				}
				if (!exactMatch) {
					numNotExactlyEqual++;
					System.err.println("\tNot exactly equal");
				}
				if (!appoxMatch) {
					numNotApproxEqual++;
					System.err.println("\tNot approximately equal");
					write = true;
				}
			}
			if (write) {
				Feature origFeature = Feature.getWithID(origReg.toFeature(), "Original "+id);
				origFeature.properties.set("ID", id);
				origFeature.properties.set("Modified", false);
				origFeature.properties.set(FeatureProperties.STROKE_COLOR_PROP, Color.BLUE);
				origFeature.properties.set(FeatureProperties.STROKE_WIDTH_PROP, 4);
				origFeature.properties.set(FeatureProperties.FILL_OPACITY_PROP, 0.05d);
				Feature modFeature = Feature.getWithID(reg.toFeature(), "Modified "+id);
				modFeature.properties.set("ID", id);
				modFeature.properties.set("Modified", false);
				modFeature.properties.set(FeatureProperties.STROKE_COLOR_PROP, Color.RED);
				modFeature.properties.set(FeatureProperties.STROKE_WIDTH_PROP, 2);
				modFeature.properties.set(FeatureProperties.FILL_OPACITY_PROP, 0.05d);
				features.add(origFeature);
				features.add(modFeature);
			}
			numSects++;
			System.out.flush();
			System.err.flush();
//			break;
//			List<String> refLines = Files.readLines(file, Charset.defaultCharset());
//			List<String> compLines = Files.readLines(new File(compDir, file.getName()), Charset.defaultCharset());
//			System.out.println("Comparing "+file.getName()+" = "+id);
//			numSects++;
//			if (refLines.size() != compLines.size()) {
//				System.err.println("\tSIZE MISMATCH! "+refLines.size()+" != "+compLines.size());
//				numSizeMismatch++;
//			} else {
//				int exactlyEqual = 0;
//				int floatEqual = 0;
//				for (int i=0; i<refLines.size(); i++) {
//					double[] refPt = load(refLines.get(i));
//					double[] compPt = load(compLines.get(i));
//					if (refPt[0] == compPt[0] && refPt[1] == compPt[1]) {
//						exactlyEqual++;
//						floatEqual++;
//					} else if ((float)refPt[0] == (float)compPt[0] && (float)refPt[1] == (float)compPt[1]) {
//						floatEqual++;
//					}
//				}
//				int num = refLines.size();
//				if (exactlyEqual != num) {
//					System.err.println("\t"+exactlyEqual+"/"+num+" are exactly equal");
//					System.err.println("\t"+floatEqual+"/"+num+" are approximately equal");
//					numNotExactlyEqual++;
//					if (floatEqual != num)
//						numNotApproxEqual++;
//				}
//			}
		}
		System.out.println(numSizeMismatch+"/"+numSects+" have a different number of points");
		System.out.println("Of those that have the same number of points:");
		System.out.println("\t"+numNotExactlyEqual+"/"+(numSects-numSizeMismatch)+" are not exactly equal (double precision)");
		System.out.println("\t"+numNotApproxEqual+"/"+(numSects-numSizeMismatch)+" are not approximately equal (float precision)");
		FeatureCollection.write(new FeatureCollection(features), new File(baseDir, "mismatched_polys.geojson"));
	}
	
	private static Region throughJSON(Region region) throws IOException {
		String json = region.toFeature().toJSON();
		Feature feature = Feature.fromJSON(json);
		return Region.fromFeature(feature);
	}
	
	private static double[] load(String line) {
		String[] split = line.split("\t");
		Preconditions.checkState(split.length == 2);
		double[] ret = new double[2];
		ret[0] = Double.parseDouble(split[0]);
		ret[1] = Double.parseDouble(split[1]);
		return ret;
	}

}
