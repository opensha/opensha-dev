package scratch.kevin.nshm29;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;
import org.opensha.commons.geo.json.Geometry.GeometryCollection;
import org.opensha.commons.geo.json.Geometry.LineString;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.geo.json.Geometry.Polygon;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBySectDetailPlots;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.PatchGenerationStrategy;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

public class VariableDepthFaultModelTests {
	
	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		File dataDir = new File("/home/kevin/OpenSHA/nshm29/fault_models/2026_initial_using_23_faults");
		Map<Integer, GeoJSONFaultSection> origSects = NSHM23_DeformationModels.getGeolFullSects(fm);
		
		DecimalFormat dipNumberDF = new DecimalFormat("000");
		int dipNumber = 0;
		File dataFile = new File(dataDir, "fault_traces_3d_final_single_r"+dipNumberDF.format(dipNumber)+".geojson");
		File outputDir = new File(dataDir, "results_r"+dipNumberDF.format(dipNumber));
		
		List<GeoJSONFaultSection> mappedSects = mapTraces(dataFile, origSects);
		Preconditions.checkState(outputDir.isDirectory() || outputDir.mkdirs(),
				"Couldn't create output directory: %s", outputDir);
		File sectOutputFile = new File(outputDir, "fault_sections.geojson");
		GeoJSONFaultReader.writeFaultSections(sectOutputFile, mappedSects);
		System.out.println("Wrote "+mappedSects.size()+" fault sections to "+sectOutputFile.getAbsolutePath());
		
		// now attempt to build subsections
		System.out.println("Building subsections");
		List<FaultSection> subSects = SubSectionBuilder.buildSubSects(mappedSects, 2, 0.5, 0);
		File subsectOutputFile = new File(outputDir, "fault_subsections.geojson");
		GeoJSONFaultReader.writeFaultSections(subsectOutputFile, subSects);
		System.out.println("Wrote "+subSects.size()+" fault subsections to "+subsectOutputFile.getAbsolutePath());
		
		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();

		GeographicMapMaker.GEOJSON_DEFAULT_SURF_OPACITY = 0.5;
		SectBySectDetailPlots.writeParentSectsMap(outputDir, "fault_sections_map", mappedSects, region);
		SectBySectDetailPlots.writeParentSectsMap(outputDir, "fault_subsections_map", subSects, region);
		
		// try building patches
		System.out.println("Building patches for "+subSects.size()+" subsections");
		SubSectStiffnessCalculator calc = new SubSectStiffnessCalculator(subSects, 2d, 30000, 30000, 0.5d, PatchGenerationStrategy.AUTO, 1d);
		Stopwatch watch = Stopwatch.createStarted();
		for (FaultSection subSect : subSects)
			calc.getPatches(subSect);
		watch.stop();
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Took "+(float)secs+" to build patches for all "+subSects.size()+" subsections");
	}

	/** Join the separately supplied traces without modifying the cached original sections. */
	static List<GeoJSONFaultSection> mapTraces(File dataFile, Map<Integer, GeoJSONFaultSection> origSects)
			throws IOException {
		FeatureCollection features;
		// The colleague's file declares EPSG:4326 and elevations in meters (negative
		// below datum). The normal OpenSHA GeoJSON reader instead assumes depth in km.
		try (Reader reader = new FileReader(dataFile)) {
			features = FeatureCollection.buildGson(DepthSerializationType.ELEVATION_M)
					.fromJson(reader, FeatureCollection.class);
		}
		Map<Integer, String> skipped = new TreeMap<>();
		Map<Integer, FaultTrace> upperTraces = new TreeMap<>();
		Map<Integer, FaultTrace> lowerTraces = new TreeMap<>();
		for (Feature feature : features) {
			int id = feature.properties.getInt("FaultID", -1);
			Preconditions.checkState(id >= 0, "Missing/invalid FaultID: %s", id);
			if (!origSects.containsKey(id)) {
				// Input can include faults outside the selected model (e.g. CEUS
				// Cheraw in a WUS run). Never invent the missing section metadata.
				skipped.put(id, feature.properties.get("FaultName", "unnamed"));
				continue;
			}
			String traceType = feature.properties.get("trace_type", null);
			Preconditions.checkState("top".equals(traceType) || "bottom".equals(traceType),
					"FaultID %s has unexpected trace_type: %s", id, traceType);
			Preconditions.checkState(feature.geometry instanceof LineString,
					"FaultID %s %s trace must be a LineString", id, traceType);
			LocationList locations = ((LineString)feature.geometry).line;
			Preconditions.checkState(locations.size() >= 2, "FaultID %s %s trace has fewer than 2 points", id, traceType);
			FaultTrace trace = new FaultTrace(origSects.get(id).getSectionName());
			for (Location loc : locations) {
				double depth = loc.getDepth();
				// Truncation at datum leaves a few sub-picometer positive elevations.
				if (depth < 0d && depth > -1e-9) depth = 0d;
				Preconditions.checkState(Double.isFinite(depth) && depth >= 0d,
						"FaultID %s %s has invalid depth: %s km", id, traceType, depth);
				trace.add(new Location(loc.getLatitude(), loc.getLongitude(), depth));
			}
			Map<Integer, FaultTrace> traces = traceType.equals("top") ? upperTraces : lowerTraces;
			Preconditions.checkState(traces.putIfAbsent(id, trace) == null,
					"Duplicate %s trace for FaultID %s", traceType, id);
		}
		for (var entry : skipped.entrySet())
			System.out.println("Skipping input fault outside selected model: "+entry.getKey()+" ("+entry.getValue()+")");
		List<GeoJSONFaultSection> mappedSects = new ArrayList<>();
		for (GeoJSONFaultSection original : new TreeMap<>(origSects).values()) {
			int id = original.getSectionId();
			String name = original.getSectionName();
			FaultTrace upper = upperTraces.get(id), lower = lowerTraces.get(id);
			Preconditions.checkState(upper != null && lower != null,
					"Missing top/bottom trace for FaultID %s (%s)", id, name);
			normalizeTraceOrder(id, name, upper, lower);
			FeatureProperties properties = new FeatureProperties(original.toFeature().properties);
			// Geometry is authoritative; retain rake/slip/other metadata, but replace
			// the original ribbon geometry's averages. Resampling weights by length.
			int samples = Math.max(10, (int)Math.max(upper.getTraceLength(), lower.getTraceLength()));
			properties.set(GeoJSONFaultSection.UPPER_DEPTH, averageDepth(upper, samples));
			properties.set(GeoJSONFaultSection.LOW_DEPTH, averageDepth(lower, samples));
			properties.remove(GeoJSONFaultSection.DIP);
			properties.remove(GeoJSONFaultSection.DIP_DIR);
			Geometry geometry = new MultiLineString(List.of(upper, lower));
			geometry.setSerializeZeroDepths(true);
			if (original.getZonePolygon() != null)
				geometry = new GeometryCollection(List.of(geometry, new Polygon(original.getZonePolygon())));
			try {
				// Validates trace ordering/right-hand rule and derives dip/dip direction.
				mappedSects.add(GeoJSONFaultSection.fromFeature(new Feature(id, geometry, properties)));
			} catch (IllegalArgumentException | IllegalStateException e) {
				throw new IllegalStateException("Invalid mapped geometry for FaultID "+id
						+" ("+name+"): "+e.getMessage(), e);
			}
		}
		return mappedSects;
	}

	/** Same orientation convention used internally by ApproxEvenlyGriddedSurface. */
	private static void normalizeTraceOrder(int id, String name, FaultTrace upper, FaultTrace lower) {
		double strike = upper.getAveStrike();
		if (FaultUtils.getAbsAngleDiff(strike, lower.getAveStrike()) > 90d) {
			lower.reverse();
			System.out.println("Reversed lower trace to match upper trace for "+id+". "+name);
		}
		double dipDirection = FaultUtils.getAngleAverage(List.of(
				LocationUtils.azimuth(upper.first(), lower.first()),
				LocationUtils.azimuth(upper.last(), lower.last())));
		double rightHandDirection = (strike+90d)%360d;
		if (FaultUtils.getAbsAngleDiff(rightHandDirection, dipDirection) > 90d) {
			// Reverse both lists, preserving the supplied vertices and their pairing.
			// The section itself must obey the same convention as its gridded surface.
			upper.reverse();
			lower.reverse();
			System.out.println("Reversed both traces for right-hand rule for "+id+". "+name);
		}
	}

	private static double averageDepth(FaultTrace trace, int samples) {
		return FaultUtils.resampleTrace(trace, samples).stream().mapToDouble(Location::getDepth).average().getAsDouble();
	}
}
