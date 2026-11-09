package scratch.kevin.segmentation.wasatchGeometrical;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.util.IDPairing;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.NormalFaultGeometricSegmentationProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.NormalFaultGeometricSegmentationProb.CoefficientSet;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.NormalFaultGeometricSegmentationProb.Probabilities;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.FaultTraceSegmentationGeometry.Measurements;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.FaultTraceSegmentationGeometry.TraceRelation;

/** Adds feature probabilities for both draft fits, keeping trace properties subsection-specific. */
public class WasatchSegmentationProbabilityDiagnostics {

	public static void write(File csvFile, File geometryGeoJSON, List<Measurements> measurements) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("first_id", "second_id", "first_name", "second_name", "relation", "coefficient_set", "status",
				"bend_applied", "gap_applied", "stepover_or_offset_applied", "bend_deg", "gap_km", "stepover_or_offset_km",
				"bend_probability", "gap_probability", "stepover_or_offset_probability", "combined_probability");
		Map<IDPairing, FeatureProperties> pairProperties = new TreeMap<>();
		for (Measurements m : measurements) {
			FeatureProperties properties = new FeatureProperties();
			for (CoefficientSet coefficients : CoefficientSet.values()) {
				String prefix = coefficients.name().toLowerCase(java.util.Locale.ROOT)+"_";
				List<String> row = new ArrayList<>(Arrays.asList(""+m.first.sectionID, ""+m.second.sectionID,
						m.first.name, m.second.name, m.relation.name(), coefficients.name()));
				if (m.first.parentID != m.second.parentID && m.relation == TraceRelation.AMBIGUOUS_STRIKE_AXIS) {
					row.add("AMBIGUOUS_STRIKE_AXIS");
					for (int i=0; i<10; i++) row.add("");
					properties.put(prefix+"status", "AMBIGUOUS_STRIKE_AXIS");
				} else {
					Probabilities p = new NormalFaultGeometricSegmentationProb(coefficients).calcProbabilities(m);
					String status = m.first.parentID == m.second.parentID && !p.bendApplied ? "BELOW_BEND_THRESHOLD" : "APPLIED";
					row.addAll(Arrays.asList(status, ""+p.bendApplied, ""+p.gapApplied, ""+p.stepoverApplied,
							number(m.bendAngleDegrees), number(p.gapLengthKm), number(p.stepoverDistanceKm),
							number(p.bendProbability), number(p.gapProbability), number(p.stepoverProbability), number(p.combinedProbability)));
					properties.put(prefix+"status", status);
					properties.put(prefix+"bend_probability", p.bendProbability);
					properties.put(prefix+"gap_probability", p.gapProbability);
					properties.put(prefix+"stepover_or_offset_probability", p.stepoverProbability);
					properties.put(prefix+"combined_probability", p.combinedProbability);
				}
				csv.addLine(row);
			}
			pairProperties.put(new IDPairing(m.first.sectionID, m.second.sectionID), properties);
		}
		csv.writeToFile(csvFile);
		FeatureCollection collection = FeatureCollection.read(geometryGeoJSON);
		for (Feature feature : collection.features) {
			if ("subsection_trace".equals(feature.properties.get("kind"))) continue;
			int first = ((Number)feature.properties.get("first_id")).intValue();
			int second = ((Number)feature.properties.get("second_id")).intValue();
			feature.properties.putAll(pairProperties.get(new IDPairing(first, second)));
		}
		FeatureCollection.write(collection, geometryGeoJSON);
	}

	private static String number(double value) {
		return Double.isFinite(value) ? Double.toString(value) : "";
	}
}
