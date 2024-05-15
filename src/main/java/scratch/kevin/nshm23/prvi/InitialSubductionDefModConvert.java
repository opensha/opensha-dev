package scratch.kevin.nshm23.prvi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.GeoJSON_Type;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.util.minisections.MinisectionSlipRecord;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;
import com.google.gson.Gson;

public class InitialSubductionDefModConvert {

	public static void main(String[] args) throws IOException {
		File fmDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/fault_models/initial");
		File dmDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/def_models/subduction");
		File inputFile = new File(fmDir, "PRVI_sub_v2_13May2024.geojson");
		File fmOutputFile = new File(fmDir, "PRVI_sub_v2_fault_model.geojson");
		File fullOutputFile = new File(dmDir, "PRVI_sub_v2_full_minisections.txt");
		File partOutputFile = new File(dmDir, "PRVI_sub_v2_partial_minisections.txt");
		
		Gson inGSON = FeatureCollection.buildGson(DepthSerializationType.DEPTH_M);
		BufferedReader read = new BufferedReader(new FileReader(inputFile));
		List<Feature> inFeatures = inGSON.fromJson(read, FeatureCollection.class).features;
		read.close();
		
		Map<String, int[]> stitchInIDs = new HashMap<>();
		Map<String, Integer> stitchOutIDs = new HashMap<>();
		
		String name1 = "Northern Hispaniola - PRVI - Lesser Antilles Subduction";
		int[] inIDs1 = {5550, 5501, 5502, 5503, 5504, 5505, 5506, 5560};
		int outID1 = 7500;
		stitchInIDs.put(name1, inIDs1);
		stitchOutIDs.put(name1, outID1);
		
		String name2 = "Muertos Thrust";
		int[] inIDs2 = {5572, 5571, 5570};
		int outID2 = 7550;
		stitchInIDs.put(name2, inIDs2);
		stitchOutIDs.put(name2, outID2);
		
		Map<String, Feature[]> stitchFeatureSets = new HashMap<>();
		for (String name : stitchInIDs.keySet())
			stitchFeatureSets.put(name, new Feature[stitchInIDs.get(name).length]);
		
		for (Feature feature : inFeatures) {
			int id = feature.properties.getInt("id", -1);
			Preconditions.checkState(id >= 0);
			boolean found = false;
			for (String name : stitchInIDs.keySet()) {
				int[] stitchIDs = stitchInIDs.get(name);
				for (int i=0; i<stitchIDs.length; i++) {
					if (id == stitchIDs[i]) {
						stitchFeatureSets.get(name)[i] = feature;
						found = true;
						break;
					}
				}
			}
			Preconditions.checkState(found);
		}
		
		List<Feature> outFeatures = new ArrayList<>();
		
		Map<Integer, List<MinisectionSlipRecord>> fullMinisectionRecsMap = new HashMap<>();
		Map<Integer, List<MinisectionSlipRecord>> partMinisectionRecsMap = new HashMap<>();
		
		for (String name : stitchInIDs.keySet()) {
			int outID = stitchOutIDs.get(name);
			Feature[] features = stitchFeatureSets.get(name);
			
			FaultTrace stitchedUpperTrace = new FaultTrace(name);
			FaultTrace stitchedLowerTrace = new FaultTrace(name);
			
			System.out.println("Building "+name);
			
			List<MinisectionSlipRecord> fullMinisectionRecs = new ArrayList<>();
			List<MinisectionSlipRecord> partMinisectionRecs = new ArrayList<>();
			
			for (Feature feature : features) {
				Preconditions.checkState(feature.geometry.type == GeoJSON_Type.MultiLineString);
				MultiLineString geometry = (MultiLineString)feature.geometry;
				Preconditions.checkState(geometry.lines.size() == 2);
				FaultTrace upper = new FaultTrace(null);
				upper.addAll(geometry.lines.get(0));
				FaultTrace lower = new FaultTrace(null);
				lower.addAll(geometry.lines.get(1));
				
				// check that lower is in the same direction
				double upperStrike = upper.getAveStrike();
				double lowerSrike = lower.getAveStrike();
				System.out.println("\tUpper trace: depth="+depthStr(upper)+"; strike="+oneDigit.format(upperStrike));
				System.out.println("\tLower trace: depth="+depthStr(lower)+"; strike="+oneDigit.format(lowerSrike));
				double diff = FaultUtils.getAbsAngleDiff(upperStrike, lowerSrike);
				if (diff > 90d) {
					System.out.println("\tREVERSING LOWER (diff="+oneDigit.format(diff)+")");
					lower.reverse();
					geometry = new MultiLineString(List.of(upper, lower));
				}
				
				// now check Aki & Richards
				Location avgUpper = avgTraceLoc(upper);
				Location avgLower = avgTraceLoc(lower);
				double azUpperToLower = LocationUtils.azimuth(avgUpper, avgLower);
				double perfectRight = upperStrike + 90d;
				while (perfectRight > 360d)
					perfectRight -= 360d;
				double arDiff = FaultUtils.getAbsAngleDiff(azUpperToLower, perfectRight);
//				System.out.println("\tARdiff="+oneDigit.format(arDiff));
				if (arDiff > 90d) {
					System.out.println("REVERSING BOTH due to Aki & Richards violation; upper + 90 is "
							+oneDigit.format(perfectRight)+", upper to lower is "
							+oneDigit.format(azUpperToLower)+", diff is "+oneDigit.format(arDiff));
					upper.reverse();
					lower.reverse();
					geometry = new MultiLineString(List.of(upper, lower));
				}
				
				if (!stitchedUpperTrace.isEmpty()) {
					// make sure we're congiguous
					double upperDist = LocationUtils.linearDistance(stitchedUpperTrace.last(), upper.first());
					Preconditions.checkState(upperDist < 0.01d, "Upper trace not contiguous: %s, %s, dist=%s",
							stitchedUpperTrace.last(), upper.first(), upperDist);
					double lowerDist = LocationUtils.linearDistance(stitchedLowerTrace.last(), lower.first());
					Preconditions.checkState(lowerDist < 0.01d, "Lower trace not contiguous: %s, %s, dist=%s",
							stitchedLowerTrace.last(), lower.first(), lowerDist);
					
					// add all but the duplicated first point
					for (int i=1; i<upper.size(); i++)
						stitchedUpperTrace.add(upper.get(i));
					for (int i=1; i<lower.size(); i++)
						stitchedLowerTrace.add(lower.get(i));
				} else {
					// add all
					stitchedUpperTrace.addAll(upper);
					stitchedLowerTrace.addAll(lower);
				}
				
				FeatureProperties props = feature.properties;
				double slipFull = props.getDouble("RateFull", Double.NaN);
				double slipUncertFull = props.getDouble("RateF_Unc", Double.NaN);
				double rakeFull = props.getDouble("Rake", Double.NaN);
				
				double slipPart = props.getDouble("RatePart", Double.NaN);
				double slipUncertPart = props.getDouble("RateP_Unc", Double.NaN);
				double rakePart = props.getDouble("RakeRateP", Double.NaN);
				
				// TODO convert to plane rates
				for (int i=1; i<upper.size(); i++) {
					Location startLoc = upper.get(i-1);
					Location endLoc = upper.get(i);
					int minisectionID = fullMinisectionRecs.size();
					fullMinisectionRecs.add(new MinisectionSlipRecord(outID, minisectionID, startLoc, endLoc, rakeFull, slipFull, slipUncertFull));
					partMinisectionRecs.add(new MinisectionSlipRecord(outID, minisectionID, startLoc, endLoc, rakePart, slipPart, slipUncertPart));
				}
			}
			
			// build stitched feature
			MultiLineString geometry = new MultiLineString(List.of(stitchedUpperTrace, stitchedLowerTrace));
			FeatureProperties props = new FeatureProperties();
			props.set(GeoJSONFaultSection.FAULT_ID, outID);
			props.set(GeoJSONFaultSection.FAULT_NAME, name);
			props.set("PrimState", "PR");
			outFeatures.add(new Feature(outID, geometry, props));
			
			fullMinisectionRecsMap.put(outID, fullMinisectionRecs);
			partMinisectionRecsMap.put(outID, partMinisectionRecs);
		}

		FeatureCollection.write(new FeatureCollection(outFeatures), fmOutputFile);
		
		MinisectionSlipRecord.writeMinisectionsFile(fullOutputFile, fullMinisectionRecsMap);
		MinisectionSlipRecord.writeMinisectionsFile(partOutputFile, partMinisectionRecsMap);
	}
	
	private static String depthStr(FaultTrace trace) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (Location loc : trace)
			track.addValue(loc.depth);
		return "avg="+oneDigit.format(track.getAverage())+" ["+oneDigit.format(track.getMin())+", "+oneDigit.format(track.getMax())+"]";
	}
	private static DecimalFormat oneDigit = new DecimalFormat("0.0");
	
	private static Location avgTraceLoc(FaultTrace trace) {
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		MinMaxAveTracker depTrack = new MinMaxAveTracker();
		for (Location loc : trace) {
			latTrack.addValue(loc.lat);
			lonTrack.addValue(loc.lon);
			depTrack.addValue(loc.depth);
		}
		return new Location(latTrack.getAverage(), lonTrack.getAverage(), depTrack.getAverage());
	}

}
