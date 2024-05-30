package scratch.kevin.nshm23.prvi;

import java.awt.Color;
import java.awt.Font;
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

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.GeoJSON_Type;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.util.minisections.MinisectionSlipRecord;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;
import com.google.gson.Gson;

public class InitialSubductionDefModConvert {

	public static void main(String[] args) throws IOException {
		File fmDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/fault_models/subduction");
		File dmDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/def_models/subduction");
		File inputsDir = new File(fmDir, "inputs");
		
//		File inputFile = new File(fmDir, "PRVI_sub_v2_13May2024.geojson");
//		File fmOutputFile = new File(fmDir, "PRVI_sub_v2_fault_model.geojson");
//		File fullOutputFile = new File(dmDir, "PRVI_sub_v2_full_minisections.txt");
//		File partOutputFile = new File(dmDir, "PRVI_sub_v2_partial_minisections.txt");
		
		
		boolean fullRate = true;
		boolean largePolys = true;
		
		String inRatePrefix = fullRate ? "FullRate" : "PartRate";
		String outRatePrefix = fullRate ? "full" : "part";
		String inGeomPrefix = largePolys ? "LargePolys" : "SmallPolys";
		String outGeomPrefix = largePolys ? "large" : "small";
		File inputFile = new File(inputsDir, "PRVI_subduction_"+inRatePrefix+"_"+inGeomPrefix+"_drape_removeFields.geojson");
		File fmOutputFile = new File(fmDir, "PRVI_sub_v3_fault_model_"+outGeomPrefix+".geojson");
		File dmOutputFile = new File(dmDir, "PRVI_sub_v3_"+outGeomPrefix+"_"+outRatePrefix+"_rate_minisections.txt");
		String ratePropName = fullRate ? "RateF_Proj" : "RateP_Proj";
		String rateUncertPropName = fullRate ? "RateF_Unc" : "RateP_Unc";
		String rakePropName = fullRate ? "RakeRateF" : "RakeRateU";
		
		Gson inGSON = FeatureCollection.buildGson(DepthSerializationType.DEPTH_M);
		BufferedReader read = new BufferedReader(new FileReader(inputFile));
		List<Feature> inFeatures = inGSON.fromJson(read, FeatureCollection.class).features;
		read.close();
		
		debugWriteSectOrders(inFeatures, new File("/tmp"),
				"geom_"+inGeomPrefix+"_"+inRatePrefix,
				inputFile.getName().replace(".geojson", ""));
		
		Map<String, int[]> stitchInIDs = new HashMap<>();
		Map<String, Integer> stitchOutIDs = new HashMap<>();
		
		String name1 = "Northern Hispaniola - PRVI - Lesser Antilles Subduction";
		int[] inIDs1 = {7501, 7502, 7503, 7504, 7505, 7506, 7507, 7508};
		int outID1 = 7500;
		stitchInIDs.put(name1, inIDs1);
		stitchOutIDs.put(name1, outID1);
		
		String name2 = "Muertos Thrust";

		int[] inIDs2 = {7553, 7552, 7551};
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
		
		Map<Integer, List<MinisectionSlipRecord>> minisectionRecsMap = new HashMap<>();
		
		for (String name : stitchInIDs.keySet()) {
			int outID = stitchOutIDs.get(name);
			Feature[] features = stitchFeatureSets.get(name);
			
			FaultTrace stitchedUpperTrace = new FaultTrace(name);
			FaultTrace stitchedLowerTrace = new FaultTrace(name);
			
			System.out.println("Building "+name);
			
			List<MinisectionSlipRecord> minisectionRecs = new ArrayList<>();
			
			for (Feature feature : features) {
				Preconditions.checkState(feature.geometry.type == GeoJSON_Type.MultiLineString);
				MultiLineString geometry = (MultiLineString)feature.geometry;
				Preconditions.checkState(geometry.lines.size() == 2);
				FaultTrace upper = new FaultTrace(null);
				upper.addAll(geometry.lines.get(0));
				FaultTrace lower = new FaultTrace(null);
				lower.addAll(geometry.lines.get(1));
				
				int sectID = feature.properties.getInt("id", -1);
				String sectName = feature.properties.getString("FaultName")+" ("+feature.properties.getString("FaultDesc")+")";
				System.out.println("Processing "+sectID+". "+sectName);
				
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
					Preconditions.checkState(upperDist < 0.01d, "Upper trace not contiguous: %s, %s, dist=%s, %s. %s",
							stitchedUpperTrace.last(), upper.first(), upperDist, sectID, sectName);
					double lowerDist = LocationUtils.linearDistance(stitchedLowerTrace.last(), lower.first());
					Preconditions.checkState(lowerDist < 0.01d, "Lower trace not contiguous: %s, %s, dist=%s, %s. %s",
							stitchedLowerTrace.last(), lower.first(), lowerDist, sectID, sectName);
					
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
				double slip = props.getDouble(ratePropName, Double.NaN);
				double slipUncert = props.getDouble(rateUncertPropName, Double.NaN);
				double rake = props.getDouble(rakePropName, Double.NaN);
				
				// TODO convert to plane rates
				for (int i=1; i<upper.size(); i++) {
					Location startLoc = upper.get(i-1);
					Location endLoc = upper.get(i);
					int minisectionID = minisectionRecs.size();
					minisectionRecs.add(new MinisectionSlipRecord(outID, minisectionID, startLoc, endLoc, rake, slip, slipUncert));
				}
			}
			
			// build stitched feature
			MultiLineString geometry = new MultiLineString(List.of(stitchedUpperTrace, stitchedLowerTrace));
			FeatureProperties props = new FeatureProperties();
			props.set(GeoJSONFaultSection.FAULT_ID, outID);
			props.set(GeoJSONFaultSection.FAULT_NAME, name);
			props.set("PrimState", "PR");
			outFeatures.add(new Feature(outID, geometry, props));
			
			minisectionRecsMap.put(outID, minisectionRecs);
		}

		if (fmOutputFile != null)
			FeatureCollection.write(new FeatureCollection(outFeatures), fmOutputFile);
		
		MinisectionSlipRecord.writeMinisectionsFile(dmOutputFile, minisectionRecsMap);
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
	
	private static void debugWriteSectOrders(List<Feature> features, File outputDir, String prefix, String title)
			throws IOException {
		List<FaultSection> sects = new ArrayList<>();
		List<Location> directionLocs = new ArrayList<>();
		List<Double> directionFracts = new ArrayList<>();
		List<Double> depths = new ArrayList<>();
		CPT orderCPT = new CPT(0d, 1d, Color.RED, Color.GREEN);
		for (Feature feature : features) {
			Feature copy = new Feature(feature.properties.getNumber("id"), feature.geometry, feature.properties);
			sects.add(GeoJSONFaultSection.fromFeature(copy));
			Preconditions.checkState(feature.geometry.type == GeoJSON_Type.MultiLineString);
			MultiLineString multi = (MultiLineString)feature.geometry;
			for (LocationList line : multi.lines) {
				FaultTrace trace = new FaultTrace(null);
				trace.addAll(line);
				trace = FaultUtils.resampleTrace(trace, (int)(trace.getTraceLength()+0.5));
				for (int l=0; l<trace.size(); l++) {
					directionLocs.add(trace.get(l));
					directionFracts.add((double)l/(double)(trace.size()-1));
					depths.add(trace.get(l).depth);
				}
			}
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 0.3f);
		mapMaker.plotScatterScalars(directionLocs, directionFracts, orderCPT, "First ----------> Last");
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		
		for (FaultSection sect : sects) {
			LocationList locs = sect.getFaultSurface(1d).getEvenlyDiscritizedListOfLocsOnSurface();
			double avgLat = locs.stream().mapToDouble(L->L.lat).average().getAsDouble();
			double avgLon = locs.stream().mapToDouble(L->L.lon).average().getAsDouble();
			XYTextAnnotation idAnn = new XYTextAnnotation(sect.getSectionId()+"", avgLon, avgLat);
			idAnn.setFont(font);;
			idAnn.setTextAnchor(TextAnchor.CENTER);
			mapMaker.addAnnotation(idAnn);
		}
		
		mapMaker.plot(outputDir, prefix, title);
		
		CPT depthCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance();
		depthCPT = depthCPT.rescale(0d, 50d);
		mapMaker.plotScatterScalars(directionLocs, depths, depthCPT, "Depth (km)");
		mapMaker.plot(outputDir, prefix+"_depth", title);
	}

}
