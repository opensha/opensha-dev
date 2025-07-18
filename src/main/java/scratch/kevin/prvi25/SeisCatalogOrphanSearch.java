package scratch.kevin.prvi25;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.CSVReader;
import org.opensha.commons.data.CSVReader.Row;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;

import com.google.common.base.Preconditions;

public class SeisCatalogOrphanSearch {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/prvi25/seis_catalogs");
		
		File dir = new File(mainDir, "2025_07_17");
		boolean filterByCrustalReg = true;
		boolean remapORtoSlab = true;
		File catFile = new File(dir, "pmmx_071725_c2-separatedAveraged_v2.csv");
		
		List<EventRecord> fullCatalog = loadCatalog(catFile);
		Region plotRegion = PRVI25_SeismicityRegions.CRUSTAL.load();
		if (filterByCrustalReg) {
			for (int i=fullCatalog.size(); --i>=0;)
				if (!plotRegion.contains(fullCatalog.get(i).loc))
					fullCatalog.remove(i);
		}
		if (remapORtoSlab) {
			double remappedToMue = 0;
			double remappedToCar = 0;
			
			for (int i=0; i<fullCatalog.size(); i++) {
				EventRecord event = fullCatalog.get(i);
				if (event.probCrust > 0d && event.locStr.equals("or")) {
					if (event.slabReg.equals("mue"))
						remappedToMue += event.probCrust/100d;
					else
						remappedToCar += event.probCrust/100d;
					fullCatalog.set(i, new EventRecord(event.eventID, event.loc, event.year, event.magnitude, 0d,
							event.probInt, event.probCrust+event.probSlab, event.locStr, event.slabReg));
				}
			}
			System.out.println("Remapped "+(float)remappedToCar+" crustal OR events to CAR slab");
			System.out.println("Remapped "+(float)remappedToMue+" crustal OR events to MUE slab");
		}
		
		Map<PRVI25_SeismicityRegions, List<EventRecord>> regRecords = new HashMap<>();
		for (PRVI25_SeismicityRegions reg : PRVI25_SeismicityRegions.values())
			regRecords.put(reg, new ArrayList<>());
		HashMap<EventRecord, Collection<PRVI25_SeismicityRegions>> eventsToRegions = new HashMap<>();

		HashMap<String, EventRecord> idsToEvent = new HashMap<>();
		for (EventRecord rec : fullCatalog) {
			Preconditions.checkState(!(idsToEvent.containsKey(rec.eventID)), "Duplicate event: %s", rec);
			idsToEvent.put(rec.eventID, rec);
			Map<PRVI25_SeismicityRegions, Double> associations = rec.getAssociatedRegions();
			Preconditions.checkState(!associations.isEmpty(), "No mappings for %s", rec);
			for (PRVI25_SeismicityRegions reg : associations.keySet())
				regRecords.get(reg).add(rec);
			eventsToRegions.put(rec, associations.keySet());
		}
		
		System.out.println("Loaded "+eventsToRegions.size()+" events");
		
		System.out.println("Looking for mapping inconsistencies");
		double sumTol = 1e-2;
		int numBadSumTo100 = 0;
		int numBadSumInFiles = 0;
		
		CSVFile<String> badMappings100 = new CSVFile<>(true);
		badMappings100.addLine("EventID", "Depth (km)", "p_crust", "p_int", "p_slab", "p_sum");
		CSVFile<String> badMappingsInFiles = new CSVFile<>(false);
		badMappingsInFiles.addLine("EventID", "Probability mapped (%)", "Probability missing (%)", "Missing region", "Distance to region (km)");

		float[] magThresholds = {0f, 5f, 6f, 7f};
		
		Map<PRVI25_SeismicityRegions, double[]> eventsOutside = new EnumMap<>(PRVI25_SeismicityRegions.class);
		Map<PRVI25_SeismicityRegions, MinMaxAveTracker[]> eventsOutsideDistTrack = new EnumMap<>(PRVI25_SeismicityRegions.class);
		Map<PRVI25_SeismicityRegions, List<EventRecord>> eventsInsideLocs = new EnumMap<>(PRVI25_SeismicityRegions.class);
		Map<PRVI25_SeismicityRegions, List<Double>> eventsInsideProbs = new EnumMap<>(PRVI25_SeismicityRegions.class);
		Map<PRVI25_SeismicityRegions, List<EventRecord>> eventsOutsideLocs = new EnumMap<>(PRVI25_SeismicityRegions.class);
		Map<PRVI25_SeismicityRegions, List<Double>> eventsOutsideProbs = new EnumMap<>(PRVI25_SeismicityRegions.class);
		Map<PRVI25_SeismicityRegions, double[]> regCounts = new HashMap<>();
		for (PRVI25_SeismicityRegions reg : PRVI25_SeismicityRegions.values()) {
			eventsOutside.put(reg, new double[magThresholds.length]);
			regCounts.put(reg, new double[magThresholds.length]);
			eventsOutsideDistTrack.put(reg, new MinMaxAveTracker[magThresholds.length]);
			eventsInsideLocs.put(reg, new ArrayList<>());
			eventsInsideProbs.put(reg, new ArrayList<>());
			eventsOutsideLocs.put(reg, new ArrayList<>());
			eventsOutsideProbs.put(reg, new ArrayList<>());
			for (int m=0; m<magThresholds.length; m++) {
				eventsOutsideDistTrack.get(reg)[m] = new MinMaxAveTracker();
			}
		}
		
		for (EventRecord rec : fullCatalog) {
			Map<PRVI25_SeismicityRegions, Double> mappedRegions = rec.getAssociatedRegions();
			double sumMapped = mappedRegions.values().stream().mapToDouble(D->D).sum();
			double sumInsidePolys = 0d;
			List<PRVI25_SeismicityRegions> insideRegions = new ArrayList<>();
			for (PRVI25_SeismicityRegions region : mappedRegions.keySet()) {
				for (int m=0; m<magThresholds.length; m++)
					if ((float)rec.magnitude >= magThresholds[m])
						regCounts.get(region)[m] += mappedRegions.get(region)/100d;
				if (region.load().contains(rec.loc)) {
					sumInsidePolys += mappedRegions.get(region);
					eventsInsideLocs.get(region).add(rec);
					eventsInsideProbs.get(region).add(mappedRegions.get(region)/100d);
					insideRegions.add(region);
				}
			}
			boolean sumsTo100 = Precision.equals(sumMapped, 100d, sumTol);
			boolean sumsInPolys = Precision.equals(sumMapped, sumInsidePolys, sumTol);
			if (!sumsTo100 || !sumsInPolys) {
//				System.out.println("Problem(s) with "+rec);
				if (!sumsTo100) {
					System.out.println("\tMapped probability doesn't sum to 1: sum="+(float)sumMapped+"\trec: "+rec);
					numBadSumTo100++;
					badMappings100.addLine(rec.eventID, (float)rec.loc.depth+"", (float)rec.probCrust+"", (float)rec.probInt+"", (float)rec.probSlab+"", (float)sumMapped+"");
				}
				if (!sumsInPolys) {
					List<String> line = new ArrayList<>();
					line.add(rec.eventID);
					line.add((float)sumInsidePolys+"");
					line.add((float)(sumMapped - sumInsidePolys)+"");
					for (PRVI25_SeismicityRegions region : mappedRegions.keySet()) {
						if (!insideRegions.contains(region)) {
							line.add(region.getShortName());
							double prob = mappedRegions.get(region)/100d;
							Location surfLoc = new Location(rec.loc.lat, rec.loc.lon);
							double dist = region.load().distanceToLocation(surfLoc);
//							if (dist > 100d)
//								System.out.println("Far field mapping issue;\n\tEvent: "+rec+"\n\t"+region.getShortName()+" is "+(float)dist+" km away");
							line.add((float)dist+"");
							eventsOutsideLocs.get(region).add(rec);
							eventsOutsideProbs.get(region).add(prob);
							for (int m=0; m<magThresholds.length; m++) {
								if ((float)rec.magnitude >= magThresholds[m]) {
									eventsOutside.get(region)[m] += prob;
									eventsOutsideDistTrack.get(region)[m].addValue(dist);
								}
							}
						}
					}
//					System.out.println("\tMapped sum in region files: "+(float)sumMappedInFiles);
					numBadSumInFiles++;
					badMappingsInFiles.addLine(line);
				}
			}
		}
		System.out.println(numBadSumTo100+"/"+eventsToRegions.size()+" events don't sum to 100%");
		System.out.println(numBadSumInFiles+"/"+eventsToRegions.size()+" events are mapped to a region they're not contained in");
		System.out.println("Regional missing event counts:");
		CPT probCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, 1d).trim(0d, 0.9).rescale(0d, 1d);
		DecimalFormat pDF = new DecimalFormat("0.0%");
		DecimalFormat oDF = new DecimalFormat("0.#");
		DecimalFormat twoDF = new DecimalFormat("0.00");
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(plotRegion);
		mapMaker.setWriteGeoJSON(false);
		
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		
		// plot all events
		for (PRVI25_SeismicityRegions reg : PRVI25_SeismicityRegions.values()) {
			System.out.println("Plotting all events for "+reg.getName());
			for (boolean m5 : new boolean[] {false,true}) {
				List<EventRecord> allEvents = new ArrayList<>();
				List<Double> allProbs = new ArrayList<>();

				allEvents.addAll(eventsInsideLocs.get(reg));
				allProbs.addAll(eventsInsideProbs.get(reg));
				allEvents.addAll(eventsOutsideLocs.get(reg));
				allProbs.addAll(eventsOutsideProbs.get(reg));
				
				List<Location> allScatterLocs = new ArrayList<>();
				List<PlotCurveCharacterstics> allScatterChars = new ArrayList<>();
				
				mapMaker.clearInsetRegions();
				
				Region testReg = reg.load();
				if (!testReg.equals(plotRegion))
					mapMaker.plotInsetRegions(List.of(testReg), new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK), null, 0);
				
				addEventProbChars(allEvents, allProbs, m5 ? 5f : Float.NEGATIVE_INFINITY, allScatterLocs, allScatterChars, probCPT, true);
				
				mapMaker.plotScatters(allScatterLocs, allScatterChars, null, probCPT, "Association probability");
				
				String prefix = reg.name();
				if (m5)
					prefix += "_m5";
				mapMaker.plot(dir, prefix, reg.getShortName());
				
				mapMaker.clearAnnotations();
			}
		}
		
//		GriddedGeoDataSet carDepths = null;
//		GriddedGeoDataSet mueDepths = null;
		double carCutDepth = 50d;
//		double mueCutDepth = 40d;
		double mueCutDepth = 50d;
		GriddedGeoDataSet carDepths = PRVI25_GridSourceBuilder.loadSubductionDepths(PRVI25_SeismicityRegions.CAR_INTRASLAB);
		GriddedGeoDataSet mueDepths = PRVI25_GridSourceBuilder.loadSubductionDepths(PRVI25_SeismicityRegions.MUE_INTRASLAB);
		
		mapMaker.clearInsetRegions();
		mapMaker.clearScatters();
		
		double maxDepth = 150d;
		CPT depthCPTbelowCut = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, maxDepth).trim(0d, carCutDepth);
		CPT depthCPTaboveCut = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(0d, maxDepth).trim(carCutDepth, maxDepth);
		CPT depthCPT = new CPT();
		depthCPT.addAll(depthCPTbelowCut);
		depthCPT.addAll(depthCPTaboveCut);
		depthCPT.setBelowMinColor(depthCPT.getMinColor());
		depthCPT.setAboveMaxColor(depthCPT.getMaxColor());
		
		mapMaker.plotInsetRegions(List.of(PRVI25_SeismicityRegions.CAR_INTERFACE.load()),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY), null, 0);
		mapMaker.plotXYZData(carDepths, depthCPT, "Depth (km)");
		mapMaker.plot(dir, "CAR_depths", "Caribbean");
		
		depthCPTbelowCut = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, maxDepth).trim(0d, mueCutDepth);
		depthCPTaboveCut = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(0d, maxDepth).trim(mueCutDepth, maxDepth);
		depthCPT = new CPT();
		depthCPT.addAll(depthCPTbelowCut);
		depthCPT.addAll(depthCPTaboveCut);
		depthCPT.setBelowMinColor(depthCPT.getMinColor());
		depthCPT.setAboveMaxColor(depthCPT.getMaxColor());

		mapMaker.clearInsetRegions();
		mapMaker.plotInsetRegions(List.of(PRVI25_SeismicityRegions.MUE_INTERFACE.load()),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY), null, 0);
		mapMaker.plotXYZData(mueDepths, depthCPT, "Depth (km)");
		mapMaker.plot(dir, "MUE_depths", "Muertos");
		
		mapMaker.clearXYZData();
		
		// plot outside events
		for (PRVI25_SeismicityRegions reg : PRVI25_SeismicityRegions.values()) {
			double missingTot = eventsOutside.get(reg)[0];
			if (missingTot == 0d)
				continue;
			System.out.println("\t"+reg.getShortName()+":\t"+oDF.format(missingTot)+" across "+eventsOutsideLocs.get(reg).size()+" events");
			
			mapMaker.clearInsetRegions();
			
			Region testReg = reg.load();
			if (!testReg.equals(plotRegion))
				mapMaker.plotInsetRegions(List.of(testReg), new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK), null, 0);
			
			List<EventRecord> insideEvents = eventsInsideLocs.get(reg);
			List<EventRecord> outsideEvents = eventsOutsideLocs.get(reg);
			List<Double> outsideProbs = eventsOutsideProbs.get(reg);
			
			GriddedGeoDataSet depths = null;
			if (reg == PRVI25_SeismicityRegions.CAR_INTERFACE || reg == PRVI25_SeismicityRegions.CAR_INTRASLAB)
				depths = carDepths;
			else if (reg == PRVI25_SeismicityRegions.MUE_INTERFACE || reg == PRVI25_SeismicityRegions.MUE_INTRASLAB)
				depths = mueDepths;
			
			for (int m=0; m<magThresholds.length; m++) {
				MinMaxAveTracker track = eventsOutsideDistTrack.get(reg)[m];
				boolean anyForMag = false;
				int numOutsideForMag = 0;
				EventRecord largest = null;
				for (EventRecord event : outsideEvents) {
					if ((float)event.magnitude >= magThresholds[m]) {
						anyForMag = true;
						numOutsideForMag++;
						if (largest == null || event.magnitude > largest.magnitude || (event.magnitude == largest.magnitude && event.year > largest.year))
							largest = event;
					}
				}
				if (!anyForMag) {
					System.out.println("\t\tNo M>="+magThresholds[m]+" events");
					continue;
				}
				double missing = eventsOutside.get(reg)[m];
				double tot = regCounts.get(reg)[m];
				String magStr, magPrefix;
				if (magThresholds[m] > 0f) {
					magStr = "M≥"+oDF.format(magThresholds[m])+" ";
					magPrefix = "_m"+oDF.format(magThresholds[m]);
				} else {
					magStr = "";
					magPrefix = "";
				}
				
				String label = "∑(P)="+oDF.format(missing)+" outside "+magStr+"across "+numOutsideForMag+" events ("
						+pDF.format(missing/tot)+" of "+oDF.format(tot)+" total), largest: M"+oDF.format(largest.magnitude)+" ("+largest.year+")";
				System.out.println("\t\t"+label);
				String distLabel = "Outside horizontal distances: avg="+twoDF.format(track.getAverage())
						+", range=["+twoDF.format(track.getMin())+", "+twoDF.format(track.getMax())+"] (km)";
				System.out.println("\t\t"+distLabel);
				
				List<Location> allScatterLocs = new ArrayList<>();
				List<PlotCurveCharacterstics> allScatterChars = new ArrayList<>();
				
				PlotCurveCharacterstics insideChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 1f, Color.GRAY);
				// first add locs inside
				for (EventRecord event : insideEvents) {
					allScatterLocs.add(event.loc);
					allScatterChars.add(insideChar);
				}
				
				addEventProbChars(outsideEvents, outsideProbs, magThresholds[m], allScatterLocs, allScatterChars, probCPT, false);
				
				mapMaker.plotScatters(allScatterLocs, allScatterChars, null, probCPT, "Association probability");
				
				double minLon = plotRegion.getMinLon();
				double minLat = plotRegion.getMinLat();
				
				List<String> annLabels = new ArrayList<>();
				annLabels.add(label);
				annLabels.add(distLabel);
				
				if (depths != null) {
					for (boolean inside : new boolean[] {true,false}) {
						MinMaxAveTracker depthTrackCatalog = new MinMaxAveTracker();
						MinMaxAveTracker depthTrackMapped = new MinMaxAveTracker();
						
						List<EventRecord> events = inside ? insideEvents : outsideEvents;
						
						for (EventRecord event : events) {
							if (!isFixedDepth(event.loc.depth))
								depthTrackCatalog.addValue(event.loc.depth);
							int mappedIndex = depths.indexOf(event.loc);
							if (mappedIndex >= 0)
								depthTrackMapped.addValue(depths.get(mappedIndex));
						}
						
						String depthLabel = inside ? "Inside" : "Outside";
						depthLabel += " depths: ";
						
						depthLabel += "cat avg="+oDF.format(depthTrackCatalog.getAverage())
								+", ["+oDF.format(depthTrackCatalog.getMin())+", "+oDF.format(depthTrackCatalog.getMax())+"]";
						depthLabel += "; mapped avg="+oDF.format(depthTrackMapped.getAverage())
						+", ["+oDF.format(depthTrackMapped.getMin())+", "+oDF.format(depthTrackMapped.getMax())+"]";
						
						annLabels.add(depthLabel);
					}
				}
				
				double lat = minLat + 0.3;
				
				Collections.reverse(annLabels);
				for (String annLabel : annLabels) {
					XYTextAnnotation ann = new XYTextAnnotation(annLabel, minLon+0.3, lat);
					ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
					ann.setFont(annFont);
					mapMaker.addAnnotation(ann);
					lat += 0.3;
				}
				
				mapMaker.plot(dir, reg.name()+"_outside_events"+magPrefix, reg.getShortName());
				
				mapMaker.clearAnnotations();
			}
		}
		if (numBadSumTo100 > 0)
			badMappings100.writeToFile(new File(dir, "bad_sums_to_100.csv"));
		if (numBadSumInFiles > 0)
			badMappingsInFiles.writeToFile(new File(dir, "bad_regional_poly_mappings.csv"));
		
//		System.out.println("Interface associations below "+cutDepth+" km");
//		double maxBelowDepth = 0d;
//		for (EventRecord event : fullCatalog) {
//			if (event.probInt > 0d && event.loc.depth > cutDepth) {
////				System.out.println("\t"+(float)event.loc.depth+";\t"+event);
//				maxBelowDepth = Math.max(maxBelowDepth, event.loc.depth);
//			}
//		}
//		System.out.println("Lowest was "+(float)maxBelowDepth);
		
		// now do psuedo rate calculations to double check Andy's output
		System.out.println("Pseudo rate calculations");
		float[] rateMags = {5f, 6f, 7f, 7.5f};
		int[] minYears = {1900, 1973};
		int endYear = 2023;
		List<PRVI25_SeismicityRegions> rateEpochs = new ArrayList<>();
		for (PRVI25_SeismicityRegions reg : PRVI25_SeismicityRegions.values())
			rateEpochs.add(reg);
		rateEpochs.add(null);
		for (PRVI25_SeismicityRegions reg : rateEpochs) {
			if (reg == null)
				System.out.println("Full catalog");
			else
				System.out.println(reg.getName());
			List<EventRecord> events = reg == null ? fullCatalog : regRecords.get(reg);
			for (int minYear : minYears) {
				int years = 1 + endYear - minYear;
				System.out.println(minYear+" onward ("+years+" years):");
				for (float minMag : rateMags) {
					int rawCount = 0;
					double count = 0d;
					for (EventRecord event : events) {
						if ((float)event.magnitude >= minMag && event.year >= minYear) {
//						if ((float)event.magnitude > minMag && event.year >= minYear) {
							rawCount++;
							double probRate = reg == null ? 1d : event.getAssociatedRegions().get(reg)/100d;
							count += probRate;
						}
					}
					double rate = count/(double)years;
					System.out.println("\tM>"+oDF.format(minMag)+":\t"+(float)rate
							+" /yr;\t"+twoDF.format(count)+" partial events ("+rawCount+" with prob>0)");
				}
			}
		}
	}
	
	private static boolean isFixedDepth(double depth) {
		return depth == 10d || depth == 15d || depth == 25d || depth == 33d || depth == 35d;
	}
	
	private static List<EventRecord> loadCatalog(File file) throws IOException {
		System.out.println("Loading catalog from "+file.getName());
		
		BufferedInputStream is = new BufferedInputStream(new FileInputStream(file));
		CSVReader csv = new CSVReader(is);
		
		List<EventRecord> recs = new ArrayList<>();
		
		csv.read(); // skip header
		for (Row row : csv) {
			EventRecord rec = new EventRecord(row);
			if (rec.magnitude > 7d)
				System.out.println("\t"+rec);
			recs.add(rec);
		}
		
		csv.close();
		is.close();
		
		System.out.println("Loaded "+recs.size()+" events");
		
		return recs;
	}
	
	private static class EventRecord {
		
		public final String eventID;
		public final Location loc;
		public final int year;
		public final double magnitude;
		public final double probCrust;
		public final double probInt;
		public final double probSlab;
		public final String locStr;
		public final String slabReg;
		
		public EventRecord(Row row) {
			loc = new Location(row.getDouble(0), row.getDouble(1), row.getDouble(2));
			eventID = row.get(27);
			year = row.getInt(18);
			magnitude = row.getDouble(6);
			probCrust = row.getDouble(29);
			probInt = row.getDouble(30);
			probSlab = row.getDouble(31);
			locStr = row.get(32);
			slabReg = row.get(34);
			Preconditions.checkState(slabReg.equals("car") || slabReg.equals("mue"));
		}
		
		public EventRecord(String eventID, Location loc, int year, double magnitude, double probCrust, double probInt,
				double probSlab, String locStr, String slabReg) {
			super();
			this.eventID = eventID;
			this.loc = loc;
			this.year = year;
			this.magnitude = magnitude;
			this.probCrust = probCrust;
			this.probInt = probInt;
			this.probSlab = probSlab;
			this.locStr = locStr;
			this.slabReg = slabReg;
		}

		public Map<PRVI25_SeismicityRegions, Double> getAssociatedRegions() {
			Map<PRVI25_SeismicityRegions, Double> ret = new EnumMap<>(PRVI25_SeismicityRegions.class);
			if (probCrust > 0d)
				ret.put(PRVI25_SeismicityRegions.CRUSTAL, probCrust);
			if (slabReg.equals("car")) {
				if (probInt > 0d)
					ret.put(PRVI25_SeismicityRegions.CAR_INTERFACE, probInt);
				if (probSlab > 0d)
					ret.put(PRVI25_SeismicityRegions.CAR_INTRASLAB, probSlab);
			} else {
				if (probInt > 0d)
					ret.put(PRVI25_SeismicityRegions.MUE_INTERFACE, probInt);
				if (probSlab > 0d)
					ret.put(PRVI25_SeismicityRegions.MUE_INTRASLAB, probSlab);
			}
			return ret;
		}

		@Override
		public int hashCode() {
			return Objects.hash(eventID, loc, magnitude, probCrust, probInt, probSlab, slabReg, year);
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			EventRecord other = (EventRecord) obj;
			return Objects.equals(eventID, other.eventID) && Objects.equals(loc, other.loc)
					&& Double.doubleToLongBits(magnitude) == Double.doubleToLongBits(other.magnitude)
					&& Double.doubleToLongBits(probCrust) == Double.doubleToLongBits(other.probCrust)
					&& Double.doubleToLongBits(probInt) == Double.doubleToLongBits(other.probInt)
					&& Double.doubleToLongBits(probSlab) == Double.doubleToLongBits(other.probSlab)
					&& Objects.equals(slabReg, other.slabReg) && year == other.year;
		}

		@Override
		public String toString() {
			return "EventRecord [eventID=" + eventID + ", loc=" + loc + ", year=" + year + ", magnitude=" + magnitude
					+ ", probCrust=" + probCrust + ", probInt=" + probInt + ", probSlab=" + probSlab + ", slabReg="
					+ slabReg + ", locStr=" + locStr + "]";
		}
	}
	
	private static void addEventProbChars(List<EventRecord> events, List<Double> probs, float magThreshold,
			List<Location> scatterLocs, List<PlotCurveCharacterstics> scatterChars, CPT probCPT, boolean sortByMag) {
		List<ComparablePairing<Double, EventRecord>> sortables = new ArrayList<>();
		Map<EventRecord, Double> eventsToProbs = new HashMap<>(events.size());
		for (int l=0; l<events.size(); l++)
			eventsToProbs.put(events.get(l), probs.get(l));
		for (int l=0; l<events.size(); l++) {
			EventRecord event = events.get(l);
			double sortable = sortByMag ? event.magnitude : probs.get(l);
			if ((float)event.magnitude >= magThreshold) {
				sortables.add(new ComparablePairing<>(sortable, event));
			}
		}
		Collections.sort(sortables);
		
		for (ComparablePairing<Double, EventRecord> pair : sortables) {
			EventRecord event = pair.getData();
			scatterLocs.add(event.loc);
			Color color = probCPT.getColor(eventsToProbs.get(event));
			float thickness = 3f + Math.max(0, 3f*((float)event.magnitude-4f));
			scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, thickness, color));
			if ((float)event.magnitude >= 5f) {
				// add outline
				scatterLocs.add(event.loc);
				scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, thickness, Color.BLACK));
			}
		}
	}

}
