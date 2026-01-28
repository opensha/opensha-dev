package scratch.kevin.nshm26;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.DownDipSubsectionBuilder;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;
import com.google.common.io.Files;

import net.mahdilamb.colormap.Colors;

public class DownDipInterfaceSubSectTests {

	public static void main(String[] args) throws IOException {
		File inFile = new File("/tmp/ker_slab2_dep_10km_contours.xyz");
		String prefix = "ker_slab2";
		Range<Double> depthRange = Range.closed(0d, 60d);
		Range<Double> lonFilter = null;
		Range<Double> latFilter = Range.atLeast(-30d);
		
//		File inFile = new File("/tmp/izu_slab2_dep_10km_contours.xyz");
//		String prefix = "izu_slab2";
//		Range<Double> depthRange = Range.closed(0d, 60d);
//		Range<Double> lonFilter = null;
//		Range<Double> latFilter = Range.atMost(27d);
		
		File outputDir = inFile.getParentFile();
		
		double scaleLength = 20d;
		boolean scaleIsMax = false;
		boolean constantCount = false;
		
		List<FaultTrace> rawContours = loadDepthContours(inFile);
		
		if (lonFilter != null || latFilter != null)
			rawContours = filterContours(rawContours, latFilter, lonFilter);
		
		MinMaxAveTracker latRange = new MinMaxAveTracker();
		MinMaxAveTracker lonRange = new MinMaxAveTracker();
		for (FaultTrace trace : rawContours) {
			for (Location loc : trace) {
				latRange.addValue(loc.lat);
				lonRange.addValue(loc.lon);
			}
		}
		
		Region mapReg = new Region(new Location(latRange.getMin()-1d, lonRange.getMin()-1d),
				new Location(latRange.getMax()+1d, lonRange.getMax()+1d));
		GeographicMapMaker mapMaker = new GeographicMapMaker(mapReg);
		
		mapMaker.setWritePDFs(false);
		
		List<LocationList> lines = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		PlotCurveCharacterstics rawChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK);
		for (FaultTrace raw : rawContours) {
			lines.add(raw);
			chars.add(rawChar);
		}
		mapMaker.plotLines(lines, chars);
		
		mapMaker.plot(outputDir, prefix+"_raw", " ");
		List<List<FaultTrace>> stitchedIndividual = new ArrayList<>();
		List<FaultTrace> depthContours = processContours(rawContours, stitchedIndividual);
//		for (int d=0; d<rawContours.size(); d++)
//			System.out.println("\t"+d+". "+traceStr(depthContours.get(d)));
//		Range<Double> depthRange = Range.closedOpen(0d, 30d);
		
		GeoJSONFaultSection refSect = new GeoJSONFaultSection.Builder(0, "Test Fault", depthContours.get(0))
				.dip(90d).upperDepth(depthRange.lowerEndpoint()).lowerDepth(depthRange.upperEndpoint())
				.rake(90).build();
		
		lines.clear();
		chars.clear();
		PlotCurveCharacterstics processedChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK);
		for (FaultTrace processed : depthContours) {
			lines.add(processed);
			chars.add(processedChar);
		}
		PlotCurveCharacterstics connectorChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_red);
		for (List<FaultTrace> stitched : stitchedIndividual) {
			for (int i=1; i<stitched.size(); i++) {
				FaultTrace prev = stitched.get(i-1);
				FaultTrace cur = stitched.get(i);
				LocationList connector = LocationList.of(prev.last(), cur.first());
				lines.add(connector);
				chars.add(connectorChar);
			}
		}
		mapMaker.plotLines(lines, chars);
		LocationList startEnds = new LocationList();
		List<PlotCurveCharacterstics> startEndChars = new ArrayList<>();
		PlotCurveCharacterstics startChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 2f, Colors.tab_blue);
		PlotCurveCharacterstics endChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, Colors.tab_green);
		for (FaultTrace trace : depthContours) {
			startEnds.add(trace.first());
			startEndChars.add(startChar);
		}
		for (FaultTrace trace : depthContours) {
			startEnds.add(trace.last());
			startEndChars.add(endChar);
		}
		mapMaker.plotScatters(startEnds, startEndChars);
		
		mapMaker.plot(outputDir, prefix+"_contour_stitch", " ");
		mapMaker.clearScatters();
		
		lines.clear();
		chars.clear();
		
		rawChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.LIGHT_GRAY);
		for (FaultTrace raw : rawContours) {
			lines.add(raw);
			chars.add(rawChar);
		}
		
		List<FaultTrace> interpTraces = DownDipSubsectionBuilder.interpolateDepthTraces(
				depthContours, 2, scaleLength, scaleIsMax, depthRange);
		System.out.println("Built "+interpTraces.size()+" interpolated traces");
		for (FaultTrace trace : interpTraces)
			System.out.println("\t"+traceStr(trace));

		PlotCurveCharacterstics interpChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1.2f, Color.BLACK);
		for (FaultTrace interp : interpTraces) {
			lines.add(interp);
			chars.add(interpChar);
		}

		mapMaker.plot(outputDir, prefix+"_contour_interp", " ");
		
		GeoJSONFaultSection[][] subSects;
		if (constantCount)
			subSects = DownDipSubsectionBuilder.buildForConstantCountAlong(
					refSect, interpTraces, 2, scaleLength, 0);
		else
			subSects = DownDipSubsectionBuilder.buildForConstantLength(
					refSect, interpTraces, 2, scaleLength, scaleIsMax, 0);
		
		mapMaker.setFillSurfaces(true);
		
		List<GeoJSONFaultSection> allSects = new ArrayList<>();
		MinMaxAveTracker lenRange = new MinMaxAveTracker();
		MinMaxAveTracker ddwRange = new MinMaxAveTracker();
		MinMaxAveTracker areaRange = new MinMaxAveTracker();
		MinMaxAveTracker magRange = new MinMaxAveTracker();
		MinMaxAveTracker dipRange = new MinMaxAveTracker();
		MinMaxAveTracker perRowRange = new MinMaxAveTracker();
		for (GeoJSONFaultSection[] row : subSects) {
			perRowRange.addValue(row.length);
			for (GeoJSONFaultSection sect : row) {
				double length = sect.getTraceLength();
				double ddw = sect.getOrigDownDipWidth();
				double area = sect.getArea(false)*1e-6;
				lenRange.addValue(length);
				ddwRange.addValue(ddw);
				areaRange.addValue(area);
				magRange.addValue(PRVI25_SubductionScalingRelationships.LOGA_C4p0.getMag(area*1e6, length*1e3, ddw*1e3, ddw*1e3, 90d));
				dipRange.addValue(sect.getAveDip());
				allSects.add(sect);
			}
		}
		mapMaker.setFaultSections(allSects);
		
		int sectMod = 20;
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, sectMod-1d);
		
		mapMaker.plotSectScalars(s -> {return (double)(s.getSubSectionIndex() % sectMod);},
				cpt, null);
//		mapMaker.plotSectScalarsByIndex(I -> {return (double)(I % sectMod);},
//				cpt, "Subsection Index");
		
		lines.clear();
		chars.clear();
		
		rawChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.LIGHT_GRAY);
		for (FaultTrace raw : rawContours) {
			lines.add(raw);
			chars.add(rawChar);
		}
		mapMaker.setPlotSectsOnTop(true);
		
//		System.out.println("First sect:\n"
//				+allSects.get(0).toFeature().toJSON()+"\n");
//		System.out.println("Second sect:\n"
//				+allSects.get(1).toFeature().toJSON()+"\n");
//		System.out.println("Last sect:\n"
//				+allSects.get(allSects.size()-1).toFeature().toJSON()+"\n");
		
		System.out.println("Sub-Sect Lenghts:\t"+lenRange);
		System.out.println("Sub-Sect DDWs:\t"+ddwRange);
		System.out.println("Sub-Sect Areas:\t"+areaRange);
		System.out.println("Sub-Sect Mmin:\t"+magRange);
		System.out.println("Sub-Sect Dip:\t"+dipRange);
		System.out.println("Sub-Sect Row counts:\t"+perRowRange);
		
		mapMaker.plot(outputDir, prefix+"_sub_sects", " ");
	}
	
	private static List<FaultTrace> processContours(List<FaultTrace> contours) {
		return processContours(contours, null);
	}
	
	private static List<FaultTrace> processContours(List<FaultTrace> contours, List<List<FaultTrace>> stitchedIndividual) {
		// group by depth
		List<List<FaultTrace>> depthGrouped = new ArrayList<>();
		double curDepth = Double.NaN;
		List<FaultTrace> curTraces = null;
		for (FaultTrace trace : contours) {
			double depth = trace.first().depth;
			if (curTraces == null || !Precision.equals(depth, curDepth, 0.1)) {
				curTraces = new ArrayList<>();
				curDepth = depth;
				depthGrouped.add(curTraces);
			}
			curTraces.add(trace);
		}
		
		// sort by depth increasing
		Collections.sort(depthGrouped, (o1, o2) -> Double.compare(o1.get(0).first().depth, o2.get(0).first().depth));
		
		// figure out overall strike direction using the middle trace
		List<FaultTrace> middleGroup = depthGrouped.get(depthGrouped.size()/2);
		Location[] furthestPair = null;
		double furthestDist = 0d;
		if (middleGroup.size() == 1) {
			FaultTrace trace = middleGroup.get(0);
			furthestDist = LocationUtils.horzDistanceFast(trace.first(), trace.last());
			furthestPair = new Location[] { trace.first(), trace.last() };
		} else {
			for (int i=0; i<middleGroup.size()-1; i++) {
				FaultTrace trace1 = middleGroup.get(i);
				Location[] locs1 = {trace1.first(), trace1.last()};
				for (int j=i+1; j<middleGroup.size(); j++) {
					FaultTrace trace2 = middleGroup.get(j);
					Location[] locs2 = {trace2.first(), trace2.last()};
					for (Location l1 : locs1) {
						for (Location l2 : locs2) {
							double dist = LocationUtils.horzDistanceFast(l1, l2);
							if (dist > furthestDist) {
								furthestDist = dist;
								furthestPair = new Location[] { l1, l2 };
							}
						}
					}
				}
			}
		}
		// that gives us an overall strike direction, but could be reversed relative to the RHR
		double strike = LocationUtils.azimuth(furthestPair[0], furthestPair[1]);
		double[] downAzimuths = new double[depthGrouped.size()-1];
		Location upperMiddle = avgLocation(depthGrouped.get(0));
		for (int i=1; i<depthGrouped.size(); i++) {
			Location lowerMiddle = avgLocation(depthGrouped.get(i));
			
			double dipAz = LocationUtils.azimuth(upperMiddle, lowerMiddle);
			downAzimuths[i-1] = dipAz;
			upperMiddle = lowerMiddle;
		}
		double dipDir = DataUtils.median(downAzimuths);
		// this handles wrapping and returns a value in [0, 360]
		System.out.println("Detected middle azimuth="+oDF.format(strike)+" and dip-dir azimuth="+oDF.format(dipDir));
		if (shouldReverse(strike+90d, dipDir))
			strike = LocationUtils.azimuth(furthestPair[1], furthestPair[0]);
		System.out.println("Detected middle strike="+oDF.format(strike)+" with dip direction of "
				+oDF.format(dipDir)+" (delta="+FaultUtils.getAbsAngleDiff(strike, dipDir)+")");
		
		// now stitch
		List<FaultTrace> stitched = new ArrayList<>();
		for (List<FaultTrace> depthTraces : depthGrouped) {
			double depth = depthTraces.get(0).first().depth;
			System.out.println(oDF.format(depth)+" km: "+depthTraces.size()+" contours");
			
			if (depthTraces.size() == 1) {
				// simple case
				FaultTrace trace = depthTraces.get(0);
				double az = trace.getStrikeDirection();
				boolean reverse = shouldReverse(strike, az, 70d);
				System.out.println("\taz="+oDF.format(az)+", reverse="+reverse);
				if (reverse)
					trace = getReversed(trace);
				stitched.add(trace);
				continue;
			}
			
			// find the most central trace and work out from there
			Location centerLoc = avgLocation(depthTraces);
			
			int centralIndex = -1;
			double centerDist = Double.POSITIVE_INFINITY;
			
			for (int i=0; i<depthTraces.size(); i++) {
				for (Location loc : depthTraces.get(i)) {
					double dist = LocationUtils.horzDistanceFast(loc, centerLoc);
					if (dist < centerDist) {
						centerDist = dist;
						centralIndex = i;
					}
				}
			}
			
			// now build it out
			List<FaultTrace> candidates = new ArrayList<>(depthTraces);
			FaultTrace central = candidates.remove(centralIndex);
			List<FaultTrace> stitchedList = new ArrayList<>(depthTraces.size());
			stitchedList.add(central);
			Location curFirst = central.first();
			Location curLast = central.last();
			while (!candidates.isEmpty()) {
				// find the closest location to either end
				int closestIndex = -1;
				boolean closestReversed = false;
				boolean closestToFirst = false;
				double minDist = Double.POSITIVE_INFINITY;
				for (int i=0; i<candidates.size(); i++) {
					FaultTrace candidate = candidates.get(i);
					Location first = candidate.first();
					Location last = candidate.last();
					double distFF = LocationUtils.horzDistanceFast(curFirst, first);
					double distFL = LocationUtils.horzDistanceFast(curFirst, last);
					double distLF = LocationUtils.horzDistanceFast(curLast, first);
					double distLL = LocationUtils.horzDistanceFast(curLast, last);
					double dist = Math.min(Math.min(distFF, distFL), Math.min(distLF, distLL));
					if (dist < minDist) {
						minDist = dist;
						closestIndex = i;
						if (distFF == dist || distFL == dist)
							closestToFirst = true;
						if (distFF == dist || distLL == dist)
							closestReversed = true;
					}
				}
				FaultTrace closest = candidates.remove(closestIndex);
				if (closestReversed)
					closest = getReversed(closest);
				if (closestToFirst) {
					stitchedList.add(0, closest);
					curFirst = closest.first();
				} else {
					stitchedList.add(closest);
					curLast = closest.last();
				}
			}
			
			double stitchedAz = LocationUtils.azimuth(curFirst, curLast);
			if (shouldReverse(strike, stitchedAz, 70d)) {
				// flip everything
				Collections.reverse(stitchedList);
				for (int i=0; i<stitchedList.size(); i++)
					stitchedList.set(i, getReversed(stitchedList.get(i)));
			}
			
			if (stitchedIndividual != null)
				stitchedIndividual.add(stitchedList);
			
			// build final trace and print distances
			FaultTrace stitchedTrace = new FaultTrace();
			for (int i=0; i<stitchedList.size(); i++) {
				FaultTrace trace = stitchedList.get(i);
				if (i == 0)
					System.out.println("\t0. "+traceStr(trace));
				else
					System.out.println("\t"+i+". "+oDF.format(LocationUtils.horzDistanceFast(stitchedTrace.last(), trace.first()))
							+" km away;\t"+traceStr(trace));
				stitchedTrace.addAll(trace);
			}
			stitched.add(stitchedTrace);
		}
		
		return stitched;
	}
	
	private static List<FaultTrace> filterContours(List<FaultTrace> contours, Range<Double> latRange, Range<Double> lonRange) {
		List<FaultTrace> ret = new ArrayList<>();
		
		for (FaultTrace trace : contours) {
			int firstInside = -1;
			int lastInside = -1;
			boolean allInside = true;
			boolean noneInside = true;
			for (int i=0; i<trace.size(); i++) {
				Location loc = trace.get(i);
				boolean inside = (latRange == null || latRange.contains(loc.lat))
						&& (lonRange == null || lonRange.contains(loc.lon));
				if (inside) {
					if (firstInside < 0)
						firstInside = i;
					lastInside = i;
					noneInside = false;
				} else {
					allInside = false;
				}
			}
			if (noneInside)
				continue;
			if (allInside) {
				ret.add(trace);
				continue;
			}
			// need to cut it; could interpolate, but this is good enough for now
			FaultTrace cutTrace = new FaultTrace(null, 1+lastInside-firstInside);
			for (int i=firstInside; i<=lastInside; i++)
				cutTrace.add(trace.get(i));
			ret.add(cutTrace);
		}
		return ret;
	}
	
	private static List<FaultTrace> loadDepthContours(File inFile) throws IOException {
		List<FaultTrace> ret = new ArrayList<>();
		FaultTrace curTrace = null;
		for (String line : Files.readLines(inFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith(">")) {
				// new contour
				curTrace = new FaultTrace();
				ret.add(curTrace);
				continue;
			}
			StringTokenizer tok = new StringTokenizer(line);
			Preconditions.checkState(tok.countTokens() == 3);
			double lon = Double.parseDouble(tok.nextToken());
			double lat = Double.parseDouble(tok.nextToken());
			double depth = -Double.parseDouble(tok.nextToken());
			curTrace.add(new Location(lat, lon, depth));;
		}
		Collections.reverse(ret);
		return ret;
	}
	
	private static final DecimalFormat twoDF = new DecimalFormat("0.00");
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static String traceStr(FaultTrace trace) {
		boolean resampled = false;
		if (trace.size() > 5) {
			trace = FaultUtils.resampleTrace(trace, 4);
			resampled = true;
		}
		StringBuilder str = null;
		for (Location loc : trace) {
			if (str == null)
				str = new StringBuilder();
			else
				str.append(" ");
			str.append("[").append(twoDF.format(loc.lat)).append(", ").append(twoDF.format(loc.lon))
					.append(", ").append(oDF.format(loc.depth)).append("]");
		}
		str.append("; strike=").append(oDF.format(LocationUtils.azimuth(trace.first(), trace.last())));
		if (resampled)
			str.append(" (resampled)");
		return str.toString();
	}
	
	private static Location avgLocation(LocationList locs) {
		return avgLocation(List.of(locs));
	}
	
	private static Location avgLocation(List<? extends LocationList> locLists) {
		double sumLat = 0d;
		double sumLon = 0d;
		double sumDepth = 0d;
		int count = 0;
		for (LocationList list : locLists) {
			for (Location loc : list) {
				count++;
				sumLat += loc.lat;
				sumLon += loc.lon;
				sumDepth += loc.depth;
			}
		}
		return new Location(sumLat/(double)count, sumLon/(double)count, sumDepth/(double)count);
	}
	
	private static final double AZ_TOL_DEFAULT = 50d;
	
	private static FaultTrace getOriented(FaultTrace trace, double overallStrike) {
		double traceStrike = trace.getStrikeDirection();
		if (shouldReverse(overallStrike, traceStrike, AZ_TOL_DEFAULT))
			return getReversed(trace);
		return trace;
	}
	
	private static boolean shouldReverse(double refStrike, double testStrike) {
		return shouldReverse(refStrike, testStrike, AZ_TOL_DEFAULT);
	}
	
	private static boolean shouldReverse(double refStrike, double testStrike, double tolerance) {
		double diff = FaultUtils.getAbsAngleDiff(refStrike, testStrike);
		Preconditions.checkState(diff >= 0d && diff <= 180d, "bad diff=%s for (%s - %s)", diff, refStrike, testStrike);
		if (diff <= tolerance)
			return false;
		if (diff > 180d-tolerance)
			return true;
		throw new IllegalStateException("Couldn't determine if the trace should be reversed: ref="
				+oDF.format(refStrike)+", dest="+oDF.format(testStrike)+", diff="+oDF.format(diff));
	}
	
	private static FaultTrace getReversed(FaultTrace trace) {
		FaultTrace reversed = new FaultTrace(trace.getName(), trace.size());
		for (int i=trace.size(); --i>=0;)
			reversed.add(trace.get(i));
		return reversed;
	}

}