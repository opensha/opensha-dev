package scratch.kevin.segmentation.wasatchGeometrical;

import java.awt.Color;
import java.awt.Font;
import java.awt.BasicStroke;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.HashSet;
import java.util.Set;
import java.util.LinkedHashSet;

import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.chart.ui.RectangleAnchor;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.FaultTraceSegmentationGeometry;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.FaultTraceSegmentationGeometryDiagnostics;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.FaultTraceSegmentationGeometry.Connector;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.FaultTraceSegmentationGeometry.Measurements;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

/**
 * Geometry-only Wasatch diagnostics. Usage:
 * OUTPUT_DIR MIN_BEND_DEG [MAX_PERP_DEVIATION_DEG [SUBSECTIONS_GEOJSON]]
 * Angular restriction defaults to positive infinity (none).
 *
 * Without a subsection file, uses FM v3 Wasatch sections, minimum two subsections
 * per parent, maximum length 0.5 DDW. These locally generated IDs are diagnostic
 * IDs, not guaranteed to match a particular inversion rupture set. Supplying the
 * exact inversion's subsection GeoJSON preserves that geometry and its IDs.
 * Inter-parent candidate pairs are retained within 15 km, independent of rupture
 * plausibility. These are measurements to review, not an inferred feature catalog.
 */
public class WasatchSegmentationGeometry {

	public static void main(String[] args) throws Exception {
		if (args.length == 0) {
			// for Kevin's convenience within eclipse
			System.err.println("No command line arguments, using hardcoded defaults; use --help to see available options instead.");
			args = new String[] {
					"/tmp/wasatch_geometry_diagnostics",	// output dir
					"10",	// minimum bend in degrees
			};
		}
		Preconditions.checkArgument(args.length >= 2 && args.length <= 4,
				"Usage: OUTPUT_DIR MIN_BEND_DEG [MAX_PERP_DEVIATION_DEG [SUBSECTIONS_GEOJSON]]");
		File outputDir = new File(args[0]);
		Preconditions.checkState(outputDir.isDirectory() || outputDir.mkdirs(), "Cannot create output directory");
		FaultTraceSegmentationGeometry geometry = new FaultTraceSegmentationGeometry(
				Double.parseDouble(args[1]), args.length > 2 ? Double.parseDouble(args[2]) : Double.POSITIVE_INFINITY);
		List<FaultSection> sections = new ArrayList<>();
		String source;
		if (args.length == 4) {
			for (FaultSection section : GeoJSONFaultReader.readFaultSections(new File(args[3])))
				if (wasatch(section)) sections.add(section);
			source = new File(args[3]).getAbsolutePath();
		} else {
			List<FaultSection> parents = new ArrayList<>();
			for (FaultSection parent : NSHM23_FaultModels.WUS_FM_v3.getFaultSections())
				if (wasatch(parent)) parents.add(parent);
			sections.addAll(SubSectionBuilder.buildSubSects(parents, 2, 0.5, Double.NaN));
			source = "NSHM23 FM v3 Wasatch; minimum 2 subsections/parent, maximum 0.5 DDW; diagnostic IDs";
		}
		Preconditions.checkState(!sections.isEmpty(), "No Wasatch sections found");
		List<Measurements> all = new ArrayList<>(), malad = new ArrayList<>(), provoNephi = new ArrayList<>(),
				virginiaStreet = new ArrayList<>();
		for (int i=0; i<sections.size(); i++) {
			FaultSection first = sections.get(i);
			for (int j=i+1; j<sections.size(); j++) {
				FaultSection second = sections.get(j);
				boolean sameParent = first.getParentSectionId() == second.getParentSectionId();
				// Input subsection lists are expected to retain their parent/along-strike order.
				if (sameParent && j != i+1) continue;
				Measurements m = geometry.measure(first, second);
				if (!sameParent && m.closestTraceConnector.lengthKm > 15d) continue;
				all.add(m);
				if (first.getSectionName().contains("Malad") && second.getSectionName().contains("Malad")) malad.add(m);
				if ((first.getSectionName().contains("Provo") && second.getSectionName().contains("Nephi"))
						|| (first.getSectionName().contains("Nephi") && second.getSectionName().contains("Provo"))) provoNephi.add(m);
				if (first.getSectionName().contains("Virginia Street") || second.getSectionName().contains("Virginia Street"))
					virginiaStreet.add(m);
			}
		}
		write(outputDir, "wasatch", sections, all);
		write(outputDir, "malad", matchingSections(sections, "Malad"), malad);
		write(outputDir, "provo_nephi", matchingSections(sections, "Provo", "Nephi"), provoNephi);
		writeDetailedPairs(outputDir, "provo_nephi", sections, provoNephi);
		write(outputDir, "virginia_street", sections, virginiaStreet);
		writeDetailedPairs(outputDir, "virginia_street", sections, virginiaStreet);
		double maxStrikeDifference = writeStrikeComparison(outputDir, sections);
		long internal = all.stream().filter(m -> m.first.parentID == m.second.parentID).count();
		long bends = all.stream().filter(m -> m.first.parentID == m.second.parentID && m.bendTransitionDefined).count();
		String description = "Source: "+source+"\nSubsections: "+sections.size()+"\n"
				+"Minimum bend angle: "+geometry.getMinBendAngleDegrees()+" degrees (inclusive; zero excluded)\n"
				+"Maximum perpendicular deviation: "+geometry.getMaxPerpendicularDeviationDegrees()+" degrees\n"
				+"Internal subsection boundaries: "+internal+"; bend transitions retained: "+bends+"\n"
				+"Inter-parent candidate pairs: "+(all.size()-internal)+" (closest trace distance <= 15 km)\n"
				+"Geometry only: no gap/stepover feature classification, probabilities, or inversion configuration changes.\n"
				+"Geometry and bends use FaultTrace.getStrikeDirection(); traces must be consistently oriented. Axial differences are included for orientation review.\n"
				+"CSV/GeoJSON include five connector alternatives per candidate; missing connectors are explicit in CSV.\n"
				+"Maps: gray=closest trace, orange=closest ends, purple=exact perpendicular, red=overlap stepover, brown dashed=facing ends.\n";
		description += "Internal boundary labels: angle in degrees; red=retained bend, gray=below threshold.\n";
		description += "White squares with black outlines mark subsection endpoints and shared boundaries.\n";
		description += "GeoJSON contains one trace per subsection with its own properties; pair measurements are attached to connector features.\n";
		description += "Pair maps show subsection endpoint-strike axes (blue/green dashed), their bisector (cyan dashed), and a right-angle marker.\n";
		description += "Virginia Street overview and individual pair maps include its internal boundary and all nearby Wasatch subsection pairs within the 15 km candidate cutoff.\n";
		description += "Gap pair maps decompose the facing-end connector: cyan=along-strike gap, brown dotted=across-strike offset.\n"
				+"Facing ends are selected by projected position, not nearest map-view distance. Interior projection extrema on folded traces are flagged separately.\n"
				+"Facing-end gap/offset are measurements; applying gap and/or offset probabilities remains a separate decision.\n";
		description += "Strike comparison CSV: FaultTrace.getAveStrike() versus getStrikeDirection(); geometry uses getStrikeDirection().\n"
				+"Maximum absolute difference: "+maxStrikeDifference+" degrees.\n";
		Files.writeString(new File(outputDir, "README.txt").toPath(), description, StandardCharsets.UTF_8);
		System.out.println(description);
	}

	private static boolean wasatch(FaultSection section) {
		return section.getSectionName().toLowerCase(Locale.ROOT).contains("wasatch");
	}

	private static List<FaultSection> matchingSections(List<FaultSection> sections, String... matches) {
		List<FaultSection> filtered = new ArrayList<>();
		for (FaultSection section : sections)
			for (String match : matches)
				if (section.getSectionName().contains(match)) { filtered.add(section); break; }
		return filtered;
	}

	private static void write(File outputDir, String prefix, List<FaultSection> sections, List<Measurements> measurements) throws Exception {
		write(outputDir, prefix, sections, measurements, false);
	}

	private static void writeDetailedPairs(File outputDir, String prefix, List<FaultSection> sections,
			List<Measurements> measurements) throws Exception {
		StringBuilder index = new StringBuilder("# "+prefix.replace('_', ' ')+" diagnostics\n\n");
		index.append("[Overview](").append(new File(outputDir, prefix+".png").getAbsolutePath()).append(") | ")
				.append("[Measurements](").append(new File(outputDir, prefix+".csv").getAbsolutePath()).append(") | ")
				.append("[GeoJSON](").append(new File(outputDir, prefix+".geojson").getAbsolutePath()).append(")\n\n")
				.append("Subsection IDs refer to the geometry described in README.txt. Nearby pairs use a 15 km closest-trace cutoff.\n\n")
				.append("| Detailed plot | First subsection | Second subsection | Projected relation |\n")
				.append("| --- | --- | --- | --- |\n");
		for (Measurements m : measurements) {
			String suffix = m.first.parentID == m.second.parentID ? "_bend_pair_"
					: m.exactPerpendicularConnector != null ? "_pair_" : "_gap_pair_";
			String plotPrefix = prefix+suffix+m.first.sectionID+"_"+m.second.sectionID;
			write(outputDir, plotPrefix, sections, Arrays.asList(m),
					Double.isFinite(m.averageStrikeAxisDegrees));
			index.append("| [").append(m.first.sectionID).append(" / ").append(m.second.sectionID).append("](")
					.append(new File(outputDir, plotPrefix+".png").getAbsolutePath()).append(") | ")
					.append(m.first.name).append(" | ").append(m.second.name).append(" | ")
					.append(m.relation).append(" |\n");
		}
		Files.writeString(new File(outputDir, prefix+"_plots.md").toPath(), index.toString(), StandardCharsets.UTF_8);
	}

	private static void write(File outputDir, String prefix, List<FaultSection> sections,
			List<Measurements> measurements, boolean showStrikeGuides) throws Exception {
		if (measurements.isEmpty()) return;
		FaultTraceSegmentationGeometryDiagnostics.write(new File(outputDir, prefix+".csv"),
				new File(outputDir, prefix+".geojson"), measurements);
		Set<Integer> usedIDs = new HashSet<>();
		for (Measurements m : measurements) { usedIDs.add(m.first.sectionID); usedIDs.add(m.second.sectionID); }
		double minLat = Double.POSITIVE_INFINITY, minLon = Double.POSITIVE_INFINITY;
		double maxLat = Double.NEGATIVE_INFINITY, maxLon = Double.NEGATIVE_INFINITY;
		List<FaultSection> usedSections = new ArrayList<>();
		for (FaultSection section : sections) {
			if (!usedIDs.contains(section.getSectionId())) continue;
			usedSections.add(section);
			for (Location location : section.getFaultTrace()) {
				minLat = Math.min(minLat, location.lat); maxLat = Math.max(maxLat, location.lat);
				minLon = Math.min(minLon, location.lon); maxLon = Math.max(maxLon, location.lon);
			}
		}
		GeographicMapMaker map = new GeographicMapMaker(new Region(
				new Location(minLat-0.03, minLon-0.05), new Location(maxLat+0.03, maxLon+0.05)), usedSections);
		map.setWritePDFs(false); map.setWriteGeoJSON(false);
		map.setFillSurfaces(false); map.setSectOutlineChar(null);
		map.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		Set<Location> boundaries = new LinkedHashSet<>();
		for (FaultSection section : usedSections) {
			for (Location endpoint : Arrays.asList(section.getFaultTrace().first(), section.getFaultTrace().last()))
				boundaries.add(new Location(endpoint.lat, endpoint.lon));
		}
		map.plotScatters(new ArrayList<>(boundaries), Color.WHITE);
		map.setScatterSymbol(PlotSymbol.FILLED_SQUARE, 7f, PlotSymbol.SQUARE, Color.BLACK);
		List<LocationList> lines = new ArrayList<>();
		List<PlotCurveCharacterstics> styles = new ArrayList<>();
		Color[] colors = {Color.GRAY, Color.ORANGE, new Color(153, 51, 170), Color.RED, new Color(140, 86, 75)};
		for (Measurements m : measurements) {
			if (m.first.parentID == m.second.parentID && !prefix.equals("wasatch")) {
				Location boundary = m.closestEndConnector.firstLocation;
				Color color = m.bendTransitionDefined ? Color.RED : Color.GRAY;
				map.addAnnotation(new XYShapeAnnotation(new Ellipse2D.Double(boundary.lon-0.0015,
						boundary.lat-0.0015, 0.003, 0.003), new BasicStroke(1f), color, color));
				XYTextAnnotation label = new XYTextAnnotation(String.format(Locale.ROOT, "%.1f deg", m.bendAngleDegrees),
						boundary.lon-0.005, boundary.lat);
				label.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 15));
				label.setPaint(color); label.setTextAnchor(TextAnchor.CENTER_RIGHT);
				map.addAnnotation(label);
			}
			List<Connector> connectors = Arrays.asList(m.closestTraceConnector, m.closestEndConnector,
					m.exactPerpendicularConnector, m.stepoverConnector, m.facingEndConnector);
			for (int i=0; i<connectors.size(); i++) {
				Connector c = connectors.get(i);
				if (c == null || c.lengthKm < 1e-6) continue;
				LocationList line = new LocationList(); line.add(c.firstLocation); line.add(c.secondLocation);
				lines.add(line);
				styles.add(new PlotCurveCharacterstics(i == 4 ? PlotLineType.DASHED : PlotLineType.SOLID,
						i >= 3 ? 3.5f : 2.5f, colors[i]));
			}
		}
		map.plotLines(lines, styles);
		List<PlotCurveCharacterstics> legend = new ArrayList<>();
		List<String> legendLabels = new ArrayList<>(Arrays.asList("Closest trace", "Closest ends", "Exact perpendicular",
				"Overlap stepover", "Facing ends", "Subsection boundary/end"));
		for (int i=0; i<colors.length; i++)
			legend.add(new PlotCurveCharacterstics(i == 4 ? PlotLineType.DASHED : PlotLineType.SOLID,
					i >= 3 ? 3.5f : 2.5f, colors[i]));
		legend.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 7f, Color.BLACK));
		if (showStrikeGuides) {
			Measurements m = measurements.get(0);
			boolean internalBoundary = m.first.parentID == m.second.parentID;
			boolean gapView = m.exactPerpendicularConnector == null && !internalBoundary;
			Connector c = m.exactPerpendicularConnector != null ? m.exactPerpendicularConnector
					: m.facingEndConnector != null ? m.facingEndConnector : m.closestTraceConnector;
			Color blue = new Color(30, 100, 210), green = new Color(20, 140, 60), cyan = new Color(0, 150, 170);
			Location center = LocationUtils.location(c.firstLocation,
					LocationUtils.azimuthRad(c.firstLocation, c.secondLocation), c.horizontalDistanceKm*0.5);
			double guideHalfLength = Math.max(1.5, Math.min(3d, c.lengthKm*0.3));
			double[] strikes = {m.first.strikeDegrees, m.second.strikeDegrees, m.averageStrikeAxisDegrees};
			Location[] centers = {c.firstLocation, c.secondLocation, center};
			Color[] guideColors = {blue, green, cyan};
			String[] labels = {"Section "+m.first.sectionID+" endpoint strike", "Section "+m.second.sectionID+" endpoint strike", "Bisector"};
			for (int i=0; i<3; i++) {
				LocationList guide = new LocationList();
				guide.add(LocationUtils.location(centers[i], Math.toRadians(strikes[i]+180d), guideHalfLength));
				guide.add(LocationUtils.location(centers[i], Math.toRadians(strikes[i]), guideHalfLength));
				PlotCurveCharacterstics style = new PlotCurveCharacterstics(PlotLineType.DASHED, 2.5f, guideColors[i]);
				lines.add(guide); styles.add(style); legend.add(style);
				legendLabels.add(String.format(Locale.ROOT, "%s: %.1f deg", labels[i], strikes[i]%180d));
			}
			if (gapView) {
				addGapComponents(map, m, lines, styles, legend, legendLabels);
			} else if (!internalBoundary) {
				// Mark the right angle between the purple connector and the cyan bisector.
				double cornerSize = Math.min(0.6, c.lengthKm*0.12);
				double strikeRadians = Math.toRadians(m.averageStrikeAxisDegrees);
				Location along = LocationUtils.location(center, strikeRadians, cornerSize);
				Location across = LocationUtils.location(center, strikeRadians+Math.PI/2d, cornerSize);
				Location corner = LocationUtils.location(along, strikeRadians+Math.PI/2d, cornerSize);
				LocationList rightAngle = new LocationList(); rightAngle.add(along); rightAngle.add(corner); rightAngle.add(across);
				lines.add(rightAngle); styles.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, cyan));
			}
			map.plotLines(lines, styles);
		}
		map.setCustomLegendItems(legendLabels, legend);
		if (showStrikeGuides)
			map.setLegendInset(false);
		else
			map.setLegendInset(RectangleAnchor.TOP_LEFT);
		String title = showStrikeGuides ? "Subsections "+measurements.get(0).first.sectionID+" / "
				+measurements.get(0).second.sectionID+(measurements.get(0).first.parentID == measurements.get(0).second.parentID
						? String.format(Locale.ROOT, ": internal bend %.1f deg", measurements.get(0).bendAngleDegrees)
						: measurements.get(0).exactPerpendicularConnector == null
								? ": facing-end gap and offset" : ": strike axes and perpendicular")
				: prefix.startsWith("provo_nephi") || prefix.startsWith("virginia_street")
						? (prefix.startsWith("virginia_street") ? "Virginia Street" : "Provo-Nephi")+" geometry ("
						+(Double.isInfinite(measurements.get(0).maxPerpendicularDeviationDegrees) ? "no angular limit"
								: "perpendicular tolerance "+measurements.get(0).maxPerpendicularDeviationDegrees+" deg")+")"
				: prefix.replace('_', ' ')+" geometry (bend threshold "+measurements.get(0).minBendAngleDegrees+" deg)";
		if (showStrikeGuides) {
			// Component corners and strike guides can extend beyond the trace bounds.
			for (LocationList line : lines) {
				for (Location location : line) {
					minLat = Math.min(minLat, location.lat); maxLat = Math.max(maxLat, location.lat);
					minLon = Math.min(minLon, location.lon); maxLon = Math.max(maxLon, location.lon);
				}
			}
			map.setRegion(new Region(new Location(minLat-0.03, minLon-0.05), new Location(maxLat+0.03, maxLon+0.05)));
		}
		map.plot(outputDir, prefix, title, 1200);
	}

	/** Draws approximate geographic guides; reported component lengths come from the common planar frame. */
	private static void addGapComponents(GeographicMapMaker map, Measurements m, List<LocationList> lines,
			List<PlotCurveCharacterstics> styles, List<PlotCurveCharacterstics> legend, List<String> legendLabels) {
		Connector c = m.facingEndConnector;
		Color cyan = new Color(0, 150, 170), brown = new Color(140, 86, 75);
		double strike = Math.toRadians(m.averageStrikeAxisDegrees);
		double connectorBearing = LocationUtils.azimuthRad(c.firstLocation, c.secondLocation);
		double alongBearing = Math.cos(connectorBearing-strike) >= 0d ? strike : strike+Math.PI;
		Location corner = LocationUtils.location(c.firstLocation, alongBearing, m.facingEndGapKm);
		Location[] starts = {c.firstLocation, corner}, ends = {corner, c.secondLocation};
		double[] lengths = {m.facingEndGapKm, m.facingEndOffsetKm};
		String[] labels = {"Along-strike gap", "Across-strike offset"};
		Color[] colors = {cyan, brown};
		for (int i=0; i<2; i++) {
			LocationList line = new LocationList(); line.add(starts[i]); line.add(ends[i]); lines.add(line);
			PlotCurveCharacterstics style = new PlotCurveCharacterstics(i == 0 ? PlotLineType.SOLID : PlotLineType.DOTTED,
					2.5f, colors[i]);
			styles.add(style); legend.add(style);
			legendLabels.add(String.format(Locale.ROOT, "%s: %.2f km", labels[i], lengths[i]));
		}
		if (m.facingEndGapKm > 0d && m.facingEndOffsetKm > 0d) {
			double size = Math.min(0.4, Math.min(m.facingEndGapKm, m.facingEndOffsetKm)*0.15);
			double towardFirst = LocationUtils.azimuthRad(corner, c.firstLocation);
			double towardSecond = LocationUtils.azimuthRad(corner, c.secondLocation);
			Location along = LocationUtils.location(corner, towardFirst, size);
			Location across = LocationUtils.location(corner, towardSecond, size);
			Location diagonal = LocationUtils.location(along, towardSecond, size);
			LocationList rightAngle = new LocationList(); rightAngle.add(along); rightAngle.add(diagonal); rightAngle.add(across);
			lines.add(rightAngle); styles.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, cyan));
		}
	}

	private static double writeStrikeComparison(File outputDir, List<FaultSection> sections) throws Exception {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("section_id", "section_name", "length_weighted_strike_deg", "endpoint_direction_deg", "absolute_difference_deg");
		double maxDifference = 0d;
		for (FaultSection section : sections) {
			double average = section.getFaultTrace().getAveStrike();
			double endpoint = section.getFaultTrace().getStrikeDirection();
			double difference = Math.abs(Math.IEEEremainder(average-endpoint, 360d));
			maxDifference = Math.max(maxDifference, difference);
			csv.addLine(""+section.getSectionId(), section.getSectionName(), ""+average, ""+endpoint, ""+difference);
		}
		csv.writeToFile(new File(outputDir, "strike_comparison.csv"));
		return maxDifference;
	}
}
