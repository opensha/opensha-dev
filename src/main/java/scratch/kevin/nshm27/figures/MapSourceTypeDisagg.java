package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.INV_DIR;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

/**
 * Builds maps of the fractional hazard contribution from each NSHM27 tectonic
 * regime. Contributions are evaluated at the IML of the total-model hazard map.
 */
public class MapSourceTypeDisagg {

	private static final boolean MASK_TO_LAND = true;
	private static final double MASK_DELTA = 0.05;
	private static final double[] PERIODS = { 0d, 0.2d, 1d, 5d };
	private static final ReturnPeriods[] RETURN_PERIODS = SolHazardMapCalc.MAP_RPS;
	private static final double CPT_MAX = 100d;
	private static final double CPT_DELTA = 5d;

	private enum SourceType {
		CRUSTAL("crustal", "Crustal", TectonicRegionType.ACTIVE_SHALLOW),
		INTERFACE("interface", "Subduction interface", TectonicRegionType.SUBDUCTION_INTERFACE),
		SLAB("slab", "Subduction intraslab and outer rise", TectonicRegionType.SUBDUCTION_SLAB);

		private final String directoryToken;
		private final String label;
		private final TectonicRegionType trt;

		private SourceType(String directoryToken, String label, TectonicRegionType trt) {
			this.directoryToken = directoryToken;
			this.label = label;
			this.trt = trt;
		}
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			// for Kevin's convenience within eclipse
			System.err.println("No command line arguments, using hardcoded defaults.");
			NSHM27_SeismicityRegions region = NSHM27_SeismicityRegions.AMSAM;
//			NSHM27_SeismicityRegions region = NSHM27_SeismicityRegions.GNMI;
			int[] vs30s = { 760 };
//			int[] vs30s = { 260, 760 };
			for (int vs30 : vs30s)
				plotRegion(region, vs30);
			return;
		}

		Preconditions.checkArgument(args.length >= 2,
				"USAGE: %s <region> <vs30> [<vs30> ...]", MapSourceTypeDisagg.class.getSimpleName());
		NSHM27_SeismicityRegions region = NSHM27_SeismicityRegions.valueOf(args[0]);
		for (int i=1; i<args.length; i++)
			plotRegion(region, Integer.parseInt(args[i]));
	}

	private static void plotRegion(NSHM27_SeismicityRegions seisReg, int vs30) throws IOException {
		String inversionName = NSHM27_PaperPaths.getSolDir(seisReg).getName();
		File combinedDir = new File(INV_DIR, inversionName+"-ba_hazard-vs"+vs30);
		EnumMap<SourceType, File> sourceDirs = new EnumMap<>(SourceType.class);
		for (SourceType type : SourceType.values())
			sourceDirs.put(type, new File(INV_DIR,
					inversionName+"-ba_hazard-"+type.directoryToken+"-vs"+vs30));

		File regionFile = new File(combinedDir, "gridded_region.geojson");
		Preconditions.checkState(regionFile.exists(), "Missing gridded region: %s", regionFile.getAbsolutePath());
		GriddedRegion gridReg = GriddedRegion.fromFeature(Feature.read(regionFile));
		System.out.println("Plotting "+seisReg+" on "+gridReg.getNodeCount()+" nodes at "
				+(float)gridReg.getSpacing()+" degree spacing");

		File reportDir = new File(combinedDir, "hazard_map_disagg");
		Preconditions.checkState(reportDir.exists() || reportDir.mkdir(), "Couldn't create %s", reportDir);
		File resourcesDir = new File(reportDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir(), "Couldn't create %s", resourcesDir);

		GriddedGeoDataSet[][] referenceMaps = loadReferenceMaps(combinedDir, gridReg);
		GriddedGeoDataSet[][] totalRates = loadExceedanceRates(combinedDir, gridReg, referenceMaps);
		EnumMap<SourceType, GriddedGeoDataSet[][]> sourceRates = new EnumMap<>(SourceType.class);
		for (SourceType type : SourceType.values())
			sourceRates.put(type, loadExceedanceRates(sourceDirs.get(type), gridReg, referenceMaps));

		GriddedGeoDataSet landMask = MASK_TO_LAND ? buildLandMask(gridReg) : null;
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setWritePDFs(true);
		mapMaker.setWriteGeoJSON(false);

		double annX = gridReg.getMinLon() + 0.95*(gridReg.getMaxLon() - gridReg.getMinLon());
		double annY1 = gridReg.getMinLat() + 0.95*(gridReg.getMaxLat() - gridReg.getMinLat());
		double annY2 = gridReg.getMinLat() + 0.85*(gridReg.getMaxLat() - gridReg.getMinLat());
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		DecimalFormat periodDF = new DecimalFormat("0.##");

		List<String> lines = new ArrayList<>();
		lines.add("# "+seisReg.getName()+" Source-Type Hazard Contributions, Vs30="+vs30+" m/s");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";

		for (int p=0; p<PERIODS.length; p++) {
			String periodLabel = PERIODS[p] == 0d ? "PGA" : periodDF.format(PERIODS[p])+"s SA";
			String periodPrefix = PERIODS[p] == 0d ? "pga" : periodDF.format(PERIODS[p])+"s";
			for (int r=0; r<RETURN_PERIODS.length; r++) {
				String label = periodLabel+", "+RETURN_PERIODS[r].label;
				String prefix = periodPrefix+"_"+RETURN_PERIODS[r].name();
				lines.add("## "+label);
				lines.add(topLink);
				lines.add("");

				TableBuilder table = MarkdownUtils.tableBuilder();
				List<String> headings = new ArrayList<>();
				for (SourceType type : SourceType.values())
					headings.add(type.label);
				table.addLine(headings);
				table.initNewLine();

				GriddedGeoDataSet sumMap = new GriddedGeoDataSet(gridReg);
				for (SourceType type : SourceType.values()) {
					GriddedGeoDataSet contribution = calcContributionMap(sourceRates.get(type)[p][r],
							totalRates[p][r], landMask, sumMap);
					Stats stats = calcStats(contribution);
					System.out.println(seisReg+", "+label+", "+type+": "+stats);

					CPT cpt = buildCPT(NSHM27_PaperPaths.getColor(type.trt, IncludeBackgroundOption.INCLUDE));
					mapMaker.clearAnnotations();
					addStatsAnnotations(mapMaker, stats, annX, annY1, annY2, annFont);
					mapMaker.plotXYZData(contribution, cpt, "% Contribution, "+label);
					String mapPrefix = type.directoryToken+"_"+prefix;
					mapMaker.plot(resourcesDir, mapPrefix, type.label);
					table.addColumn("![Map]("+resourcesDir.getName()+"/"+mapPrefix+".png)");
				}
				table.finalizeLine();
				lines.addAll(table.build());
				Stats sumStats = calcStats(sumMap);
				lines.add("");
				lines.add("Contribution sum range: "+sumStats.rangeString()+"; average: "
						+formatPercent(sumStats.average)+".");
				lines.add("");
				System.out.println(seisReg+", "+label+", SUM: "+sumStats);
			}
		}

		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "");
		lines.add(tocIndex, "## Table Of Contents");
		MarkdownUtils.writeReadmeAndHTML(lines, reportDir);
		System.out.println("Wrote "+reportDir.getAbsolutePath());
	}

	private static GriddedGeoDataSet[][] loadReferenceMaps(File combinedDir, GriddedRegion gridReg)
			throws IOException {
		GriddedGeoDataSet[][] maps = new GriddedGeoDataSet[PERIODS.length][RETURN_PERIODS.length];
		File zipFile = new File(combinedDir, "results_hazard_INCLUDE.zip");
		Preconditions.checkState(zipFile.exists(), "Missing reference hazard archive: %s", zipFile);
		try (ZipFile zip = new ZipFile(zipFile)) {
			for (int p=0; p<PERIODS.length; p++) {
				for (int r=0; r<RETURN_PERIODS.length; r++) {
					String entryName = MPJ_LogicTreeHazardCalc.mapPrefix(PERIODS[p], RETURN_PERIODS[r])+".txt";
					ZipEntry entry = zip.getEntry(entryName);
					Preconditions.checkNotNull(entry, "%s doesn't exist in %s", entryName, zipFile);
					maps[p][r] = loadMap(zip, entry, gridReg);
				}
			}
		}
		return maps;
	}

	private static GriddedGeoDataSet loadMap(ZipFile zip, ZipEntry entry, GriddedRegion gridReg)
			throws IOException {
		GriddedGeoDataSet map = new GriddedGeoDataSet(gridReg, false);
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)))) {
			String line;
			int index = 0;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#"))
					continue;
				StringTokenizer tok = new StringTokenizer(line);
				double lon = Double.parseDouble(tok.nextToken());
				double lat = Double.parseDouble(tok.nextToken());
				double value = Double.parseDouble(tok.nextToken());
				Preconditions.checkState(index < gridReg.getNodeCount(), "Too many values in %s", entry.getName());
				Preconditions.checkState(LocationUtils.areSimilar(new Location(lat, lon), gridReg.getLocation(index)),
						"Location mismatch at index %s in %s", index, entry.getName());
				map.set(index++, value);
			}
			Preconditions.checkState(index == gridReg.getNodeCount(), "Loaded %s/%s values from %s",
					index, gridReg.getNodeCount(), entry.getName());
		}
		return map;
	}

	private static GriddedGeoDataSet[][] loadExceedanceRates(File hazardDir, GriddedRegion gridReg,
			GriddedGeoDataSet[][] referenceMaps) throws IOException {
		Preconditions.checkState(hazardDir.exists(), "Missing hazard directory: %s", hazardDir.getAbsolutePath());
		File curvesArchive = new File(hazardDir, "results_hazard_INCLUDE_curves.zip");
		Preconditions.checkState(curvesArchive.isFile(), "Missing INCLUDE curves archive: %s", curvesArchive);
		GriddedGeoDataSet[][] rates = new GriddedGeoDataSet[PERIODS.length][RETURN_PERIODS.length];
		try (ZipFile zip = new ZipFile(curvesArchive)) {
			for (int p=0; p<PERIODS.length; p++) {
				String curvesName = SolHazardMapCalc.getCSV_FileName("curves", PERIODS[p]);
				CSVFile<String> csv = loadCurvesCSV(zip, curvesName);
				DiscretizedFunc[] curves = SolHazardMapCalc.loadCurvesCSV(csv, gridReg, true);
				for (int r=0; r<RETURN_PERIODS.length; r++) {
					GriddedGeoDataSet map = new GriddedGeoDataSet(gridReg);
					for (int i=0; i<curves.length; i++) {
						if (curves[i] == null) {
							map.set(i, 0d);
						} else {
							double iml = referenceMaps[p][r].get(i);
							double probability = curves[i].getInterpolatedY_inLogXLogYDomain(iml);
							map.set(i, -Math.log1p(-probability));
						}
					}
					rates[p][r] = map;
				}
			}
		}
		return rates;
	}

	private static CSVFile<String> loadCurvesCSV(ZipFile zip, String curvesName) throws IOException {
		ZipEntry entry = zip.getEntry(curvesName);
		if (entry == null)
			entry = zip.getEntry(curvesName+".gz");
		Preconditions.checkNotNull(entry, "Missing %s[.gz] in %s", curvesName, zip.getName());
		InputStream input = zip.getInputStream(entry);
		if (entry.getName().endsWith(".gz"))
			input = new GZIPInputStream(input);
		try (InputStream closeableInput = input) {
			return CSVFile.readStream(closeableInput, true);
		}
	}

	private static GriddedGeoDataSet calcContributionMap(GriddedGeoDataSet rateMap,
			GriddedGeoDataSet totalRateMap, GriddedGeoDataSet mask, GriddedGeoDataSet sumMap) {
		GriddedGeoDataSet contribution = new GriddedGeoDataSet(rateMap.getRegion());
		for (int i=0; i<contribution.size(); i++) {
			if (mask != null && mask.get(i) != 1d) {
				contribution.set(i, Double.NaN);
				sumMap.set(i, Double.NaN);
				continue;
			}
			double totalRate = totalRateMap.get(i);
			double value = totalRate > 0d && Double.isFinite(totalRate)
					? 100d*rateMap.get(i)/totalRate : Double.NaN;
			contribution.set(i, value);
			if (Double.isFinite(value) && Double.isFinite(sumMap.get(i)))
				sumMap.add(i, value);
			else
				sumMap.set(i, Double.NaN);
		}
		return contribution;
	}

	private static CPT buildCPT(Color maxColor) {
		CPT cpt = new CPT(0d, CPT_MAX, Color.WHITE, maxColor.darker());
		cpt.setNanColor(Color.WHITE);
		cpt.setAboveMaxColor(cpt.getMaxColor());
		if (CPT_DELTA > 0d)
			cpt = cpt.asDiscrete(CPT_DELTA, false);
		return cpt;
	}

	private static Stats calcStats(GriddedGeoDataSet map) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		double sum = 0d;
		int count = 0;
		for (int i=0; i<map.size(); i++) {
			double value = map.get(i);
			if (Double.isFinite(value)) {
				min = Math.min(min, value);
				max = Math.max(max, value);
				sum += value;
				count++;
			}
		}
		Preconditions.checkState(count > 0, "Map has no finite values");
		return new Stats(min, max, sum/count, count);
	}

	private static void addStatsAnnotations(GeographicMapMaker mapMaker, Stats stats, double x,
			double y1, double y2, Font font) {
		XYTextAnnotation rangeAnn = new XYTextAnnotation("Range: "+stats.rangeString(), x, y1);
		rangeAnn.setFont(font);
		rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		mapMaker.addAnnotation(rangeAnn);
		XYTextAnnotation averageAnn = new XYTextAnnotation("Average: "+formatPercent(stats.average), x, y2);
		averageAnn.setFont(font);
		averageAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		mapMaker.addAnnotation(averageAnn);
	}

	private static String formatPercent(double value) {
		return value < 0.5 ? "<1%" : Math.round(value)+"%";
	}

	private record Stats(double min, double max, double average, int count) {
		String rangeString() {
			return "["+formatPercent(min)+", "+formatPercent(max)+"]";
		}

		@Override
		public String toString() {
			return "range="+rangeString()+", average="+formatPercent(average)+", count="+count;
		}
	}

	public static GriddedGeoDataSet buildLandMask(GriddedRegion gridReg) {
		GriddedGeoDataSet mask = new GriddedGeoDataSet(gridReg);
		double latDelta = Math.max(MASK_DELTA, gridReg.getLatSpacing()*0.5);
		double lonDelta = Math.max(MASK_DELTA, gridReg.getLonSpacing()*0.5);
		XY_DataSet[] outlines = PoliticalBoundariesData.loadDefaultOutlines(gridReg);
		List<Region> outlineRegions = new ArrayList<>();
		for (XY_DataSet outline : outlines) {
			if (outline.getMaxY() < gridReg.getMinLat()-latDelta
					|| outline.getMinY() > gridReg.getMaxLat()+latDelta
					|| outline.getMaxX() < gridReg.getMinLon()-lonDelta
					|| outline.getMinX() > gridReg.getMaxLon()+lonDelta)
				continue;
			LocationList border = new LocationList();
			for (Point2D point : outline)
				border.add(new Location(point.getY(), point.getX()));
//			System.out.println("Masking for region with "+border.size()+" pts:\n"+border);
			Region region = new Region(border, BorderType.MERCATOR_LINEAR);
			outlineRegions.add(region);
		}
		int numInside = 0;
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			Location center = gridReg.getLocation(i);
			GriddedRegion microGrid = new GriddedRegion(
					new Location(center.lat-latDelta, center.lon-lonDelta),
					new Location(center.lat+latDelta, center.lon+lonDelta),
					MASK_DELTA/10d, center);
			LocationList testLocations = new LocationList(microGrid.getNodeCount()+1);
			testLocations.add(center);
			testLocations.addAll(microGrid.getNodeList());
//			Location[] testLocations = {
//					center,
//					new Location(center.lat+latDelta, center.lon+lonDelta),
//					new Location(center.lat+latDelta, center.lon-lonDelta),
//					new Location(center.lat-latDelta, center.lon+lonDelta),
//					new Location(center.lat-latDelta, center.lon-lonDelta)
//			};
			boolean anyInside = false;
			for (Location location : testLocations) {
				for (Region region : outlineRegions) {
					if (region.contains(location)) {
						anyInside = true;
						break;
					}
				}
				if (anyInside) {
					mask.set(i, 1d);
					numInside++;
					break;
				}
			}
		}
		System.out.println("Land mask kept "+numInside+"/"+gridReg.getNodeCount()+" nodes ("+(float)(100d*mask.getSumZ()/gridReg.getNodeCount())+"%)");
		return mask;
	}
}
