package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.function.Consumer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.comcat.ComcatAccessor;
import org.opensha.commons.data.comcat.ComcatRegion;
import org.opensha.commons.data.comcat.ComcatRegionAdapter;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.latex.LaTeXUtils;
import scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.Models;

public class HazardMapFigures {

	public static void main(String[] args) throws IOException {
		double[] periods = { 0d, 1d };
		ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50 };
		
		Map<ReturnPeriods, String> texRPSuffixes = new HashMap<>();
		texRPSuffixes.put(ReturnPeriods.TWO_IN_50, "");
		texRPSuffixes.put(ReturnPeriods.TEN_IN_50, "TenFifty");
		
		boolean replot = false;
		
		String mapZipName = "results_hazard_INCLUDE.zip";

		File hazardDir = new File(FIGURES_DIR, "hazard_calcs");
		Preconditions.checkState(hazardDir.exists() || hazardDir.mkdir());
		File hazardRawDir = new File(PAPER_DIR.getParentFile(), "hazard_calcs_raw");
		Preconditions.checkState(hazardRawDir.exists() || hazardRawDir.mkdir());
		
		Map<Models, List<Models>> comparisons = new HashMap<>();

		Table<Models, Models, String> compTitles = HashBasedTable.create();
		Table<Models, Models, Boolean> compForceTights = HashBasedTable.create();
		
		Models[] models = Models.values();
		Models prevModel = null;
		for (Models model : models) {
			if (model.getCustomGridLocationAnchor() != null)
				continue;
			if (model != REF_FINITE_MODEL && model.ordinal() <= PROPOSED_DIST_CORR_MODEL.ordinal())
				// this is a distance correction/point surface representation test
				// include a comparison with our reference model
				compAdd(comparisons, model, REF_FINITE_MODEL);
			if (model.ordinal() >= PROPOSED_DIST_CORR_MODEL.ordinal()) {
				// this is a further grid property test
				if (model.ordinal() < PROPOSED_FULL_MODEL.ordinal())
					// compare with proposed full model
					compAdd(comparisons, model, PROPOSED_FULL_MODEL);
				if (prevModel != PROPOSED_FULL_MODEL)
					// also add the incremental change from just this compared to the previous
					compAdd(comparisons, model, prevModel);
			}
			prevModel = model;
		}
		
		// published vs full proposed
		compAdd(comparisons, Models.AS_PUBLISHED, PROPOSED_FULL_MODEL);
		
		// mag threshold only
		compAdd(comparisons, Models.SPINNING_AVG_CENTERED_M5, Models.SPINNING_AVG_CENTERED_M6);
		
		// alt random and realization count tests
		compAdd(comparisons, Models.FINITE_1X_UNCENTERED, Models.FINITE_1X_UNCENTERED_ALT_RAND);
		compTitles.put(Models.FINITE_1X_UNCENTERED, Models.FINITE_1X_UNCENTERED_ALT_RAND, "Virtual faults (1x)");
		compForceTights.put(Models.FINITE_1X_UNCENTERED, Models.FINITE_1X_UNCENTERED_ALT_RAND, true);
		
		compAdd(comparisons, Models.FINITE_2X_UNCENTERED, Models.FINITE_2X_UNCENTERED_ALT_RAND);
		compTitles.put(Models.FINITE_2X_UNCENTERED, Models.FINITE_2X_UNCENTERED_ALT_RAND, "Virtual faults (2x)");
		compForceTights.put(Models.FINITE_2X_UNCENTERED, Models.FINITE_2X_UNCENTERED_ALT_RAND, true);
		
		compAdd(comparisons, Models.FINITE_5X_UNCENTERED, Models.FINITE_5X_UNCENTERED_ALT_RAND);
		compTitles.put(Models.FINITE_5X_UNCENTERED, Models.FINITE_5X_UNCENTERED_ALT_RAND, "Virtual faults (5x)");
		compForceTights.put(Models.FINITE_5X_UNCENTERED, Models.FINITE_5X_UNCENTERED_ALT_RAND, true);
		
		compAdd(comparisons, Models.FINITE_10X_UNCENTERED, Models.FINITE_10X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_20X_UNCENTERED, Models.FINITE_20X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_50X_UNCENTERED, Models.FINITE_50X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_100X_UNCENTERED, Models.FINITE_100X_UNCENTERED_ALT_RAND);
		
		// fixed (openquake) vs 100 random
		compAdd(comparisons, Models.FINITE_FIXED_1X, Models.FINITE_100X_CENTERED);
		compForceTights.put(Models.FINITE_FIXED_1X, Models.FINITE_100X_CENTERED, true);
		compTitles.put(Models.FINITE_FIXED_1X, Models.FINITE_100X_CENTERED, "1-2 fixed strikes");
		compAdd(comparisons, Models.FINITE_FIXED_1X, Models.FINITE_100X_UNCENTERED);
		compForceTights.put(Models.FINITE_FIXED_1X, Models.FINITE_100X_UNCENTERED, true);
		compTitles.put(Models.FINITE_FIXED_1X, Models.FINITE_100X_UNCENTERED, "1-2 fixed strikes");
		compAdd(comparisons, Models.FINITE_FIXED_2X, Models.FINITE_100X_CENTERED);
		compForceTights.put(Models.FINITE_FIXED_2X, Models.FINITE_100X_CENTERED, true);
		compTitles.put(Models.FINITE_FIXED_2X, Models.FINITE_100X_CENTERED, "2-4 fixed strikes");
		compAdd(comparisons, Models.FINITE_FIXED_2X, Models.FINITE_100X_UNCENTERED);
		compForceTights.put(Models.FINITE_FIXED_2X, Models.FINITE_100X_UNCENTERED, true);
		compTitles.put(Models.FINITE_FIXED_2X, Models.FINITE_100X_UNCENTERED, "2-4 fixed strikes");
		compAdd(comparisons, Models.OPENQUAKE_FINITE, Models.FINITE_100X_CENTERED);
		compForceTights.put(Models.OPENQUAKE_FINITE, Models.FINITE_100X_CENTERED, true);
		compTitles.put(Models.OPENQUAKE_FINITE, Models.FINITE_100X_CENTERED, "4-8 fixed strikes");
		compAdd(comparisons, Models.OPENQUAKE_FINITE, Models.FINITE_100X_UNCENTERED);
		compForceTights.put(Models.OPENQUAKE_FINITE, Models.FINITE_100X_UNCENTERED, true);
		compTitles.put(Models.OPENQUAKE_FINITE, Models.FINITE_100X_UNCENTERED, "4-8 fixed strikes");
//		compAdd(comparisons, Models.OPENQUAKE_FINITE_UNCENTERED, Models.FINITE_100X_UNCENTERED);
		
		compAdd(comparisons, PROPOSED_DIST_CORR_MODEL, Models.AS_PUBLISHED);
		compTitles.put(PROPOSED_DIST_CORR_MODEL, Models.AS_PUBLISHED, "Proposed corrections");
		compAdd(comparisons, PROPOSED_FULL_MODEL, Models.AS_PUBLISHED);
		compTitles.put(PROPOSED_FULL_MODEL, Models.AS_PUBLISHED, "Proposed properties & corrections");
		
		compAdd(comparisons, PROPOSED_FULL_MODEL, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M4);
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3, PROPOSED_FULL_MODEL);
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_CORR_M5, PROPOSED_FULL_MODEL);
		
		// incremental comps, table 1
		compAdd(comparisons, Models.SPINNING_AVG_CENTERED_M6, Models.AS_PUBLISHED);
		compTitles.put(Models.SPINNING_AVG_CENTERED_M6, Models.AS_PUBLISHED, "Improved Rrup and hanging wall");
		compForceTights.put(Models.SPINNING_AVG_CENTERED_M6, Models.AS_PUBLISHED, true);
		
		compAdd(comparisons, Models.SPINNING_AVG_CENTERED_M5, Models.SPINNING_AVG_CENTERED_M6);
		compTitles.put(Models.SPINNING_AVG_CENTERED_M5, Models.SPINNING_AVG_CENTERED_M6, "Lower correction M>5");
		compForceTights.put(Models.SPINNING_AVG_CENTERED_M5, Models.SPINNING_AVG_CENTERED_M6, true);
		
		compAdd(comparisons, Models.FINITE_100X_CENTERED, Models.SPINNING_AVG_CENTERED_M5);
		compTitles.put(Models.FINITE_100X_CENTERED, Models.SPINNING_AVG_CENTERED_M5, "Virtual faults, centered");
		compForceTights.put(Models.FINITE_100X_CENTERED, Models.SPINNING_AVG_CENTERED_M5, true);
		
		compAdd(comparisons, Models.FINITE_100X_UNCENTERED, Models.FINITE_100X_CENTERED);
		compTitles.put(Models.FINITE_100X_UNCENTERED, Models.FINITE_100X_CENTERED, "Virtual faults, uncentered");
		compForceTights.put(Models.FINITE_100X_UNCENTERED, Models.FINITE_100X_CENTERED, true);
		
		compAdd(comparisons, PROPOSED_DIST_CORR_MODEL, Models.FINITE_100X_UNCENTERED);
		compTitles.put(PROPOSED_DIST_CORR_MODEL, Models.FINITE_100X_UNCENTERED, "Distribution-based correction");
		compForceTights.put(PROPOSED_DIST_CORR_MODEL, Models.FINITE_100X_UNCENTERED, true);

		// incremental comps, table 3
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, PROPOSED_DIST_CORR_MODEL);
		compTitles.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, PROPOSED_DIST_CORR_MODEL, "Smooth Ztor transition");
		compForceTights.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, PROPOSED_DIST_CORR_MODEL, true);
		
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR);
		compTitles.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, "Leonard (2010) M~L");
		compForceTights.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, true);
		
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN);
		compTitles.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, "Reduced Mmin=3.5");
		compForceTights.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, true);
		
		// supsersampling/offset comps
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS_OFFSET,
				Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_OFFSET);
		compTitles.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS_OFFSET,
				Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_OFFSET, "No supersampling, offset sites");
		compForceTights.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS_OFFSET,
				Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_OFFSET, true);
		
		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS,
				Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5);
		compTitles.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS,
				Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, "No supersampling");
		compForceTights.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS,
				Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, true);
		
		EnumSet<Models> seismicityMapModels = EnumSet.of(Models.FINITE_100X_CENTERED);
		
		// ones that we'll always plot even for 10 in 50
		Map<Models, List<Models>> alwaysPlotComps = new HashMap<>();
		alwaysPlotComps.put(Models.AS_PUBLISHED, new ArrayList<>()); // plot it as is
		compAdd(alwaysPlotComps, PROPOSED_DIST_CORR_MODEL, Models.AS_PUBLISHED);
		compAdd(alwaysPlotComps, PROPOSED_FULL_MODEL, Models.AS_PUBLISHED);
		
		// boolean here is true:1/3 width, false 1/2 width
		Table<Models, Models, Boolean> paperFigs = HashBasedTable.create();
		// figure 1
		paperFigs.put(Models.AS_PUBLISHED, Models.AS_PUBLISHED, true);
		paperFigs.put(PROPOSED_DIST_CORR_MODEL, Models.AS_PUBLISHED, true);
		paperFigs.put(PROPOSED_FULL_MODEL, Models.AS_PUBLISHED, true);
		
		// incremental 1
		paperFigs.put(Models.SPINNING_AVG_CENTERED_M6, Models.AS_PUBLISHED, true);
		paperFigs.put(Models.SPINNING_AVG_CENTERED_M5, Models.SPINNING_AVG_CENTERED_M6, true);
		paperFigs.put(Models.FINITE_100X_CENTERED, Models.SPINNING_AVG_CENTERED_M5, true);
		paperFigs.put(Models.FINITE_100X_UNCENTERED, Models.FINITE_100X_CENTERED, true);
		paperFigs.put(PROPOSED_DIST_CORR_MODEL, Models.FINITE_100X_UNCENTERED, true);
		
		// incremental 2
		paperFigs.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, PROPOSED_DIST_CORR_MODEL, true);
		paperFigs.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR, true);
		paperFigs.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN, true);
		
		paperFigs.put(Models.FINITE_1X_UNCENTERED, Models.FINITE_1X_UNCENTERED_ALT_RAND, true);
		paperFigs.put(Models.FINITE_2X_UNCENTERED, Models.FINITE_2X_UNCENTERED_ALT_RAND, true);
		paperFigs.put(Models.FINITE_5X_UNCENTERED, Models.FINITE_5X_UNCENTERED_ALT_RAND, true);
		paperFigs.put(Models.FINITE_FIXED_1X, Models.FINITE_100X_CENTERED, true);
		paperFigs.put(Models.FINITE_FIXED_2X, Models.FINITE_100X_CENTERED, true);
		paperFigs.put(Models.OPENQUAKE_FINITE, Models.FINITE_100X_CENTERED, true);
		paperFigs.put(Models.FINITE_FIXED_1X, Models.FINITE_100X_UNCENTERED, true);
		paperFigs.put(Models.FINITE_FIXED_2X, Models.FINITE_100X_UNCENTERED, true);
		paperFigs.put(Models.OPENQUAKE_FINITE, Models.FINITE_100X_UNCENTERED, true);
		paperFigs.put(Models.AS_PUBLISHED, Models.FINITE_100X_UNCENTERED, false); // half width
		paperFigs.put(Models.FINITE_100X_CENTERED, Models.FINITE_100X_UNCENTERED, false); // half width
		
		// offset plots
		paperFigs.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS_OFFSET, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_OFFSET, true);
		paperFigs.put(Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5_NO_SS, Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5, true);
		
//		// + improved Rrup and Rx
//		compAdd(comparisons, Models.SPINNING_AVG_M6, Models.AS_PUBLISHED);
//		// + lower magnitude-threshold
//		compAdd(comparisons, Models.SPINNING_AVG, Models.SPINNING_AVG_M6);
//		// + uncentered
//		compAdd(comparisons, Models.SPINNING_AVG_UNCENTERED, Models.SPINNING_AVG);
//
//		// effect of multiple realizations
//		compAdd(comparisons, Models.FINITE_20X_CENTERED, Models.SPINNING_AVG);
//		// effect of uncentering
//		compAdd(comparisons, Models.FINITE_20X_UNCENTERED, Models.FINITE_20X_CENTERED);
//		
//		// new model vs higher resolution
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED, Models.FINITE_50X_UNCENTERED);
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED, Models.FINITE_100X_UNCENTERED);
//		
//		// modeling choice comparisons
//		// M3 vs M4
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_M3, Models.SPINNING_DIST_5X_UNCENTERED_M4);
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_M3, Models.SPINNING_DIST_5X_UNCENTERED);
//		// M4 vs M5
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_M4, Models.SPINNING_DIST_5X_UNCENTERED);
//		// alt depth profile (interpolated, not shelved)
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_ALT_DEPTH, Models.SPINNING_DIST_5X_UNCENTERED);
//		// alt WC94 lengths
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_ALT_LENGTH, Models.SPINNING_DIST_5X_UNCENTERED);
//		// no SS
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_NO_SS, Models.SPINNING_DIST_5X_UNCENTERED);
		
		ZipFile[] modelZips = new ZipFile[models.length];
		ZipFile[] modelZoomZips = new ZipFile[models.length];
		for (int i=0; i<models.length; i++) {
			Models model = models[i];
			File mapFile = new File(model.getMapDir(), mapZipName);
			if (mapFile.exists())
				modelZips[i] = new ZipFile(mapFile);
			File zoomFile = new File(model.getZoomMapDir(), mapZipName);
			if (zoomFile.exists())
				modelZoomZips[i] = new ZipFile(zoomFile);
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(FULL_GRID_REG);
		GeographicMapMaker zoomMapMaker = new GeographicMapMaker(ZOOM_GRID_REG);
		
		List<? extends FaultSection> allSects = NSHM23_FaultModels.WUS_FM_v3.getFaultSections();

		double thirdWidth = PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/3d;
		double defaultWidth = PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d;
		PlotPreferences plotPrefs = GeographicMapMaker.PLOT_PREFS_PRINT_DEFAULT.clone();
		plotPrefs.setSizeScalar(1d/3d);
		plotPrefs.setPlotLabelFontSize(10);
		PlotPreferences thirdPlotPrefs = plotPrefs.clone();
		thirdPlotPrefs.setTickLabelFontSize(6);
		thirdPlotPrefs.setAxisLabelFontSize(10);
		
		for (GeographicMapMaker map : List.of(mapMaker, zoomMapMaker)) {
			map.setFaultSections(allSects);
			map.setSectOutlineChar(null);
			map.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.7f, Color.DARK_GRAY));
			map.setPoliticalBoundaryChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.7f, Color.GRAY));
//			map.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(127, 127, 127, 100)));
			map.setWriteGeoJSON(false);
			map.setPrintPlotPrefs(plotPrefs);
			map.setDefaultPlotWidthInches(defaultWidth);
		}
//		mapMaker.setSectOutlineChar(null);
//		zoomMapMaker.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(127, 127, 127, 100)));
		
		LocationList zoomScatterLocs = new LocationList();
		List<PlotCurveCharacterstics> zoomScatterChars = new ArrayList<>();
		
		if (ZOOM_CITIES != null && !ZOOM_CITIES.isEmpty()) {
			Font cityFont = new Font(Font.SANS_SERIF, Font.BOLD, 10);
			for (String name : ZOOM_CITIES.keySet()) {
				Location loc = ZOOM_CITIES.get(name);
				XYTextAnnotation ann = new XYTextAnnotation("  "+name, loc.lon, loc.lat);
				ann.setFont(cityFont);
				ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
				zoomMapMaker.addAnnotation(ann);
				
				zoomScatterLocs.add(loc);
				zoomScatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_INV_TRIANGLE, 10f, Colors.tab_brown));
				zoomScatterLocs.add(loc);
				zoomScatterChars.add(new PlotCurveCharacterstics(PlotSymbol.INV_TRIANGLE, 10f, Color.BLACK));
			}
		}
		zoomMapMaker.plotScatters(zoomScatterLocs, zoomScatterChars);
		
		LocationList seisZoomScatterLocs = null;
		List<PlotCurveCharacterstics> seisZoomScatterChars = null;
		
		Color trans = new Color(255, 255, 255, 0);
		
		CPT hazardCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		hazardCPT.setLog10(true);
		hazardCPT.setNanColor(trans);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-20, 20);
		pDiffCPT.setNanColor(trans);
		CPT pDiffTightCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-5, 5);
		pDiffTightCPT.setNanColor(trans);
		double maxForTight = 6d;
		
		GriddedGeoDataSet landMask = buildLandMask(FULL_GRID_REG);
		Map<Location, GriddedRegion> customAnchorFullRegs = new HashMap<>();
		Map<Location, GriddedRegion> customAnchorZoomRegs = new HashMap<>();
		Map<Location, GriddedGeoDataSet> customAnchorLandMasks = new HashMap<>();
		for (Models model : models) {
			Location anchor = model.getCustomGridLocationAnchor();
			if (anchor != null && !customAnchorFullRegs.containsKey(anchor)) {
				GriddedRegion fullReg = new GriddedRegion(FULL_GRID_REG, FULL_GRID_REG.getSpacing(), anchor);
				GriddedRegion zoomReg = new GriddedRegion(ZOOM_GRID_REG, ZOOM_GRID_REG.getSpacing(), anchor);
				customAnchorFullRegs.put(anchor, fullReg);
				customAnchorZoomRegs.put(anchor, zoomReg);
				customAnchorLandMasks.put(anchor, buildLandMask(fullReg));
			}
		}
		
		DecimalFormat oDF = new DecimalFormat("0.#");
		
		for (int p=0; p<periods.length; p++) {
			String perPrefix, perLabel, perUnits;
			String perTexPrefix;
			if (periods[p] == 0d) {
				perPrefix = "pga";
				perLabel = "PGA";
				perUnits = "(g)";
				perTexPrefix = "PGA";
			} else {
				perPrefix = (float)periods[p]+"s";
				perLabel = oDF.format(periods[p])+"s SA";
				perUnits = "(g)";
				Preconditions.checkState(periods[p] == 1d);
				perTexPrefix = "SAOne";
			}
			
			for (ReturnPeriods rp : rps) {
				GriddedGeoDataSet[] maps = new GriddedGeoDataSet[models.length];
				GriddedGeoDataSet[] zoomMaps = new GriddedGeoDataSet[models.length];
				
				String entryName = "map_"+perPrefix+"_"+rp.name()+".txt";
				
				System.out.println("Loading maps for "+perLabel+", "+rp.label);
				
				for (int i=0; i<models.length; i++) {
					GriddedRegion fullReg = FULL_GRID_REG;
					GriddedRegion zoomReg = ZOOM_GRID_REG;
					GriddedGeoDataSet modelLandMask = landMask;
					Location anchor = models[i].getCustomGridLocationAnchor();
					if (anchor != null) {
						fullReg = customAnchorFullRegs.get(anchor);
						zoomReg = customAnchorZoomRegs.get(anchor);
						modelLandMask = customAnchorLandMasks.get(anchor);
					}
					if (modelZips[i] != null)
						maps[i] = loadXYZ(modelZips[i], fullReg, entryName, modelLandMask);
					if (modelZoomZips[i] != null)
						zoomMaps[i] = loadXYZ(modelZoomZips[i], zoomReg, entryName, null);
				}
				
				String mapLabel = perLabel+" "+perUnits+", "+rp.label;
				String mapDiffLabel = "% change, "+perLabel+", "+rp.label;
				
				for (boolean paper : new boolean[] {false,true}) {
					File subDir = new File(paper ? hazardDir : hazardRawDir, perPrefix+"_"+rp.name());
					Preconditions.checkState(subDir.exists() || subDir.mkdir());
					
					String texRPSuffix = texRPSuffixes.get(rp);
					FileWriter texFW = !paper || texRPSuffix == null ? null : new FileWriter(new File(subDir, "hazard_change_stats.tex"));
					
					File mapsDir = new File(subDir, "maps");
					if (alwaysPlotComps != null && !alwaysPlotComps.isEmpty())
						Preconditions.checkState(mapsDir.exists() || mapsDir.mkdir());
					File mapCompDir = new File(subDir, "comparisons");
					if (alwaysPlotComps != null && !alwaysPlotComps.isEmpty())
						Preconditions.checkState(mapCompDir.exists() || mapCompDir.mkdir());
					
					File zoomMapsDir = new File(subDir, "maps_zoom");
					if (alwaysPlotComps != null && !alwaysPlotComps.isEmpty())
						Preconditions.checkState(zoomMapsDir.exists() || zoomMapsDir.mkdir());
					File zoomMapCompDir = new File(subDir, "comparisons_zoom");
					if (alwaysPlotComps != null && !alwaysPlotComps.isEmpty())
						Preconditions.checkState(zoomMapCompDir.exists() || zoomMapCompDir.mkdir());
					
					File cptsDir = new File(subDir, "cpts");
					Preconditions.checkState(cptsDir.exists() || cptsDir.mkdir());
					
					if (replot || !new File(cptsDir, "hazard_cpt.pdf").exists())
						PlotUtils.writeScaleLegendOnly(cptsDir, "hazard_cpt",
								GeographicMapMaker.buildCPTLegend(hazardCPT, mapLabel, mapMaker.getPrintPlotPrefs()),
								defaultWidth, 300, true, true);
					if (replot || !new File(cptsDir, "pdiff_cpt.pdf").exists())
						PlotUtils.writeScaleLegendOnly(cptsDir, "pdiff_cpt",
								GeographicMapMaker.buildCPTLegend(pDiffCPT, mapDiffLabel, mapMaker.getPrintPlotPrefs()),
								defaultWidth, 300, true, true);
					if (replot || !new File(cptsDir, "pdiff_tight_cpt.pdf").exists())
						PlotUtils.writeScaleLegendOnly(cptsDir, "pdiff_tight_cpt",
								GeographicMapMaker.buildCPTLegend(pDiffTightCPT, mapDiffLabel, mapMaker.getPrintPlotPrefs()),
								defaultWidth, 300, true, true);
					
					CSVFile<String> compCSV = new CSVFile<>(true);
					compCSV.addLine("Model", "Comparison Model", "Mean Change", "Mean Absolute Change", "Median Absolute Change", "Maximum Change");
					
					for (int i=0; i<models.length; i++) {
						Models model = models[i];
						if (maps[i] == null) {
							System.err.println("Skipping "+models[i].getName()+", calculations not found");
							continue;
						}
						boolean plotModelMaps = !paper || paperFigs.contains(model, model);
						boolean doSeis = seismicityMapModels != null && seismicityMapModels.contains(model) && zoomMaps[i] != null;
						boolean doThird = plotModelMaps && paper && paperFigs.get(model, model);
						String myMapLabel = mapLabel;
						if (doThird) {
							mapMaker.setDefaultPlotWidthInches(thirdWidth);
							mapMaker.setPrintPlotPrefs(thirdPlotPrefs);
							zoomMapMaker.setDefaultPlotWidthInches(thirdWidth);
							zoomMapMaker.setPrintPlotPrefs(thirdPlotPrefs);
							myMapLabel = null;
						}
						System.out.println("Plotting maps for "+models[i].getName());
						if (plotModelMaps) {
							if (replot || !new File(mapsDir, model.name()+".pdf").exists()) {
								mapMaker.plotXYZData(maps[i], hazardCPT, myMapLabel);
								mapMaker.plot(mapsDir, model.name(), model.getName());
							}
							if (zoomMaps[i] != null && (replot || !new File(zoomMapsDir, model.name()+".pdf").exists())) {
								zoomMapMaker.plotXYZData(zoomMaps[i], hazardCPT, myMapLabel);
								zoomMapMaker.plot(zoomMapsDir, model.name(), model.getName());
							}
						}
						
						if (plotModelMaps && doSeis && (replot || !new File(zoomMapsDir, model.name()+"_seis.pdf").exists())) {
							if (seisZoomScatterLocs == null) {
								seisZoomScatterLocs = new LocationList(zoomScatterLocs);
								seisZoomScatterChars = new ArrayList<>(zoomScatterChars);
								fetchComcatEvents(seisZoomScatterLocs, seisZoomScatterChars);
							}
							
							zoomMapMaker.plotScatters(seisZoomScatterLocs, seisZoomScatterChars);
							zoomMapMaker.plotXYZData(zoomMaps[i], hazardCPT, myMapLabel);
							zoomMapMaker.plot(zoomMapsDir, model.name()+"_seis", model.getName());
							zoomMapMaker.plotScatters(zoomScatterLocs, zoomScatterChars);
						}
						
						if (doThird) {
							mapMaker.setDefaultPlotWidthInches(defaultWidth);
							mapMaker.setPrintPlotPrefs(plotPrefs);
							zoomMapMaker.setDefaultPlotWidthInches(defaultWidth);
							zoomMapMaker.setPrintPlotPrefs(plotPrefs);
						}
						
						List<Models> comps = comparisons.get(model);
						if (comps != null) {
							for (Models comp : comps) {
								int compIndex = comp.ordinal();
								System.out.println("\tVs "+comp.getName());
								if (maps[compIndex] == null) {
									System.err.println("\tSkipping comparison with "+models[compIndex].getName()+", calculations not found");
									continue;
								}
								GriddedGeoDataSet pDiff = pDiff(maps[i], maps[compIndex]);
								DiffStats stats = new DiffStats(pDiff);
								
								boolean plotCompModelMaps = !paper || paperFigs.contains(model, comp);
								doThird = plotCompModelMaps && paper && paperFigs.get(model, comp);
								
								GriddedGeoDataSet zoomPDiff = null;
								double mean = stats.mean;
								double meanAbs = stats.meanAbs;
								double medianAbs = stats.medianAbs;
								double maxSigned = stats.maxSigned;
								if (zoomMaps[i] != null && zoomMaps[compIndex] != null) {
									zoomPDiff = pDiff(zoomMaps[i], zoomMaps[compIndex]);
									DiffStats zoomStats = new DiffStats(zoomPDiff);
									if (Math.abs(zoomStats.maxSigned) > Math.abs(maxSigned))
										maxSigned = zoomStats.maxSigned;
								}
								
								Boolean tight = compForceTights.get(model, comp);
								if (tight == null)
									tight = compForceTights.get(comp, model);
								if (tight == null)
									tight = Math.abs(maxSigned) <= maxForTight;
								
								CPT cpt = tight ? pDiffTightCPT : pDiffCPT;
								
//								System.out.println("\t\t"+stats);
								System.out.println("\t\tmean="+twoDF.format(mean)
									+"%;\tmeanAbs="+twoDF.format(meanAbs)
									+"%;\tmedianAbs="+twoDF.format(medianAbs)
									+"%;\tmaxSigned="+twoDF.format(maxSigned)+"%");

								String prefix = model.name()+"_vs_"+comp.name();
								String label = mapDiffLabel;
								String title = model.getName();
								if (compTitles.contains(model, comp))
									title = compTitles.get(model, comp);
								
								if (doThird) {
									mapMaker.setDefaultPlotWidthInches(thirdWidth);
									mapMaker.setPrintPlotPrefs(thirdPlotPrefs);
									zoomMapMaker.setDefaultPlotWidthInches(thirdWidth);
									zoomMapMaker.setPrintPlotPrefs(thirdPlotPrefs);
									label = null;
								}
								
								if (plotCompModelMaps) {
//									String label = model.getName()+" vs "+comp.getName()+", "+mapDiffLabel;
									if (replot || !new File(mapCompDir, prefix+".pdf").exists()) {
										mapMaker.plotXYZData(pDiff, cpt, label);
										mapMaker.plot(mapCompDir, prefix, title);
									}
									if (zoomPDiff != null && (replot || !new File(zoomMapCompDir, prefix+".pdf").exists())) {
										zoomMapMaker.plotXYZData(zoomPDiff, cpt, label);
										zoomMapMaker.plot(zoomMapCompDir, prefix, title);
									}
								}
								
								if (plotCompModelMaps && doSeis && zoomPDiff != null && (replot || !new File(zoomMapCompDir, prefix+"_seis.pdf").exists())) {
									if (seisZoomScatterLocs == null) {
										seisZoomScatterLocs = new LocationList(zoomScatterLocs);
										seisZoomScatterChars = new ArrayList<>(zoomScatterChars);
										fetchComcatEvents(seisZoomScatterLocs, seisZoomScatterChars);
									}
									
									zoomMapMaker.plotScatters(seisZoomScatterLocs, seisZoomScatterChars);
									zoomMapMaker.plotXYZData(zoomPDiff, cpt, label);
									zoomMapMaker.plot(zoomMapCompDir, prefix+"_seis", title);
									zoomMapMaker.plotScatters(zoomScatterLocs, zoomScatterChars);
								}
								
								if (doThird) {
									mapMaker.setDefaultPlotWidthInches(defaultWidth);
									mapMaker.setPrintPlotPrefs(plotPrefs);
									zoomMapMaker.setDefaultPlotWidthInches(defaultWidth);
									zoomMapMaker.setPrintPlotPrefs(plotPrefs);
								}
								
								compCSV.addLine(model.getName(), comp.getName(), twoDF.format(mean)+"%",
										twoDF.format(meanAbs)+"%", twoDF.format(medianAbs)+"%", twoDF.format(maxSigned)+"%");
								
								if (texFW != null) {
									String texPrefix = model.texName+"Vs"+comp.texName+perTexPrefix+texRPSuffix;
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"Mean", noZeroChangeFormat(oneDF, mean)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"AbsMean", noZeroChangeFormat(oneDF, Math.abs(mean))+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MeanAbs", noZeroChangeFormat(oneDF, meanAbs)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MedianAbs", noZeroChangeFormat(oneDF, medianAbs)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"Max", noZeroChangeFormat(oneDF, maxSigned)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"AbsMax", noZeroChangeFormat(oneDF, Math.abs(maxSigned))+"%")+"\n");
									
									if (model == Models.AS_PUBLISHED && (comp == REF_FINITE_MODEL || comp == PROPOSED_FULL_MODEL)) {
										// add extra rounded values
										// always absolute value for these (used in the text)
										texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MeanRounded", noZeroChangeFormat(roundedDF, Math.abs(mean))+"%")+"\n");
										texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MeanAbsRounded", noZeroChangeFormat(roundedDF, meanAbs)+"%")+"\n");
										texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MedianAbsRounded", noZeroChangeFormat(roundedDF, Math.abs(medianAbs))+"%")+"\n");
										texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MaxRounded", noZeroChangeFormat(roundedDF, Math.abs(maxSigned))+"%")+"\n");
									}
								}
							}
						}
					}
					
					compCSV.writeToFile(new File(subDir, "hazard_change_stats.csv"));
					
					if (texFW != null)
						texFW.close();
				}
			}
		}
	}
	
	private static String noZeroChangeFormat(DecimalFormat df, double change) {
		int digits = df.getMaximumFractionDigits();
		Preconditions.checkState(digits >= 0 && digits < 10);
		double minThreshold = 0.5*Math.pow(10, -digits);
		if (Math.abs(change) < minThreshold) {
			if (change < 0)
				return "-<"+df.format(minThreshold);
			return "<"+df.format(minThreshold);
		}
		return df.format(change);
	}
	
	static GriddedGeoDataSet loadXYZ(ZipFile zip, GriddedRegion gridReg, String entryName, GriddedGeoDataSet mask) throws IOException {
		System.out.println("Loading "+entryName+" from "+zip.getName());
		ZipEntry mapEntry = zip.getEntry(entryName);
		InputStream is = zip.getInputStream(mapEntry);
		Preconditions.checkNotNull(is, "IS is null for %s", entryName);
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(mapEntry)));
		String line = bRead.readLine();
		int index = 0;
		while (line != null) {
			line = line.trim();
			if (!line.startsWith("#")) {
				StringTokenizer tok = new StringTokenizer(line);
				double lon = Double.parseDouble(tok.nextToken());
				double lat = Double.parseDouble(tok.nextToken());
				double val = Double.parseDouble(tok.nextToken());
				Location loc = new Location(lat, lon);
				Preconditions.checkState(LocationUtils.areSimilar(loc, gridReg.getLocation(index)));
				xyz.set(index++, val);
			}
			line = bRead.readLine();
		}
		Preconditions.checkState(index == gridReg.getNodeCount());
		
		if (mask != null) {
			Preconditions.checkState(mask.size() == xyz.size());
			for (int i=0; i<mask.size(); i++)
				xyz.set(i, xyz.get(i)*mask.get(i));
		}
		
		return xyz;
	}
	
	private static GriddedGeoDataSet pDiff(GriddedGeoDataSet map, GriddedGeoDataSet refMap) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map.getRegion());
		for (int i=0; i<ret.size(); i++) {
			double val1 = map.get(i);
			double val2 = refMap.get(i);
			if (!Double.isFinite(val1) || !Double.isFinite(val2)) {
				ret.set(i, Double.NaN);
			} else if (val1 == 0d && val2 == 0d) {
				ret.set(i, 0d);	
			} else {
				ret.set(i, 100d*(val1 - val2)/val2);
			}
		}
		return ret;
	}

	private static DecimalFormat roundedDF = new DecimalFormat("0");
	private static DecimalFormat oneDF = new DecimalFormat("0.0");
	private static DecimalFormat twoDF = new DecimalFormat("0.00");
	
	private static class DiffStats {
		
		public final double mean;
		public final double meanAbs;
		public final double medianAbs;
		public final double maxSigned;
		public final int numFinite;
		
		public DiffStats(GriddedGeoDataSet pDiff) {
			double sum = 0d;
			int numFinite = 0;
			double sumAbs = 0d;
			double largest = 0;
			double largestAbs = 0d;
			
			List<Double> absVals = new ArrayList<>(pDiff.size());
			
			for (int i=0; i<pDiff.size(); i++) {
				double val = pDiff.get(i);
				if (Double.isFinite(val)) {
					numFinite++;
					sum += val;
					double absVal = Math.abs(val);
					absVals.add(absVal);
					sumAbs += absVal;
					if (absVal > largestAbs) {
						largestAbs = absVal;
						largest = val;
					}
					
				}
			}
			
			mean = sum/(double)numFinite;
			meanAbs = sumAbs/(double)numFinite;
			medianAbs = DataUtils.median(Doubles.toArray(absVals));
			maxSigned = largest;
			this.numFinite = numFinite;
		}
		
		public String toString() {
			return "mean="+twoDF.format(mean)+"%;\tmeanAbs="+twoDF.format(meanAbs)+"%;\tmaxSigned="+twoDF.format(maxSigned)+"%";
		}
	}
	
	private static void compAdd(Map<Models, List<Models>> comparisons, Models from, Models to) {
		List<Models> comps = comparisons.get(from);
		if (comps == null) {
			comps = new ArrayList<>(1);
			comparisons.put(from, comps);
		}
		if (!comps.contains(to))
			comps.add(to);
	}
	
	private static void fetchComcatEvents(LocationList locs, List<PlotCurveCharacterstics> chars) {
		double minMag = 5d;
		ComcatRegion cReg = new ComcatRegionAdapter(ZOOM_GRID_REG);
		long startTime = new GregorianCalendar(1900, 0, 1).getTimeInMillis();
		long endTime = new GregorianCalendar(2023, 0, 1).getTimeInMillis();
		ObsEqkRupList events = new ComcatAccessor().fetchEventList(
				null, startTime, endTime, -1, 100d, cReg, false, false, minMag);
		for (ObsEqkRupture event : events) {
			locs.add(event.getHypocenterLocation());
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, Colors.tab_purple));
		}
	}

}
