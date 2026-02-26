package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.FIGURES_DIR;
import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.FULL_GRID_REG;
import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.ORIG_SOL_FILE;

import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;


public class SupersamplingGridFigure {

	public static void main(String[] args) throws IOException {
		double middleLat = Math.round(0.5*(FULL_GRID_REG.getMaxGridLat() + FULL_GRID_REG.getMinLat())-0.5)+0.5;
		double middleLon = Math.round(0.5*(FULL_GRID_REG.getMaxGridLon() + FULL_GRID_REG.getMinLon())-0.5)+0.5;
//		Location middleLoc = FULL_GRID_REG.getLocation(FULL_GRID_REG.indexForLocation(new Location(middleLat, middleLon)));\
		System.out.println("Middle: "+(float)middleLat+", "+(float)middleLon);
		
		Range latRange = new Range(middleLat-0.55, middleLat+0.55);
		Range lonRange = new Range(middleLon-0.55, middleLon+0.55);
		
		PlotPreferences prefs = PlotPreferences.getDefaultPrintFigurePrefs();
		HeadlessGraphPanel gp = PlotUtils.initHeadless(prefs);
		prefs.setLegendFontSize(8);
		prefs.setLegendItemGraphicPadding(new RectangleInsets(0, 0, 0, 0));
		prefs.setLegendItemLabelPadding(new RectangleInsets(0, 0, 0, 0));
		
		GridSourceList gridList = FaultSystemSolution.load(ORIG_SOL_FILE).requireModule(GridSourceList.class);
		List<PointSource> sources = new ArrayList<>();
		GriddedSeismicitySettings gridSettings = GriddedSeismicitySettings.DEFAULT.forDistanceCorrection(
				PointSourceDistanceCorrections.NONE.get()).forSupersamplingSettings(GridCellSupersamplingSettings.DEFAULT);
		for (int l=0; l<gridList.getNumLocations(); l++) {
			Location gridLoc = gridList.getLocation(l);
			if (latRange.contains(gridLoc.lat) && lonRange.contains(gridLoc.lon)) {
				sources.add(gridList.getSource(TectonicRegionType.ACTIVE_SHALLOW, l, 1d, null, gridSettings));
			}
		}
		
		for (boolean colocated : new boolean[] {true,false}) {
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			// draw grid lines
			PlotCurveCharacterstics gridChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.GRAY);
			for (double lat : FULL_GRID_REG.getLatNodes()) {
				lat += 0.05; // edges, not centers
				if (latRange.contains(lat)) {
					funcs.add(new DefaultXY_DataSet(lonRange.getLowerBound(), lat, lonRange.getUpperBound(), lat));
					chars.add(gridChar);
				}
			}
			for (double lon : FULL_GRID_REG.getLonNodes()) {
				lon += 0.05; // edges, not centers
				if (lonRange.contains(lon)) {
					funcs.add(new DefaultXY_DataSet(lon, latRange.getLowerBound(), lon, latRange.getUpperBound()));
					chars.add(gridChar);
				}
			}
			DefaultXY_DataSet fakeCell = new DefaultXY_DataSet(-100d, -100d);
			fakeCell.setName("Grid cells");
			funcs.add(fakeCell);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 6f, gridChar.getColor()));
			
			DefaultXY_DataSet gridNodes = new DefaultXY_DataSet();
			for (Location gridLoc : FULL_GRID_REG.getNodeList())
				if (latRange.contains(gridLoc.lat) && lonRange.contains(gridLoc.lon))
					gridNodes.set(gridLoc.lon, gridLoc.lat);
			gridNodes.setName("Grid nodes");
			funcs.add(gridNodes);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.BLACK));
			
			Location siteLoc = new Location(latRange.getCentralValue(), lonRange.getCentralValue());
			if (!colocated)
				siteLoc = new Location(siteLoc.lat+0.05, siteLoc.lon+0.05);
			
			DefaultXY_DataSet ssLocs = new DefaultXY_DataSet();
			for (PointSource source : sources) {
				HashSet<Location> uniqueLocs = new HashSet<>();
				for (ProbEqkRupture rup : source.getForSite(new Site(siteLoc))) {
					PointSurface surf = (PointSurface) rup.getRuptureSurface();
					Location loc = surf.getLocation();
					uniqueLocs.add(new Location(loc.lat, loc.lon));
				}
				for (Location loc : uniqueLocs)
					ssLocs.set(loc.lon, loc.lat);
			}
			ssLocs.setName("Supersampled");
			funcs.add(ssLocs);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 1.25f, Colors.tab_blue));
			
//			DefaultXY_DataSet siteNodes = new DefaultXY_DataSet();
//			for (Location gridLoc : FULL_GRID_REG.getNodeList()) {
//				if (!colocated)
//					gridLoc = new Location(gridLoc.lat+0.05, gridLoc.lon+0.05);
//				if (latRange.contains(gridLoc.lat) && lonRange.contains(gridLoc.lon))
//					siteNodes.set(gridLoc.lon, gridLoc.lat);
//			}
			DefaultXY_DataSet siteXY = new DefaultXY_DataSet(siteLoc.lon, siteLoc.lat);
			siteXY.setName("Site location");
			funcs.add(siteXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_INV_TRIANGLE, 6f, Colors.tab_red));
			funcs.add(new DefaultXY_DataSet(siteLoc.lon, siteLoc.lat));
			chars.add(new PlotCurveCharacterstics(PlotSymbol.INV_TRIANGLE, 6f, Colors.tab_red.darker().darker()));
//			for (PlotSymbol sym : PlotSymbol.values()) {
//				double offset = (1+sym.ordinal())*0.025;
//				Location loc = new Location(siteLoc.lat + offset, siteLoc.lon);
//				funcs.add(new DefaultXY_DataSet(loc.lon, loc.lat));
//				chars.add(new PlotCurveCharacterstics(sym, 3f, Colors.tab_red));
//			}
			
			double[] dists =  { gridSettings.supersamplingSettings.fullDist, gridSettings.supersamplingSettings.borderDist };
//			PlotLineType[] distLTs = { PlotLineType.SHORT_DASHED, PlotLineType.DOTTED };
			PlotLineType[] distLTs = { PlotLineType.SOLID, PlotLineType.SOLID };
			List<XYAnnotation> anns = new ArrayList<>();
			for (int d=0; d<dists.length; d++) {
				DefaultXY_DataSet circle = new DefaultXY_DataSet();
				double maxAngle = 2d*Math.PI;
				int numAngles = 360;
				double deltaAngle = maxAngle/(double)numAngles;
				for (int i=0; i<numAngles; i++) {
					Location loc = LocationUtils.location(siteLoc, i*deltaAngle, dists[d]);
					circle.set(loc.lon, loc.lat);
				}
				circle.set(circle.get(0));
				funcs.add(circle);
				chars.add(new PlotCurveCharacterstics(distLTs[d], 0.25f, Color.BLACK));
				
				XYTextAnnotation ann = new XYTextAnnotation((int)dists[d]+" km", siteLoc.getLongitude(), circle.getMinY()+0.02);
				ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 8));
				ann.setBackgroundPaint(new Color(255, 255, 255, 200));
				anns.add(ann);
			}
			
			String title = colocated ? "Colocated site" : "Offset site";
			String prefix = colocated ? "grid_supersample_colocated" : "grid_supersample_offset";
			
			Preconditions.checkState(funcs.size() == chars.size());
			PlotSpec plot = new PlotSpec(funcs, chars, title, "Latitude", "Longitude");
			plot.setPlotAnnotations(anns);
			plot.setLegendVisible(true);
			
			gp.drawGraphPanel(plot, false, false, lonRange, latRange);
			
			PlotUtils.setGridLinesVisible(gp, false, false);
			PlotUtils.setXTick(gp, 0.2);
			PlotUtils.setYTick(gp, 0.2);
			
			PlotUtils.writePrintPlots(FIGURES_DIR, prefix, gp, 0.5*PlotUtils.DEFAULT_USABLE_PAGE_WIDTH, -1d, 300, true, true, true, false);
		}
	}

}
