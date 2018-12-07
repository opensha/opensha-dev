package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;
import scratch.UCERF3.utils.FaultSystemIO;

public class RegionSubtractPlot {
	
	private static final PlotCurveCharacterstics minuendChar = new PlotCurveCharacterstics(
			PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK);
	private static final PlotCurveCharacterstics subtrahendChar = new PlotCurveCharacterstics(
			PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 4f, Color.RED);
	private static final Color resultMinColor = Color.BLUE;
	private static final Color resultMaxColor = Color.GREEN.darker();
	
	private static void plot(Region minuend, Region subtrahend) {
		Region[] result;
		System.out.println("=====================");
		try {
			result = Region.subtract(minuend, subtrahend);
		} catch (Exception e) {
			e.printStackTrace();
			result = null;
		}
		System.out.println("=====================\n");
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (XY_DataSet xy : toXY(minuend, "Minuend")) {
			funcs.add(xy);
			chars.add(minuendChar);
		}
		
		for (XY_DataSet xy : toXY(subtrahend, "Subtrahend")) {
			funcs.add(xy);
			chars.add(subtrahendChar);
		}
		
		if (result != null) {
			CPT cpt = new CPT(0d, result.length > 1 ? result.length-1 : 1, resultMinColor, resultMaxColor);
			for (int i=0; i<result.length; i++) {
				PlotCurveCharacterstics resultChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f,
						PlotSymbol.CIRCLE, 1.5f, cpt.getColor((float)i));
				for (XY_DataSet xy : toXY(result[i], "Result "+i)) {
					funcs.add(xy);
					chars.add(resultChar);
				}
			}
		}
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (XY_DataSet func : funcs) {
			latTrack.addValue(func.getMinY());
			latTrack.addValue(func.getMaxY());
			lonTrack.addValue(func.getMinX());
			lonTrack.addValue(func.getMaxX());
		}
		double centerLat = 0.5*(latTrack.getMin() + latTrack.getMax());
		double centerLon = 0.5*(lonTrack.getMin() + lonTrack.getMax());
		double maxSpan = Math.max(latTrack.getMax() - latTrack.getMin(), lonTrack.getMax() - lonTrack.getMin());
		double centerBuffer = 0.5*maxSpan + 0.2;
		Region plotReg = new Region(new Location(centerLat - centerBuffer, centerLon - centerBuffer),
				new Location(centerLat + centerBuffer, centerLon + centerBuffer));
		
		GriddedRegion gridReg = new GriddedRegion(plotReg, maxSpan/100d, new Location(centerLat, centerLon));
		DefaultXY_DataSet insideXY = new DefaultXY_DataSet();
		for (Location loc : gridReg.getNodeList()) {
			boolean inside = false;
			if (result != null)
				for (Region r : result)
					inside = inside || r.contains(loc);
			if (inside)
				insideXY.set(loc.getLongitude(), loc.getLatitude());
		}
		if (insideXY.size() > 0) {
			funcs.add(0, insideXY);
			chars.add(0, new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, new Color(230, 230, 230)));
		}
		
		Range xRange = new Range(plotReg.getMinLon(), plotReg.getMaxLon());
		Range yRange = new Range(plotReg.getMinLat(), plotReg.getMaxLat());
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Region Subtraction", "Longitude", "Latitude");
		spec.setLegendVisible(true);
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		gw.setAxisRange(xRange, yRange);
	}
	
	private static List<XY_DataSet> toXY(Region region, String name) {
		List<XY_DataSet> ret = new ArrayList<>();
		ret.add(toXY(region.getBorder()));
		ret.get(0).setName(name);
		List<LocationList> interiors = region.getInteriors();
		if (interiors != null)
			for (LocationList interior : interiors)
				ret.add(toXY(interior));
		return ret;
	}
	
	private static XY_DataSet toXY(LocationList locs) {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		for (Location loc : locs)
			xy.set(loc.getLongitude(), loc.getLatitude());
		xy.set(xy.get(0));
		return xy;
	}

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		plot(new Region(new Location(0, 0.1), new Location(2, 1.9)), new Region(new Location(0.8, 0), new Location(1.2, 2)));
		plot(new Region(new Location(0, -0.1), new Location(2, 2.1)), new Region(new Location(0.8, 0), new Location(1.2, 2)));
		plot(new Region(new Location(0, -0.1), new Location(2, 1.9)), new Region(new Location(0.8, 0), new Location(1.2, 2)));
		plot(new Region(new Location(0, 0), new Location(2, 2)), new Region(new Location(-0.1, -0.1), new Location(2.1, 2.1)));
		plot(new Region(new Location(0, 0), new Location(0.8, 2)), new Region(new Location(1.2, 0), new Location(2, 2)));
		
//		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(new File("/home/kevin/git/ucerf3-etas-launcher/inputs/"
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip"));
//		FaultPolyMgr polyManager = FaultPolyMgr.create(rupSet.getFaultSectionDataList(), 12d);
//
//		plot(polyManager.getPoly(1944), polyManager.getPoly(1851));
//		plot(polyManager.getPoly(1949), polyManager.getPoly(324));
//		plot(polyManager.getPoly(324), polyManager.getPoly(1949));
		
		LocationList minuendBorder = new LocationList();
		minuendBorder.add(new Location(37.46166436009728, -122.22562712732683));
		minuendBorder.add(new Location(37.460260316141316, -122.2232239458848));
		minuendBorder.add(new Location(37.45658278241739, -122.22971334120768));
		Region minuend = new Region(minuendBorder, null);
		LocationList subtrahendBorder = new LocationList();
		subtrahendBorder.add(new Location(37.3676151335851, -122.38364795753529));
		subtrahendBorder.add(new Location(37.36810172870189, -122.38412484292463));
		subtrahendBorder.add(new Location(37.411842857831566, -122.42703216149287));
		subtrahendBorder.add(new Location(37.47564754881593, -122.31743787530517));
		subtrahendBorder.add(new Location(37.538367207355925, -122.20669450856569));
		subtrahendBorder.add(new Location(37.494634195855156, -122.16392394398072));
		subtrahendBorder.add(new Location(37.494139494768504, -122.16344056935307));
		subtrahendBorder.add(new Location(37.431419019746954, -122.27411744897063));
		Region subtrahend = new Region(subtrahendBorder, null);
		plot(minuend, subtrahend);
	}

}
