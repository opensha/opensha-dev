package scratch.kevin.simulators.multiFault;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;

public class SlipVectorTests {
	
	public static void main(String[] args) throws IOException {
		// fault 1
		FaultTrace trace1 = new FaultTrace("test trace");
		trace1.add(new Location(0d, 0d));
		trace1.add(new Location(1d, 0d));
		double rake1 = 180d;
		double dip1 = 90d;
		
		
		// fault 2
		FaultTrace trace2 = new FaultTrace("test trace");
//		trace2.add(new Location(1.1d, 0.1d));
//		trace2.add(new Location(2.1d, 0.1d));
		trace2.add(new Location(1d, 0d));
		trace2.add(new Location(0.5d, 2d));
		double rake2 = 0d;
		double dip2 = 90d;
		
		
		double upperDepth = 0d;
		double lowerDepth = 10d;
		
		FaultSectionPrefData sect1 = new FaultSectionPrefData();
		sect1.setAveDip(dip1);
		sect1.setDipDirection((float)(90d + trace1.getAveStrike()));
		sect1.setAveRake(rake1);
		sect1.setFaultTrace(trace1);
		sect1.setAveUpperDepth(upperDepth);
		sect1.setAveLowerDepth(lowerDepth);
		
		FaultSectionPrefData sect2 = new FaultSectionPrefData();
		sect2.setAveDip(dip2);
		sect2.setDipDirection((float)(90d + trace2.getAveStrike()));
		sect2.setAveRake(rake2);
		sect2.setFaultTrace(trace2);
		sect2.setAveUpperDepth(upperDepth);
		sect2.setAveLowerDepth(lowerDepth);
		
		plot(sect1, sect2);
		
		plotCompatibility(sect1, sect2);
		
//		for (double dip : dips) {
//			for (double rake : rakes) {
//				FaultSectionPrefData sect = new FaultSectionPrefData();
//				sect.setAveDip(dip);
//				sect.setDipDirection((float)(90d + trace.getAveStrike()));
//				sect.setAveRake(rake);
//				sect.setFaultTrace(trace);
//				sect.setAveUpperDepth(upperDepth);
//				sect.setAveLowerDepth(lowerDepth);
//				
//				System.out.println("Dip: "+dip);
//				System.out.println("Rake: "+rake);
//				
//				plot(sect);
//			}
//		}
	}
	
	private static void plot(FaultSection... sects) {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		boolean first = true;
		
		for (FaultSection sect : sects) {
			RuptureSurface surf = sect.getFaultSurface(1d, false, false);
			
			plotSurf(surf, funcs, chars);
			
			LocationVector hwVect = RSQSimU3RupturePageGen.calcSlipVector(sect, false);
			System.out.println("\tHanging wall motion: "+hwVect);
			LocationVector fwVect = RSQSimU3RupturePageGen.calcSlipVector(sect, true);
			System.out.println("\tFoot wall motion: "+fwVect);
			
			double scale = 25d;
			hwVect = new LocationVector(hwVect.getAzimuth(), scale*hwVect.getHorzDistance(),
					scale*hwVect.getVertDistance());
			fwVect = new LocationVector(fwVect.getAzimuth(), scale*fwVect.getHorzDistance(),
					scale*fwVect.getVertDistance());
			
			double strike = sect.getFaultTrace().getStrikeDirection();
			
			double offset = sect.getFaultTrace().getTraceLength()*0.05;
			
			Location midPt = GriddedSurfaceUtils.getSurfaceMiddleLoc(surf);
			Location fwPt = LocationUtils.location(midPt, new LocationVector(strike-90d, offset, 0d));
			Location hwPt = LocationUtils.location(midPt, new LocationVector(strike+90d, offset, 0d));
			
			Location hwEnd = LocationUtils.location(hwPt, stripVert(hwVect));
//			System.out.println("hw moved "+LocationUtils.linearDistance(midPt, hwEnd));
			Location fwEnd = LocationUtils.location(fwPt, stripVert(fwVect));
//			System.out.println("fw moved "+LocationUtils.linearDistance(midPt, fwEnd));
			
			DefaultXY_DataSet hwLine = new DefaultXY_DataSet();
			hwLine.set(hwPt.getLongitude(), hwPt.getLatitude());
			hwLine.set(hwEnd.getLongitude(), hwEnd.getLatitude());
			funcs.add(hwLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
			if (first)
				hwLine.setName("Hanging Wall Motion");
			hwLine.setInfo(hwVect.toString());
			System.out.println(hwLine);
			
			DefaultXY_DataSet fwLine = new DefaultXY_DataSet();
			fwLine.set(fwPt.getLongitude(), fwPt.getLatitude());
			fwLine.set(fwEnd.getLongitude(), fwEnd.getLatitude());
			funcs.add(fwLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.GREEN));
			if (first)
				fwLine.setName("Footwall Motion");
			fwLine.setInfo(fwVect.toString());
			System.out.println(fwLine);
			
			first = false;
		}
		
		String title;
		if (sects.length == 1) {
			title = "Strike="+df.format(sects[0].getFaultTrace().getAveStrike());
			title += ", Dip="+df.format(sects[0].getAveDip());
			title += ", Rake="+df.format(sects[0].getAveRake());
		} else {
			title = " ";
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Longitude", "Latitude");
		spec.setLegendVisible(true);
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		
		for (XY_DataSet func : funcs) {
			for (Point2D pt : func) {
				latTrack.addValue(pt.getY());
				lonTrack.addValue(pt.getX());
			}
		}
		
		double maxSpan = Math.max(latTrack.getMax()-latTrack.getMin(), lonTrack.getMax() - lonTrack.getMin());
		double midLat = latTrack.getMin() + 0.5*(latTrack.getMax()-latTrack.getMin());
		double midLon = lonTrack.getMin() + 0.5*(lonTrack.getMax()-lonTrack.getMin());
		
		Range yRange = new Range(midLat-0.6*maxSpan, midLat+0.6*maxSpan);
		Range xRange = new Range(midLon-0.6*maxSpan, midLon+0.6*maxSpan);
		
		GraphWindow gw = new GraphWindow(spec, true);
		gw.setAxisRange(xRange, yRange);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}
	
	private static void plotSurf(RuptureSurface surf, List<XY_DataSet> funcs,
			List<PlotCurveCharacterstics> chars) {
		if (surf.getAveDip() != 90d) {
			DefaultXY_DataSet outline = new DefaultXY_DataSet();
			for (Location loc : surf.getEvenlyDiscritizedPerimeter())
				outline.set(loc.getLongitude(), loc.getLatitude());
			
			funcs.add(outline);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		}
		
		DefaultXY_DataSet traceXY = new DefaultXY_DataSet();
		for (Location loc : surf.getUpperEdge())
			traceXY.set(loc.getLongitude(), loc.getLatitude());
		
		funcs.add(traceXY);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
	}
	
	private static void plotCompatibility(FaultSection sect1, FaultSection sect2) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		RuptureSurface surf1 = sect1.getFaultSurface(1d, false, false);
		RuptureSurface surf2 = sect2.getFaultSurface(1d, false, false);
		
		plotSurf(surf1, funcs, chars);
		plotSurf(surf2, funcs, chars);
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		
		for (XY_DataSet func : funcs) {
			for (Point2D pt : func) {
				latTrack.addValue(pt.getY());
				lonTrack.addValue(pt.getX());
			}
		}
		
		double maxSpan = Math.max(latTrack.getMax()-latTrack.getMin(), lonTrack.getMax() - lonTrack.getMin());
		double midLat = latTrack.getMin() + 0.5*(latTrack.getMax()-latTrack.getMin());
		double midLon = lonTrack.getMin() + 0.5*(lonTrack.getMax()-lonTrack.getMin());
		
		Range yRange = new Range(midLat-0.6*maxSpan, midLat+0.6*maxSpan);
		Range xRange = new Range(midLon-0.6*maxSpan, midLon+0.6*maxSpan);

		LocationVector hwVect1 = RSQSimU3RupturePageGen.calcSlipVector(sect1, false);
		LocationVector fwVect1 = RSQSimU3RupturePageGen.calcSlipVector(sect1, true);
		LocationVector hwVect2 = RSQSimU3RupturePageGen.calcSlipVector(sect2, false);
		LocationVector fwVect2 = RSQSimU3RupturePageGen.calcSlipVector(sect2, true);
		
		Region reg = new Region(new Location(yRange.getLowerBound(), xRange.getLowerBound()),
				new Location(yRange.getUpperBound(), xRange.getUpperBound()));
		GriddedRegion gridReg = new GriddedRegion(reg, 0.05, null);
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		double[] vals = new double[xyz.size()];
		for (int i=0; i<xyz.size(); i++) {
			Location loc = xyz.getLocation(i);
			boolean fw1 = RSQSimU3RupturePageGen.isOnFootwall(surf1, loc);
			boolean fw2 = RSQSimU3RupturePageGen.isOnFootwall(surf2, loc);
			Location l2 = LocationUtils.location(loc, fw1 ? fwVect1 : hwVect1);
			l2 = LocationUtils.location(l2, fw2 ? fwVect2 : hwVect2);
			double dist = LocationUtils.linearDistance(loc, l2);
			xyz.set(i, dist);
			vals[i] = dist;
		}
		
		System.out.println("Z Range: "+(float)StatUtils.min(vals)+" "+(float)StatUtils.max(vals));
		System.out.println("Mean: "+(float)StatUtils.mean(vals));
		System.out.println("Median: "+(float)DataUtils.median(vals));
		
		CPT cpt = GMT_CPT_Files.GMT_POLAR.instance().rescale(0d, 2d);
		
		XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, "", "Longitude",
				"Latitude", "Compatibility");
		spec.setXYElems(funcs);
		spec.setXYChars(chars);
		
		GraphWindow xyzWind = new GraphWindow(spec, xRange, yRange);
		xyzWind.setVisible(true);
		xyzWind.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}
	
	private static final DecimalFormat df = new DecimalFormat("0.0");
	
	private static LocationVector stripVert(LocationVector vect) {
		return new LocationVector(vect.getAzimuth(), vect.getHorzDistance(), 0d);
	}

}
