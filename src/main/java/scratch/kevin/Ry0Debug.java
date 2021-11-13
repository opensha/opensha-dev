package scratch.kevin;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;

import com.google.common.base.Preconditions;

public class Ry0Debug {
	
	private static FaultTrace buildTrace(double... vals) {
		Preconditions.checkState(vals.length % 2 == 0);
		FaultTrace trace = new FaultTrace("Trace");
		
		for (int i=0; i<vals.length; i+=2)
			trace.add(new Location(vals[i], vals[i+1]));
		
		return trace;
	}

	public static void main(String[] args) throws IOException {
		List<FaultTrace> traces = new ArrayList<>();
		
		traces.add(buildTrace(34, -118, 34, -117));
		traces.add(buildTrace(34, -118, 34.03, -117.9, 34.03, -117.1, 34, -117));
		traces.add(buildTrace(34, -118, 33.97, -117.9, 34.03, -117.1, 34, -117));
		traces.add(buildTrace(34, -118, 34.25, -117.5, 34, -117));
		traces.add(buildTrace(34, -118, 34.25, -117.25, 34, -117));
		
		for (FaultTrace trace : traces) {
			SummaryStatistics latTrack = new SummaryStatistics();
			SummaryStatistics lonTrack = new SummaryStatistics();
			
			for (Location loc : trace) {
				latTrack.addValue(loc.getLatitude());
				lonTrack.addValue(loc.getLongitude());
			}
			
			double buffer = 2d;
			double spacing = 0.01;
			
			Region reg = new Region(new Location(latTrack.getMin()-buffer, lonTrack.getMin()-buffer),
					new Location(latTrack.getMax()+buffer, lonTrack.getMax()+buffer));
			GriddedRegion gridReg = new GriddedRegion(reg, spacing, null);
			
			GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
			
			for (int i=0; i<xyz.size(); i++) {
				Location siteLoc = gridReg.getLocation(i);
				xyz.set(i, GriddedSurfaceUtils.getDistanceY0(trace, siteLoc));
			}
			
			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0.01d, xyz.getMaxZ());
			cpt.setBelowMinColor(Color.WHITE);
//			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(xyz.getMinZ(), xyz.getMaxZ());
			
			XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, "Ry0 Test", "Longitude", "Latitude", "Ry0 (km)");
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DefaultXY_DataSet traceXY = new DefaultXY_DataSet();
			for (Location loc : trace)
				traceXY.set(loc.getLongitude(), loc.getLatitude());
			
			funcs.add(traceXY);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
			
			spec.setXYElems(funcs);
			spec.setXYChars(chars);
			
			GraphWindow window = new GraphWindow(spec, new Range(reg.getMinLon(), reg.getMaxLon()), new Range(reg.getMinLat(), reg.getMaxLat()));
			window.setVisible(true);
			window.setSize(450, 450);
			window.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		}
	}

}
