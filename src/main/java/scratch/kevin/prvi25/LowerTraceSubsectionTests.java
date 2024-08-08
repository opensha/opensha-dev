package scratch.kevin.prvi25;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.faultSurface.ApproxEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class LowerTraceSubsectionTests {

	public static void main(String[] args) throws IOException {
		FaultTrace upper = new FaultTrace("upper");
		FaultTrace lower = new FaultTrace("upper");
		
		upper.add(new Location(0, 0, 1d));
		lower.add(new Location(-0.5, 0, 40d));
		
		upper.add(new Location(0.15, 0.35, 1.5d));
		lower.add(new Location(-0.42, 0.4, 42d));
		
		upper.add(new Location(0.2, 1, 0d));
		lower.add(new Location(-0.4, 1.05, 40d));
		
		double upperDepth = Double.POSITIVE_INFINITY;
		for (Location loc : upper)
			upperDepth = Math.min(upperDepth, loc.depth);
		double lowerDepth = 0d;
		for (Location loc : lower)
			lowerDepth = Math.max(lowerDepth, loc.depth);
		
		double rake = 90d;
		ApproxEvenlyGriddedSurface approxSurf = new ApproxEvenlyGriddedSurface(
				upper, lower, 1d);
		double dip = approxSurf.getAveDip();
		double dipDir = approxSurf.getAveDipDirection();
		String sectJSON =
				"    {\n"+
				"      \"type\": \"Feature\",\n"+
				"      \"id\": 0,\n"+
				"      \"properties\": {\n"+
				"        \"FaultID\": 0,\n"+
				"        \"FaultName\": \"Test Fault\",\n"+
				"        \"DipDeg\": "+(float)dip+",\n"+
				"        \"DipDir\": "+(float)dipDir+",\n"+
				"        \"Rake\": "+(float)rake+",\n"+
				"        \"LowDepth\": "+(float)lowerDepth+",\n"+
//				"        \"UpDepth\": "+(float)upperDepth+",\n"+
				"        \"SlipRate\": 10\n"+
				"      },\n"+
				"      \"geometry\": {\n"+
				"        \"type\": \"LineString\",\n"+
				"        \"coordinates\": [\n";
		for (int i=0; i<upper.size(); i++) {
			Location loc = upper.get(i);
			sectJSON += 
				"          [\n"+
				"            "+loc.getLongitude()+",\n"+
				"            "+loc.getLatitude()+",\n"+
				"            "+loc.getDepth()+"\n"+
				"          ]";
			if (i < upper.size()-1)
				sectJSON += ",";
			sectJSON += "\n";
		}
		sectJSON +=
				"        ]\n"+
				"      }\n"+
				"    }";
		System.out.println("Input JSON:\n"+sectJSON);
		GeoJSONFaultSection origSect = GeoJSONFaultSection.fromFeature(Feature.fromJSON(sectJSON));
		System.out.println("Parsed JSON:\n"+origSect.toFeature().toJSON());
		Feature feature = origSect.toFeature();
		feature = new Feature((Number)feature.id, new MultiLineString(List.of(upper, lower)), new FeatureProperties(feature.properties));
		GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(feature);
		System.out.println("JSON with lower trace:\n"+sect.toFeature().toJSON());
		Preconditions.checkNotNull(sect.getLowerFaultTrace());
		
		double ddw = approxSurf.getAveWidth();
		double halfDDW = 0.5*ddw;
		
		List<FaultSection> origSectAndSubs = new ArrayList<>();
		origSectAndSubs.add(origSect);
		origSectAndSubs.addAll(origSect.getSubSectionsList(halfDDW, 1, 2));
		GeographicMapMaker mapMaker = new GeographicMapMaker(origSectAndSubs);
		mapMaker.plot(new File("/tmp"), "orig_subsects", "Simplified (no lower trace)");
		plotDepth(new File("/tmp"), "orig_subsects_depth", "Simplified (no lower trace)", origSectAndSubs);
		
		List<FaultSection> sectAndSubs = new ArrayList<>();
		sectAndSubs.add(sect);
		sectAndSubs.addAll(sect.getSubSectionsList(halfDDW, 1, 2));
		mapMaker = new GeographicMapMaker(sectAndSubs);
		mapMaker.plot(new File("/tmp"), "mod_subsects", "Complex (with lower trace)");
		plotDepth(new File("/tmp"), "mod_subsects_depth", "Complex (with lower trace)", sectAndSubs);
		
		for (FaultSection subSect : sectAndSubs) {
			System.out.println(subSect.getSectionId()+". "+subSect.getName());
			System.out.println("\tDip:\t"+(float)subSect.getAveDip());
			System.out.println("\tDip Direction:\t"+(float)subSect.getDipDirection());
			System.out.print("\tLower Trace:");
			for (Location loc : subSect.getLowerFaultTrace())
				System.out.print("\t"+(float)loc.lat+","+(float)loc.lon+","+(float)loc.depth);
			System.out.println();
		}
	}
	
	private static void plotDepth(File outputDir, String prefix, String title, List<FaultSection> sects) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY);
		PlotCurveCharacterstics traceChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK);
		
		for (FaultSection sect : sects) {
			DefaultXY_DataSet outlineXY = new DefaultXY_DataSet();
			for (Location loc : sect.getFaultSurface(1d).getPerimeter())
				outlineXY.set(loc.getLongitude(), loc.getDepth());
			outlineXY.set(outlineXY.get(0));
			funcs.add(outlineXY);
			chars.add(outlineChar);
			
			DefaultXY_DataSet traceXY = new DefaultXY_DataSet();
			for (Location loc : sect.getFaultTrace())
				traceXY.set(loc.getLongitude(), loc.getDepth());
			funcs.add(traceXY);
			chars.add(traceChar);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Longitude", "Depth (km)");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		double minX = Double.POSITIVE_INFINITY;
		double maxX = Double.NEGATIVE_INFINITY;
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		for (XY_DataSet xy : funcs) {
			for (Point2D pt : xy) {
				minX = Math.min(minX, pt.getX());
				maxX = Math.max(maxX, pt.getX());
				minY = Math.min(minY, pt.getY());
				maxY = Math.max(maxY, pt.getY());
			}
		}
		Range xRange = new Range(Math.floor((minX-0.25)*2d)/2d, Math.ceil((maxX+0.25)*2d)/2d);
		Range yRange = new Range(Math.floor((minY-5d)*2d)/2d, Math.ceil((maxY+5d)*2d)/2d);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().getChart().getXYPlot().getRangeAxis().setInverted(true);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 700, true, false, false);
	}

}
