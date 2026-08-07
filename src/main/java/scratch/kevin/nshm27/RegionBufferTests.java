package scratch.kevin.nshm27;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_MapRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.nshm27.figures.NSHM27_PaperPaths;

public class RegionBufferTests {

	public static void main(String[] args) throws IOException {
		Region containsRegion = NSHM27_RegionLoader.NSHM27_MapRegions.AMSAM.load();
//		Region region = containsRegion;
		Region region = new Region(new Location(-17, 185), new Location(-11, 193));
		
		XY_DataSet[] outlines = PoliticalBoundariesData.loadDefaultOutlines(region);
		
		double minDist = Double.POSITIVE_INFINITY;
		
		LocationList border = new LocationList(region.getBorder());
		// close it
		border.add(border.first());
		
		for (XY_DataSet xy : outlines) {
			for (Point2D pt : xy) {
				Location loc = new Location(pt.getY(), pt.getX());
				if (containsRegion.contains(loc))
					minDist = Math.min(minDist, border.minDistToLine(loc));
			}
		}
		System.out.println("Nearest coastline is "+(float)minDist+" km from border");
		
		List<Region> regions = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		List<Color> colors = new ArrayList<>();
		
		Color color = Colors.tab_blue;
		
		regions.add(NSHM27_SeismicityRegions.AMSAM.load());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		colors.add(color);
		
		regions.add(region);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLACK));
		colors.add(color);
		
		regions.add(NSHM27_MapRegions.AMSAM.load());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLACK));
		colors.add(color);

		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (Region reg : regions) {
			latTrack.addValue(reg.getMinLat());
			latTrack.addValue(reg.getMaxLat());
			lonTrack.addValue(reg.getMinLon());
			lonTrack.addValue(reg.getMaxLon());
		}
		List<? extends FaultSection> sects = NSHM27_PaperPaths.getInterfaceSolution(NSHM27_SeismicityRegions.AMSAM)
				.getRupSet().getFaultSectionDataList();
		for (FaultSection sect : sects) {
			for (Location loc : sect.getFaultTrace()) {
				latTrack.addValue(loc.lat);
				lonTrack.addValue(loc.lon);
			}
		}
		
		double minLat = Math.floor(latTrack.getMin()-0.5);
		double maxLat = Math.ceil(latTrack.getMax()+0.5);
		double minLon = Math.floor(lonTrack.getMin()-0.5);
		double maxLon = Math.ceil(lonTrack.getMax()+0.5);
		GeographicMapMaker mapMaker = new GeographicMapMaker(
				new Region(new Location(minLat, minLon), new Location(maxLat, maxLon)));
		
		mapMaker.setFaultSections(sects);
		
		mapMaker.plotInsetRegions(regions, chars, colors, 0.3);
		
		mapMaker.plot(new File("/tmp"), "amsam_seismicity_regions", " ");
	}

}
