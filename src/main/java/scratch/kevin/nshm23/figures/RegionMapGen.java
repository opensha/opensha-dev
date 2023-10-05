package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.*;

public class RegionMapGen {

	public static void main(String[] args) throws IOException {
		File dir = Regional_MFD_Plots.NSHM23_PLOTS_DIR;
		
		List<Region> regions = new ArrayList<>();
		
		for (SeismicityRegions seisReg : SeismicityRegions.values())
			regions.add(seisReg.load());
		
		for (AnalysisRegions aReg : AnalysisRegions.values())
			regions.add(aReg.load());
		
		for (LocalRegions lReg : LocalRegions.values())
			regions.add(lReg.load());
		
		double minLat = Double.POSITIVE_INFINITY;
		double maxLat = Double.NEGATIVE_INFINITY;
		double minLon = Double.POSITIVE_INFINITY;
		double maxLon = Double.NEGATIVE_INFINITY;
		
		for (Region region : regions) {
			minLat = Math.min(minLat, region.getMinLat());
			maxLat = Math.max(maxLat, region.getMaxLat());
			minLon = Math.min(minLon, region.getMinLon());
			maxLon = Math.max(maxLon, region.getMaxLon());
		}
		minLat -= 0.1;
		minLon -= 0.1;
		maxLat += 0.1;
		maxLon += 0.1;
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(
				new Region(new Location(minLat, minLon), new Location(maxLat, maxLon)));
		
		mapMaker.plotInsetRegions(regions,
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK), new Color(100, 100, 200), 0.1);
		
		mapMaker.setWritePDFs(true);
		mapMaker.setWriteGeoJSON(true);
		mapMaker.setDefaultPlotWidth(1200);
		mapMaker.plot(dir, "conus_regions", " ");
	}

}
