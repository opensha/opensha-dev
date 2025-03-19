package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class SiteLocationFigures {
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/prvi_sites");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Region reg = PRVI25_RegionLoader.loadPRVI_MapExtents();
		
		List<Site> sites = PRVI25_RegionLoader.loadHazardSites();
		GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
		mapMaker.setWritePDFs(true);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_INV_TRIANGLE, 10f, PlotSymbol.INV_TRIANGLE, Color.BLACK);
		for (Site site : sites) {
			mapMaker.plotScatters(List.of(site.getLocation()), Colors.tab_green.darker());
			String fname = FileNameUtils.simplify(site.getName());
			mapMaker.plot(outputDir, fname, " ");
		}
	}

}
