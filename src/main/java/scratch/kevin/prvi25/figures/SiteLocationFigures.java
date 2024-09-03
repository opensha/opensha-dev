package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class SiteLocationFigures {
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/prvi_sites");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Region reg = PRVI25_RegionLoader.loadPRVI_MapExtents();
		
		CSVFile<String> csv = CSVFile.readStream(PRVI25_CrustalFaultModels.class.getResourceAsStream(
				"/data/erf/prvi25/sites/prvi_sites.csv"), true);
		GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
		mapMaker.setWritePDFs(true);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_INV_TRIANGLE, 10f, PlotSymbol.INV_TRIANGLE, Color.BLACK);
		for (int row=1; row<csv.getNumRows(); row++) {
			String name = csv.get(row, 0);
			Location loc = new Location(csv.getDouble(row, 1), csv.getDouble(row, 2));
			
			mapMaker.plotScatters(List.of(loc), Colors.tab_green.darker());
			String fname = name.replaceAll("\\W+", "_");
			mapMaker.plot(outputDir, fname, " ");
		}
	}

}
