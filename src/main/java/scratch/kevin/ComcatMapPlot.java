package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.comcat.ComcatAccessor;
import org.opensha.commons.data.comcat.ComcatRegion;
import org.opensha.commons.data.comcat.ComcatRegionAdapter;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import net.mahdilamb.colormap.Colors;

public class ComcatMapPlot {

	public static void main(String[] args) throws IOException {
		double minMag = 6d;
		double maxMag = 10;
		int startYear = 1900;
//		int startYear = 1981;
		int endYear = 2025;
		Date startDate = new GregorianCalendar(startYear, 0, 1).getTime();
		Date endDate = new GregorianCalendar(endYear, 0, 1).getTime();
		
		Region region = AnalysisRegions.CONUS_U3_RELM.load();
		String title = "ComCat Events, CA, "+startYear+"-"+endYear;
		String prefix = "ca_comcat_map_plot_"+startYear+"_"+endYear;
		
//		Region region = new CaliforniaRegions.RELM_SOCAL();
//		String title = "ComCat Events, Southern CA, "+startYear+"-"+endYear;
//		String prefix = "socal_comcat_map_plot_"+startYear+"_"+endYear;
		
//		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
//		String title = "ComCat Events, Western US, "+startYear+"-"+endYear;
//		String prefix = "wus_comcat_map_plot_"+startYear+"_"+endYear;
		
		File outputDir = new File("/tmp");
		
		double minDepth = -10;
		double maxDepth = 200;
		
		ComcatRegion cReg = new ComcatRegionAdapter(region);
		
		ComcatAccessor access = new ComcatAccessor();
		
		ObsEqkRupList events = access.fetchEventList(null, startDate.getTime(), endDate.getTime(), minDepth, maxDepth,
				cReg, false, false, minMag);
		
		System.out.println("Loaded "+events.size()+" events");

		GeographicMapMaker mapMaker = new GeographicMapMaker(region);
		
		List<Location> locs = new ArrayList<>(events.size());
		List<PlotCurveCharacterstics> chars = new ArrayList<>(events.size());
		CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(startYear, endYear);
		cpt.setPreferredTickInterval(10d);
		
		for (ObsEqkRupture rup : events) {
			locs.add(rup.getHypocenterLocation());
			int year = rup.getOriginTimeCal().get(Calendar.YEAR);
			
			Color color = cpt.getColor(year);
			
			double mag = rup.getMag();
			double size = 3d + 4*(mag-minMag);
			PlotSymbol symbol;
			if ((float)mag <= 6.5f)
				symbol = PlotSymbol.FILLED_CIRCLE;
			else if ((float)mag < 7f)
				symbol = PlotSymbol.FILLED_SQUARE;
			else if ((float)mag < 7.5f)
				symbol = PlotSymbol.FILLED_TRIANGLE;
			else
				symbol = PlotSymbol.FILLED_DIAMOND;
			chars.add(new PlotCurveCharacterstics(symbol, (float)size, color));
		}
		
		mapMaker.plotScatters(locs, chars, cpt, "Year");
		
		mapMaker.plot(outputDir, prefix, title);
	}

}