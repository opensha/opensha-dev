package scratch.kevin.prvi25;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;

public class SeisCatalogRegionTests {
	
	public static void main(String[] args) throws IOException {
		Region testReg = PRVI25_SeismicityRegions.CAR_INTERFACE.load();
		
//		String label = "iscgem899944 1943 7.76";
//		Location testLoc = new Location(18.907, -67.118);

//		String label = "iscgem898498 1946 7.7";
//		Location testLoc = new Location(19.124, -69.103);

		String label = "usp00037fr 1987 2.8";
		Location testLoc = new Location(15.11100, -62.06500);
		
		System.out.println("Testing location: "+testLoc);
		if (testReg.contains(testLoc)) {
			System.out.println("\tcontained!");
			label += " (contained)";
		} else {
			String away = new DecimalFormat("0.0").format(testReg.distanceToLocation(testLoc))+" km outside";
			System.out.println("\t"+away);
			label += "("+away+")";
		}
		
		double minLat = testReg.getMinLat();
		double maxLat = testReg.getMaxLat();
		double minLon = testReg.getMinLon();
		double maxLon = testReg.getMaxLon();
		
		minLat = Math.min(minLat, testLoc.getLatitude());
		maxLat = Math.max(maxLat, testLoc.getLatitude());
		minLon = Math.min(minLon, testLoc.getLongitude());
		maxLon = Math.max(maxLon, testLoc.getLongitude());
		
		Region buffered = new Region(new Location(minLat-0.5, minLon-0.5),
				new Location(maxLat+0.5, maxLon+0.5));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(buffered);
		
		mapMaker.plotInsetRegions(List.of(testReg), new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK), null, 0);
		
		mapMaker.plotScatters(List.of(testLoc), List.of(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.BLUE)));
		
		mapMaker.plot(new File("/tmp"), "prvi_reg_test", label);
	}

}
