package scratch.kevin.nshm26;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.NEHRP_TestCity;

public class GeoJSON3DExample {

	public static void main(String[] args) throws IOException {
		List<? extends FaultSection> allSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.WUS_FM_v3);
		Region region = new Region(new Location(33, -118), new Location(35, -120));
		List<FaultSection> plotSects = new ArrayList<>();
		for (FaultSection sect : allSects)
			if (FaultSectionUtils.sectInRegion(region, sect, false))
				plotSects.add(sect);
		System.out.println("Keeping "+plotSects.size()+" sections");
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(region, plotSects);
		
		mapMaker.setFillSurfaces(true);
		mapMaker.setFillVerticalSurfaces(true);
		
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 10d);
		mapMaker.plotSectScalars(s->s.getOrigAveSlipRate(), slipCPT, "Slip Rate (mm/yr)");
		
		CPT colorCPT = GMT_CPT_Files.CATEGORICAL_TAB10.instance().rescale(0d, 1d);
		
		Random r = new Random(123456l);
		
		// add some scatter data
		LocationList scatters = new LocationList();
		List<PlotCurveCharacterstics> scatterChars = new ArrayList<>();
		for (NEHRP_TestCity city : NEHRP_TestCity.values()) {
			Location loc = city.location();
			if (region.contains(loc)) {
				scatters.add(loc);
				scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, (float)(12d*r.nextDouble()), colorCPT.getColor(r.nextDouble())));
			}
		}
		mapMaker.plotScatters(scatters, scatterChars);
		
		mapMaker.setWriteGeoJSON(true);
		mapMaker.plot(new File("/tmp"), "geo3d_test", " ");
	}

}
