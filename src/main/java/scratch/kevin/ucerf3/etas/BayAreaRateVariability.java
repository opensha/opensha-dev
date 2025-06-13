package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.LocalRegions;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class BayAreaRateVariability {

	public static void main(String[] args) throws IOException {
		File etasDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2024_05_24-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog/");
		File configFile = new File(etasDir, "config.json");
		File binFile = new File(etasDir, "results_m5_preserve_chain.bin");
		
		ETAS_Config config = ETAS_Config.readJSON(configFile);
		
//		Region reg = new CaliforniaRegions.CYBERSHAKE_BAY_AREA_MAP_REGION();
		Region reg = LocalRegions.CONUS_SF_BAY.load();
		double windowLen = 28d;
		List<Double> counts = new ArrayList<>();
		
		long durationMillis = (long) (windowLen*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		long delta = (long) ProbabilityModelsCalc.MILLISEC_PER_YEAR;
//		long delta = durationMillis;
		long overallStartTime = config.getSimulationStartTimeMillis();
		long overallEndTime = overallStartTime + (long)(ProbabilityModelsCalc.MILLISEC_PER_YEAR*config.getDuration());
		
		int windows = 0;
		for (ETAS_Catalog cat : ETAS_CatalogIO.getBinaryCatalogsIterable(binFile, 5d)) {
			long startTime = config.getSimulationStartTimeMillis();
			while (startTime + durationMillis < overallEndTime) {
				windows++;
				long endTime = startTime + durationMillis;
				
				int count = 0;
				for (ETAS_EqkRupture event : cat) {
					long time = event.getOriginTime();
					if (time < startTime)
						continue;
					if (time > endTime)
						break;
					Preconditions.checkState(event.getMag() >= 5d);
					if (reg.contains(event.getHypocenterLocation()))
						count++;
				}
				
				counts.add((double)count);
				
				startTime += delta;
			}
		}
		System.out.println("Processed "+windows+" "+(int)windowLen+"-year windows");
		
		double[] countsArray = Doubles.toArray(counts);
		double mean = StatUtils.mean(countsArray);
		double median = DataUtils.median(countsArray);
		double min = StatUtils.min(countsArray);
		double max = StatUtils.max(countsArray);
		System.out.println("M>5 SF box nucleation counts:");
		System.out.println("Mean: "+(float)mean);
		System.out.println("Median: "+(int)median);
		System.out.println("Min: "+(int)min);
		System.out.println("Max: "+(int)max);
		
		FeatureCollection features = new FeatureCollection(
				List.of(LocalRegions.CONUS_SF_BAY.load().toFeature(),
						new CaliforniaRegions.CYBERSHAKE_BAY_AREA_MAP_REGION().toFeature()));
		FeatureCollection.write(features, new File("/tmp/sf_boxes.geojson"));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(new Region(new Location(35, -125), new Location(40, -118)));
		List<PlotCurveCharacterstics> outlineChars = List.of(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		List<Color> fillColors = List.of(Colors.tab_orange, Colors.tab_blue);
		mapMaker.plotInsetRegions(List.of(LocalRegions.CONUS_SF_BAY.load(),
						new CaliforniaRegions.CYBERSHAKE_BAY_AREA_MAP_REGION()), outlineChars, fillColors, 0.4);
		mapMaker.setWriteGeoJSON(true);
		mapMaker.plot(new File("/tmp"), "sf_box_map", " ");
	}

}
