package scratch.kevin.simulators.hazard;

import java.awt.Color;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;

public class U3_U2_Compare {

	public static void main(String[] args) throws Exception {
//		File u3CurveFile = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_07_14-ucerf3-gridded-tests/"
//		File u3CurveFile = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_08_07-ucerf3-full-ba-gridded-tests/"
//				+ "supra_plus_sub_only/curves/imrs1.bin");
		File u3CurveFile = new File("/home/kevin/Simulators/hazard/ucerf-comparisons/ucerf3-m6.5-8xPoints/curves_pga/imrs1.bin");
		File u2CurveFile = new File("/home/kevin/Simulators/hazard/ucerf-comparisons/ucerf2-faults/curves_pga/imrs1.bin");
		
		String imtLabel = "PGA (g)";
		
		File outputDir = new File("/home/kevin/Simulators/hazard/u3_vs_u2");
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		int[] mapRPs = { 1000, 2500, 10000};
		int[] histRPs = { 1000, 2500, 5000, 10000 };
		int[] nehrpRPs = { 2500 };
		int histHighlightIndex = 1;
		
		BinaryHazardCurveReader u3Reader = new BinaryHazardCurveReader(u3CurveFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> u3Curves = u3Reader.getCurveMap();
		
		BinaryHazardCurveReader u2Reader = new BinaryHazardCurveReader(u2CurveFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> u2Curves = u2Reader.getCurveMap();
		
		HazardMapComparePlotter.plotHists(u3Curves, u2Curves, "UCERF2", gridReg, histRPs, histHighlightIndex, outputDir, true, imtLabel);
		HazardMapComparePlotter.plotNEHRP_Hists(u3Curves, u2Curves, "UCERF2", gridReg, nehrpRPs, outputDir, imtLabel);
		HazardMapComparePlotter.plotMeanStdDevTrend(1e-6, u3Curves, u2Curves, "UCERF2", gridReg, outputDir, imtLabel);
		
		System.out.println("Plotting maps");
		String catalogName = "UCERF2";
		String catalogFileName = "ucerf2";
		CPT hazardCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		hazardCPT = hazardCPT.rescale(0d, 1.2d);
		hazardCPT.setNanColor(Color.WHITE);
		for (int mapRP : mapRPs) {
			boolean isProbAtIML = false;
			double level = 1d/(double)mapRP;
			
			String durationLabel = mapRP+"yr";
			
			GriddedGeoDataSet rsqsimData = HazardMapComparePlotter.loadFromBinary(gridReg, u2Curves, isProbAtIML, level);
			GriddedGeoDataSet u3Data = HazardMapComparePlotter.loadFromBinary(gridReg, u3Curves, isProbAtIML, level);
			
			File csvFile = new File(outputDir, "map_"+durationLabel+".csv.gz");
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Index", "Longitude", "Latitude", catalogName, "UCERF3");
			for (int i=0; i<gridReg.getNodeCount(); i++) {
				Location loc = rsqsimData.getLocation(i);
				csv.addLine(i+"", (float)loc.getLongitude()+"", (float)loc.getLatitude()+"",
						(float)rsqsimData.get(i)+"", (float)u3Data.get(i)+"");
			}
			System.out.println("Writing CSV to "+csvFile.getAbsolutePath());
			OutputStream gzFileOut = new GZIPOutputStream(new FileOutputStream(csvFile));
			csv.writeToStream(gzFileOut);
			
			HazardMapComparePlotter.plotMaps(outputDir, "map_"+durationLabel+"_"+catalogFileName, rsqsimData, region,
					(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), catalogName+", "+durationLabel+", "+imtLabel, hazardCPT, false);
			HazardMapComparePlotter.plotMaps(outputDir, "map_"+durationLabel+"_u3", u3Data, region,
					(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), "UCERF3, "+durationLabel+", "+imtLabel, hazardCPT, false);
			
			GriddedGeoDataSet ratioData = new GriddedGeoDataSet(gridReg, false);
			for (int i=0; i<gridReg.getNodeCount(); i++)
				ratioData.set(i, rsqsimData.get(i) / u3Data.get(i));
			ratioData.log();
			
			CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance();
			ratioCPT.setNanColor(Color.WHITE);
			ratioCPT = ratioCPT.rescale(-0.2d, 0.2d);
			HazardMapComparePlotter.plotMaps(outputDir, "map_"+durationLabel+"_ratio_log_tight", ratioData, region, (double)ratioCPT.getMinValue(),
					(double)ratioCPT.getMaxValue(), "Ln("+catalogName+" / UCERF3), "+durationLabel+", "+imtLabel, ratioCPT, false);
			ratioCPT = ratioCPT.rescale(-0.5d, 0.5d);
			HazardMapComparePlotter.plotMaps(outputDir, "map_"+durationLabel+"_ratio_log", ratioData, region, (double)ratioCPT.getMinValue(),
					(double)ratioCPT.getMaxValue(), "Ln("+catalogName+" / UCERF3), "+durationLabel+", "+imtLabel, ratioCPT, false);
		}
		HazardMapComparePlotter.waitOnFutures(true);
		
		
//		File curveDir = new File(outputDir, "curves");
//		Preconditions.checkState(curveDir.exists() || curveDir.mkdir());
//		
//		Map<String, Location> curveLocs = new HashMap<>();
//		curveLocs.put("Pasadena", new Location(34.148426, -118.17119));
//		curveLocs.put("San Bernardino", new Location(34.064986, -117.29201));
//		curveLocs.put("USC", new Location(34.0192, -118.286));
//		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
//			curveLocs.put(city.toString(), city.location());
//		}
//		
//		for (String siteName : curveLocs.keySet()) {
//			Location siteLoc = curveLocs.get(siteName);
////			Location gridLoc = rsqsimData.getLocation(rsqsimData.indexOf(siteLoc));
//			Location gridLoc = gridReg.locationForIndex(gridReg.indexForLocation(siteLoc));
//			double dist = LocationUtils.horzDistance(siteLoc, gridLoc);
//			System.out.println(siteName+" is "+(float)dist+" km away from nearest grid point");
//			
//			DiscretizedFunc u3Curve = u3Curves.get(gridLoc);
//			DiscretizedFunc rsqsimCurve = u2Curves.get(gridLoc);
//			DiscretizedFunc u3FullCurve = null;
//			
//			plotCurves(u3Curve, u3FullCurve, rsqsimCurve, "UCERF2", curveDir, siteName, imtLabel, mapRPs);
//		}
	}

}
