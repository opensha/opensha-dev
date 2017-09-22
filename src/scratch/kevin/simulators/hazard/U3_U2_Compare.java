package scratch.kevin.simulators.hazard;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;

import com.google.common.base.Preconditions;

public class U3_U2_Compare {

	public static void main(String[] args) throws Exception {
//		File u3CurveFile = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_07_14-ucerf3-gridded-tests/"
		File u3CurveFile = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_08_07-ucerf3-full-ba-gridded-tests/"
				+ "supra_plus_sub_only/curves/imrs1.bin");
		File u2CurveFile = new File("/home/kevin/Simulators/hazard/ucerf-comparisons/ucerf2-faults/curves/imrs1.bin");
		
		File outputDir = new File("/home/kevin/Simulators/hazard/u3_vs_u2");
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		int[] histRPs = { 1000, 2500, 5000, 10000 };
		int[] nehrpRPs = { 2500 };
		int histHighlightIndex = 1;
		
		BinaryHazardCurveReader u3Reader = new BinaryHazardCurveReader(u3CurveFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> u3Curves = u3Reader.getCurveMap();
		
		BinaryHazardCurveReader u2Reader = new BinaryHazardCurveReader(u2CurveFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> u2Curves = u2Reader.getCurveMap();
		
		HazardMapComparePlotter.plotHists(u3Curves, u2Curves, "UCERF2", gridReg, histRPs, histHighlightIndex, outputDir, true);
		HazardMapComparePlotter.plotNEHRP_Hists(u3Curves, u2Curves, "UCERF2", gridReg, nehrpRPs, outputDir);
		HazardMapComparePlotter.plotMeanStdDevTrend(1e-6, u3Curves, u2Curves, "UCERF2", gridReg, outputDir);
		
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
