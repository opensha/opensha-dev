package scratch.kevin.simulators;

import java.io.File;
import java.util.Map;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;

public class U3_U2_Compare {

	public static void main(String[] args) throws Exception {
		File u3CurveFile = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_07_14-ucerf3-gridded-tests/"
				+ "supra_plus_sub_only/curves/imrs1.bin");
		File u2CurveFile = new File("/home/kevin/Simulators/hazard/ucerf-comparisons/ucerf2-faults/curves/imrs1.bin");
		
		File outputDir = new File("/tmp");
		
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
	}

}
