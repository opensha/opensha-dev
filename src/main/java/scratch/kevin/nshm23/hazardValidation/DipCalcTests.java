package scratch.kevin.nshm23.hazardValidation;

import static org.opensha.commons.geo.GeoTools.TO_RAD;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.util.TectonicRegionType;

import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.nshmp.shaded.model.NshmErf;
import scratch.kevin.nshm23.hazardValidation.WrapperRupSetRupMapper.WrapperMatch;

public class DipCalcTests {

	public static void main(String[] args) throws IOException {
		File inputSolFile = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
//				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip");
				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_matchPrev.zip");
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
//		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		
		NshmpHazardModel model = NshmpHazardModel.load(modelDir.toPath());
		
		NshmErf wrapper = new NshmErf(model, Set.of(TectonicRegionType.ACTIVE_SHALLOW), IncludeBackgroundOption.EXCLUDE);
		wrapper.getTimeSpan().setDuration(1d);
		wrapper.updateForecast();
		
		WrapperMatch[] mappings = WrapperRupSetRupMapper.map(sol, wrapper);
		
		int rupIndex = 229940;
		ProbEqkRupture wrapperRup = mappings[rupIndex].wrapperRup();
		RuptureSurface wrapperSurf = wrapperRup.getRuptureSurface();
		RuptureSurface fssSurf = sol.getRupSet().getSurfaceForRupture(rupIndex, 1d);
		
		System.out.println("Wrapper rup:");
		System.out.println("\tDepth range: ["+(float)wrapperSurf.getAveRupTopDepth()+", "+(float)wrapperSurf.getAveRupBottomDepth()+"]");
		System.out.println("\tDDW: "+(float)wrapperSurf.getAveWidth());
		System.out.println("\tDip: "+(float)wrapperSurf.getAveDip());
		System.out.println("\tZHyp-from-dip: "+(float)nhazZHypeCalc(wrapperSurf));
		System.out.println("\tZHyp-from-middle: "+(float)middleZHypeCalc(wrapperSurf));
		System.out.println("FSS rup:");
		System.out.println("\tDepth range: ["+(float)fssSurf.getAveRupTopDepth()+", "+(float)fssSurf.getAveRupBottomDepth()+"]");
		System.out.println("\tDDW: "+(float)fssSurf.getAveWidth());
		System.out.println("\tDip: "+(float)fssSurf.getAveDip());
		System.out.println("\tZHyp: "+(float)nhazZHypeCalc(fssSurf));
		System.out.println("\tZHyp-from-middle: "+(float)middleZHypeCalc(fssSurf));
		
		// calculate raw area-weighted dip
		CompoundSurface cSurf = (CompoundSurface)fssSurf;
		double sumArea = 0d;
		double sumDipArea = 0d;
		for (RuptureSurface surf : cSurf.getSurfaceList()) {
			double area = surf.getArea();
			double dip = surf.getAveDip();
			sumArea += area;
			sumDipArea += dip*area;
		}
		double rawAvgDip = sumDipArea/sumArea;
		System.out.println("Area-averaged-dip:\t"+(float)rawAvgDip);
	}
	
	private static double nhazZHypeCalc(RuptureSurface surf) {
		return surf.getAveRupTopDepth() +
				Math.sin(surf.getAveDip() * TO_RAD) * surf.getAveWidth()/2.0;
	}
	
	private static double middleZHypeCalc(RuptureSurface surf) {
		return 0.5*(surf.getAveRupTopDepth() + surf.getAveRupBottomDepth());
	}

}
