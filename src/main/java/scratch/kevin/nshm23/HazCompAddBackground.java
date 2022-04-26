package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.HazardMapPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

public class HazCompAddBackground {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/");
//				+ "2022_03_31-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrTap-TotNuclRate-SubB1-JumpProb-2000ip/"
//				+ "hazard_maps_comp_ucerf3_new_anneal/resources/");
		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;
		double spacing = 0.1;
		double period = 0d;
		IncludeBackgroundOption bgOp = IncludeBackgroundOption.INCLUDE;
//		IncludeBackgroundOption bgOp = IncludeBackgroundOption.EXCLUDE;
		GriddedRegion reg = new GriddedRegion(new CaliforniaRegions.RELM_TESTING(), spacing, GriddedRegion.ANCHOR_0_0);
		
		FaultSystemSolution primarySol = FaultSystemSolution.load(new File(invDir,
				"2022_03_24-u3_branches-FM3_1-2000ip/results_FM3_1_branch_averaged.zip"));
//				"2022_03_31-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrTap-TotNuclRate-SubB1-JumpProb-2000ip"
//				+ "/results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		FaultSystemSolution compSol = FaultSystemSolution.load(new File(invDir,
				"2021_11_30-u3_branches-orig_calcs-5h/results_FM3_1_branch_averaged.zip"));
		
		Preconditions.checkState(compSol.hasModule(GridSourceProvider.class));
		primarySol.addModule(compSol.getGridSourceProvider());
		Preconditions.checkState(primarySol.hasModule(GridSourceProvider.class));
		
		HazardMapPlot plot = new HazardMapPlot(gmpeRef, spacing, 0d);
		File outputDir = new File("/tmp/hazard_with_gridded");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		SolHazardMapCalc calc = new SolHazardMapCalc(primarySol, gmpeRef, reg, bgOp, period);
		SolHazardMapCalc compCalc = new SolHazardMapCalc(compSol, gmpeRef, reg, bgOp, period);
		
		calc.calcHazardCurves(32);
		calc.writeCurvesCSV(new File(outputDir, "curves.csv"), period);
		compCalc.calcHazardCurves(32);
		compCalc.writeCurvesCSV(new File(outputDir, "curves_comp.csv"), period);
		
		plot.plot(outputDir, "", "", reg, calc, compCalc);
		
	}

}
