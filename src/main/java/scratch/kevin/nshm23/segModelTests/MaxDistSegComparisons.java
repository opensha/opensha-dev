package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateSegmentationConstraint.RateCombiner;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.MaxJumpDistModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class MaxDistSegComparisons {

	public static void main(String[] args) throws IOException {
		File inputFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
//				+ "2022_02_08-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip/"
//				+ "node_branch_averaged/SegModel_ShawR0_3.zip");
//				+ "results_FM3_1_CoulombRupSet_branch_averaged.zip");
				+ "2022_01_28-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip/"
				+ "results_FM3_1_CoulombRupSet_branch_averaged_reweight_r0_3.0.zip");
//				+ "node_branch_averaged/MaxDist_MaxDist3km.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(inputFile);
		ClusterRuptures cRups = ClusterRuptures.singleStranged(sol.getRupSet());
		PlausibilityConfiguration config = sol.getRupSet().getModule(PlausibilityConfiguration.class);
		ClusterConnectionStrategy connStrat = config.getConnectionStrategy();
		SegmentationCalculator calc = new SegmentationCalculator(sol, cRups.getAll(),
				connStrat, config.getDistAzCalc(), new double[] { 0d });
		calc = calc.combineMultiJumps(true);
		
		File parent = inputFile.getParentFile();
		if (parent.getName().startsWith("node_branch"))
			parent = parent.getParentFile();
		File outputDir = new File(parent, "seg_comparisons");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		calc.plotDistDependComparison(outputDir, "shaw_test", true, RateCombiner.MIN);
		double a = inputFile.getParentFile().getName().contains("-max_dist") ? MaxJumpDistModels.invertForWeights() : 1d;
		for (MaxJumpDistModels maxDist : MaxJumpDistModels.values())
			System.out.println(maxDist.getName()+": "+maxDist.getNodeWeight(null));
		System.out.println("a="+a);
		calc.plotMaxDistModelsComparison(outputDir, "max_dist_compare", false,
				MaxJumpDistModels.WEIGHT_TARGET_R0, a, RateCombiner.MIN);
		
		// now map view
		RupSetMapMaker mapMaker = new RupSetMapMaker(sol.getRupSet(), new CaliforniaRegions.RELM_TESTING());
		CPT minDistCPT = new CPT();
		double prevDist = 0d;
		MaxJumpDistModels[] values = MaxJumpDistModels.values();
		CPT colorIndexCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, values.length-1d);
		for (int i = 0; i < values.length; i++) {
			MaxJumpDistModels maxDist = values[i];
			Color color = colorIndexCPT.getColor((float)i);;
			double dist = maxDist.getMaxDist();
			minDistCPT.add(new CPTVal((float)prevDist, color, (float)dist, color));
			prevDist = dist;
		}
		minDistCPT.setAboveMaxColor(Color.GRAY);
		double[] sectScalars = new double[sol.getRupSet().getNumSections()];
		for (int s=0; s<sectScalars.length; s++)
			sectScalars[s] = 9999999d;
		
		for (ClusterRupture rup : cRups) {
			for (Jump jump : rup.getJumpsIterable()) {
				for (FaultSection sect : jump.fromCluster.subSects)
					sectScalars[sect.getSectionId()] = Math.min(sectScalars[sect.getSectionId()], jump.distance);
				for (FaultSection sect : jump.toCluster.subSects)
					sectScalars[sect.getSectionId()] = Math.min(sectScalars[sect.getSectionId()], jump.distance);
			}
		}
		
		mapMaker.plotSectScalars(sectScalars, minDistCPT, "Closest Neighbor Jump Distance (km)");
		mapMaker.plot(outputDir, "sect_closest_dists", " ");
		
		mapMaker.clearSectScalars();
		Map<Jump, Double> jumpDists = new HashMap<>();
		Map<Jump, Double> jumpNormDists = new HashMap<>();
		for (Jump jump : calc.getNonZeroJumps()) {
			jumpDists.put(jump, jump.distance);
			for (int i=0; i<values.length-1; i++) {
				double dist1 = values[i].getMaxDist();
				double dist2 = values[i+1].getMaxDist();
				if (jump.distance >= dist1 && jump.distance < dist2) {
					jumpNormDists.put(jump, (jump.distance - dist1) / (dist2 - dist1));
					break;
				}
			}
		}
		
		mapMaker.plotJumpScalars(jumpDists, minDistCPT, "Jump Distance (km)");
		mapMaker.plot(outputDir, "sect_jump_dists", " ");
		
		mapMaker.clearJumpScalars();
		CPT normCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(0d, 1d);
		mapMaker.plotJumpScalars(jumpNormDists, normCPT, "Jump Distance In Bin");
		mapMaker.plot(outputDir, "sect_jump_norm_dists", " ");
	}

}
