package scratch.kevin.prvi25;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysHazardCalcSettings;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenRelSupplier;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.gmm.GmmInput;
import gov.usgs.earthquake.nshmp.gmm.GroundMotion;

class GMMFilterTests {
	
	public static void main(String[] args) {
		List<LogicTreeLevel<? extends LogicTreeNode>> allLevels = new ArrayList<>();
		allLevels.addAll(PRVI25_LogicTree.levelsCrustalGMM);
		allLevels.addAll(PRVI25_LogicTree.levelsInterfaceGMM);
		allLevels.addAll(PRVI25_LogicTree.levelsSlabGMM);
		List<LogicTree<?>> logicTrees = List.of(
				LogicTree.buildExhaustive(PRVI25_LogicTree.levelsCrustalGMM, true),
				LogicTree.buildExhaustive(PRVI25_LogicTree.levelsInterfaceGMM, true),
				LogicTree.buildExhaustive(PRVI25_LogicTree.levelsSlabGMM, true),
				LogicTree.buildExhaustive(allLevels, true)
				);
		
		Options ops = new Options();
		FaultSysHazardCalcSettings.addCommonOptions(ops, false);
		String[] hazardArgs = {
				"--gmm-sigma-trunc-one-sided",
				"3.0",
				"--vs30",
				"260",
		};
		CommandLine cmd = FaultSysTools.parseOptions(ops, hazardArgs, GMMFilterTests.class);
		Map<TectonicRegionType, AttenRelSupplier> upstreamGMMs = FaultSysHazardCalcSettings.getGMMs(cmd);
		
		for (LogicTree<?> tree : logicTrees) {
			System.out.println("Testing LogicTree with "+tree.size()+" levels");
			for (LogicTreeBranch<?> branch : tree) {
				System.out.println("\tBranch: "+branch);
				Map<TectonicRegionType, ? extends Supplier<ScalarIMR>> gmms = FaultSysHazardCalcSettings.getGMM_Suppliers(branch, upstreamGMMs, true);
				for (TectonicRegionType trt : gmms.keySet()) {
					NSHMP_GMM_Wrapper gmm = (NSHMP_GMM_Wrapper)gmms.get(trt).get();
					gmm.setIntensityMeasure(PGA_Param.NAME);
					System.out.println("\t\t"+trt.name()+" GMM: "+gmm.getName());
					System.out.println("\t\t\tVs30 is "+gmm.getSiteParams().getValue(Vs30_Param.NAME));
					System.out.println("\t\t\tTrunc is "+gmm.getParameter(SigmaTruncTypeParam.NAME).getValue());
					if (gmm.getOtherParams().containsParameter(SigmaTruncLevelParam.NAME))
						System.out.println("\t\t\tTrunc level is "+gmm.getParameter(SigmaTruncLevelParam.NAME).getValue());
					gmm.setCurrentGmmInput(GmmInput.builder().withDefaults().build());
					gov.usgs.earthquake.nshmp.tree.LogicTree<GroundMotion> gmmTree = gmm.getGroundMotionTree();
					System.out.println("\t\t\tTree has "+gmmTree.size()+" values: "+gmmTree);
				}
			}
		}
	}

}
