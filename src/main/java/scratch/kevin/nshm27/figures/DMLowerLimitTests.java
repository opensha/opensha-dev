package scratch.kevin.nshm27.figures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeLevel.SamplingMethod;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.dmSampling.DeformationModelDistSampler.FixedFractileSampler;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM27_InvConfigFactory;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceHingedBValue;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;

public class DMLowerLimitTests {

	public static void main(String[] args) throws IOException {
		ModuleContainer.VERBOSE_DEFAULT = false;
		NSHM27_SeismicityRegions seisReg = NSHM27_SeismicityRegions.GNMI;
		int numSamples = 50000;
		NSHM27_LogicTree.INTERFACE_B_HINGED_WEIGHT = 1d;
		LogicTree<LogicTreeNode> tree = NSHM27_LogicTree.buildLogicTree(seisReg,
				TectonicRegionType.SUBDUCTION_INTERFACE, numSamples, true,
				SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE);
		int fixedIndex = 27551;
		double searchThreshold = 0d;
		
		NSHM27_InvConfigFactory factory = new NSHM27_InvConfigFactory();
		FaultSystemRupSet rs = factory.buildGenericRupSet(tree.getBranch(0), 16);
		SectSlipRates slips1 = rs.requireModule(SectSlipRates.class);
		FaultSystemRupSet genericRS = rs;
		
		double minDM = Double.POSITIVE_INFINITY;
		LogicTreeBranch<LogicTreeNode> minDMBranch = null;
		if (fixedIndex >= 0) {
			LogicTreeBranch<LogicTreeNode> branch = tree.getBranch(fixedIndex);
			minDM = ((FixedFractileSampler)branch.requireValue(NSHM27_InterfaceDeformationModels.class).getValue()).getFixedFractile();
			minDMBranch = branch;
		} else {
			for (LogicTreeBranch<LogicTreeNode> branch : tree) {
				double fract = ((FixedFractileSampler)branch.requireValue(NSHM27_InterfaceDeformationModels.class).getValue()).getFixedFractile();
				if (fract < minDM) {
					minDM = fract;
					minDMBranch = branch;
				}
			}
		}
		System.out.println("Min DM fractile: "+minDM);
		System.out.println("Branch: "+minDMBranch);
		rs = factory.updateRuptureSetForBranch(rs, minDMBranch);
		Preconditions.checkState(rs != genericRS);
		SectSlipRates slips2 = rs.requireModule(SectSlipRates.class);
		Preconditions.checkState(slips1 != slips2);
		for (int s=0; s<rs.getNumSections(); s++) {
			FaultSection sect1 = rs.getFaultSectionData(s);
			FaultSection sect2 = genericRS.getFaultSectionData(s);
			Preconditions.checkState(sect1 != sect2);
			Preconditions.checkState(sect1.getOrigAveSlipRate() != sect2.getOrigAveSlipRate());
			Preconditions.checkState(slips1.getSlipRate(s) == 0d || slips1.getSlipRate(s) != slips2.getSlipRate(s));
		}
		FaultGridAssociations assoc = rs.requireModule(FaultGridAssociations.class);
		MinMaxAveTracker slipTrack = new MinMaxAveTracker();
		MinMaxAveTracker momentTrack = new MinMaxAveTracker();
		MinMaxAveTracker assocMomentTrack = new MinMaxAveTracker();
		for (FaultSection sect : rs.getFaultSectionDataList()) {
			slipTrack.addValue(sect.getReducedAveSlipRate());
			double moment = sect.calcMomentRate(true);
			momentTrack.addValue(moment);
			assocMomentTrack.addValue(moment*assoc.getSectionFractInRegion(sect.getSectionId()));
		}
		System.out.println("Slips:\t"+slipTrack);
		System.out.println("Moments:\t"+momentTrack);
		System.out.println("Assoc moments:\t"+momentTrack);
		double b = NSHM27_InterfaceHingedBValue.calcRawInterfaceHingedBValue(rs, minDMBranch);
		System.out.println(b);
		System.out.println("Running inversion");
		FaultSystemSolution sol = Inversions.run(factory, minDMBranch, FaultSysTools.defaultNumThreads());
		System.out.println("DONE with inversion");
		System.out.println("Total rate: "+sol.getTotalRateForAllFaultSystemRups());
		
		if (searchThreshold > 0d) {
			List<LogicTreeBranch<LogicTreeNode>> lowBranches = new ArrayList<>();
			for (LogicTreeBranch<LogicTreeNode> branch : tree) {
				double fract = ((FixedFractileSampler)branch.requireValue(NSHM27_InterfaceDeformationModels.class).getValue()).getFixedFractile();
				if (fract < searchThreshold)
					lowBranches.add(branch);
			}
			System.out.println("Now trying tree for "+lowBranches.size()+" low branches");
			List<CompletableFuture<Double>> bFutures = new ArrayList<>(lowBranches.size());
			for (LogicTreeBranch<LogicTreeNode> branch : lowBranches) {
				bFutures.add(CompletableFuture.supplyAsync(()->{
					FaultSystemRupSet branchRS;
					try {
						branchRS = factory.updateRuptureSetForBranch(genericRS, branch);
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					return NSHM27_InterfaceHingedBValue.calcRawInterfaceHingedBValue(branchRS, branch);
				}));
			}
			
			for (int i=0; i<bFutures.size(); i++) {
				LogicTreeBranch<LogicTreeNode> branch = lowBranches.get(i);
				try {
					bFutures.get(i).get();
				} catch (Exception e) {
					System.err.println("Failed on "+branch);
					e.printStackTrace();
					System.err.flush();
					System.exit(1);
				}
			}
			System.out.println("All "+lowBranches.size()+" succeeded for searchThreshold="+(float)searchThreshold);
			System.exit(0);
		}
	}

}
