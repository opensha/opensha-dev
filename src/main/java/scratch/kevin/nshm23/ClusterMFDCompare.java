package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

public class ClusterMFDCompare {

	public static void main(String[] args) throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_08_08-nshm23_branches-wide_seg_branches-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-SubB1-"
//				+ "ThreshAvgIterRelGR/node_branch_averaged/SegModel_Max2km.zip");
				+ "ThreshAvgIterRelGR/node_branch_averaged/SegModel_HighSeg.zip");
		LogicTreeBranch<LogicTreeNode> defaultBranch = NSHM23_LogicTreeBranch.DEFAULT;
		InversionConfigurationFactory factory = new NSHM23_InvConfigFactory();
		
		int includeSubSect = 4366; // Sangre de Cristo (San Luis)
		boolean invert = true;
		
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		ConnectivityClusters clusters = rupSet.getModule(ConnectivityClusters.class);
		if (clusters == null)
			clusters = ConnectivityClusters.build(rupSet);
		
		ConnectivityCluster cluster = null;
		for (ConnectivityCluster c : clusters) {
			if (c.containsSect(includeSubSect)) {
				cluster = c;
				break;
			}
		}
		
		Preconditions.checkNotNull(cluster);
		
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(rupSet);
		
		SummedMagFreqDist summedTarget = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		SummedMagFreqDist summedSolution = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		
		InversionTargetMFDs targets = rupSet.requireModule(InversionTargetMFDs.class);
		List<? extends IncrementalMagFreqDist> sectTargets = targets.getOnFaultSupraSeisNucleationMFDs();
		
		for (int sectID : cluster.getSectIDs()) {
			summedTarget.addIncrementalMagFreqDist(sectTargets.get(sectID));
			summedSolution.addIncrementalMagFreqDist(sol.calcNucleationMFD_forSect(
					sectID, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size()));
		}
		
		int minNonZeroIndex = refMFD.size();
		int maxNonZeroIndex = 0;
		for (int i=0; i<refMFD.size(); i++) {
			if (summedTarget.getY(i) > 0 || summedSolution.getY(i) > 0) {
				minNonZeroIndex = Integer.min(i, minNonZeroIndex);
				maxNonZeroIndex = Integer.max(i, maxNonZeroIndex);
			}
		}
		
		System.out.println("Mag\tTarget\tSolution");
		for (int i=minNonZeroIndex; i<=maxNonZeroIndex; i++)
			printComp((float)refMFD.getX(i)+"", summedTarget.getY(i), summedSolution.getY(i));
		printComp("TOTAL", summedTarget.calcSumOfY_Vals(), summedSolution.calcSumOfY_Vals());
		
		if (invert) {
			FaultSystemRupSet rupSubSet = rupSet.getForSectionSubSet(cluster.getSectIDs());
			
			LogicTreeBranch<LogicTreeNode> branch = rupSet.requireModule(LogicTreeBranch.class).copy();
			for (int i=0; i<branch.size(); i++) {
				if (branch.getValue(i) == null) {
					LogicTreeNode val;
					if (NSHM23_DeformationModels.class.isAssignableFrom(branch.getLevel(i).getType()))
						val = NSHM23_DeformationModels.AVERAGE;
					else if (NSHM23_ScalingRelationships.class.isAssignableFrom(branch.getLevel(i).getType()))
						val = NSHM23_ScalingRelationships.AVERAGE;
					else
						val = defaultBranch.getValue(i);
					System.err.println("setting branch value to default: "+val.getName());
					branch.setValue(val);
				}
			}
			
			InversionConfiguration config = factory.buildInversionConfig(rupSubSet, branch, FaultSysTools.defaultNumThreads());
			
			FaultSystemSolution subsetSol = Inversions.run(rupSubSet, config);
			
			subsetSol.write(new File("/tmp/subset_"+solFile.getName()));
		}
	}
	
	private static final DecimalFormat eDF = new DecimalFormat("0.000E0");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static void printComp(String prefix, double targetVal, double solVal) {
		String line = prefix+"\t"+eDF.format(targetVal)+"\t"+eDF.format(solVal)+"\t";
		if (solVal > targetVal)
			line += "+";
		line += pDF.format((solVal - targetVal)/targetVal);
		System.out.println(line);
	}

}
