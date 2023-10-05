package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.geo.Region;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAveragingOrder;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

public class MFD_ExceedObsStats {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File modelDir = new File(invsDir,
				"2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		FaultSystemSolution baSol = FaultSystemSolution.load(new File(modelDir,
				"results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		LogicTree<?> tree = LogicTree.read(new File(modelDir, "logic_tree.json"));
		
		double minMag = 6d;
		
		BranchAveragingOrder order = baSol.requireModule(BranchAveragingOrder.class);
		BranchRegionalMFDs regMFDs = baSol.requireModule(BranchRegionalMFDs.class);
		RegionsOfInterest roi = baSol.getRupSet().getModule(RegionsOfInterest.class);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.45);
		
//		boolean contermWUS = true;
//		Region fullRegion = NSHM23_RegionLoader.loadFullConterminousWUS();
//		IncrementalMagFreqDist obsMFD = NSHM23_RegionalSeismicity.getRemapped(fullRegion,
//				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, 7.95);
		
		boolean contermWUS = false;
		IncrementalMagFreqDist obsMFD = NSHM23_RegionalSeismicity.PREFFERRED.build(
				SeismicityRegions.CONUS_WEST, refMFD, 7.95);
		
		ExceedTracker allTrack = new ExceedTracker();
		Map<LogicTreeNode, ExceedTracker> nodeTracks = new HashMap<>();
		
		int startBinIndex = obsMFD.getClosestXIndex(minMag+0.01);
		
		IncrementalMagFreqDist[] branchMFDs;
		if (contermWUS) {
			SummedMagFreqDist[] sumMFDs = new SummedMagFreqDist[tree.size()];
			
			Region[] sumRegions = {
					AnalysisRegions.CONUS_U3_RELM.load(),
					AnalysisRegions.CONUS_IMW.load(),
					AnalysisRegions.CONUS_PNW.load()
			};
			
			for (Region region : sumRegions) {
				int regionIndex = -1;
				List<Region> allRegions = roi.getRegions();
				for (int i=0; i<allRegions.size(); i++) {
					Region testReg = allRegions.get(i);
//					System.out.println("Checking "+testReg.getName()+" against "+region.getName()
//						+": "+testReg.equalsRegion(region));
					if (testReg.equalsRegion(region)) {
						regionIndex = i;
						break;
					}
				}
				System.out.println(region.getName()+" is index "+regionIndex);
				Preconditions.checkState(regionIndex >= 0);
				IncrementalMagFreqDist[] mfds = regMFDs.getRegionalBranchMFDs(MFDType.SUPRA_ONLY, regionIndex);
				
				for (int i=0; i<mfds.length; i++) {
					if (sumMFDs[i] == null) {
						sumMFDs[i] = new SummedMagFreqDist(mfds[i].getMinX(), mfds[i].getMaxX(), mfds[i].size());
					} else {
						Preconditions.checkState(sumMFDs[i].getMinX() == mfds[i].getMinX());
						Preconditions.checkState(sumMFDs[i].size() == mfds[i].size());
					}
					sumMFDs[i].addIncrementalMagFreqDist(mfds[i]);
				}
			}
			
			branchMFDs = sumMFDs;
		} else {
			branchMFDs = regMFDs.getTotalBranchMFDs(MFDType.SUPRA_ONLY);
		}
		
		System.out.println("Obs MFD: "+obsMFD);
		
		for (int t=0; t<tree.size(); t++) {
			LogicTreeBranch<?> branch = tree.getBranch(t);
			int index = order.getBranchAveragingIndex(branch);
			IncrementalMagFreqDist mfd = branchMFDs[index];
			
			int numExceeds = 0;
			for (int i=startBinIndex; i<obsMFD.size(); i++) {
				double obsRate = obsMFD.getY(i);
				int mfdXIndex = mfd.getClosestXIndex(obsMFD.getX(i));
				double supraRate = mfd.getY(mfdXIndex);
				boolean exceeds = obsRate > 0 && supraRate > obsRate;
				if (exceeds)
					numExceeds++;
			}
			
			System.out.println("Branch "+t+"->"+index+" is "+branch);
			System.out.println("\tExceeds: "+numExceeds);
			if (t == 0 || numExceeds > 0) {
				for (int i=startBinIndex; i<obsMFD.size(); i++) {
					double obsRate = obsMFD.getY(i);
					int mfdXIndex = mfd.getClosestXIndex(obsMFD.getX(i));
					double supraRate = mfd.getY(mfdXIndex);
					boolean exceeds = obsRate > 0 && supraRate > obsRate;
					if (exceeds)
						System.out.println("\t"+(float)mfd.getX(mfdXIndex)+"->"+(float)obsMFD.getX(i)
							+":\tobs="+(float)obsRate+"\tsupra="+(float)supraRate+"\texeeds="+exceeds);
				}
			}
			
			allTrack.add(numExceeds);
			for (LogicTreeNode node : branch) {
				ExceedTracker nodeTrack = nodeTracks.get(node);
				if (nodeTrack == null) {
					nodeTrack = new ExceedTracker();
					nodeTracks.put(node, nodeTrack);
				}
				nodeTrack.add(numExceeds);
			}
		}
		
		FileWriter fw = new FileWriter(new File(outputDir, "mfd_exceed_stats.txt"));
		
		logPrint(fw, "All Branchches:\t"+allTrack);
		
		for (int l=0; l<tree.getLevels().size(); l++) {
			LogicTreeLevel<?> level = tree.getLevels().get(l);
			List<LogicTreeNode> nodes = new ArrayList<>();
			for (LogicTreeNode node : level.getNodes())
				if (nodeTracks.containsKey(node))
					nodes.add(node);
			if (nodes.size() > 1) {
				logPrint(fw, level.getName());
				for (LogicTreeNode node : nodes)
					logPrint(fw, "\t"+node.getShortName()+":\t"+nodeTracks.get(node));
			}
		}
		
		fw.close();
	}
	
	private static void logPrint(FileWriter fw, String str) throws IOException {
		fw.write(str+"\n");
		System.out.println(str);
	}
	
	private static DecimalFormat pDF = new DecimalFormat("0.00%");
	private static DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static class ExceedTracker {
		private int numBranches = 0;
		private int numAnyExceed = 0;
		private MinMaxAveTracker binExceedsTrack = new MinMaxAveTracker();
		
		public void add(int numExceeds) {
			numBranches++;
			if(numExceeds > 0) {
				numAnyExceed++;
				binExceedsTrack.addValue((double)numExceeds);
			}
		}
		
		public String toString() {
			return numAnyExceed+"/"+numBranches+" exceedances ("
					+pDF.format((double)numAnyExceed/(double)numBranches)+"), "
					+oDF.format(binExceedsTrack.getAverage())+" bins per exceed";
		}
	}

}
