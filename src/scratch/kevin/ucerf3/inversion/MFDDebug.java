package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileUtils;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class MFDDebug {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException {
		int parentID = 206;
		
		File dir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions");
		File fssFile = new File(dir, "2013_01_14-stampede_3p2_production_runs_combined_COMPOUND_SOL.zip");
		
		CompoundFaultSystemSolution fss = CompoundFaultSystemSolution.fromZipFile(fssFile);
		
		Map<FaultModels, List<Integer>> rupsListMap = Maps.newHashMap();
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		
		LogicTreeBranch minBranch = null;
		
		List<BranchVal> valList = Lists.newArrayList();
		
		for (LogicTreeBranch branch : fss.getBranches()) {
			FaultModels fm = branch.getValue(FaultModels.class);
			if (track.getNum() % 50 == 0)
				System.out.println("Workng on solution "+(track.getNum()+1));
			List<Integer> rupList = rupsListMap.get(fm);
			if (rupList == null) {
				FaultSystemSolution sol = fss.getSolution(branch);
				rupList = sol.getRupSet().getRupturesForParentSection(parentID);
				rupsListMap.put(fm, rupList);
			}
			
			double[] rates = fss.getRates(branch);
			
			double sum = 0;
			List<Double> ratesList = Lists.newArrayList();
			for (int rupID : rupList) {
				sum += rates[rupID];
				ratesList.add(rates[rupID]);
			}
			
			if (sum < 1e-20)
//				System.out.println(fm+" rups: "+Joiner.on(",").join(rupList));
				System.out.println(fm+" rates: "+Joiner.on(",").join(ratesList));
			
			if (sum < track.getMin())
				minBranch = branch;
			
			valList.add(new BranchVal(sum, branch));
			
			track.addValue(sum);
		}
		
		Collections.sort(valList);
		
		System.out.println(track);
		System.out.println("Min branch: "+minBranch.buildFileName());
		System.out.println("Min branches:");
		for (int i=0; i<200; i++)
			System.out.println(i+". "+valList.get(i));
		
		listBranchCorrelations(valList, 38);
		
		File binDir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_12_27-ucerf3p2_prod_runs_1/bins/" +
				"2012_12_27-ucerf3p2_prod_runs_1/FM3_1_ABM_EllBsqrtLen_DsrUni_CharConst_M5Rate7.6_MMaxOff8.0_" +
				"NoFix_SpatSeisU2_run2");
		File ratesBin = new File(binDir, "FM3_1_ABM_EllBsqrtLen_DsrUni_CharConst_M5Rate7.6_MMaxOff8.0_NoFix_SpatSeisU2_run2.bin");
		File noMinRatesBin = new File(binDir, "FM3_1_ABM_EllBsqrtLen_DsrUni_CharConst_M5Rate7.6_MMaxOff8.0_NoFix_SpatSeisU2_run2_noMinRates.bin");
		
		double[] theRates = MatrixIO.doubleArrayFromFile(ratesBin);
		double[] theRatesNoMin = MatrixIO.doubleArrayFromFile(noMinRatesBin);
		
		for (int rupID : rupsListMap.get(FaultModels.FM3_1)) {
			double waterlevel = theRates[rupID] - theRatesNoMin[rupID];
			System.out.println("Rup: "+rupID+". Waterlevel: "+waterlevel);
		}
	}
	
	private static class BranchVal implements Comparable<BranchVal> {
		private double val;
		private LogicTreeBranch branch;
		public BranchVal(double val, LogicTreeBranch branch) {
			this.val = val;
			this.branch = branch;
		}
		@Override
		public int compareTo(BranchVal o) {
			return Double.compare(val, o.val);
		}
		
		public String toString() {
			return val+": "+branch.buildFileName();
		}
	}
	
	private static void listBranchCorrelations(List<BranchVal> branchVals, int num) {
		Map<Class<? extends LogicTreeBranchNode<?>>, int[]> countsMap = Maps.newHashMap();
		
		for (Class<? extends LogicTreeBranchNode<?>> clazz : LogicTreeBranch.getLogicTreeNodeClasses()) {
			int[] counts = new int[clazz.getEnumConstants().length];
			countsMap.put(clazz, counts);
		}
		
		for (int i=0; i<num; i++) {
			LogicTreeBranch branch = branchVals.get(i).branch;
			for (LogicTreeBranchNode<?> node : branch) {
				Class<? extends LogicTreeBranchNode> enclClass = LogicTreeBranch.getEnumEnclosingClass(node.getClass());
				int[] counts = countsMap.get(enclClass);
				LogicTreeBranchNode<?>[] consts = enclClass.getEnumConstants();
				int ind = -1;
				for (int j=0; j<consts.length; j++)
					if (consts[j].name().equals(node.name()))
						ind = j;
				counts[ind] = counts[ind] +1;
			}
		}
		
		for (Class<? extends LogicTreeBranchNode<?>> clazz : LogicTreeBranch.getLogicTreeNodeClasses()) {
			LogicTreeBranchNode<?>[] consts = clazz.getEnumConstants();
			int[] counts = countsMap.get(clazz);
			System.out.println("\nBranch Choice: "+ClassUtils.getClassNameWithoutPackage(clazz));
			for (int i=0; i<counts.length; i++)
				if (counts[i] > 0)
					System.out.println("\t"+consts[i].getShortName()+": "+counts[i]);
		}
	}
	
	private static HashSet<Integer> loadMinimizedFM3_1Rups() {
		HashSet<Integer> ids = new HashSet<Integer>();
		try {
			for (String line : FileUtils.loadFile("/tmp/fm3_1_minimized.txt")) {
				line = line.trim();
				if (line.isEmpty())
					continue;
				String[] split = line.split(" ");
				int id = Integer.parseInt(split[split.length-1]);
				ids.add(id);
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return ids;
	}

}
