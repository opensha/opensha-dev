package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;

public class BranchAverageAverageBuild {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File("/tmp/branch_avg_avg");
		InversionFaultSystemRupSet rupSet = null;
		
		List<double[]> ratesList = Lists.newArrayList();
		
		for (File file : dir.listFiles()) {
			if (!file.isFile())
				continue;
			if (file.getName().endsWith("mean.zip"))
				continue;
			InversionFaultSystemSolution invSol = FaultSystemIO.loadInvSol(file);
			if (rupSet == null)
				rupSet = invSol.getRupSet();
			ratesList.add(invSol.getRateForAllRups());
		}
		AverageFaultSystemSolution avgSol = new AverageFaultSystemSolution(rupSet, ratesList);
		FaultSystemIO.writeSol(avgSol, new File(dir, "mean.zip"));
		
		// now make average no means
		File binsDir = new File("/home/kevin/OpenSHA/UCERF3/inversions/" +
				"2013_05_10-ucerf3p3-production-10runs/bins");
		
		APrioriBranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		
		Map<String, double[]> map = Maps.newHashMap();
		double totWt = loadNoMins(binsDir, weightProv, map);
		totWt /= (double)map.size();
		System.out.println("Tot wt: "+totWt);
		System.out.println("# runs: "+map.size());
		
		// now build averages
		ratesList = Lists.newArrayList();
		for (double[] runningRates : map.values()) {
			// now scale to total rate
			for (int r=0; r<runningRates.length; r++)
				runningRates[r] = runningRates[r]/totWt;
			ratesList.add(runningRates);
		}
		avgSol = new AverageFaultSystemSolution(rupSet, ratesList);
		FaultSystemIO.writeSol(avgSol, new File(dir, "mean_noMins.zip"));
	}
	
	private static double loadNoMins(File dir, APrioriBranchWeightProvider weightProv,
			Map<String, double[]> map) throws IOException {
		double totWt = 0;
		for (File file : dir.listFiles()) {
			if (file.isDirectory()) {
				totWt += loadNoMins(file, weightProv, map);
			} else {
				String name = file.getName();
				if (!name.endsWith("_noMinRates.bin"))
					continue;
				String runStr = name.substring(name.indexOf("_run")+1);
				runStr = runStr.substring(0, runStr.indexOf("_"));
				name = name.substring(0, name.indexOf("_run"));
				LogicTreeBranch branch = LogicTreeBranch.fromFileName(name);
				Preconditions.checkState(branch.isFullySpecified());
				double wt = weightProv.getWeight(branch);
				totWt += wt;
				double[] rates = MatrixIO.doubleArrayFromFile(file);
				double[] runningRates = map.get(runStr);
				if (runningRates == null) {
					runningRates = new double[rates.length];
					map.put(runStr, runningRates);
				}
				for (int r=0; r<rates.length; r++)
					runningRates[r] = runningRates[r] + wt*rates[r];
			}
		}
		return totWt;
	}

}
