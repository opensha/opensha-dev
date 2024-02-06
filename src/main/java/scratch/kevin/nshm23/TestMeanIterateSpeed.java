package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.TrueMeanRuptureMappings;

import com.google.common.base.Stopwatch;

public class TestMeanIterateSpeed {

	public static void main(String[] args) throws IOException {
		File dir = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/");
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File(dir, "true_mean_solution.zip"));
		SolutionLogicTree slt = SolutionLogicTree.load(new File(dir, "results.zip"));
		
		boolean loadRates = true;
		boolean loadProps = true;
		boolean loadMappings = true;
		boolean loadGridded = false;
		
		LogicTree<?> tree = slt.getLogicTree();
		TrueMeanRuptureMappings mappings = sol.getRupSet().requireModule(TrueMeanRuptureMappings.class);
		
		Stopwatch rateWatch = Stopwatch.createUnstarted();
		Stopwatch propWatch = Stopwatch.createUnstarted();
		Stopwatch mapWatch = Stopwatch.createUnstarted();
		Stopwatch gridWatch = Stopwatch.createUnstarted();
		
		Stopwatch totalWatch = Stopwatch.createStarted();
		for (int i=0; i<tree.size(); i++) {
			LogicTreeBranch<?> branch = tree.getBranch(i);
			System.out.println(i+". "+branch);
			
			if (loadRates) {
				rateWatch.start();
				slt.loadRatesForBranch(branch);
				rateWatch.stop();
			}
			
			if (loadProps) {
				propWatch.start();
				slt.loadPropsForBranch(branch);
				propWatch.stop();
			}
			
			if (loadMappings) {
				mapWatch.start();
				mappings.getRuptureMappings(branch);
				mapWatch.stop();
			}
			
			if (loadGridded) {
				gridWatch.start();
				slt.loadGridProvForBranch(branch);
				gridWatch.stop();
			}
		}
		totalWatch.stop();
		
		System.out.println("Done in "+timeStr(totalWatch, 1d)+"\t("+timeStr(totalWatch, tree.size())+" each)");
		
		if (loadRates)
			System.out.println("Rates in "+timeStr(rateWatch, 1d)+"\t("+timeStr(rateWatch, tree.size())+" each)");
		if (loadProps)
			System.out.println("Props in "+timeStr(propWatch, 1d)+"\t("+timeStr(propWatch, tree.size())+" each)");
		if (loadMappings)
			System.out.println("Mappings in "+timeStr(mapWatch, 1d)+"\t("+timeStr(mapWatch, tree.size())+" each)");
		if (loadGridded)
			System.out.println("Grid Provs in "+timeStr(gridWatch, 1d)+"\t("+timeStr(gridWatch, tree.size())+" each)");
	}
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	private static String timeStr(Stopwatch watch, double divide) {
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		secs /= divide;
		if (secs < 60d)
			return twoDigits.format(secs)+" s";
		double mins = secs/60d;
		if (mins < 60d)
			return twoDigits.format(mins)+" m";
		double hours = secs/60d;
		return twoDigits.format(hours)+" h";
	}

}
