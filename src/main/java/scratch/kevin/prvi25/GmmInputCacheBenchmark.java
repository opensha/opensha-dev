package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.AbstractSitewiseThreadedLogicTreeCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Stopwatch;

class GmmInputCacheBenchmark {

	public static void main(String[] args) throws IOException {
		boolean verbose = false;
//		ModuleContainer.VERBOSE_DEFAULT = verbose;
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		File treeFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_07_01-prvi25_crustal_branches-dmSample5x-gmTreeCalcs/logic_tree.json");
		
//		File resultFsile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_01-prvi25_crustal_branches-dmSample5x/results.zip");
//		IncludeBackgroundOption bgOp = IncludeBackgroundOption.EXCLUDE;
		File resultsFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_01-prvi25_crustal_branches-dmSample5x/results_avg_gridded.zip");
		IncludeBackgroundOption bgOp = IncludeBackgroundOption.INCLUDE;
		
		int reuseMod = 4;
		
		LogicTree<?> tree = LogicTree.read(treeFile);
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile, tree);
		
		Location loc = new Location(18.47,-66.72);
		double[] periods = {0d, 1d};
		ExecutorService exec = Executors.newSingleThreadExecutor();
		AbstractSitewiseThreadedLogicTreeCalc calc = new AbstractSitewiseThreadedLogicTreeCalc(exec, 1, slt,
				AttenRelRef.ASK_2014, periods, bgOp, new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS)) {
			
			@Override
			public Site siteForIndex(int siteIndex, Map<TectonicRegionType, ScalarIMR> gmms) {
				Site site = new Site(loc);
				site.addParameterList(gmms.values().iterator().next().getSiteParams());
				return site;
			}
			
			@Override
			public void debug(String message) {
				if (verbose)
					System.out.println(message);
			}
		};
		
		calc.setDoGmmInputCache(true);
		calc.setCacheGridSources(true);
		
		int printDelta = 10;
		Stopwatch watch = Stopwatch.createStarted();
		int reuseCount = 0;
		Stopwatch reuseWatch = Stopwatch.createUnstarted();
		int newCount = 0;
		Stopwatch newWatch = Stopwatch.createUnstarted();
		for (int i=0; i<tree.size(); i++) {
			if (!verbose)
				System.out.print(".");
			boolean reuse = i % reuseMod != 0;
			if (reuse)
				reuseWatch.start();
			else
				newWatch.start();
			DiscretizedFunc[][] curves = calc.calcForBranch(i);
			if (reuse) {
				reuseCount++;
				reuseWatch.stop();
			} else {
				newCount++;
				newWatch.stop();
			}
			if (verbose) {
				System.out.print("\tCURVE[0][0]:");
				System.out.print("\tsumY="+curves[0][0].calcSumOfY_Vals()+"; vals: ");
				for (int x=0; x<curves[0][0].size(); x++)
					System.out.print(" "+(float)curves[0][0].getY(x));
				System.out.println();
			}
			int numDone = i+1;
			if (numDone % printDelta == 0) {
				System.out.println(numDone+"/"+tree.size()+" done, "+getRate(numDone, watch)
						+" [reuseERF: "+getRate(reuseCount, reuseWatch)+", newERF: "+getRate(newCount, newWatch)+"]");
			}
		}
		watch.stop();
		System.out.println("DONE: "+getRate(tree.size(), watch)
				+" [reuseERF: "+getRate(reuseCount, reuseWatch)+", newERF: "+getRate(reuseCount, reuseWatch)+"]");
	}
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	private static String getRate(int numDone, Stopwatch watch) {
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		double perSec = (double)numDone/secs;
		return twoDigits.format(perSec)+" /sec";
	}

}
