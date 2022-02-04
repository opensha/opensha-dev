package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

public class MisfitComparisonOverTime {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
//		File mainDir = new File(invDir, "2022_01_27-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-SubB1-152_samples-10000ip");
//		File mainDir = new File(invDir, "2022_01_27-nshm23_u3_hybrid_branches-FM3_1-U3RedRupSet-SubB1-152_samples-10000ip");
//		File mainDir = new File(invDir, "2022_01_28-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-SubB1-5000ip");
		File mainDir = new File(invDir, "2022_01_28-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
		int maxCount = 2000;
		int deltaCount = 50;
		int scatterDelta = 500;
		
		List<InversionMisfitProgress> progresses = new ArrayList<>();
		
		Preconditions.checkState(mainDir.exists(), "%s doesn't exist", mainDir.getAbsolutePath());
		
		File resultsFile = new File(mainDir, "results.zip");
		Preconditions.checkState(resultsFile.exists(), "%s doesn't exist", resultsFile.getAbsolutePath());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
		
		ZipFile	zip = new ZipFile(resultsFile);
		
		int numRups = -1;
		
		for (LogicTreeBranch<?> branch : tree) {
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(InversionMisfitProgress.MISFIT_PROGRESS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += InversionMisfitProgress.MISFIT_PROGRESS_FILE_NAME;
			System.out.println("Loading "+entryName);
			ZipEntry entry = zip.getEntry(entryName);
			Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			InversionMisfitProgress stats = new InversionMisfitProgress(csv);
			
			progresses.add(stats);
			
			if (numRups < 0)
				numRups = slt.forBranch(branch).getRupSet().getNumRuptures();
		}
		
		zip.close();

		ArbitrarilyDiscretizedFunc averages = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc targets = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc maxs = new ArbitrarilyDiscretizedFunc();
		
		List<Integer> counts = new ArrayList<>();
		for (int count=deltaCount; count<=maxCount; count+=deltaCount)
			counts.add(count);
		if (counts.get(counts.size()-1) != maxCount)
			counts.add(maxCount);
		
		List<DefaultXY_DataSet> scatters = new ArrayList<>();
		List<Integer> scatterCounts = new ArrayList<>();
		
		Quantity targetQuantity = progresses.get(0).getTargetQuantity();
		
		double maxAvg = 0d;
		for (int c=0; c<counts.size(); c++) {
			int count = counts.get(c);

			MinMaxAveTracker maxToTargetRatiosTrack = new MinMaxAveTracker();
			MinMaxAveTracker avgTrack = new MinMaxAveTracker();
			MinMaxAveTracker targetTrack = new MinMaxAveTracker();
			MinMaxAveTracker maxTrack = new MinMaxAveTracker();
			
			DefaultXY_DataSet scatter = null;
			if (count % scatterDelta == 0) {
				scatter = new DefaultXY_DataSet();
				scatters.add(scatter);
				scatterCounts.add(count);
			}
			
			int matchingIndex = -1;
			long closestDiff = Long.MAX_VALUE;
			long closest = -1l;
			List<Long> iters = progresses.get(0).getIterations();
			long iterations = (long)count*(long)numRups;
			for (int i=0; i<iters.size(); i++) {
				long myIters = iters.get(i);
				if (myIters == iterations) {
					matchingIndex = i;
					break;
				} else {
					long diff = myIters - iterations;
					if (diff < 0l)
						diff = -diff;
					if (diff < closestDiff) {
						closestDiff = diff;
						closest = myIters;
					}
				}
			}
			Preconditions.checkState(matchingIndex >= 0, "No match found for %s iterations. Closest was %s (%s away)",
					iterations, closest, closestDiff);
			System.out.println(count+" iteartions/rup = "+iterations+" iters, index "+matchingIndex);
			
			int refIndex = iters.size()-1;
			
			for (InversionMisfitProgress progress : progresses) {
				InversionMisfitStats myStats = progress.getStats().get(matchingIndex);
				InversionMisfitStats refStats = progress.getStats().get(refIndex);
				
				double myTarget = progress.getTargetVals().get(matchingIndex);
				targetTrack.addValue(myTarget);
				
				double myAvg = avgUncertWeightVal(myStats, targetQuantity);
				double refAvg = avgUncertWeightVal(refStats, targetQuantity);
				
				avgTrack.addValue(myAvg);
				
				double myMax = maxUncertWeightVal(myStats, targetQuantity);
				maxTrack.addValue(myMax);
				maxToTargetRatiosTrack.addValue(myMax/myTarget);
				
				if (scatter != null)
					scatter.set(refAvg, myAvg);
				
				maxAvg = Math.max(maxAvg, myAvg);
			}
			
			System.out.println(count+" stats:");
			System.out.println("\tAverages: "+avgTrack);
			System.out.println("\tTargets: "+targetTrack);
			System.out.println("\tMaxs: "+maxTrack);
			System.out.println("\tMax-To-Target ratios: "+maxToTargetRatiosTrack);

			averages.set((double)count, avgTrack.getAverage());
			targets.set((double)count, targetTrack.getAverage());
			maxs.set((double)count, maxTrack.getAverage());
		}
		
		for (int c=0; c<scatters.size(); c++) {
			int count = scatterCounts.get(c);
			if (count == maxCount)
				continue;
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(scatters.get(c));
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
			
			Range range = new Range(0d, 0.1*Math.ceil(10*maxAvg));
			
			DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
			oneToOne.set(0d, 0d);
			oneToOne.set(range.getUpperBound(), range.getUpperBound());
			
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Convergence Test",
					maxCount+" Iterations Per Rupture Misfits", count+" Iterations Per Rupture Misfits");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, false, range, range);
			
			PlotUtils.writePlots(mainDir, "convergence_scatter_vs_"+count, gp, 800, -1, true, false, false);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		maxs.setName("Max");
		funcs.add(maxs);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		
		averages.setName("Average");
		funcs.add(averages);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		targets.setName("Target");
		funcs.add(targets);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Convergence Test",
				"Iterations Per Rupture", "Average Constraint Misfit");
		spec.setLegendVisible(true);
		
		MinMaxAveTracker yTrack = new MinMaxAveTracker();
		for (DiscretizedFunc func : funcs) {
			yTrack.addValue(func.getMinY());
			yTrack.addValue(func.getMaxY());
		}
		
		for (boolean logY : new boolean[] {false, true}) {
			Range yRange = new Range(0.1*Math.floor(10*yTrack.getMin()), 0.1*Math.ceil(10*yTrack.getMax()));
			if (logY)
				yRange = new Range(Math.max(yRange.getLowerBound(), Math.pow(10, Math.floor(Math.log10(yTrack.getMin())))),
						Math.min(yRange.getUpperBound(), Math.pow(10, Math.ceil(Math.log10(yTrack.getMax())))));
				
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, logY, null, yRange);
			
			String prefix = "constraint_convergence";
			if (logY)
				prefix += "_log";
			PlotUtils.writePlots(mainDir, prefix, gp, 800, 650, true, false, false);
		}
	}
	
	private static double avgUncertWeightVal(InversionMisfitStats stats, Quantity quantity) {
		int num = 0;
		double sum = 0d;
		for (MisfitStats misfits : stats.getStats()) {
			if (misfits.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
				num++;
				sum += misfits.get(quantity);
			}
		}
		return sum / (double)num;
	}
	
	private static double maxUncertWeightVal(InversionMisfitStats stats, Quantity quantity) {
		double max = 0d;
		for (MisfitStats misfits : stats.getStats())
			if (misfits.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
				max = Math.max(max, misfits.get(quantity));
		return max;
	}

}
