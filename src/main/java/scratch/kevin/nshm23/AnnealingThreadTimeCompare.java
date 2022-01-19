package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

public class AnnealingThreadTimeCompare {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");

//		File mainDir = new File(invDir, "2022_01_18-reproduce-ucerf3-ref_branch-long_reweight_test-NEOK-EllB-DsrTap-SupraB0.0-NuclMFD-reweight_sqrt_conserve_phased-initial_parkfield");
//		File mainDir = new File(invDir, "2022_01_18-reproduce-ucerf3-ref_branch-long_reweight_test-NEOK-EllB-DsrTap-SupraB0.0-TotNuclRate-reweight_sqrt_conserve_phased-initial_parkfield");
//		File mainDir = new File(invDir, "2022_01_18-reproduce-ucerf3-ref_branch-long_reweight_test-ZENGBB-Shaw09Mod-DsrUni-SupraB1.0-NuclMFD-reweight_sqrt_conserve_phased-initial_parkfield");
		File mainDir = new File(invDir, "2022_01_18-reproduce-ucerf3-ref_branch-long_reweight_test-ZENGBB-Shaw09Mod-DsrUni-SupraB1.0-TotNuclRate-reweight_sqrt_conserve_phased-initial_parkfield");
		
		Table<Double, Double, DiscretizedFunc> avgVsRounds = HashBasedTable.create();
		Table<Double, Double, DiscretizedFunc> avgVsTimes = HashBasedTable.create();
		Table<Double, Double, DiscretizedFunc> spreadVsRounds = HashBasedTable.create();
		Table<Double, Double, DiscretizedFunc> spreadVsTimes = HashBasedTable.create();
		
		Quantity quantity = Quantity.MAD;
		int maxPointsWithoutDecimation = 50;
		boolean skipAvg0 = true;
		
		int numRuptures = FaultSystemRupSet.load(new File(mainDir, "rupture_set.zip")).getNumRuptures();
		System.out.println("Num ruptures: "+numRuptures);
		
		double bestFinalAvg = Double.POSITIVE_INFINITY;
		double worstFinalAvg = Double.NEGATIVE_INFINITY;
		double bestFinalSpread = Double.POSITIVE_INFINITY;
		double worstFinalSpread = Double.NEGATIVE_INFINITY;
		
		for (File subDir : mainDir.listFiles()) {
			if (!subDir.isDirectory() || !subDir.getName().startsWith("sub_"))
				continue;
			
			String dirName = subDir.getName();
			String[] split = dirName.split("_");
			double subDelta = Double.parseDouble(split[1]);
			double avgDelta = split.length == 5 ? Double.parseDouble(split[3]) : 0d;
			
			if (skipAvg0 && avgDelta == 0d)
				continue;
			
			System.out.println(dirName+": subDelta="+(float)subDelta+", avgDelta="+(float)avgDelta);
			
			File solFile = new File(subDir, "solution.zip");
			
			if (!solFile.exists()) {
				System.out.println("\tmissing");
				continue;
			}
			
			System.out.println("\tLoading from "+solFile.getAbsolutePath());
			ZipFile zip = new ZipFile(solFile);
			
			ZipEntry entry = zip.getEntry("solution/"+InversionMisfitProgress.MISFIT_PROGRESS_FILE_NAME);
			Preconditions.checkNotNull(entry);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			
			InversionMisfitProgress progress = new InversionMisfitProgress(csv);
			
			List<Long> iters = progress.getIterations();
			List<Long> times = progress.getTimes();
			List<InversionMisfitStats> stats = progress.getStats();
			
			DiscretizedFunc avgVsRound = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc avgVsTime = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc spreadVsRound = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc spreadVsTime = new ArbitrarilyDiscretizedFunc();
			
			for (int i=0; i<iters.size(); i++) {
				// ms -> min
				double time = (double)times.get(i)/(1000d*60d);
				double rounds = (double)iters.get(i)/(double)numRuptures;
				
				double avgMisfit = 0d;
				List<Double> vals = new ArrayList<>();
				for (MisfitStats misfits : stats.get(i).getStats()) {
					double val = misfits.get(quantity);
					vals.add(val);
					avgMisfit += val;
				}
				avgMisfit /= vals.size();
				double spread = 0d;
				for (double val : vals)
					spread += Math.abs(val - avgMisfit)/avgMisfit;
				spread /= vals.size();
				
				avgVsRound.set(rounds, avgMisfit);
				avgVsTime.set(time, avgMisfit);
				spreadVsRound.set(rounds, spread);
				spreadVsTime.set(time, spread);
				
				if (i == iters.size()-1) {
					bestFinalAvg = Math.min(bestFinalAvg, avgMisfit);
					worstFinalAvg = Math.max(worstFinalAvg, avgMisfit);
					bestFinalSpread = Math.min(bestFinalSpread, spread);
					worstFinalSpread = Math.max(worstFinalSpread, spread);
				}
			}
			
			spreadVsTimes.put(subDelta, avgDelta, spreadVsTime);
			avgVsTimes.put(subDelta, avgDelta, avgVsTime);
			spreadVsRounds.put(subDelta, avgDelta, spreadVsRound);
			avgVsRounds.put(subDelta, avgDelta, avgVsRound);
			
			zip.close();
		}
		
		if (avgVsRounds.isEmpty()) {
			System.out.println("No solutions found");
			System.exit(1);
		}
		
		System.out.println("Plotting "+avgVsRounds.size()+" solution misfits");
		
		System.out.println("Final misfit range: ["+(float)bestFinalAvg+", "+(float)worstFinalAvg+"]");
		System.out.println("Final spread range: ["+(float)bestFinalSpread+", "+(float)worstFinalSpread+"]");
		
		CPT rainbow = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
		Map<Double, Color> subDeltaColors = getValColors(spreadVsTimes.rowKeySet(), rainbow);
		Map<Double, Color> avgDeltaColors = getValColors(spreadVsTimes.columnKeySet(), rainbow);
		
		File outputDir = new File(mainDir, "misfits_vs_thread_params");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		for (boolean isSpread : new boolean[] {false, true}) {
			for (boolean isTime : new boolean[] {false, true}) {
				for (boolean isAvg : new boolean[] {false, true}) {
					Table<Double, Double, DiscretizedFunc> table;
					String title = " ", prefix, yAxisLabel, xAxisLabel;
					if (isSpread) {
						prefix = "spread";
						yAxisLabel = "Average Constraint Fractional Spread";
						if (isTime)
							table = spreadVsTimes;
						else
							table = spreadVsRounds;
					} else {
						prefix = "avg";
						yAxisLabel = "Average Constraint "+quantity;
						if (isTime)
							table = avgVsTimes;
						else
							table = avgVsRounds;
					}
					
					if (isTime) {
						prefix += "_vs_time";
						xAxisLabel = "Time (minutes)";
					} else {
						prefix += "_vs_rounds";
						xAxisLabel = "Iteration Rounds (per rupture)";
					}
					
					if (isAvg)
						prefix += "_avg";

					List<DiscretizedFunc> subFuncs = new ArrayList<>();
					List<DiscretizedFunc> avgFuncs = new ArrayList<>();
					List<PlotCurveCharacterstics> subChars = new ArrayList<>();
					List<PlotCurveCharacterstics> avgChars = new ArrayList<>();
					
					double minY = isSpread ? bestFinalSpread : bestFinalAvg;
					double maxY = isSpread ? worstFinalSpread : worstFinalAvg;
					if (isSpread) {
						minY = Math.pow(10, Math.floor(Math.log10(minY)));
						maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
					} else {
						minY *= 0.9;
						maxY *= 1.1;
					}
					double minX = 0d;
					double maxX = 0d;
					
					double maxAny = 0d;
					double minAny = Double.POSITIVE_INFINITY;
					
					if (isAvg) {
						List<Double> subVals = new ArrayList<>(subDeltaColors.keySet());
						Collections.sort(subVals);
						for (double subVal : subVals) {
							DiscretizedFunc func = averageFuncs(table, subVal, true);
							func.setName("Sub="+oDF.format(subVal));
							
							Color color = subDeltaColors.get(subVal);
							color = new Color(color.getRed(), color.getGreen(), color.getBlue());
							
							subFuncs.add(func);
							subChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
							
							maxX = Math.max(maxX, func.getMaxX());
							double halfway = 0.5*func.getMaxX();
							if (halfway >= func.getMinX())
								maxY = Math.max(maxY, func.getInterpolatedY(halfway));
						}
						
						List<Double> avgVals = new ArrayList<>(avgDeltaColors.keySet());
						Collections.sort(avgVals);
						for (double avgVal : avgVals) {
							DiscretizedFunc func = averageFuncs(table, avgVal, false);
							func.setName("Avg="+oDF.format(avgVal));
							
							Color color = avgDeltaColors.get(avgVal);
							color = new Color(color.getRed(), color.getGreen(), color.getBlue());
							
							avgFuncs.add(func);
							avgChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
							
							maxX = Math.max(maxX, func.getMaxX());
							double halfway = 0.5*func.getMaxX();
							if (halfway >= func.getMinX())
								maxY = Math.max(maxY, func.getInterpolatedY(halfway));
						}
					} else {
						HashSet<Double> processedSubs = new HashSet<>();
						HashSet<Double> processedAvgs = new HashSet<>();
						for (Cell<Double, Double, DiscretizedFunc> cell : cellsSorted(table)) {
							Double sub = cell.getRowKey();
							Double avg = cell.getColumnKey();
							Color subColor = subDeltaColors.get(sub);
							Color avgColor = avgDeltaColors.get(avg);
							Preconditions.checkNotNull(subColor);
							Preconditions.checkNotNull(avgColor);
							DiscretizedFunc func = cell.getValue();
							
							if (func.size() > maxPointsWithoutDecimation) {
								int mod = (int)(0.5*(double)func.size()/(double)maxPointsWithoutDecimation);
								if (mod > 1) {
									DiscretizedFunc decimated = new ArbitrarilyDiscretizedFunc();
									for (int i=0; i<func.size(); i++) {
										if (i % mod == 0 || i == func.size()-1)
											decimated.set(func.get(i));
									}
									func = decimated;
								}
							}
							
							double halfway = 0.5*func.getMaxX();
							if (halfway >= func.getMinX())
								maxY = Math.max(maxY, func.getInterpolatedY(halfway));
							
							maxX = Math.max(maxX, func.getMaxX());
							
							maxAny = Math.max(maxAny, func.getMaxY());
							if (func.getMinY() > 0d)
								minAny = Math.min(minAny, func.getMinY());
							
							DiscretizedFunc subFunc = func;
							if (!processedSubs.contains(sub)) {
								subFunc = subFunc.deepClone();
								subFunc.setName("Sub="+oDF.format(sub));
								processedSubs.add(sub);
							}
							
							DiscretizedFunc avgFunc = func;
							if (!processedAvgs.contains(avg)) {
								avgFunc = avgFunc.deepClone();
								avgFunc.setName("Avg="+oDF.format(avg));
								processedAvgs.add(avg);
							}
							
							subFuncs.add(subFunc);
							subChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, subColor));
							avgFuncs.add(avgFunc);
							avgChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, avgColor));
						}
						
						if (minAny*0.95 > minY)
							minY = minAny*0.95;
						if (maxY > maxAny*1.05)
							maxY = maxAny*1.05;
					}
					
					System.out.println("Spread? "+isSpread+", yRange=["+minY+", "+maxY+"]");
					Range xRange = new Range(minX, maxX);
					Range yRange = new Range(minY, maxY);

					PlotSpec subSpec = new PlotSpec(subFuncs, subChars, title, xAxisLabel, yAxisLabel);
					subSpec.setLegendInset(true);
					PlotSpec avgSpec = new PlotSpec(avgFuncs, avgChars, title, xAxisLabel, yAxisLabel);
					avgSpec.setLegendInset(true);
					
					List<PlotSpec> specs = List.of(subSpec, avgSpec);
					List<Range> xRanges = List.of(xRange);
					List<Range> yRanges = List.of(yRange, yRange);
					
					HeadlessGraphPanel gp = PlotUtils.initHeadless();
					
					gp.drawGraphPanel(specs, false, isSpread, xRanges, yRanges);
					
					PlotUtils.writePlots(outputDir, prefix, gp, 1000, 1400, true, true, false);
				}
			}
		}
	}
	
	private static List<Cell<Double, Double, DiscretizedFunc>> cellsSorted(Table<Double, Double, DiscretizedFunc> table) {
		List<Cell<Double, Double, DiscretizedFunc>> cells = new ArrayList<>(table.cellSet());
		cells.sort(new Comparator<Cell<Double, Double, DiscretizedFunc>>() {

			@Override
			public int compare(Cell<Double, Double, DiscretizedFunc> o1, Cell<Double, Double, DiscretizedFunc> o2) {
				double sub1 = o1.getRowKey();
				double sub2 =o2.getRowKey();
				if (sub1 != sub2)
					return Double.compare(sub1, sub2);
				double avg1 = o1.getColumnKey();
				double avg2 = o2.getColumnKey();
				if (avg1 != avg2)
					return Double.compare(avg1, avg2);
				return 0;
			}
		});
		
		return cells;
	}
	
	private static DiscretizedFunc averageFuncs(Table<Double, Double, DiscretizedFunc> table, double val, boolean isRow) {
		List<DiscretizedFunc> matching = new ArrayList<>();
		if (isRow)
			matching.addAll(table.row(val).values());
		else
			matching.addAll(table.column(val).values());
		double maxMinX = Double.NEGATIVE_INFINITY;
		double minMaxX = Double.POSITIVE_INFINITY;
		for (DiscretizedFunc func : matching) {
			maxMinX = Math.max(maxMinX, func.getMinX());
			minMaxX = Math.min(minMaxX, func.getMaxX());
		}
		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(maxMinX, minMaxX, 100);
		for (int i=0; i<ret.size(); i++) {
			double x = ret.getX(i);
			double y = 0d;
			for (DiscretizedFunc func : matching) {
				if ((float)x == (float)func.getMinX())
					y += func.getY(0);
				else if ((float)x == (float)func.getMaxX())
					y += func.getY(func.size()-1);
				else
					y += func.getInterpolatedY(x);
			}
			y /= (double)matching.size();
			ret.set(i, y);
		}
		return ret;
	}
	
	private static Map<Double, Color> getValColors(Set<Double> valSet, CPT rainbow) {
		List<Double> vals = new ArrayList<>(valSet);
		Collections.sort(vals);
		rainbow = rainbow.rescale(0d, Double.max(1d, vals.size()-1d));
		Map<Double, Color> ret = new HashMap<>();
		for (int i=0; i<vals.size(); i++) {
			Color color = rainbow.getColor((float)i);
			color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 127);
			ret.put(vals.get(i), color);
		}
		return ret;
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");

}
