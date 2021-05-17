package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import javax.swing.SwingUtilities;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.collect.Lists;

public class ResultPlotter {
	
	private static final int max_curve_pts = 800;
	private static final double pDiffMins = 30;
	private static final double pDiffIters = 180000*10;
	private static final boolean do_extrap = true;
	private static final boolean include_min_max = true;
	private static final boolean sort_by_size = true;
	
	private static final boolean thread_speed_csv = true;
	
	static final String time_label = "Time (minutes)";
	static final String serial_time_label = "Serial Time (minutes)";
	static final String parallel_time_label = "Parallel Time (minutes)";
	static final String energy_label = "Energy";
	static final String iterations_label = "Iterations";
	static final String parallel_iterations_label = "Tot. Parallel Iterations";
	static final String time_speedup_label = "Time (serial) / Time (parallel)";
	static final String energy_speedup_label = "Energy (serial) / Energy (parallel)";
	static final String std_dev_label = "Std. Dev. (Energy)";
	static final String improvement_label = "% Improvement";
	static final String threads_label = "# Threads";
	
	static final String avg_energy_vs_time_title = "Averaged Energy Vs Time (m)";
	static final String improvement_vs_time_title = "% Improvement (over "+pDiffMins+" min invervals)";
	static final String improvement_vs_iterations = "% Improvement (over previous "+pDiffIters+" iterations)";
	static final String time_speedup_vs_energy_title = "Time Speedup Vs Energy";
	static final String time_speedup_vs_time_title = "Parallel Speedup Vs Time";
	static final String time_comparison_title = "Serial Time Vs Parallel Time";
	static final String energy_speedup_vs_time_title = "Energy Speedup Vs Time";
	static final String avg_energy_vs_iterations_title = "Averaged Energy Vs Iterations";
	static final String energy_vs_iterations_title = "Energy Vs Iterations";
	static final String energy_vs_time_title = "Energy Vs Time (m)";
	static final String std_dev_vs_time_title = "Std. Dev. Vs Time (m)";
	static final String iterations_vs_time_title = "Iterations Vs Time (m)";
	static final String parallel_iterations_vs_time_title = "Parallel Iterations Vs Time (m)";
	static final String energy_vs_parallel_iterations_title = "Energy vs Parallel Iterations";
	static final String speedup_vs_threads_title = "Speedup vs Threads";
	
	private static ArbitrarilyDiscretizedFunc[] loadCSV (
			File file, int targetPPM) throws IOException {
		String name = file.getName();
		CSVFile<String> csv = CSVFile.readFile(file, true);
		ArbitrarilyDiscretizedFunc energyVsIter = new ArbitrarilyDiscretizedFunc();
		energyVsIter.setName("Energy Vs Iterations ("+name+")");
		energyVsIter.setYAxisName("Energy");
		energyVsIter.setXAxisName("Iterations");
		ArbitrarilyDiscretizedFunc energyVsTime = new ArbitrarilyDiscretizedFunc();
		energyVsTime.setName("Energy Over Time ("+name+")");
		energyVsTime.setYAxisName("Energy");
		energyVsTime.setXAxisName("Time (m)");
		ArbitrarilyDiscretizedFunc iterVsTime = new ArbitrarilyDiscretizedFunc();
		iterVsTime.setName("Iterations Over Time ("+name+")");
		iterVsTime.setYAxisName("Iterations");
		iterVsTime.setXAxisName("Time (m)");
		
		System.out.print("Loading: "+name+" ...");
		
		int rows = csv.getNumRows();
		long totMilis = Long.parseLong(csv.getLine(rows-1).get(1));
		double totMins = totMilis / 1000d / 60d;
		double targetPts = targetPPM * totMins;
		if (targetPts > max_curve_pts)
			targetPts = max_curve_pts;
		double newMod = (double)rows/targetPts;
		
		int mod = (int)(newMod + 0.5);
		if (mod < 1)
			mod = 1;
		
		for (int i=1; i<rows; i++) {
			if (i % mod > 0)
				continue;
			
			List<String> line = csv.getLine(i);
			
//			for (String val : line)
//				System.out.println(val);
			
			long iter = Long.parseLong(line.get(0));
			long millis = Long.parseLong(line.get(1));
			double energy = Double.parseDouble(line.get(2));
			
			double secs = millis / 1000d;
			double mins = secs / 60d;
			
			energyVsIter.set((double)iter, energy);
			energyVsTime.set((double)mins, energy);
			iterVsTime.set((double)mins, (double)iter);
		}
		ArbitrarilyDiscretizedFunc[] ret = { energyVsIter, energyVsTime, iterVsTime };
		
		System.out.println("DONE (size="+energyVsTime.size()+")");
		
		return ret;
	}
	
	private static int getIntForKey(String name, String key) {
		String[] splits = name.split("_");
		for (String split : splits) {
			if (split.contains(key)) {
				split = split.substring(0, split.indexOf(key));
				try {
					return Integer.parseInt(split);
				} catch (NumberFormatException e) {
					return 0;
				}
			}
		}
		return 0;
	}
	
	private static int getSizeScore(String name) {
		if (name == null || name.isEmpty()) {
			return 0;
		} else {
			if (name.contains("dsa")) {
				int score =  10 + getIntForKey(name, "nodes");
				if (name.contains("EXTRAPOLATED"))
					score--;
				if (name.contains("minimum") || name.contains("maximum"))
					score--;
				return score;
			} else {
				return getIntForKey(name, "threads");
			}
		}
	}
	
	private static ArrayList<DiscretizedFunc> getParallelIters(ArrayList<DiscretizedFunc> funcs, boolean itersIsX) {
		ArrayList<DiscretizedFunc> ret = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc func : funcs) {
			String name = func.getName();
			if (name == null || name.length()<10 || !(name.contains("dsa") || name.contains("tsa"))) {
				ret.add(new ArbitrarilyDiscretizedFunc());
				continue;
			}
			int threads = getIntForKey(name, "threads");
			if (threads < 1) {
				ret.add(new ArbitrarilyDiscretizedFunc());
				continue;
			}
			int nodes = getIntForKey(name, "nodes");
			if (nodes == 0)
				nodes = 1;
			
			threads *= nodes;
			
			ArbitrarilyDiscretizedFunc newFunc = new ArbitrarilyDiscretizedFunc();
			newFunc.setName(name);
			
			for (int i=0; i<func.size(); i++) {
				double iters, other;
				if (itersIsX) {
					iters = func.getX(i);
					other = func.getY(i);
				} else {
					iters = func.getY(i);
					other = func.getX(i);
				}
				
				iters *= (double)threads;
				
				if (itersIsX)
					newFunc.set(iters, other);
				else
					newFunc.set(other, iters);
			}
			
			ret.add(newFunc);
		}
		
		return ret;
	}
	
	protected static GraphWindow getGraphWindow(ArrayList<? extends DiscretizedFunc> funcs, String title,
			ArrayList<PlotCurveCharacterstics> chars, String xAxisName, String yAxisName, boolean visible) {
		if (sort_by_size) {
			ArrayList<DiscretizedFunc> sortedFuncs = new ArrayList<DiscretizedFunc>();
			ArrayList<PlotCurveCharacterstics> sortedChars = new ArrayList<PlotCurveCharacterstics>();
			ArrayList<Integer> sortedSizes = new ArrayList<Integer>();
			for (int i=0; i<funcs.size(); i++) {
				String name = funcs.get(i).getName();
				int size = getSizeScore(name);
				int insertion = -1;
				for (int j=0; j<sortedSizes.size(); j++) {
					if (sortedSizes.get(j) < size) {
						insertion = j;
						break;
					}
				}
				if (insertion < 0)
					insertion = sortedSizes.size();
				sortedFuncs.add(insertion, funcs.get(i));
				sortedChars.add(insertion, chars.get(i));
				sortedSizes.add(insertion, size);
			}
			funcs = sortedFuncs;
			chars = sortedChars;
		}
		
		GraphWindow gwAPI = new GraphWindow(funcs, title, chars);
		GraphWidget gw = gwAPI.getGraphWidget();
		gw.setPlotLabelFontSize(30);
		gw.setAxisLabelFontSize(18);
		gw.setTickLabelFontSize(14);
		gw.setPlottingOrder(DatasetRenderingOrder.REVERSE);
//		gw.getGraphPanel().setBackgroundColor(Color.WHITE);
		gw.setXAxisLabel(xAxisName);
		gw.setYAxisLabel(yAxisName);
		gw.setSize(1000, 1000);
		return gwAPI;
	}
	
	private static ArbitrarilyDiscretizedFunc[] asArray(List<ArbitrarilyDiscretizedFunc> funcs) {
		return funcs.toArray(new ArbitrarilyDiscretizedFunc[funcs.size()]);
	}
	
	public static EvenlyDiscretizedFunc calcAvg(List<ArbitrarilyDiscretizedFunc> funcs, int numX) {
		return calcAvgAndStdDev(funcs, numX, false)[0];
	}
	
	public static EvenlyDiscretizedFunc[] calcAvgAndStdDev(List<ArbitrarilyDiscretizedFunc> funcs, int numX,
			boolean includeStdDev) {
		ArbitrarilyDiscretizedFunc[] funcsArray = asArray(funcs);
		double largestMin = Double.MIN_VALUE;
		double smallestMax = Double.MAX_VALUE;
		
		for (ArbitrarilyDiscretizedFunc func : funcsArray) {
			double min = func.getMinX();
			double max = func.getMaxX();
			if (min > largestMin)
				largestMin = min;
			if (max < smallestMax)
				smallestMax = max;
		}
		
		largestMin += 0.001;
		smallestMax -= 0.001;
		
		int num = funcsArray.length;
		
		EvenlyDiscretizedFunc avg = new EvenlyDiscretizedFunc(largestMin, smallestMax, numX);
		EvenlyDiscretizedFunc stdDevs = null;
		if (includeStdDev)
			stdDevs = new EvenlyDiscretizedFunc(largestMin, smallestMax, numX);
		
		EvenlyDiscretizedFunc minFunc = null, maxFunc = null;
		if (include_min_max) {
			minFunc = new EvenlyDiscretizedFunc(largestMin, smallestMax, numX);
			maxFunc = new EvenlyDiscretizedFunc(largestMin, smallestMax, numX);
		}
		
		double[] values = new double[num];
		for (int i=0; i<numX; i++) {
			double x = avg.getX(i);
			for (int j=0; j<num; j++) {
//				values[j] = funcsArray[j].getInterpolatedY_inLogYDomain(x);
				values[j] = funcsArray[j].getInterpolatedY(x);
			}
			double mean = StatUtils.mean(values);
			
			avg.set(i, mean);
			if (includeStdDev) {
				double stdDev = Math.sqrt(StatUtils.variance(values, mean));
				stdDevs.set(i, stdDev);
			}
			if (include_min_max) {
				minFunc.set(i, StatUtils.min(values));
				maxFunc.set(i, StatUtils.max(values));
			}
		}
		EvenlyDiscretizedFunc[] ret = { avg, stdDevs, minFunc, maxFunc };
		return ret;
	}
	
	private static ArbitrarilyDiscretizedFunc getSwapped(DiscretizedFunc func) {
		ArbitrarilyDiscretizedFunc swapped = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<func.size(); i++)
			swapped.set(func.getY(i), func.getX(i));
		return swapped;
	}
	
	private static ArrayList<DiscretizedFunc> generateEnergySpeedup(
			List<DiscretizedFunc> funcs, DiscretizedFunc ref, int numX) {
		double refMin = ref.getMinY();
		double refMax = ref.getMaxY();
		
		DiscretizedFunc refSwapped = getSwapped(ref);
		
		ArrayList<DiscretizedFunc> speedups = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc func : funcs) {
			DiscretizedFunc funcSwapped = getSwapped(func);
			
			double min = funcSwapped.getMinX();
			double max = funcSwapped.getMaxX();
			if (min < refMin)
				min = refMin;
			min += 0.001;
			if (max > refMax)
				max = refMax;
			max -= 0.001;
			
			EvenlyDiscretizedFunc speedupFunc = new EvenlyDiscretizedFunc(min, max, numX);
			speedupFunc.setName(func.getName());
			
			for (int i=0; i<numX; i++) {
				double energy = speedupFunc.getX(i);
				double refTime = refSwapped.getInterpolatedY(energy);
				double myTime = funcSwapped.getInterpolatedY(energy);
				
				double speedup = refTime / myTime;
				speedupFunc.set(energy, speedup);
			}
			
			speedups.add(speedupFunc);
		}
		
		return speedups;
	}
	
	private static double getLowestEnergy(List<DiscretizedFunc> funcs) {
		double min = Double.MAX_VALUE;
		for (DiscretizedFunc func : funcs) {
			double minE = func.getY(func.size()-1);
			if (minE < min)
				min = minE;
		}
		return min;
	}
	
//	private static DiscretizedFunc getExtrapolatedRef(DiscretizedFunc ref, double targetEnergy) {
//		double energy = ref.getY(ref.getNum()-1);
//		if (targetEnergy > energy)
//			// ref func already goes low enough!
//			return null;
//		
//		ArbitrarilyDiscretizedFunc arbRef;
//		if (ref instanceof ArbitrarilyDiscretizedFunc) {
//			arbRef = (ArbitrarilyDiscretizedFunc)ref;
//		} else {
//			arbRef = new ArbitrarilyDiscretizedFunc();
//			for (int i=ref.getNum()-100; i<ref.getNum(); i++) {
//				if (i < 0)
//					continue;
//				arbRef.set(ref.get(i));
//			}
//		}
//		
//		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
//		double startTime = ref.getX(ref.getNum()-1);
//		double timeDelta = startTime - ref.getX(ref.getNum()-2);
//		
//		double maxTime = startTime * 4;
//		for (double time=startTime; energy>targetEnergy && time < maxTime; time += timeDelta) {
//			energy = arbRef.getInterpExterpY_inLogXLogYDomain(time);
//			ret.set(time, energy);
//		}
//		return ret;
//	}
	
	private static final double extrap_slope_pt_fraction = 0.90;
	
	private static DiscretizedFunc getExtrapolatedRef(DiscretizedFunc ref, double targetEnergy) {
		double energy = ref.getY(ref.size()-1);
		if (targetEnergy > energy)
			// ref func already goes low enough!
			return null;
		
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		ret.setName(ref.getName()+" (EXTRAPOLATED)");
		double startTime = ref.getX(ref.size()-1);
		double timeDelta = startTime - ref.getX(ref.size()-2);
		
		int i1 = ref.size()-1;
		int i0 = -1;
		double extrapPtTime = startTime * extrap_slope_pt_fraction;
		for (int i=0; i<ref.size(); i++) {
			double time = ref.getX(i);
			if (time > extrapPtTime) {
				i0 = i;
				break;
			}
		}
		
		double x0 = Math.log(ref.getX(i0));
		double x1 = Math.log(ref.getX(i1));
		
		double y0 = Math.log(ref.getY(i0));
		double y1 = Math.log(ref.getY(i1));
		
		double logSlope = (y1 - y0) / (x1 - x0);
		double logIntercept = y1 - (logSlope * x1);
		
		double maxTime = startTime * 4;
		for (double time=startTime; energy>targetEnergy && time < maxTime; time += timeDelta) {
//			double logEnergy = (logSlope * Math.log(time)) + logIntercept;
			energy = Math.exp(logSlope * Math.log(time) + logIntercept);
//			energy = arbRef.getInterpExterpY_inLogXLogYDomain(time);
			ret.set(time, energy);
		}
		return ret;
	}
	
	private static ArrayList<DiscretizedFunc> generateEnergyTimeSpeedup(
			List<DiscretizedFunc> energyTimeComparisons) {
		
		ArrayList<DiscretizedFunc> speedups = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc comp : energyTimeComparisons) {
			ArbitrarilyDiscretizedFunc speedupFunc = new ArbitrarilyDiscretizedFunc();
			speedupFunc.setName(comp.getName());
			
			for (int i=0; i<comp.size(); i++) {
				double time = comp.getX(i);
				double refTime = comp.getY(i);
				
				double speedup = refTime / time;
				speedupFunc.set(time, speedup);
			}
			
			speedups.add(speedupFunc);
		}
		
		return speedups;
	}
	
	private static ArrayList<DiscretizedFunc> generateEnergyTimeComparison(
			List<DiscretizedFunc> funcs, DiscretizedFunc ref, int numX) {
		double refMinEnergy = ref.getMinY();
		double refMaxEnergy = ref.getMaxY();
		
		ArrayList<DiscretizedFunc> speedups = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc func : funcs) {
			double funcMinEnergy = func.getMinY();
			double funcMaxEnergy = func.getMaxY();
			
			if (funcMinEnergy > refMaxEnergy) {
				speedups.add(null);
				continue;
			}
			
			double maxEnergy = funcMaxEnergy;
			if (maxEnergy > refMaxEnergy)
				maxEnergy = refMaxEnergy;
			double minEnergy = refMinEnergy;
			if (minEnergy < funcMinEnergy)
				minEnergy = funcMinEnergy;
			
			if (Double.isNaN(minEnergy) || Double.isNaN(maxEnergy)) {
				speedups.add(null);
				continue;
			}
			
			double minTime = func.getFirstInterpolatedX(maxEnergy);
			minTime += 0.001;
			double maxTime = func.getFirstInterpolatedX(minEnergy);
			maxTime -= 0.001;
			
			if (minTime >= maxTime) {
				speedups.add(null);
				continue;
			}
			
			EvenlyDiscretizedFunc speedupFunc = new EvenlyDiscretizedFunc(minTime, maxTime, numX);
			speedupFunc.setName(func.getName());
			
			for (int i=0; i<numX; i++) {
				double time = speedupFunc.getX(i);
				double energy = func.getInterpolatedY(time);
				double refTime = ref.getFirstInterpolatedX(energy);
				speedupFunc.set(time, refTime);
				
//				double speedup = refTime / time;
//				speedupFunc.set(time, speedup);
			}
			
			speedups.add(speedupFunc);
		}
		
		return speedups;
	}
	
	private static ArrayList<DiscretizedFunc> generateTimeSpeedup(
			List<DiscretizedFunc> funcs, DiscretizedFunc ref, int numX) {
		double refMin = ref.getMinX();
		double refMax = ref.getMaxX();
		
		ArrayList<DiscretizedFunc> speedups = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc func : funcs) {
			double min = func.getMinX();
			double max = func.getMaxX();
			if (min < refMin)
				min = refMin;
			if (max > refMax)
				max = refMax;
			
			EvenlyDiscretizedFunc speedupFunc = new EvenlyDiscretizedFunc(min, max, numX);
			speedupFunc.setName(func.getName());
			
			for (int i=0; i<numX; i++) {
				double time = speedupFunc.getX(i);
				double refEnergy = ref.getInterpolatedY(time);
				double myEnergy = func.getInterpolatedY(time);
				
				double speedup = refEnergy / myEnergy;
				speedupFunc.set(time, speedup);
			}
			
			speedups.add(speedupFunc);
		}
		
		return speedups;
	}
	
	private static ArrayList<DiscretizedFunc> generateSpeedupVsThreads(File dir, ArrayList<DiscretizedFunc> speedupVsTime) throws IOException {
		double smallestMaxTime = Double.MAX_VALUE;
		
		for (DiscretizedFunc func : speedupVsTime) {
			if (func.getMaxX() < smallestMaxTime)
				smallestMaxTime = func.getMaxX();
		}
		
		ArbitrarilyDiscretizedFunc spd = new ArbitrarilyDiscretizedFunc();
		spd.setName(speedup_vs_threads_title);
		
		int maxThreads = 0;
		
		for (DiscretizedFunc func : speedupVsTime) {
			String name = func.getName();
			System.out.println("Processing threads for func: "+name);
			if (name == null || name.length()<10 || !(name.contains("dsa") || name.contains("tsa") || name.startsWith("FM"))) {
				continue;
			}
			if (name.contains("min") || name.contains("max"))
				continue;
			int threads = getIntForKey(name, "threads");
			if (threads < 1)
				continue;
				
			int nodes = getIntForKey(name, "nodes");
			if (nodes == 0)
				nodes = 1;
			
			threads *= nodes;
			
			if (threads > maxThreads)
				maxThreads = threads;
			
			double maxSpeedup = func.getInterpolatedY(smallestMaxTime);
			
			try {
				double prev = func.getY(smallestMaxTime);
				if (prev > maxSpeedup) {
					System.out.println("already have a better speedup for "+threads+" threads, keeping that.");
					continue;
				} else {
					System.out.println("already have a worse speedup for "+threads+" threads, replacing.");
				}
			} catch (Exception e) {}
			
			System.out.println("Adding thread speedup for "+threads+" threads: "+name);
			
			spd.set((double)threads, maxSpeedup);
		}
		
		ArrayList<DiscretizedFunc> spds = new ArrayList<DiscretizedFunc>();
		spds.add(spd);

		ArbitrarilyDiscretizedFunc linearFunc = new ArbitrarilyDiscretizedFunc();
		linearFunc.setName("Linear Speedup");
		ArbitrarilyDiscretizedFunc sqrtFunc = new ArbitrarilyDiscretizedFunc();
		sqrtFunc.setName("Sqrt Speedup");
		for (int threads=1; threads<maxThreads; threads++) {
			double dthreads = (double)threads;
			linearFunc.set(dthreads, dthreads);
			
			sqrtFunc.set(dthreads, Math.sqrt(dthreads));
		}
		
		spds.add(linearFunc);
		spds.add(sqrtFunc);
		
		// write speedup file
		if (thread_speed_csv && dir != null) {
			CSVFile<String> csv = new CSVFile<String>(true);
			csv.addLine(Lists.newArrayList("Threads", "Speedup"));
			for (int i=0; i<spd.size(); i++) {
				Integer threads = (int)spd.getX(i);
				Double speedup = spd.getY(i);
				
				csv.addLine(Lists.newArrayList(threads.toString(), speedup.toString()));
			}
			csv.writeToFile(new File(dir, "spd_vs_threads.csv"));
		}
		
		return spds;
	}
	
	private static ArrayList<DiscretizedFunc> generatePercentImprovementOverTime(
			List<DiscretizedFunc> funcs, double mins) {
		ArrayList<DiscretizedFunc> ret = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc func : funcs) {
			// x is time in m
			// y is energy
			
			DiscretizedFunc retFunc = new ArbitrarilyDiscretizedFunc();
			
			double start = func.getMinX()+mins + 0.000001;
			double delta = (func.getMaxX() - start) / 1000d;
			
			for (double time=start; time<func.getMaxX(); time += delta) {
				double energy = func.getInterpolatedY(time);
				
				double prevEnergy = func.getInterpolatedY(time - mins);
				double deltaE = prevEnergy - energy;
				double improvement = deltaE / prevEnergy;
				double percent = improvement * 100d;
				
				prevEnergy = energy;
				
				retFunc.set(time, percent);
			}
			
			ret.add(retFunc);
		}
		
		return ret;
	}
	
	private static ArrayList<DiscretizedFunc> generatePercentImprovementVsIterations(
			List<DiscretizedFunc> funcs, double iterations) {
		ArrayList<DiscretizedFunc> ret = new ArrayList<DiscretizedFunc>();
		
		for (DiscretizedFunc func : funcs) {
			// x is time in m
			// y is energy
			
			DiscretizedFunc retFunc = new ArbitrarilyDiscretizedFunc();
			
			double start = func.getMinX()+iterations;
			double delta = (func.getMaxX() - start) / 1000d;
			
			for (double iter=start; iter<func.getMaxX(); iter += delta) {
				double energy = func.getInterpolatedY(iter);
				
				double prevEnergy = func.getInterpolatedY(iter - iterations);
				double deltaE = prevEnergy - energy;
				double improvement = deltaE / prevEnergy;
				double percent = improvement * 100d;
				
				prevEnergy = energy;
				
				retFunc.set(iter, percent);
			}
			
			ret.add(retFunc);
		}
		
		return ret;
	}
	
	private static Color getSaturated(Color c, float factor) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		
		if ( r == 0 && g == 0 && b == 0) {
			return Color.GRAY;
//			float val = (float)r * factor / 255f;
//			int adjusted = (int)(val + 0.5f);
//			if (adjusted < 0)
//				adjusted = 0;
//			if (adjusted > 255)
//				adjusted = 255;
//			return new Color(adjusted, adjusted, adjusted);
		} else {
			float[] hsb = Color.RGBtoHSB(r, g, b, null);
			int rgb = Color.HSBtoRGB(hsb[0], hsb[1]*factor, hsb[2]);
			
			return new Color(rgb);
		}
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws InvocationTargetException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException, InvocationTargetException {
		
		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/test_inversion/bench/");
//		File mainDir = new File("D:\\Documents\\temp\\Inversion Results");
		
		final File tsaDir = null;
//		final File dsaDir = null;
		
//		tsaDir = new File(mainDir, "results_5");
//		tsaDir = new File(mainDir, "results_6");
//		dsaDir = new File(mainDir, "dsa_results_8");
//		dsaDir = new File(mainDir, "mult_state_1_7hrs");
//		dsaDir = new File(mainDir, "mult_ncal_1");
//		dsaDir = new File(mainDir, "mult_ncal_2");
//		dsaDir = new File(mainDir, "multi/ncal_1");
//		dsaDir = new File(mainDir, "mult_state_2_comb_3");
//		dsaDir = new File(mainDir, "multi/ranger_ncal_1");
//		dsaDir = new File(mainDir, "2011_10_19-threads_test");
//		dsaDir = new File(mainDir, "2011_10_20-ncal-bench-orig_single");
//		dsaDir = new File(mainDir, "2011_10_21-ncal-bench-sub200");
//		dsaDir = new File(mainDir, "2011_10_21-ncal-bench-sub200-orig_single");
//		dsaDir = new File(mainDir, "2011_10_27-ncal-bench-sub-secs-test");
//		dsaDir = new File(mainDir, "2011_10_31-allcal-bench-sub-secs-test");
//		dsaDir = new File(mainDir, "2012_02_22-model2-bench/8threads");
//		dsaDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2012_02_27-unconstrained");
//		final File dsaDir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_04-threaded-bench/csvs/char");
		final File dsaDir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_04_19-threaded-bench/csvs");
		
//		List<String> bundleKeys = null;
		final List<String> bundleKeys = Lists.newArrayList("VarNormal", "VarKeepAsBest", "VarSingleThreadKeep",
				"VarSingleThreadNoKeep", "VarSub5_1", "VarSub5_0.3", "VarSub10_1");
		
//		dsaDir = new File(mainDir, "agu/ncal_constrained");
//		dsaDir = new File(mainDir, "agu/ncal_unconstrained");
//		dsaDir = new File(mainDir, "agu/allcal_constrained");
//		dsaDir = new File(mainDir, "agu/allcal_unconstrained");
		
//		dsaDir = new File(mainDir, "2011_09_07_morgan_NoCS_UCERF2MagDist");
//		dsaDir = new File(mainDir, "2011_09_16_genetic_test");
		
		final ArrayList<String> plots = new ArrayList<String>();
		plots.add(energy_vs_time_title);
//		plots.add(improvement_vs_iterations);
		plots.add(avg_energy_vs_time_title);
//		plots.add(std_dev_vs_time_title);
//		plots.add(improvement_vs_time_title);
//		plots.add(time_speedup_vs_energy_title);
//		plots.add(time_comparison_title);
		plots.add(time_speedup_vs_time_title);
		plots.add(speedup_vs_threads_title);
//		plots.add(speedup_vs_time_title);
//		plots.add(energy_vs_iterations_title);
//		plots.add(avg_energy_vs_iterations_title);
//		plots.add(energy_vs_parallel_iterations_title);
//		plots.add(parallel_iterations_vs_time_title);
//		plots.add(iterations_vs_time_title);
		
//		String highlight = "Unconst";
		final String highlight = null;
		
//		highlight = "dsa_8threads_10nodes_FAST_SA_dSub200_sub100";
		
//		String coolType = "VERYFAST";
		final String coolType = null;
		final int threads = -1;
//		int threads = 6;
		final int nodes = -1;
		final boolean includeStartSubZero = false;
		final boolean plotAvg = true;
		final boolean bundleDsaBySubs = false;
		final boolean bundleTsaBySubs = false;
		
		CPT tempCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		CPTVal last = tempCPT.remove(tempCPT.size()-1);
		final CPT threadCPT = tempCPT.rescale(1, 8);
		last.start = 8f;
		last.end = 23f;
		threadCPT.add(last);
		
		final int avgNumX = 400;
		final int targetPPM = 2;
		
		SwingUtilities.invokeAndWait(new Runnable() {
			
			@Override
			public void run() {
				try {
					generatePlots(tsaDir, dsaDir, highlight, coolType, threads, nodes,
							includeStartSubZero, plotAvg, bundleDsaBySubs, bundleTsaBySubs,
							avgNumX, targetPPM, true, plots, threadCPT, bundleKeys);
				} catch (Exception e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
			}
		});
	}

	protected static HashMap<String, GraphWindow> generatePlots(File tsaDir, File dsaDir,
			String highlight, String coolType, int threads, int nodes,
			boolean includeStartSubZero, boolean plotAvg,
			boolean bundleDsaBySubs, boolean bundleTsaBySubs, int avgNumX, int targetPPM,
			boolean visible, Collection<String> plots, CPT threadCPT, List<String> bundleKeys)
			throws IOException {
		File[] tsaFiles;
		if (tsaDir == null)
			tsaFiles = new File[0];
		else
			tsaFiles = tsaDir.listFiles();
		File[] dsaFiles;
		if (dsaDir == null)
			dsaFiles = new File[0];
		else
			dsaFiles = dsaDir.listFiles();
		
		File[] files = new File[tsaFiles.length+dsaFiles.length];
		System.arraycopy(tsaFiles, 0, files, 0, tsaFiles.length);
		System.arraycopy(dsaFiles, 0, files, tsaFiles.length, dsaFiles.length);
		
		ArrayList<DiscretizedFunc> energyVsIter = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> energyVsTime = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> iterVsTime = new ArrayList<DiscretizedFunc>();
		
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<PlotCurveCharacterstics>();
		
		HashMap<String, ArrayList<ArbitrarilyDiscretizedFunc>> runsEnergyTimeMap =
			new HashMap<String, ArrayList<ArbitrarilyDiscretizedFunc>>();
		HashMap<String, ArrayList<ArbitrarilyDiscretizedFunc>> runsIterTimeMap =
			new HashMap<String, ArrayList<ArbitrarilyDiscretizedFunc>>();
		HashMap<String, ArrayList<ArbitrarilyDiscretizedFunc>> runsEnergyIterMap =
			new HashMap<String, ArrayList<ArbitrarilyDiscretizedFunc>>();
		HashMap<String, PlotCurveCharacterstics> runsChars = new HashMap<String, PlotCurveCharacterstics>();
		
		DiscretizedFunc refFunc = null;
		CPT bundleCPT = null;
		
		for (File file : files) {
			String name = file.getName();
			if (!name.endsWith(".csv"))
				continue;
			if (name.contains("spd_vs_threads"))
				continue;
			if (highlight != null && !name.contains(highlight))
				continue;
			
//			if (name.contains("dSub600"))
//				continue;
			
			if (coolType != null && !name.contains(coolType))
				continue;
			if (threads > 0 && !name.contains(threads+"thread"))
				continue;
			if (nodes == 1 && name.contains("nodes"))
				continue;
			if (nodes > 1 && !name.contains(nodes+"nodes"))
				continue;
			if (!includeStartSubZero && name.contains("startSub"))
				continue;
			if (tsaDir != null && name.startsWith("tsa") && file.getAbsolutePath().contains("mult_"))
				continue;
			System.out.println("Loading CSV: "+name);
			ArbitrarilyDiscretizedFunc[] funcs = loadCSV(file, targetPPM);
			
			energyVsIter.add(funcs[0]);
			energyVsTime.add(funcs[1]);
			iterVsTime.add(funcs[2]);
			
			PlotLineType type;
			if (name.contains("CLASSICAL_SA"))
				type = PlotLineType.DOTTED;
			else if (name.contains("VERYFAST_SA"))
				type = PlotLineType.DASHED;
			else
				type = PlotLineType.SOLID;
			
			float size = 1f;
			
			Color c;
			if (name.contains("dsa")) {
				if (bundleDsaBySubs && name.contains("dSub")) {
					if (name.contains("dSub25_") || name.contains("dSub500mi_"))
						c = Color.BLACK;
					else if (name.contains("dSub50_") || name.contains("dSub1s_"))
						c = Color.BLUE;
					else if (name.contains("dSub100_") || name.contains("dSub2500mi_"))
						c = Color.GREEN;
					else if (name.contains("dSub200_") || name.contains("dSub5s_"))
						c = Color.YELLOW;
					else if (name.contains("dSub500_") || name.contains("dSub10s_"))
						c = Color.ORANGE;
					else
						c = Color.MAGENTA;
					
					if (name.contains("dSub500mi_") || name.contains("dSub1s_")
							|| name.contains("dSub2500mi_") || name.contains("dSub5s_")
							|| name.contains("dSub10s_"))
						type = PlotLineType.DASHED;
					
//					if (name.contains("dSub25_") || name.contains("dSub500mi_"))
//						c = Color.BLACK;
//					else if (name.contains("dSub50_") || name.contains("dSub1s_"))
//						c = Color.BLUE;
//					else if (name.contains("dSub100_") || name.contains("dSub2s_"))
//						c = Color.GREEN;
//					else if (name.contains("dSub200_") || name.contains("dSub3s_"))
//						c = Color.YELLOW;
//					else if (name.contains("dSub500_") || name.contains("dSub4s_"))
//						c = Color.ORANGE;
//					else if (name.contains("dSub5s_"))
//						c = Color.RED;
//					else
//						c = Color.MAGENTA;
//					
//					if (name.contains("dSub500mi_") || name.contains("dSub1s_")
//							|| name.contains("dSub2s_") || name.contains("dSub3s_")
//							|| name.contains("dSub4s_") ||name.contains("dSub5s_"))
//						type = PlotLineType.DASHED;
					
//					if (name.contains("5nodes"))
//						size = 2f;
//					else if (name.contains("10nodes"))
//						size = 3f;
//					else if (name.contains("20nodes"))
//						size = 4f;
					if (name.contains("50nodes"))
						size = 3f;
				} else {
					if (name.contains("2nodes"))
						c = Color.BLACK;
					else if (name.contains("5nodes"))
						c = Color.BLUE;
					else if (name.contains("10nodes"))
						c = Color.GREEN;
					else if (name.contains("20nodes"))
						c = Color.RED;
					else if (name.contains("50nodes"))
						c = Color.ORANGE;
					else if (name.contains("100nodes"))
						c = Color.MAGENTA;
					else
						c = Color.PINK;
					
					if (name.contains("1thread"))
						type = PlotLineType.DASHED;
//					else if (name.contains("4threads"))
//						c = Color.GREEN;
//					else if (name.contains("6threads"))
//						size += 1;
//					else if (name.contains("8threads"))
//						size += 1;
				}
			} else {
				if (bundleTsaBySubs && name.contains("sub")) {
					if (name.contains("sub50_"))
						c = Color.BLACK;
					else if (name.contains("sub100_"))
						c = Color.DARK_GRAY;
					else if (name.contains("sub200_"))
						c = Color.BLUE;
					else if (name.contains("sub400_"))
						c = Color.CYAN;
					else if (name.contains("sub600_"))
						c = Color.GREEN;
//					else if (name.contains("dSub5000_"))
//						c = Color.YELLOW;
//					else if (name.contains("dSub10000_"))
//						c = Color.ORANGE;
//					else if (name.contains("dSub15000_"))
//						c = Color.RED;
					else
						c = Color.MAGENTA;
					
					if (name.contains("4threads"))
						size = 2f;
					else if (name.contains("6threads"))
						size = 3f;
					else if (name.contains("8threads"))
						size = 4f;
				} else if (bundleKeys != null && !bundleKeys.isEmpty()) {
					if (bundleCPT == null) {
						bundleCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
						bundleCPT = bundleCPT.rescale(0, bundleKeys.size()-1);
//						bundleCPT = threadCPT.rescale(0, bundleKeys.size());
						bundleCPT.setNanColor(Color.BLACK);
					}
					float val = Float.NaN;
					for (int i=0; i<bundleKeys.size(); i++) {
						String key = bundleKeys.get(i);
						if (name.contains(key)) {
							val = (float)i;
							break;
						}
					}
					System.out.println("VAL: "+val);
					c = bundleCPT.getColor(val);
					size = 1f;
				} else {
					if (threadCPT != null && name.contains("thread")) {
						String sub = name.substring(0, name.indexOf("thread"));
						sub = sub.substring(sub.lastIndexOf("_")+1);
						int numThreads = Integer.parseInt(sub);
						c = threadCPT.getColor((float)numThreads);
					} else {
						if (name.contains("1thread"))
							c = Color.BLACK;
						else if (name.contains("_2threads"))
							c = Color.BLUE;
						else if (name.contains("4threads"))
							c = Color.GREEN;
						else if (name.contains("6threads"))
							c = Color.MAGENTA;
						else if (name.contains("8threads"))
							c = Color.RED;
						else
							c = Color.ORANGE;
					}
				}
			}
			
			if (highlight != null && name.startsWith(highlight)) {
				c = Color.PINK;
				size *= 2;
			}
			
			if (includeStartSubZero && name.contains("startSubIterationsAtZero"))
				size += 1f;
			
			if (nodes == 1 || name.contains("dsa"))
				size += 1f;
			
			if (plotAvg && name.contains("run")) {
				String shortName = name.substring(0, name.indexOf("run"));
				if (!runsEnergyTimeMap.containsKey(shortName)) {
					runsEnergyTimeMap.put(shortName, new ArrayList<ArbitrarilyDiscretizedFunc>());
					runsIterTimeMap.put(shortName, new ArrayList<ArbitrarilyDiscretizedFunc>());
					runsEnergyIterMap.put(shortName, new ArrayList<ArbitrarilyDiscretizedFunc>());
					runsChars.put(shortName, new PlotCurveCharacterstics(type, size, c));
				}
				
				runsEnergyTimeMap.get(shortName).add(funcs[1]);
				runsIterTimeMap.get(shortName).add(funcs[2]);
				runsEnergyIterMap.get(shortName).add(funcs[0]);
			}
			
			chars.add(new PlotCurveCharacterstics(type, size, c));
		}
		System.gc();
		
		System.out.println("Averaging");
		ArrayList<DiscretizedFunc> averages = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> stdDevs = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> minFuncs = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> maxFuncs = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> iterTimeAvgs = new ArrayList<DiscretizedFunc>();
		ArrayList<DiscretizedFunc> energyIterAvgs = new ArrayList<DiscretizedFunc>();
		ArrayList<PlotCurveCharacterstics> avgChars = new ArrayList<PlotCurveCharacterstics>();
		if (plots.contains(avg_energy_vs_time_title)
				|| plots.contains(std_dev_vs_time_title)
				|| plots.contains(time_speedup_vs_energy_title)
				|| plots.contains(time_comparison_title)
				|| plots.contains(time_speedup_vs_time_title)
				|| plots.contains(avg_energy_vs_iterations_title)
				|| plots.contains(energy_vs_parallel_iterations_title)) {
			for (String name : runsEnergyTimeMap.keySet()) {
				ArrayList<ArbitrarilyDiscretizedFunc> runs = runsEnergyTimeMap.get(name);
				if (runs == null || runs.size() <= 1)
					continue;
				
				String avgName = name + " (average of "+runs.size()+" curves)";
				String maxName = name + " (maximum of "+runs.size()+" curves)";
				String minName = name + " (minimum of "+runs.size()+" curves)";
				String stdDevName = name + " (std dev of "+runs.size()+" curves)";
				
				System.out.println("Averaging: "+name);
				
				DiscretizedFunc[] ret = calcAvgAndStdDev(runs, avgNumX, true);
				DiscretizedFunc avg = ret[0];
				avg.setName(avgName);
				averages.add(avg);
				
				DiscretizedFunc stdDev = ret[1];
				if (stdDev != null) {
					stdDev.setName(stdDevName);
					stdDevs.add(stdDev);
				}
				
				DiscretizedFunc minFunc = ret[2];
				DiscretizedFunc maxFunc = ret[3];
				if (minFunc != null && maxFunc != null) {
					minFuncs.add(minFunc);
					minFunc.setName(minName);
					maxFuncs.add(maxFunc);
					maxFunc.setName(maxName);
				}
				
				DiscretizedFunc iterTimeAvg = calcAvg(runsIterTimeMap.get(name), avgNumX);
				iterTimeAvg.setName(avgName);
				iterTimeAvgs.add(iterTimeAvg);
				
				DiscretizedFunc energyIterAvg = calcAvg(runsEnergyIterMap.get(name), avgNumX);
				energyIterAvg.setName(avgName);
				energyIterAvgs.add(energyIterAvg);
				
				avgChars.add(runsChars.get(name));
				
				if (refFunc == null) {
					if (coolType == null && name.startsWith("tsa_1threads_FAST")
							|| name.startsWith("tsa_1threads_"+coolType))
						refFunc = avg;
					else if (name.startsWith("FM") && name.contains("_1thread"))
						refFunc = avg;
				}
//				if (refFunc == null && name.startsWith("tsa_1threads_VERYFAST"))
//					refFunc = avg;
			}
		}
		
		ArrayList<DiscretizedFunc> timeComparisons = null;
		ArrayList<PlotCurveCharacterstics> timeCompChars = null;
		ArrayList<DiscretizedFunc> timeSpeedups = null;
		ArrayList<DiscretizedFunc> combinedSpeedups = null;
		if (refFunc != null && 
				(plots.contains(time_comparison_title)
						|| plots.contains(time_speedup_vs_time_title)
						|| plots.contains(speedup_vs_threads_title))) {
			
			timeCompChars = new ArrayList<PlotCurveCharacterstics>();
			for (PlotCurveCharacterstics pltChar : avgChars)
				timeCompChars.add((PlotCurveCharacterstics)pltChar.clone());
			
			ArrayList<DiscretizedFunc> plotFuncs = averages;
			
			if (minFuncs.size() > 0) {
				plotFuncs = new ArrayList<DiscretizedFunc>();
				plotFuncs.addAll(averages);
				for (int i=0; i<averages.size(); i++) {
					DiscretizedFunc func = averages.get(i);
					if (func.getName().contains("dsa_")) {
						plotFuncs.add(minFuncs.get(i));
						plotFuncs.add(maxFuncs.get(i));
						
						PlotCurveCharacterstics minChar = (PlotCurveCharacterstics)avgChars.get(i).clone();
						minChar.setLineWidth(1f);
						minChar.setColor(getSaturated(minChar.getColor(), 0.5f));
						timeCompChars.add(minChar);
						
						PlotCurveCharacterstics maxChar = (PlotCurveCharacterstics)avgChars.get(i).clone();
						maxChar.setLineWidth(1f);
						maxChar.setColor(getSaturated(maxChar.getColor(), 0.5f));
						timeCompChars.add(maxChar);
					}
				}
			}
			
			System.out.println("generating time comparisons");
			timeComparisons = generateEnergyTimeComparison(plotFuncs, refFunc, avgNumX);
			ArrayList<DiscretizedFunc> combTimeComps = null;
			if (plots.contains(speedup_vs_threads_title)) {
				combTimeComps = new ArrayList<DiscretizedFunc>();
				for (DiscretizedFunc func : timeComparisons) {
					String name = func.getName().toLowerCase();
					if (name.contains("min") || name.contains("max"))
						combTimeComps.add(null);
					else
						combTimeComps.add(func);
				}
			}
			
			int numOrig = timeComparisons.size();
			if (do_extrap) {
				DiscretizedFunc extrapRef = getExtrapolatedRef(refFunc, getLowestEnergy(plotFuncs));
				if (extrapRef != null) {
					if (visible) {
						ArrayList<DiscretizedFunc> extrapFuncs = new ArrayList<DiscretizedFunc>();
						extrapFuncs.add(refFunc);
						extrapFuncs.add(extrapRef);
						ArrayList<PlotCurveCharacterstics> extrapChars = new ArrayList<PlotCurveCharacterstics>();
						extrapChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
						extrapChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
						getGraphWindow(extrapFuncs, "Extrapolated Curve", extrapChars, time_label, energy_label, visible);
					}
					timeComparisons.addAll(generateEnergyTimeComparison(plotFuncs, extrapRef, avgNumX/2));
					for (int i=numOrig-1; i>=0; i--) {
						DiscretizedFunc origFunc = timeComparisons.get(i);
						int extrapI = i+numOrig;
						if (timeComparisons.get(extrapI) == null) {
							timeComparisons.remove(extrapI);
						} else {
							DiscretizedFunc func = timeComparisons.get(extrapI);
							func.setName(func.getName()+" (EXTRAPOLATED)");
							PlotCurveCharacterstics cloned = (PlotCurveCharacterstics)timeCompChars.get(i).clone();
							cloned.setLineType(PlotLineType.DOTTED);
							cloned.setColor(getSaturated(cloned.getColor(), 0.5f));
							timeCompChars.add(numOrig, cloned);
							if (combTimeComps != null && combTimeComps.get(i) !=  null) {
								DiscretizedFunc stitched = new ArbitrarilyDiscretizedFunc();
								stitched.setName(origFunc.getName()+" (STICHED WITH EXTRAP)");
								for (int j=0; j<origFunc.size(); j++) {
									stitched.set(origFunc.get(j));
								}
								for (int j=0; j<func.size(); j++) {
									stitched.set(func.get(j));
								}
								combTimeComps.set(i, stitched);
							}
						}
					}
				}
			}
			
			if (plots.contains(time_speedup_vs_time_title) || plots.contains(speedup_vs_threads_title)) {
				System.out.println("generating "+time_speedup_vs_time_title);
				timeSpeedups = generateEnergyTimeSpeedup(timeComparisons);
			}
			if (combTimeComps != null) {
				for (int i=combTimeComps.size()-1; i>=0; i--)
					if (combTimeComps.get(i) == null)
						combTimeComps.remove(i);
				combinedSpeedups = generateEnergyTimeSpeedup(combTimeComps);
			}
		}
		
		HashMap<String, GraphWindow> windows = new HashMap<String, GraphWindow>();
		
		for (String plot : plots) {
			if (plot.equals(avg_energy_vs_time_title)) {
				if (averages.size() > 0) {
					System.out.println("displaying "+avg_energy_vs_time_title);
					windows.put(avg_energy_vs_time_title,
							getGraphWindow(averages, avg_energy_vs_time_title, avgChars, time_label, energy_label, visible));
				}
			} else if (plot.equals(std_dev_vs_time_title)) {
				if (stdDevs.size() > 0) {
					System.out.println("displaying "+std_dev_vs_time_title);
					windows.put(std_dev_vs_time_title,
							getGraphWindow(stdDevs, std_dev_vs_time_title, avgChars, time_label, std_dev_label, visible));
				}
			} else if (plot.equals(improvement_vs_time_title)) {
				System.out.println("generating percent improvements over "+pDiffMins+" mins");
				ArrayList<DiscretizedFunc> pImpFuncs = generatePercentImprovementOverTime(energyVsTime, pDiffMins);
				System.out.println("displaying "+improvement_vs_time_title);
				windows.put(improvement_vs_time_title,
						getGraphWindow(pImpFuncs, improvement_vs_time_title, chars,
						time_label, improvement_label, visible));
			} else if (plot.equals(improvement_vs_iterations)) {
				System.out.println("generating percent improvements over "+pDiffIters+" iters");
				ArrayList<DiscretizedFunc> pImpFuncs = generatePercentImprovementVsIterations(energyVsIter, pDiffIters);
				System.out.println("displaying "+improvement_vs_iterations);
				windows.put(improvement_vs_time_title,
						getGraphWindow(pImpFuncs, improvement_vs_iterations, chars,
						iterations_label, improvement_label, visible));
			} else if (plot.equals(time_speedup_vs_energy_title)) {
				if (refFunc != null) {
					System.out.println("generating energy speedup");
					ArrayList<DiscretizedFunc> energySpeedups = generateEnergySpeedup(averages, refFunc, avgNumX);
					System.out.println("displaying "+time_speedup_vs_energy_title);
					windows.put(time_speedup_vs_energy_title,
							getGraphWindow(energySpeedups, time_speedup_vs_energy_title, avgChars, energy_label, time_speedup_label, visible));
				}
			} else if (plot.equals(time_comparison_title)) {
				if (timeComparisons != null) {
					System.out.println("displaying "+time_comparison_title);
					windows.put(time_comparison_title,
							getGraphWindow(timeComparisons, time_comparison_title, timeCompChars,
									parallel_time_label, serial_time_label, visible));
				}
			} else if (plot.equals(time_speedup_vs_time_title)) {
				if (timeSpeedups != null) {
					System.out.println("displaying "+time_speedup_vs_time_title);
					windows.put(time_speedup_vs_time_title,
							getGraphWindow(timeSpeedups, time_speedup_vs_time_title, timeCompChars,
									parallel_time_label, time_speedup_label, visible));
				} else
					System.out.println("can't display "+time_speedup_vs_time_title+"...no funcs!");
			} else if (plot.equals(speedup_vs_threads_title)) {
				if (combinedSpeedups != null) {
					System.out.println("generating "+speedup_vs_threads_title);
					ArrayList<DiscretizedFunc> spdThreadFuncs = generateSpeedupVsThreads(dsaDir, combinedSpeedups);
					System.out.println("displaying "+speedup_vs_threads_title);
					ArrayList<PlotCurveCharacterstics> spdThreadChars = new ArrayList<PlotCurveCharacterstics>();
					spdThreadChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 8f, Color.BLACK));
					spdThreadChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
					spdThreadChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GREEN));
					windows.put(speedup_vs_threads_title,
							getGraphWindow(spdThreadFuncs, speedup_vs_threads_title, spdThreadChars, threads_label, time_speedup_label, visible));
				}
			} else if (plot.equals(energy_speedup_vs_time_title)) {
				if (refFunc != null) {
					System.out.println("generating energy speedup");
					ArrayList<DiscretizedFunc> energySpeedups = generateTimeSpeedup(averages, refFunc, avgNumX);
					System.out.println("displaying "+energy_speedup_vs_time_title);
					windows.put(energy_speedup_vs_time_title,
							getGraphWindow(energySpeedups, energy_speedup_vs_time_title, avgChars, time_label, energy_speedup_label, visible));
				}
			} else if (plot.equals(avg_energy_vs_iterations_title)) {
				if (energyIterAvgs.size() > 0) {
					System.out.println("displaying "+avg_energy_vs_iterations_title);
					windows.put(avg_energy_vs_iterations_title,
							getGraphWindow(energyIterAvgs, avg_energy_vs_iterations_title, avgChars, iterations_label, energy_label, visible));
				}
			} else if (plot.equals(energy_vs_iterations_title)) {
				System.out.println("displaying "+energy_vs_iterations_title);
				windows.put(energy_vs_iterations_title,
						getGraphWindow(energyVsIter, energy_vs_iterations_title, chars, iterations_label, energy_label, visible));
			} else if (plot.equals(energy_vs_time_title)) {
				System.out.println("displaying "+energy_vs_time_title);
				windows.put(energy_vs_time_title,
						getGraphWindow(energyVsTime, energy_vs_time_title, chars, time_label, energy_label, visible));
			} else if (plot.equals(iterations_vs_time_title)) {
				System.out.println("displaying "+iterations_vs_time_title);
				if (iterTimeAvgs.size() > 0)
					windows.put(iterations_vs_time_title,
							getGraphWindow(iterTimeAvgs, iterations_vs_time_title, avgChars, time_label, iterations_label, visible));
				else
					windows.put(iterations_vs_time_title,
							getGraphWindow(iterVsTime, iterations_vs_time_title, chars, time_label, iterations_label, visible));
			} else if (plot.equals(parallel_iterations_vs_time_title)) {
				System.out.println("displaying "+parallel_iterations_vs_time_title);
				ArrayList<DiscretizedFunc> regIterTimeCurves;
				ArrayList<PlotCurveCharacterstics> iterTimeChars;
				if (iterTimeAvgs.size() > 0) {
					regIterTimeCurves = iterTimeAvgs;
					iterTimeChars = avgChars;
				} else {
					regIterTimeCurves = iterVsTime;
					iterTimeChars = chars;
				}
				ArrayList<DiscretizedFunc> piters = getParallelIters(regIterTimeCurves, false);
				windows.put(parallel_iterations_vs_time_title,
							getGraphWindow(piters, parallel_iterations_vs_time_title, iterTimeChars, time_label, parallel_iterations_label, visible));
			} else if (plot.equals(energy_vs_parallel_iterations_title)) {
				System.out.println("displaying "+energy_vs_parallel_iterations_title);
				ArrayList<DiscretizedFunc> regEnergyIterCurves;
				ArrayList<PlotCurveCharacterstics> iterTimeChars;
				if (energyIterAvgs.size() > 0) {
					regEnergyIterCurves = energyIterAvgs;
					iterTimeChars = avgChars;
				} else {
					regEnergyIterCurves = energyVsIter;
					iterTimeChars = chars;
				}
				ArrayList<DiscretizedFunc> piters = getParallelIters(regEnergyIterCurves, true);
				windows.put(energy_vs_parallel_iterations_title,
							getGraphWindow(piters, energy_vs_parallel_iterations_title, iterTimeChars, parallel_iterations_label, energy_label, visible));
			} else {
				System.out.println("Unknown plot type: "+plot);
			}
		}
		
		return windows;
	}

}
