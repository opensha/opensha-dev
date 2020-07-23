package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

public class Stampede2CoresPerNodePlot {

	public static void main(String[] args) throws IOException {
		File dataDir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2020_07_13-stampede2-benchmark/results");
		Map<Integer, List<DiscretizedFunc>> funcMap = new HashMap<>();
		
		for (File file : dataDir.listFiles()) {
			String name = file.getName();
			if (!name.endsWith(".csv") || !name.contains("_thread"))
				continue;
			String threadStr = name.substring(name.lastIndexOf("_")+1);
			threadStr = threadStr.substring(0, threadStr.indexOf(".csv"));
			int threads = Integer.parseInt(threadStr);
			System.out.println("Loading "+threads+" threads from "+name);
			CSVFile<String> csv = CSVFile.readFile(file, true);
			
			ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
			for (int row=1; row<csv.getNumRows(); row++) {
				long millis = csv.getLong(row, 1);
				double secs = millis/1000d;
				double mins = secs/60d;
				double energy = csv.getDouble(row, 2);
				func.set(mins, energy);
			}
			System.out.println("\tmin energy: "+func.getMinY());
			List<DiscretizedFunc> funcs = funcMap.get(threads);
			if (funcs == null) {
				funcs = new ArrayList<>();
				funcMap.put(threads, funcs);
			}
			funcs.add(func);
		}
		
		List<Integer> threads = new ArrayList<>(funcMap.keySet());
		Collections.sort(threads);
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(1, threads.size()-1);
		List<Color> colors = new ArrayList<>();
		for (int i=0; i<threads.size(); i++) {
			if (i == 0)
				colors.add(Color.BLACK);
			else
				colors.add(cpt.getColor((float)i).darker());
		}
		
		List<DiscretizedFunc> meanFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		double largestMinEnergy = 0d;
		double smallestMinEnergy = Double.POSITIVE_INFINITY;
		for (int i=0; i<threads.size(); i++) {
			List<DiscretizedFunc> tFuncs = funcMap.get(threads.get(i));
			double largestMin = 0d;
			double smallestMax = Double.POSITIVE_INFINITY;
			for (DiscretizedFunc func : tFuncs) {
				largestMin = Math.max(largestMin, func.getMinX());
				smallestMax = Math.min(smallestMax, func.getMaxX());
			}
			EvenlyDiscretizedFunc meanFunc = new EvenlyDiscretizedFunc(largestMin, smallestMax, 1000);
			for (int j=0; j<meanFunc.size(); j++) {
				double avgVal = 0d;
				double x = meanFunc.getX(j);
				for (DiscretizedFunc func : tFuncs)
					avgVal += func.getInterpolatedY(x);
				avgVal /= (double)tFuncs.size();
				meanFunc.set(j, avgVal);
			}
			largestMinEnergy = Math.max(largestMinEnergy, meanFunc.getMinY());
			smallestMinEnergy = Math.min(smallestMinEnergy, meanFunc.getMinY());
			
			meanFunc.setName(threads.get(i)+" Threads");
			meanFuncs.add(meanFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(i)));
		}
		
		PlotSpec spec = new PlotSpec(meanFuncs, chars, "Convergence Thread Test",
				"Time (minutes)", "Energy");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(16);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		Range xRange = null;
		Range yRange = new Range(00.9*smallestMinEnergy, 4*largestMinEnergy);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File file = new File(dataDir, "energy_vs_time");
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		// now time to 150 plot
		double[] targets = { 150d, 130d };
		
		for (double target : targets) {
			DiscretizedFunc threadTimeFunc = new ArbitrarilyDiscretizedFunc();
			double oneThreadTime = Double.NaN;
			for (int i=0; i<meanFuncs.size(); i++) {
				DiscretizedFunc meanFunc = meanFuncs.get(i);
				if (meanFunc.getMinY() > target)
					continue;
				double time = meanFunc.getFirstInterpolatedX(target);
				int myThreads = threads.get(i);
				if (myThreads == 1)
					oneThreadTime = time;
				threadTimeFunc.set((double)myThreads, time);
			}
			EvenlyDiscretizedFunc idealScaling = new EvenlyDiscretizedFunc(
					threadTimeFunc.getMinX(), threadTimeFunc.getMaxX(), 500);
			EvenlyDiscretizedFunc sqrtScaling = new EvenlyDiscretizedFunc(
					threadTimeFunc.getMinX(), threadTimeFunc.getMaxX(), 500);
			for (int i=0; i<idealScaling.size(); i++) {
				double myThreads = idealScaling.getX(i);
				double ideal = oneThreadTime/myThreads;
				idealScaling.set(i, ideal);
				double sqrt = oneThreadTime/Math.sqrt(myThreads);
				sqrtScaling.set(i, sqrt);
			}
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			threadTimeFunc.setName("Actual Scaling");
			funcs.add(threadTimeFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
			
			idealScaling.setName("Ideal Scaling");
			funcs.add(idealScaling);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			
			sqrtScaling.setName("Sqrt(N) Scaling");
			funcs.add(sqrtScaling);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			spec = new PlotSpec(funcs, chars, "Scaling Test",
					"Threads", "Time to E="+(int)target+" (m)");
			spec.setLegendVisible(true);
			
			xRange = null;
			yRange = null;

			gp.setLegendFontSize(22);
			gp.drawGraphPanel(spec, true, false, xRange, yRange);
			
			file = new File(dataDir, "scaling_to_"+(int)target+"_log");
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			file = new File(dataDir, "scaling_to_"+(int)target);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			
			// strong scaling
			
			DiscretizedFunc strongScalingFunc = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<meanFuncs.size(); i++) {
				int myThreads = threads.get(i);
				double myTime = threadTimeFunc.getY((double)myThreads);
				double speedup = oneThreadTime/myTime;
				strongScalingFunc.set((double)myThreads, speedup);
			}
			for (int i=0; i<idealScaling.size(); i++) {
				double myThreads = idealScaling.getX(i);
				double ideal = myThreads;
				idealScaling.set(i, ideal);
				double sqrt = Math.sqrt(myThreads);
				sqrtScaling.set(i, sqrt);
			}
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			threadTimeFunc.setName("Actual Scaling");
			funcs.add(strongScalingFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
			
			idealScaling.setName("Ideal Scaling");
			funcs.add(idealScaling);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			
			sqrtScaling.setName("Sqrt(N) Scaling");
			funcs.add(sqrtScaling);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			spec = new PlotSpec(funcs, chars, "Strong Scaling",
					"Threads", "Speedup");
			spec.setLegendVisible(true);
			
			double maxThreads = idealScaling.getMaxX();
			xRange = new Range(0d, maxThreads);
			yRange = xRange;
			
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			file = new File(dataDir, "strong_scaling_"+(int)target);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			
			// scaling fraction
			
			DiscretizedFunc strongScalingFract = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<meanFuncs.size(); i++) {
				strongScalingFract.set(strongScalingFunc.getX(i), strongScalingFunc.getY(i)/strongScalingFunc.getX(i));
			}
			for (int i=0; i<idealScaling.size(); i++) {
				double myThreads = idealScaling.getX(i);
				double sqrt = Math.sqrt(myThreads);
				sqrtScaling.set(i, sqrt/myThreads);
			}
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			strongScalingFract.setName("Actual Scaling");
			funcs.add(strongScalingFract);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
			
			sqrtScaling.setName("Sqrt(N) Scaling");
			funcs.add(sqrtScaling);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			spec = new PlotSpec(funcs, chars, "Strong Scaling",
					"Threads", "Scaling Efficiency");
			spec.setLegendVisible(true);
			
			xRange = new Range(0d, maxThreads);
			yRange = new Range(0d, 1d);
			
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			file = new File(dataDir, "strong_scaling_fract_"+(int)target);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		}
	}

}
