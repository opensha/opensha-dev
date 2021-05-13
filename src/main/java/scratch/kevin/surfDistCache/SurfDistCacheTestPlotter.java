package scratch.kevin.surfDistCache;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.io.Files;

public class SurfDistCacheTestPlotter {

	public static void main(String[] args) throws IOException {
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/direct_1");
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/hazard_1");
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/hazard_2_java7");
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/direct_2_java7");
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/direct_3_java7_no_exp");
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/hazard_3_java7_no_exp");/
//		File dir = new File("/home/kevin/OpenSHA/dist_cache/direct_5_hybrid_no_exp");
		File dir = new File("/home/kevin/OpenSHA/dist_cache/hazard_5_hybrid_no_exp");
		
		boolean memoryDebug = false;
		boolean memoryDebugTot = true;
		
		int maxCacheSize = 0;
		int totCount = 0;
		int maxThreads = 8;
		
		Map<Config, List<Double>> valsMap = Maps.newHashMap();
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (!name.contains("pbs.o"))
				continue;
			name = name.substring(0, name.indexOf(".pbs"));
			String[] split = name.split("_");
			int threads = Integer.parseInt(split[1]);
			if (threads > maxThreads)
				continue;
			int size = Integer.parseInt(split[3]);
			if (size > maxCacheSize)
				maxCacheSize = size;
			boolean forceMulti = Boolean.parseBoolean(split[5]);
			
			double vals = Double.NaN;
			
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				line = line.trim();
				if (memoryDebug) {
					if (line.startsWith("mem ")) {
						split = line.split("/");
						if (memoryDebugTot)
							vals = Double.parseDouble(split[split.length-3].substring(2)); // tot mem
						else
							vals = Double.parseDouble(split[split.length-2]); // used mem
					}
				} else {
					if (line.startsWith("Took ")) {
						vals = Double.parseDouble(line.split(" ")[1]); // seconds
					}
				}
			}
			if (Double.isNaN(vals)) {
				System.out.println("Job didn't finish: "+file.getName());
				continue;
			}
			
			Config conf = new Config(threads, size, forceMulti);
			List<Double> times = valsMap.get(conf);
			if (times == null) {
				times = Lists.newArrayList();
				valsMap.put(conf, times);
			}
			times.add(vals);
			
			totCount++;
		}
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (List<Double> times : valsMap.values())
			track.addValue(times.size());
		System.out.println("Loaded "+valsMap.size()+" configs (runs: "+track.toString()+")");
		Preconditions.checkState(!valsMap.isEmpty());
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(1d, 12d);
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		Table<Integer, Boolean, DiscretizedFunc> funcTable = HashBasedTable.create();
		
		List<Config> confs = Lists.newArrayList(valsMap.keySet());
		Collections.sort(confs);
		
		DiscretizedFunc refFunc = null;
		
		double maxX = 0;
		
		for (Config conf : confs) {
			List<Double> vals = valsMap.get(conf);
			double meanVal = 0;
			for (double val : vals)
				meanVal += val;
			meanVal /= vals.size();
			
			DiscretizedFunc func = funcTable.get(conf.size, conf.forceMulti);
			if (func == null) {
				func = new ArbitrarilyDiscretizedFunc();
				func.setName("Size="+conf.size+", Force="+conf.forceMulti);
				funcs.add(func);
				PlotLineType type;
				if (conf.forceMulti)
					type = PlotLineType.DASHED;
				else
					type = PlotLineType.SOLID;
				Color color;
				if (conf.size == 1)
					color = Color.BLACK;
				else if (conf.size == 0)
					color = Color.GRAY;
				else
					color = cpt.getColor((float)conf.size);
				if (!conf.forceMulti && conf.size == 1)
					refFunc = func;
				chars.add(new PlotCurveCharacterstics(type, 2f, PlotSymbol.FILLED_CIRCLE, 4f, color));
				funcTable.put(conf.size, conf.forceMulti, func);
			}
			func.set((double)conf.threads, meanVal);
			if (conf.threads > maxX)
				maxX = conf.threads;
		}
		
		String yAxisLabel, fname;
		if (memoryDebug) {
			yAxisLabel = "Used Memory (MB)";
			if (memoryDebugTot)
				fname = dir.getName()+"_totmem.png";
			else
				fname = dir.getName()+"_usedmem.png";
		} else {
			yAxisLabel = "Time (s)";
			fname = dir.getName()+"_results.png";
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Cache Speed Tests", "# Threads", yAxisLabel);
		spec.setLegendVisible(true);
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gw.setAxisRange(0.8, maxX + 0.2, 0, 2000);
		try {
			Thread.sleep(100);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		gw.saveAsPNG(new File(dir, fname).getAbsolutePath());
		
		if (refFunc != null) {
			// plot % improvement from reference
			List<DiscretizedFunc> impFuncs = Lists.newArrayList();
			List<PlotCurveCharacterstics> impChars = Lists.newArrayList();
			
			for (int i=0; i<funcs.size(); i++) {
				DiscretizedFunc func = funcs.get(i);
				PlotCurveCharacterstics plotChars = chars.get(i);
				
				if (func == refFunc)
					continue;
				
				ArbitrarilyDiscretizedFunc impFunc = new ArbitrarilyDiscretizedFunc();
				for (Point2D pt : func) {
					double x = pt.getX();
					double y = pt.getY();
					double refY;
					try {
						refY = refFunc.getY(x);
					} catch (Exception e) {
						continue;
					}
					impFunc.set(x, (refY - y)/refY);
				}
				if (impFunc.size() > 0) {
					impFunc.setName(func.getName());
					impFuncs.add(impFunc);
					impChars.add(plotChars);
				}
			}
			
			if (impFuncs.size() > 0) {
				spec = new PlotSpec(impFuncs, impChars, "Dist Cache Improvement", "# Threads", "Fract. Improvement");
				spec.setLegendVisible(true);
				
				new GraphWindow(spec);
			}
		}
		// now plot scaling of each one
		List<DiscretizedFunc> scaleFuncs = Lists.newArrayList();
		List<PlotCurveCharacterstics> scaleChars = Lists.newArrayList();
		
		for (int i=0; i<funcs.size(); i++) {
			DiscretizedFunc func = funcs.get(i);
			PlotCurveCharacterstics plotChars = chars.get(i);
			
			if (!func.hasX(1d) || func.size() < 2)
				continue;
			
			double serialTime = func.getY(1d);
			
			ArbitrarilyDiscretizedFunc speedup = new ArbitrarilyDiscretizedFunc();
			speedup.setName(func.getName());
			
			for (Point2D pt : func) {
				double x = pt.getX();
				double y = serialTime/pt.getY();
				speedup.set(x, y);
			}
			
			scaleFuncs.add(speedup);
			scaleChars.add(plotChars);
		}
		
		if (!scaleFuncs.isEmpty()) {
			// add ideal scaling
			ArbitrarilyDiscretizedFunc ideal = new ArbitrarilyDiscretizedFunc();
			ideal.setName("Ideal Scaling");
			ideal.set(1d, 1d);
			ideal.set(maxX, maxX);
			
			scaleFuncs.add(0, ideal);
			scaleChars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			spec = new PlotSpec(scaleFuncs, scaleChars, "Dist Cache Scaling", "# Threads", "Scaling");
			spec.setLegendVisible(true);
			
			new GraphWindow(spec);
		}
		
		// now plot relative scaling of each one
		scaleFuncs = Lists.newArrayList();
		scaleChars = Lists.newArrayList();

		for (int i=0; i<funcs.size(); i++) {
			DiscretizedFunc func = funcs.get(i);
			PlotCurveCharacterstics plotChars = chars.get(i);

			if (func.size() < 2)
				continue;

			ArbitrarilyDiscretizedFunc relScaling = new ArbitrarilyDiscretizedFunc();
			relScaling.setName(func.getName());
			
			for (int j=1; j<func.size(); j++) {
				double prevThreads = func.getX(j-1);
				double prevTime = func.getY(j-1);
				double threads = func.getX(j);
				double time = func.getY(j);
				
				double threadFactor = threads/prevThreads;
				
				double idealTime = prevTime/threadFactor;
				double eff = (prevTime - time)/(prevTime - idealTime);
				
				relScaling.set(threads, eff);
			}

			scaleFuncs.add(relScaling);
			scaleChars.add(plotChars);
		}

		if (!scaleFuncs.isEmpty()) {
			// add ideal scaling
			ArbitrarilyDiscretizedFunc ideal = new ArbitrarilyDiscretizedFunc();
			ideal.setName("Ideal Scaling");
			ideal.set(1d, 1d);
			ideal.set(maxX, 1d);

			scaleFuncs.add(0, ideal);
			scaleChars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));

			spec = new PlotSpec(scaleFuncs, scaleChars, "Dist Cache Scaling", "# Threads", "Relative Efficiency");
			spec.setLegendVisible(true);

			new GraphWindow(spec);
		}
	}
	
	private static class Config implements Comparable<Config> {
		private int threads, size;
		private boolean forceMulti;
		public Config(int threads, int size, boolean forceMulti) {
			super();
			this.threads = threads;
			this.size = size;
			this.forceMulti = forceMulti;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (forceMulti ? 1231 : 1237);
			result = prime * result + size;
			result = prime * result + threads;
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Config other = (Config) obj;
			if (forceMulti != other.forceMulti)
				return false;
			if (size != other.size)
				return false;
			if (threads != other.threads)
				return false;
			return true;
		}
		@Override
		public int compareTo(Config o) {
			int sizeComp = new Integer(size).compareTo(o.size);
			if (sizeComp != 0)
				return sizeComp;
			return new Boolean(forceMulti).compareTo(o.forceMulti);
		}
	}

}
