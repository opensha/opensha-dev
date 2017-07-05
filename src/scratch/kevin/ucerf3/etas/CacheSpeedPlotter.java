package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

public class CacheSpeedPlotter {
	
	private static String timeStr = "Looping over events took ";

	public static void main(String[] args) throws IOException {
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/cache_tests/hard_values");
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/cache_tests/soft_values");
		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/cache_tests/soft_10Gjvm");
		
		Map<Double, List<Double>> timesMap = Maps.newHashMap();
		
		for (File subDir : dir.listFiles()) {
			if (!subDir.isDirectory())
				continue;
			String name = subDir.getName();
			if (!name.startsWith("cache_"))
				continue;
			
			name = name.substring(6);
			name = name.substring(0, name.indexOf("gb"));
			double size = Double.parseDouble(name);
			
//			System.out.println("Size: "+size);
			
			for (File runDir : subDir.listFiles()) {
				String runName = runDir.getName();
				if (!runDir.isDirectory() || !runName.startsWith("run_"))
					continue;
				File infoFile = new File(runDir, "infoString.txt");
				if (!infoFile.exists())
					continue;
				for (String line : Files.readLines(infoFile, Charset.defaultCharset())) {
					if (line.contains(timeStr)) {
						line = line.substring(timeStr.length());
						line = line.substring(0, line.indexOf("secs"));
						double time = Double.parseDouble(line);
						
						System.out.println("Size: "+(float)size+"\tTime: "+(float)time);
						
						List<Double> times = timesMap.get(size);
						if (times == null) {
							times = Lists.newArrayList();
							timesMap.put(size, times);
						}
						times.add(time);
					}
				}
			}
		}
		
		ArbitrarilyDiscretizedFunc minFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc maxFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
		
		for (double size : timesMap.keySet()) {
			List<Double> times = timesMap.get(size);
			
			double[] timesArray = Doubles.toArray(times);
			
			minFunc.set(size, StatUtils.min(timesArray));
			maxFunc.set(size, StatUtils.max(timesArray));
			meanFunc.set(size, StatUtils.mean(timesArray));
		}
		
		UncertainArbDiscDataset uncertFunc = new UncertainArbDiscDataset(meanFunc, minFunc, maxFunc);
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(uncertFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.BLUE));
		
		funcs.add(meanFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Cache Speed Tests", "Max Cache Memory (GB)", "Time (sec)");
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		gw.setY_AxisRange(0, 1800);
	}

}
