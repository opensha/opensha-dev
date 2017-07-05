package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.swing.JFrame;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import scratch.UCERF3.inversion.CommandLineInversionRunner;

import com.google.common.collect.Lists;
public class LegendTests {
	
	public static void main(String[] args) throws IOException {
		List<PlotSpec> specs = Lists.newArrayList();
		int numSpecs = 1;
		boolean hist = true;
		
		for (int s=0; s<numSpecs; s++) {
			DiscretizedFunc data1, data2;
			
			if (hist) {
				data1 = new HistogramFunction(0d, 1d, 10);
				data1.setName("Histogram 1 P"+s);
				for (int i=0; i<10; i++)
					((HistogramFunction)data1).add(Math.random(), Math.random());
				
				data2 = new HistogramFunction(0d, 1d, 10);
				data2.setName("Histogram 2 P"+s);
				for (int i=0; i<10; i++)
					((HistogramFunction)data2).add(Math.random(), Math.random());
			} else {
				data1 = new ArbitrarilyDiscretizedFunc();
				data1.setName("Super Long Title one P"+s);
				data1.set(Math.random(), Math.random());
				data1.set(Math.random(), Math.random());
				data1.set(Math.random(), Math.random());
				data1.set(Math.random(), Math.random());
				data1.set(Math.random(), Math.random());
				data1.set(Math.random(), Math.random());
				
				data2 = new ArbitrarilyDiscretizedFunc();
				data2.setName("Asdf2 P"+s);
				data2.set(Math.random(), Math.random());
				data2.set(Math.random(), Math.random());
				data2.set(Math.random(), Math.random());
				data2.set(Math.random(), Math.random());
				data2.set(Math.random(), Math.random());
				data2.set(Math.random(), Math.random());
			}
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			funcs.add(data1);
			funcs.add(data2);
			List<PlotCurveCharacterstics> chars;
			if (hist) {
				chars = Lists.newArrayList(
						new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED),
						new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
			} else {
				chars = Lists.newArrayList(
						new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED),
						new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			}
			PlotSpec spec = new PlotSpec(funcs, chars, "Title P"+s, "X", "Y");
			if (s == 0)
				spec.setLegendVisible(true);
			specs.add(spec);
		}
//		GraphWindow gw = new GraphWindow(spec);
//		gw.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		
//		List<Range> xRanges = Lists.newArrayList();
//		List<Range> yRanges = Lists.newArrayList();
		gp.drawGraphPanel(specs, false, false, null, null);

		File file = new File("/tmp/legend_test.png");
		gp.getChartPanel().setSize(1000, 1000);
		gp.saveAsPNG(file.getAbsolutePath());
	}

}
