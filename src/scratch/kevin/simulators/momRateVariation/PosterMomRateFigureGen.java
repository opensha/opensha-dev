package scratch.kevin.simulators.momRateVariation;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class PosterMomRateFigureGen {

	public static void main(String[] args) throws IOException {
		int windowLen = 25;
		int length = 3000;
//		int length = 900;
		
		File inputDir = new File("/home/kevin/Simulators/time_series");
		File outputDir = new File(inputDir, "ts_figures");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
//		String[] prefixes = {"actual", "random", "ucerf3_mid"};
//		boolean[] yAxisLabels = {false, true, false};
//		String[] names = {"RSQSim", "RSQSim Poisson", "UCERF3-TD"};
//		Color[] colors = {Color.RED, Color.BLUE, new Color(64, 128, 64)};
		String[] prefixes = {"actual", "random", "u3sighi", "ucerf3_mid", "ucerf3_etas"};
		boolean[] yAxisLabels = {false, false, true, false, false};
		String[] names = {"RSQSim", "RSQSim Poisson", "RSQSim-U3 (Prelim.)", "UCERF3-TD", "UCERF3-ETAS"};
		Color[] colors = {new Color(64, 128, 64), Color.BLUE, Color.RED, Color.CYAN.darker(), Color.MAGENTA.darker()};
		
		Range xRange = new Range(0, length);
//		Range yRange = new Range(1e18, 1e20);
		Range yRange = new Range(5e17, 1e20);
		
		List<double[]> series = Lists.newArrayList();
		
		for (int i=0; i<prefixes.length; i++) {
			String prefix = prefixes[i];
			File file = new File(inputDir, prefix+"_"+windowLen+"yr.bin");
			Preconditions.checkState(file.exists(), "File doesn't exist: %s", file.getName());
			
			double[] mySeries = MatrixIO.doubleArrayFromFile(file);
			System.out.println(file.getName()+" length: "+mySeries.length);
			Preconditions.checkState(mySeries.length >= length);
			
			int startIndex = 0;
			if (!prefix.startsWith("ucerf3"))
				// simulator
				startIndex = 10000;
			
			mySeries = Arrays.copyOfRange(mySeries, startIndex, startIndex+length);
			
			series.add(mySeries);
		}
		
		List<PlotSpec> specs = Lists.newArrayList();
		
		String title = "Tapered Seismic Moment Release";
		
		for (int i=0; i<series.size(); i++) {
			double[] vals = series.get(i);
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0, vals.length, 1d);
			for (int n=0; n<vals.length; n++)
				func.set(n, vals[n]);
			
			List<EvenlyDiscretizedFunc> funcs = Lists.newArrayList();
			funcs.add(func);
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors[i]));
			
			String yAxisLabel = "Moment Rate (N-m/yr)";
			if (!yAxisLabels[i])
				yAxisLabel = " ";
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Years", yAxisLabel);
			
			double x = xRange.getLowerBound()+(xRange.getUpperBound()-xRange.getLowerBound())*0.1;
			double y = yRange.getLowerBound()+(yRange.getUpperBound()-yRange.getLowerBound())*0.9;
			String label = names[i];
			XYTextAnnotation ann = new XYTextAnnotation(label, x, y);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
			ann.setTextAnchor(TextAnchor.TOP_LEFT);
			List<XYTextAnnotation> anns = Lists.newArrayList(ann);
			spec.setPlotAnnotations(anns);
			
			specs.add(spec);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		List<Range> xRanges = Lists.newArrayList(xRange);
		List<Range> yRanges = Lists.newArrayList();
		for (int i=0; i<specs.size(); i++)
			yRanges.add(yRange);
		gp.drawGraphPanel(specs, false, true, xRanges, yRanges);
		gp.getChartPanel().setSize(1200, 900);
		gp.saveAsPNG(new File(outputDir, "combined_time_series.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "combined_time_series.pdf").getAbsolutePath());
	}

}
