package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.block.ColumnArrangement;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

public class EmbeddedXYLegendTest {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "legend_test";
		
		Range xRange = new Range(0.1d, 100d);
		Range yRange = new Range(0.1d, 100d);
		
		DiscretizedFunc curve1 = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), 20);
		DiscretizedFunc curve2 = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), 20);
		
		for (int i=0; i<curve1.size(); i++) {
			curve1.set(i, curve1.getX(i) + 0.25*Math.random());
			curve2.set(i, curve2.getX(i) + 0.4*Math.random());
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		curve1.setName("Curve 1");
		funcs.add(curve1);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		curve2.setName("Curve 2");
		funcs.add(curve2);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Legend Test", "X", "Y");
//		spec.setLegendInset(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
//		double x = xRange.getCentralValue();
//		double y = yRange.getCentralValue();
		double x = 0.95;
		double y = 0.95;
		
		System.out.println("Plotting legend at "+x+", "+y);
		XYPlot plot = gp.getPlot();
		LegendTitle lt = new LegendTitle(plot, new ColumnArrangement(), new ColumnArrangement());
//		LegendTitle lt = new OnePerLineLegendTitle(plot);
		lt.setItemFont(new Font(Font.SANS_SERIF, Font.PLAIN, gp.getPlotPrefs().getLegendFontSize()));
		lt.setBackgroundPaint(new Color(255, 255, 255, 200));
		lt.setFrame(new BlockBorder(Color.BLACK));
		lt.setPosition(RectangleEdge.TOP);
//		lt.setHorizontalAlignment(HorizontalAlignment.RIGHT);
//		lt.set
//		lt.setVerticalAlignment(VerticalAlignment.TOP);
		XYTitleAnnotation ann = new XYTitleAnnotation(x, y, lt, RectangleAnchor.TOP_RIGHT);
		ann.setMaxWidth(0.7d);
//		ann.set
		plot.addAnnotation(ann);
		gp.getChartPanel().getChart().fireChartChanged();
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
}
