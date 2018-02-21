package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;

import com.google.common.collect.Lists;

public class GroundMotionScatterPlot {
	
	public static int PLOT_WIDTH = 900;
	public static boolean WRITE_PDF = true;
	public static String SCATTER_QUANTITY_NAME = "Ruptures";
	public static boolean YELLOW_REGION = true;
	public static int MAX_SCATTER_POINTS = 100000;
	
	public static void plot(XY_DataSet xy, String xAxisLabel, String yAxisLabel, List<String> binDescriptions,
			String title, File outputDir, String prefix) throws IOException {
		Range range = new Range(Math.min(xy.getMinX(), xy.getMinY()),
				Math.max(xy.getMaxX(), xy.getMaxY()));
		
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		// shaded yellow area
		ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
		oneToOne.set(xy.getMinX(), xy.getMinX());
		oneToOne.set(xy.getMaxX(), xy.getMaxX());
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		oneToOne.set(range.getLowerBound()*2, range.getLowerBound()*2);
		oneToOne.set(range.getLowerBound()/2, range.getLowerBound()/2);
		oneToOne.set(range.getUpperBound()*2, range.getUpperBound()*2);
		oneToOne.set(range.getUpperBound()/2, range.getUpperBound()/2);
		ArbitrarilyDiscretizedFunc upper = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : oneToOne)
			upper.set(pt.getX(), pt.getY()*2);
		ArbitrarilyDiscretizedFunc lower = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : oneToOne)
			lower.set(pt.getX(), pt.getY()/2);
		UncertainArbDiscDataset shaded = new UncertainArbDiscDataset(oneToOne, lower, upper);
		
		if (YELLOW_REGION) {
			funcs.add(shaded);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.YELLOW));
			
			funcs.add(upper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
			
			funcs.add(lower);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		}
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		
		XY_DataSet origXY = xy;
		if (xy.size() > MAX_SCATTER_POINTS) {
			System.out.println("Filtering scatter points from "+xy.size()+" to ~"+MAX_SCATTER_POINTS);
			// use fixed seed for reproducibility of downsampled plots
			Random r = new Random(xy.size());
			double rand = (double)MAX_SCATTER_POINTS/(double)xy.size();
			XY_DataSet filtered = new DefaultXY_DataSet();
			for (Point2D pt : xy)
				if (r.nextDouble() < rand)
					filtered.set(pt);
			System.out.println("\tNew size: "+filtered.size());
			xy = filtered;
		}
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.RED));
		
		// regression in log space
		SimpleRegression regression = new SimpleRegression();
		for (Point2D pt : origXY)
			regression.addData(Math.log(pt.getX()), Math.log(pt.getY()));
		double b = regression.getIntercept();
		double m = regression.getSlope();
		DefaultXY_DataSet fit = new DefaultXY_DataSet();
		// use one to one for x values
		for (int i=0; i<oneToOne.size(); i++) {
			double origX = oneToOne.getX(i);
			if (origX < xy.getMinX() || origX > xy.getMaxX())
				continue;
			double x = Math.log(origX);
			double y = m*x + b;
			fit.set(Math.exp(x), Math.exp(y));
		}
		funcs.add(fit);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		
		// bottom up
		List<String> annLines = new ArrayList<>(binDescriptions);
		Collections.reverse(annLines);
		annLines.add(0, xy.size()+" "+SCATTER_QUANTITY_NAME);
//		annLines.add(str(distRange.getLowerBound())+" km < Dist < "+str(distRange.getUpperBound())+" km");
//		annLines.add(str(magRange.getLowerBound())+" < Mw < "+str(magRange.getUpperBound()));
		
		// this lines up annotations in log space
		EvenlyDiscretizedFunc lineSpacings = new EvenlyDiscretizedFunc(
				Math.log(range.getLowerBound()), Math.log(range.getUpperBound()), 17);
		
		List<XYTextAnnotation> anns = Lists.newArrayList();
		for (int i=0; i<annLines.size(); i++) {
			double x = Math.exp(lineSpacings.getX(lineSpacings.size()-2));
			double y = Math.exp(lineSpacings.getX(1+i));
			String annText = annLines.get(i);
			XYTextAnnotation ann = new XYTextAnnotation(annText, x, y);
			ann.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
			if (PLOT_WIDTH >= 800)
				ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));
			else
				ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
			anns.add(ann);
		}
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, true, true, range, range);
		gp.getChartPanel().setSize(PLOT_WIDTH, PLOT_WIDTH-100);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		if (WRITE_PDF)
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

}
