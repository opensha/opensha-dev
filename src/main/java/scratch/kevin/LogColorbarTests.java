package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

public class LogColorbarTests {

	public static void main(String[] args) throws IOException {
		logXYdemo(new File("C:\\Users\\Kevin Milner\\Downloads"), "log_xy");
		logXYZdemo(new File("C:\\Users\\Kevin Milner\\Downloads"), "log_xyz");
	}
	
	public static void logXYdemo(File outputDir, String outputPrefix) throws IOException {
		DiscretizedFunc data = new ArbitrarilyDiscretizedFunc();
		for (double x=0d; x<2.001; x+=0.01)
			data.set(x, Math.pow(10, x));
		
		PlotSpec plot = new PlotSpec(
				List.of(data),
				List.of(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK)),
				"XY Plot Demo", "X", "Y");
		
		Range xRange = new Range(data.getMinX(), data.getMaxX());
		Range yRange = new Range(1d, data.getMaxY());
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		// create a multi plot with the same data, but the top one in
		// linear-linear and the bottom linear-log
		gp.drawGraphPanel(List.of(plot, plot),
				List.of(false), List.of(false, true),
				List.of(xRange), List.of(yRange, yRange));
		
		PlotUtils.writePlots(outputDir, outputPrefix, gp, 800, 800, true, false, false);
	}
	
	public static void logXYZdemo(File outputDir, String outputPrefix) throws IOException {
		EvenlyDiscrXYZ_DataSet data = new EvenlyDiscrXYZ_DataSet(100, 100, 1d, 1d, 0.1);
		for (int xInd=0; xInd<data.getNumX(); xInd++) {
			double x = data.getX(xInd);
			for (int yInd=0; yInd<data.getNumY(); yInd++) {
				double y = data.getY(yInd);
				data.set(xInd, yInd, x*y);
			}
		}
		// convert the data to log10
		data.log10();
		CPT colorscale = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(
				data.getMinZ(), data.getMaxZ());
		
		XYZPlotSpec plot = new XYZPlotSpec(data, colorscale, "XYZ Plot Demo", "X", "Y", "Log10 Z");
		plot.setCPTPosition(RectangleEdge.BOTTOM);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();

		gp.drawGraphPanel(plot, false, false, new Range(data.getMinX(), data.getMaxX()),
				new Range(data.getMinY(), data.getMaxY()));
		PlotUtils.setXTick(gp, 1d);
		PlotUtils.setYTick(gp, 1d);
		
		PlotUtils.writePlots(outputDir, outputPrefix, gp, 800, 800, true, false, false);
		
		// now write out a plot with the x-axis on a log scale to make the mockup of what we wanted
		// that to actually look like
		
		DiscretizedFunc emptyFunc = new LightFixedXFunc(new double[] {0d, 0d}, new double[] {0d, 0d});
		gp.drawGraphPanel(new PlotSpec(List.of(emptyFunc),
				List.of(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, Color.WHITE)),
				" ", "Z", null), true, false,
				new Range(Math.pow(10, data.getMinZ()), Math.pow(10, data.getMaxZ())),
				new Range(-0.5d, 1d));
		PlotUtils.setYTick(gp, 1d);
		
		PlotUtils.writePlots(outputDir, outputPrefix+"_for_mockup", gp, 800, 800, true, false, false);
	}

}
