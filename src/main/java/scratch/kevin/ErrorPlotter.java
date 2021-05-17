package scratch.kevin;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

import scratch.UCERF3.inversion.CommandLineInversionRunner;

import com.google.common.collect.Lists;

public class ErrorPlotter {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		List<double[]> dtVals = Lists.newArrayList();
		
		// dt for explicit
//		dtVals.add(toArray(22222, 0.1845, 0.015, 0.0789));
//		dtVals.add(toArray(20000, 0.1609, 0.0147, 0.079));
//		dtVals.add(toArray(1.67E+004, 0.1241, 0.0143, 0.0791));
//		dtVals.add(toArray(12500, 0.0792, 0.0138, 0.0792));
//		dtVals.add(toArray(8.33E+003, 0.0794, 0.0134, 0.0794));
//		dtVals.add(toArray(5000, 0.0797, 0.0133, 0.0797));
		
		// dt for implicit
		dtVals.add(toArray(22222, 0.1845, 0.015, 0.0789));
		dtVals.add(toArray(20000, 0.1609, 0.0147, 0.079));
		dtVals.add(toArray(1.67E+004, 0.1241, 0.0143, 0.0791));
		dtVals.add(toArray(12500, 0.0792, 0.0138, 0.0792));
		dtVals.add(toArray(8.33E+003, 0.0794, 0.0134, 0.0794));
		dtVals.add(toArray(5000, 0.0797, 0.0133, 0.0797));
		
		ArbitrarilyDiscretizedFunc maxFunc = new ArbitrarilyDiscretizedFunc("Max");
		ArbitrarilyDiscretizedFunc middleFunc = new ArbitrarilyDiscretizedFunc("Middle");
		ArbitrarilyDiscretizedFunc endFunc = new ArbitrarilyDiscretizedFunc("End");
		
		for (int i=0; i<dtVals.size(); i++) {
			double[] vals = dtVals.get(i);
			double x = vals[0];
			maxFunc.set(x, vals[1]);
			middleFunc.set(x, vals[2]);
			endFunc.set(x, vals[3]);
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(maxFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_SQUARE, 4f, Color.RED));
		funcs.add(middleFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_SQUARE, 4f, Color.GREEN));
		funcs.add(endFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_SQUARE, 4f, Color.BLUE));
		
		MinMaxAveTracker xTrack = new MinMaxAveTracker();
		MinMaxAveTracker yTrack = new MinMaxAveTracker();
		
		for (DiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				xTrack.addValue(pt.getX());
				yTrack.addValue(pt.getY());
			}
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setUserBounds(xTrack.getMin()*0.5, xTrack.getMax()*1.5, yTrack.getMin()*0.5, yTrack.getMax()*1.5);
		
		CommandLineInversionRunner.setFontSizes(gp);
		gp.setXLog(true);
		gp.setYLog(true);
		
		gp.drawGraphPanel("DT (years)", "Relative Error",
				funcs, chars, "Error Vs DT");
		
		gp.getChartPanel().setSize(500, 500);
		gp.saveAsPNG("/home/kevin/Documents/Geol 557/prob_set_06/explicit_error_dt.png");
	}
	
	public static double[] toArray(double... vals) {
		return vals;
	}

}
