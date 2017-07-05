package scratch.kevin;

import java.awt.Color;

import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

public class GraphRedrawTest {
	
	public static void main(String[] args) throws InterruptedException {
		HistogramFunction hist = new HistogramFunction(0.05, 0.95, 10);
		
		GraphWindow gw = new GraphWindow(hist, "Random Number Hist",
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		
		while (true) {
			Thread.sleep(100);
			hist.add(Math.random(), 1d);
			gw.redrawGraph();
		}
	}

}
