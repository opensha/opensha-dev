package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;

import java.util.List;

public class InsetLegendTests {

	public static void main(String[] args) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		double relX = 0.975;
		double relY = 0.975;
		
		Range xRange = new Range(1e-8, 1e4);
		Range yRange = new Range(1e-6, 0.5);
		
		boolean xLog = true;
		boolean yLog = true;
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), 5);
		for (int i=0; i<func.size(); i++)
			func.set(i, Math.random()*yRange.getLength() + yRange.getLowerBound());
		
		func.setName("This is my function");
		funcs.add(func);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3, Color.BLUE));
		
		ArbitrarilyDiscretizedFunc pt = new ArbitrarilyDiscretizedFunc();
		double actualX = getActual(relX, xRange, xLog);
		double actualY = getActual(relY, yRange, yLog);
		pt.set(actualX, actualY);
		pt.setName("Location");
		funcs.add(pt);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 8f, Color.GREEN));
		System.out.println("Actual: "+actualX+"\t"+actualY);
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Title", "X", "Y");
		spec.setLegendInset(RectangleAnchor.TOP_RIGHT, relX, relY, 0.35, true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, xLog, yLog, xRange, yRange);
		
		PlotUtils.writePlots(new File("/tmp"), "inset_test", gp, 1000, 800, true, false, false);
	}
	
	private static double getActual(double rel, Range range, boolean log) {
		if (log)
			return Math.pow(10, rel * (Math.log10(range.getUpperBound())
					- Math.log10(range.getLowerBound())) + Math.log10(range.getLowerBound()));
		return rel*range.getLength() + range.getLowerBound();
	}

}
