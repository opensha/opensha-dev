package scratch.kevin;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;

public class GPLegendTest {

	public static void main(String[] args) throws IOException {
		PlotLineType[] types = {
				PlotLineType.HISTOGRAM,
				null,
				PlotLineType.SOLID,
				PlotLineType.DASHED,
				PlotLineType.DOTTED,
				PlotLineType.DOTTED_AND_DASHED
		};
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<types.length; i++) {
			DiscretizedFunc xy = new EvenlyDiscretizedFunc(0d, 1d, 10);
			for (int j=0; j<xy.size(); j++)
				xy.set(j, i+1d);
			funcs.add(xy);
			if (types[i] == null) {
				xy.setName("Circles");
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
//			} else if (types[i] == PlotLineType.HISTOGRAM) {
//				xy.setName(types[i].name());
//				chars.add(new PlotCurveCharacterstics(types[i], 2f, Color.BLACK));
			} else {
				xy.setName(types[i].name());
				chars.add(new PlotCurveCharacterstics(types[i], 2f, Color.BLACK));
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Test", "x", "y");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, new Range(0d, 1d), new Range(0d, types.length+1));

		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG("/tmp/legend.png");
	}

}
