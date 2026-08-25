package scratch.kevin.mfdInversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleFunction;

import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator;

import net.mahdilamb.colormap.Colors;

public class LikelihoodDistPlot {

	public static void main(String[] args) throws IOException {
		Range xRange = new Range(0d, 5d);
		
		List<DoubleFunction<Double>> functions = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		functions.add(PaleoBValueEstimator.GAUSSIAN_LIKELIHOOD);
		names.add("Gaussian");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		functions.add(PaleoBValueEstimator.getStudentsTLogLikelihood(4, false));
		names.add("Student's t, nu=4 (orig)");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, Colors.tab_orange));
		
		functions.add(PaleoBValueEstimator.getStudentsTLogLikelihood(4, true));
		names.add("Student's t, nu=4 (scaled)");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
		
		functions.add(PaleoBValueEstimator.getStudentsTLogLikelihood(5, false));
		names.add("Student's t, nu=5 (orig)");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, Colors.tab_blue));
		
		functions.add(PaleoBValueEstimator.getStudentsTLogLikelihood(5, true));
		names.add("Student's t, nu=5 (scaled)");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
		
		List<EvenlyDiscretizedFunc> funcs = new ArrayList<>();
		for (int i=0; i<functions.size(); i++) {
			DoubleFunction<Double> f = functions.get(i);
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), 100);
			func.setName(names.get(i));
			for (int x=0; x<func.size(); x++)
				func.set(x, f.apply(func.getX(x)));
			funcs.add(func);
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, "Log-likelihood Functions", "|z-score|", "Log-likelihood");
		plot.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		
		gp.drawGraphPanel(plot, false, false, xRange, null);
		
		PlotUtils.writePlots(new File("/tmp"), "likelihoods", gp, 800, 800, true, true, false);
	}

}
