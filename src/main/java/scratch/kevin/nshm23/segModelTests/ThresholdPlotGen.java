package scratch.kevin.nshm23.segModelTests;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;

public class ThresholdPlotGen {

	public static void main(String[] args) throws IOException {
		Shaw07JumpDistProb calc = new Shaw07JumpDistProb(1, 3);
		
		double maxJumpDist = 15d;
		
		double maxY = 1.2;
		double maxX = 20;
		double minY = calc.calcJumpProbability(maxJumpDist);
		
		double[] binEdges = { 0d, 1d, 3d, 5d, 7d, 9d, 11d, 13d };
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		// add shaw line
		EvenlyDiscretizedFunc modelLine = new EvenlyDiscretizedFunc(0d, maxJumpDist, 1000);
		for (int i=0; i<modelLine.size(); i++)
			modelLine.set(i, calc.calcJumpProbability(modelLine.getX(i)));
		
		modelLine.setName(calc.getName());
		funcs.add(modelLine);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED.darker()));
		
		XY_DataSet binPoints = new DefaultXY_DataSet();
		binPoints.setName("Probability Levels");
		
		double[] binProbs = new double[binEdges.length];
		for (int i=0; i<binEdges.length; i++) {
			double x = binEdges[i];
			double y = calc.calcJumpProbability(x);
			binPoints.set(x, y);
			binProbs[i] = y;
			
			DefaultXY_DataSet line = new DefaultXY_DataSet();
			line.set(x, y);
			line.set(maxX, y);
			funcs.add(line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
		}
		
		funcs.add(binPoints);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 5f, Color.BLACK));
		
		List<XYTextAnnotation> anns = new ArrayList<>();
		
		for (int i=0; i<binProbs.length; i++) {
			double x = binEdges[i];
			double y = binProbs[i];
			
			XYTextAnnotation ptAnn = new XYTextAnnotation(" P"+(i+1), x, y);
			ptAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 15));
			ptAnn.setTextAnchor(TextAnchor.BOTTOM_LEFT);
			anns.add(ptAnn);
			
			double plotY2, weight;
			String weightStr = "W"+(i+1)+" = P"+(i+1)+" - ";
			if (i == binProbs.length-1) {
				plotY2 = minY;
				weightStr += "0";
				weight = y;
			} else {
				plotY2 = binProbs[i+1];
				weightStr += "P"+(i+2);
				weight = y-binProbs[i+1];
			}
			double plotY = Math.pow(10, 0.5*(Math.log10(plotY2)+Math.log10(y)));
			weightStr += " = "+wtDF.format(weight)+" ";
			XYTextAnnotation weightAnn = new XYTextAnnotation(weightStr, maxX, plotY);
			weightAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 15));
			weightAnn.setTextAnchor(TextAnchor.CENTER_RIGHT);
			anns.add(weightAnn);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Threshold Averaging", "Distance (km)", "Probability");
		spec.setLegendVisible(true);
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, true, new Range(-0.2, maxX), new Range(minY, maxY));
		
		PlotUtils.writePlots(new File("/tmp"), "threshold_avg_demp", gp, 800, 750, true, false, false);
	}
	
	private static DecimalFormat wtDF = new DecimalFormat("0.00E0");

}
