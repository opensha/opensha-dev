package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.BiasiWesnouskyJumpProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc.DistDependentJumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;

public class SegWeightSearch {

	public static void main(String[] args) throws IOException {
		List<DistDependentJumpProbabilityCalc> calcs = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		
//		DistDependentJumpProbabilityCalc target = new Shaw07JumpDistProb(1, 3d);
		
//		calcs.add(new Shaw07JumpDistProb(1, 1d));
//		weights.add(0.15);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 3d));
//		weights.add(0.7);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 6d));
//		weights.add(0.15);
		
//		calcs.add(new Shaw07JumpDistProb(1, 2d));
//		weights.add(0.25);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 3d));
//		weights.add(0.6);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 4d));
//		weights.add(0.15);
		
//		calcs.add(new Shaw07JumpDistProb(1, 2d));
//		weights.add(0.3);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 3d));
//		weights.add(0.6);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 5d));
//		weights.add(0.10);
		
//		calcs.add(new Shaw07JumpDistProb(1, 2d));
//		weights.add(0.25);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 3d));
//		weights.add(0.6);
//		
//		calcs.add(new Shaw07JumpDistProb(1, 5d));
//		weights.add(0.15);
//		
//		calcs.add(new BiasiWesnouskyJumpProb.BiasiWesnousky2016SSJumpProb(1d));
//		weights.add(0.2);
		
		DistDependentJumpProbabilityCalc target = Shaw07JumpDistProb.forHorzOffset(1d, 3d, 2d);
		
		calcs.add(Shaw07JumpDistProb.forHorzOffset(1d, 2d, 1d));
//		calcs.add(new Shaw07JumpDistProb(1d, 2d));
//		weights.add(0.3);
		weights.add(1d);
		
		calcs.add(Shaw07JumpDistProb.forHorzOffset(1d, 3d, 2d));
//		weights.add(0.5);
		weights.add(1d);
		
		calcs.add(Shaw07JumpDistProb.forHorzOffset(1d, 4d, 3d));
//		weights.add(0.2);
		weights.add(1d);
		
		plotWeights(calcs, weights, target);
	}
	
	private static final EvenlyDiscretizedFunc xVals = new EvenlyDiscretizedFunc(0d, 15d, 5000);
	
	private static EvenlyDiscretizedFunc calcWeighted(List<DistDependentJumpProbabilityCalc> calcs, List<Double> weights) {
		EvenlyDiscretizedFunc weightFunc = new EvenlyDiscretizedFunc(xVals.getMinX(), xVals.getMaxX(), xVals.size());
		
		double totWeight = 0d;
		
		for (int i=0; i<calcs.size(); i++) {
			DistDependentJumpProbabilityCalc calc = calcs.get(i);
			double weight = weights.get(i);
			
			totWeight += weight;
			for (int j=0; j<weightFunc.size(); j++)
				weightFunc.add(j, calc.calcJumpProbability(weightFunc.getX(j))*weight);
		}
		if (totWeight != 1d)
			weightFunc.scale(1d/totWeight);
		return weightFunc;
	}
	
	private static void plotWeights(List<DistDependentJumpProbabilityCalc> calcs, List<Double> weights,
			DistDependentJumpProbabilityCalc target) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		EvenlyDiscretizedFunc weightFunc = calcWeighted(calcs, weights);
		
		for (DistDependentJumpProbabilityCalc calc : calcs) {
			EvenlyDiscretizedFunc subFunc = new EvenlyDiscretizedFunc(xVals.getMinX(), xVals.getMaxX(), xVals.size());
			for (int j=0; j<subFunc.size(); j++)
				subFunc.set(j, calc.calcJumpProbability(subFunc.getX(j)));
			funcs.add(subFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		}
		
		EvenlyDiscretizedFunc targetFunc = new EvenlyDiscretizedFunc(xVals.getMinX(), xVals.getMaxX(), xVals.size());
		
		for (int j=0; j<targetFunc.size(); j++)
			targetFunc.set(j, target.calcJumpProbability(targetFunc.getX(j)));
		
		funcs.add(targetFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		funcs.add(weightFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Jump Distance Models", "Jump Distance (km)", "Jump Probability");
		spec.setLegendInset(true);
		
		Range xRange = new Range(xVals.getMinX(), xVals.getMaxX());
		Range yRange = new Range(0d, 1d);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File outputDir = new File("/tmp");
		String prefix = "jump_weight_test";
		
		PlotUtils.writePlots(outputDir, prefix, gp, 1000, 850, true, true, false);
		
		yRange = new Range(1e-4, 1);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_log", gp, 1000, 850, true, true, false);
	}

}
