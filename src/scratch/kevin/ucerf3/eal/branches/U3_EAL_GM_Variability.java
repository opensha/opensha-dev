package scratch.kevin.ucerf3.eal.branches;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

public enum U3_EAL_GM_Variability implements LogicTreeBranchNode<U3_EAL_GM_Variability> {
	NEG2("-2", "Neg2", 0.1, -2d),
	NEG1("-1", "Neg1", 0.1, -1d),
	ZERO("0", "Zero", 0.6, 0d),
	POS1("+1", "Pos1", 0.1, 1d),
	POS2("+2", "Pos2", 0.1, 2d);

	private String name;
	private String shortName;
	private double weight;
	private double val;
	
	private U3_EAL_GM_Variability(String name, String shortName, double weight, double val) {
		this.name = name;
		this.shortName = shortName;
		this.weight = weight;
		this.val = val;
	}
	
	@Override
	public String getShortName() {
		return shortName;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public double getRelativeWeight(InversionModels im) {
		return weight;
	}

	@Override
	public String encodeChoiceString() {
		return getShortName();
	}

	@Override
	public String getBranchLevelName() {
		return "Combined GM Variability";
	}
	
	public double calcGMVarLoss(double lambda) {
		return calcGMVarLoss(lambda, val);
	}
	
	// constants for GM variability
	private static final double a_sub_L = 0.1992;
	private static final double b_sub_L = 1.0696;
	private static final double c_sub_L = 3.1055;
	private static final double d_sub_L = -0.0935;
	private static final double e_sub_L = 0.091;

	private static double calcGMVarLoss(double lambda, double epsilon_sub_L) {
		if (lambda == 0d)
			return 0d;
		// eqn 35
		double theta_sub_L = a_sub_L * Math.pow(lambda, b_sub_L);
		// eqn 37
		double beta_sub_L0 = c_sub_L * Math.pow(theta_sub_L, d_sub_L);
		// eqn 38
		double beta_sub_L1 = e_sub_L;
		// eqn 36
		double beta_sub_L = Math.sqrt(beta_sub_L0*beta_sub_L0 + beta_sub_L1*beta_sub_L1);
		// eqn 34
		double A_sub_L = (1d/lambda) * theta_sub_L * Math.exp(epsilon_sub_L * beta_sub_L);
		// eqn 33
		return A_sub_L * lambda;
	}
	
	public static void main(String[] args) throws IOException {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(1000d, 1000000d, 1001);
		
		for (int i=0; i<func.size(); i++) {
			double lambda = func.getX(i);
//			System.out.println("Nominal loss: "+lambda);
			double wtVal = 0d;
			for (U3_EAL_GM_Variability var : values()) {
				double varLoss = var.calcGMVarLoss(lambda);
				double ratio = varLoss/lambda;
//				System.out.println(var.name+": varLoss = "+(float)varLoss+",\tratio="+(float)ratio);
				wtVal += varLoss*var.weight;
			}
//			System.out.println("Weight average: "+wtVal);
			func.set(i, wtVal);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = new Range(Math.min(func.getMinX(), func.getMinY()),
				Math.max(func.getMaxX(), func.getMaxY()));
		Range yRange = xRange;
		
		DiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
		oneToOne.set(xRange.getLowerBound(), xRange.getLowerBound());
		oneToOne.set(xRange.getUpperBound(), xRange.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, Color.GRAY));
		
		funcs.add(func);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "GM Var Bias", "Nominal Loss", "Weighted-Mean GMVar Loss");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File file = new File("/tmp/gm_var_bias");
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		file = new File("/tmp/gm_var_bias_log");
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}

}
