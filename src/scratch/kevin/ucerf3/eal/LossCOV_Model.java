package scratch.kevin.ucerf3.eal;

import java.awt.Color;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

import com.google.common.base.Preconditions;

public enum LossCOV_Model {
	
	PORTER_POWER_LAW_2020_09_01 {
		@Override
		public double getCOV(double meanLoss) {
			return 4.4751*Math.pow(Math.max(1000, meanLoss), -0.1152);
		}
	};
	
	public abstract double getCOV(double meanLoss);
	
	public LogNormalDistribution getDistribution(double meanLoss) {
		Preconditions.checkState(meanLoss > 0d);
		return new LogNormalDistribution(Math.log(meanLoss), getCOV(meanLoss));
	}
	
	public DiscretizedFunc calcLossExceedanceProbs(DiscretizedFunc xVals, double meanLoss) {
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		LogNormalDistribution dist = getDistribution(meanLoss);
		for (int i=0; i<xVals.size(); i++) {
			double x = xVals.getX(i);
			double y = 1d - dist.cumulativeProbability(x);
			ret.set(x, y);
		}
		return ret;
	}
	
	public double calcLossExceedanceProb(double x, double meanLoss) {
		LogNormalDistribution dist = getDistribution(meanLoss);
		return 1d - dist.cumulativeProbability(x);
	}
	
	public static void main(String[] args) {
		LossCOV_Model model = PORTER_POWER_LAW_2020_09_01;
		
		EvenlyDiscretizedFunc logFunc = new EvenlyDiscretizedFunc(0d, 8d, 500);
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		
		for (int i=0; i<logFunc.size(); i++) {
			double x = Math.pow(10, logFunc.getX(i));
			func.set(x, model.getCOV(x));
		}
		
		GraphWindow gw = new GraphWindow(func, "Power-Law COV Model",
				new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		gw.setAxisRange(new Range(func.getMinX(), func.getMaxX()), null);
		gw.setXLog(true);
		gw.setYLog(false);
		gw.setX_AxisLabel("Mean loss (thousands)");
		gw.setY_AxisLabel("COV");
		
		double[] lossVals = { 1d, 10d, 100d, 1000d, 10000d };
		for (double meanLoss : lossVals) {
			System.out.println("Mean loss: "+(float)meanLoss);
			for (double x : lossVals)
				System.out.println("\tP(>"+(float)x+")="+(float)model.calcLossExceedanceProb(x, meanLoss));
			DiscretizedFunc exceeds = model.calcLossExceedanceProbs(func, meanLoss);
			
			gw = new GraphWindow(exceeds, "Power-Law COV Model, meanLoss="+(float)meanLoss,
					new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			gw.setAxisRange(new Range(func.getMinX(), func.getMaxX()), null);
			gw.setXLog(true);
			gw.setYLog(false);
			gw.setX_AxisLabel("Mean loss (thousands)");
			gw.setY_AxisLabel("Exceedance Prob");
		}
	}

}
