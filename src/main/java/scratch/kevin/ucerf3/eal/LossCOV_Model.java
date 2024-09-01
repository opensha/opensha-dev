package scratch.kevin.ucerf3.eal;

import java.awt.Color;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

import com.google.common.base.Preconditions;

public enum LossCOV_Model {
	
//	// not implementing this because it requires the x-axis values of the func supplied in calcLossExceedanceProbs to be modified too
//	PORTER_POWER_LAW_2024_08_09 { // meanLoss is loss (in dollars) divided by total portfolio value in units of $1000
//		@Override
//		public double getCOV(double meanLoss) {
//			return Math.min(2.026, 0.9832*Math.pow(meanLoss, -0.117));
//		}
//	},
	PORTER_POWER_LAW_2020_09_01 { // Here, meanLoss is loss in units of $1000
		@Override
		public double getCOV(double meanLoss) {
			return 4.4751*Math.pow(Math.max(1000, meanLoss), -0.1152);
		}
	},
	PORTER_POWER_LAW_2020_09_01_fixed { // fixed to be consistent with Keith's manuscripts; Here, meanLoss is loss in units of $1000
		@Override
		public double getCOV(double meanLoss) {
			return 4.546*Math.pow(Math.max(1000, meanLoss), -0.117);
		}
	},
	PORTER_POWER_LAW_2020_09_01_DOLLARS { // here, meanLoss is in dollars (not units of $1000)
		@Override
		public double getCOV(double meanLoss) {
			return 4.546*Math.pow(Math.max(1000, meanLoss/1000d), -0.117);
		}
	};

	
	public abstract double getCOV(double meanLoss);
	
	public LogNormalDistribution getDistribution(double meanLoss) {
		return getDistribution(meanLoss, null);
	}
	
	public LogNormalDistribution getDistribution(double meanLoss, RandomGenerator randomGen) {
		Preconditions.checkState(meanLoss > 0d);
		double cov = getCOV(meanLoss);
		double sigma = Math.sqrt(Math.log(cov*cov+1));
		double mu = Math.log(meanLoss)-(sigma*sigma/2);
		if(randomGen == null)
			return new LogNormalDistribution(mu, sigma);
		else
			return new LogNormalDistribution(randomGen, mu, sigma);
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
		LossCOV_Model model = PORTER_POWER_LAW_2020_09_01_fixed;
		
		// test lognormal distribution
        double mean = 1e5;
        double cov = model.getCOV(mean);
        double samples[] =      model.getDistribution(mean).sample(100000000);
        DescriptiveStatistics stats = new DescriptiveStatistics(samples);
        double meanFromDist = stats.getMean();
        double covFromDist = stats.getStandardDeviation()/meanFromDist;
       System.out.println("mean="+mean+";\tmeanFromDist="+meanFromDist+";\tcov="+cov+";\tcovFromDist="+
    		   covFromDist+"; meanRatio="+(float)(meanFromDist/mean)+"; covRatio="+(float)(covFromDist/cov));
		
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
