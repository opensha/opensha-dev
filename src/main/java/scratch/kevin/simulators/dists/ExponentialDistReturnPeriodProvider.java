package scratch.kevin.simulators.dists;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.stat.StatUtils;

public class ExponentialDistReturnPeriodProvider implements RandomReturnPeriodProvider {
	
	private ExponentialDistribution n;
	
	public ExponentialDistReturnPeriodProvider(double[] rps) {
		double mean = StatUtils.mean(rps);
		double sd = Math.sqrt(StatUtils.variance(rps, mean));
		n = new ExponentialDistribution(mean);
	}

	@Override
	public double getReturnPeriod() {
		return n.sample();
	}
}