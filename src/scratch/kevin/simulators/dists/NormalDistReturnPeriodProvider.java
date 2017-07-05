package scratch.kevin.simulators.dists;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

public class NormalDistReturnPeriodProvider implements RandomReturnPeriodProvider {
	
	private NormalDistribution n;
	
	public NormalDistReturnPeriodProvider(double[] rps) {
		double mean = StatUtils.mean(rps);
		double sd = Math.sqrt(StatUtils.variance(rps, mean));
		n = new NormalDistribution(mean, sd);
	}

	@Override
	public double getReturnPeriod() {
		return n.sample();
	}
}