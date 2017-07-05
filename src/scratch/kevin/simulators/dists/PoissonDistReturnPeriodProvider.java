package scratch.kevin.simulators.dists;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.StatUtils;

public class PoissonDistReturnPeriodProvider implements RandomReturnPeriodProvider {
	
	private PoissonDistribution n;
	
	public PoissonDistReturnPeriodProvider(double[] rps) {
		double mean = StatUtils.mean(rps);
		double sd = Math.sqrt(StatUtils.variance(rps, mean));
		n = new PoissonDistribution(mean);
//		System.out.println("Poisson Distribution. mean="+mean+", std dev="+sd);
//		for (int i=0; i<10; i++)
//			System.out.println("\t"+i+". "+getReturnPeriod());
	}

	@Override
	public double getReturnPeriod() {
		return n.sample();
	}
}