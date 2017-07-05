package scratch.kevin.simulators.dists;

import java.util.List;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class LogNormalDistReturnPeriodProvider implements RandomReturnPeriodProvider {
		
		private LogNormalDistribution n;
		
		public static double[] getTrimmedRPs(double[] rps, double mean) {
			double trimAmount = mean*0.4;
			double trimMin = mean-trimAmount;
			double trimMax = mean+trimAmount;
			List<Double> trimmedRPs = Lists.newArrayList();
			for (double rp : rps)
				if (rp>=trimMin && rp<=trimMax)
					trimmedRPs.add(rp);
			return Doubles.toArray(trimmedRPs);
		}
		
		public LogNormalDistReturnPeriodProvider(double[] rps) {
			double mean = StatUtils.mean(rps);
			double var = StatUtils.variance(rps, mean);
			double sd = Math.sqrt(var);
			System.out.println("ORIG Mean: "+mean);
			System.out.println("ORIG Variance: "+StatUtils.variance(rps, mean));
			System.out.println("ORIG SD: "+sd);
			// trim dist
			rps = getTrimmedRPs(rps, mean);
			mean = StatUtils.mean(rps);
			var = StatUtils.variance(rps, mean);
			sd = Math.sqrt(var);
//			n = new LogNormalDistribution(mean, sd);
//			n = new LogNormalDistribution(Math.log(mean), Math.log(sd));
			System.out.println("Mean: "+mean);
			System.out.println("Variance: "+var);
			System.out.println("SD: "+sd);
			double shape = sd / mean;
//			double shape = 0.2;
//			mean = mean-10;
			System.out.println("Shape: "+shape);
			n = new LogNormalDistribution(Math.log(mean), shape);
//			n = new LogNormalDistribution(mean, 1d);
			
//			System.out.println("Log-Normal Distribution. mean="+mean+", std dev="+sd+", num_mean="+n.getNumericalMean());
//			for (int i=0; i<10; i++)
//				System.out.println("\t"+i+". "+getReturnPeriod());
		}
		
		public LogNormalDistReturnPeriodProvider(double scale, double shape) {
			n = new LogNormalDistribution(scale, shape);
		}
		
		public void setSeed(long seed) {
			n.reseedRandomGenerator(seed);
		}

		@Override
		public double getReturnPeriod() {
			return n.sample();
		}
		
		@Override
		public String toString() {
			return "LogNormal[scale="+n.getScale()+", mean="+Math.exp(n.getScale())+", shape="+n.getShape()+"]";
		}
		
		public LogNormalDistribution getDist() {
			return n;
		}
	}