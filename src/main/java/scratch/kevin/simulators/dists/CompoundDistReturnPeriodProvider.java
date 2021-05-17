package scratch.kevin.simulators.dists;

import java.util.List;

import com.google.common.base.Preconditions;

public class CompoundDistReturnPeriodProvider implements RandomReturnPeriodProvider {
	
	private List<RandomReturnPeriodProvider> provs;
	private List<Double> weights;
	private double[] weightEnds;
	
	public CompoundDistReturnPeriodProvider(
			List<RandomReturnPeriodProvider> provs, List<Double> weights) {
		Preconditions.checkState(provs.size() == weights.size());
		double tot = 0;
		for (double weight : weights)
			tot += weight;
		
		weightEnds = new double[provs.size()];
		double running = 0;
		for (int i=0; i<weights.size(); i++) {
			double weight = weights.get(i);
			running += weight;
			weightEnds[i] = running/tot;
		}
		this.provs = provs;
		this.weights = weights;
	}

	@Override
	public double getReturnPeriod() {
		double r = Math.random();
		int i;
		for (i=0; i<weightEnds.length-1; i++) {
			if (r < weightEnds[i])
				break;
		}
		return provs.get(i).getReturnPeriod();
	}

	public List<RandomReturnPeriodProvider> getProvs() {
		return provs;
	}

	public List<Double> getWeights() {
		return weights;
	}
	
}