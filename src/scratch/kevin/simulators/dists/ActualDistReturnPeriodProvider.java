package scratch.kevin.simulators.dists;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.google.common.primitives.Doubles;

public class ActualDistReturnPeriodProvider implements RandomReturnPeriodProvider {
	
	private int index = 0;
	private double[] random_rps;
	
	public ActualDistReturnPeriodProvider(double[] rps) {
		random_rps = Arrays.copyOf(rps, rps.length);
		List<Double> randomized = Doubles.asList(random_rps);
		Collections.shuffle(randomized);
		random_rps = Doubles.toArray(randomized);
	}

	@Override
	public double getReturnPeriod() {
		if (index == random_rps.length)
			index = 0;
		return random_rps[index++];
	}
	
	public double[] getRPs() {
		return random_rps;
	}
}