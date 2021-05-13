package scratch.kevin.simulators.dists;

import java.util.List;

import org.opensha.sha.simulators.SimulatorEvent;

public interface ProbabilisticReturnPeriodProvider extends
		RandomReturnPeriodProvider {
	
	public PossibleRupture getPossibleRupture(List<SimulatorEvent> prevEvents, double windowStart, double windowEnd);
	
	public double getPreferredWindowLength();

}
