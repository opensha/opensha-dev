package scratch.kevin.simulators.catBuild;

import java.util.List;

import org.opensha.sha.simulators.SimulatorEvent;

import scratch.kevin.simulators.dists.RandomReturnPeriodProvider;

public interface CatalogBuilder {
	public List<SimulatorEvent> buildCatalog(List<? extends SimulatorEvent> events,
			List<RandomReturnPeriodProvider> randomRPsList,
			List<List<? extends SimulatorEvent>> eventListsToResample, boolean trim);
}