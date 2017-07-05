package scratch.kevin.ucerf3.inversion;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class FaultModelCacheTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		for (FaultModels fm : FaultModels.values())
			fm.fetchFaultSections();
		System.exit(0);
	}

}
