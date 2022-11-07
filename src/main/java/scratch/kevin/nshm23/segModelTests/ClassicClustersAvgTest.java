package scratch.kevin.nshm23.segModelTests;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;

public class ClassicClustersAvgTest {

	public static void main(String[] args) throws IOException {
		File resultsDir = new File(args[0]);
		
		AveragingAccumulator<ConnectivityClusters> accumulator = null;
		for (File dir : resultsDir.listFiles()) {
			if (!dir.getName().contains("Classic"))
				continue;
			File solFile = new File(dir, "solution.zip");
			FaultSystemRupSet rupSet = FaultSystemRupSet.load(solFile);
			
			ConnectivityClusters clusters = rupSet.requireModule(ConnectivityClusters.class);
			
			if (accumulator == null)
				accumulator = clusters.averagingAccumulator();
			
			accumulator.process(clusters, 1d);
		}
		
		accumulator.getAverage();
		
		System.out.println("DONE");
	}

}
