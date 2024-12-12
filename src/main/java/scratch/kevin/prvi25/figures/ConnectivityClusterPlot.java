package scratch.kevin.prvi25.figures;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.FaultSectionConnectionsPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

class ConnectivityClusterPlot {
	
	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(CRUSTAL_SOL_SUPRA_ONLY);
		
		File outputDir = new File(FIGURES_DIR, "crustal_fm");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		if (!rupSet.hasModule(ConnectivityClusters.class)) {
			System.out.println("Calculating connection clusters");
			Stopwatch watch = Stopwatch.createStarted();
			ConnectivityClusters clusters = ConnectivityClusters.build(rupSet);
			watch.stop();
			System.out.println("Found "+clusters.size()+" connectivity clusters in "+(float)(watch.elapsed(TimeUnit.MILLISECONDS)/1000)+" s");
			rupSet.addModule(clusters);
		}
		List<ConnectivityCluster> clustersUnsorted = rupSet.requireModule(ConnectivityClusters.class).get();
		
		// see if any clusters rupture completely
		System.out.println("Looking for clusters that rupture fully...");
		boolean anyFound = false;
		for (ConnectivityCluster cluster : clustersUnsorted) {
			HashSet<Integer> parents = new HashSet<>(cluster.getParentSectIDs());
			if (parents.size() == 1)
				continue;
			Set<Integer> sects = cluster.getSectIDs();
			for (List<Integer> rupSects : rupSet.getSectionIndicesForAllRups()) {
				if (rupSects.size() == sects.size() && rupSects.containsAll(sects)) {
					System.out.println("Cluster ruptures fully: "+sects);
					anyFound = true;
					System.out.println("Parents:");
					for (int sectID: rupSects) {
						FaultSection sect = rupSet.getFaultSectionData(sectID);
						if (parents.contains(sect.getParentSectionId())) {
							System.out.println("\t"+sect.getParentSectionId()+". "+sect.getParentSectionName());
							parents.remove(sect.getParentSectionId());
						}
					}
				}
			}
		}
		if (!anyFound)
			System.out.println("No clusters rupture fully");
		
		GeographicMapMaker plotter = FaultSectionConnectionsPlot.buildConnectedClustersPlot(rupSet, sol,
				SlipRateFigures.CRUSTAL_FAULT_MAP_REG, null, clustersUnsorted);
		
		plotter.setWritePDFs(true);
		plotter.setWriteGeoJSON(false);
		
		plotter.plot(outputDir, "connectivity_clusters", " ");
		
		CrustalFaultNamesFigure.addStandardFaultLabels(plotter, rupSet.getFaultSectionDataList());
		
		plotter.plot(outputDir, "connectivity_clusters_with_fault_names", " ");
	}

}
