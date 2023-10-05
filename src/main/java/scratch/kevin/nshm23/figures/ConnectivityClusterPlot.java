package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.FaultSectionConnectionsPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Stopwatch;

class ConnectivityClusterPlot {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		
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
		
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		reg = new Region(new Location(reg.getMinLat(), reg.getMinLon()), new Location(reg.getMaxLat(), reg.getMaxLon()));
		GeographicMapMaker plotter = FaultSectionConnectionsPlot.buildConnectedClustersPlot(rupSet, sol, reg, null, clustersUnsorted);
		
		plotter.setWritePDFs(true);
		plotter.setWriteGeoJSON(false);
		
		plotter.plot(outputDir, "connectivity_clusters", " ", 1200);
	}

}
