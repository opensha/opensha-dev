package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class StateByStateMFDCarveouts {

	public static void main(String[] args) throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip");
		File outputDir = new File("/tmp/nshm23_mfd_carveout");
//				+ "2024_02_02-nshm23_branches-WUS_FM_v3/node_branch_averaged/SegModel_Classic.zip");
//		File outputDir = new File("/tmp/nshm23_mfd_carveout_classic");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<NSHM23_SingleStates> states = new ArrayList<>();
		states.add(null);
		for (NSHM23_SingleStates state : NSHM23_SingleStates.values())
			states.add(state);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 8.01);
		Range xRange = new Range(5d, 8d);
		Range yRange = new Range(1e-5, 3e0);
		
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		
		for (NSHM23_SingleStates state : states) {
			Region region = state == null ? NSHM23_RegionLoader.SeismicityRegions.CONUS_WEST.load() : state.loadRegion();
			
			IncrementalMagFreqDist onFaultMFD = sol.calcNucleationMFD_forRegion(
					region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size(), false);
			SummedMagFreqDist griddedMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			
			for (int l=0; l<gridProv.getNumLocations(); l++) {
				Location loc = gridProv.getLocation(l);
				if (region.contains(loc)) {
					IncrementalMagFreqDist gridMFD = gridProv.getMFD(l, refMFD.getMinX());
					if (gridMFD != null)
						griddedMFD.addIncrementalMagFreqDist(gridMFD);
				}
			}
			SummedMagFreqDist totalMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			totalMFD.addIncrementalMagFreqDist(onFaultMFD);
			totalMFD.addIncrementalMagFreqDist(griddedMFD);
			
			List<IncrementalMagFreqDist> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			onFaultMFD.setName("On-fault supra-seis");
			funcs.add(onFaultMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
			
			griddedMFD.setName("Gridded seismicity");
			funcs.add(griddedMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
			
			totalMFD.setName("Total model");
			funcs.add(totalMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
			
			String title, prefix;
			if (state == null) {
				title = "Full Collection Region";
				prefix = "full";
			} else {
				title = state.getName();
				prefix = state.getFilePrefix();
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", "Incremental Rate (1/yr)");
			plot.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plot, false, true, xRange, yRange);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 800, true, true, false);
		}
	}

}
