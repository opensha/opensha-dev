package scratch.kevin.nshm23.prvi.figures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import net.mahdilamb.colormap.Colors;

public class CombinedMFDsPlot {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution crustalSol = FaultSystemSolution.load(new File(""
				+ "/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_07_26-prvi25_crustal_branches-dmSample5x/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		FaultSystemSolution subductionSol = FaultSystemSolution.load(new File(""
				+ "/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_08_01-prvi25_subduction_branches/results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip"));
		
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(9.45);
		
		List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
		List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Region crustalReg = SeismicityRegions.CRUSTAL.load();
		
		IncrementalMagFreqDist crustalFault = calcFaultMFD(crustalReg, crustalSol, refMFD);
		crustalFault.setName("Crustal, Fault");
		incrFuncs.add(crustalFault);
		cmlFuncs.add(crustalFault.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
		
		IncrementalMagFreqDist crustalGrid = calcGriddedMFD(crustalReg, TectonicRegionType.ACTIVE_SHALLOW, crustalSol, refMFD);
		crustalGrid.setName("Crustal, Gridded");
		incrFuncs.add(crustalGrid);
		cmlFuncs.add(crustalGrid.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightblue));
		
		IncrementalMagFreqDist slabMuertos = calcGriddedMFD(SeismicityRegions.MUE_INTRASLAB.load(),
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol, refMFD);
		slabMuertos.setName("Muertos, Slab");
		incrFuncs.add(slabMuertos);
		cmlFuncs.add(slabMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightgreen));
		
		IncrementalMagFreqDist griddedMuertos = calcGriddedMFD(SeismicityRegions.MUE_INTERFACE.load(),
				TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol, refMFD);
		griddedMuertos.setName("Muertos, Interface Gridded");
		incrFuncs.add(griddedMuertos);
		cmlFuncs.add(griddedMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightgreen));
		
		IncrementalMagFreqDist faultMuertos = calcFaultMFD(SeismicityRegions.MUE_INTERFACE.load(), subductionSol, refMFD);
		faultMuertos.setName("Muertos, Interface Fault");
		incrFuncs.add(faultMuertos);
		cmlFuncs.add(faultMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
		
		
		IncrementalMagFreqDist slabCar = calcGriddedMFD(SeismicityRegions.CAR_INTRASLAB.load(),
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol, refMFD);
		slabCar.setName("Caribbean, Slab");
		incrFuncs.add(slabCar);
		cmlFuncs.add(slabCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightorange));
		
		IncrementalMagFreqDist griddedCar = calcGriddedMFD(SeismicityRegions.CAR_INTERFACE.load(),
				TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol, refMFD);
		griddedCar.setName("Caribbean, Interface Gridded");
		incrFuncs.add(griddedCar);
		cmlFuncs.add(griddedCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightorange));
		
		IncrementalMagFreqDist faultCar = calcFaultMFD(SeismicityRegions.CAR_INTERFACE.load(), subductionSol, refMFD);
		faultCar.setName("Caribbean, Interface Fault");
		incrFuncs.add(faultCar);
		cmlFuncs.add(faultCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
		
		PlotSpec incrSpec = new PlotSpec(incrFuncs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
		incrSpec.setLegendInset(RectangleAnchor.TOP_RIGHT);
		PlotSpec cmlSpec = new PlotSpec(cmlFuncs, chars, " ", "Magnitude", "Cumulative Rate (1/yr)");
		cmlSpec.setLegendInset(RectangleAnchor.TOP_RIGHT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();

		File outputDir = new File("/tmp");
		Range yRange = new Range(1e-6, 1e1);
		Range xRange = new Range(5d, 9.5d);

		gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
		PlotUtils.writePlots(outputDir, "combined_mfds", gp, 800, 750, true, false, false);
		gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
		PlotUtils.writePlots(outputDir, "combined_mfds_cml", gp, 800, 750, true, false, false);
	}
	
	private static IncrementalMagFreqDist calcFaultMFD(Region region, FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) {
		return sol.calcNucleationMFD_forRegion(region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
	}
	
	private static IncrementalMagFreqDist calcGriddedMFD(Region region, TectonicRegionType trt,
			FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) {
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		SummedMagFreqDist ret = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int i=0; i<gridProv.getNumLocations(); i++) {
			if (region.contains(gridProv.getLocation(i))) {
				IncrementalMagFreqDist mfd = gridProv.getMFD(trt, i);
				ret.addIncrementalMagFreqDist(mfd);
			}
		}
		return ret;
	}
	
	private static IncrementalMagFreqDist average(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(mfd1.getMinX(), mfd1.size(), mfd1.getDelta());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, 0.5*mfd1.getY(i) + 0.5*mfd2.getY(i));
		return ret;
	}

}
