package scratch.kevin.prvi25.figures;

import java.awt.Color;
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
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class CombinedMFDsPlot {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution crustalSol = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED);
		
		FaultSystemSolution subductionSol1 = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		FaultSystemSolution subductionSol2 = FaultSystemSolution.load(SUBDUCTION_SOL_SMALL);
		
		File outputDir = new File(FIGURES_DIR, "combined_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(9.45);
		
		List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
		List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Region crustalReg = PRVI25_SeismicityRegions.CRUSTAL.load();
		
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
		
		IncrementalMagFreqDist faultMuertos = average(
				calcFaultMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), subductionSol1, refMFD),
				calcFaultMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), subductionSol2, refMFD));
		faultMuertos.setName("Muertos, Interface Fault");
		incrFuncs.add(faultMuertos);
		cmlFuncs.add(faultMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
		
		IncrementalMagFreqDist griddedMuertos = average(
				calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol1, refMFD),
				calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol2, refMFD));
		griddedMuertos.setName("Muertos, Interface Gridded");
		incrFuncs.add(griddedMuertos);
		cmlFuncs.add(griddedMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightgreen));
		
		IncrementalMagFreqDist slabMuertos = calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTRASLAB.load(),
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol1, refMFD);
		slabMuertos.setName("Muertos, Slab");
		incrFuncs.add(slabMuertos);
		cmlFuncs.add(slabMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightgreen));
		
		IncrementalMagFreqDist faultCar = average(
				calcFaultMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), subductionSol1, refMFD),
				calcFaultMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), subductionSol2, refMFD));
		faultCar.setName("Caribbean, Interface Fault");
		incrFuncs.add(faultCar);
		cmlFuncs.add(faultCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
		
		IncrementalMagFreqDist griddedCar = average(
				calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol1, refMFD),
				calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol2, refMFD));
		griddedCar.setName("Caribbean, Interface Gridded");
		incrFuncs.add(griddedCar);
		cmlFuncs.add(griddedCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightorange));
		
		IncrementalMagFreqDist slabCar = calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTRASLAB.load(),
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol1, refMFD);
		slabCar.setName("Caribbean, Slab");
		incrFuncs.add(slabCar);
		cmlFuncs.add(slabCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightorange));
		
		SummedMagFreqDist sum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (IncrementalMagFreqDist mfd : incrFuncs)
			sum.addIncrementalMagFreqDist(mfd);
		sum.setName("Total");
		incrFuncs.add(0, sum);
		cmlFuncs.add(0, sum.getCumRateDistWithOffset());
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
		
		PlotSpec incrSpec = new PlotSpec(incrFuncs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
		incrSpec.setLegendInset(RectangleAnchor.TOP_RIGHT);
		PlotSpec cmlSpec = new PlotSpec(cmlFuncs, chars, " ", "Magnitude", "Cumulative Rate (1/yr)");
		cmlSpec.setLegendInset(RectangleAnchor.TOP_RIGHT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();

		Range yRange = new Range(1e-6, 1e1);
		Range xRange = new Range(5d, 9.5d);

		gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
		PlotUtils.writePlots(outputDir, "combined_mfds", gp, 800, 750, true, false, false);
		gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
		PlotUtils.writePlots(outputDir, "combined_mfds_cml", gp, 800, 750, true, false, false);
	}
	
	static IncrementalMagFreqDist calcFaultMFD(Region region, FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) {
		return sol.calcNucleationMFD_forRegion(region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
	}
	
	static IncrementalMagFreqDist calcGriddedMFD(Region region, TectonicRegionType trt,
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
	
	static IncrementalMagFreqDist average(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(mfd1.getMinX(), mfd1.size(), mfd1.getDelta());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, 0.5*mfd1.getY(i) + 0.5*mfd2.getY(i));
		return ret;
	}

}
