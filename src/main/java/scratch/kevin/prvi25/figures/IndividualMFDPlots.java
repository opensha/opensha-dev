package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import net.mahdilamb.colormap.Colors;

public class IndividualMFDPlots {

	public static void main(String[] args) throws IOException {
		File crustalDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_11_19-prvi25_crustal_branches-dmSample5x");
		FaultSystemSolution crustalSol = FaultSystemSolution.load(new File(crustalDir,
				"results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(4.05, 9.55);
		
		File subductionDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_11_19-prvi25_subduction_branches");
		FaultSystemSolution subductionSol1 = FaultSystemSolution.load(new File(subductionDir,
				"results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip"));
		FaultSystemSolution subductionSol2 = FaultSystemSolution.load(new File(subductionDir,
				"results_PRVI_SUB_FM_SMALL_branch_averaged_gridded.zip"));
		
		File crustalOutputDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures/crustal_sol");
		Preconditions.checkState(crustalOutputDir.exists() || crustalOutputDir.mkdir());
		File subductionOutputDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures/sub_sol");
		Preconditions.checkState(subductionOutputDir.exists() || subductionOutputDir.mkdir());
		
		UncertainBoundedIncrMagFreqDist[] crustalOnFaultDists = getOnFaultMFDs(crustalSol);
		IncrementalMagFreqDist crustalMean = CombinedMFDsPlot.calcFaultMFD(null, crustalSol, refMFD);
		UncertainBoundedIncrMagFreqDist crustalObs = PRVI25_RegionalSeismicity.getBounded(PRVI25_SeismicityRegions.CRUSTAL, refMFD, 8.6);
		IncrementalMagFreqDist crustalGridded = CombinedMFDsPlot.calcGriddedMFD(
				PRVI25_SeismicityRegions.CRUSTAL.load(), TectonicRegionType.ACTIVE_SHALLOW, crustalSol, refMFD);
		
		plot(crustalOutputDir, "crustal_mfds", new Range(5d, 8.5d), crustalOnFaultDists, crustalMean, null, crustalGridded, crustalObs);
	}
	
	private static void plot(File outputDir, String prefix, Range xRange,
			UncertainBoundedIncrMagFreqDist[] onFaultDists, IncrementalMagFreqDist onFaultMean,
			UncertainBoundedIncrMagFreqDist[] griddedDists, IncrementalMagFreqDist gridded,
			UncertainBoundedIncrMagFreqDist obs) throws IOException {
		Color onFaultColor = Colors.tab_red;
		Color onFaultTransColor = new Color(onFaultColor.getRed(), onFaultColor.getGreen(), onFaultColor.getBlue(), 60);
		Color obsColor = Colors.tab_green;
		Color griddedColor = Colors.tab_blue;
		Color griddedTransColor = new Color(griddedColor.getRed(), griddedColor.getGreen(), griddedColor.getBlue(), 60);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (obs != null) {
			funcs.add(obs);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, obsColor));
			
			IncrementalMagFreqDist obsLower = obs.getLower();
			obsLower.setName(obs.getBoundName());
			funcs.add(obsLower);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, obsColor));
			
			IncrementalMagFreqDist obsUpper = obs.getUpper();
			obsUpper.setName(null);
			funcs.add(obsUpper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, obsColor));
		}
		
		onFaultMean.setName("On-Fault");
		funcs.add(onFaultMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, onFaultColor));
		
		if (onFaultDists != null) {
			for (int i=0; i<onFaultDists.length; i++) {
				if (i == 0)
					onFaultDists[i].setName(fractileLabel);
				else
					onFaultDists[i].setName(null);
				funcs.add(onFaultDists[i]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, onFaultTransColor));
			}
		}
		
		if (gridded != null) {
			gridded.setName("Gridded");
			funcs.add(gridded);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, griddedColor));
		}
		
		if (griddedDists != null) {
			for (int i=0; i<griddedDists.length; i++) {
				if (i == 0)
					griddedDists[i].setName(fractileLabel);
				else
					griddedDists[i].setName(null);
				funcs.add(griddedDists[i]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, griddedTransColor));
			}
		}
		
		// again on top
		onFaultMean = onFaultMean.deepClone();
		onFaultMean.setName(null);
		funcs.add(onFaultMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, onFaultColor));
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
		plot.setLegendInset(true);
		
		Range yRange = new Range(1e-6, 1e1);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
	}
	
	private static double[] fractiles = {0d, 0.025, 0.16, 0.84, 0.975d, 1d};
	private static String fractileLabel = "Sample p[0,2.5,16,84,97.5,100]";
	
	private static UncertainBoundedIncrMagFreqDist[] getOnFaultMFDs(FaultSystemSolution sol, int... parents) {
		BranchSectNuclMFDs branchSectMFDs = sol.requireModule(BranchSectNuclMFDs.class);
		
		List<Integer> sectIDs = new ArrayList<>();
		if (parents != null && parents.length > 0) {
			for (FaultSection sect : sol.getRupSet().getFaultSectionDataList())
				if (Ints.contains(parents, sect.getParentSectionId()))
					sectIDs.add(sect.getSectionId());
		} else {
			for (int i=0; i<sol.getRupSet().getNumSections(); i++)
				sectIDs.add(i);
		}
		IncrementalMagFreqDist[] fractileMFDs = branchSectMFDs.calcIncrementalSectFractiles(sectIDs, fractiles);
		
		int numRet = fractiles.length/2;
		UncertainBoundedIncrMagFreqDist[] ret = new UncertainBoundedIncrMagFreqDist[numRet];
		for (int i=0; i<numRet; i++) {
			IncrementalMagFreqDist lower = fractileMFDs[i];
			IncrementalMagFreqDist upper = fractileMFDs[fractileMFDs.length - (1 + i)];
			ret[i] = new UncertainBoundedIncrMagFreqDist(CombinedMFDsPlot.average(lower, upper), lower, upper, null);
		}
		
		return ret;
	}

}
