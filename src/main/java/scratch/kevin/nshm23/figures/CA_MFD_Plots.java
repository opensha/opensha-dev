package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

class CA_MFD_Plots {

	public static void main(String[] args) throws IOException {
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(8.95);
		
		File outputDir = new File("/tmp/u3_reg_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		AnalysisRegions analysis = AnalysisRegions.CONUS_U3_RELM;
		Region region = analysis.load();
		
		UncertainBoundedIncrMagFreqDist dataBounds = NSHM23_RegionalSeismicity.getRemapped(region,
				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
		dataBounds.setName("Observed");
		dataBounds.setBoundName("95% Bounds");
		
		FaultSystemSolution u3Sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/branch_avgs_combined.zip"));
		
		FaultSystemSolution methodsSol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_06-nshm23_u3_hybrid_branches-no_paleo_slip-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		
		FaultSystemSolution modelSol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_07-nshm23_branches-no_paleo_slip-mod_dm_weights-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-"
				+ "ThreshAvgIterRelGR/results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		List<XY_DataSet> incrFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
		List<XY_DataSet> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> cmlChars = new ArrayList<>();
		
		// add observed bounds
		Color obsColor = new Color(125, 80, 145); // "indigo"
		incrFuncs.add(dataBounds);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, obsColor));
		EvenlyDiscretizedFunc dataCumulative = dataBounds.getCumRateDistWithOffset();
		cmlFuncs.add(dataBounds.getCumRateDistWithOffset());
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, obsColor));
		dataBounds = dataBounds.deepClone();
		dataBounds.setName(dataBounds.getBoundName());
		
		incrFuncs.add(dataBounds);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
				new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60)));
		
		EvenlyDiscretizedFunc upperCumulative = dataBounds.getUpper().getCumRateDistWithOffset();
		EvenlyDiscretizedFunc lowerCumulative = dataBounds.getLower().getCumRateDistWithOffset();
		Preconditions.checkState(dataCumulative.size() == upperCumulative.size());
		for (int i=0; i<dataCumulative.size(); i++) {
			upperCumulative.set(i, Math.max(dataCumulative.getY(i), upperCumulative.getY(i)));
			lowerCumulative.set(i, Math.max(0, Math.min(dataCumulative.getY(i), lowerCumulative.getY(i))));
		}
		
		UncertainArbDiscFunc cmlBounded = new UncertainArbDiscFunc(dataCumulative, lowerCumulative, upperCumulative);
		cmlBounded.setName(dataBounds.getName());
		cmlFuncs.add(cmlBounded);
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
				new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60)));
		
		// add models
		addMFDs(u3Sol, region, refMFD, Color.BLUE, "UCERF3", incrFuncs, incrChars, cmlFuncs, cmlChars);
		addMFDs(methodsSol, region, refMFD, Color.GREEN.darker(), "NSHM23 Methodology", incrFuncs, incrChars, cmlFuncs, cmlChars);
		addMFDs(modelSol, region, refMFD, Color.RED, "NSHM23 Model", incrFuncs, incrChars, cmlFuncs, cmlChars);
		
		Range xRange = new Range(6d, 8.5d);
		Range yRange = new Range(1e-6, 2e0);
		
		// add full model bounds
		BranchRegionalMFDs branchMFDs = modelSol.requireModule(BranchRegionalMFDs.class);
		List<Region> rois = modelSol.getRupSet().requireModule(RegionsOfInterest.class).getRegions();
		for (int i=0; i<rois.size(); i++) {
			Region testReg = rois.get(i);
			if (region.equalsRegion(testReg)) {
				System.out.println("Found region match!");
				
				IncrementalMagFreqDist[] incrPercentiles = branchMFDs.calcRegionalIncrementalFractiles(
						MFDType.SUPRA_ONLY, i, SolMFDPlot.standardFractiles);
				
				Color transColor = new Color(255, 0, 0, 30);
				PlotCurveCharacterstics minMaxChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, transColor);
				
				for (IncrementalMagFreqDist bounds : SolMFDPlot.processIncrFractiles(incrPercentiles)) {
					incrFuncs.add(bounds);
					incrChars.add(minMaxChar);
				}
				
				EvenlyDiscretizedFunc[] cmlPercentiles = branchMFDs.calcRegionalCumulativeFractiles(
						MFDType.SUPRA_ONLY, i, SolMFDPlot.standardFractiles);
				for (UncertainArbDiscFunc cmlBounds : SolMFDPlot.processCmlFractiles(cmlPercentiles, xRange.getLowerBound())) {
					cmlFuncs.add(cmlBounds);
					cmlChars.add(minMaxChar);
				}
			}
		}
		
		// add model target
		FaultSystemRupSet modelRupSet = modelSol.getRupSet();
		List<? extends IncrementalMagFreqDist> targetMFDs = modelRupSet.requireModule(InversionTargetMFDs.class).getOnFaultSupraSeisNucleationMFDs();
		IncrementalMagFreqDist target = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		double[] modelSectFracts = modelRupSet.getFractSectsInsideRegion(region, false);
		for (int s=0; s<modelRupSet.getNumSections(); s++) {
			if (modelSectFracts[s] > 0) {
				IncrementalMagFreqDist mfd = targetMFDs.get(s);
				for (int i=0; i<mfd.size(); i++)
					target.set(i, target.getY(i)+mfd.getY(i)*modelSectFracts[s]);
			}
		}
		target.setName("NSHM23 Model Target");
		incrFuncs.add(target);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.GRAY));
		cmlFuncs.add(target.getCumRateDistWithOffset());
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.GRAY));
		
		// now add again on top without names
		addMFDs(u3Sol, region, refMFD, Color.BLUE, null, incrFuncs, incrChars, cmlFuncs, cmlChars);
		addMFDs(methodsSol, region, refMFD, Color.GREEN.darker(), null, incrFuncs, incrChars, cmlFuncs, cmlChars);
		addMFDs(modelSol, region, refMFD, Color.RED, null, incrFuncs, incrChars, cmlFuncs, cmlChars);
		
		PlotSpec incrSpec = new PlotSpec(incrFuncs, incrChars, " ", "Magnitude", "Incremental Rate (/yr)");
		incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		PlotSpec cmlSpec = new PlotSpec(cmlFuncs, cmlChars, " ", "Magnitude", "Cumulative Rate (/yr)");
		cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setAxisLabelFontSize(26);
		gp.setTickLabelFontSize(22);
		
		gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "u3_reg_mfds_incr", gp, 800, 750, true, true, false);
		
		gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "u3_reg_mfds_cml", gp, 800, 750, true, true, false);
		
		for (XY_DataSet func : cmlFuncs) {
			if (func.getName() != null) {
				for (Point2D pt : func) {
					if ((float)pt.getX() == (float)xRange.getLowerBound()) {
						System.out.println(func.getName()+": rate M>"+(float)pt.getX()+"="+(float)pt.getY());
					}
				}
			}
		}
	}
	
	private static void addMFDs(FaultSystemSolution sol, Region region, EvenlyDiscretizedFunc refMFD, Color color,
			String name, List<XY_DataSet> incrFuncs, List<PlotCurveCharacterstics> incrChars,
			List<XY_DataSet> cmlFuncs, List<PlotCurveCharacterstics> cmlChars) {
		IncrementalMagFreqDist incrMFD = sol.calcNucleationMFD_forRegion(region,
				refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
		
		incrMFD.setName(name);
		incrFuncs.add(incrMFD);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
		
		EvenlyDiscretizedFunc cmlMFD = incrMFD.getCumRateDistWithOffset();
		cmlMFD.setName(name);
		cmlFuncs.add(cmlMFD);
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
	}

}
