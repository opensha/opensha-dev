package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

public class SubNodeBAFigures {

	public static void main(String[] args) throws IOException {
		File resultsZip = new File(SUBDUCTION_DIR, "results.zip");
		File outputDir = new File(FIGURES_DIR, "sub_sol");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsZip);
		
		LogicTree<?> tree = slt.getLogicTree();
		
		List<IncrementalMagFreqDist> carBranchMFDs = new ArrayList<>(tree.size());
		List<IncrementalMagFreqDist> mueBranchMFDs = new ArrayList<>(tree.size());
		List<IncrementalMagFreqDist> totBranchMFDs = new ArrayList<>(tree.size());
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 9.49);
		Region carReg = PRVI25_SeismicityRegions.CAR_INTERFACE.load();
		Region mueReg = PRVI25_SeismicityRegions.MUE_INTERFACE.load();
		
		for (LogicTreeBranch<?> branch : tree) {
			FaultSystemSolution sol = slt.forBranch(branch);
			
			totBranchMFDs.add(sol.calcNucleationMFD_forRegion(
					null, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false));
			carBranchMFDs.add(sol.calcNucleationMFD_forRegion(carReg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false));
			mueBranchMFDs.add(sol.calcNucleationMFD_forRegion(mueReg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false));
		}
		
		for (LogicTreeLevel<?> level : tree.getLevels()) {
			String prefix = "branch_mfds_"+level.getFilePrefix();
			plotLevelMFDs(tree, totBranchMFDs, level, outputDir, prefix, level.getName());
			plotLevelMFDs(tree, carBranchMFDs, level, outputDir, prefix+"_car", level.getName()+" (Caribbean Trench)");
			plotLevelMFDs(tree, mueBranchMFDs, level, outputDir, prefix+"_mue", level.getName()+" (Muertos Trough)");
		}
	}
	
	private static void plotLevelMFDs(LogicTree<?> tree, List<IncrementalMagFreqDist> mfds, LogicTreeLevel<?> level,
			File outputDir, String prefix, String title) throws IOException {
		Range xRange = new Range(7d, 9.5d);
		Range yRange = new Range(1e-6, 1e-3);
		Range cmlYRange = new Range(1e-6, 1e-2);
		
		IncrementalMagFreqDist totMFD = mfdForBranchNode(tree, mfds, null);
		
		List<IncrementalMagFreqDist> funcs = new ArrayList<>();
		for (LogicTreeNode node : level.getNodes()) {
			IncrementalMagFreqDist choiceMFD = mfdForBranchNode(tree, mfds, node);
			if (choiceMFD == null)
				continue;
			choiceMFD.setName(node.getShortName());
			funcs.add(choiceMFD);
		}
		if (funcs.size() < 2)
			return;
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0, funcs.size()-1);
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		for (int i=0; i<funcs.size(); i++)
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt.getColor((float)i)));
		
		totMFD.setName("Branch-averaged");
		funcs.add(totMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 5f, Color.BLACK));
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", "Incremental Rate (1/yr)");
		plot.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
		
		List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>(funcs.size());
		for (IncrementalMagFreqDist mfd : funcs)
			cmlFuncs.add(mfd.getCumRateDistWithOffset());
		
		plot = new PlotSpec(cmlFuncs, chars, title, "Magnitude", "Cumulative Rate (1/yr)");
		plot.setLegendInset(true);
		
		gp.drawGraphPanel(plot, false, true, xRange, cmlYRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, true, false);
	}
	
	private static IncrementalMagFreqDist mfdForBranchNode(LogicTree<?> tree, List<IncrementalMagFreqDist> mfds, LogicTreeNode node) {
		double totWeight = 0d;
		SummedMagFreqDist totMFD = null;
		
		for (int b=0; b<mfds.size(); b++) {
			LogicTreeBranch<?> branch = tree.getBranch(b);
			if (node != null && !branch.hasValue(node))
				continue;
			double weight = tree.getBranchWeight(b);
			totWeight += weight;
			IncrementalMagFreqDist mfd = mfds.get(b);
			if (totMFD == null)
				totMFD = new SummedMagFreqDist(mfd.getMinX(), mfd.getMaxX(), mfd.size());
			totMFD.addIncrementalMagFreqDist(mfd, weight);
		}
		if (totMFD != null && (float)totWeight != 1f)
			totMFD.scale(1d/totWeight);
		return totMFD;
	}

}
