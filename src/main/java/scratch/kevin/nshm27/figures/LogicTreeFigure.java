package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.SamplingMethod;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class LogicTreeFigure {
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "logic_trees");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int samples = 10000;
		boolean useLevelWeights = true;
		
		TectonicRegionType[] trts = {TectonicRegionType.SUBDUCTION_INTERFACE, TectonicRegionType.SUBDUCTION_SLAB,
				TectonicRegionType.ACTIVE_SHALLOW};
		
		SamplingMethod samplingMethod = SamplingMethod.MONTE_CARLO;
		
		for (NSHM27_SeismicityRegions seisReg : NSHM27_SeismicityRegions.values()) {
			for (TectonicRegionType trt : trts) {
				LogicTree<LogicTreeNode> tree = NSHM27_LogicTree.buildLogicTree(seisReg, trt, samples, true, samplingMethod);
				tree = stripFaultModels(tree);
				
				LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
				ltFig.write(outputDir, seisReg.name()+"_"+trt.name(), true, true);
				
				boolean doSeparate = trt == TectonicRegionType.SUBDUCTION_INTERFACE
						|| (trt == TectonicRegionType.ACTIVE_SHALLOW && seisReg == NSHM27_SeismicityRegions.GNMI);
				if (doSeparate) {
					List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM27_LogicTree.buildLevels(seisReg, trt, useLevelWeights, true, false, false);
					levels = stripFaultModels(levels);
					tree = LogicTree.buildSampled(levels, samples, 123456l, NSHM27_InterfaceFaultModels.regionDefault(seisReg));
					
					ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
					ltFig.write(outputDir, seisReg.name()+"_"+trt.name()+"_inversion", true, true);
					
					// include common with gridded
					levels = NSHM27_LogicTree.buildLevels(seisReg, trt, useLevelWeights, false, true, true);
					levels = stripFaultModels(levels);
					tree = LogicTree.buildSampled(levels, samples, 123456l);
					
					ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
					ltFig.write(outputDir, seisReg.name()+"_"+trt.name()+"_gridded", true, true);
				}
			}
			
			LogicTree<LogicTreeNode> multiTree = NSHM27_LogicTree.buildMultiRegimeTree(seisReg, samples, true, samplingMethod);
			LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(stripFaultModels(LogicTree.unrollTRTs(multiTree)), false, useLevelWeights);
			ltFig.write(outputDir, seisReg.name()+"_combined", true, true);
		}
	}
	
	private static List<LogicTreeLevel<? extends LogicTreeNode>> stripFaultModels(List<LogicTreeLevel<? extends LogicTreeNode>> levels) {
		List<LogicTreeLevel<? extends LogicTreeNode>> ret = new ArrayList<>();
		
		for (LogicTreeLevel<? extends LogicTreeNode> level : levels) {
			if (RupSetFaultModel.class.isAssignableFrom(level.getType()))
				continue;
			ret.add(level);
		}
		
		return ret;
	}
	
	private static LogicTree<LogicTreeNode> stripFaultModels(LogicTree<LogicTreeNode> tree) {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = stripFaultModels(tree.getLevels());
		List<LogicTreeBranch<LogicTreeNode>> branches = new ArrayList<>(tree.size());
		for (LogicTreeBranch<LogicTreeNode> branch : tree) {
			List<LogicTreeNode> values = new ArrayList<>(levels.size());
			for (LogicTreeNode value : branch) {
				if (value instanceof RupSetFaultModel)
					continue;
				values.add(value);
			}
			Preconditions.checkState(values.size() == levels.size());
			LogicTreeBranch<LogicTreeNode> modBranch = new LogicTreeBranch<>(levels, values);
			modBranch.setOrigBranchWeight(branch.getOrigBranchWeight());
			branches.add(modBranch);
		}
		return LogicTree.fromExisting(levels, branches);
	}

}
