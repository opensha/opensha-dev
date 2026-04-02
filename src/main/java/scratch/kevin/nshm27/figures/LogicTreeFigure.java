package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.SamplingMethod;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceFaultModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;

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
				
				LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
				ltFig.write(outputDir, seisReg.name()+"_"+trt.name(), true, true);
				
				boolean doSeparate = trt == TectonicRegionType.SUBDUCTION_INTERFACE
						|| (trt == TectonicRegionType.ACTIVE_SHALLOW && seisReg == NSHM27_SeismicityRegions.GNMI);
				if (doSeparate) {
					List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM27_LogicTree.buildLevels(seisReg, trt, useLevelWeights, true, false);
					tree = LogicTree.buildSampled(levels, samples, 123456l, NSHM27_InterfaceFaultModels.regionDefault(seisReg));
					
					ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
					ltFig.write(outputDir, seisReg.name()+"_"+trt.name()+"_inversion", true, true);
					
					levels = NSHM27_LogicTree.buildLevels(seisReg, trt, useLevelWeights, false, true);
					tree = LogicTree.buildSampled(levels, samples, 123456l);
					
					ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
					ltFig.write(outputDir, seisReg.name()+"_"+trt.name()+"_gridded", true, true);
				}
			}
			
			LogicTree<LogicTreeNode> multiTree = NSHM27_LogicTree.buildMultiRegimeTree(seisReg, samples, true, samplingMethod);
			LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(LogicTree.unrollTRTs(multiTree), false, useLevelWeights);
			ltFig.write(outputDir, seisReg.name()+"_combined", true, true);
		}
	}

}
