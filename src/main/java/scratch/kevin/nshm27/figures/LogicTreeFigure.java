package scratch.kevin.nshm27.figures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM26_RegionLoader.NSHM26_SeismicityRegions;

public class LogicTreeFigure {
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/nshm26_logic_trees");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int samples = 10000;
		boolean useLevelWeights = true;
		
		TectonicRegionType[] trts = {TectonicRegionType.SUBDUCTION_INTERFACE, TectonicRegionType.SUBDUCTION_SLAB,
				TectonicRegionType.ACTIVE_SHALLOW};
		
		for (NSHM26_SeismicityRegions seisReg : NSHM26_SeismicityRegions.values()) {
			for (TectonicRegionType trt : trts) {
				LogicTree<LogicTreeNode> tree = NSHM26_LogicTree.buildLogicTree(seisReg, trt, samples, true);
				
				LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(tree, false, useLevelWeights);
				ltFig.write(outputDir, seisReg.name()+"_"+trt.name(), true, true);
			}
			
			LogicTree<LogicTreeNode> multiTree = NSHM26_LogicTree.buildMultiRegimeTree(seisReg, samples, true);
			LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(LogicTree.unrollTRTs(multiTree), false, useLevelWeights);
			ltFig.write(outputDir, seisReg.name()+"_combined", true, true);
		}
	}

}
