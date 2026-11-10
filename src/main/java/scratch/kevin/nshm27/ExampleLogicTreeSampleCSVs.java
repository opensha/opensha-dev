package scratch.kevin.nshm27;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.ValuedLogicTreeNode;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.sectDistSampling.SectDistributionSampler.FixedFractileSampler;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SectionSupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceHingedBValue;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceHingedBValue.CombinedSampledType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_SeisRateModel.ClassificationDependentGR;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import scratch.kevin.nshm27.figures.LogicTreeFigure;

public class ExampleLogicTreeSampleCSVs {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/OpenSHA/nshm27/sampling/example_trees");
		NSHM27_SeismicityRegions seisReg = NSHM27_SeismicityRegions.AMSAM;
		String prefix = "nshm27_amsam_interface";
		
		boolean inversion = true;
		boolean gridded = false;
		
//		boolean inversion = false;
//		boolean gridded = true;
		
//		boolean inversion = true;
//		boolean gridded = true;
		
		int[] sampleCounts = {
				256,
				512,
				1024,
				2048,
				4096,
				8192,
				16384
		};
		long seed = 123456l;
		SamplingMethod method = SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE;
		
		if (inversion && gridded)
			prefix += "_full";
		else if (inversion)
			prefix += "_inversion_only";
		else if (gridded)
			prefix += "_gridded_only";
		else
			throw new IllegalStateException();
		
		boolean common = gridded;
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM27_LogicTree.buildLevels(
				seisReg, TectonicRegionType.SUBDUCTION_INTERFACE,
				true, inversion, gridded, common);
		
		for (int s=0; s<sampleCounts.length; s++) {
			int samples = sampleCounts[s];
			System.out.println("=========================");
			System.out.println("Doing "+samples+" samples");
			String samplePrefix = prefix+"_"+samples+"samples_"+method.name();
			
			LogicTree<LogicTreeNode> logicTree = LogicTree.buildSampled(levels, samples, seed, method,
					NSHM27_InterfaceFaultModels.regionDefault(seisReg));

			// remove the unnecessary fault model level (it's fixed)
			logicTree = LogicTreeFigure.stripFaultModels(logicTree);
			// remove the unnecessary overall model/regime level (it's fixed)
			logicTree = LogicTreeFigure.stripModelLevel(logicTree);
			
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>();
			for (LogicTreeLevel<?> level : logicTree.getLevels())
				header.add(level.getName());
			csv.addLine(header);
			
			for (LogicTreeBranch<LogicTreeNode> branch : logicTree) {
				List<String> line = new ArrayList<>(header.size());
				
				for (LogicTreeNode node : branch) {
					if (node instanceof NSHM27_InterfaceHingedBValue.CombinedSampledType combBSample) {
						if (combBSample.isHinged())
							line.add(NSHM27_InterfaceHingedBValue.SHORT_NAME);
						else
							line.add((float)combBSample.getB(null, null)+"");
					} else if (node instanceof ValuedLogicTreeNode<?> valued) {
						Object value = valued.getValue();
						if (value instanceof FixedFractileSampler fractiles) {
							line.add((float)fractiles.getFixedFractile()+"");
						} else if (value instanceof ClassificationDependentGR gr) {
							line.add((float)gr.getSampleFractile()+"");
						} else {
//							if (!(value instanceof Number))
//								System.out.println("UNKNOWN value "+value+" of type "+value.getClass().getName()
//										+" with node "+node.getShortName()+" of type "+node.getClass());
							line.add(value.toString());
						}
					} else {
						line.add(node.getShortName());
					}
				}
				
				csv.addLine(line);
			}
			
			csv.writeToFile(new File(outputDir, samplePrefix+".csv"));
			
			if (s == 0) {
				LogicTreeFigureWriter ltFig = new LogicTreeFigureWriter(logicTree, false, true);
				ltFig.write(outputDir, prefix, true, true);
			}
			System.out.println("=========================");
			System.out.println();
		}
	}

}
