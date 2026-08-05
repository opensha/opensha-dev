package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.dmSampling.RupSetDeformationModelDistribution.BinnedUniformSamplingLevel;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceCouplingDepthModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceHingedBValue;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_InterfaceObsSeisDMAdjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

class InterfaceParticipationRateFigures {

	public static void main(String[] args) throws IOException {
		ModuleContainer.VERBOSE_DEFAULT = false;
		File outputDir = new File(FIGURES_DIR, "interface_partic");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		boolean doSLT = true;
		
		double[] minMags = {0d, 7d, 8d, 9d};
		String[] magLabels = new String[minMags.length];
		String[] magPrefixes = new String[minMags.length];
		for (int m=0; m<minMags.length; m++) {
			if (minMags[m] > 0) {
				magLabels[m] = "On-fault M>"+oDF.format(minMags[m]);
				magPrefixes[m] = "m"+oDF.format(minMags[m]);
			} else {
				magLabels[m] = "On-fault";
				magPrefixes[m] = "supra";
			}
		}
		
		double[] sltMinMags = {0d};
		String[] sltMagLabels = new String[minMags.length];
		String[] sltMagPrefixes = new String[minMags.length];
		for (int m=0; m<sltMinMags.length; m++) {
			if (sltMinMags[m] > 0) {
				sltMagLabels[m] = "On-fault M>"+oDF.format(sltMinMags[m]);
				sltMagPrefixes[m] = "m"+oDF.format(sltMinMags[m]);
			} else {
				sltMagLabels[m] = "On-fault";
				sltMagPrefixes[m] = "supra";
			}
		}
		
		List<LogicTreeNode> fixedBatchNodes = List.of(
				NSHM27_InterfaceDeformationModels.Aggregated.HIGH_COUPLING,
				NSHM27_InterfaceDeformationModels.Aggregated.LOW_COUPLING,
				NSHM27_InterfaceHingedBValue.HINGED_SINGLE_NODE,
				NSHM27_InterfaceObsSeisDMAdjustment.EXTRAPOLATE,
				NSHM27_InterfaceCouplingDepthModels.DEEP_TAPER,
				NSHM27_InterfaceCouplingDepthModels.DOUBLE_TAPER,
				NSHM27_InterfaceCouplingDepthModels.NONE);
		List<? extends Class<? extends LogicTreeLevel<?>>> fixedBatchLevelClasses = List.of(
				BinnedUniformSamplingLevel.class
				);
		
		for (NSHM27_SeismicityRegions seisReg : NSHM27_SeismicityRegions.values()) {
			CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
			if (seisReg == NSHM27_SeismicityRegions.AMSAM)
				cpt = cpt.rescale(-4, -1);
			else
				cpt = cpt.rescale(-5, -2);
			cpt.setLog10(true);
			File dir = getSolDir(seisReg);
			
			File solFile = getSolFile(seisReg);
			FaultSystemSolution baSol = FaultSystemSolution.load(solFile);
			
			List<FaultSection> interfaceSects = new ArrayList<>();
			List<Integer> interfaceIndexes = new ArrayList<>();
			for (FaultSection sect : baSol.getRupSet().getFaultSectionDataList()) {
				if (sect.getTectonicRegionType() == TectonicRegionType.SUBDUCTION_INTERFACE) {
					interfaceSects.add(sect);
					interfaceIndexes.add(sect.getSectionId());
				}
			}
			int[] interfaceRemaps = interfaceIndexes.size() == baSol.getRupSet().getNumSections() ? null : Ints.toArray(interfaceIndexes);
			
			GeographicMapMaker mapMaker = new GeographicMapMaker(interfaceSects);
			
			String baPrefix = seisReg.name()+"_ba_";
			
			for (int m=0; m<minMags.length; m++) {
				double[] rates = getRemapRates(interfaceRemaps, baSol.calcParticRateForAllSects(minMags[m], Double.POSITIVE_INFINITY));
				mapMaker.plotSectScalars(rates, cpt, magLabels[m]+" participation rate (1/yr)");
				
				mapMaker.plot(outputDir, baPrefix+magPrefixes[m], " ");
			}
			
			if (doSLT) {
				LogicTree<LogicTreeNode> analysisTree = getAnalysisLogicTree(seisReg);
				
				File sltFile = new File(dir, "results_nogrid.zip");
				SolutionLogicTree slt = SolutionLogicTree.load(sltFile);
				LogicTree<?> tree = slt.getLogicTree();
				Preconditions.checkState(tree.size() == analysisTree.size());
				
				List<LogicTreeNode> batchNodesList = new ArrayList<>();
				List<LogicTreeLevel<?>> batchNodesLevels = new ArrayList<>();
				for (LogicTreeNode node : fixedBatchNodes) {
					LogicTreeLevel<?> match = null;
					for (LogicTreeLevel<?> level : tree.getLevels()) {
						if (level.isMember(node)) {
							match = level;
							break;
						}
					}
					if (match == null) {
						// check analysis tree
						for (LogicTreeLevel<?> level : analysisTree.getLevels()) {
							if (level.isMember(node)) {
								match = level;
								break;
							}
						}
					}
					if (match != null) {
						batchNodesList.add(node);
						batchNodesLevels.add(match);
					}
				}
				for (LogicTree<?> theTree : List.of(tree, analysisTree)) {
					for (LogicTreeLevel<?> level : theTree.getLevels()) {
						for (Class<? extends LogicTreeLevel<?>> clazz : fixedBatchLevelClasses) {
							if (clazz.isAssignableFrom(level.getClass())) {
								for (LogicTreeNode node : level.getNodes()) {
									if (!batchNodesList.contains(node)) {
										batchNodesList.add(node);
										batchNodesLevels.add(level);
									}
								}
							}
						}
					}
				}
				
				LogicTreeNode[] batchNodes = new LogicTreeNode[batchNodesList.size()];
				for (int n=0; n<batchNodes.length; n++) {
					LogicTreeNode node = batchNodesList.get(n);
					LogicTreeLevel<?> level = batchNodesLevels.get(n);
					System.out.println("Will average for node: "+level.getShortName()+" "+node.getFilePrefix());
					batchNodes[n] = node;
				}
				
				double[][][] batchValues = new double[batchNodes.length][][];
				double[] batchWeightSums = new double[batchNodes.length];
				int[] batchCounts = new int[batchNodes.length];
				
				CompletableFuture<Void> processFuture = null;
				for (int b=0; b<tree.size(); b++) {
					LogicTreeBranch<?> branch = tree.getBranch(b);
					LogicTreeBranch<LogicTreeNode> analysisBranch = analysisTree.getBranch(b);
					
					System.out.println("Branch "+b+":\t"+analysisBranch);
					
					List<String> matches = new ArrayList<>();
					for (LogicTreeNode node : batchNodes)
						if (branch.hasValue(node) || analysisBranch.hasValue(node))
							matches.add(node.getFilePrefix());
					if (matches.isEmpty()) {
						System.out.println("\tNo matches, skipping");
						continue;
					}
					System.out.println("\tMatches:\t"+matches);
					
					FaultSystemSolution sol = slt.forBranch(branch);
					
					if (processFuture != null)
						processFuture.join();
					
					double weight = tree.getBranchWeight(b);
					
					processFuture = CompletableFuture.runAsync(() -> {
						int numSects = interfaceRemaps == null ? sol.getRupSet().getNumSections() : interfaceRemaps.length;
						
						double[][] magRates = new double[sltMinMags.length][];
						for (int m=0; m<sltMinMags.length; m++)
							magRates[m] = getRemapRates(interfaceRemaps, sol.calcParticRateForAllSects(sltMinMags[m], Double.POSITIVE_INFINITY));
						
						for (int n=0; n<batchNodes.length; n++) {
							LogicTreeNode node = batchNodes[n];
							if (branch.hasValue(node) || analysisBranch.hasValue(node)) {
								if (batchValues[n] == null)
									batchValues[n] = new double[sltMinMags.length][numSects];
								else
									Preconditions.checkState(batchValues[n][0].length == numSects);
								batchCounts[n]++;
								batchWeightSums[n] += weight;
								for (int m=0; m<sltMinMags.length; m++)
									for (int s=0; s<numSects; s++)
										batchValues[n][m][s] = Math.fma(magRates[m][s], weight, batchValues[n][m][s]);
							}
						}
					});
				}
				if (processFuture != null)
					processFuture.join();
				for (int n=0; n<batchNodes.length; n++) {
					System.out.println("Processing node "+batchNodes[n]+" with "+batchCounts[n]+" matches");
					if (batchCounts[n] == 0) {
						System.out.flush();
						System.err.println("\tNo matches!");
						System.err.flush();
						continue;
					}
					double weightScale = 1d/batchWeightSums[n];
					
					String nodePrefix = seisReg.name()+"_indv_"+batchNodesLevels.get(n).getFilePrefix()+"_"+batchNodes[n].getFilePrefix()+"_";
					
					for (int m=0; m<sltMinMags.length; m++) {
						double[] rates = batchValues[n][m];
						for (int s=0; s<rates.length; s++)
							rates[s] *= weightScale;
						mapMaker.plotSectScalars(rates, cpt, sltMagLabels[m]+" participation rate (1/yr)");
						
						mapMaker.plot(outputDir, nodePrefix+sltMagPrefixes[m], batchNodes[n].getName());
					}
				}
			}
		}
	}
	
	private static double[] getRemapRates(int[] interfaceRemaps, double[] rates) {
		if (interfaceRemaps == null)
			return rates;
		double[] ret = new double[interfaceRemaps.length];
		for (int i=0; i<interfaceRemaps.length; i++)
			ret[i] = rates[interfaceRemaps[i]];
		return ret;
	}

}
