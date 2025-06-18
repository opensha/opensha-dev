package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.Supplier;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeLevel.FileBackedLevel;
import org.opensha.commons.logicTree.LogicTreeNode.FileBackedNode;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlySampledNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SingleSolHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultCubeAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolLogicTreeSampler;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_MaxMagOffFault;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class GriddedRateDistributionSolutionWriter {

	public static void main(String[] args) throws IOException {
		File ratesDir = new File("/home/kevin/OpenSHA/nshm23/prvi/rate_raw_data");
		File invsDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/");
		FaultSystemSolution baSol = FaultSystemSolution.load(new File(invsDir,
//				"2025_01_17-prvi25_subduction_branches/results_PRVI_SUB_FM_LARGE_branch_averaged.zip"));
				"2025_01_17-prvi25_crustal_subduction_combined_branches/combined_branch_averaged_solution.zip"));
		FaultSystemSolution subLargeSol = FaultSystemSolution.load(new File(invsDir,
				"2025_01_17-prvi25_subduction_branches/results_PRVI_SUB_FM_LARGE_branch_averaged.zip"));
		FaultSystemSolution crustalBASol = FaultSystemSolution.load(new File(invsDir,
				"2025_01_17-prvi25_crustal_branches-dmSample10x/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		
//		File fullRandOutDir = new File(invsDir, "2025_05_13-prvi_gridded_variability-full_random");
		File fullRandOutDir = null;
//		File threeBranchOutDir = new File(invsDir, "2025_05_13-prvi_gridded_variability-three_branch");
//		File threeBranchOutDir = null;
		boolean threeBranchCorrelateMueAndCar = true;
		File threeBranchOutDir = new File(invsDir, "2025_05_13-prvi_gridded_variability-three_branch-corr_mue_car");
		
		Preconditions.checkState(fullRandOutDir == null || fullRandOutDir.exists() || fullRandOutDir.mkdir());
		Preconditions.checkState(threeBranchOutDir == null || threeBranchOutDir.exists() || threeBranchOutDir.mkdir());
		
		List<double[]> crustalPairs = loadRates(new File(ratesDir, "rbpairs-Crustal-Full-v3.csv"));
		List<double[]> carInterfacePairs = loadRates(new File(ratesDir, "rbpairs-CAR Interface-Full-v3.csv"));
		List<double[]> mueInterfacePairs = loadRates(new File(ratesDir, "rbpairs-MUE Interface-Full-v3.csv"));
		List<double[]> carSlabPairs = loadRates(new File(ratesDir, "rbpairs-CAR Intraslab-Full-v3.csv"));
		List<double[]> mueSlabPairs = loadRates(new File(ratesDir, "rbpairs-MUE Intraslab-Full-v3.csv"));
		int numSamplesPer = 100;
		int numTotalSamples = numSamplesPer*3*2*3*3*2; // times declustering and smoothing branches
		Preconditions.checkState(numTotalSamples <= crustalPairs.size(), "Have %s samples but want %s sampling nodes",
				crustalPairs.size(), numTotalSamples);
		System.out.println(numTotalSamples+" total samples");
//		System.exit(0);
		Random rand = new Random(numTotalSamples);
		
		List<PRVI25_DeclusteringAlgorithms> origDeclusterNodes = new ArrayList<>();
		LogicTreeLevel<PRVI25_DeclusteringAlgorithms> origDeclusterLevel = PRVI25_LogicTreeBranch.SEIS_DECLUSTER;
		List<FileBackedNode> crustalDeclusterNodes = new ArrayList<>();
		List<FileBackedNode> interfaceDeclusterNodes = new ArrayList<>();
		for (PRVI25_DeclusteringAlgorithms node : origDeclusterLevel.getNodes()) {
			if (node.getNodeWeight(null) == 0d)
				continue;
			origDeclusterNodes.add(node);
			crustalDeclusterNodes.add(new FileBackedNode("Crustal "+node.getName(), node.getShortName(),
					node.getNodeWeight(null), "Crustal"+node.getFilePrefix()));
			interfaceDeclusterNodes.add(new FileBackedNode("Interface "+node.getName(), node.getShortName(),
					node.getNodeWeight(null), "Interface"+node.getFilePrefix()));
		}
		FileBackedLevel crustalDeclusterLevel = new FileBackedLevel("Crustal "+origDeclusterLevel.getName(),
				"Crustal"+origDeclusterLevel.getShortName(), crustalDeclusterNodes);
		crustalDeclusterLevel.setAffected(origDeclusterLevel.getAffected(), origDeclusterLevel.getNotAffected(), false);
		FileBackedLevel interfaceDeclusterLevel = new FileBackedLevel("Interface "+origDeclusterLevel.getName(),
				"Interface"+origDeclusterLevel.getShortName(), interfaceDeclusterNodes);
		interfaceDeclusterLevel.setAffected(origDeclusterLevel.getAffected(), origDeclusterLevel.getNotAffected(), false);
		
		List<PRVI25_SeisSmoothingAlgorithms> origSmoothingNodes = new ArrayList<>();
		LogicTreeLevel<PRVI25_SeisSmoothingAlgorithms> origSmoothingLevel = PRVI25_LogicTreeBranch.SEIS_SMOOTH;
		List<FileBackedNode> crustalSmoothingNodes = new ArrayList<>();
		List<FileBackedNode> interfaceSmoothingNodes = new ArrayList<>();
		for (PRVI25_SeisSmoothingAlgorithms node : origSmoothingLevel.getNodes()) {
			if (node.getNodeWeight(null) == 0d)
				continue;
			origSmoothingNodes.add(node);
			crustalSmoothingNodes.add(new FileBackedNode("Crustal "+node.getName(), node.getShortName(),
					node.getNodeWeight(null), "Crustal"+node.getFilePrefix()));
			interfaceSmoothingNodes.add(new FileBackedNode("Interface "+node.getName(), node.getShortName(),
					node.getNodeWeight(null), "Interface"+node.getFilePrefix()));
		}
		FileBackedLevel crustalSmoothingLevel = new FileBackedLevel("Crustal "+origSmoothingLevel.getName(),
				"Crustal"+origSmoothingLevel.getShortName(), crustalSmoothingNodes);
		crustalSmoothingLevel.setAffected(origSmoothingLevel.getAffected(), origSmoothingLevel.getNotAffected(), false);
		FileBackedLevel interfaceSmoothingLevel = new FileBackedLevel("Interface "+origSmoothingLevel.getName(),
				"Interface"+origSmoothingLevel.getShortName(), interfaceSmoothingNodes);
		interfaceSmoothingLevel.setAffected(origSmoothingLevel.getAffected(), origSmoothingLevel.getNotAffected(), false);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> fullRandLevels = new ArrayList<>();
		fullRandLevels.add(crustalDeclusterLevel);
		fullRandLevels.add(crustalSmoothingLevel);
		fullRandLevels.add(PRVI25_LogicTreeBranch.MMAX_OFF);
		CrustalRateSamplingLevel crustalSampler = new CrustalRateSamplingLevel(crustalPairs);
		crustalSampler.buildNodes(rand, numTotalSamples);
		fullRandLevels.add(crustalSampler);
		fullRandLevels.add(interfaceDeclusterLevel);
		fullRandLevels.add(interfaceSmoothingLevel);
		fullRandLevels.add(PRVI25_LogicTreeBranch.SUB_SCALE);
		CarSlabRateSamplingLevel carSampler = new CarSlabRateSamplingLevel(carSlabPairs, carInterfacePairs);
		carSampler.buildNodes(rand, numTotalSamples);
		MueRateSamplingLevel mueSampler = new MueRateSamplingLevel(mueSlabPairs, mueInterfacePairs);
		mueSampler.buildNodes(rand, numTotalSamples);
		fullRandLevels.add(carSampler);
		fullRandLevels.add(mueSampler);
		List<LogicTreeLevel<? extends LogicTreeNode>> threeBranchLevels = new ArrayList<>();
		threeBranchLevels.add(crustalDeclusterLevel);
		threeBranchLevels.add(crustalSmoothingLevel);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.MMAX_OFF);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.CRUSTAL_SEIS_RATE);
		threeBranchLevels.add(interfaceDeclusterLevel);
		threeBranchLevels.add(interfaceSmoothingLevel);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.SUB_SCALE);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.CAR_SEIS_RATE);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.MUE_SEIS_RATE);
		
		Region reg = PRVI25_RegionLoader.loadPRVI_Tight();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.05, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		SolutionLogicTree.FileBuilder threeBranchBuilder = null;
		if (threeBranchOutDir != null) {
			baSol.write(new File(threeBranchOutDir, "full_branch_averaged.zip"));
			writeHazardScripts(threeBranchOutDir, gridReg);
			threeBranchBuilder = new SolutionLogicTree.FileBuilder(
//					new ArchiveOutput.AsynchronousZipFileOutput(new File(threeBranchOutDir, "results.zip")));
					new ArchiveOutput.ParallelZipFileOutput(new File(threeBranchOutDir, "results.zip"), 20));
			threeBranchBuilder.setSerializeGridded(true);
		}
		SolutionLogicTree.FileBuilder fullRandBuilder = null;
		if (fullRandOutDir != null) {
			baSol.write(new File(fullRandOutDir, "full_branch_averaged.zip"));
			writeHazardScripts(fullRandOutDir, gridReg);
			fullRandBuilder = new SolutionLogicTree.FileBuilder(
//					new ArchiveOutput.AsynchronousZipFileOutput(new File(fullRandOutDir, "results.zip")));
					new ArchiveOutput.ParallelZipFileOutput(new File(fullRandOutDir, "results.zip"), 20));
			fullRandBuilder.setSerializeGridded(true);
		}

		List<CrustalSamplingNode> crustalSamples = crustalSampler.getNodes();
		List<CarSlabSamplingNode> carSamples = carSampler.getNodes();
		List<MueSamplingNode> mueSamples = mueSampler.getNodes();
		
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(PRVI25_GridSourceBuilder.OVERALL_MMIN, 8.6);
		
		baSol.setVerbose(false);
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.LOGA_C4p0;
		List<LogicTreeLevel<? extends LogicTreeNode>> levelsForInterface = new ArrayList<>();
		levelsForInterface.add(PRVI25_LogicTreeBranch.SUB_SCALE);
		levelsForInterface.addAll(PRVI25_LogicTreeBranch.levelsSubductionGridded);
		
		int sampleIndex = 0;
		int threeBranchWriteCount = 0;
		int sampleWriteCount = 0;
		for (int i1=0; i1<crustalDeclusterNodes.size(); i1++) {
			FileBackedNode crustalDeclusterNode = crustalDeclusterNodes.get(i1);
			PRVI25_DeclusteringAlgorithms crustalOrigDeclusterNode = origDeclusterNodes.get(i1);
			for (int i2=0; i2<crustalSmoothingNodes.size(); i2++) {
				FileBackedNode crustalSmoothingNode = crustalSmoothingNodes.get(i2);
				PRVI25_SeisSmoothingAlgorithms crustalOrigSmoothingNode = origSmoothingNodes.get(i2);
				for (NSHM23_MaxMagOffFault mmaxOff : NSHM23_MaxMagOffFault.values()) {
					if (mmaxOff.getNodeWeight(null) == 0d)
						continue;
					for (int i3=0; i3<interfaceDeclusterNodes.size(); i3++) {
						FileBackedNode interfaceDeclusterNode = interfaceDeclusterNodes.get(i3);
						PRVI25_DeclusteringAlgorithms interfaceOrigDeclusterNode = origDeclusterNodes.get(i3);
						for (int i4=0; i4<interfaceSmoothingNodes.size(); i4++) {
							FileBackedNode interfaceSmoothingNode = interfaceSmoothingNodes.get(i4);
							PRVI25_SeisSmoothingAlgorithms interfaceOrigSmoothingNode = origSmoothingNodes.get(i4);

							if (threeBranchBuilder != null) {
								System.out.println("Building 3-branch for "+crustalOrigDeclusterNode+", "+crustalOrigSmoothingNode
										+", "+interfaceOrigDeclusterNode+", "+interfaceOrigSmoothingNode);
								List<CompletableFuture<GridResult>> futures = new ArrayList<>();
								for (PRVI25_CrustalSeismicityRate crustalRate : PRVI25_CrustalSeismicityRate.values()) {
									if (crustalRate.getNodeWeight(null) == 0d)
										continue;
									for (PRVI25_SubductionCaribbeanSeismicityRate carRate : PRVI25_SubductionCaribbeanSeismicityRate.values()) {
										if (carRate.getNodeWeight(null) == 0d)
											continue;
										for (PRVI25_SubductionMuertosSeismicityRate mueRate : PRVI25_SubductionMuertosSeismicityRate.values()) {
											if (mueRate.getNodeWeight(null) == 0d)
												continue;
											
											if (threeBranchCorrelateMueAndCar && !carRate.getShortName().equals(mueRate.getShortName()))
												continue;
											
											LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(threeBranchLevels);
											branch.setValue(crustalDeclusterNode);
											branch.setValue(crustalSmoothingNode);
											branch.setValue(mmaxOff);
											branch.setValue(crustalRate);
											branch.setValue(interfaceDeclusterNode);
											branch.setValue(interfaceSmoothingNode);
											branch.setValue(carRate);
											branch.setValue(mueRate);
											branch.setValue(scale);
											
											if (threeBranchCorrelateMueAndCar) {
												// need to correct the weight to remove the influence of the muertos
												// rate branch
												double curWeight = branch.getBranchWeight();
												// remove meurtos
												curWeight /= mueRate.getNodeWeight(branch);
												branch.setOrigBranchWeight(curWeight);
											}
											
											LogicTreeBranch<LogicTreeNode> branchForCrustal = new LogicTreeBranch<>(PRVI25_LogicTreeBranch.levelsCrustalOffFault);
											branchForCrustal.setValue(crustalOrigDeclusterNode);
											branchForCrustal.setValue(crustalOrigSmoothingNode);
											branchForCrustal.setValue(crustalRate);
											branchForCrustal.setValue(mmaxOff);
											
											LogicTreeBranch<LogicTreeNode> branchForInterface = new LogicTreeBranch<>(levelsForInterface);
											branchForInterface.setValue(interfaceOrigDeclusterNode);
											branchForInterface.setValue(interfaceOrigSmoothingNode);
											branchForInterface.setValue(mueRate);
											branchForInterface.setValue(carRate);
											branchForInterface.setValue(scale);
											
											futures.add(CompletableFuture.supplyAsync(new Supplier<GridResult>() {

												@Override
												public GridResult get() {
													try {
														GridSourceList crustal = PRVI25_GridSourceBuilder.buildCrustalGridSourceProv(crustalBASol, branchForCrustal);
														GridSourceList sub = PRVI25_GridSourceBuilder.buildCombinedSubductionGridSourceList(subLargeSol, branchForInterface);
														GridSourceList comb = GridSourceList.combine(crustal, sub);
														return new GridResult(comb, branch);
													} catch (IOException e) {
														e.printStackTrace();
														System.exit(1);
														return null;
													}
												}
											}));
										}
									}
								}
								
								for (CompletableFuture<GridResult> future : futures) {
									GridResult result = future.join();
									baSol.setGridSourceProvider(result.gridList);
									threeBranchBuilder.solution(baSol, result.branch);
									System.out.println("Writing 3-branch "+(threeBranchWriteCount++));
								}
							}
							
							if (fullRandBuilder != null) {
								System.out.println("Building "+numSamplesPer+" for "+crustalOrigDeclusterNode+", "+crustalOrigSmoothingNode
										+", "+interfaceOrigDeclusterNode+", "+interfaceOrigSmoothingNode);
								
//								CompletableFuture<Void> prevFuture = null;
								
								Stopwatch totalWatch = Stopwatch.createStarted();
								Stopwatch ioWatch = Stopwatch.createUnstarted();
								List<CompletableFuture<GridResult>> futures = new ArrayList<>();
								for (int s=0; s<numSamplesPer; s++) {
									if (s % 10 == 0) {
										double ioTime = ioWatch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
										double totTime = totalWatch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
										System.out.println("\tBuilding sample "+s+"/"+numSamplesPer+" ("+pDF.format(ioTime/totTime)+" waiting on blocking I/O)");
									}
									CrustalSamplingNode crustalSample = crustalSamples.get(sampleIndex);
									CarSlabSamplingNode carSample = carSamples.get(sampleIndex);
									MueSamplingNode mueSample = mueSamples.get(sampleIndex);
									sampleIndex++;
									LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(fullRandLevels);
									branch.setValue(crustalDeclusterNode);
									branch.setValue(crustalSmoothingNode);
									branch.setValue(crustalSample);
									branch.setValue(mmaxOff);
									branch.setValue(interfaceDeclusterNode);
									branch.setValue(interfaceSmoothingNode);
									branch.setValue(carSample);
									branch.setValue(mueSample);
									branch.setValue(scale);
									
									LogicTreeBranch<LogicTreeNode> branchForCrustal = new LogicTreeBranch<>(PRVI25_LogicTreeBranch.levelsCrustalOffFault);
									branchForCrustal.setValue(crustalOrigDeclusterNode);
									branchForCrustal.setValue(crustalOrigSmoothingNode);
									branchForCrustal.setValue(mmaxOff);
									
									LogicTreeBranch<LogicTreeNode> branchForInterface = new LogicTreeBranch<>(levelsForInterface);
									branchForInterface.setValue(interfaceOrigDeclusterNode);
									branchForInterface.setValue(interfaceOrigSmoothingNode);
									branchForInterface.setValue(scale);
									
									futures.add(CompletableFuture.supplyAsync(new Supplier<GridResult>() {

										@Override
										public GridResult get() {
											try {
												GutenbergRichterMagFreqDist carSlabGR = buildGR(carSample.slabRate, carSample.slabB,
														PRVI25_GridSourceBuilder.SLAB_MMAX, 5d, refMFD);
												GutenbergRichterMagFreqDist mueSlabGR = buildGR(mueSample.slabRate, mueSample.slabB,
														PRVI25_GridSourceBuilder.SLAB_MMAX, 5d, refMFD);
												
												GridSourceList carSlabList = PRVI25_GridSourceBuilder.buildSlabGridSourceList(
														branchForInterface, PRVI25_SeismicityRegions.CAR_INTRASLAB, carSlabGR);
												GridSourceList mueSlabList = PRVI25_GridSourceBuilder.buildSlabGridSourceList(
														branchForInterface, PRVI25_SeismicityRegions.MUE_INTRASLAB, mueSlabGR);
												
												GridSourceList combSlabList = GridSourceList.combine(carSlabList, mueSlabList);

												Function<Double, IncrementalMagFreqDist> carMFDBuilderFunc = new Function<Double, IncrementalMagFreqDist>() {
													
													@Override
													public IncrementalMagFreqDist apply(Double mmax) {
														return buildGR(carSample.interfaceRate, carSample.interfaceB,
																mmax, 5d, refMFD);
													}
												};
												Function<Double, IncrementalMagFreqDist> mueMFDBuilderFunc = new Function<Double, IncrementalMagFreqDist>() {
													
													@Override
													public IncrementalMagFreqDist apply(Double mmax) {
														return buildGR(mueSample.interfaceRate, mueSample.interfaceB,
																mmax, 5d, refMFD);
													}
												};
												
												GridSourceList carInterfaceList = PRVI25_GridSourceBuilder.buildInterfaceGridSourceList(
														subLargeSol, branchForInterface, PRVI25_SeismicityRegions.CAR_INTERFACE,
														scale.getMagAreaRelationship(), carMFDBuilderFunc);
												GridSourceList mueInterfaceList = PRVI25_GridSourceBuilder.buildInterfaceGridSourceList(
														subLargeSol, branchForInterface, PRVI25_SeismicityRegions.MUE_INTERFACE,
														scale.getMagAreaRelationship(), mueMFDBuilderFunc);
												
												GridSourceList combInterfaceList = GridSourceList.combine(carInterfaceList, mueInterfaceList);
												
												GridSourceList combSubList = GridSourceList.combine(combSlabList, combInterfaceList);
												
												GutenbergRichterMagFreqDist crustalGR = buildGR(crustalSample.rate, crustalSample.b,
														mmaxOff.getMaxMagOffFault(), 5d, refMFD);
												NSHM23_SingleRegionGridSourceProvider crustalMFD = PRVI25_GridSourceBuilder.buildCrustalGridSourceProv(
														crustalBASol, branchForCrustal, crustalBASol.getRupSet().requireModule(FaultCubeAssociations.class), crustalGR);
												GridSourceList crustalList = crustalMFD.convertToGridSourceList(5d);
												
												GridSourceList combList = GridSourceList.combine(crustalList, combSubList);
												return new GridResult(combList, branch);
											} catch (IllegalStateException | IOException e) {
												e.printStackTrace();
												System.exit(1);
												return null;
											}
										}
									}));
									
//									if (prevFuture != null) {
//										ioWatch.start();
//										prevFuture.join();
//										ioWatch.stop();
//									}
//									
//									prevFuture = CompletableFuture.runAsync(new Runnable() {
//										
//										@Override
//										public void run() {
//											baSol.setGridSourceProvider(combList);
//											try {
//												fullRandBuilder.solution(baSol, branch);
//											} catch (IOException e) {
//												e.printStackTrace();
//												System.exit(1);
//											}
//										}
//									});
								}
//								prevFuture.join();
								
								for (CompletableFuture<GridResult> future : futures) {
									GridResult result = future.join();
									baSol.setGridSourceProvider(result.gridList);
									fullRandBuilder.solution(baSol, result.branch);
									System.out.println("Writing sampled "+(sampleWriteCount++));
								}
								
								totalWatch.stop();
							}
						}
					}
				}
			}
		}
		
		if (fullRandBuilder != null)
			fullRandBuilder.close();
		if (threeBranchBuilder != null) {
			if (threeBranchCorrelateMueAndCar)
				threeBranchBuilder.setWeightProv(new BranchWeightProvider.OriginalWeights());
			threeBranchBuilder.close();
		}
		Preconditions.checkState(fullRandBuilder == null || sampleIndex == numTotalSamples,
				"Only used %s samples but expected %s", sampleIndex, numTotalSamples);
	}
	
	private static class GridResult {
		public final GridSourceList gridList;
		public final LogicTreeBranch<LogicTreeNode> branch;
		
		private GridResult(GridSourceList gridList, LogicTreeBranch<LogicTreeNode> branch) {
			super();
			this.gridList = gridList;
			this.branch = branch;
		}
	}
	
	private static GutenbergRichterMagFreqDist buildGR(double rateAboveM1, double b, double mMax, double m1, IncrementalMagFreqDist refMFD) {
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		// this sets shape, min/max
		// subtract a tiny amount from mMax so that if it's exactly at a bin edge, e.g. 7.9, it rounds down, e.g. to 7.85
		gr.setAllButTotCumRate(refMFD.getX(0), refMFD.getX(refMFD.getClosestXIndex(mMax-0.001)), 1e16, b);
		// this scales it to match
		// similarly, add a tiny amount to M1 so that if it's exactly at a bin edge (which it should be as it's determined
		// using cumulative binning), it rounds up to the incremental bin for that cumulative edge
		gr.scaleToCumRate(refMFD.getClosestXIndex(m1+0.001), rateAboveM1);
		return gr;
	}
	
	private static List<double[]> loadRates(File csvFile) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, false);
		
		boolean reading = false;
		List<double[]> ret = new ArrayList<>();
		
		MinMaxAveTracker rateTrack = new MinMaxAveTracker();
		MinMaxAveTracker bTrack = new MinMaxAveTracker();
		for (int row=0; row<csv.getNumRows(); row++) {
			if (reading) {
				double b = csv.getDouble(row, 0);
				double rate = csv.getDouble(row, 1);
				rateTrack.addValue(rate);
				bTrack.addValue(b);
				ret.add(new double[] {rate, b});
			} else if (csv.get(row, 0).startsWith("next ") && csv.get(row, 0).contains("lines are b")) {
				reading = true;
			}
		}
		System.out.println("Loaded "+ret.size()+" rate and b-value pairs from "+csvFile.getName());
		System.out.println("\trates: "+rateTrack);
		System.out.println("\tbs: "+bTrack);
		return ret;
	}
	
	@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
	@DoesNotAffect(FaultSystemSolution.RATES_FILE_NAME)
	@DoesNotAffect(GridSourceProvider.ARCHIVE_GRID_REGION_FILE_NAME)
	@DoesNotAffect(MFDGridSourceProvider.ARCHIVE_MECH_WEIGHT_FILE_NAME)
	@DoesNotAffect(GridSourceList.ARCHIVE_GRID_LOCS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_SUB_SEIS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_UNASSOCIATED_FILE_NAME)
	@Affects(GridSourceList.ARCHIVE_GRID_SOURCES_FILE_NAME)
	private static class CrustalSamplingNode implements RandomlySampledNode {

		private String name;
		private String shortName;
		private String prefix;
		private double weight;
		private long seed;
		private double rate;
		private double b;

		private CrustalSamplingNode() {
			
		}
		
		private CrustalSamplingNode(String name, String shortName, String prefix, double weight, long seed, double rate, double b) {
			super();
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
			this.rate = rate;
			this.b = b;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return weight;
		}

		@Override
		public String getFilePrefix() {
			return prefix;
		}

		@Override
		public String getShortName() {
			return shortName;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public long getSeed() {
			return seed;
		}

		@Override
		public void init(String name, String shortName, String prefix, double weight, long seed) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
		}
		
	}
	
	private static class CrustalRateSamplingLevel extends LogicTreeLevel.RandomlySampledLevel<CrustalSamplingNode> {
		
		private List<double[]> samples;
		private List<double[]> randomizedSamples;
		
		private CrustalRateSamplingLevel() {}

		public CrustalRateSamplingLevel(List<double[]> samples) {
			this.samples = samples;
		}

		@Override
		public String getShortName() {
			return "Crustal-Sampling";
		}

		@Override
		public String getName() {
			return "Crustal Rate/b Distribution Sampling";
		}

		@Override
		public CrustalSamplingNode buildNodeInstance(int index, long seed, double weight) {
			if (randomizedSamples == null) {
				randomizedSamples = new ArrayList<>(samples);
				Collections.shuffle(randomizedSamples, new Random(seed));
			}
			double[] sample = randomizedSamples.get(index);
			return new CrustalSamplingNode("Crustal Sample "+index, "Crustal-Sample"+index, "crustal_sample_"+index, weight, seed,
					sample[0], sample[1]);
		}

		@Override
		public Class<? extends CrustalSamplingNode> getType() {
			return CrustalSamplingNode.class;
		}
		
	}
	
	@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
	@DoesNotAffect(FaultSystemSolution.RATES_FILE_NAME)
	@DoesNotAffect(GridSourceProvider.ARCHIVE_GRID_REGION_FILE_NAME)
	@DoesNotAffect(MFDGridSourceProvider.ARCHIVE_MECH_WEIGHT_FILE_NAME)
	@DoesNotAffect(GridSourceList.ARCHIVE_GRID_LOCS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_SUB_SEIS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_UNASSOCIATED_FILE_NAME)
	@Affects(GridSourceList.ARCHIVE_GRID_SOURCES_FILE_NAME)
	private static class CarSlabSamplingNode implements RandomlySampledNode {

		private String name;
		private String shortName;
		private String prefix;
		private double weight;
		private long seed;
		private double slabRate;
		private double slabB;
		private double interfaceRate;
		private double interfaceB;


		private CarSlabSamplingNode() {
			
		}
		
		private CarSlabSamplingNode(String name, String shortName, String prefix, double weight, long seed,
				double slabRate, double slabB, double interfaceRate, double interfaceB) {
			super();
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
			this.slabRate = slabRate;
			this.slabB = slabB;
			this.interfaceRate = interfaceRate;
			this.interfaceB = interfaceB;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return weight;
		}

		@Override
		public String getFilePrefix() {
			return prefix;
		}

		@Override
		public String getShortName() {
			return shortName;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public long getSeed() {
			return seed;
		}

		@Override
		public void init(String name, String shortName, String prefix, double weight, long seed) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
		}
		
	}
	
	private static class CarSlabRateSamplingLevel extends LogicTreeLevel.RandomlySampledLevel<CarSlabSamplingNode> {

		private List<double[]> slabSamples;
		private List<double[]> interfaceSamples;
		private List<Integer> randomizedIndexes;
		
		private CarSlabRateSamplingLevel() {}

		public CarSlabRateSamplingLevel(List<double[]> slabSamples, List<double[]> interfaceSamples) {
			TotalRateComparator comp = new TotalRateComparator(8d);
			Collections.sort(slabSamples, comp);
			Collections.sort(interfaceSamples, comp);
			this.slabSamples = slabSamples;
			this.interfaceSamples = interfaceSamples;
		}

		@Override
		public String getShortName() {
			return "CAR-Sampling";
		}

		@Override
		public String getName() {
			return "CAR Rate/b Distribution Sampling";
		}

		@Override
		public CarSlabSamplingNode buildNodeInstance(int index, long seed, double weight) {
			if (randomizedIndexes == null) {
				randomizedIndexes = new ArrayList<>(slabSamples.size());
				for (int i=0; i<slabSamples.size(); i++)
					randomizedIndexes.add(i);
				Collections.shuffle(randomizedIndexes, new Random(seed));
			}
			int randIndex = randomizedIndexes.get(index);
			double[] slabSample = slabSamples.get(randIndex);
			double[] interfaceSample = interfaceSamples.get(randIndex);
			return new CarSlabSamplingNode("CAR Sample "+index, "CAR-Sample"+index, "car_sample_"+index, weight, seed,
					slabSample[0], slabSample[1], interfaceSample[0], interfaceSample[1]);
		}

		@Override
		public Class<? extends CarSlabSamplingNode> getType() {
			return CarSlabSamplingNode.class;
		}
		
	}
	
	@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
	@DoesNotAffect(FaultSystemSolution.RATES_FILE_NAME)
	@DoesNotAffect(GridSourceProvider.ARCHIVE_GRID_REGION_FILE_NAME)
	@DoesNotAffect(MFDGridSourceProvider.ARCHIVE_MECH_WEIGHT_FILE_NAME)
	@DoesNotAffect(GridSourceList.ARCHIVE_GRID_LOCS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_SUB_SEIS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_UNASSOCIATED_FILE_NAME)
	@Affects(GridSourceList.ARCHIVE_GRID_SOURCES_FILE_NAME)
	private static class MueSamplingNode implements RandomlySampledNode {

		private String name;
		private String shortName;
		private String prefix;
		private double weight;
		private long seed;
		private double slabRate;
		private double slabB;
		private double interfaceRate;
		private double interfaceB;

		private MueSamplingNode() {
			
		}
		
		private MueSamplingNode(String name, String shortName, String prefix, double weight, long seed,
				double slabRate, double slabB, double interfaceRate, double interfaceB) {
			super();
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
			this.slabRate = slabRate;
			this.slabB = slabB;
			this.interfaceRate = interfaceRate;
			this.interfaceB = interfaceB;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return weight;
		}

		@Override
		public String getFilePrefix() {
			return prefix;
		}

		@Override
		public String getShortName() {
			return shortName;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public long getSeed() {
			return seed;
		}

		@Override
		public void init(String name, String shortName, String prefix, double weight, long seed) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
		}
		
	}
	
	private static class MueRateSamplingLevel extends LogicTreeLevel.RandomlySampledLevel<MueSamplingNode> {

		private List<double[]> slabSamples;
		private List<double[]> interfaceSamples;
		private List<Integer> randomizedIndexes;
		
		private MueRateSamplingLevel() {}

		public MueRateSamplingLevel(List<double[]> slabSamples, List<double[]> interfaceSamples) {
			TotalRateComparator comp = new TotalRateComparator(8d);
			Collections.sort(slabSamples, comp);
			Collections.sort(interfaceSamples, comp);
			this.slabSamples = slabSamples;
			this.interfaceSamples = interfaceSamples;
		}

		@Override
		public String getShortName() {
			return "MUE-Sampling";
		}

		@Override
		public String getName() {
			return "MUE Rate/b Distribution Sampling";
		}

		@Override
		public MueSamplingNode buildNodeInstance(int index, long seed, double weight) {
			if (randomizedIndexes == null) {
				randomizedIndexes = new ArrayList<>(slabSamples.size());
				for (int i=0; i<slabSamples.size(); i++)
					randomizedIndexes.add(i);
				Collections.shuffle(randomizedIndexes, new Random(seed));
			}
			int randIndex = randomizedIndexes.get(index);
			double[] slabSample = slabSamples.get(randIndex);
			double[] interfaceSample = interfaceSamples.get(randIndex);
			return new MueSamplingNode("MUE Sample "+index, "MUE-Sample"+index, "mue_sample_"+index, weight, seed,
					slabSample[0], slabSample[1], interfaceSample[0], interfaceSample[1]);
		}

		@Override
		public Class<? extends MueSamplingNode> getType() {
			return MueSamplingNode.class;
		}
		
	}
	
	private static class TotalRateComparator implements Comparator<double[]> {
		
		private double mmax;

		public TotalRateComparator(double mmax) {
			this.mmax = mmax;
		}

		@Override
		public int compare(double[] o1, double[] o2) {
			GutenbergRichterMagFreqDist mfd1 = new GutenbergRichterMagFreqDist(5d, mmax, 10);
			mfd1.setAllButTotMoRate(mfd1.getMinX(), mfd1.getMaxX(), o1[0], o1[1]);
			GutenbergRichterMagFreqDist mfd2 = new GutenbergRichterMagFreqDist(5d, mmax, 10);
			mfd2.setAllButTotMoRate(mfd2.getMinX(), mfd2.getMaxX(), o2[0], o2[1]);
			return Double.compare(mfd1.calcSumOfY_Vals(), mfd2.calcSumOfY_Vals());
		}
		
	}
	
	private static void writeHazardScripts(File outputDir, GriddedRegion gridReg) throws IOException {
		String dirName = outputDir.getName();
		File localDir = outputDir;
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteTotalThreads = 20;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 36;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
//		JavaShellScriptWriter mpjWrite = new NoMPJSingleNodeShellScriptWriter(USC_CARC_ScriptWriter.JAVA_BIN,
//				remoteTotalMemGB*1024, null); nodes = 1; remoteInversionsPerBundle = 2;
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
		mpjWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		mpjWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		String dirPath = "$DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		if (mpjWrite instanceof NoMPJSingleNodeShellScriptWriter)
			classpath.add(new File("/project/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
		
		mpjWrite.setClasspath(classpath);
		if (mpjWrite instanceof MPJExpressShellScriptWriter)
			((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		else if (mpjWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		File gridRegFile = new File(outputDir, "gridded_region.geojson");
		Feature.write(gridReg.toFeature(), gridRegFile);
		
		Double sigmaTrunc = 3d;
		
		for (IncludeBackgroundOption bgOp : IncludeBackgroundOption.values()) {
			String mapScriptName = "batch_hazard_"+bgOp.name()+".slurm";
			String argz = "--input-file "+dirPath+"/results.zip";
			argz += " --output-dir "+dirPath+"/results";
			argz += " --output-file "+dirPath+"/results_hazard_"+bgOp.name()+".zip";
			argz += " --gridded-seis "+bgOp.name();
			if (bgOp != IncludeBackgroundOption.ONLY)
				argz += " --external-fss "+dirPath+"/full_branch_averaged.zip";
			argz += " --quick-grid-calc";
			argz += " --region "+dirPath+"/"+gridRegFile.getName();
			argz += " --gmpe "+AttenRelRef.USGS_PRVI_ACTIVE.name();
			argz += " --gmpe "+AttenRelRef.USGS_PRVI_INTERFACE.name();
			argz += " --gmpe "+AttenRelRef.USGS_PRVI_SLAB.name();
			if (sigmaTrunc != null)
				argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
			argz += " "+MPJTaskCalculator.argumentBuilder().minDispatch(1).maxDispatch(100).threads(remoteTotalThreads).build();
			List<String> script = mpjWrite.buildScript(MPJ_LogicTreeHazardCalc.class.getName(), argz);
			pbsWrite.writeScript(new File(localDir, mapScriptName), script, 1440, nodes, remoteTotalThreads, queue);
		}
	}

}
