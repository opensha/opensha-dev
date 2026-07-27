package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletableFuture;

import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.jfree.data.Range;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.RandomLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.SamplingMethod;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.dmSampling.DeformationModelDistSampler.FixedFractileSampler;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.MaxRuptureLengthBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SectionSupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM27_InvConfigFactory;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceCouplingDepthModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceFaultModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceHingedBValue;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceHingedBValue.CombinedSampledType;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceMinSubSects;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceObsSeisDMAdjustment;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisClassificationMethod;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisRateModel;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import gov.usgs.earthquake.nshmp.erf.seismicity.SeismicityRateFileLoader.PureGR;
import net.mahdilamb.colormap.Colors;

public class InterfaceLogicTreeMFDExploration {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(NSHM27_PaperPaths.FIGURES_DIR, "interface_mfd_exploration");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		ModuleContainer.VERBOSE_DEFAULT = false;
		NSHM27_SeismicityRegions seisReg = NSHM27_SeismicityRegions.AMSAM;
		
		outputDir = new File(outputDir, seisReg.name());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		TectonicRegionType trt = TectonicRegionType.SUBDUCTION_INTERFACE;

		boolean includeObs = true;
		boolean includeTotalObs = true;

		LogicTreeBranch<LogicTreeNode> branch = NSHM27_LogicTree.buildDefault(seisReg, trt, false);
		LogicTreeBranch<LogicTreeNode> sampledBranch = NSHM27_LogicTree.buildDefault(seisReg, trt, true);
//		int samples = 50;
//		int samples = 500;
//		int samples = 1000;
		int samples = 10000;
		
//		boolean skipSampled = true;
		boolean skipSampled = false;
		
		for (int l=0; l<branch.size(); l++) {
			LogicTreeLevel<? extends LogicTreeNode> level = branch.getLevel(l);
			System.out.println(l+". "+level.getName());
			System.out.println("\tClass:\t"+level.getClass());
			LogicTreeNode value = branch.getValue(l);
			System.out.println("\tValue:\t"+value);
			System.out.println("\tValue type:\t"+level.getType());
			System.out.println("\tNum nodes:\t"+level.getNodes().size());
		}
		
		System.out.println(branch);
		
		int[] levelIndexes = {
				branch.getLevelTypeIndex(NSHM27_InterfaceCouplingDepthModels.class),
				branch.getLevelTypeIndex(NSHM27_InterfaceDeformationModels.Aggregated.class),
				branch.getLevelTypeIndex(PRVI25_SubductionScalingRelationships.class),
				branch.getLevelTypeIndex(SectionSupraSeisBValues.class),
				branch.getLevelTypeIndex(NSHM27_InterfaceObsSeisDMAdjustment.class),
				branch.getLevelTypeIndex(NSHM27_InterfaceMinSubSects.class),
				branch.getLevelTypeIndex(MaxRuptureLengthBranchNode.class),
		};
		
		WeightedList<NSHM27_SeisClassificationMethod> classificationChoices = new WeightedList<>();
		for (NSHM27_SeisClassificationMethod classification : NSHM27_SeisClassificationMethod.values())
			if (classification.getNodeWeight() > 0d)
				classificationChoices.add(classification, classification.getNodeWeight());
		List<NSHM27_SeisClassificationMethod> classificationSamples = classificationChoices.sampleEvenly(samples, new Random(12345l));
		LogicTreeLevel<?> rateLevel = sampledBranch.getLevel(sampledBranch.getLevelTypeIndex(NSHM27_SeisRateModel.class));
		((RandomLevel<?, ?>)rateLevel).build(12345l, samples, SamplingMethod.LATIN_HYPERCUBE);
		List<? extends LogicTreeNode> rateModelSamples = rateLevel.getNodes();
		
		NSHM27_InvConfigFactory factory = new NSHM27_InvConfigFactory();
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 9.99);
		
		FaultGridAssociations assoc = rupSet.requireModule(FaultGridAssociations.class);
		
		IncrementalMagFreqDist defaultMFD = calculateMFD(factory, rupSet, branch, refMFD, assoc);
		
		double defaultMmax = 0d;
		for (int m=0; m<defaultMFD.size(); m++)
			if (defaultMFD.getY(m) > 0)
				defaultMmax = defaultMFD.getX(m);
		
		CPT tab10CPT = GMT_CPT_Files.CATEGORICAL_TAB10_NOGRAY.instance();
		Color[] tab10 = new Color[tab10CPT.size()];
		for (int i=0; i<tab10.length; i++)
			tab10[i] = tab10CPT.get(i).minColor;
		
		Range magRange = new Range(6d, 9.5d);
		Range rateRange = new Range(1e-5, 1e0);
		
		List<String> avgStrings = new ArrayList<>();

		IncrementalMagFreqDist obsMFD = null;
		EvenlyDiscretizedFunc obsCmlMFD = null;
		if (includeObs) {
			NSHM27_SeisRateModel rateModel = branch.requireValue(NSHM27_SeisRateModel.class);
			obsMFD = rateModel.build(seisReg, branch.requireValue(NSHM27_SeisClassificationMethod.class), trt, refMFD, defaultMmax);
			obsMFD.setName("Observed (interface)");
			System.out.println("Observed rate M>5: "+obsMFD.getCumRate(obsMFD.getClosestXIndex(5.01)));
			System.out.println("Observed rate M>6: "+obsMFD.getCumRate(obsMFD.getClosestXIndex(6.01)));
			obsCmlMFD = obsMFD.getCumRateDistWithOffset();
			// this could be used to make it appear flat to Mmax, but it's not realistic and won't match the extrap branch
//			IncrementalMagFreqDist obsTmp = obsMFD.deepClone();
//			obsTmp.scaleToIncrRate(0, obsCmlMFD.getY(0));
//			obsCmlMFD = obsTmp;
		}
		IncrementalMagFreqDist obsTotalMFD = null;
		EvenlyDiscretizedFunc obsTotalCmlMFD = null;
		if (includeTotalObs) {
			NSHM27_SeisRateModel rateModel = branch.requireValue(NSHM27_SeisRateModel.class);
			obsTotalMFD = rateModel.build(seisReg, branch.requireValue(NSHM27_SeisClassificationMethod.class), null, refMFD, defaultMmax);
			obsTotalMFD.setName("Observed (total)");
			System.out.println("Observed total rate M>5: "+obsTotalMFD.getCumRate(obsTotalMFD.getClosestXIndex(5.01)));
			System.out.println("Observed total rate M>6: "+obsTotalMFD.getCumRate(obsTotalMFD.getClosestXIndex(6.01)));
			obsTotalCmlMFD = obsTotalMFD.getCumRateDistWithOffset();
		}
		
		for (int l : levelIndexes) {
			Preconditions.checkState(l >= 0);
			LogicTreeLevel<? extends LogicTreeNode> level = branch.getLevel(l);
			LogicTreeNode defaultValue = branch.getValue(l);
			LogicTreeLevel<? extends LogicTreeNode> sampledLevel = sampledBranch.getLevel(l);
			System.out.println("Processing level "+l+". "+level.getName());
			List<DiscretizedFunc> incrFuncs = new ArrayList<>();
			List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (includeTotalObs) {
				incrFuncs.add(obsTotalMFD);
				cmlFuncs.add(obsTotalCmlMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.GRAY));
			}
			
			if (includeObs) {
				incrFuncs.add(obsMFD);
				cmlFuncs.add(obsCmlMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, OBS_RATE_COLOR));
			}
			
			if (sampledLevel == level) {
				System.out.println("\tNormal level");
				SummedMagFreqDist avgMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				double sumWeight = 0d;
				int index = 0;
				String avgStr = level.getName()+" RIs:";
				for (LogicTreeNode node : level.getNodes()) {
					double weight = node.getNodeWeight(branch);
					if (weight > 0d) {
						branch.setValue(node);
						IncrementalMagFreqDist mfd = calculateMFD(factory, rupSet, branch, refMFD, assoc);
						avgMFD.addIncrementalMagFreqDist(mfd, weight);
						sumWeight += weight;
						Color color = tab10[index++ % tab10.length];
						mfd.setName(node.getShortName());
						incrFuncs.add(mfd);
						EvenlyDiscretizedFunc cml = mfd.getCumRateDistWithOffset();
						cmlFuncs.add(cml);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
						avgStr += "\n\t"+node.getShortName()+":\t"+(float)(1d/cml.getY(0));
					}
				}
				if ((float)sumWeight != 1f)
					avgMFD.scale(1d/sumWeight);
				avgMFD.setName("Weighted Average");
				incrFuncs.add(avgMFD);
				EvenlyDiscretizedFunc cml = avgMFD.getCumRateDistWithOffset();
				cmlFuncs.add(avgMFD.getCumRateDistWithOffset());
				avgStr += "\n\tAverage:\t"+(float)(1d/cml.getY(0));
				avgStrings.add(avgStr);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 4f, Color.DARK_GRAY));
				
				branch.setValue(defaultValue);
			} else {
//				baselineMFD = defaultMFD.deepClone();
//				baselineMFD.setName(defaultValue.getName());
				System.out.println("\tSampled version: "+sampledLevel.getName());
				if (skipSampled) {
					System.out.println("Skipping sampling this time");
					continue;
				}
				boolean doRateSamples = SectionSupraSeisBValues.class.isAssignableFrom(sampledLevel.getType());
				List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
				List<LogicTreeNode> values = new ArrayList<>();
				for (int i=0; i<branch.size(); i++) {
					if (i == l) {
						levels.add(sampledLevel);
						values.add(null);
					} else if (doRateSamples && branch.getValue(i) instanceof NSHM27_SeisRateModel) {
						levels.add(rateLevel);
						values.add(null);
					} else {
						levels.add(branch.getLevel(i));
						values.add(branch.getValue(i));
					}
				}
				Preconditions.checkNotNull(rateLevel);
				LogicTreeBranch<LogicTreeNode> myBranch = new LogicTreeBranch<>(levels, values);
				((RandomLevel<?, ?>)sampledLevel).build(12345l, samples, SamplingMethod.LATIN_HYPERCUBE);
				Preconditions.checkState(sampledLevel.getNodes().size() == samples);
				double weightEach = 1d/samples;
				SummedMagFreqDist avgMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				List<CompletableFuture<IncrementalMagFreqDist>> mfdFutures = new ArrayList<>(samples);
				List<? extends LogicTreeNode> sampledNodes = sampledLevel.getNodes();
				Preconditions.checkState(sampledNodes.size() == samples);
				for (int n=0; n<samples; n++) {
					LogicTreeNode node = sampledNodes.get(n);
					LogicTreeNode classNode = classificationSamples.get(n);
					LogicTreeNode rateNode = rateModelSamples.get(n);
					Preconditions.checkNotNull(rateNode);
					mfdFutures.add(CompletableFuture.supplyAsync(()->{
						LogicTreeBranch<LogicTreeNode> myBranch2 = myBranch.copy();
						myBranch2.setValue(classNode);
						myBranch2.setValue(l, node);
						if (doRateSamples)
							myBranch2.setValue(rateNode);
						try {
							return calculateMFD(factory, rupSet, myBranch2, refMFD, assoc);
						} catch (IOException e) {
							e.printStackTrace();
							System.err.flush();
							System.exit(1);
							return null;
						}
					}));
				}

				ArbDiscrEmpiricalDistFunc[] incrPDFs = new ArbDiscrEmpiricalDistFunc[refMFD.size()];
				ArbDiscrEmpiricalDistFunc[] cmlPDFs = new ArbDiscrEmpiricalDistFunc[refMFD.size()];
				for (int i=0; i<samples; i++) {
					CompletableFuture<IncrementalMagFreqDist> mfdFuture = mfdFutures.get(i);
					IncrementalMagFreqDist mfd = mfdFuture.join();
					EvenlyDiscretizedFunc cml = mfd.getCumRateDistWithOffset();
					avgMFD.addIncrementalMagFreqDist(mfd, weightEach);
					for (int m=0; m<mfd.size(); m++) {
						double y = mfd.getY(m);
						if (y > 0 || incrPDFs[m] != null) {
							if (incrPDFs[m] == null) {
								incrPDFs[m] = new ArbDiscrEmpiricalDistFunc();
								if (i > 0) {
									// add weight for prior ones that were zero
									incrPDFs[m].set(0d, weightEach*(i-1));
								}
							}
							incrPDFs[m].set(y, weightEach);
						}
						double cmlY = cml.getY(m);
						if (cmlY > 0 || cmlPDFs[m] != null) {
							if (cmlPDFs[m] == null) {
								cmlPDFs[m] = new ArbDiscrEmpiricalDistFunc();
								if (i > 0) {
									// add weight for prior ones that were zero
									cmlPDFs[m].set(0d, weightEach*(i-1));
								}
							}
							cmlPDFs[m].set(cmlY, weightEach);
						}
					}
				}
				double[] fractiles = {0d, 0.025, 0.16, 0.5, 0.84, 0.975, 1d};
				String fractileNames = "p[0, 2.5, 16, 84, 97.5, 100]";
//				Color transColor = new Color(0, 0, 0, 60);
				Color base = Colors.tab_blue;
				IncrementalMagFreqDist[] incrFractiles = new IncrementalMagFreqDist[fractiles.length];
				EvenlyDiscretizedFunc[] cmlFractiles = new EvenlyDiscretizedFunc[fractiles.length];
				for (int f=0; f<fractiles.length; f++) {
					incrFractiles[f] = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
					cmlFractiles[f] = incrFractiles[f].getCumRateDistWithOffset();
					for (int m=0; m<refMFD.size(); m++) {
						if (incrPDFs[m] != null)
							incrFractiles[f].set(m, incrPDFs[m].getInterpolatedFractile(fractiles[f]));
						if (cmlPDFs[m] != null)
							cmlFractiles[f].set(m, cmlPDFs[m].getInterpolatedFractile(fractiles[f]));
					}
				}
				IncrementalMagFreqDist incrMedian = incrFractiles[3];
				EvenlyDiscretizedFunc cmlMedian = cmlFractiles[3];
				
				if (NSHM27_InterfaceDeformationModels.class.isAssignableFrom(sampledLevel.getType())) {
					// expand extrema to actual bounds
					FixedFractileSampler fixedMin = new FixedFractileSampler(NSHM27_InterfaceDeformationModels.MIN_DM_FRACTILE);
					LogicTreeBranch<LogicTreeNode> myBranch2 = myBranch.copy();
					myBranch2.setValue(l, new NSHM27_InterfaceDeformationModels("Lower", "Lower", "Lower", 0d, fixedMin));
					incrFractiles[0] = calculateMFD(factory, rupSet, myBranch2, refMFD, assoc);
					cmlFractiles[0] = incrFractiles[0].getCumRateDistWithOffset();
					
					FixedFractileSampler fixedMax = new FixedFractileSampler(1d);
					myBranch2.setValue(l, new NSHM27_InterfaceDeformationModels("Upper", "Upper", "Upper", 0d, fixedMax));
					incrFractiles[6] = calculateMFD(factory, rupSet, myBranch2, refMFD, assoc);
					cmlFractiles[6] = incrFractiles[6].getCumRateDistWithOffset();
				}
				
				UncertainBoundedIncrMagFreqDist incrExtrema = new UncertainBoundedIncrMagFreqDist(
						incrMedian, incrFractiles[0], incrFractiles[6], null);
				UncertainArbDiscFunc cmlExtrema = new UncertainArbDiscFunc(cmlMedian, cmlFractiles[0], cmlFractiles[6]);
				incrExtrema.setName(null);
				cmlExtrema.setName(null);
				incrFuncs.add(incrExtrema);
				cmlFuncs.add(cmlExtrema);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(base.getRed(), base.getGreen(), base.getBlue(), 50)));
				
				UncertainBoundedIncrMagFreqDist incr95 = new UncertainBoundedIncrMagFreqDist(
						incrMedian, incrFractiles[1], incrFractiles[5], null);
				UncertainArbDiscFunc cml95 = new UncertainArbDiscFunc(cmlMedian, cmlFractiles[1], cmlFractiles[5]);
				incr95.setName(null);
				cml95.setName(null);
				incrFuncs.add(incr95);
				cmlFuncs.add(cml95);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(base.getRed(), base.getGreen(), base.getBlue(), 70)));
				
				UncertainBoundedIncrMagFreqDist incr68 = new UncertainBoundedIncrMagFreqDist(
						incrMedian, incrFractiles[2], incrFractiles[4], null);
				UncertainArbDiscFunc cml68 = new UncertainArbDiscFunc(cmlMedian, cmlFractiles[2], cmlFractiles[4]);
				incr68.setName(fractileNames);
				cml68.setName(fractileNames);
				incrFuncs.add(incr68);
				cmlFuncs.add(cml68);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(base.getRed(), base.getGreen(), base.getBlue(), 100)));
				
				if (sampledLevel instanceof MaxRuptureLengthBranchNode.DistributionSamplingLevel) {
					MaxRuptureLengthBranchNode.DistributionSamplingLevel distSampleLevel =
							(MaxRuptureLengthBranchNode.DistributionSamplingLevel)sampledLevel;
					MaxRuptureLengthBranchNode.FixedValueLevel fixedValueLevel = (MaxRuptureLengthBranchNode.FixedValueLevel)level;
					ContinuousDistribution dist = distSampleLevel.getDistribution();
					double origFixedValue = fixedValueLevel.getValue();
					
					fixedValueLevel.setValue(dist.getSupportLowerBound());
					branch.setValue(l, fixedValueLevel.getNodes().get(0));
					IncrementalMagFreqDist mfdLow = calculateMFD(factory, rupSet, branch, refMFD, assoc);
					mfdLow.setName("Lmax="+(int)dist.getSupportLowerBound()+" km");
					incrFuncs.add(mfdLow);
					cmlFuncs.add(mfdLow.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_green));
					
					IncrementalMagFreqDist baseline = defaultMFD.deepClone();
					baseline.setName("Lmax="+((MaxRuptureLengthBranchNode)defaultValue).getValue().intValue()+" km");
					incrFuncs.add(baseline);
					cmlFuncs.add(baseline.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
					
					fixedValueLevel.setValue(dist.getSupportUpperBound());
					branch.setValue(l, fixedValueLevel.getNodes().get(0));
					IncrementalMagFreqDist mfdHigh = calculateMFD(factory, rupSet, branch, refMFD, assoc);
					mfdHigh.setName("Lmax="+(int)dist.getSupportUpperBound()+" km");
					incrFuncs.add(mfdHigh);
					cmlFuncs.add(mfdHigh.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_orange));
					
					// again on top
					baseline = baseline.deepClone();
					baseline.setName(null);
					incrFuncs.add(baseline);
					cmlFuncs.add(baseline.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
					
					fixedValueLevel.setValue(origFixedValue);
					branch.setValue(l, fixedValueLevel.getNodes().get(0));
				} else {
					IncrementalMagFreqDist baseline = defaultMFD.deepClone();
					baseline.setName(defaultValue.getName());
					incrFuncs.add(baseline);
					cmlFuncs.add(baseline.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				}
				
				if (SectionSupraSeisBValues.class.isAssignableFrom(level.getType())) {
					SummedMagFreqDist avgHingedMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
					SummedMagFreqDist avgSampledMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
					double sumHingedWeight = 0d;
					double sumSampledWeight = 0d;
					for (int n=0; n<samples; n++) {
						LogicTreeNode node = sampledNodes.get(n);
						Preconditions.checkState(node instanceof CombinedSampledType);
						CombinedSampledType bNode = (CombinedSampledType)node;
						IncrementalMagFreqDist mfd = mfdFutures.get(n).join();
//						System.out.println("b node "+n+": "+node+" ("+node.getClass()+")");
						if (bNode.isHinged()) {
							sumHingedWeight += weightEach;
							avgHingedMFD.addIncrementalMagFreqDist(mfd, weightEach);
						} else {
							sumSampledWeight += weightEach;
							avgSampledMFD.addIncrementalMagFreqDist(mfd, weightEach);
						}
					}
					avgHingedMFD.scale(1d/sumHingedWeight);
					avgSampledMFD.scale(1d/sumSampledWeight);
					
					avgHingedMFD.setName("Average hinged b");
					incrFuncs.add(avgHingedMFD);
					cmlFuncs.add(avgHingedMFD.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_orange));
					
					avgSampledMFD.setName("Average sampled b");
					incrFuncs.add(avgSampledMFD);
					cmlFuncs.add(avgSampledMFD.getCumRateDistWithOffset());
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_green));
					
					System.out.println("sumHingedWeight="+sumHingedWeight);
					System.out.println("sumSampledWeight="+sumSampledWeight);
//					System.exit(0);
				}
				
				avgMFD.setName("Average");
				incrFuncs.add(avgMFD);
				cmlFuncs.add(avgMFD.getCumRateDistWithOffset());
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 4f, Color.DARK_GRAY));
			}
			
			String prefix = l+"_"+level.getFilePrefix();
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			PlotSpec incrPlot = new PlotSpec(incrFuncs, chars, level.getName(), "Magnitude", "Incremental Rate (1/yr)");
			incrPlot.setLegendInset(true);
			gp.drawGraphPanel(incrPlot, false, true, magRange, rateRange);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 800, true, true, false);
			
			PlotSpec cmlPlot = new PlotSpec(cmlFuncs, chars, level.getName(), "Magnitude", "Cumulative Rate (1/yr)");
			cmlPlot.setLegendInset(true);
			gp.drawGraphPanel(cmlPlot, false, true, magRange, rateRange);
			
			PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 800, true, true, true);
			
			if (includeObs && level.getType().isAssignableFrom(NSHM27_InterfaceMinSubSects.class)) {
				int indexOffset = includeTotalObs ? 2 : 1;
				// do special version with sub-seis to mMin lines
				List<DiscretizedFunc> origIncrs = new ArrayList<>();
				List<? extends LogicTreeNode> nodes = level.getNodes();
				for (int i=0; i<nodes.size(); i++)
					origIncrs.add(incrFuncs.get(i+indexOffset));
				List<IncrementalMagFreqDist> addFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> addChars = new ArrayList<>();
				CPT tab10light = GMT_CPT_Files.CATEGORICAL_TAB10_LIGHT_NOGRAY.instance();
				for (int i=0; i<origIncrs.size(); i++) {
					Color color = tab10light.get(i).minColor;
					DiscretizedFunc fault = origIncrs.get(i);
					IncrementalMagFreqDist subSeis = obsMFD.deepClone();
					int faultMminIndex = 0;
					for (faultMminIndex=0; faultMminIndex<fault.size(); faultMminIndex++)
						if (fault.getY(faultMminIndex) > 0)
							break;
					for (int m=faultMminIndex+1; m<subSeis.size(); m++)
						subSeis.set(m, 0d);
					subSeis.setName(null);
					addFuncs.add(subSeis);
					addChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
					DiscretizedFunc faultClone = new ArbitrarilyDiscretizedFunc();
					// prevent the down spikes
					for (int m=0; m<faultMminIndex; m++)
						faultClone.set(fault.getX(m), Double.NaN);
					// stitch it to fault
					faultClone.set(fault.getX(faultMminIndex-1), subSeis.getY(faultMminIndex-1));
					for (int m=faultMminIndex; m<fault.size(); m++)
						faultClone.set(fault.getX(m), fault.getY(m));
					faultClone.setName(fault.getName());
					incrFuncs.set(i+indexOffset, faultClone);
					// prevent sub-seis down spike
					for (int m=faultMminIndex; m<subSeis.size(); m++)
						subSeis.set(m, Double.NaN);
				}
				Collections.reverse(addFuncs);
				Collections.reverse(addChars);
				incrFuncs.addAll(indexOffset, addFuncs);
				chars.addAll(indexOffset, addChars);
				// remove average
				incrFuncs.remove(incrFuncs.size()-1);
				chars.remove(chars.size()-1);
				incrPlot = new PlotSpec(incrFuncs, chars, level.getName(), "Magnitude", "Incremental Rate (1/yr)");
				incrPlot.setLegendInset(true);
				gp.drawGraphPanel(incrPlot, false, true, magRange, rateRange);
				
				PlotUtils.writePlots(outputDir, prefix+"_sub_seis", gp, 800, 800, true, true, false);
			}
		}
		
		for (String str : avgStrings)
			System.out.println(str);
	}
	
	private static IncrementalMagFreqDist calculateMFD(NSHM27_InvConfigFactory factory, FaultSystemRupSet rupSet,
			LogicTreeBranch<LogicTreeNode> branch, EvenlyDiscretizedFunc refMFD, FaultGridAssociations assoc) throws IOException {
		System.out.println("Calculating MFD for: "+branch);
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		rupSet = factory.updateRuptureSetForBranch(rupSet, branch);
		BinaryRuptureProbabilityCalc exclusionModel = NSHM27_InvConfigFactory.getExclusionModel(rupSet, branch, cRups);
		
		BitSet includedRups = new BitSet(rupSet.getNumRuptures());
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++)
			if (exclusionModel.isRupAllowed(cRups.get(rupIndex), false))
				includedRups.set(rupIndex);
		
		// calculate average mMin and mMax for this simple demonstration
		MinMaxAveTracker mMinTrack = new MinMaxAveTracker();
		MinMaxAveTracker mMaxTrack = new MinMaxAveTracker();
		for (int s=0; s<rupSet.getNumSections(); s++) {
			double mMin = Double.POSITIVE_INFINITY;
			double mMax = Double.NEGATIVE_INFINITY;
			for (int rupIndex : rupSet.getRupturesForSection(s)) {
				if (includedRups.get(rupIndex)) {
					double mag = rupSet.getMagForRup(rupIndex);
					mMin = Math.min(mag, mMin);
					mMax = Math.max(mag, mMax);
				}
			}
			mMinTrack.addValue(mMin);
			mMaxTrack.addValue(mMax);
		}
		// take average and snap to binning
		double mMin = refMFD.getX(refMFD.getClosestXIndex(mMinTrack.getAverage()));
		double mMax = refMFD.getX(refMFD.getClosestXIndex(mMaxTrack.getAverage()));
		
		double b = branch.requireValue(SectionSupraSeisBValues.class).getB(rupSet, branch);
		if (branch.hasValue(NSHM27_InterfaceObsSeisDMAdjustment.EXTRAPOLATE))
			b = ((PureGR)branch.requireValue(NSHM27_SeisRateModel.class).getRateRecord(
					branch.requireValue(NSHM27_InterfaceFaultModels.class).getSeismicityRegion(),
					branch.requireValue(NSHM27_SeisClassificationMethod.class), TectonicRegionType.SUBDUCTION_INTERFACE)).b;
		Preconditions.checkState(Double.isFinite(b));
		
		SectSlipRates slipRates = rupSet.requireModule(SectSlipRates.class);
		
		SummedMagFreqDist mfd = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int s=0; s<rupSet.getNumSections(); s++) {
			double moRate = slipRates.calcMomentRate(s);
			moRate *= assoc.getSectionFractInRegion(s);
			if (moRate > 0d) {
				GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				gr.setAllButTotCumRate(mMin, mMax, moRate, b);
				mfd.addIncrementalMagFreqDist(gr);
			}
		}
		
		return mfd;
	}

}
