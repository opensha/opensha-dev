package scratch.kevin.nshm27.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletableFuture;

import org.jfree.data.Range;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
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
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
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

		LogicTreeBranch<LogicTreeNode> branch = NSHM27_LogicTree.buildDefault(seisReg, trt, false);
		LogicTreeBranch<LogicTreeNode> sampledBranch = NSHM27_LogicTree.buildDefault(seisReg, trt, true);
//		int samples = 50;
//		int samples = 500;
		int samples = 1000;
		
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
		
		NSHM27_InvConfigFactory factory = new NSHM27_InvConfigFactory();
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(6.01, 9.99);
		
		IncrementalMagFreqDist defaultMFD = calculateMFD(factory, rupSet, branch, refMFD);
		
		CPT tab10CPT = GMT_CPT_Files.CATEGORICAL_TAB10_NOGRAY.instance();
		Color[] tab10 = new Color[tab10CPT.size()];
		for (int i=0; i<tab10.length; i++)
			tab10[i] = tab10CPT.get(i).minColor;
		
		Range magRange = new Range(6d, 9.5d);
		Range rateRange = new Range(1e-5, 1e0);
		
		List<String> avgStrings = new ArrayList<>();
		
		for (int l : levelIndexes) {
			Preconditions.checkState(l >= 0);
			LogicTreeLevel<? extends LogicTreeNode> level = branch.getLevel(l);
			LogicTreeNode defaultValue = branch.getValue(l);
			LogicTreeLevel<? extends LogicTreeNode> sampledLevel = sampledBranch.getLevel(l);
			System.out.println("Processing level "+l+". "+level.getName());
			List<DiscretizedFunc> incrFuncs = new ArrayList<>();
			List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
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
						IncrementalMagFreqDist mfd = calculateMFD(factory, rupSet, branch, refMFD);
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
				List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
				List<LogicTreeNode> values = new ArrayList<>();
				for (int i=0; i<branch.size(); i++) {
					if (i == l) {
						levels.add(sampledLevel);
						values.add(null);
					} else {
						levels.add(branch.getLevel(i));
						values.add(branch.getValue(i));
					}
				}
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
					mfdFutures.add(CompletableFuture.supplyAsync(()->{
						LogicTreeBranch<LogicTreeNode> myBranch2 = myBranch.copy();
						myBranch2.setValue(classNode);
						myBranch2.setValue(l, node);
						try {
							return calculateMFD(factory, rupSet, myBranch2, refMFD);
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
				
				IncrementalMagFreqDist baseline = defaultMFD.deepClone();
				baseline.setName(defaultValue.getName());
				incrFuncs.add(baseline);
				cmlFuncs.add(baseline.getCumRateDistWithOffset());
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
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
		}
		
		for (String str : avgStrings)
			System.out.println(str);
	}
	
	private static IncrementalMagFreqDist calculateMFD(NSHM27_InvConfigFactory factory, FaultSystemRupSet rupSet,
			LogicTreeBranch<LogicTreeNode> branch, EvenlyDiscretizedFunc refMFD) throws IOException {
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
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			gr.setAllButTotCumRate(mMin, mMax, moRate, b);
			mfd.addIncrementalMagFreqDist(gr);
		}
		
		return mfd;
	}

}
