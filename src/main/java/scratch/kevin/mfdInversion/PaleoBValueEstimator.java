package scratch.kevin.mfdInversion;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.stream.Collectors;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.UniformContinuousDistribution;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.WeightedContinuousDistribution;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscrFuncContinuousDistribution;
import org.opensha.commons.data.function.EvenlyDiscrFuncContinuousDistribution.DiscretizationType;
import org.opensha.commons.data.uncertainty.BoundedUncertainty;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.commons.util.Interpolate;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBySectDetailPlots;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBySectDetailPlots.AlongStrikePlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoDataLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoProbabilityModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SectionSupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class PaleoBValueEstimator {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/paleo_b_value");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
		branch.setValue(NSHM23_DeformationModels.AVERAGE);
//		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		System.out.println(branch);
		
		NSHM23_InvConfigFactory.APPLY_DEF_MODEL_UNCERTAINTIES_DEFAULT = false;
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		FaultSystemRupSet rs = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
		// if true, sections use weighted combinations of distributions
		// if false, sections use weighted combinations of likelihoods
		boolean sectUseWeightedDistribution = false;
		
		ClusterRuptures cRups = rs.requireModule(ClusterRuptures.class);
		BinaryRuptureProbabilityCalc exclusionModel = NSHM23_InvConfigFactory.getExclusionModel(rs, branch, cRups);
		BitSet includedRups = exclusionModel == null ? null : new BitSet(rs.getNumRuptures());
		if (exclusionModel != null) {
			for (int r=0; r<rs.getNumRuptures(); r++)
				if (exclusionModel.isRupAllowed(cRups.get(r), false))
					includedRups.set(r);
		}
		
		ContinuousDistribution priorDist = UniformContinuousDistribution.of(0d, 1d);
//		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 11);
		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 21);
		double connectivityWeightB = priorDist.getMean();
//		double connectivityWeightB = 1d;
		
		PaleoseismicConstraintData paleoData = NSHM23_PaleoDataLoader.load(rs);
		
		List<? extends SectMappedUncertainDataConstraint> paleoConstraints =
				paleoData.getPaleoRateConstraints();
		List<EvenlyDiscretizedFunc> paleoBValMisfits = new ArrayList<>(paleoConstraints.size());
		List<EvenlyDiscretizedFunc> paleoBValEstRates = new ArrayList<>(paleoConstraints.size());
		for (int s=0; s<paleoConstraints.size(); s++) {
			paleoBValMisfits.add(bVals.deepClone());
			paleoBValEstRates.add(bVals.deepClone());
		}
		
		// calculate for each b-value
		List<CompletableFuture<InversionTargetMFDs>> targetFutures = new ArrayList<>();
		for (int i=0; i<bVals.size(); i++) {
			double b = bVals.getX(i);
			targetFutures.add(CompletableFuture.supplyAsync(()->buildTargetsForB(factory, rs, branch, b)));
		}
		InversionTargetMFDs connectivityWeightTargetMFDs = buildTargetsForB(factory, rs, branch, connectivityWeightB);
		InversionTargetMFDs[] bValTargetMFDs = new InversionTargetMFDs[bVals.size()];
		for (int i=0; i<bValTargetMFDs.length; i++)
			bValTargetMFDs[i] = targetFutures.get(i).join();
		System.out.println("DONE calculating for "+bVals.size()+" b-values");
		
		NSHM23_PaleoProbabilityModel paleoProb = new NSHM23_PaleoProbabilityModel();
		
		for (int i=0; i<bValTargetMFDs.length; i++) {
			rs.removeModuleInstances(InversionTargetMFDs.class);
			double b = bVals.getX(i);
			System.out.println("====================================");
			System.out.println("Site results for b="+(float)b);
			for (int s=0; s<paleoConstraints.size(); s++) {
				SectMappedUncertainDataConstraint paleoConstr = paleoConstraints.get(s);
				IncrementalMagFreqDist mfd = bValTargetMFDs[i].
						getOnFaultSupraSeisNucleationMFDs().get(paleoConstr.sectionIndex);
				double[] sumParticScalars = new double[mfd.size()];
				double[] sumPaleoScalars = new double[mfd.size()];
				int[] binCounts = new int[mfd.size()];
				double sectArea = rs.getAreaForSection(paleoConstr.sectionIndex);
				for (int rupIndex : rs.getRupturesForSection(paleoConstr.sectionIndex)) {
					int magIndex = mfd.getClosestXIndex(rs.getMagForRup(rupIndex));
					binCounts[magIndex]++;
					double rupArea = rs.getAreaForRup(rupIndex);
					sumParticScalars[magIndex] += rupArea/sectArea;
					double paleoVisibleProb =  paleoProb.getProbPaleoVisible(rs, rupIndex, paleoConstr.sectionIndex);
					sumPaleoScalars[magIndex] += paleoVisibleProb;
				}
				double paleoVisibleRate = 0d;
				for (int m=0; m<mfd.size(); m++) {
					double totNuclRate = mfd.getY(m);
					if (totNuclRate > 0d) {
						double particRate = totNuclRate * sumParticScalars[m]/binCounts[m];
						paleoVisibleRate += particRate * sumPaleoScalars[m]/binCounts[m];
					}
				}
//				double misfit = (paleoVisibleRate - paleoConstr.bestEstimate)/paleoConstr.getPreferredStdDev();
				BoundedUncertainty uncertainty = (BoundedUncertainty)paleoConstr.uncertainties[0];
				double misfit = paleoConstr.estimateDataZ(paleoVisibleRate);
				System.out.println(paleoConstr.name+" ("+paleoConstr.sectionName+"):");
				System.out.println("\tpaleoRate="+(float)paleoConstr.bestEstimate+"\testSolRate="+(float)paleoVisibleRate+"\tmisfit="+(float)misfit);
				System.out.println("\tpaleoRI="+(float)(1d/paleoConstr.bestEstimate)+"\testSolRI="+(float)(1d/paleoVisibleRate));
				paleoBValMisfits.get(s).set(i, misfit);
				paleoBValEstRates.get(s).set(i, paleoVisibleRate);
			}
			System.out.println("====================================");
		}
		
		DecimalFormat df = new DecimalFormat("0.000");
		NamedFaults faults = rs.requireModule(NamedFaults.class);
		Map<String, List<EvenlyDiscretizedFunc>> namedFaultResults = new HashMap<>();
		List<EvenlyDiscretizedFunc> ssafResults = new ArrayList<>();
		List<EvenlyDiscretizedFunc> nsafResults = new ArrayList<>();
		namedFaultResults.put("San Andreas (Southern)", nsafResults);
		namedFaultResults.put("San Andreas (Northern)", ssafResults);
		for (int s=0; s<paleoConstraints.size(); s++) {
			SectMappedUncertainDataConstraint paleoConstr = paleoConstraints.get(s);
			System.out.println(paleoConstr.name+" ("+paleoConstr.sectionName+"):");
			EvenlyDiscretizedFunc result = paleoBValMisfits.get(s);
			for (int b=0; b<bVals.size(); b++)
				System.out.print("\t"+df.format(result.getY(b)));
			System.out.println();
			System.out.println();
			int parentID = rs.getFaultSectionData(paleoConstr.sectionIndex).getParentSectionId();
			String faultName = faults.getFaultName(parentID);
			if (faultName != null) {
				List<EvenlyDiscretizedFunc> faultResults =  namedFaultResults.get(faultName);
				if (faultResults == null) {
					faultResults = new ArrayList<>();
					namedFaultResults.put(faultName, faultResults);
				}
				faultResults.add(result);
				if (faultName.toLowerCase().startsWith("san andreas")) {
					if (paleoConstr.dataLocation.lat > 36.5)
						nsafResults.add(result);
					else
						ssafResults.add(result);
				}
			}
		}
		System.out.println("=====================");
		System.out.println("Named Faults:");
		System.out.println("=====================");
		List<String> faultNames = new ArrayList<>(namedFaultResults.keySet());
		Collections.sort(faultNames);
		for (String faultName : faultNames) {
			List<EvenlyDiscretizedFunc> faultResults = namedFaultResults.get(faultName);
			EvenlyDiscretizedFunc avgResult = new EvenlyDiscretizedFunc(bVals.getMinX(), bVals.getMaxX(), bVals.size());
			for (EvenlyDiscretizedFunc result : faultResults)
				for (int i=0; i<result.size(); i++)
					avgResult.add(i, result.getY(i));
			avgResult.scale(1d/faultResults.size());
			System.out.println(faultName+":");
			for (int b=0; b<bVals.size(); b++)
				System.out.print("\t"+df.format(avgResult.getY(b)));
			System.out.println();
			System.out.println();
		}
		
		// now build posteriori distributions
		System.out.println("Building site posteriors");
		Map<Integer, List<Integer>> parentPaleoIndexes = new HashMap<>();
		Map<String, List<Integer>> faultPaleoIndexes = new HashMap<>();
		List<ContinuousDistribution> sitePosteriorDists = new ArrayList<>(paleoConstraints.size());
		List<double[]> siteLogLikelihoods = new ArrayList<>(paleoConstraints.size());
		int numRups = rs.getNumRuptures();
		BitSet[] siteRups = new BitSet[paleoConstraints.size()];
		for (int s=0; s<paleoConstraints.size(); s++) {
			SectMappedUncertainDataConstraint paleoConstr = paleoConstraints.get(s);
			int parentID = rs.getFaultSectionData(paleoConstr.sectionIndex).getParentSectionId();
			if (!parentPaleoIndexes.containsKey(parentID))
				parentPaleoIndexes.put(parentID, new ArrayList<>());
			parentPaleoIndexes.get(parentID).add(s);
			String faultName = faults.getFaultName(parentID);
			if (faultName != null) {
				if (!faultPaleoIndexes.containsKey(faultName))
					faultPaleoIndexes.put(faultName, new ArrayList<>());
				faultPaleoIndexes.get(faultName).add(s);
			}

			double[] logLikelihoods = new double[bVals.size()];
			double[] posteriorWeights = new double[bVals.size()];
			double sumWeights = 0d;
			for (int i=0; i<posteriorWeights.length; i++) {
				double b = bVals.getX(i);
				double priorDensity = priorDist.density(b);
				Preconditions.checkState(Double.isFinite(priorDensity) && priorDensity > 0d);
				
				double misfit = paleoBValMisfits.get(s).getY(i);
				Preconditions.checkState(Double.isFinite(misfit));

				// Unnormalized posterior density:
				// prior density * Gaussian likelihood
				double logLikelihood = -0.5 * misfit * misfit;
				logLikelihoods[i] = logLikelihood;
				posteriorWeights[i] = priorDensity * Math.exp(logLikelihood);
				sumWeights += posteriorWeights[i];
			}

			Preconditions.checkState(Double.isFinite(sumWeights) && sumWeights > 0d);

			EvenlyDiscretizedFunc posteriorPDF = new EvenlyDiscretizedFunc(
					bVals.getMinX(), bVals.getMaxX(), bVals.size());

			// Normalize such that sum(pdf[i] * delta) = 1
			double normalization = sumWeights * posteriorPDF.getDelta();

			for (int i=0; i<posteriorWeights.length; i++)
				posteriorPDF.set(i, posteriorWeights[i] / normalization);
			
			siteLogLikelihoods.add(logLikelihoods);
			sitePosteriorDists.add(new EvenlyDiscrFuncContinuousDistribution(posteriorPDF, DiscretizationType.INTERPOLATE));
			
			siteRups[s] = new BitSet(numRups);
			for (int r : rs.getRupturesForSection(paleoConstr.sectionIndex))
				siteRups[s].set(r);
		}
		
		System.out.println("Building section weighted posteriors");
		int numSects = rs.getNumSections();
		List<CompletableFuture<SectResult>> sectBDistFutures = new ArrayList<>(numSects);
		
		for (int s=0; s<numSects; s++) {
			int sectIndex = s;
			int parentID = rs.getFaultSectionData(sectIndex).getParentSectionId();
			String faultName = faults.getFaultName(parentID);
			
			sectBDistFutures.add(CompletableFuture.supplyAsync(() -> {
				List<Integer> connectedPaleoIndexes;
				if (faultName != null && faultPaleoIndexes.containsKey(faultName)) {
					connectedPaleoIndexes = faultPaleoIndexes.get(faultName);
				} else if (parentPaleoIndexes.containsKey(parentID)) {
					connectedPaleoIndexes = parentPaleoIndexes.get(parentID);
				} else {
					// no connected paleo data, use prior
					return new SectResult(priorDist, null, 1d);
				}
				
				IncrementalMagFreqDist mfd = connectivityWeightTargetMFDs.getOnFaultSupraSeisNucleationMFDs().get(sectIndex);
				
				if (mfd.calcSumOfY_Vals() == 0d)
					return new SectResult(priorDist, null, 1d);
				
				/*
				 * Weighting scheme:
				 * 
				 * figure out estimated fractional rate of each rupture by converting the MFD to estimated participation
				 * rates and dividing the total rate for each magnitude bin evenly among them
				 * 
				 * for each rupture, give that fractional rate as weight to each paleo PDF that it hits (divided evenly
				 * among them, or to the prior if none)
				 */
				
				// bin ruptures by magnitude
				int mMinIndex = mfd.getClosestXIndex(rs.getMinMagForSection(sectIndex));
				int mMaxIndex = mfd.getClosestXIndex(rs.getMaxMagForSection(sectIndex));
				int numMag = 1 + mMaxIndex - mMinIndex;
				List<List<Integer>> magRupIndexes = new ArrayList<>(numMag);
				for (int m=0; m<numMag; m++)
					magRupIndexes.add(new ArrayList<>());
				double[] rupAreaSums = new double[numMag];
				for (int rupIndex : rs.getRupturesForSection(sectIndex)) {
					if (includedRups == null || includedRups.get(rupIndex)) {
						double mag = rs.getMagForRup(rupIndex);
						int magIndex = mfd.getClosestXIndex(mag) - mMinIndex;
						rupAreaSums[magIndex] += rs.getAreaForRup(rupIndex);
						magRupIndexes.get(magIndex).add(rupIndex);
					}
				}
				
				double sectArea = rs.getAreaForSection(sectIndex);
				
				double[] weightPerPaleo = new double[paleoConstraints.size()];
				double weightNoPaleo = 0d;
				BitSet rupPaleoIndexes = new BitSet(paleoConstraints.size());
				double sectParticRate = 0d;
				for (int m=0; m<numMag; m++) {
					List<Integer> rups = magRupIndexes.get(m);
					if (rups.isEmpty())
						continue;
					double nuclRate = mfd.getY(m+mMinIndex);
					if (nuclRate == 0d)
						continue;
					double particScalar = rupAreaSums[m] / (rups.size() * sectArea);
					double particRate = nuclRate * particScalar;
					double particRateEach = particRate / rups.size();
					Preconditions.checkState(Double.isFinite(particRateEach) && particRateEach > 0d);
					for (int rupIndex : rups) {
						rupPaleoIndexes.clear();
						for (int paleoIndex : connectedPaleoIndexes)
							if (siteRups[paleoIndex].get(rupIndex))
								rupPaleoIndexes.set(paleoIndex);
						int numPaleo = rupPaleoIndexes.cardinality();
						if (numPaleo == 0) {
							// rupture hits no paleo sites
							weightNoPaleo += particRateEach;
						} else {
							double rupRatePerPaleo = particRateEach / (double)numPaleo;
							for (int i = rupPaleoIndexes.nextSetBit(0); i >= 0; i = rupPaleoIndexes.nextSetBit(i + 1))
								weightPerPaleo[i] += rupRatePerPaleo;
						}
					}
					sectParticRate += particRate;
				}
				
				Preconditions.checkState(sectParticRate > 0d);
				
				// normalize rates to weights
				weightNoPaleo /= sectParticRate;
				for (int p=0; p<weightPerPaleo.length; p++)
					weightPerPaleo[p] /= sectParticRate;
				
				ContinuousDistribution distribution;
				if (sectUseWeightedDistribution) {
					WeightedList<ContinuousDistribution> distWeights = new WeightedList<>();
					
					if (weightNoPaleo > 0d) {
						// weight for ruptures that hit no paleo sites
						distWeights.add(priorDist, weightNoPaleo);
					}
					
					for (int paleoIndex : connectedPaleoIndexes) {
						double paleoWeight = weightPerPaleo[paleoIndex];
						if (paleoWeight > 0d)
							distWeights.add(sitePosteriorDists.get(paleoIndex), paleoWeight);
					}
					Preconditions.checkState(!distWeights.isEmpty());
					Preconditions.checkState(distWeights.isNormalized());
					distribution = new WeightedContinuousDistribution(distWeights);
				} else {
					EvenlyDiscretizedFunc posterior = bVals.deepClone();

					double maxLogPosterior = Double.NEGATIVE_INFINITY;

					for (int bIndex=0; bIndex<posterior.size(); bIndex++) {
						double b = bVals.getX(bIndex);
						double logPosterior = Math.log(priorDist.density(b));

						for (int paleoIndex : connectedPaleoIndexes) {
							double paleoWeight = weightPerPaleo[paleoIndex];
							if (paleoWeight > 0d) {
								double logLikelihood = siteLogLikelihoods.get(paleoIndex)[bIndex];

								logPosterior += paleoWeight * logLikelihood;
							}
						}

						posterior.set(bIndex, logPosterior);
						maxLogPosterior = Math.max(maxLogPosterior, logPosterior);
					}

					// exponentiate and normalize
					double sum = 0d;
					for (int bIndex=0; bIndex<posterior.size(); bIndex++) {
						double density = Math.exp(posterior.getY(bIndex) - maxLogPosterior);
						posterior.set(bIndex, density);
						sum += density;
					}

					double scalar = 1d / (sum * posterior.getDelta());
					posterior.scale(scalar);
					
					distribution = new EvenlyDiscrFuncContinuousDistribution(posterior, DiscretizationType.INTERPOLATE);
				}
				return new SectResult(distribution, weightPerPaleo, weightNoPaleo);
			}));
		}

		List<SectResult> sectBDists = new ArrayList<>(numSects);
		for (CompletableFuture<SectResult> future : sectBDistFutures)
			sectBDists.add(future.join());
		
		System.out.println("DONE, plotting individual dists");
		// now plots
		Map<Integer, List<FaultSection>> parentMappedSects = rs.getFaultSectionDataList().stream().collect(
				Collectors.groupingBy(S->S.getParentSectionId()));
		EvenlyDiscretizedFunc priorFunc = bVals.deepClone();
		for (int i=0; i<priorFunc.size(); i++)
			priorFunc.set(i, priorDist.density(bVals.getX(i)));
		priorFunc.setName("Prior");
		PlotCurveCharacterstics priorDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_green);
		PlotCurveCharacterstics parentDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_blue);
		PlotCurveCharacterstics otherDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY);
		PlotCurveCharacterstics posteriorDistChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Colors.tab_orange);
		PlotCurveCharacterstics posteriorAvgDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK);
		
		DiscretizedFunc thicknessForWeights = new ArbitrarilyDiscretizedFunc();
		thicknessForWeights.set(0d, 1d);
		thicknessForWeights.set(0.5, 4d);
		thicknessForWeights.set(1d, 5d);
		for (int parentID : parentMappedSects.keySet()) {
			String faultName = faults.getFaultName(parentID);
			List<Integer> connectedPaleoIndexes;
			if (faultName != null && faultPaleoIndexes.containsKey(faultName))
				connectedPaleoIndexes = faultPaleoIndexes.get(faultName);
			else if (parentPaleoIndexes.containsKey(parentID))
				connectedPaleoIndexes = parentPaleoIndexes.get(parentID);
			else
				continue;
			List<FaultSection> sects = parentMappedSects.get(parentID);
			
			List<EvenlyDiscretizedFunc> pdfFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> pdfChars = new ArrayList<>();
			List<EvenlyDiscretizedFunc> misfitFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> misfitChars = new ArrayList<>();
			
			WeightedList<ContinuousDistribution> mySectDists = new WeightedList<>();
			double[] paleoSiteWeights = new double[paleoConstraints.size()];
			double noPaleoWeight = 0d;
			for (FaultSection sect : sects) {
				SectResult sectDist = sectBDists.get(sect.getSectionId());
				if (sectDist.paleoSiteWeights != null)
					for (int p=0; p<paleoSiteWeights.length; p++)
						paleoSiteWeights[p] += sectDist.paleoSiteWeights[p];
				mySectDists.add(sectDist.distribution, 1d);
				noPaleoWeight += sectDist.noPaleoWeight;
			}
			for (int p=0; p<paleoSiteWeights.length; p++)
				paleoSiteWeights[p] /= (double)sects.size();
			noPaleoWeight /= (double)sects.size();
			
			pdfFuncs.add(priorFunc);
			pdfChars.add(getForThickness(priorDistChar, (float)thicknessForWeights.getInterpolatedY(noPaleoWeight)));
			
			boolean firstSame = true;
			boolean firstOther = true;
			for (int paleoIndex : connectedPaleoIndexes) {
				int paleoSectIndex = paleoConstraints.get(paleoIndex).sectionIndex;
				if (paleoSiteWeights[paleoIndex] == 0d)
					continue;
				float thickness = (float)thicknessForWeights.getInterpolatedY(paleoSiteWeights[paleoIndex]);
				boolean sameParent = parentID == rs.getFaultSectionData(paleoSectIndex).getParentSectionId();
				EvenlyDiscretizedFunc misfits = paleoBValMisfits.get(paleoIndex);
				ContinuousDistribution dist = sitePosteriorDists.get(paleoIndex);
				EvenlyDiscretizedFunc pdf = bVals.deepClone();
				for (int i=0; i<pdf.size(); i++)
					pdf.set(i, dist.density(pdf.getX(i)));
				if (sameParent) {
					if (firstSame) {
						misfits = misfits.deepClone();
						pdf = pdf.deepClone();
						misfits.setName("Paleo site (this section)");
//						pdf.setName(misfits.getName());
						firstSame = false;
					}
					misfitFuncs.add(misfits);
					misfitChars.add(getForThickness(parentDistChar, thickness));
					pdfFuncs.add(pdf);
					pdfChars.add(getForThickness(parentDistChar, thickness));
				} else {
					if (firstOther) {
						misfits = misfits.deepClone();
						pdf = pdf.deepClone();
						misfits.setName("Paleo site (other connected sections)");
//						pdf.setName(misfits.getName());
						firstOther = false;
					}
					misfitFuncs.add(0, misfits);
					misfitChars.add(0, getForThickness(otherDistChar, thickness));
					pdfFuncs.add(0, pdf);
					pdfChars.add(0, getForThickness(otherDistChar, thickness));
				}
			}

			boolean firstSectDist = true;
			for (int s=0; s<mySectDists.size(); s++) {
				ContinuousDistribution sectDist = mySectDists.getValue(s);
				EvenlyDiscretizedFunc sectPDF = bVals.deepClone();
				for (int i=0; i<sectPDF.size(); i++)
					sectPDF.set(i, sectDist.density(bVals.getX(i)));
				if (firstSectDist) {
					sectPDF.setName("Subsection posteriors");
					firstSectDist = false;
				}
				pdfFuncs.add(sectPDF);
				pdfChars.add(posteriorDistChar);
			}
			ContinuousDistribution poisteriorDist = new WeightedContinuousDistribution(mySectDists);
			EvenlyDiscretizedFunc posteriorFunc = bVals.deepClone();
			for (int i=0; i<posteriorFunc.size(); i++)
				posteriorFunc.set(i, poisteriorDist.density(bVals.getX(i)));
			posteriorFunc.setName("Section average posterior");
			pdfFuncs.add(posteriorFunc);
			pdfChars.add(posteriorAvgDistChar);
			
			String parentName = sects.get(0).getParentSectionName();
			PlotSpec pdfPlot = new PlotSpec(pdfFuncs, pdfChars, parentName, "b-value", "Density");
			pdfPlot.setLegendInset(true);
			
//			PlotSpec misfitsPlot = new PlotSpec(misfitFuncs, misfitChars, parentName, "b-value", "Misfit (z-score)");
//			misfitsPlot.setLegendInset(true);
			
			List<EvenlyDiscretizedFunc> absMisfitFuncs = new ArrayList<>();
			for (EvenlyDiscretizedFunc func : misfitFuncs) {
				EvenlyDiscretizedFunc absFunc = func.deepClone();
				for (int i=0; i<func.size(); i++)
					absFunc.set(i, Math.abs(func.getY(i)));
				absMisfitFuncs.add(absFunc);
			}
			PlotSpec misfitsPlot = new PlotSpec(absMisfitFuncs, misfitChars, parentName, "b-value", "|Misfit (z-score)|");
			misfitsPlot.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			gp.drawGraphPanel(List.of(pdfPlot, misfitsPlot), false, false, List.of(new Range(bVals.getMinX(), bVals.getMaxX())), null);
			
			String prefix = FileNameUtils.simplify(parentName);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 1200, true, true, false);
		}
		
		System.out.println("Plotting along-strike");
		
		for (String faultName : faults.getFaultNames()) {
			List<AlongStrikePlot> plots = new ArrayList<>();
			
			List<FaultSection> sects = new ArrayList<>();
			MinMaxAveTracker latTrack = new MinMaxAveTracker();
			MinMaxAveTracker lonTrack = new MinMaxAveTracker();
			Map<Integer, List<FaultSection>> parentsMap = new HashMap<>();
			for (FaultSection sect : rs.getFaultSectionDataList()) {
				int parentID = sect.getParentSectionId();
				String sectFaultName = faults.getFaultName(parentID);
				if (sectFaultName != null && sectFaultName.equals(faultName)) {
					sects.add(sect);
					if (!parentsMap.containsKey(parentID))
						parentsMap.put(parentID, new ArrayList<>());
					parentsMap.get(parentID).add(sect);
					for (Location loc : sect.getFaultTrace()) {
						latTrack.addValue(loc.lat);
						lonTrack.addValue(loc.lon);
					}
				}
			}
			Range latRange = new Range(latTrack.getMin(), latTrack.getMax());
			Range lonRange = new Range(lonTrack.getMin(), lonTrack.getMax());
			
			boolean latX = latRange.getLength() > 0.6*lonRange.getLength() || faultName.contains("San Andreas");
			String xLabel;
			Range xRange;
			if (latX) {
				xLabel = "Latitude";
				xRange = latRange;
			} else {
				xLabel = "Longitude";
				xRange = lonRange;
			}
			
			List<XY_DataSet> emptySectFuncs = new ArrayList<>();
			for (FaultSection sect : sects) {
				XY_DataSet func = new DefaultXY_DataSet();
				for (Location loc : sect.getFaultTrace()) {
					if (latX)
						func.set(loc.getLatitude(), 0d);
					else
						func.set(loc.getLongitude(), 0d);
				}
				emptySectFuncs.add(func);
			}
			
			// top plot: paleo fits
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DefaultXY_DataSet dataXY = new DefaultXY_DataSet();
			double halfWhisker = 0.005*xRange.getLength();
			
			double minY = Double.POSITIVE_INFINITY;
			double maxY = 0d;
			for (int c=0; c<paleoConstraints.size(); c++) {
				SectMappedUncertainDataConstraint constraint = paleoConstraints.get(c);
				int parentID = rs.getFaultSectionData(constraint.sectionIndex).getParentSectionId();
				if (parentsMap.containsKey(parentID)) {
					double x = latX ? constraint.dataLocation.getLatitude() : constraint.dataLocation.getLongitude();
					dataXY.set(x, constraint.bestEstimate);
					
					BoundedUncertainty range95 = constraint.estimateUncertaintyBounds(UncertaintyBoundType.CONF_95);
					
					funcs.add(line(x-halfWhisker, range95.upperBound, x+halfWhisker, range95.upperBound));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
					
					funcs.add(line(x, range95.lowerBound, x, range95.upperBound));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
					
					funcs.add(line(x-halfWhisker, range95.lowerBound, x+halfWhisker, range95.lowerBound));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
					
					BoundedUncertainty range68 = constraint.estimateUncertaintyBounds(UncertaintyBoundType.CONF_68);
					// check against a smaller range
					// but 68 could even be negative (determined from std dev, whereas 95 might be direct), so guard against that
					double lowerCheck = Math.max(range68.lowerBound, 0.5*(constraint.bestEstimate + range95.lowerBound));
					double upperCheck = Math.min(range68.upperBound, 0.5*(constraint.bestEstimate + range95.upperBound));
					minY = Math.min(minY, lowerCheck);
					maxY = Math.max(maxY, upperCheck);
				}
			}
			if (dataXY.size() > 0) {
				dataXY.setName("Paleo Constraints");
				funcs.add(dataXY);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
			} else {
				continue;
			}
			
			List<CompletableFuture<double[]>> sectPaleoFutures = new ArrayList<>();
			for (int i=0; i<sects.size(); i++) {
				int sectIndex = sects.get(i).getSectionId();
				sectPaleoFutures.add(CompletableFuture.supplyAsync(()->{
					double[] ret = new double[bVals.size()];
					for (int b=0; b<bVals.size(); b++)
						ret[b] = estimatePaleoRate(rs, sectIndex,
								bValTargetMFDs[b].getOnFaultSupraSeisNucleationMFDs().get(sectIndex), paleoProb, includedRups);
					return ret;
				}));
			}
			List<double[]> sectPaleoRates = new ArrayList<>();
			for (int i=0; i<sects.size(); i++)
				sectPaleoRates.add(sectPaleoFutures.get(i).join());
			
			PlotCurveCharacterstics charBounds = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY);
			DecimalFormat oDF = new DecimalFormat("0.#");
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				
				double rate0 = sectPaleoRates.get(i)[0];
				double rate1 = sectPaleoRates.get(i)[bVals.size()-1];
				
				funcs.add(copyAtY(emptyFunc, rate0));
				chars.add(charBounds);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("b={"+oDF.format(bVals.getMinX())+", "+oDF.format(bVals.getMaxX())+"}");
				funcs.add(copyAtY(emptyFunc, rate1));
				chars.add(charBounds);
				minY = Math.min(minY, rate0);
				minY = Math.min(minY, rate1);
				maxY = Math.max(maxY, rate0);
				maxY = Math.max(maxY, rate1);
			}
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				
				double sumWeight = 0d;
				double sumRateWeight = 0d;
				for (int b=0; b<bVals.size(); b++) {
					double weight = priorDist.density(bVals.getX(b));
					sumWeight += weight;
					sumRateWeight += sectPaleoRates.get(i)[b]*weight;
				}
				
				funcs.add(copyAtY(emptyFunc, sumRateWeight/sumWeight));
				chars.add(priorDistChar);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("Average prior");
			}
			PlotCurveCharacterstics postModeChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue);
			double[] modes = new double[sects.size()];
			for (int i=0; i<sects.size(); i++) {
				FaultSection sect = sects.get(i);
				int sectIndex = sect.getSectionId();
				
				ContinuousDistribution dist = sectBDists.get(sectIndex).distribution;
				double maxDensity = 0d;
				double mode = Double.NaN;
				for (int b=0; b<bVals.size(); b++) {
					double density = dist.density(bVals.getX(b));
					if (density > maxDensity) {
						maxDensity = density;
						mode = bVals.getX(b);
					}
				}
				Preconditions.checkState(Double.isFinite(mode));
				modes[i] = mode;
			}
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				int bIndex = bVals.getClosestXIndex(modes[i]);
				
				funcs.add(copyAtY(emptyFunc, sectPaleoRates.get(i)[bIndex]));
				chars.add(postModeChar);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("Modal posterior");
			}
			PlotCurveCharacterstics postAvgChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_orange);
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				FaultSection sect = sects.get(i);
				int sectIndex = sect.getSectionId();
				
				ContinuousDistribution dist = sectBDists.get(sectIndex).distribution;
				double sumWeight = 0d;
				double sumRateWeight = 0d;
				for (int b=0; b<bVals.size(); b++) {
					double weight = dist.density(bVals.getX(b));
					sumWeight += weight;
					sumRateWeight += sectPaleoRates.get(i)[b]*weight;
				}
				
				funcs.add(copyAtY(emptyFunc, sumRateWeight/sumWeight));
				chars.add(postAvgChar);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("Average posterior");
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, " ", xLabel, "Paleo-visible rate estimate");
			plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT, 0.025, 0.025, 0.95, false);
			
			Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minY))), Math.pow(10, Math.ceil(Math.log10(maxY))));
			
			plots.add(new AlongStrikePlot(plot, funcs, chars, yRange, true));
			
			// b-value plot
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			double priorAvg = priorDist.getMean();
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				
				funcs.add(copyAtY(emptyFunc, priorAvg));
				chars.add(priorDistChar);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("Average prior");
			}
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				
				funcs.add(copyAtY(emptyFunc, modes[i]));
				chars.add(postModeChar);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("Modal posterior");
			}
			for (int i=0; i<sects.size(); i++) {
				XY_DataSet emptyFunc = emptySectFuncs.get(i);
				FaultSection sect = sects.get(i);
				int sectIndex = sect.getSectionId();
				
				double avgB = sectBDists.get(sectIndex).distribution.getMean();
				
				funcs.add(copyAtY(emptyFunc, avgB));
				chars.add(postAvgChar);
				if (i == 0)
					funcs.get(funcs.size()-1).setName("Average posterior");
			}
			
			plot = new PlotSpec(funcs, chars, " ", xLabel, "b-value");
			plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT, 0.025, 0.025, 0.95, false);
			
			yRange = new Range(bVals.getMinX()-0.02, bVals.getMaxX()+0.02);
			
			plots.add(new AlongStrikePlot(plot, funcs, chars, yRange, false));
			
			String prefix = "along_strike_"+FileNameUtils.simplify(faultName);
			
			SectBySectDetailPlots.writeAlongStrikePlots(outputDir, prefix, plots, parentsMap, latX, xLabel, xRange, faultName);
		}
		
		System.out.println("DONE");
	}
	
	private static PlotCurveCharacterstics getForThickness(PlotCurveCharacterstics pChar, float thickness) {
		return new PlotCurveCharacterstics(pChar.getLineType(), thickness, pChar.getSymbol(), pChar.getSymbolWidth(), pChar.getColor());
	}
	
	private static record SectResult(ContinuousDistribution distribution, double[] paleoSiteWeights, double noPaleoWeight) {}; 
	
	private static InversionTargetMFDs buildTargetsForB(NSHM23_InvConfigFactory factory, FaultSystemRupSet rupSet,
			LogicTreeBranch<LogicTreeNode> branch, double b) {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		List<LogicTreeNode> nodes = new ArrayList<>();
		SectionSupraSeisBValues.FixedValueLevel fixedB = null;
		boolean replaced = false;
		for (int l=0; l<branch.size(); l++) {
			LogicTreeLevel<? extends LogicTreeNode> level = branch.getLevel(l);
			if (SupraSeisBValues.class.isAssignableFrom(level.getType())) {
				Preconditions.checkState(!replaced);
				System.out.println("Replacing level "+l+": level");
				fixedB = new SectionSupraSeisBValues.FixedValueLevel("Fixed b", "FixedB", b);
				levels.add(fixedB);
				nodes.add(fixedB.getNodes().get(0));
				replaced = true;
			} else {
				levels.add(level);
				nodes.add(branch.getValue(l));
			}
		}
		Preconditions.checkState(replaced);
		branch = new LogicTreeBranch<>(levels, nodes);
		try {
			return factory.updateRuptureSetForBranch(rupSet, branch).requireModule(InversionTargetMFDs.class);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	private static double estimatePaleoRate(FaultSystemRupSet rs, int sectIndex,
			IncrementalMagFreqDist nuclMFD, NSHM23_PaleoProbabilityModel paleoProb, BitSet includedRups) {
		int mMinIndex = nuclMFD.getClosestXIndex(rs.getMinMagForSection(sectIndex));
		int mMaxIndex = nuclMFD.getClosestXIndex(rs.getMaxMagForSection(sectIndex));
		int numMag = 1 + mMaxIndex - mMinIndex;
		List<List<Integer>> magRupIndexes = new ArrayList<>(numMag);
		for (int m=0; m<numMag; m++)
			magRupIndexes.add(new ArrayList<>());
		double[] rupAreaSums = new double[numMag];
		for (int rupIndex : rs.getRupturesForSection(sectIndex)) {
			if (includedRups == null || includedRups.get(rupIndex)) {
				double mag = rs.getMagForRup(rupIndex);
				int magIndex = nuclMFD.getClosestXIndex(mag) - mMinIndex;
				rupAreaSums[magIndex] += rs.getAreaForRup(rupIndex);
				magRupIndexes.get(magIndex).add(rupIndex);
			}
		}
		
		double sectArea = rs.getAreaForSection(sectIndex);
		
		double ret = 0d;
		
		for (int m=0; m<numMag; m++) {
			List<Integer> rups = magRupIndexes.get(m);
			if (rups.isEmpty())
				continue;
			double nuclRate = nuclMFD.getY(m+mMinIndex);
			if (nuclRate == 0d)
				continue;
			double particScalar = rupAreaSums[m] / (rups.size() * sectArea);
			double particRate = nuclRate * particScalar;
			double particRateEach = particRate / rups.size();
			Preconditions.checkState(Double.isFinite(particRateEach) && particRateEach > 0d);
			for (int rupIndex : rups)
				ret += particRateEach * paleoProb.getProbPaleoVisible(rs, rupIndex, sectIndex);
		}
		return ret;
	}
	
	static XY_DataSet copyAtY(XY_DataSet func, double y) {
		double[] xVals = new double[func.size()];
		double[] yVals = new double[xVals.length];
		for (int i=0; i<xVals.length; i++) {
			xVals[i] = func.getX(i);
			yVals[i] = y;
		}
		return new DefaultXY_DataSet(xVals, yVals);
	}

	private static XY_DataSet line(double x1, double y1, double x2, double y2) {
		return new DefaultXY_DataSet(new double[] { x1, x2 }, new double[] { y1, y2 });
	}
	
}