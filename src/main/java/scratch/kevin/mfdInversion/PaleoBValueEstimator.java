package scratch.kevin.mfdInversion;

import java.awt.Color;
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
import org.jfree.data.Range;
import org.opensha.commons.data.WeightedContinuousDistribution;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscrFuncContinuousDistribution;
import org.opensha.commons.data.function.EvenlyDiscrFuncContinuousDistribution.DiscretizationType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoDataLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoProbabilityModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SectionSupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class PaleoBValueEstimator {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/sect_b_val_est");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
		branch.setValue(NSHM23_DeformationModels.AVERAGE);
		System.out.println(branch);
		
		NSHM23_InvConfigFactory.APPLY_DEF_MODEL_UNCERTAINTIES_DEFAULT = false;
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		FaultSystemRupSet rs = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
		ContinuousDistribution priorDist = UniformContinuousDistribution.of(0d, 1d);
		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 11);
		double priorMeanB = priorDist.getMean();
		
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
		InversionTargetMFDs priorMeanTargetMFDs = buildTargetsForB(factory, rs, branch, priorMeanB);
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
				double misfit = (paleoVisibleRate - paleoConstr.bestEstimate)/paleoConstr.getPreferredStdDev();
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
				double likelihood = Math.exp(-0.5 * misfit * misfit);
				posteriorWeights[i] = priorDensity * likelihood;
				sumWeights += posteriorWeights[i];
			}

			Preconditions.checkState(Double.isFinite(sumWeights) && sumWeights > 0d);

			EvenlyDiscretizedFunc posteriorPDF = new EvenlyDiscretizedFunc(
					bVals.getMinX(), bVals.getMaxX(), bVals.size());

			// Normalize such that sum(pdf[i] * delta) = 1
			double normalization = sumWeights * posteriorPDF.getDelta();

			for (int i=0; i<posteriorWeights.length; i++)
				posteriorPDF.set(i, posteriorWeights[i] / normalization);
			
			sitePosteriorDists.add(new EvenlyDiscrFuncContinuousDistribution(posteriorPDF, DiscretizationType.INTERPOLATE));
			
			siteRups[s] = new BitSet(numRups);
			for (int r : rs.getRupturesForSection(paleoConstr.sectionIndex))
				siteRups[s].set(r);
		}
		
		System.out.println("Building section weighted posteriors");
		int numSects = rs.getNumSections();
		List<CompletableFuture<ContinuousDistribution>> sectBDistFutures = new ArrayList<>(numSects);
		
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
					// no paleo data, use prior
					return priorDist;
				}
				
				Preconditions.checkState(!connectedPaleoIndexes.isEmpty());
				// for each paleo site, estimate fraction of paleo-visible co-ruptures with that site
				IncrementalMagFreqDist mfd = priorMeanTargetMFDs.getOnFaultSupraSeisNucleationMFDs().get(sectIndex);
				WeightedList<ContinuousDistribution> distWeights = new WeightedList<>();
				for (int paleoIndex : connectedPaleoIndexes) {
					int paleoSectIndex = paleoConstraints.get(paleoIndex).sectionIndex;
					
					double[] sumParticScalars = new double[mfd.size()];
					double[] sumPaleoScalars = new double[mfd.size()];
					int[] binCounts = new int[mfd.size()];
					double sectArea = rs.getAreaForSection(sectIndex);
					for (int rupIndex : rs.getRupturesForSection(sectIndex)) {
						int magIndex = mfd.getClosestXIndex(rs.getMagForRup(rupIndex));
						binCounts[magIndex]++;
						double rupArea = rs.getAreaForRup(rupIndex);
						sumParticScalars[magIndex] += rupArea/sectArea;
						if (siteRups[paleoIndex].get(rupIndex)) {
							// this rupture includes this paleo site
							double paleoVisibleProb =  paleoProb.getProbPaleoVisible(rs, rupIndex, paleoSectIndex);
							sumPaleoScalars[magIndex] += paleoVisibleProb;
						}
					}
					double corupPaleoVisibleRate = 0d;
					for (int m=0; m<mfd.size(); m++) {
						double totNuclRate = mfd.getY(m);
						if (totNuclRate > 0d) {
							double particRate = totNuclRate * sumParticScalars[m]/binCounts[m];
							corupPaleoVisibleRate += particRate * sumPaleoScalars[m]/binCounts[m];
						}
					}
					if (corupPaleoVisibleRate > 0d)
						distWeights.add(sitePosteriorDists.get(paleoIndex), corupPaleoVisibleRate);
				}
				Preconditions.checkState(!distWeights.isEmpty());
				distWeights.normalize();
				return new WeightedContinuousDistribution(distWeights);
			}));
		}

		List<ContinuousDistribution> sectBDists = new ArrayList<>(numSects);
		for (CompletableFuture<ContinuousDistribution> future : sectBDistFutures)
			sectBDists.add(future.join());
		
		System.out.println("DONE, plotting");
		// now plots
		Map<Integer, List<FaultSection>> parentMappedSects = rs.getFaultSectionDataList().stream().collect(
				Collectors.groupingBy(S->S.getParentSectionId()));
		EvenlyDiscretizedFunc priorFunc = bVals.deepClone();
		for (int i=0; i<priorFunc.size(); i++)
			priorFunc.set(i, priorDist.density(bVals.getX(i)));
		priorFunc.setName("Prior");
		PlotCurveCharacterstics priorDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_green);
		PlotCurveCharacterstics parentDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_blue);
		PlotCurveCharacterstics otherDistChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.DARK_GRAY);
		PlotCurveCharacterstics posteriorDistChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Colors.tab_orange);
		PlotCurveCharacterstics posteriorAvgDistChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK);
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
			
			pdfFuncs.add(priorFunc);
			pdfChars.add(priorDistChar);
			
			boolean firstSame = true;
			boolean firstOther = true;
			for (int paleoIndex : connectedPaleoIndexes) {
				int paleoSectIndex = paleoConstraints.get(paleoIndex).sectionIndex;
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
						misfits.setName("Same parent");
						pdf.setName("Same parent");
						firstSame = false;
					}
					misfitFuncs.add(misfits);
					misfitChars.add(parentDistChar);
					pdfFuncs.add(pdf);
					pdfChars.add(parentDistChar);
				} else {
					if (firstOther) {
						misfits = misfits.deepClone();
						pdf = pdf.deepClone();
						misfits.setName("Other connected parents");
						pdf.setName("Other connected parents");
						firstOther = false;
					}
					misfitFuncs.add(0, misfits);
					misfitChars.add(0, otherDistChar);
					pdfFuncs.add(0, pdf);
					pdfChars.add(0, otherDistChar);
				}
			}
			
			WeightedList<ContinuousDistribution> mySectDists = new WeightedList<>();
			boolean firstSectDist = true;
			for (FaultSection sect : sects) {
				ContinuousDistribution sectDist = sectBDists.get(sect.getSectionId());
				mySectDists.add(sectDist, 1d);
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
			posteriorFunc.setName("Subsection average posterior");
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
		System.out.println("DONE");
	}
	
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

}