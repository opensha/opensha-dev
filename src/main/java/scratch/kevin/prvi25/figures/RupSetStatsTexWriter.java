package scratch.kevin.prvi25.figures;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.latex.LaTeXUtils;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class RupSetStatsTexWriter {

	public static void main(String[] args) throws IOException {
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		
		File outputDir = FIGURES_DIR;
		
		FileWriter fw = new FileWriter(new File(outputDir, "rup_set_stats.tex"));
		
		LogicTreeBranch<LogicTreeNode> avgCrustalBranch = PRVI25_LogicTree.DEFAULT_CRUSTAL_ON_FAULT.copy();
		avgCrustalBranch.setValue(NSHM23_ScalingRelationships.AVERAGE);
		List<LogicTreeBranch<LogicTreeNode>> crustalScaleBranches = new ArrayList<>();
		for (NSHM23_ScalingRelationships scale : NSHM23_ScalingRelationships.values())
			if (scale.getNodeWeight(null) > 0d)
				crustalScaleBranches.add(copyAndSet(avgCrustalBranch, scale));
		
		fw.write("% Crustal\n");
		write(fw, "Crustal", factory, List.of(avgCrustalBranch), crustalScaleBranches);
		
		LogicTreeBranch<LogicTreeNode> avgSubLargeBranch = PRVI25_LogicTree.DEFAULT_SUBDUCTION_INTERFACE.copy();
		avgSubLargeBranch.setValue(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE);
		avgSubLargeBranch.setValue(PRVI25_SubductionScalingRelationships.AVERAGE);
		List<LogicTreeBranch<LogicTreeNode>> subLargeScaleBranches = new ArrayList<>();
		for (PRVI25_SubductionScalingRelationships scale : PRVI25_SubductionScalingRelationships.values())
			if (scale.getNodeWeight(null) > 0d)
				subLargeScaleBranches.add(copyAndSet(avgSubLargeBranch, scale));
		
		fw.write("\n% Subduction, Large FM\n");
		write(fw, "SubLarge", factory, List.of(avgSubLargeBranch), subLargeScaleBranches);
		
		LogicTreeBranch<LogicTreeNode> avgSubSmallBranch = PRVI25_LogicTree.DEFAULT_SUBDUCTION_INTERFACE.copy();
		avgSubSmallBranch.setValue(PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL);
		avgSubSmallBranch.setValue(PRVI25_SubductionScalingRelationships.AVERAGE);
		List<LogicTreeBranch<LogicTreeNode>> subSmallScaleBranches = new ArrayList<>();
		for (PRVI25_SubductionScalingRelationships scale : PRVI25_SubductionScalingRelationships.values())
			if (scale.getNodeWeight(null) > 0d)
				subSmallScaleBranches.add(copyAndSet(avgSubSmallBranch, scale));
		
		fw.write("\n% Subduction, Small FM\n");
		write(fw, "SubSmall", factory, List.of(avgSubSmallBranch), subSmallScaleBranches);
		
		fw.write("\n% Subduction, Combined FM\n");
		List<LogicTreeBranch<LogicTreeNode>> subAllBranches = new ArrayList<>();
		subAllBranches.addAll(subSmallScaleBranches);
		subAllBranches.addAll(subLargeScaleBranches);
		write(fw, "SubComb", factory, List.of(avgSubLargeBranch, avgSubSmallBranch), subAllBranches);
		
		fw.close();
	}
	
	private static LogicTreeBranch<LogicTreeNode> copyAndSet(LogicTreeBranch<LogicTreeNode> branch, LogicTreeNode node) {
		branch = branch.copy();
		branch.setValue(node);
		return branch;
	}
	
	static final DecimalFormat magDF = new DecimalFormat("0.0");
	static final DecimalFormat groupedIntDF = new DecimalFormat("0");
	static {
		groupedIntDF.setGroupingUsed(true);
		groupedIntDF.setGroupingSize(3);
	}
	
	private static void write(FileWriter fw, String texPrefix, PRVI25_InvConfigFactory factory,
			List<? extends LogicTreeBranch<?>> avgScaleBranches,
			List<? extends LogicTreeBranch<?>> scaleBranches) throws IOException {
		MinMaxAveTracker sectMminTrack = new MinMaxAveTracker();
		MinMaxAveTracker sectMmaxTrack = new MinMaxAveTracker();
		double momentWeightedSectMmin = 0d;
		double momentWeightedSectMmax = 0d;
		double momentSum = 0d;
		
		double weightSum = 0d;
		double avgMmin = 0d;
		double avgMmax = 0d;
		
		for (LogicTreeBranch<?> avgBranch : avgScaleBranches) {
			double weight = avgBranch.requireValue(RupSetFaultModel.class).getNodeWeight(null);
			Preconditions.checkState(weight > 0d);
			FaultSystemRupSet avgRupSet = factory.buildRuptureSet(avgBranch, FaultSysTools.defaultNumThreads());
			avgMmin += weight*avgRupSet.getMinMag();
			avgMmax += weight*avgRupSet.getMaxMag();
			weightSum += weight;
			ClusterRuptures cRups = avgRupSet.requireModule(ClusterRuptures.class);
			BinaryRuptureProbabilityCalc exclusion = PRVI25_InvConfigFactory.getExclusionModel(
					avgRupSet, avgBranch, cRups);
			double myMaxMin = 0d;
			double myMinMin = Double.POSITIVE_INFINITY;
			double mueMaxMin = 0d;
			double mueMinMin = Double.POSITIVE_INFINITY;
			double mueMax = 0d;
			for (int s=0; s<avgRupSet.getNumSections(); s++) {
				FaultSection sect = avgRupSet.getFaultSectionData(s);
				double mMin, mMax;
				if (exclusion == null) {
					mMin = avgRupSet.getMinMagForSection(s);
					mMax = avgRupSet.getMaxMagForSection(s);
				} else {
					mMin = Double.POSITIVE_INFINITY;
					mMax = 0d;
					for (int rupIndex : avgRupSet.getRupturesForSection(s)) {
						if (exclusion.isRupAllowed(cRups.get(rupIndex), false)) {
							double mag = avgRupSet.getMagForRup(rupIndex);
							mMin = Math.min(mMin, mag);
							mMax = Math.max(mMax, mag);
						}
					}
				}
				sectMminTrack.addValue(mMin);
				sectMmaxTrack.addValue(mMax);
				double moment = FaultMomentCalc.getMoment(avgRupSet.getAreaForSection(s), avgRupSet.getSlipRateForSection(s))*weight;
				momentSum += moment;
				momentWeightedSectMmin += moment*mMin;
				momentWeightedSectMmax += moment*mMax;
				myMaxMin = Math.max(myMaxMin, mMin);
				myMinMin = Math.min(myMinMin, mMin);
				if (sect.getParentSectionName().contains("Muertos")) {
					mueMaxMin = Math.max(myMaxMin, mMin);
					mueMinMin = Math.min(myMinMin, mMin);
					mueMax = Math.max(mueMax, mMax);
				}
			}
			if (avgScaleBranches.size() == 1) {
				fw.write(LaTeXUtils.defineValueCommand(texPrefix+"NumRups", groupedIntDF.format(avgRupSet.getNumRuptures()))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(texPrefix+"NumSubsects", groupedIntDF.format(avgRupSet.getNumSections()))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgSectMinMmin", magDF.format(myMinMin))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgSectMaxMmin", magDF.format(myMaxMin))+"\n");
				if (avgBranch.hasValue(PRVI25_SubductionFaultModels.class)) {
					fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgGridMinMmax", magDF.format(myMinMin-0.1))+"\n");
					fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgGridMaxMmax", magDF.format(myMaxMin-0.1))+"\n");
					fw.write(LaTeXUtils.defineValueCommand(texPrefix+"MueAvgSectMinMmin", magDF.format(mueMinMin))+"\n");
					fw.write(LaTeXUtils.defineValueCommand(texPrefix+"MueAvgSectMaxMmin", magDF.format(mueMaxMin))+"\n");
					fw.write(LaTeXUtils.defineValueCommand(texPrefix+"MueAvgSectMmax", magDF.format(mueMax))+"\n");
				}
			}
		}
		momentWeightedSectMmin /= momentSum;
		momentWeightedSectMmax /= momentSum;
		avgMmin /= weightSum;
		avgMmax /= weightSum;
		
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgScaleMmin", magDF.format(avgMmin))+"\n");
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgScaleMmax", magDF.format(avgMmax))+"\n");
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgSectMmin", magDF.format(sectMminTrack.getAverage()))+"\n");
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgSectMmax", magDF.format(sectMmaxTrack.getAverage()))+"\n");
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgSectWeightedMmin", magDF.format(momentWeightedSectMmin))+"\n");
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"AvgSectWeightedMmax", magDF.format(momentWeightedSectMmax))+"\n");
		
		double overallMmin = Double.POSITIVE_INFINITY;
		double overallMmax = 0d;
		for (LogicTreeBranch<?> branch : scaleBranches) {
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
			ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
			BinaryRuptureProbabilityCalc exclusion = PRVI25_InvConfigFactory.getExclusionModel(
					rupSet, branch, cRups);
			double mMin, mMax;
			if (exclusion == null) {
				mMin = rupSet.getMinMag();
				mMax = rupSet.getMaxMag();
			} else {
				mMin = Double.POSITIVE_INFINITY;
				mMax = 0d;
				for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
					if (exclusion.isRupAllowed(cRups.get(rupIndex), false)) {
						double mag = rupSet.getMagForRup(rupIndex);
						mMin = Math.min(mMin, mag);
						mMax = Math.max(mMax, mag);
					}
				}
			}
			overallMmin = Math.min(overallMmin, mMin);
			overallMmax = Math.max(overallMmax, mMax);
		}
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"OverallMmin", magDF.format(overallMmin))+"\n");
		fw.write(LaTeXUtils.defineValueCommand(texPrefix+"OverallMmax", magDF.format(overallMmax))+"\n");
	}

}
