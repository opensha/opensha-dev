package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.gui.plot.GraphWindow;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.logicTree.TreeTrimmer;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class CarrizoRateHist {
	
	public static void main(String[] args) throws IOException {
		File file = new File("/tmp/carrizo_paleo_obs_rates.csv");
		
		CSVFile<String> csv = CSVFile.readFile(file, true);
		
		HistogramFunction hist = new HistogramFunction(0d, 40, 10d);
		HistogramFunction branchHist = new HistogramFunction(0d, 40, 10d);
		HistogramFunction subHist = new HistogramFunction(0d, 40, 10d);
		
		LogicTreeBranchNode<?> branchNode = TotalMag5Rate.RATE_9p6;
		
		double subWeightTot = 0d;
		double subRate = 0d;

		double fm31WeightTot = 0d;
		double fm31Rate = 0d;
		
		double weightTot = 0d;
		double rateTot = 0;
		
		TreeTrimmer subsetTrim = LogicTreePBSWriter.getCustomTrimmer();
		List<Class<? extends LogicTreeBranchNode<?>>> nodeClasses = LogicTreeBranch.getLogicTreeNodeClasses();
		
		Map<LogicTreeBranch, Double> ratesMap = Maps.newHashMap();
		Map<LogicTreeBranch, Double> weightsMap = Maps.newHashMap();
		
		//9, 10
		for (int i=1; i<csv.getNumRows(); i++) {
			List<String> line = csv.getLine(i);
			double weight = Double.parseDouble(line.get(9));
			double rate = Double.parseDouble(line.get(10));
			double ri = 1d/rate;
			
			weightTot += weight;
			rateTot += rate*weight;
			
			hist.add(ri, weight);
			
			List<String> branchVals = line.subList(0, 9);
			List<LogicTreeBranchNode<?>> vals = Lists.newArrayList();
			for (int j=0; j<nodeClasses.size(); j++) {
				for (LogicTreeBranchNode<?> node : nodeClasses.get(j).getEnumConstants()) {
					if (node.getShortName().equals(branchVals.get(j))) {
						vals.add(node);
						break;
					}
				}
			}
			LogicTreeBranch branch = LogicTreeBranch.fromValues(vals);
			Preconditions.checkState(branch.isFullySpecified(), branch);
			
			if (subsetTrim.isTreeValid(branch)) {
				subRate += rate*weight;
				subWeightTot += weight;
				subHist.add(ri, weight);
			}
			
			if (branch.getValueUnchecked((Class<? extends LogicTreeBranchNode<?>>) branchNode.getClass()) == branchNode) {
				fm31Rate += rate*weight;
				fm31WeightTot += weight;
				branchHist.add(ri, weight);
			}
			
			ratesMap.put(branch, rate);
			weightsMap.put(branch, weight);
		}
		
		rateTot /= weightTot;
		double meanRI = 1d/rateTot;
		System.out.println("Mean RI: "+meanRI);
		
		subRate /= subWeightTot;
		double subRI = 1d/subRate;
		System.out.println("Subset Mean RI: "+subRI);
		
		fm31Rate /= fm31WeightTot;
		double fm31RI = 1d/fm31Rate;
		System.out.println(branchNode+" Mean RI: "+fm31RI);
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(hist);
		funcs.add(branchHist);
		funcs.add(subHist);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		new GraphWindow(funcs, "Carrizo RI Dist", chars);
		
		for (Class<? extends LogicTreeBranchNode<?>> clazz : nodeClasses) {
			System.out.println(ClassUtils.getClassNameWithoutPackage(clazz));
			for (LogicTreeBranchNode<?> val : clazz.getEnumConstants()) {
				if (val.getRelativeWeight(InversionModels.CHAR_CONSTRAINED) > 0) {
					double totWt = 0;
					double totRate = 0;
					for (LogicTreeBranch branch : ratesMap.keySet()) {
						if (branch.getValueUnchecked((Class<? extends LogicTreeBranchNode<?>>)clazz) == val) {
							double weight = weightsMap.get(branch);
							totRate += ratesMap.get(branch)*weight;
							totWt += weight;
						}
					}
					if (totRate == 0)
						continue;
					totRate /= totWt;
					double ri = 1d/totRate;
					System.out.println("\t"+val.getShortName()+":\t"+(float)ri);
				}
			}
		}
	}

}
