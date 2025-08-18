package scratch.kevin.prvi25.figures;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.latex.LaTeXUtils;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class FaultTableWriter {
	
	static Map<String, String> getNameRemappings() {
		Map<String, String> prefixRemappings = new HashMap<>();
		prefixRemappings.put("Anegada Passage", "Anegada Passage (2 zones)");
		prefixRemappings.put("Bunce", "Bunce (8 sections)");
		prefixRemappings.put("Lesser Antilles", "Lesser Antilles (2 sections)");
		prefixRemappings.put("Main Ridge", "Main Ridge (2 sections)");
		prefixRemappings.put("Cerro Goden PROXY", "Cerro Goden (zone)");
		prefixRemappings.put("Mayaguez Bay PROXY", "Mayaguez Bay (zone)");
		prefixRemappings.put("Mona Passage", "Mona Passage (4 sections \\& 1 zone)");
		prefixRemappings.put("SW Puerto Rico PROXY", "SW Puerto Rico (zone)");
		prefixRemappings.put("Septentrional", "Septentrional (2 sections)");
		prefixRemappings.put("Yabon", "Yabon (7 sections)");
		return prefixRemappings;
	}

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "crustal_sol");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File outputFile = new File(outputDir, "fault_rate_table.tex");
		
		FaultSystemSolution sol = FaultSystemSolution.load(CRUSTAL_SOL_SUPRA_ONLY);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		SolutionLogicTree slt = SolutionLogicTree.load(CRUSTAL_SLT);
		LogicTree<?> logicTree = slt.getLogicTree();
		List<double[]> branchRates = new ArrayList<>(logicTree.size());
		List<Double> branchWeights = new ArrayList<>(logicTree.size());
		for (int i=0; i<logicTree.size(); i++) {
			LogicTreeBranch<?> branch = logicTree.getBranch(i);
			double[] rates = slt.loadRatesForBranch(branch);
			Preconditions.checkState(rates.length == rupSet.getNumRuptures());
			branchRates.add(rates);
			double weight = logicTree.getBranchWeight(i);
			branchWeights.add(weight);
		}
		
		Map<String, String> prefixRemappings = getNameRemappings();
		
		Map<String, Set<Integer>> namesToSectIDs = new HashMap<>();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			String name = sect.getParentSectionName();
			for (String prefix : prefixRemappings.keySet()) {
				if (name.startsWith(prefix)) {
					name = prefixRemappings.get(prefix);
					break;
				}
			}
			Set<Integer> ids = namesToSectIDs.get(name);
			if (ids == null) {
				ids = new HashSet<>();
				namesToSectIDs.put(name, ids);
			}
			ids.add(sect.getSectionId());
		}
		
		List<String> namesSorted = new ArrayList<>(namesToSectIDs.keySet());
		Collections.sort(namesSorted);
		
		FileWriter fw = new FileWriter(outputFile);
		
		fw.write("\\begin{table}[h!]\n");
		fw.write("\\tbl{Participation recurrence intervals (branch-averaged and 95\\% bounds across all "
				+ "\\dynvalCrustalFaultBranches{} branches) and supra-seismogenic magnitude ranges (using average scaling) for crustal faults. "
				+ "Faults with multiple sections or zones are grouped, and the listed recurrence intervals apply to any "
				+ "rupture involving any part of any section.\\label{tbl:crustal_fault_ris}}\n");
		fw.write("{\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}lccc@{}}\n");
		fw.write("\\textbf{Name} "
				+ "& \\multicolumn{1}{c}{\\textbf{Magnitude range}} "
				+ "& \\multicolumn{1}{c}{\\textbf{Recurrence interval (years)}} "
				+ "& \\multicolumn{1}{c}{\\textbf{Branch recurrence 95\\% bounds (years)}} "
				+ "\\\\\n");
		fw.write("\\colrule\n");
		DecimalFormat magDF = new DecimalFormat("0.0");
		for (String name : namesSorted) {
			Set<Integer> ids = namesToSectIDs.get(name);
			double mMin = Double.POSITIVE_INFINITY;
			double mMax = 0d;
			
			List<Integer> rupIDs = new ArrayList<>();
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				for (int sectID : rupSet.getSectionsIndicesForRup(r)) {
					if (ids.contains(sectID)) {
						rupIDs.add(r);
						break;
					}
				}
			}
			for (int rupIndex : rupIDs) {
				double rupRate = sol.getRateForRup(rupIndex);
				if (rupRate == 0d)
					continue;
				double mag = rupSet.getMagForRup(rupIndex);
				mMin = Math.min(mMin, mag);
				mMax = Math.max(mMax, mag);
			}
			double rate = calcRate(rupIDs, sol.getRateForAllRups());
			double ri = 1d/rate;
			double minRate = Double.POSITIVE_INFINITY;
			double maxRate = 0d;
			List<Double> faultRates = new ArrayList<>(branchRates.size());
			for (double[] rates : branchRates) {
				double branchRate = calcRate(rupIDs, rates);
				faultRates.add(branchRate);
				minRate = Math.min(minRate, branchRate);
				maxRate = Math.max(maxRate, branchRate);
			}
			LightFixedXFunc ncdf = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(faultRates, branchWeights);
			double rate2p5 = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(ncdf, 0.025);
			double rate97p5 = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(ncdf, 0.975);
			
			System.out.println(name);
			System.out.println("\tBA rate="+(float)rate);
			System.out.println("\t95% rates=["+(float)rate2p5+", "+(float)rate97p5+"]");
			System.out.println("\tExtrema rates=["+(float)minRate+", "+(float)maxRate+"]");
			System.out.println("\tBA RI="+LaTeXUtils.groupedIntNumber(ri));
			System.out.println("\t95% RIs=["+LaTeXUtils.groupedIntNumber(1d/rate97p5)+", "+LaTeXUtils.groupedIntNumber(1d/rate2p5)+"]");
			System.out.println("\tExtrema RIs=["+LaTeXUtils.groupedIntNumber(1d/maxRate)+", "+LaTeXUtils.groupedIntNumber(1d/minRate)+"]");
			double minRI = 1d/maxRate;
			double maxRI = 1d/minRate;
			Preconditions.checkState(ri > minRI && ri < maxRI, "RI=%s for %s outside extrema: [%s, %s]", (int)ri, name, (int)minRI, (int)maxRI);
			double ri2p5 = 1d/rate97p5;
			double ri97p5 = 1d/rate2p5;
			Preconditions.checkState(ri > ri2p5 && ri < ri97p5, "RI=%s for %s outside 95% bounds: [%s, %s]", (int)ri, name, (int)ri2p5, (int)ri97p5);
			String line = name+" & ["+magDF.format(mMin)+", "+magDF.format(mMax)+"] "
					+ "& "+LaTeXUtils.groupedIntNumber(ri)+" "
					+ "& ["+(int)(ri2p5+0.5)+", "+(int)(ri97p5+0.5)+"]"+" \\\\";
			System.out.println(line);
			fw.write(line+"\n");
		}
		fw.write("\\botrule\n");
		fw.write("\\end{tabular*}}\n");
		fw.write("{}\n");
		fw.write("\\end{table}\n");
		
		fw.close();
	}
	
	private static double calcRate(List<Integer> rupIDs, double[] rates) {
		double rate = 0d;
		for (int rupIndex : rupIDs)
			rate += rates[rupIndex];
		return rate;
	}

}
