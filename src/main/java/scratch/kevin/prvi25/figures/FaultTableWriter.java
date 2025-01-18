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

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
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
		fw.write("\\tbl{Participation recurrence intervals and magnitude ranges (using average scaling) for crustal faults. Faults consisting of multiple sections or zones are grouped.\\label{tbl:crustal_fault_ris}}\n");
		fw.write("{\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}lccc@{}}\n");
		fw.write("\\textbf{Name} & \\multicolumn{1}{c}{\\textbf{Magnitude range}} & \\multicolumn{1}{c}{\\textbf{Recurrence interval (years)}} \\\\\n");
		fw.write("\\colrule\n");
		DecimalFormat magDF = new DecimalFormat("0.0");
		for (String name : namesSorted) {
			Set<Integer> ids = namesToSectIDs.get(name);
			double mMin = Double.POSITIVE_INFINITY;
			double mMax = 0d;
			double rate = 0d;
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				boolean match = false;
				for (int sectID : rupSet.getSectionsIndicesForRup(r)) {
					if (ids.contains(sectID)) {
						match = true;
						break;
					}
				}
				if (match) {
					double rupRate = sol.getRateForRup(r);
					if (rupRate == 0d)
						continue;
					rate += rupRate;
					double mag = rupSet.getMagForRup(r);
					mMin = Math.min(mMin, mag);
					mMax = Math.max(mMax, mag);
				}
			}
			fw.write(name+" & $["+magDF.format(mMin)+", "+magDF.format(mMax)+"]$ & "
					+LaTeXUtils.groupedIntNumber(1d/rate)+" \\\\\n");
		}
		fw.write("\\botrule\n");
		fw.write("\\end{tabular*}}\n");
		fw.write("{}\n");
		fw.write("\\end{table}\n");
		
		fw.close();
	}

}
