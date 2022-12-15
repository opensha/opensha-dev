package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ReweightEvenFitSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;

import com.google.common.base.Preconditions;

public class ReweightMatchToTargetStats {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_13-nshm23_u3_hybrid_branches-full_sys_inv-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		LogicTree<?> tree = LogicTree.read(new File(dir, "logic_tree.json"));
		
		ZipFile zip = new ZipFile(new File(dir, "results.zip"));
		
		Map<String, MinMaxAveTracker> constraintFullTracks = new HashMap<>();
		Map<String, MinMaxAveTracker> constraintPhasedInTracks = new HashMap<>();
		Map<String, Integer> constraintPhasedOutCount = new HashMap<>();
		Map<String, Integer> constraintFullPhasedOutCount = new HashMap<>();
		MinMaxAveTracker fullTracks = new MinMaxAveTracker();
		MinMaxAveTracker phasedInTracks = new MinMaxAveTracker();
		for (LogicTreeBranch<?> branch : tree) {
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(InversionMisfitProgress.MISFIT_PROGRESS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += InversionMisfitProgress.MISFIT_PROGRESS_FILE_NAME;
//			System.out.println("Loading "+entryName);
			ZipEntry entry = zip.getEntry(entryName);
			Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			InversionMisfitProgress progress = new InversionMisfitProgress(csv);
			
			Quantity quantity = progress.getTargetQuantity();
			List<Double> targets = progress.getTargetVals();
			double target = targets.get(targets.size()-1);
			
			InversionMisfitStats initialStats = progress.getStats().get(0);
			Map<String, Double> initialWeights = new HashMap<>();
			for (MisfitStats stats : initialStats.getStats())
				initialWeights.put(stats.range.name, stats.range.weight);
			InversionMisfitStats finalStats = progress.getStats().get(targets.size()-1);
			
			for (MisfitStats stats : finalStats.getStats()) {
				ConstraintRange range = stats.range;
				if (range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
					continue;
				String name = range.name;
				double val = stats.get(quantity);
				double pDiff = 100d*(val-target)/target;
				track(constraintFullTracks, name, pDiff);
				fullTracks.addValue(pDiff);
				
				double weight = range.weight;
				
				if (pDiff > 10) {
					System.out.println("Poor fit encountered in "+entryName);
					System.out.println("\t"+name+" weight="+(float)weight+", val="+(float)val+", target="+(float)target+", pDiff="+(float)pDiff);
				}
				
				double initialWeight = initialWeights.get(range.name);
				boolean phased = weight < initialWeight/ReweightEvenFitSimulatedAnnealing.PHASE_OUT_START_FACTOR;
				if (phased) {
					increment(constraintPhasedOutCount, name);
					if ((float)weight <= (float)(initialWeight/ReweightEvenFitSimulatedAnnealing.PHASE_OUT_END_FACTOR))
						increment(constraintFullPhasedOutCount, name);
				} else {
					track(constraintPhasedInTracks, name, pDiff);
					phasedInTracks.addValue(pDiff);
				}
			}
		}
		
		zip.close();
		
		System.out.println("All constraints target pDiffs:");
		System.out.println("\t"+fullTracks);
		System.out.println("All phased-in constraints target pDiffs:");
		System.out.println("\t"+phasedInTracks);
		System.out.println();
		System.out.println("Constraint-specific stats:");
		for (String name : constraintFullTracks.keySet()) {
			System.out.println(name);
			System.out.println("\ttarget pDiffs:");
			System.out.println("\t"+constraintFullTracks.get(name));
			System.out.println("\tphased-in target pDiffs:");
			System.out.println("\t"+constraintPhasedInTracks.get(name));
			System.out.println("\tphase out stats:");
			System.out.println("\tcount: "+constraintPhasedOutCount.get(name)+", full: "+constraintFullPhasedOutCount.get(name));
		}
	}
	
	private static void increment(Map<String, Integer> map, String name) {
		int prev = map.containsKey(name) ? map.get(name) : 0;
		map.put(name, prev+1);
	}
	
	private static void track(Map<String, MinMaxAveTracker> map, String name, double value) {
		MinMaxAveTracker track = map.get(name);
		if (track == null) {
			track = new MinMaxAveTracker();
			map.put(name, track);
		}
		track.addValue(value);
	}

}
