package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ReweightEvenFitSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;

import com.google.common.base.Preconditions;

public class ConstraintWeightPhaseOutStats {

	public static void main(String[] args) throws IOException {
		File resultsDir = new File(args[0]);
		
		int countThreshold = -1;
		boolean largestOnly = false;
		if (args.length > 1) {
			if (args[1].toLowerCase().equals("true")) {
				largestOnly = true;
				System.out.println("Will only apply to the largest cluster");
			} else if (args[1].toLowerCase().equals("false")) {
				largestOnly = false;
			} else {
				countThreshold = Integer.parseInt(args[1]);
				System.out.println("Will only apply to clusters with at least "+countThreshold+" subsections");
			}
		}

		Map<String, Double> origWeights = new LinkedHashMap<>();
		origWeights.put("Uncertain Slip Rate", 1d);
		origWeights.put("Uncertain MFD Equality", 10d);
		origWeights.put("Paleoseismic Event Rate", 5d);
		origWeights.put("Paleoseismic Average Slip", 5d);
		origWeights.put("Uncertain Parkfield", 10d);
		origWeights.put("Uncertain Subsection Total Nucleation Rates", 0.5d);
		
		Map<String, Integer> counts = new HashMap<>();
		Map<String, Double> sumWeights = new HashMap<>();
		Map<String, Double> maxWeights = new HashMap<>();
		Map<String, Double> minWeights = new HashMap<>();
		Map<String, Integer> partialPhaseCounts = new HashMap<>();
		Map<String, Integer> fullPhaseCounts = new HashMap<>();
		Map<String, Integer> maxCounts = new HashMap<>();
		
		int numWithFullPhase = 0;
		int numWithPartialPhase = 0;
		int numWithMaxAdjustment = 0;
		
		int branchCount = 0;
		for (File runDir : resultsDir.listFiles()) {
			if (!runDir.isDirectory())
				continue;
			File solFile = new File(runDir, "solution.zip");
			if (!solFile.exists())
				continue;
			
			System.out.println("Loading branch "+branchCount+": "+solFile.getAbsolutePath());
			
			ZipFile zip = new ZipFile(solFile);
			
			ZipEntry entry = zip.getEntry("solution/"+ConnectivityClusters.CLUSTER_MISFITS_FILE_NAME);
			Preconditions.checkNotNull(entry);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			final int nameCol = 3;
			final int weightCol = 5;
			final int countCol = 1;
			zip.close();
			
			int largestCount = -1;
			if (largestOnly)
				for (int row=1; row<csv.getNumRows(); row++)
					largestCount = Integer.max(largestCount, csv.getInt(row, countCol));
			
			boolean anyPhased = false;
			boolean anyFullyPhased = false;
			boolean anyMax = false;
			
			for (String name : origWeights.keySet()) {
				double origWeight = origWeights.get(name);
				
				boolean myPhased = false;
				boolean myFullyPhased = false;
				boolean myMax = false;
				for (int row=1; row<csv.getNumRows(); row++) {
					if (largestOnly && csv.getInt(row, countCol) < largestCount)
						continue;
					if (countThreshold > 0 && csv.getInt(row, countCol) < countThreshold)
						continue;
					if (name.equals(csv.get(row, nameCol))) {
						// match
						double weight = csv.getDouble(row, weightCol);
						double factor = origWeight/weight;
						Preconditions.checkState((float)factor >= (float)(1d/ReweightEvenFitSimulatedAnnealing.MAX_ADJUSTMENT_FACTOR)
								&& (float)factor <= (float)ReweightEvenFitSimulatedAnnealing.MAX_ADJUSTMENT_FACTOR,
								"Bad orig weight? factor=%s with orig=%s, final=%s, for %s", factor, origWeight, weight, name);
						
						increment(counts, name);
						
						if (maxWeights.containsKey(name)) {
							sumWeights.put(name, sumWeights.get(name)+weight);
							maxWeights.put(name, Math.max(weight, maxWeights.get(name)));
							minWeights.put(name, Math.min(weight, minWeights.get(name)));
						} else {
							sumWeights.put(name, weight);
							maxWeights.put(name, weight);
							minWeights.put(name, weight);
						}
						
						if ((float)factor > (float)ReweightEvenFitSimulatedAnnealing.PHASE_OUT_START_FACTOR) {
							myPhased = true;
							
							if ((float)factor == (float)ReweightEvenFitSimulatedAnnealing.PHASE_OUT_END_FACTOR) {
								myFullyPhased = true;
								System.out.println(name+" was fully phased out ("+(float)factor
										+"x, orig="+(float)origWeight+", final="+(float)weight+")");
							} else {
								System.out.println(name+" was partially phased out ("+(float)factor
										+"x, orig="+(float)origWeight+", final="+(float)weight+")");
							}
						} else if ((float)factor == (float)(1d/ReweightEvenFitSimulatedAnnealing.MAX_ADJUSTMENT_FACTOR)) {
							myMax = true;
							System.out.println(name+" was set to max factor of "+(float)factor
									+"x, orig="+(float)origWeight+", final="+(float)weight);
						}
					}
				}
				if (myPhased) {
					increment(partialPhaseCounts, name);
					anyPhased = true;
					if (myFullyPhased) {
						anyFullyPhased = true;
						increment(fullPhaseCounts, name);
					}
				}
				if (myMax) {
					increment(maxCounts, name);
					anyMax = true;
				}
			}
			
			if (anyPhased)
				numWithPartialPhase++;
			if (anyFullyPhased)
				numWithFullPhase++;
			if (anyMax)
				numWithMaxAdjustment++;
			
			branchCount++;
		}
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		System.out.println(numWithPartialPhase+"/"+branchCount+" ("
				+pDF.format((double)numWithPartialPhase/(double)branchCount)
				+") had a phased out (fully or partially) constraint");
		System.out.println(numWithFullPhase+"/"+branchCount+" ("
				+pDF.format((double)numWithFullPhase/(double)branchCount)
				+") had a fully phased out constraint");
		System.out.println(numWithMaxAdjustment+"/"+branchCount+" ("
				+pDF.format((double)numWithMaxAdjustment/(double)branchCount)
				+") had a max weight adjustment");
		
		System.out.println("Constraint phase out stats:");
		
		for (String name : origWeights.keySet()) {
			Integer count = counts.get(name);
			if (count == null)
				continue;
			double origWeight = origWeights.get(name);
			double minWeight = minWeights.get(name);
			double maxWeight = maxWeights.get(name);
			System.out.println(name+":");
			
			double average = sumWeights.get(name)/count;
			System.out.println("\tOriginal weight: "+(float)origWeight);
			System.out.println("\tAverage weight: "+(float)average+" ("+(float)(average/origWeight)+" x)");
			System.out.println("\tMin weight: "+(float)minWeight+" ("+(float)(minWeight/origWeight)+" x)");
			System.out.println("\tMax weight: "+(float)maxWeight+" ("+(float)(maxWeight/origWeight)+" x)");
			
			Integer partialCount = partialPhaseCounts.get(name);
			Integer fullCount = fullPhaseCounts.get(name);
			Integer maxCount = maxCounts.get(name);
			if (partialCount == null)
				partialCount = 0;
			if (fullCount == null)
				fullCount = 0;
			if (maxCount == null)
				maxCount = 0;
			System.out.println("\t"+partialCount+"/"+branchCount+" at least partially phased out ("
					+pDF.format((double)partialCount/(double)branchCount)+")");
			System.out.println("\t"+fullCount+"/"+branchCount+" at fully phased out ("
					+pDF.format((double)fullCount/(double)branchCount)+")");
			System.out.println("\t"+maxCount+"/"+branchCount+" had a max weight adjustment ("
					+pDF.format((double)maxCount/(double)branchCount)+")");
		}
	}
	
	private static void increment(Map<String, Integer> counts, String name) {
		Integer prev = counts.get(name);
		if (prev == null)
			prev = 0;
		counts.put(name, prev+1);
	}

}
