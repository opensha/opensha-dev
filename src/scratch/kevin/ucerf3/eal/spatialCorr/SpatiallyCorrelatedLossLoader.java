package scratch.kevin.ucerf3.eal.spatialCorr;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

public class SpatiallyCorrelatedLossLoader {

	public static void main(String[] args) throws IOException {
		File resultsDir = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "2020_08_21-ucerf3-ngaw2-cea-100pct-spatial-randTau100/results");
		
		// do all of them
//		File[] csvFiles = resultsDir.listFiles();
		// just do one
		File[] csvFiles = { new File(resultsDir, "ASK_2014_EpiNONE_Wills2015.csv") };
		
		for (File file : csvFiles) {
			if (!file.getName().endsWith(".csv"))
				continue;
			String name = file.getName();
			name = name.substring(0, name.indexOf(".csv"));
			name = name.replaceAll("_", " ");
			
			System.out.println("Loading "+name);
			
			CSVFile<String> csv = CSVFile.readFile(file, true);
			// original mean, between-event samples, list of losses for within-event samples
			Table<Double, Double, List<Double>> resultsTable = HashBasedTable.create();
			
			int numZeros = 0;
			int numVals = 0;
			
			for (int row=1; row<csv.getNumRows(); row++) {
				double meanLoss = csv.getDouble(row, 6);
				double between = csv.getDouble(row, 8);
				Preconditions.checkState(!resultsTable.contains(meanLoss, between));
				List<Double> vals = new ArrayList<>();
				for (int col=9; col<csv.getNumCols(); col++) {
					double val = csv.getDouble(row, col);
					vals.add(val);
					if (val == 0d)
						numZeros++;
					numVals++;
				}
				resultsTable.put(meanLoss, between, vals);
			}
			double zeroPercent = 100d*(double)numZeros/(double)numVals;
			System.out.println("\t"+numZeros+" zero values in total ("+(float)zeroPercent+" %)");
			
			for (Double origMean : resultsTable.rowKeySet()) {
				// for this mean value, this is a map of:
				// <between event sample, list of loss values for within-event samples>
				Map<Double, List<Double>> betweenToLosses = resultsTable.row(origMean);
				double linearMean = 0d;
				double logMeanNonZero = 0d;
				
				List<Double> nonZeroLogVals = new ArrayList<>();
				
				int totalNum = 0;
				for (Double betweenSample : betweenToLosses.keySet()) {
					List<Double> lossVals = betweenToLosses.get(betweenSample);
					for (Double val : lossVals) {
						linearMean += val;
						double logVal = Math.log(val);
						if (val > 0) {
							logMeanNonZero += logVal;
							nonZeroLogVals.add(logVal);
						}
						totalNum++;
					}
				}
				// normalize
				linearMean /= (double)totalNum;
				logMeanNonZero /= (double)nonZeroLogVals.size();
				
				System.out.println("\tRupture with original meanLoss="+origMean+", "
						+nonZeroLogVals.size()+"/"+totalNum+" nonzero");
				System.out.println("\t\tLinear mean: "+ linearMean);
				System.out.println("\t\tNonzero log mean: "+logMeanNonZero);
				double variance = StatUtils.variance(Doubles.toArray(nonZeroLogVals));
				double stdDev = Math.sqrt(variance);
				System.out.println("\t\tNonzero log std dev: "+stdDev);
			}
		}
	}

}
