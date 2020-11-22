package scratch.ned.SpatiallyCorrelatedLossPDF;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

public class SpatiallyCorrelatedLossPDF {
	
	public static void doit() {
		// this data came from: http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2020_08_21-ucerf3-ngaw2-cea-100pct-spatial-randTau100/results.zip
		File resultsDir = new File("/Users/field/workspace/git/opensha-dev/src/scratch/ned/SpatiallyCorrelatedLossPDF/data");

		
		// do all of them
		File[] csvFiles = resultsDir.listFiles();
		// just do one
//		File[] csvFiles = { new File(resultsDir, "ASK_2014_EpiNONE_Wills2015.csv") };
		
		System.out.println("ModalRupInfo"+"\tMag\t"+"OrigMean"+"\t"+"NumNonZero"+"\t"
				+"LinearMean"+"\t"+"LnMeanNonZero"+"\t"+"LnStdDev"+"\t"+"linearMeanCalc"+"\t"+"linearCOVcalc"+"\t"+"linearCOV");

		
		for (File file : csvFiles) {
			if (!file.getName().endsWith(".csv"))
				continue;
			String name = file.getName();
			name = name.substring(0, name.indexOf(".csv"));
//			name = name.replaceAll("_", " ");
			
//			System.out.println("Loading "+name);
			
			CSVFile<String> csv=null;
			try {
				csv = CSVFile.readFile(file, true);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// original mean, between-event samples, list of losses for within-event samples
			Table<Double, Double, List<Double>> resultsTable = HashBasedTable.create();
			
			int numZeros = 0;
			int numVals = 0;
			
//			System.out.println("Working on "+name);
//			ArrayList<Double> magList = new ArrayList<Double>();
			HashMap<Double,Double> magHashMap = new HashMap<Double,Double>();
			
			for (int row=1; row<csv.getNumRows(); row++) {
				double mag = csv.getDouble(row, 5);
				double meanLoss = csv.getDouble(row, 6);
				magHashMap.put(meanLoss, mag);
				double betweenEventIndex = csv.getDouble(row, 7);
//				String epistemicBranch = csv.get(row, 1);
//				if(!epistemicBranch.equals("NONE"))
//					continue;
				Preconditions.checkState(!resultsTable.contains(meanLoss, betweenEventIndex), "Error on "+name+" row="+row+" meanLoss="+meanLoss+" betweenEventIndex="+betweenEventIndex);
				List<Double> vals = new ArrayList<>();
				for (int col=9; col<csv.getNumCols(); col++) {
					double val = csv.getDouble(row, col);
					vals.add(val);
					if (val == 0d)
						numZeros++;
					numVals++;
				}
				resultsTable.put(meanLoss, betweenEventIndex, vals);
			}
			double zeroPercent = 100d*(double)numZeros/(double)numVals;
//			System.out.println("\t"+numZeros+" zero values in total ("+(float)zeroPercent+" %)");
			
			int lossBinIndex = -1;
			for (Double origMean : resultsTable.rowKeySet()) {
				lossBinIndex += 1;
				String fullName = name+"_LossBin"+lossBinIndex;
				// for this mean value, this is a map of:
				// <between event sample, list of loss values for within-event samples>
				Map<Double, List<Double>> betweenToLosses = resultsTable.row(origMean);
				double linearMean = 0d;
				double logMeanNonZero = 0d;
				
				double mag = magHashMap.get(origMean);
				
				List<Double> nonZeroLogVals = new ArrayList<>();
				List<Double> allLinearVals = new ArrayList<>();
				
				int totalNum = 0;
				for (Double betweenSample : betweenToLosses.keySet()) {
					List<Double> lossVals = betweenToLosses.get(betweenSample);
					for (Double val : lossVals) {
						linearMean += val;
						allLinearVals.add(val);
						double logVal = Math.log(val);
						if (val > 0) {
							logMeanNonZero += logVal;
							nonZeroLogVals.add(logVal);
						}
						totalNum++;
					}
				}
				
//				if(nonZeroLogVals.size()==10000) {
//				if(fullName.equals("ASK_2014_EpiLOWER_WaldAllen_LossBin43")) {
//					System.out.println("\nData For ASK_2014_EpiLOWER_WaldAllen_LossBin43");
//					for(double lnVal:nonZeroLogVals)
//						System.out.println(lnVal);
//					System.exit(-1);
//				}
						
						
				// normalize
				linearMean /= (double)totalNum;
				logMeanNonZero /= (double)nonZeroLogVals.size();
				
//				System.out.println("\tRupture with original meanLoss="+origMean+", "
//						+nonZeroLogVals.size()+"/"+totalNum+" nonzero");
//				System.out.println("\t\tLinear mean: "+ linearMean);
//				System.out.println("\t\tNonzero log mean: "+logMeanNonZero);
				double stdDevLogVals = Math.sqrt(StatUtils.variance(Doubles.toArray(nonZeroLogVals)));				
				double linearCOVcalc = Math.sqrt(Math.exp(stdDevLogVals*stdDevLogVals)-1.0);
				
				double linearCOV = Math.sqrt(StatUtils.variance(Doubles.toArray(allLinearVals)))/linearMean;
				
				double linearMeanCalc = Math.exp(logMeanNonZero+stdDevLogVals*stdDevLogVals/2.0);
				
				System.out.println(fullName+"\t"+(float)mag+"\t"+origMean.floatValue()+"\t"+nonZeroLogVals.size()+"\t"
						+(float)linearMean+"\t"+(float)logMeanNonZero+"\t"+(float)stdDevLogVals+"\t"+(float)linearMeanCalc+"\t"+
						(float)linearCOVcalc+"\t"+(float)linearCOV);
			}
		}

	}
	
	
	public static void computeUncertaintiesAndBiases(double mu, double sigma) {

		double targetMean = Math.exp(mu+sigma*sigma/2.0);
		double targetCOV = Math.sqrt(Math.exp(sigma*sigma)-1.0);

		System.out.println("Target mean = "+(float)targetMean);
		System.out.println("Target COV = "+(float)targetCOV);

		int numLossSamples = 10000;
		int numEstimateSamples = 10000;
		
		
		double[] muEstArray = new double[numEstimateSamples];
		double[] sigmaEstArray = new double[numEstimateSamples];
		double[] meanFromLnEstArray = new double[numEstimateSamples];
		double[] covFromLnEstArray = new double[numEstimateSamples];
		
		double[] meanArray = new double[numEstimateSamples];
		double[] covArray = new double[numEstimateSamples];
		
		Random r = new Random();
		for(int i=0;i<numEstimateSamples;i++) {
			double[] lnLossSamples = new double[numLossSamples];
			double[] lossSamples = new double[numLossSamples];
			for(int j=0;j<numLossSamples;j++) {
				double randValue = r.nextGaussian()*sigma+mu;
				lnLossSamples[j] = randValue;
				lossSamples[j] = Math.exp(randValue);
			}
			muEstArray[i] = StatUtils.mean(lnLossSamples);
			sigmaEstArray[i] = Math.sqrt(StatUtils.variance(lnLossSamples));

			meanFromLnEstArray[i] = Math.exp(muEstArray[i]+sigmaEstArray[i]*sigmaEstArray[i]/2.0);
			covFromLnEstArray[i] = Math.sqrt(Math.exp(sigmaEstArray[i]*sigmaEstArray[i])-1.0);
			
			meanArray[i] = StatUtils.mean(lossSamples);
			covArray[i] = Math.sqrt(StatUtils.variance(lossSamples))/meanArray[i];
		}
		
		double meanFromLnEst = StatUtils.mean(meanFromLnEstArray);
		double meanStdDevFromLnEst = Math.sqrt(StatUtils.variance(meanFromLnEstArray));
		double meanStdomFromLnEst = meanStdDevFromLnEst/Math.sqrt(numEstimateSamples);
		double meanSpreadFromLnEst = meanStdDevFromLnEst/meanFromLnEst;
		double meanBiasFromLnEst = (meanFromLnEst-targetMean)/meanStdomFromLnEst;
		
		double covFromLnEst = StatUtils.mean(covFromLnEstArray);
		double covStdDevFromLnEst = Math.sqrt(StatUtils.variance(covFromLnEstArray));
		double covStdomFromLnEst = covStdDevFromLnEst/Math.sqrt(numEstimateSamples);
		
		System.out.println("\nmean from Ln: "+(float)meanFromLnEst+"\t"+(float)meanStdDevFromLnEst+
				"\t"+(float)meanStdomFromLnEst+"\tmeanSpread: "+(float)meanSpreadFromLnEst+"\tmeanBias: "+(float)meanBiasFromLnEst);
		System.out.println("cov from Ln: "+(float)covFromLnEst+"\t"+(float)covStdDevFromLnEst+"\t"+(float)covStdomFromLnEst);
		
		
		
		double meanEst = StatUtils.mean(meanArray);
		double meanStdDevEst = Math.sqrt(StatUtils.variance(meanArray));
		double meanStdomEst = meanStdDevEst/Math.sqrt(numEstimateSamples);
		double meanSpreadEst = meanStdDevEst/meanEst;
		double meanBiasEst = (meanEst-targetMean)/meanStdomEst;
		
		double covEst = StatUtils.mean(covArray);
		double covStdDevEst = Math.sqrt(StatUtils.variance(covArray));
		double covStdomEst = covStdDevEst/Math.sqrt(numEstimateSamples);
		
		System.out.println("\nmean: "+(float)meanEst+"\t"+(float)meanStdDevEst+
				"\t"+(float)meanStdomEst+"\tmeanSpread: "+(float)meanSpreadEst+"\tmeanBias: "+(float)meanBiasEst);
		System.out.println("cov: "+(float)covEst+"\t"+(float)covStdDevEst+"\t"+(float)covStdomEst);

	}

	public static void main(String[] args) throws IOException {
		
//		computeUncertaintiesAndBiases(9.28,1.22);
		
		doit();
	}

}
