package scratch.ned.U3_TreeValuation;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.math.IEEE754rUtils;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.Precision;
import org.jfree.ui.RectangleEdge;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.calc.recurInterval.LognormalDistCalc;
import org.opensha.sha.gui.infoTools.CalcProgressBar;

import com.google.common.io.Files;

import scratch.UCERF3.inversion.CommandLineInversionRunner;

public class Analysis {
	
	/* Data files came from:
	*  http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2018_11_20-ucerf3-ngaw2-cea-consolidate/
	*  
	*  
	*/
	
	final static String ROOT_DIR = "src/scratch//ned/U3_TreeValuation/";
	
	String inputFileName, outDirName;

	String infoString = "";
	double[] branchWt, branchEAL, branchNormEAL;
	HashMap<String, String[]> allBranchValuesMap;
	HashMap<String, ArrayList<String>> optionsForBranchHashMap;
	
	int totNumBranches;
	double totMeanEAL, totMedianNormEAL, totModalNormEAL;
	double totStdDevEAL, totCOV_EAL, totFactor95perc, totWithin10perc;
	
	int dataCol=2;
	
	double totOrigWeight;
	HashMap<String, String> closestValuesForBranchMap; // this stores the branch value that has a mean closest to the totalMeanEAL
	
	// for branch name only
	HashMap<String, Double> covForBrMap;
	HashMap<String, Double> wtAbsValMeanForBrMap;

	
	// for branch name/value
	HashMap<String, Double> weightForBrValMap;
	HashMap<String, Double> meanForBrValMap;
	HashMap<String, Double> meanDiffForBrValMap;
	HashMap<String, Double> meanDiffWtedForBrValMap;
	HashMap<String, Double> covForBrValMap;
	HashMap<String, Double> covDiffForBrValMap;
	HashMap<String, Double> fracWithIn10percForBrValMap;
	HashMap<String, Double> fracWithIn10percDiffForBrValMap;
	HashMap<String, Double> factor95percForBrValMap;
	HashMap<String, Double> factor95percRatioForBrValMap;
	
	HashMap<String, Double> meanIfBrRemovedMap;
	HashMap<String, Double> meanDiffIfBrRemovedMap;
	HashMap<String, Double> meanDiffWtedIfBrRemovedMap;
	HashMap<String, Double> covIfBrRemovedMap;
	HashMap<String, Double> covDiffIfBrRemovedMap;
	HashMap<String, Double> fracWithIn10percIfBrRemovedMap;
	HashMap<String, Double> fracWithIn10percDiffIfBrRemovedMap;
	HashMap<String, Double> factor95percIfBrRemovedMap;
	HashMap<String, Double> mfactor95percRatioDiffIfBrRemovedMap;

	HashMap<String, String> nameChangeMap=null;

	ArrayList<String> erf_branches;
	ArrayList<String> gmm_branches;
	
	
	private double computeWeightedAverage(double[] wtArray, double[] valArray) {
		double wtAve = 0;
		double wtSum = 0;
		for(int i=0;i<wtArray.length;i++) {
			wtAve += wtArray[i]*valArray[i];
			wtSum += wtArray[i];
		}
		return wtAve/wtSum;
	}
	
	
	private double computeWeightedStdDev(double[] wtArray, double[] valArray, double wtAve) {
		double valSum = 0;
		double wtSum = 0;
		double numNonZero=0;
		for(int i=0;i<wtArray.length;i++) {
			if(wtArray[i] > 0) {
				valSum += wtArray[i]*(valArray[i]-wtAve)*(valArray[i]-wtAve);
				wtSum += wtArray[i];
				numNonZero += 1.0;
			}
		}
		double stdDevSquared = (numNonZero*valSum)/((numNonZero-1)*wtSum);
		return Math.sqrt(stdDevSquared);
	}

	/**
	 * This is for AAL data
	 * @param fileInputName
	 * @param outDirName
	 * @param branchesToRemove
	 * @param popupWindows
	 * @param mkDiffForEachBranchPlots
	 */
	public Analysis(String fileInputName, String outDirName, ArrayList<String> branchesToRemove, boolean popupWindows, 
			boolean mkDiffForEachBranchPlots) {
		this(fileInputName, outDirName, branchesToRemove, popupWindows, mkDiffForEachBranchPlots, 2);
	}
	
	/**
	 * Note that output files and plots will say "AAL" even though the analysis may be for the 
	 * other points on the curve (as set by dataCol).
	 * 
	 * @param fileInputName
	 * @param outDirName
	 * @param branchesToRemove
	 * @param popupWindows
	 * @param mkDiffForEachBranchPlots
	 * @param dataCol - 2 for AAL, 18 for Loss@0.01; 19 for	Loss@0.004; 20 for Loss@0.0025; 21 for Loss@0.0018; and 22 for Loss@4.0E-4
	 */
	public Analysis(String fileInputName, String outDirName, ArrayList<String> branchesToRemove, boolean popupWindows, 
			boolean mkDiffForEachBranchPlots, int dataCol) {
		
		this.dataCol = dataCol;
		
		boolean savePlots = true;
		
		this.inputFileName = fileInputName;
		this.outDirName = outDirName;
		
		if(outDirName != null) {
			// make directory if it doesn't exist
			File dir = new File(ROOT_DIR+outDirName);
			if(!dir.exists())
				dir.mkdirs();			
		}

		infoString += "Input File Name: "+inputFileName+"\n";
		infoString += "Output Dir Name: "+outDirName+"\n";
		
		if(dataCol != 2)
			infoString += "\nNOTE: EAL BELOW IS REALLY THE LOSS AT ANOTHER POINT ON THE CURVE, AS INDICATED IN THE FILENAME\n";

		
//		readBranchLevelSummaryDataFromFile();
		readAllBranchDataFromFile();
		
		if(branchesToRemove != null) {
			for(String toRemove:branchesToRemove) {
				String[] nameValArray = toRemove.split(" = ");
				removeBranchFromData(nameValArray[0],nameValArray[1]);
			}			
		}
		
		makeEAL_Historgram(popupWindows, savePlots);
	
		generateBranchValueResults(popupWindows, outDirName != null);
		
		if(mkDiffForEachBranchPlots) {
			
			makeDiffForEachBranchPlot(popupWindows, outDirName != null, wtAbsValMeanForBrMap, "wtAbsValMeanForBrMap", "Expected Fractional Mean Change", 0d, -0.3, 0.3);
			makeDiffForEachBranchPlot(popupWindows, outDirName != null, covForBrMap, "covForBrMap", "Expected COV Change", totCOV_EAL, 0.1, 0.43);
			
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, meanDiffForBrValMap, "meanDiffForBrValMap", "Mean Fractional Change", 0d, -1d, 1d);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, meanDiffWtedForBrValMap, "meanDiffWtedForBrValMap", "Mean Fractional Change", 0d, -0.2, 0.2);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, meanDiffIfBrRemovedMap, "meanDiffIfBrRemovedMap", "Mean Fractional Change", 0d, -0.2, 0.2);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, meanDiffWtedIfBrRemovedMap, "meanDiffWtedIfBrRemovedMap", "Mean Fractional Change", 0d, -0.2, 0.2);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, covForBrValMap, "covForBrValMap", "COV Change", totCOV_EAL, 0.1, 0.601);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, covIfBrRemovedMap, "covIfBrRemovedMap", "COV Change", totCOV_EAL, 0.1, 0.601);
			
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, fracWithIn10percForBrValMap, "fracWithIn10percForBrValMap", "Fraction Within 10%", totWithin10perc, 0d, 0.4);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, fracWithIn10percIfBrRemovedMap, "fracWithIn10percIfBrRemovedMap", "Fraction Within 10%", totWithin10perc, 0d, 0.4);

			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, factor95percForBrValMap, "factor95percForBrValMap", "95% Conf Factor", totFactor95perc, 1d, 3d);
			makeDiffForEachBranchOptionPlot(popupWindows, outDirName != null, factor95percIfBrRemovedMap, "factor95percIfBrRemovedMap", "95% Conf Factor", totFactor95perc, 1d, 3d);
			
//			//  the following tests factor95 and fract10 are close to that for lognormal distribution with associated COV
//			HashMap<String, Double> test_fracWithIn10percForBrValMap = new HashMap<String, Double>();
//			HashMap<String, Double> test_factor95percForBrValMap = new HashMap<String, Double>();
//			for(String key:fracWithIn10percForBrValMap.keySet()) {
//				double cov = new Double(covForBrValMap.get(key));
//				double[] valArray = computeFact95_fract10_fromLogNorm_COV(cov);
//				test_factor95percForBrValMap.put(key, valArray[0]);
//				test_fracWithIn10percForBrValMap.put(key, valArray[1]);
//			}
//			makeDiffForEachBranchPlot(popupWindows, outDirName != null, test_fracWithIn10percForBrValMap, "test_fracWithIn10percForBrValMap", "Fraction Within 10%", totWithin10perc, 0d, 0.4);
//			makeDiffForEachBranchPlot(popupWindows, outDirName != null, test_factor95percForBrValMap, "test_factor95percForBrValMap", "95% Conf Factor", totFactor95perc, 1d, 3d);
	//
//			HashMap<String, Double> test_fracWithIn10percIfBrRemovedMap = new HashMap<String, Double>();
//			HashMap<String, Double> test_factor95percIfBrRemovedMap = new HashMap<String, Double>();
//			for(String key:fracWithIn10percIfBrRemovedMap.keySet()) {
//				double cov = new Double(covIfBrRemovedMap.get(key));
//				double[] valArray = computeFact95_fract10_fromLogNorm_COV(cov);
//				test_factor95percIfBrRemovedMap.put(key, valArray[0]);
//				test_fracWithIn10percIfBrRemovedMap.put(key, valArray[1]);
//			}
//			makeDiffForEachBranchPlot(popupWindows, outDirName != null, test_fracWithIn10percIfBrRemovedMap, "test_fracWithIn10percIfBrRemovedMap", "Fraction Within 10%", totWithin10perc, 0d, 0.4);
//			makeDiffForEachBranchPlot(popupWindows, outDirName != null, test_factor95percIfBrRemovedMap, "test_factor95percIfBrRemovedMap", "95% Conf Factor", totFactor95perc, 1d, 3d);
		}
	


//		System.out.println(infoString);
		
		//		makeLog10_EAL_Historgram();

		
//		// hard coded assignment of branch type:
//		erf_branches = new ArrayList<String>();
//		gmm_branches = new ArrayList<String>();
//
//		erf_branches.add("Fault Model");
//		erf_branches.add("Deformation Model");
//		erf_branches.add("Scaling Relationship");
//		erf_branches.add("Slip Along Rupture (Dsr)");
////		erf_branches.add("Inversion Model");
//		erf_branches.add("Total Mag>5 Rate");
//		erf_branches.add("Mmax Off Fault");
////		erf_branches.add("Moment Rate Fixes");
//		erf_branches.add("Spatial Seismicity PDF");
//		erf_branches.add("Probability Model");
//
//		gmm_branches.add("Ground Motion Model");
//		gmm_branches.add("GMM Added Uncertainty");
//		gmm_branches.add("Vs30 Model");

//		make_ERF_vs_GMM_UncertHists(true, true);
		
		if(outDirName != null) {
			try{
				FileWriter fw = new FileWriter(ROOT_DIR+outDirName+"/"+outDirName+"_Info.txt");
				fw.write(infoString+"\n");
				fw.close();
			}catch(Exception e) {
				e.printStackTrace();
			}			
		}

	}
	
	/**
	 * This reads the branch level data
	 * @return
	 */
	private void old_readBranchLevelSummaryDataFromFile() {
		String fileName = "src/scratch//ned/U3_TreeValuation/branch_level_summary.csv";
		
		String[] branchLevel, branchChoice;
		double[] totalWeight, weightedTotalMeanEAL, weightedFaultMeanEAL, weightedGriddedMeanEAL;
		double totalMeanEAL, faultMeanEAL, griddedMeanEAL;

		File file = new File(fileName);
		List<String> fileLines;
		branchLevel = new String[37];
		branchChoice = new String[37];
		totalWeight = new double[37]; 
		weightedTotalMeanEAL = new double[37]; 
		weightedFaultMeanEAL = new double[37]; 
		weightedGriddedMeanEAL = new double[37]; 
		
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			for(int i=1; i<38;i++ ) {
				String str = fileLines.get(i);
				String[] split = str.split(",");
				branchLevel[i-1] = split[0];
				branchChoice[i-1] = split[1];
				totalWeight[i-1] = Double.parseDouble(split[2]);
				weightedTotalMeanEAL[i-1] = Double.parseDouble(split[3]);
				weightedFaultMeanEAL[i-1] = Double.parseDouble(split[4]);
				weightedGriddedMeanEAL[i-1] = Double.parseDouble(split[5]);
//				System.out.println(i+"\t"+branchLevel[i-1]+"\t"+branchChoice[i-1]+"\t"+
//						totalWeight[i-1]+"\t"+weightedTotalMeanEAL[i-1]+"\t"+
//						weightedFaultMeanEAL[i-1]+"\t"+weightedGriddedMeanEAL[i-1]);
			}
			String str = fileLines.get(38);
			String[] split = str.split(",");
			totalMeanEAL= Double.parseDouble(split[3]);
			faultMeanEAL= Double.parseDouble(split[4]);
			griddedMeanEAL= Double.parseDouble(split[5]);
//			System.out.println("Totals:\t"+totalMeanEAL+"\t"+faultMeanEAL+"\t"+griddedMeanEAL);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void removeBranchFromData(String branchName, String branchValue) {
				
		HashMap<String, String[]> allBranchValuesMapNew = new HashMap<String, String[]>();
		
		ArrayList<String> optionsList = optionsForBranchHashMap.get(branchName);
		
		int oldNumOptions = optionsList.size();
		int totNumBranchesNew = totNumBranches*(oldNumOptions-1)/oldNumOptions;
		double[] branchWtNew = new double[totNumBranchesNew];
		double[] branchNormEAL_New = new double[totNumBranchesNew];
		
		optionsList.remove(branchValue);
		int newNumOptions = optionsForBranchHashMap.get(branchName).size();
//		System.out.println("oldNumOptions="+oldNumOptions+"; newNumOptions="+newNumOptions);
		
		String[] targetBranchValues = allBranchValuesMap.get(branchName);
		
		for(String brName:allBranchValuesMap.keySet()) {
			allBranchValuesMapNew.put(brName, new String[totNumBranchesNew]);
		}
		
		int indexNew=0;
		for(int i=0;i<branchNormEAL.length;i++) {
			if(!targetBranchValues[i].equals(branchValue)) {
				branchWtNew[indexNew] = branchWt[i];
				branchNormEAL_New[indexNew] = branchNormEAL[i];
				for(String brName:allBranchValuesMap.keySet()) {
					String[] stringArray = allBranchValuesMap.get(brName);
					String[] stringArrayNew = allBranchValuesMapNew.get(brName);
					stringArrayNew[indexNew] = stringArray[i];
				}
				indexNew += 1;
			}
		}
		
		// normalize new weights and EALs
		double totWt = 0;
		for(double wt:branchWtNew) totWt += wt;
		for(int i=0; i<branchWtNew.length;i++) branchWtNew[i] /= totWt;
		
		double prevMeanEAL = totMeanEAL;
		double newNormMean = computeWeightedAverage(branchWtNew, branchNormEAL_New);
		double newNormStdDev = computeWeightedStdDev(branchWtNew, branchNormEAL_New, newNormMean);
		for(int i=0; i<branchNormEAL_New.length;i++) branchNormEAL_New[i] /= newNormMean;
		totMeanEAL = newNormMean*prevMeanEAL;
		totStdDevEAL = newNormStdDev*prevMeanEAL;
		totCOV_EAL = newNormStdDev/newNormMean;
		totNumBranches=totNumBranchesNew;
		
		infoString += "Removed branch: "+branchName+" = "+branchValue+"\n";

//		infoString += "\n\n After removing branch "+branchName+" = "+branchValue+
//				":\n\n\ttotMeanEAL="+ Precision.round(totMeanEAL, 3)+
//				"\n\ttotStdDevEAL="+Precision.round(totStdDevEAL, 3)+
//				"\n\ttotCOV_EAL="+Precision.round(totCOV_EAL, 3)+
//				"\n\ttotNumBranches="+totNumBranches+"\n";
		
		allBranchValuesMap = allBranchValuesMapNew;
		branchWt = branchWtNew;
		branchNormEAL = branchNormEAL_New;
	}
	
	
	/**
	 * This reads the branch level data
	 * @return
	 */
	private void readAllBranchDataFromFile() {
		
		File file = new File(ROOT_DIR+inputFileName);
		List<String> fileLines;
		double[] branchEAL;
		
		
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			int numData = fileLines.size()-1;
			branchWt = new double[numData];
			branchEAL = new double[numData];
			branchNormEAL = new double[numData];
			allBranchValuesMap = new HashMap<String, String[]>();
			String str = fileLines.get(0);
			String[] colNamesArray = str.split(",");
			for(int c=5;c<18; c++) {
				String keyName = colNamesArray[c];
				// change name if needed
				keyName = getNameChange(keyName);
				colNamesArray[c] = keyName;
				String[] strArray = new String[numData];
				allBranchValuesMap.put(keyName, strArray);
//				infoString += "\t"+keyName+"\n";
//				System.out.println(keyName);
			}
			
			double minEAL = Double.MAX_VALUE;
			double maxEAL = 0;
			int maxLineNum=-1;
			int minLineNum=-1;
			
			totOrigWeight = 0;
			int arrayIndex = 0;
			for(int i=1; i<fileLines.size();i++ ) {
				str = fileLines.get(i);
				String[] split = str.split(",");
				branchWt[arrayIndex] = Double.parseDouble(split[1]);
				branchEAL[arrayIndex] = Double.parseDouble(split[dataCol])*1e-6;	// convert to Billions
				totOrigWeight+=branchWt[arrayIndex];
				if(minEAL > branchEAL[arrayIndex]) {
					minEAL = branchEAL[arrayIndex];
					minLineNum = i;
				}
				if(maxEAL < branchEAL[arrayIndex]) {
					maxEAL = branchEAL[arrayIndex];
					maxLineNum = i;
				}
				// set branch values
				for(int c=5;c<18; c++) {
					String keyName = colNamesArray[c];
					String[] valArray = allBranchValuesMap.get(keyName);
					String value = split[c];
					// Change name if needed
					value = getNameChange(value);
					valArray[arrayIndex] = value;
				}
				arrayIndex +=1;
			}
			
			optionsForBranchHashMap = new HashMap<String, ArrayList<String>>();
			for(String brName:allBranchValuesMap.keySet()) {
				String[] allValues = allBranchValuesMap.get(brName);
				ArrayList<String> optionsList = new ArrayList<String>();
				for(String opt:allValues)
					if(!optionsList.contains(opt))
						optionsList.add(opt);
				optionsForBranchHashMap.put(brName,optionsList);
			}
			
			infoString += "\nBranch Names/values:\n\n";
			for(String name:optionsForBranchHashMap.keySet()) 
				for(String value:optionsForBranchHashMap.get(name))
					infoString += name+" = "+value+"\n";

			
			
			
			infoString += "\nMinEAL:\t"+(float)minEAL+"; fileline:";
			infoString += "\n\t"+fileLines.get(minLineNum);
			infoString += "\nMaxEAL:\t"+(float)maxEAL+"; fileline:";
			infoString += "\n\t"+fileLines.get(maxLineNum);
			
			totMeanEAL = computeWeightedAverage(branchWt, branchEAL);
			totStdDevEAL = computeWeightedStdDev(branchWt, branchEAL, totMeanEAL);
			totCOV_EAL = totStdDevEAL/totMeanEAL;
						
			infoString += "\ntotMeanEAL = "+ Precision.round(totMeanEAL, 3);
			infoString += "\ntotStdDevEAL = "+ Precision.round(totStdDevEAL, 3);
			infoString += "\ntotCOV_EAL = "+ Precision.round(totCOV_EAL, 3);
			infoString += "\ntotOrigWeight = "+(float)totOrigWeight;
			infoString += "\nsize = "+branchWt.length;
			
			// do normalizations
			double maxNormEAL = 0;
			for(int i=0;i<branchNormEAL.length;i++) {
				branchNormEAL[i] = branchEAL[i]/totMeanEAL;
				if(branchNormEAL[i]>maxNormEAL)
					maxNormEAL = branchNormEAL[i];
				branchWt[i] /= totOrigWeight;
			}
			System.out.println("maxNormEAL = "+maxNormEAL);
			double totWt=0;
			for(double wt:branchWt) totWt+=wt;
			infoString += "\nTotal Weight (after normalization; should = 1) = "+(float)totWt;
			
			totNumBranches = branchNormEAL.length;
			
			// find the closest branch to mean (one weighted and one not)
//			int minDiffIndex = getBranchClosestToValue(1.0, false);
//			infoString += "\nBranch Closest to Mean (non weighted):\n";
//			infoString += "\tindex = "+minDiffIndex+"\n\twt = "+(float)branchWt[minDiffIndex]+"\n\tmean = "+(float)(branchNormEAL[minDiffIndex]*totMeanEAL);
//			for(String branchName:branchValuesHashMap.keySet())
//				infoString += "\n\t"+branchName+" = "+branchValuesHashMap.get(branchName)[minDiffIndex];

			int minDiffWtedIndex = getBranchClosestToValue(1.0, true);
			infoString += "\nBranch Closest to Mean (weighted):\n";
			infoString += "\tindex = "+minDiffWtedIndex+"\n\twt = "+(float)branchWt[minDiffWtedIndex]+"\n\tEAL = "+(float)(branchNormEAL[minDiffWtedIndex]*totMeanEAL);
			closestValuesForBranchMap = new HashMap<String, String> ();
			for(String branchName:allBranchValuesMap.keySet()) {
				infoString += "\n\t"+branchName+" = "+allBranchValuesMap.get(branchName)[minDiffWtedIndex];
				closestValuesForBranchMap.put(branchName, allBranchValuesMap.get(branchName)[minDiffWtedIndex]);
			}
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	/**
	 * If there are more than one with same degree of closeness, this returns the first one encountered.
	 * @param normValue - normalized EAL value
	 * @param weighted - whether to consider the branch weight in determining closeness
	 * @return
	 */
	private int getBranchClosestToValue(double normValue, boolean weighted) {
		double minDiff = Double.MAX_VALUE;
		int minDiffIndex = -1;
		for(int i=0;i<branchNormEAL.length;i++) {
			double diff = Math.abs(normValue-branchNormEAL[i]);
//			if(diff<0.0001)
//				System.out.println(fileLines.get(i+1));
			if(weighted)
				diff /= branchWt[i];
			if(minDiff>diff) {
				minDiff = diff;
				minDiffIndex=i;
			}
		}
		return minDiffIndex;
	}
	
	int nextedLoopCounter=0;
	int numCOVsForNestedLoop = -1;
	int numBranchesForCOVcalcNestedLoop = -1;
	double[] covValsFromNestedLoop = null;
	double[] covWtsFromNestedLoop = null;
	String[] brValuesFromNestedLoop = null;
	CalcProgressBar  nestedLoopProgBar;
	public void nestedCOV_Loop(ArrayList<String> brNamesList, ArrayList<ArrayList<String>> brValuesLists, int currentLevel, String[] targetNames, String[] targetValues) {
		if(currentLevel == 0) {
			targetNames = new String[brNamesList.size()];
			targetValues = new String[brNamesList.size()];
			nextedLoopCounter=0;
			numCOVsForNestedLoop = 1;
			for(ArrayList<String> optionsList : brValuesLists)
				numCOVsForNestedLoop *= optionsList.size();
			numBranchesForCOVcalcNestedLoop = totNumBranches/numCOVsForNestedLoop;
			covValsFromNestedLoop = new double[numCOVsForNestedLoop];
			covWtsFromNestedLoop = new double[numCOVsForNestedLoop];
			brValuesFromNestedLoop = new String[numCOVsForNestedLoop];
			System.out.println("numCOVsForNestedLoop="+numCOVsForNestedLoop+"\nnumBranchesForCOVcalcNestedLoop="+numBranchesForCOVcalcNestedLoop);
			nestedLoopProgBar = new CalcProgressBar("Progress", "Nested Loop");
		}
		targetNames[currentLevel] = brNamesList.get(currentLevel);
		for(String val:brValuesLists.get(currentLevel)) {
			targetValues[currentLevel] = val;
			if(currentLevel == brNamesList.size()-1) {
				// process
				double[] wtArray = new double[numBranchesForCOVcalcNestedLoop];
				double[] valArray = new double[numBranchesForCOVcalcNestedLoop];
				int counter = 0;
				for(int i=0;i<branchNormEAL.length;i++) {
					boolean keep = true;
					for(int j=0;j<targetNames.length;j++) {
						String brVal = allBranchValuesMap.get(targetNames[j])[i];
						if(!brVal.equals(targetValues[j]))
							keep = false;
					}
					if(keep) {
						wtArray[counter] = branchWt[i];
						valArray[counter] = branchNormEAL[i];
						counter += 1;
					}
				}
				double mean = computeWeightedAverage(wtArray, valArray);
				double stdDev = computeWeightedStdDev(wtArray, valArray, mean);
				covValsFromNestedLoop[nextedLoopCounter] = stdDev/mean;
				double totWt = 0;
				for(double wt:wtArray)
					totWt += wt;
				covWtsFromNestedLoop[nextedLoopCounter] = totWt;
				String valuesString = "";
				for(String opt:targetValues)
					valuesString += opt+";  ";
				brValuesFromNestedLoop[nextedLoopCounter] = valuesString;
				
				nextedLoopCounter +=1;
				nestedLoopProgBar.updateProgress(nextedLoopCounter, numCOVsForNestedLoop);

//				System.out.print(globalCounter+"; ");
//				
//				for(int i=0;i<targetNames.length;i++)
//					System.out.print(targetValues[i]+";  ");
//				System.out.print("\n");
			}
			else {
				nestedCOV_Loop(brNamesList, brValuesLists, currentLevel+1, targetNames, targetValues);
			}
		}
	}
	
	
	
	public void doCOV_ReductionAnalysisComplete(boolean popUpWindows, boolean saveResults) {

		String[] branchesToRemoveNamesArray = {
				"GMM Added Uncertainty",
				"Total Mag>5 Rate",
				"Probability Model",
				"Scaling Relationship",
				"Ground Motion Model",
				"Spatial Seismicity PDF",
				"Deformation Model",
				"Mmax Off Fault",
				"Vs30 Model",
				"Slip Along Rupture (Dsr)" //,
//				"Fault Model"
		};
		
		DefaultXY_DataSet meanCovValues = new DefaultXY_DataSet();
		DefaultXY_DataSet minCovValues = new DefaultXY_DataSet();
		DefaultXY_DataSet maxCovValues = new DefaultXY_DataSet();
		meanCovValues.set(0.0,totCOV_EAL);
		minCovValues.set(0.0,totCOV_EAL);
		maxCovValues.set(0.0,totCOV_EAL);

//		ArrayUtils.reverse(branchesToRemoveNamesArray);
		
		int numToInclude = branchesToRemoveNamesArray.length;
		
		for(int j=1;j<=numToInclude;j++) {
			System.out.println("****** Working on "+j+" (of "+numToInclude+") ******");
			ArrayList<String> brNamesList = new ArrayList<String>();
			ArrayList<ArrayList<String>> brValuesLists = new ArrayList<ArrayList<String>>();
			for(int i=0;i<j;i++) {
				brNamesList.add(branchesToRemoveNamesArray[i]);
				brValuesLists.add(optionsForBranchHashMap.get(branchesToRemoveNamesArray[i]));
			}
			
			long timeTakenMillis = System.currentTimeMillis();
			nestedCOV_Loop(brNamesList, brValuesLists, 0, null, null);
			timeTakenMillis = System.currentTimeMillis()-timeTakenMillis;
			System.out.println("Took (sec): " + timeTakenMillis/1000);
			nestedLoopProgBar.dispose();
			
//			for(int i=0;i<covValsFromNestedLoop.length;i++)
//				System.out.println((float)covValsFromNestedLoop[i]+"\t"+(float)covWtsFromNestedLoop[i]);
			
			// Check wts
			double totWt = 0;
			for(double wt:covWtsFromNestedLoop) totWt +=wt;
			System.out.println("totWt="+(float)totWt);
			double meanCOV = computeWeightedAverage(covWtsFromNestedLoop, covValsFromNestedLoop);
//			double minCOV = IEEE754rUtils.min(covValsFromNestedLoop);
//			double maxCOV = IEEE754rUtils.max(covValsFromNestedLoop);
			double minCOV = Double.MAX_VALUE;
			int minCOV_index = -1;
			double maxCOV = 0;
			int maxCOV_index = -1;

			for(int i=0;i<covValsFromNestedLoop.length;i++) {
				if(minCOV>covValsFromNestedLoop[i]) {
					minCOV=covValsFromNestedLoop[i];
					minCOV_index = i;
				}
				if(maxCOV<covValsFromNestedLoop[i]) {
					maxCOV=covValsFromNestedLoop[i];
					maxCOV_index = i;
				}
			}
				
			System.out.println("meanCOV="+(float)meanCOV);
			System.out.println("minCOV="+(float)minCOV);
			System.out.println("\t"+brValuesFromNestedLoop[minCOV_index]);
			System.out.println("maxCOV="+(float)maxCOV);
			System.out.println("\t"+brValuesFromNestedLoop[maxCOV_index]);
			
			meanCovValues.set((double)j,meanCOV);
			minCovValues.set((double)j,minCOV);
			maxCovValues.set((double)j,maxCOV);

		}
		meanCovValues.setName("meanCovValues");
		minCovValues.setName("minCovValues");
		maxCovValues.setName("maxCovValues");
		System.out.println(meanCovValues);
		System.out.println(minCovValues);
		System.out.println(maxCovValues);
		
		ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		// convert 
		for(int i=0;i<meanCovValues.size();i++) {
			DefaultXY_DataSet xyData = new DefaultXY_DataSet();
			xyData.set(meanCovValues.getX(i), 0.0);
			xyData.set(meanCovValues.getX(i), meanCovValues.getY(i));
			funcsArray.add(xyData);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		}
		
		funcsArray.add(meanCovValues);
		funcsArray.add(minCovValues);
		funcsArray.add(maxCovValues);
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.DASH, 4f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.DASH, 4f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.DASH, 4f, Color.BLACK));

		
		PlotSpec spec = new PlotSpec(funcsArray, plotChars, "", "Branch Removed", "COV");
//		specCombinedLegend.setLegendVisible(false);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);
		
		double yMax = 0.5;
		double xMax = meanCovValues.getMaxX()+0.5;

		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(spec);
			gw.setAxisRange(-0.5, xMax, 0d, yMax);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
		}

		if(saveResults) {
			String fname = "COV_ReductionAnalysisComplete";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(-0.5, xMax, 0d, yMax);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec);
			try {
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(500, 400);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");	
				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				
				FileWriter fw = new FileWriter(ROOT_DIR+outDirName+"/"+"COV_ReductionData.txt");
				for(int i=0;i<meanCovValues.size();i++)
					fw.write(meanCovValues.getX(i)+"\t"+meanCovValues.getY(i)+"\n");
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}
	
	
	
	public void makeCOV_ReductionAnalysisCompletePlotsAgain() {
		DefaultXY_DataSet meanCovValues = new DefaultXY_DataSet();
		DefaultXY_DataSet factor95conf_Values = new DefaultXY_DataSet();
		DefaultXY_DataSet fract10perc_Values = new DefaultXY_DataSet();
		
		try {
			File file = new File(ROOT_DIR+outDirName+"/"+"COV_ReductionData.txt");
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			for(String line:fileLines) {
//				System.out.println(line);
				String[] split = line.split("\t");
				double xVal = Double.parseDouble(split[0]);
				double yVal = Double.parseDouble(split[1]);
				meanCovValues.set(xVal,yVal);
				double[] valArray = computeFact95_fract10_fromLogNorm_COV(yVal);
				factor95conf_Values.set(xVal,valArray[0]);
				fract10perc_Values.set(xVal,valArray[1]);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		String[] fileNamesArray = {"COV_ReductionAnalysisCompleteAlt", "fact95_ReductionAnalysisCompleteAlt", "fract10_ReductionAnalysisCompleteAlt"};
		String[] yAxisNamesArray = {"COV", "95% Conf Factor", "Fraction Within 10% of Mean"};
		DefaultXY_DataSet[] functionsArray = {meanCovValues, factor95conf_Values, fract10perc_Values};
		double[] ymaxArray = {0.45,2.2,1.01};
		double[] yminArray = {0.0,1.0,0.0};
		
		for(int j=0;j<3;j++) {
			ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			
			// convert to sticks
			for(int i=0;i<functionsArray[j].size();i++) {
				DefaultXY_DataSet xyData = new DefaultXY_DataSet();
				xyData.set(functionsArray[j].getX(i), 0.0);
				xyData.set(functionsArray[j].getX(i), functionsArray[j].getY(i));
				funcsArray.add(xyData);
				plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
			}
					
			PlotSpec spec = new PlotSpec(funcsArray, plotChars, "", "After Branch Removed", yAxisNamesArray[j]);
			
			double xMax = functionsArray[j].getMaxX()+0.5;

			String fname = fileNamesArray[j];
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(-0.5, xMax, yminArray[j], ymaxArray[j]);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec);
			try {
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(468, 192);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
	}

	

	public void generateBranchValueResults(boolean popUpWindows, boolean saveResults) {

		double min = 0;
		double max = 7;
		int num = 140;

		String tableHeaderLine = "Branch Name\tBranch Value\tWeight\tMean (fractional change)\tCOV\t95% Conf Factor\tFract Within 10%\t"+
				"Mean (fractional change) if removed\tCOV if removed\t95% Conf Factor if removed\tFract Within 10% if removed";	
		HashMap<String, String> tableLinesMap = new HashMap<String, String>();
		
		covForBrMap = new HashMap<String, Double> ();
		wtAbsValMeanForBrMap = new HashMap<String, Double> ();


		// for branch-value averages
		weightForBrValMap = new HashMap<String, Double>();
		meanForBrValMap = new HashMap<String, Double>();
		meanDiffForBrValMap = new HashMap<String, Double>();
		meanDiffWtedForBrValMap = new HashMap<String, Double>();
		covForBrValMap = new HashMap<String, Double>();
		covDiffForBrValMap = new HashMap<String, Double>();
		fracWithIn10percForBrValMap = new HashMap<String, Double>();
		fracWithIn10percDiffForBrValMap = new HashMap<String, Double>();
		factor95percForBrValMap = new HashMap<String, Double>();
		factor95percRatioForBrValMap = new HashMap<String, Double>();

		meanIfBrRemovedMap = new HashMap<String, Double>();
		meanDiffIfBrRemovedMap = new HashMap<String, Double>();
		meanDiffWtedIfBrRemovedMap = new HashMap<String, Double>();
		covIfBrRemovedMap = new HashMap<String, Double>();
		covDiffIfBrRemovedMap = new HashMap<String, Double>();
		fracWithIn10percIfBrRemovedMap = new HashMap<String, Double>();
		fracWithIn10percDiffIfBrRemovedMap = new HashMap<String, Double>();
		factor95percIfBrRemovedMap = new HashMap<String, Double>();
		mfactor95percRatioDiffIfBrRemovedMap = new HashMap<String, Double>();

		infoString += "\n\nExpected COV etc if Branches Removed:\n\n" ;
		infoString += "name\tcov\tcovRatio\tfact95perc\twithin10perc\twtAveAbsValMeanDiff\n";

		HashMap<String, Double> infoLinesMap = new HashMap<String, Double>();
		
		for(String branchName:allBranchValuesMap.keySet()) {
			//System.out.println("working on" + branchName);
			String[] branchValuesArray = allBranchValuesMap.get(branchName);
			ArrayList<String> optionsForBranchList = optionsForBranchHashMap.get(branchName);
			int numOptions = optionsForBranchList.size();
			if(numOptions == 1)
				continue;
			int numForBranch = totNumBranches/numOptions;
			int numOtherBranches = totNumBranches - numForBranch;

			// for plots
			HashMap<String, HistogramFunction> histogramsHashMap = new HashMap<String, HistogramFunction>();
			HashMap<String, Double> meanHashMap = new HashMap<String, Double>();		
			HashMap<String, Double> wtHashMap = new HashMap<String, Double>();	

			double wtedAbsValOfMeanDiffs = 0;
			double expCOV_ifBranchesRemoved = 0;
			double expFactor95percIfBranchesRemoved = 0;
			double expwithin10percIfBranchesRemoved = 0;

			for(String brOpt:optionsForBranchList) {
				double[] wtForBranchArray = new double[numForBranch];
				double[] normEAL_ForBranchArray = new double[numForBranch];
				double[] wtForOtherBranchesArray = new double[numOtherBranches];
				double[] normEAL_ForOtherBranchesArray = new double[numOtherBranches];
				HistogramFunction histForBr = new HistogramFunction(min+(max-min)/(num*2),max-(max-min)/(num*2), num);
				histForBr.setName(brOpt);
				HistogramFunction histForOtherBrs = new HistogramFunction(min+(max-min)/(num*2),max-(max-min)/(num*2), num);
				histForOtherBrs.setName("histForOtherBrs: "+brOpt);

				int indexForBranch = 0;
				int indexForOtherBranches = 0;
				for(int i=0; i<branchNormEAL.length;i++) {
					String brVal = branchValuesArray[i];
					if(brOpt.equals(brVal)) {
						wtForBranchArray[indexForBranch] = branchWt[i];
						normEAL_ForBranchArray[indexForBranch] = branchNormEAL[i];
						histForBr.add(branchNormEAL[i], branchWt[i]);
						indexForBranch+=1;
					}
					else {
						wtForOtherBranchesArray[indexForOtherBranches] = branchWt[i]; 
						normEAL_ForOtherBranchesArray[indexForOtherBranches] = branchNormEAL[i]; 
						histForOtherBrs.add(branchNormEAL[i], branchWt[i]);
						indexForOtherBranches+=1;

					}					
				}
				// results for branch
				double wtForBr =0;
				for(double wt:wtForBranchArray) wtForBr += wt;
				double testWt = histForBr.calcSumOfY_Vals();
				if(Math.abs(wtForBr-testWt) > 0.0001)
					throw new RuntimeException(brOpt+"\tSomething wrong: "+wtForBr+"\t"+testWt);
				double meanForBr = computeWeightedAverage(wtForBranchArray, normEAL_ForBranchArray);
				double stdDevForBr = computeWeightedStdDev(wtForBranchArray, normEAL_ForBranchArray, meanForBr);
				double covForBr = stdDevForBr/meanForBr;
				HistogramFunction histForBrCumulative = histForBr.getCumulativeDistFunctionWithHalfBinOffset();
				histForBrCumulative.scale(1.0/histForBrCumulative.getMaxY());; // so it sums to 1.0
				double factor95percForBr = get95percConfFactorForValue(meanForBr, histForBrCumulative);
				double within10percForBr = histForBrCumulative.getInterpolatedY(meanForBr*1.1) - histForBrCumulative.getInterpolatedY(meanForBr*0.9);

				wtedAbsValOfMeanDiffs += Math.abs(meanForBr-1.0)*wtForBr;
				expCOV_ifBranchesRemoved += covForBr*wtForBr;
				expFactor95percIfBranchesRemoved += factor95percForBr*wtForBr;
				expwithin10percIfBranchesRemoved += within10percForBr*wtForBr;

				// for stacked histogram plotting
				wtHashMap.put(brOpt, wtForBr);
				meanHashMap.put(brOpt, meanForBr);
				histForBr.scale(1.0/histForBr.getDelta()); // scale to density
				histogramsHashMap.put(brOpt, histForBr);

				// results for other branches
				double wtForOtherBrs =0;
				for(double wt:wtForOtherBranchesArray) wtForOtherBrs += wt;
				testWt = histForOtherBrs.calcSumOfY_Vals();
				if(Math.abs(wtForOtherBrs-testWt) > 0.0001)
					throw new RuntimeException("Something wrong");
				double meanForOtherBrs = computeWeightedAverage(wtForOtherBranchesArray, normEAL_ForOtherBranchesArray);
				double stdDevForOtherBrs = computeWeightedStdDev(wtForOtherBranchesArray, normEAL_ForOtherBranchesArray, meanForOtherBrs);
				double covForOtherBrs = stdDevForOtherBrs/meanForOtherBrs;
				HistogramFunction histForOtherBrsCumulative = histForOtherBrs.getCumulativeDistFunctionWithHalfBinOffset();
				histForOtherBrsCumulative.scale(1.0/histForOtherBrsCumulative.getMaxY());; // so it sums to 1.0
				double factor95percForOtherBrs = get95percConfFactorForValue(meanForOtherBrs, histForOtherBrsCumulative);
				double within10percForOtherBrs = histForOtherBrsCumulative.getInterpolatedY(meanForOtherBrs*1.1) - histForOtherBrsCumulative.getInterpolatedY(meanForOtherBrs*0.9);

				//				if(brOpt.equals("NONE")) {
				//					GraphWindow graph = new GraphWindow(histForOtherBrs, "");
				//				}

				//				"branchName\tbranchValue\tweight\tmean\tmeanDiff\twtedMeanDiff\tmeanIfRemoved\tmeanIfRemovedDiff";	
				String tableLine = branchName + "\t" +
						brOpt + "\t"+
						(float)wtForBr + "\t" +
//						(float)meanForBr + "\t" +
						(float)(meanForBr-1.0) + "\t" +
//						(float)((meanForBr-1.0)*wtForBr) + "\t"+
						(float)covForBr + "\t" +
						(float)factor95percForBr + "\t" +
						(float)within10percForBr + "\t" +

//						(float)meanForOtherBrs + "\t" +
						(float)(meanForOtherBrs-1.0) + "\t" +
						(float)covForOtherBrs + "\t" +
						(float)factor95percForOtherBrs + "\t" +
						(float)within10percForOtherBrs + "\t";
//				tableLinesList.add(tableLine);

				String combinedName = branchName+" = "+brOpt;
				
				tableLinesMap.put(combinedName, tableLine);
				
				weightForBrValMap.put(combinedName, wtForBr);
				meanForBrValMap.put(combinedName, meanForBr);
				meanDiffForBrValMap.put(combinedName, meanForBr-1.0);
				meanDiffWtedForBrValMap.put(combinedName, (meanForBr-1.0)*wtForBr);
				meanIfBrRemovedMap.put(combinedName, meanForOtherBrs);
				meanDiffIfBrRemovedMap.put(combinedName, (meanForOtherBrs-1.0));
				meanDiffWtedIfBrRemovedMap.put(combinedName, ((1.0-wtForBr)*(meanForOtherBrs-1.0)));
				covForBrValMap.put(combinedName, covForBr);
				covIfBrRemovedMap.put(combinedName, covForOtherBrs);
				fracWithIn10percForBrValMap.put(combinedName, within10percForBr);
				fracWithIn10percIfBrRemovedMap.put(combinedName, within10percForOtherBrs);
				factor95percForBrValMap.put(combinedName, factor95percForBr);
				factor95percIfBrRemovedMap.put(combinedName, factor95percForOtherBrs);

			}
			
			String lineSting = branchName+"\t"+(float)expCOV_ifBranchesRemoved+"\t"+(float)(expCOV_ifBranchesRemoved/totCOV_EAL)+
					"\t"+(float)expFactor95percIfBranchesRemoved+"\t"+(float)expwithin10percIfBranchesRemoved+"\t"+(float)wtedAbsValOfMeanDiffs+"\n";
			infoLinesMap.put(lineSting, expCOV_ifBranchesRemoved);
			
			covForBrMap.put(branchName, expCOV_ifBranchesRemoved);
			wtAbsValMeanForBrMap.put(branchName, wtedAbsValOfMeanDiffs);


//			infoString += branchName+"\t"+(float)expCOV_ifBranchesRemoved+"\t"+(float)(expCOV_ifBranchesRemoved/totCOV_EAL)+
//					"\t"+(float)expFactor95percIfBranchesRemoved+"\t"+(float)expwithin10percIfBranchesRemoved+"\n";

			if(saveResults || popUpWindows)
				makeStackedHistograms(branchName, histogramsHashMap, meanHashMap, wtHashMap, popUpWindows, saveResults);


			//			for(String branchValue: meanHashMap.keySet()) {
			//				String combinedName = branchName+" = "+branchValue;
			//				weightForBrValMap.put(combinedName, weightHashMap.get(branchValue));
			//				meanForBrValMap.put(combinedName, meanHashMap.get(branchValue));
			//				meanDiffForBrValMap.put(combinedName, meanDiffHashMap.get(branchValue));
			//				meanDiffWtedForBrValMap.put(combinedName, wtedMeanDiffHashMap.get(branchValue));
			//				meanIfBrRemovedMap.put(combinedName, meanIfRemovedHashMap.get(branchValue));
			//				meanDiffIfBrRemovedMap.put(combinedName, meanIfRemovedDiffHashMap.get(branchValue));
			//			}

		}

		//		for(String line:	 tableLinesList)		
		//			System.out.println(line);
		
		List<String> sortedInfoLinesList = getKeysSortedByMapValues(infoLinesMap);
		for(String line:sortedInfoLinesList)
			infoString+=line;
		
//		nextBranchesToRemoveToMinimizeCOV_List = new ArrayList<String>();
//		String[] splitLine = sortedInfoLinesList.get(0).split("\t");
//		String nextBranchName = splitLine[0];
//		double meanCOV = Double.parseDouble(splitLine[1]);
////System.out.println("nextBranchName="+nextBranchName+"\tmeanCOV="+meanCOV);
//
//		double minDiff = Double.MAX_VALUE;
//		String minCombinedName = null;
//		for(String opt:optionsForBranchHashMap.get(nextBranchName)) {
//			String combinedName = nextBranchName+" = "+opt;
//			nextBranchesToRemoveToMinimizeCOV_List.add(combinedName);
//			double cov = covForBrValMap.get(combinedName);
//			double wt = weightForBrValMap.get(combinedName);
//			double absDiff = Math.abs(cov-meanCOV)/wt;
//			if(minDiff>absDiff) {
//				minDiff = absDiff;
//				minCombinedName = combinedName;
//			}
////System.out.println(combinedName+"\t"+cov);
//
//		}
//		nextBranchesToRemoveToMinimizeCOV_List.remove(minCombinedName); // remove the one to keep.
//		
//System.out.println("branch to set for next:\t"+minCombinedName);
////System.out.println("nextBranchesToRemoveToMinimizeCOV_List:\n"+nextBranchesToRemoveToMinimizeCOV_List);
////System.exit(0);
		
		ArrayList<String> combineNameList = getLogicTreeBranchValueOrder();

		if(saveResults) {
			try{
				FileWriter fw = new FileWriter(ROOT_DIR+outDirName+"/"+outDirName+"_BrValStats.txt");
				fw.write(tableHeaderLine+"\n");

				for(String combinedName:	 combineNameList)		
					fw.write(tableLinesMap.get(combinedName)+"\n");
				fw.close();
			}catch(Exception e) {
				e.printStackTrace();
			}			
		}

	}

	
	/**
	 * This makes a list of branch options minus the one that has a COV that is closest 
	 * to the wt-average COV (the one to keep)
	 * @param branchName
	 * @return
	 */
	private ArrayList<String> getBranchesToRemoveForCOV_List(String branchName) {
		
		ArrayList<String> nextBranchesToRemoveToMinimizeCOV_List = new ArrayList<String>();

		double wtMeanCOV = 0d;
		double totWt = 0d;
		for(String opt:optionsForBranchHashMap.get(branchName)) {
			String combinedName = branchName+" = "+opt;
			double wt = weightForBrValMap.get(combinedName);
			wtMeanCOV += wt*covForBrValMap.get(combinedName);
			totWt += wt;
		}
		if(Math.abs(totWt-1.0)>0.00001)
			throw new RuntimeException("Problem");
		
		double minDiff = Double.MAX_VALUE;
		String minCombinedName = null;
		for(String opt:optionsForBranchHashMap.get(branchName)) {
			String combinedName = branchName+" = "+opt;
			nextBranchesToRemoveToMinimizeCOV_List.add(combinedName);
			double cov = covForBrValMap.get(combinedName);
			double wt = weightForBrValMap.get(combinedName);
			double absDiff = Math.abs(cov-wtMeanCOV)/wt;
			if(minDiff>absDiff) {
				minDiff = absDiff;
				minCombinedName = combinedName;
			}
System.out.println("\t"+(float)wtMeanCOV+"\t"+(float)cov+"\t"+(float)wt+"\t"+(float)absDiff+"\t"+opt);
		}
		nextBranchesToRemoveToMinimizeCOV_List.remove(minCombinedName); // remove the one to keep.
		
System.out.println("branch to set for next:\t"+minCombinedName);
//System.out.println("nextBranchesToRemoveToMinimizeCOV_List:\n"+nextBranchesToRemoveToMinimizeCOV_List);
//System.exit(0);

		return nextBranchesToRemoveToMinimizeCOV_List;
	}

	private void makeStackedHistograms(String branchName, HashMap<String, HistogramFunction> histogramsHashMap, HashMap<String, Double> meanHashMap, 
			HashMap<String, Double> weightHashMap, boolean popUpWindows, boolean saveResults) {
		
		File stackedHistDir = null;
		if(saveResults) {
			stackedHistDir = new File(ROOT_DIR+outDirName+"/stackedHistograms");
			if(!stackedHistDir.exists())
				stackedHistDir.mkdirs();
		}
				
		double xMinPlot=0;
		double xMaxPlot=3;
		
		ArrayList<HistogramFunction> histArrayList = new ArrayList<HistogramFunction>();
		for(String branchValue:meanHashMap.keySet())
			histArrayList.add(histogramsHashMap.get(branchValue));
		List<HistogramFunction> stackedHistograms = HistogramFunction.getStackedHists(histArrayList, false);
		
		HashMap<String, Color> colorMap = getColorMap(meanHashMap);
		ArrayList<PlotCurveCharacterstics> plotCharsHists = new ArrayList<PlotCurveCharacterstics>();
		ArrayList<PlotCurveCharacterstics> plotCharsLines = new ArrayList<PlotCurveCharacterstics>();
		ArrayList<PlotCurveCharacterstics> plotCharsValuePDFs = new ArrayList<PlotCurveCharacterstics>();

		for(String branchOption: meanHashMap.keySet()) {
			plotCharsHists.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colorMap.get(branchOption)));
			plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colorMap.get(branchOption)));
			plotCharsValuePDFs.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colorMap.get(branchOption)));
		}

		
		ArrayList<DefaultXY_DataSet> meanLinesList = new ArrayList<DefaultXY_DataSet>();
		double totMeanEAL_here=0;
		// the following are to plot the mean lines offset from the graph
		double maxYhist = stackedHistograms.get(0).getMaxY();
		maxYhist = Math.ceil(maxYhist*100)/100.0;
		double yPlotMax = 1.3*maxYhist;
		double yOffset=maxYhist*1.03;
		double yScale=(yPlotMax-yOffset)/0.8;	// 0.8 is the max branch weight (MMax)
		double minMean = Double.MAX_VALUE;
		double maxMean = 0;
		for(String branchValue: histogramsHashMap.keySet()) {
			HistogramFunction hist = histogramsHashMap.get(branchValue);
			DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
			double mean = meanHashMap.get(branchValue);
			double wt = weightHashMap.get(branchValue);
			meanLine.set(mean,0.0+yOffset);
			meanLine.set(mean,(0.0+yOffset)+wt*yScale);
			meanLine.setName(hist.getName());
			double meanFromHist = hist.computeMean();
			double stdevFromHist = hist.computeStdDev();
			totMeanEAL_here += mean*wt;

			meanLine.setInfo("mean="+(float)mean+"\nweight="+(float)wt+"\tmeanFromHist="+(float)meanFromHist+
					"\nstdevFromHist="+(float)stdevFromHist);
			meanLinesList.add(meanLine);
			if(minMean>mean)
				minMean=mean;
			if(maxMean<mean)
				maxMean=mean;
//			System.out.println(branchName+"\t"+branchValue+"\t"+(float)wt+"\t"+(float)mean+"\t"+meanFromHist+"\t"+stdevFromHist);
		}
		
		DefaultXY_DataSet totMeanLine = new DefaultXY_DataSet();
		totMeanLine.set(totMeanEAL_here,0.0);
		totMeanLine.set(totMeanEAL_here,yPlotMax);
		totMeanLine.setName("Mean");
		totMeanLine.setInfo("totMean="+(float)totMeanEAL_here);
		
		DefaultXY_DataSet meanLinesPlatform = new DefaultXY_DataSet();
		meanLinesPlatform.set(minMean-0.005,yOffset);
		meanLinesPlatform.set(maxMean+0.005,yOffset);
		meanLinesPlatform.setName("meanLinesPlatform");

					
		// this is for combined stack hists and means
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		funcs.addAll(stackedHistograms);
		funcs.add(totMeanLine);
		funcs.addAll(meanLinesList);
//		funcs.addAll(histArrayList);
//		funcs.add(meanLinesPlatform);

		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.addAll(plotCharsHists);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.addAll(plotCharsLines);
//		plotChars.addAll(plotCharsValuePDFs);
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));

		PlotSpec specCombined = new PlotSpec(funcs, plotChars, branchName, "Normalized Loss", "Density");
		specCombined.setLegendVisible(false);
		specCombined.setLegendLocation(RectangleEdge.RIGHT);
		
		if(stackedHistDir != null) {
			String fname = branchName.replace(" ", "_");
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xMinPlot, xMaxPlot, 0d, yPlotMax);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombined);
			try {
//				File file = new File(dirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				File file = new File(stackedHistDir, fname);
				gp.getChartPanel().setSize(400, 300);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
		// now with legend
		PlotSpec specCombinedLegend = new PlotSpec(funcs, plotChars, branchName, "Normalized Loss", "Density");
		specCombinedLegend.setLegendVisible(true);
		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);

		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(specCombinedLegend);
//			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			gw.setAxisRange(xMinPlot, xMaxPlot, 0d, yPlotMax);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
//			gw.setBackgroundColor(Color.WHITE);	
		}

		if(stackedHistDir != null) {
			String fname = branchName.replace(" ", "_")+"_wLegend";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xMinPlot, xMaxPlot, 0d, yPlotMax);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombinedLegend);
			try {
//				File file = new File(dirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				File file = new File(stackedHistDir, fname);
				gp.getChartPanel().setSize(500, 300);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		
		// now plot histograms
		ArrayList<XY_DataSet> funcs2 = new ArrayList<XY_DataSet>();
		funcs2.addAll(histArrayList);
		funcs2.add(totMeanLine);
		funcs2.addAll(meanLinesList);
		funcs2.add(stackedHistograms.get(0));


		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		plotChars2.addAll(plotCharsLines);
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.BLACK));
		plotChars2.addAll(plotCharsLines);
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));


		PlotSpec specCombined2 = new PlotSpec(funcs2, plotChars2, branchName, "Normalized AAL", "Density");
		specCombined2.setLegendVisible(false);
		specCombined2.setLegendLocation(RectangleEdge.RIGHT);
		
		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(specCombined2);
//			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			gw.setAxisRange(xMinPlot, xMaxPlot, 0d, yPlotMax);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
//			gw.setBackgroundColor(Color.WHITE);	
		}

		
		if(stackedHistDir != null) {
			String fname = branchName.replace(" ", "_");
			fname += "_PDFs";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xMinPlot, xMaxPlot, 0d, yPlotMax);
//			gp.setTickLabelFontSize(23);
//			gp.setAxisLabelFontSize(24);
//			gp.setPlotLabelFontSize(26);
			gp.setTickLabelFontSize(9);
			gp.setAxisLabelFontSize(10);
			gp.setPlotLabelFontSize(10);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombined2);
			try {
//				File file = new File(dirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				File file = new File(stackedHistDir, fname);
				gp.getChartPanel().setSize(400, 300);
				gp.getChartPanel().setSize(180, 155);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
	}


	
	
	
	
	
	
	

	

	
	public void make_ERF_vs_GMM_UncertHists(boolean popUpWindows, boolean saveResults) {
		
		double min = 0;
		double max = 6;
		int num = 60;
		
		// make ERF histogram
		double mean_erf =0;
		double totWt_erf=0;
		HistogramFunction hist_erf = new HistogramFunction(0.0+max/(num*2),max-max/(num*2), num);
		for(int i=0; i<branchNormEAL.length; i++) {
				// use only gmm branches that are the closest one to the mean
				boolean useThis = true;
				for(String gmmBranchName:gmm_branches) {
					String branchValue = allBranchValuesMap.get(gmmBranchName)[i];
					String closestValue = closestValuesForBranchMap.get(gmmBranchName);
					if(!branchValue.equals(closestValue))
						useThis = false;
				}
				if(useThis) {
					hist_erf.add(branchNormEAL[i], branchWt[i]);
					mean_erf += branchNormEAL[i]*branchWt[i];
					totWt_erf += branchWt[i];						
				}
		}
		mean_erf /= totWt_erf;
		hist_erf.scale(1.0/(hist_erf.getDelta()*hist_erf.calcSumOfY_Vals()));
		double stdDev_erf = hist_erf.computeStdDev();
		HistogramFunction cumDist_erf = hist_erf.getCumulativeDistFunctionWithHalfBinOffset();
		double low95_erf = cumDist_erf.getFirstInterpolatedX(0.025);
		double upp95_erf = cumDist_erf.getFirstInterpolatedX(0.975);
		double diff95_erf = upp95_erf-low95_erf;
		hist_erf.setName("hist_erf");
		String info_erf = "mean="+(float)+mean_erf+"\nstdDev="+(float)stdDev_erf+"\ntotWt="+(float)totWt_erf+
				"\n95% bounds (diff): "+(float)low95_erf +", "+ (float)upp95_erf +", ("+diff95_erf+")";
		hist_erf.setInfo(info_erf);

				
		// make GMM histogram
		double mean_gmm = 0;
		double totWt_gmm = 0;
		HistogramFunction hist_gmm = new HistogramFunction(0.0+max/(num*2),max-max/(num*2), num);
		for(int i=0; i<branchNormEAL.length; i++) {
			// use only erf branches that are the closest one to the mean
				boolean useThis = true;
				for(String erfBranchName:erf_branches) {
					String branchValue = allBranchValuesMap.get(erfBranchName)[i];
					String closestValue = this.closestValuesForBranchMap.get(erfBranchName);
					if(!branchValue.equals(closestValue))
						useThis = false;
				}
				if(useThis) {
					hist_gmm.add(branchNormEAL[i], branchWt[i]);
					mean_gmm += branchNormEAL[i]*branchWt[i];
					totWt_gmm += branchWt[i];						
				}
		}
		mean_gmm /= totWt_gmm;
		hist_gmm.scale(1.0/(hist_gmm.getDelta()*hist_gmm.calcSumOfY_Vals()));
		double stdDev_gmm = hist_gmm.computeStdDev();
		HistogramFunction cumDist_gmm = hist_gmm.getCumulativeDistFunctionWithHalfBinOffset();
		double low95_gmm = cumDist_gmm.getFirstInterpolatedX(0.025);
		double upp95_gmm = cumDist_gmm.getFirstInterpolatedX(0.975);
		double diff95_gmm = upp95_gmm-low95_gmm;
		
		hist_gmm.setName("hist_gmm");
		String info_gmm = "mean="+(float)+mean_gmm+"\nstdDev="+(float)stdDev_gmm+"\ntotWt="+(float)totWt_gmm+
				"\n95% bounds (diff): "+(float)low95_gmm +", "+ (float)upp95_gmm +", ("+diff95_gmm+")";
		hist_gmm.setInfo(info_gmm);

		ArrayList<HistogramFunction> histArrayList = new ArrayList<HistogramFunction>();
		histArrayList.add(hist_erf);
		histArrayList.add(hist_gmm);
		
		ArrayList<PlotCurveCharacterstics> plotCharsLines = new ArrayList<PlotCurveCharacterstics>();
		plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));

		
		PlotSpec specCombinedLegend = new PlotSpec(histArrayList, plotCharsLines, "ERF vs GMM Dist", "EAL (Billion $)", "Density");
		specCombinedLegend.setLegendVisible(true);
		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);

		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(specCombinedLegend);
//			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
//			gw.setAxisRange(xMinPlot, xMaxPlot, 0d, yPlotMax);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
//			gw.setBackgroundColor(Color.WHITE);	
		}

		if(saveResults) {
			String fname = "ERF_vs_GMM_EAL_Dists";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
//			gp.setUserBounds(xMinPlot, xMaxPlot, 0d, yPlotMax);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombinedLegend);
			try {
//				File file = new File(dirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(500, 300);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		
		
		
		
	}
	
	/**
	 * this lists the combined name-value options in the order that they are shown in the
	 * logic-tree figure
	 * @return
	 */
	private ArrayList<String> getLogicTreeBranchValueOrder() {
		ArrayList<String> list = new ArrayList<String>();
		list.add("Fault Model = FM3.1");
		list.add("Fault Model = FM3.2");
		list.add("Deformation Model = GEOL");
		list.add("Deformation Model = ABM");
		list.add("Deformation Model = NEOK");
		list.add("Deformation Model = ZENG");
		list.add("Scaling Relationship = Shaw09Mod");
		list.add("Scaling Relationship = EllB");
		list.add("Scaling Relationship = HB08");
		list.add("Scaling Relationship = EllBsqrtLen");
		list.add("Scaling Relationship = ShConStrDrp");
		list.add("Slip Along Rupture (Dsr) = Tapered");
		list.add("Slip Along Rupture (Dsr) = Uniform");
		list.add("Total Mag>5 Rate = 6.5");
		list.add("Total Mag>5 Rate = 7.9");
		list.add("Total Mag>5 Rate = 9.6");
		list.add("Mmax Off Fault = 7.3");
		list.add("Mmax Off Fault = 7.6");
		list.add("Mmax Off Fault = 7.9");
		list.add("Spatial Seismicity PDF = U2");
		list.add("Spatial Seismicity PDF = U3");
		list.add("Probability Model = LOW_VALUES");
		list.add("Probability Model = MID_VALUES");
		list.add("Probability Model = HIGH_VALUES");
		list.add("Probability Model = POISSON");
		list.add("Vs30 Model = WillsEtAl");
		list.add("Vs30 Model = WaldAllen");
		list.add("Ground Motion Model = ASK_2014");
		list.add("Ground Motion Model = BSSA_2014");
		list.add("Ground Motion Model = CB_2014");
		list.add("Ground Motion Model = CY_2014");
		list.add("Ground Motion Model = IDRISS_2014");
		list.add("GMM Added Uncertainty = LOWER");
		list.add("GMM Added Uncertainty = NONE");
		list.add("GMM Added Uncertainty = UPPER");
		return list;
	}
	
	private static HashMap<String, Color> getColorForBranchNameMap() {
		HashMap<String, Color> colorMap = new HashMap<String, Color>();
		colorMap.put("Fault Model", Color.BLACK);
		colorMap.put("Deformation Model", new Color(0, 160, 0));
		colorMap.put("Scaling Relationship", new Color(0, 0, 255));
		colorMap.put("Slip Along Rupture (Dsr)", new Color(100, 100, 255));
		colorMap.put("Total Mag>5 Rate", new Color(0, 0, 200));
		colorMap.put("Mmax Off Fault", new Color(128, 0, 255));
		colorMap.put("Spatial Seismicity PDF", new Color(0, 128, 255));
		colorMap.put("Probability Model", new Color(255, 0, 0));
		colorMap.put("Vs30 Model", new Color(200, 150, 0));
		colorMap.put("Ground Motion Model", new Color(220, 0, 220));
		colorMap.put("GMM Added Uncertainty", new Color(255, 147, 0));
		return colorMap;
	}
	
	
	
	public void makeDiffForEachBranchOptionPlot(boolean popUpWindows, boolean savePlots, HashMap<String, Double> fractChangeMap, String fileName, 
			String xAxisName, double anchorX, double minX, double maxX) {
		
		HashMap<String, Color> colorMap = getColorForBranchNameMap();
		
//		System.out.println("Unsorted list:\n");
//		for(String name:fractChangeMap.keySet())
//			System.out.println("\t"+name+"\t"+fractChangeMap.get(name));
		
		//convert map to a List
		List<Entry<String, Double>> list = new LinkedList<Map.Entry<String, Double>>(fractChangeMap.entrySet());

		//sorting the list with a comparator
		Collections.sort(list, new Comparator<Entry<String, Double>>() {
			public int compare(Map.Entry<String, Double> o1, Map.Entry<String, Double> o2) {
				return (o1.getValue()).compareTo(o2.getValue());
			}
		});		


		//convert sorted list ArrayLists
		ArrayList<String> sortedNamesList = new ArrayList<String>();
		ArrayList<Double> sortedVauesList = new ArrayList<Double>();
//		System.out.println("\nSorted list:\n");
		for (Entry<String, Double> entry : list) {
			sortedNamesList.add(entry.getKey());
			sortedVauesList.add(entry.getValue());
//			System.out.println("\t"+entry.getValue().floatValue()+"\t"+entry.getKey());
		}
		
//		double maxDiff = Math.max(-sortedVauesList.get(0), sortedVauesList.get(sortedVauesList.size()-1));
//		double xMax = Math.ceil(maxDiff*10.0)/10.0;
		
		// OVERRIDE TO ORDER IN LOGIC_TREE FIGURE
		sortedNamesList = getLogicTreeBranchValueOrder();
		Collections.reverse(sortedNamesList);
		sortedVauesList = new ArrayList<Double>();
		for(String nameVal:sortedNamesList)
			sortedVauesList.add(fractChangeMap.get(nameVal));
		
//		// OVERRIDE TO TORNADO TYPE PLOT:
//		int lowIndex = 0;
//		int highIndex = sortedNamesList.size()-1;
//		ArrayList<String> sortedNamesListTornado = new ArrayList<String>();
//		ArrayList<Double> sortedVauesListTornado = new ArrayList<Double>();
//
//		while(lowIndex<=highIndex) {
//			sortedNamesListTornado.add(sortedNamesList.get(highIndex));
//			sortedVauesListTornado.add(sortedVauesList.get(highIndex));
//			if(lowIndex != highIndex) {
//				sortedNamesListTornado.add(sortedNamesList.get(lowIndex));
//				sortedVauesListTornado.add(sortedVauesList.get(lowIndex));
//			}
//			lowIndex +=1;
//			highIndex -= 1;
//		}
//		sortedNamesList = sortedNamesListTornado;
//		sortedVauesList = sortedVauesListTornado;
//		Collections.reverse(sortedNamesList);
//		Collections.reverse(sortedVauesList);
		
		ArrayList<DefaultXY_DataSet> funcList = new ArrayList<DefaultXY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		for(int i=0;i<sortedNamesList.size();i++) {
			DefaultXY_DataSet xyFunc = new DefaultXY_DataSet();
			double yVal = i+0.5;
			xyFunc.set(sortedVauesList.get(i), yVal);
			xyFunc.set(anchorX, yVal);
			xyFunc.setName(sortedNamesList.get(i));
			funcList.add(xyFunc);
			String[] branchName = sortedNamesList.get(i).split(" = ");
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colorMap.get(branchName[0])));
		}
		
		String stringForInfo="";
		for(int i=sortedNamesList.size()-1;i>=0;i--) {
			stringForInfo += sortedNamesList.get(i)+"\n";
		}
		funcList.get(0).setInfo("PLOT ORDER:\n"+stringForInfo);

		
		double xMinPlot = minX;
		double xMaxPlot = maxX;
		double yMinPlot = 0;
		double yMaxPlot = sortedNamesList.size();

		//add 10% change lines if anchor is 0.0 (fractional change plot)
		if(anchorX == 0.0) {
			DefaultXY_DataSet tenPercLowFunc = new DefaultXY_DataSet();
			tenPercLowFunc.set(-0.1, 0.0);
			tenPercLowFunc.set(-0.1, yMaxPlot);
			DefaultXY_DataSet tenPercHighFunc = new DefaultXY_DataSet();
			tenPercHighFunc.set(0.1, 0.0);
			tenPercHighFunc.set(0.1, yMaxPlot);
			funcList.add(tenPercLowFunc);
			funcList.add(tenPercHighFunc);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));			
		}
		else {	// add anchor line
			DefaultXY_DataSet anchorFunc = new DefaultXY_DataSet();
			anchorFunc.set(anchorX, 0.0);
			anchorFunc.set(anchorX, yMaxPlot);
			anchorFunc.setName("anchorFunc");
			funcList.add(anchorFunc);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));

		}

		PlotSpec specCombined = new PlotSpec(funcList, plotChars, null, xAxisName, null);
		specCombined.setLegendVisible(false);
		specCombined.setLegendLocation(RectangleEdge.RIGHT);
		
		if(savePlots) {
			String fname = fileName;
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(22);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombined);
			try {
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(400, 500);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(specCombined);
//			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			gw.setAxisRange(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
//			gw.setBackgroundColor(Color.WHITE);	
		}

		
//		// now with legend
//		PlotSpec specCombinedLegend = new PlotSpec(funcList, plotChars, null, "Fractional Change", "Branch");
//		specCombinedLegend.setLegendVisible(true);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);
//
//		if(popUpWindows) {
//			GraphWindow gw = new GraphWindow(specCombinedLegend);
////			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
//			gw.setAxisRange(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
//			gw.setTickLabelFontSize(19);
//			gw.setAxisLabelFontSize(20);
//			gw.setPlotLabelFontSize(21);
////			gw.setBackgroundColor(Color.WHITE);	
//		}
//
//		if(savePlots) {
//			String fname = fileName+"_wLegend";
//			HeadlessGraphPanel gp = new HeadlessGraphPanel();
//			gp.setUserBounds(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
//			gp.setTickLabelFontSize(23);
//			gp.setAxisLabelFontSize(24);
//			gp.setPlotLabelFontSize(22);
//			gp.setBackgroundColor(Color.WHITE);
//			gp.drawGraphPanel(specCombinedLegend);
//			try {
//				File file = new File(ROOT_DIR+outDirName, fname);
//				gp.getChartPanel().setSize(800, 500);
//				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
////				gp.saveAsPNG(file.getAbsolutePath()+".png");				
//
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}
	}
	
	
	
	public void makeDiffForEachBranchPlot(boolean popUpWindows, boolean savePlots, HashMap<String, Double> fractChangeMap, String fileName, 
			String xAxisName, double anchorX, double minX, double maxX) {
		
		HashMap<String, Color> colorMap = getColorForBranchNameMap();
		
		//convert map to a List
		List<Entry<String, Double>> list = new LinkedList<Map.Entry<String, Double>>(fractChangeMap.entrySet());

		//sorting the list with a comparator
		Collections.sort(list, new Comparator<Entry<String, Double>>() {
			public int compare(Map.Entry<String, Double> o1, Map.Entry<String, Double> o2) {
				return (o1.getValue()).compareTo(o2.getValue());
			}
		});		


		//convert sorted list ArrayLists
		ArrayList<String> sortedNamesList = new ArrayList<String>();
		ArrayList<Double> sortedVauesList = new ArrayList<Double>();
//		System.out.println("\nSorted list:\n");
		for (Entry<String, Double> entry : list) {
			sortedNamesList.add(entry.getKey());
			sortedVauesList.add(entry.getValue());
//			System.out.println("\t"+entry.getValue().floatValue()+"\t"+entry.getKey());
		}
		
		// reverse order is it's the COV
		if(anchorX != 0.0) {
			Collections.reverse(sortedNamesList);
			Collections.reverse(sortedVauesList);
		}

		
		
		ArrayList<DefaultXY_DataSet> funcList = new ArrayList<DefaultXY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		for(int i=0;i<sortedNamesList.size();i++) {
			DefaultXY_DataSet xyFunc = new DefaultXY_DataSet();
			double yVal = i+0.5;
			double xVal = sortedVauesList.get(i);
			xyFunc.set(xVal, yVal);
			xyFunc.set(anchorX, yVal);
			xyFunc.setName(sortedNamesList.get(i));
			xyFunc.setInfo("anchorX, yVal: "+(float)anchorX+", "+(float)xVal);
			funcList.add(xyFunc);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colorMap.get(sortedNamesList.get(i))));
		}
		
		String stringForInfo="";
		for(int i=sortedNamesList.size()-1;i>=0;i--) {
			stringForInfo += sortedNamesList.get(i)+"\n";
		}
		funcList.get(0).setInfo(funcList.get(0).getInfo()+"\n\nPLOT ORDER:\n"+stringForInfo);

		double xMinPlot = minX;
		double xMaxPlot = maxX;
		double yMinPlot = 0;
		double yMaxPlot = sortedNamesList.size();

		//add 10% change lines if anchor is 0.0 (fractional change plot)
		if(anchorX == 0.0) {
			DefaultXY_DataSet tenPercLowFunc = new DefaultXY_DataSet();
			tenPercLowFunc.set(-0.1, 0.0);
			tenPercLowFunc.set(-0.1, yMaxPlot);
			DefaultXY_DataSet tenPercHighFunc = new DefaultXY_DataSet();
			tenPercHighFunc.set(0.1, 0.0);
			tenPercHighFunc.set(0.1, yMaxPlot);
			funcList.add(tenPercLowFunc);
			funcList.add(tenPercHighFunc);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));			
		}
		else {	// add anchor line
			DefaultXY_DataSet anchorFunc = new DefaultXY_DataSet();
			anchorFunc.set(anchorX, 0.0);
			anchorFunc.set(anchorX, yMaxPlot);
			anchorFunc.setName("anchorFunc");
			funcList.add(anchorFunc);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));

		}

		PlotSpec specCombined = new PlotSpec(funcList, plotChars, null, xAxisName, null);
		specCombined.setLegendVisible(false);
		specCombined.setLegendLocation(RectangleEdge.RIGHT);
		
		if(savePlots) {
			String fname = fileName;
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
			gp.setTickLabelFontSize(22);
			gp.setAxisLabelFontSize(22);
			gp.setPlotLabelFontSize(22);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombined);
			try {
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(400, 500);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(specCombined);
//			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			gw.setAxisRange(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
//			gw.setBackgroundColor(Color.WHITE);	
		}
	}


	public void makeEAL_Historgram(boolean popUpWindows, boolean savePlots) {
		double min = 0;
		double max =7;
		int num = 140;
		
		HistogramFunction hist = new HistogramFunction(0.0+(max-min)/(num*2),max-(max-min)/(num*2), num);
		for(int i=0; i<branchNormEAL.length; i++) {
			hist.add(branchNormEAL[i], branchWt[i]);
		}
		
		hist.normalizeBySumOfY_Vals();
		hist.scale(1.0/hist.getDelta());
		
		String histName = "Normalized PDF; orig mean = "+totMeanEAL+" and orig stdDev = "+totStdDevEAL;
		hist.setName(histName);
		
		double meanFromPDF = hist.computeMean();
		infoString += "\nMean from Norm PDF = "+(float)meanFromPDF+"\n";
		double stdDevFromPDF = hist.computeStdDev();
		infoString += "StdDev from Norm PDF: " + (float)stdDevFromPDF+"\n";
		infoString += "COV from data Norm PDF = " + (float)(stdDevFromPDF/meanFromPDF)+"\n";
		infoString += "COV from orig (non-norm) data = " + (float)(totStdDevEAL/totMeanEAL)+"\n";
		
		HistogramFunction histCumulative = hist.getCumulativeDistFunctionWithHalfBinOffset();
		histCumulative.scale(hist.getDelta()); // so it sums to 1.0
		histCumulative.setName("Cumulative Dist");
		
		totFactor95perc = get95percConfFactorForValue(1.0, histCumulative);
		infoString += "\n95% Conf Bounds Factor About Mean: " + (float)totFactor95perc+"\n";
		totWithin10perc = histCumulative.getInterpolatedY(1.0*1.1) - histCumulative.getInterpolatedY(1.0*0.9);
		infoString += "Fract within 10% of Mean: " + (float)totWithin10perc+" (chance it's outside: "+(float)(1.0-totWithin10perc)+")\n";

		//totFactor95perc, totWithin10perc
		
		// Median info
		totMedianNormEAL = histCumulative.getFirstInterpolatedX(0.5);
		infoString += "\nMedian (normalized) from Cumulative Dist: " + (float)+totMedianNormEAL+"\n";
		int minDiffWtedIndex = getBranchClosestToValue(totMedianNormEAL, true);
		double factor95percent = get95percConfFactorForValue(totMedianNormEAL, histCumulative);
		infoString += "95% Conf Bounds Factor About Median: " + (float)factor95percent+"\n";
		double within10percent = histCumulative.getInterpolatedY(totMedianNormEAL*1.1) - histCumulative.getInterpolatedY(totMedianNormEAL*0.9);
		infoString += "Fract within 10% of Median: " + (float)within10percent+" (chance it's outside: "+(float)(1.0-within10percent)+")\n";
		infoString += "Relative Likelihood of Median vs Mean: "+(float)(hist.getInterpolatedY(totMedianNormEAL)/hist.getInterpolatedY(1.0));
		infoString += "\nBranch Closest to Median (weighted):\n";
		infoString += "\tindex = "+minDiffWtedIndex+"\n\twt = "+(float)branchWt[minDiffWtedIndex]+"\n\tEAL = "+(float)(branchNormEAL[minDiffWtedIndex]*totMeanEAL);
		for(String branchName:allBranchValuesMap.keySet())
			infoString += "\n\t"+branchName+" = "+allBranchValuesMap.get(branchName)[minDiffWtedIndex];

		
		int modeIndex = hist.getXindexForMaxY();
		totModalNormEAL = hist.getX(modeIndex);
		infoString += "\n\nMode (normalized) from PDF: " + (float)totModalNormEAL;
		minDiffWtedIndex = getBranchClosestToValue(totModalNormEAL, true);
		factor95percent = get95percConfFactorForValue(totModalNormEAL, histCumulative);
		infoString += "\n95% Conf Bounds Factor About Mode: " + (float)factor95percent+"\n";
		within10percent = histCumulative.getInterpolatedY(totModalNormEAL*1.1) - histCumulative.getInterpolatedY(totModalNormEAL*0.9);
		infoString += "Fract within 10% of Mode: " + (float)within10percent+" (chance it's outside: "+(float)(1.0-within10percent)+")\n";
		infoString += "Relative Likelihood of Mode vs Mean: "+(float)(hist.getInterpolatedY(totModalNormEAL)/hist.getInterpolatedY(1.0));
		infoString += "\nBranch Closest to Mode (weighted):\n";
		infoString += "\tindex = "+minDiffWtedIndex+"\n\twt = "+(float)branchWt[minDiffWtedIndex]+"\n\tEAL = "+(float)(branchNormEAL[minDiffWtedIndex]*totMeanEAL);
		for(String branchName:allBranchValuesMap.keySet())
			infoString += "\n\t"+branchName+" = "+allBranchValuesMap.get(branchName)[minDiffWtedIndex];
		
		
		hist.setInfo("Mean, Median, Mode from PDF: "+(float)meanFromPDF+", "+(float)totMedianNormEAL+", "+(float)totModalNormEAL);
		histCumulative.setInfo(hist.getInfo());
		
		HistogramFunction logNormFitHist = getLogNormalFitHist(hist);
		HistogramFunction logNormFitHistCumulative = logNormFitHist.getCumulativeDistFunctionWithHalfBinOffset();
		logNormFitHistCumulative.scale(logNormFitHistCumulative.getDelta());
		logNormFitHistCumulative.setInfo("Lognormal Dist Fit Cumulative");

		ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		funcsArray.add(hist);
		funcsArray.add(logNormFitHist);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcsArray, plotChars, null, "Normalized AAL", "Density");
//		specCombinedLegend.setLegendVisible(false);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);
		
		double yVal = hist.getMaxY()*1.1;
//		yVal = 2.25;

		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(spec);
			gw.setAxisRange(0d, 3d, 0d, yVal);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
		}

		if(savePlots) {
			String fname = "EAL_Histogram";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(0d, 3d, 0d, yVal);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec);
			try {
				File file = new File(ROOT_DIR+outDirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				gp.getChartPanel().setSize(500, 400);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				gp.saveAsTXT(file.getAbsolutePath() + ".txt");

			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		
		ArrayList<XY_DataSet> funcsArray2 = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		funcsArray2.add(histCumulative);
		funcsArray2.add(logNormFitHistCumulative);
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));

		PlotSpec specCum = new PlotSpec(funcsArray2, plotChars2, null, "Normalized AAL", "Cumulative Density");
//		specCombinedLegend.setLegendVisible(false);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);

		if(popUpWindows) {
			GraphWindow gw2 = new GraphWindow(specCum);
			gw2.setAxisRange(0d, 3d, 0d, 1d);
			gw2.setTickLabelFontSize(19);
			gw2.setAxisLabelFontSize(20);
			gw2.setPlotLabelFontSize(21);
		}

		if(outDirName != null) {
			String fname = "EAL_CumDist";
			HeadlessGraphPanel gp2 = new HeadlessGraphPanel();
			gp2.setUserBounds(0d, 3d, 0d, 1d);
			gp2.setTickLabelFontSize(23);
			gp2.setAxisLabelFontSize(24);
			gp2.setPlotLabelFontSize(26);
			gp2.setBackgroundColor(Color.WHITE);
			gp2.drawGraphPanel(specCum);
			try {
				File file = new File(ROOT_DIR+outDirName, fname);
//				gp2.getChartPanel().setSize(1000, 800);
//				gp2.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp2.saveAsPNG(file.getAbsolutePath() + ".png");
				gp2.getChartPanel().setSize(500, 400);
				gp2.saveAsPDF(file.getAbsolutePath()+".pdf");
				gp2.saveAsTXT(file.getAbsolutePath() + ".txt");

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	
	}
	
	/**
	 * This returns the 95% confidence bounds about the given value, defined as a multiplicative factor, meaning
	 * there is a 95% chance the true normalized EAL lies between value/factor and value*factor.
	 * @param normEAL_value
	 * @param cumDist
	 * @return
	 */
	private  static double get95percConfFactorForValue(double normEAL_value, HistogramFunction cumDist) {
		
		double maxX = cumDist.getMaxX()/normEAL_value;
		EvenlyDiscretizedFunc fractWithinScaleFactFunc = new EvenlyDiscretizedFunc(1.01,maxX,1000);
		for(int i=0;i<fractWithinScaleFactFunc.size();i++) {
			double scaleFact = fractWithinScaleFactFunc.getX(i);
			double fractWithin = cumDist.getInterpolatedY(scaleFact*normEAL_value) - cumDist.getInterpolatedY(normEAL_value/scaleFact);
			fractWithinScaleFactFunc.set(i,fractWithin);
		}
		return fractWithinScaleFactFunc.getFirstInterpolatedX(0.95);
	}
	
	
	private HistogramFunction getLogNormalFitHist(HistogramFunction hist) {
		
		LognormalDistCalc logNormCalc = new LognormalDistCalc();
		logNormCalc.fitToThisFunction(hist, 0.5, 1.0, 100, 0.1, 1.0, 100);
		EvenlyDiscretizedFunc fitHist = logNormCalc.getPDF();
		HistogramFunction logNormHist = new HistogramFunction(fitHist.getMinX(),fitHist.getMaxX(),fitHist.size());
		for(int i=0;i<logNormHist.size();i++)
			logNormHist.set(i,fitHist.getY(i));
		logNormHist.normalizeBySumOfY_Vals();
		logNormHist.scale(1.0/logNormHist.getDelta());
		logNormHist.setName("Lognormal Dist Fit");
		return logNormHist;
	}
	
	
	/**
	 * This computes the factor 95% and fract 10% as a function of COV for a perfect lognormal distribution
	 */
	private static void computeFact95_fract10_vsLogNormCOV() {
		
		for( double cov=0.45;cov >=0.009; cov -= 0.01) {
			double[] result = computeFact95_fract10_fromLogNorm_COV(cov);
			System.out.println((float)cov+"\t"+(float)result[0]+"\t"+(float)result[1]);		
		}
	}
	
	/**
	 * This computes the factor 95% and fract 10% as a function of COV for a perfect lognormal distribution
	 */
	/**
	 * This computes the factor about the mean which contains 95% of the distribution (fact 95) and the fraction
	 * of the distribution that is within 10% of the mean for a Lognormal distributions with the given COV.
	 * @param cov
	 * @return double array where element 0 is factor 95% and 1 is fract 10%
	 */
	private static double[] computeFact95_fract10_fromLogNorm_COV(double cov) {

		double[] result = new double[2];
		LognormalDistCalc logNormCalc = new LognormalDistCalc();
		int numPoints = 10000;
		double deltaX = 5.0/numPoints;
		logNormCalc.setAll(1.0, cov, deltaX, numPoints);
		EvenlyDiscretizedFunc pdf = logNormCalc.getPDF();
		HistogramFunction pdfHist = new HistogramFunction(pdf.getMinX(),pdf.getMaxX(),pdf.size());
		for(int i=0;i<pdfHist.size();i++)
			pdfHist.set(i,pdf.getY(i));
		pdfHist.normalizeBySumOfY_Vals();
		HistogramFunction cumDist = pdfHist.getCumulativeDistFunctionWithHalfBinOffset();
		result[0] = get95percConfFactorForValue(1.0, cumDist);
		result[1] = cumDist.getInterpolatedY(1.1) - cumDist.getInterpolatedY(0.9);
		return result;
	}
	
	
	public void makeLog10_EAL_Historgram() {
		double min = Math.log10(1);
		double max = Math.log10(20);
		int num = 40;
		double delta = (max-min)/num;
		
		HistogramFunction hist = new HistogramFunction(min+delta/2,max-delta/2, num);
		for(int i=0; i<branchEAL.length; i++) {
			hist.add(Math.log10(branchEAL[i]), branchWt[i]);
		}
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plottingFuncsArray.add(hist);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		
		GraphWindow graph = new GraphWindow(plottingFuncsArray, "Log10 EAL Distribution");
//		graph.setX_AxisRange(xAxisRange.getLowerBound(),xAxisRange.getUpperBound());
//		graph.setY_AxisRange(yAxisRange.getLowerBound(),yAxisRange.getUpperBound());
//		graph.setXLog(logX);
//		graph.setYLog(logY);
		graph.setPlotChars(plotChars);
		graph.setX_AxisLabel("log10 EAL (Billion $)");
		graph.setY_AxisLabel("Density");
		graph.setTickLabelFontSize(18);
		graph.setAxisLabelFontSize(20);

	}
	
	/**
	 * This returns a list of key Strings sorted by the double values in the map
	 */
	public List<String> getKeysSortedByMapValues(HashMap<String, Double> valuesMap) {
		//convert map to a List
		List<Entry<String, Double>> list = new LinkedList<Map.Entry<String, Double>>(valuesMap.entrySet());

		//sorting the list with a comparator
		Collections.sort(list, new Comparator<Entry<String, Double>>() {
			public int compare(Map.Entry<String, Double> o1, Map.Entry<String, Double> o2) {
				return (o1.getValue()).compareTo(o2.getValue());
			}
		});

		//convert sorted list back to Map
		ArrayList<String> sortedKeys = new ArrayList<String>();
		for (Entry<String, Double> entry : list) {
			sortedKeys.add(entry.getKey());
		}
		return sortedKeys;
	}
	
	
	public HashMap<String, Color> getColorMap(HashMap<String, Double> valuesMap) {
//		Map<String, Integer> unsortedMap = new HashMap<String, Integer>();

		//convert map to a List
		List<Entry<String, Double>> list = new LinkedList<Map.Entry<String, Double>>(valuesMap.entrySet());

		//sorting the list with a comparator
		Collections.sort(list, new Comparator<Entry<String, Double>>() {
			public int compare(Map.Entry<String, Double> o1, Map.Entry<String, Double> o2) {
				return (o1.getValue()).compareTo(o2.getValue());
			}
		});

		//convert sorted list back to Map
		Map<String, Double> sortedMap = new LinkedHashMap<String, Double>();
		for (Entry<String, Double> entry : list) {
			sortedMap.put(entry.getKey(), entry.getValue());
		}
		
		HashMap<String, Color> colorMap = new HashMap<String, Color>();
		  
        // switch statement with int data type 
		Object[] keySet = sortedMap.keySet().toArray();
        switch (keySet.length) { 
        case 1: 
        		colorMap.put((String)keySet[0], Color.BLUE); 
            break; 
        case 2: 
			colorMap.put((String)keySet[0], Color.BLUE); 
			colorMap.put((String)keySet[1], Color.RED); 
            break; 
        case 3: 
			colorMap.put((String)keySet[0], Color.BLUE); 
			colorMap.put((String)keySet[1], Color.GREEN); 
			colorMap.put((String)keySet[2], Color.RED); 
            break; 
        case 4: 
			colorMap.put((String)keySet[0], Color.BLUE); 
			colorMap.put((String)keySet[1], Color.GREEN); 
			colorMap.put((String)keySet[2], Color.ORANGE); 
			colorMap.put((String)keySet[3], Color.RED); 
            break; 
        case 5: 
			colorMap.put((String)keySet[0], Color.BLUE); 
			colorMap.put((String)keySet[1], Color.GREEN); 
			colorMap.put((String)keySet[2], Color.ORANGE); 
			colorMap.put((String)keySet[3], Color.RED); 
			colorMap.put((String)keySet[4], Color.MAGENTA); 
            break; 
        default: 
            throw new RuntimeException("Too many elements"); 
        } 
		
		return colorMap;

	}
		
	
	public double getTotMeanEAL() {
		return totMeanEAL;
	}
	
	public double getTotMedianNormEAL() {
		return totMedianNormEAL;
	}
	
	public double getTotModalNormEAL() {
		return totModalNormEAL;
	}
	
	public double getTotStdDevEAL() {
		return totStdDevEAL;
	}
	
	public double getTotCOV_EAL() {
		return totCOV_EAL;
	}
	
	public double getTotFactor95perc() {
		return totFactor95perc;
	}
	
	public double getTotWithin10perc() {
		return totWithin10perc;
	}
	
	
	/**
	 * This looks at reduction in COV and 95% bounds, and increase in fraction within 10% of mean, as branches are elliminated, 
	 * where the next branches to removed is based on the previous result (and branch kept is either the highest weighted or the
	 * one that maximizes the next COV reduction).
	 * 
	 */
	public static void doAllAnalysis1() {
		
		DefaultXY_DataSet covValues = new DefaultXY_DataSet();
		DefaultXY_DataSet fact95Values = new DefaultXY_DataSet();
		DefaultXY_DataSet fractWithin10percValues = new DefaultXY_DataSet();
		
		ArrayList<String> branchesToRemove = new ArrayList<String>();
		
//		branchesToRemove.add("Ground Motion Model"+" = "+"IDRISS_2014");
		
		Analysis analysis0 = new Analysis("all_branch_results.csv", "analysis1_fullReg_allBranches", branchesToRemove, false, false);
		covValues.set(0.0,analysis0.getTotCOV_EAL());
		fact95Values.set(0.0,analysis0.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis0.getTotWithin10perc());
		
		// keep top weighted branch
		branchesToRemove.add("GMM Added Uncertainty"+" = "+"LOWER");
		branchesToRemove.add("GMM Added Uncertainty"+" = "+"UPPER");
		
		Analysis analysis1 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_01", branchesToRemove, false, false);
		covValues.set(1.0,analysis1.getTotCOV_EAL());
		fact95Values.set(0.0,analysis1.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis1.getTotWithin10perc());

		// keep top weighted branch
		branchesToRemove.add("Total Mag>5 Rate"+" = "+"6.5");
		branchesToRemove.add("Total Mag>5 Rate"+" = "+"9.6");
		
		Analysis analysis2 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_02", branchesToRemove, false, false);
		covValues.set(2.0,analysis2.getTotCOV_EAL());
		fact95Values.set(0.0,analysis2.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis2.getTotWithin10perc());
		
		// the branch kept has lowest resultant COV
		branchesToRemove.add("Scaling Relationship"+" = "+"EllBsqrtLen");
		branchesToRemove.add("Scaling Relationship"+" = "+"HB08");
		branchesToRemove.add("Scaling Relationship"+" = "+"ShConStrDrp");
		branchesToRemove.add("Scaling Relationship"+" = "+"Shaw09Mod");
		
		Analysis analysis3 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_03", branchesToRemove, false, false);
		covValues.set(3.0,analysis3.getTotCOV_EAL());
		fact95Values.set(0.0,analysis3.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis3.getTotWithin10perc());
		
		
		// the branch kept has lowest resultant COV
		branchesToRemove.add("Ground Motion Model"+" = "+"ASK_2014");
		branchesToRemove.add("Ground Motion Model"+" = "+"BSSA_2014");
		branchesToRemove.add("Ground Motion Model"+" = "+"CY_2014");
		branchesToRemove.add("Ground Motion Model"+" = "+"IDRISS_2014");
		
		Analysis analysis4 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_04", branchesToRemove, false, false);
		covValues.set(4.0,analysis4.getTotCOV_EAL());
		fact95Values.set(0.0,analysis4.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis4.getTotWithin10perc());

		
		// the branch kept has lowest resultant COV
		branchesToRemove.add("Spatial Seismicity PDF"+" = "+"U2");
		
		Analysis analysis5 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_05", branchesToRemove, false, false);
		covValues.set(5.0,analysis5.getTotCOV_EAL());
		fact95Values.set(0.0,analysis5.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis5.getTotWithin10perc());

		// keep top weighted branch
		branchesToRemove.add("Probability Model"+" = "+"POISSON");
		branchesToRemove.add("Probability Model"+" = "+"LOW_VALUES");
		branchesToRemove.add("Probability Model"+" = "+"HIGH_VALUES");
		
		Analysis analysis6 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_06", branchesToRemove, false, false);
		covValues.set(6.0,analysis6.getTotCOV_EAL());
		fact95Values.set(0.0,analysis6.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis6.getTotWithin10perc());

		// keep the branch with lowest COV among the top weighted branches
		branchesToRemove.add("Deformation Model"+" = "+"ABM");
		branchesToRemove.add("Deformation Model"+" = "+"NEOK");
		branchesToRemove.add("Deformation Model"+" = "+"GEOL");
		
		Analysis analysis7 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_07", branchesToRemove, false, false);
		covValues.set(7.0,analysis7.getTotCOV_EAL());
		fact95Values.set(0.0,analysis7.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis7.getTotWithin10perc());

		
		// the branch kept has lowest resultant COV
		branchesToRemove.add("Vs30 Model"+" = "+"WillsEtAl");
		
		Analysis analysis8 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_08", branchesToRemove, false, false);
		covValues.set(8.0,analysis8.getTotCOV_EAL());
		fact95Values.set(0.0,analysis8.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis8.getTotWithin10perc());

		// keep top weighted branch
		branchesToRemove.add("Mmax Off Fault"+" = "+"7.9");
		branchesToRemove.add("Mmax Off Fault"+" = "+"7.3");
		
		Analysis analysis9 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_09", branchesToRemove, false, false);
		covValues.set(9.0,analysis9.getTotCOV_EAL());
		fact95Values.set(0.0,analysis9.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis9.getTotWithin10perc());

		// the branch kept has lowest resultant COV
		branchesToRemove.add("Slip Along Rupture (Dsr)"+" = "+"Uniform");
		
		Analysis analysis10 = new Analysis("all_branch_results.csv", "analysis1_fullReg_BrRemoval_10", branchesToRemove, false, false);
		covValues.set(10.0,analysis10.getTotCOV_EAL());
		fact95Values.set(0.0,analysis10.getTotFactor95perc());
		fractWithin10percValues.set(0.0,analysis10.getTotWithin10perc());

		
		
		// last would be to remove either Fault Model to zero COV
		covValues.set(11.0,0.0);
		fact95Values.set(0.0,0.0);
		fractWithin10percValues.set(0.0,1.0);

		covValues.setName("covValues");
		fact95Values.setName("fact95Values");
		fractWithin10percValues.setName("fractWithin10percValues");


		System.out.println(covValues);
		System.out.println(fact95Values);
		System.out.println(fractWithin10percValues);
		
		
		ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		funcsArray.add(covValues);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		
		PlotSpec spec = new PlotSpec(funcsArray, plotChars, "", "Branch Removed", "COV");
//		specCombinedLegend.setLegendVisible(false);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);
		
		double yMax = 0.5;
		double xMax = covValues.getMaxX()-0.5;

		if(true) {
			GraphWindow gw = new GraphWindow(spec);
			gw.setAxisRange(-0.5, xMax, 0d, yMax);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
		}

			String fname = "analysis1_COV_AsBramchesRemoved";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(-0.5, xMax, 0d, yMax);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec);
			try {
				File file = new File(ROOT_DIR, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				gp.getChartPanel().setSize(500, 400);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}



	}
	
	
	

	/**
	 * This looks at reduction in COV and 95% bounds, and increase in fraction within 10% of mean, as branches are eliminated, 
	 * where the next branch to removed is based on wt-average COV reduction implied by the run with all branches (hard coded here).
	 * The option that is kept for each branch is that with a COV closest to the weight average.
	 * 
	 */
	public static void doCOV_ReductionAnalysis() {
		
		DefaultXY_DataSet covValues = new DefaultXY_DataSet();
		DefaultXY_DataSet fact95Values = new DefaultXY_DataSet();
		DefaultXY_DataSet fractWithin10percValues = new DefaultXY_DataSet();
		
		String[] branchesToRemoveNamesArray = {
				"GMM Added Uncertainty",
				"Total Mag>5 Rate",
				"Probability Model",
				"Scaling Relationship",
				"Ground Motion Model",
				"Spatial Seismicity PDF",
				"Deformation Model",
				"Mmax Off Fault",
				"Vs30 Model",
				"Slip Along Rupture (Dsr)",
				"Fault Model"
		};
		
//		String[] reversedArray = new String[branchesToRemoveNamesArray.length];
//		int j=0;
//		for(int i = branchesToRemoveNamesArray.length-1;i>=0;i--) {
//			reversedArray[j] = branchesToRemoveNamesArray[i];
//			j+=1;
//		}
//		branchesToRemoveNamesArray = reversedArray;
//		for(String brName:branchesToRemoveNamesArray)
//			System.out.println(brName);
		
		
		ArrayList<String> branchesToRemoveNameValue = new ArrayList<String>();
		
		for(int index=0;index<branchesToRemoveNamesArray.length;index++) {
			String branchName = branchesToRemoveNamesArray[index];
			String fileSuffix = "";
			if(index<10)
				fileSuffix = "0"+index;
			else
				fileSuffix = ""+index;
			String name = "COV_ReductionAnalysisBrRemoval_"+fileSuffix;
			Analysis analysis = new Analysis("all_branch_results.csv", name, branchesToRemoveNameValue, false, false);
			covValues.set((double)index, analysis.getTotCOV_EAL());
			fact95Values.set((double)index, analysis.getTotFactor95perc());
			fractWithin10percValues.set((double)index, analysis.getTotWithin10perc());
			branchesToRemoveNameValue.addAll(analysis.getBranchesToRemoveForCOV_List(branchName));
		}
		
		
		// last would be to remove either Fault Model to zero COV
		int index = covValues.size();
		System.out.println("index="+index);
		covValues.set((double)index,0.0);
		fact95Values.set((double)index,0.0);
		fractWithin10percValues.set((double)index,1.0);

		covValues.setName("covValues");
		fact95Values.setName("fact95Values");
		fractWithin10percValues.setName("fractWithin10percValues");


		System.out.println(covValues);
		System.out.println(fact95Values);
		System.out.println(fractWithin10percValues);
		
		
		ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		funcsArray.add(covValues);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		
		PlotSpec spec = new PlotSpec(funcsArray, plotChars, "", "Branch Removed", "COV");
//		specCombinedLegend.setLegendVisible(false);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);
		
		double yMax = 0.5;
		double xMax = covValues.getMaxX()+0.5;

		if(true) {
			GraphWindow gw = new GraphWindow(spec);
			gw.setAxisRange(-0.5, xMax, 0d, yMax);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
		}

			String fname = "COV_ReductionAnalysis";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(-0.5, xMax, 0d, yMax);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec);
			try {
				File file = new File(ROOT_DIR, fname);
				gp.getChartPanel().setSize(500, 400);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}



	}
	
	
	public void writeFigure1_Data() {
		String outputString = "";
		ArrayList<String> nameList = getLogicTreeBranchValueOrder();
		for(String name:nameList) {
			String elements[] = name.split(" =");
			double wt = weightForBrValMap.get(name);
			outputString += elements[0]+"\t"+elements[1]+"\t"+(float)wt+"\n";
		}
		System.out.println(outputString);
	}

	
	public String getNameChange(String name) {
		if(nameChangeMap == null) {			
			nameChangeMap = new HashMap<String, String>();
			nameChangeMap.put("Slip Along Rupture Model (Dsr)", "Slip Along Rupture (Dsr)");
			nameChangeMap.put("MMax Off Fault", "Mmax Off Fault");
			nameChangeMap.put("GMM Additional Epistemic Uncertainty", "GMM Added Uncertainty");
			nameChangeMap.put("ERF Probability Model", "Probability Model");
			nameChangeMap.put("Total Mag 5 Rate", "Total Mag>5 Rate");
			nameChangeMap.put("ZENGBB", "ZENG");
			nameChangeMap.put("Tap", "Tapered");
			nameChangeMap.put("Uni", "Uniform");
			nameChangeMap.put("FM3_1", "FM3.1");
			nameChangeMap.put("FM3_2", "FM3.2");
			nameChangeMap.put("Wills2015", "WillsEtAl");
		}
		
		if(nameChangeMap.keySet().contains(name))
			return nameChangeMap.get(name);
		else
			return name;
	}
	
	
	public static void main(String[] args) {
	
//		// for Figure 1 data:
//		Analysis analysisFig1 = new Analysis("all_branch_results.csv", "junk", null, false);
//		analysisFig1.writeFigure1_Data();
		
		// CEA PORTFOLIO ANALYSIS ******************
		// for Figure 2-7 data for CEA-specific portolio from: 
		// http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2020_09_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC-covModel/all_branch_results.csv
		// OLD: http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2020_04_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC/all_branch_results.csv
//		Analysis cea_analysisFig2to7 = new Analysis("ceaPortfolio_all_branch_results.csv", "CEA_PortfolioFigure2to7_Data", null, false, true);
	
		//dataCol: - 2 for AAL, 18 for Loss@0.01; 19 for Loss@0.004; 20 for Loss@0.0025; 21 for Loss@0.0018; and 22 for Loss@4.0E-4
//		Analysis cea_analysisFig2to7_ = new Analysis("ceaPortfolio_all_branch_results.csv", "CEA_PortfolioFigure2to7_Data_LossAt0.01", null, false, true, 18);
//		Analysis cea_analysisFig2to7_ = new Analysis("ceaPortfolio_all_branch_results.csv", "CEA_PortfolioFigure2to7_Data_LossAt0.004", null, false, true, 19);
//		Analysis cea_analysisFig2to7_ = new Analysis("ceaPortfolio_all_branch_results.csv", "CEA_PortfolioFigure2to7_Data_LossAt0.0025", null, false, true, 20);
//		Analysis cea_analysisFig2to7_ = new Analysis("ceaPortfolio_all_branch_results.csv", "CEA_PortfolioFigure2to7_Data_LossAt0.0018", null, false, true, 21);
		Analysis cea_analysisFig2to7_ = new Analysis("ceaPortfolio_all_branch_results.csv", "CEA_PortfolioFigure2to7_Data_LossAt4e-4", null, false, true, 22);

		
		
		
		
//		// for Figure 2-7 data:
//		Analysis analysisFig2to7 = new Analysis("all_branch_results.csv", "Figure2to7_Data", null, false, true);
		
		
//		computeFact95_fract10_vsLogNormCOV();
		
//		doCOV_ReductionAnalysis();
		
//		ArrayList<String> branchesToRemove = new ArrayList<String>();
//		Analysis analysisFig8 = new Analysis("all_branch_results.csv", "Figure8_Data", branchesToRemove, false, false);
//		// This takes a very long time:
////		analysisFig8.doCOV_ReductionAnalysisComplete(true, true);
//		// this recreates the plot from saved data
//		analysisFig8.makeCOV_ReductionAnalysisCompletePlotsAgain();
		
		// LA Analysis
		// from:  http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2020_03_17-ucerf3-ngaw2-cea-consolidate-los-angeles/
//		Analysis la_analysis = new Analysis("LosAngeles_all_branch_results.csv", "LosAngeles_Data", null, false, true);

//		// SF Analysis
//		// from:  http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2020_03_17-ucerf3-ngaw2-cea-consolidate-san-francisco/
//		Analysis sf_analysis = new Analysis("SanFrancisco_all_branch_results.csv", "SanFrancisco_Data", null, false, true);

		
//		// Fairfield census tracts
//		// from: http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2019_02_13-ucerf3-ngaw2-cea-consolidate-fairfield/
//		// for Figure 2-4 data:
//		Analysis fairfieldAnalysis = new Analysis("fairfield_all_branch_results.csv", "FairfieldResults", null, false, true);

// OLD:		
//		Analysis analysis0 = new Analysis("all_branch_results.csv", "fullReg_allBranches", branchesToRemove, false);
//
//		branchesToRemove.add("GMM Added Uncertainty"+" = "+"LOWER");
//		branchesToRemove.add("GMM Added Uncertainty"+" = "+"UPPER");
//		
//		Analysis analysis1 = new Analysis("all_branch_results.csv", "fullReg_BrRemoval_1", branchesToRemove, false);
		
	}

}
