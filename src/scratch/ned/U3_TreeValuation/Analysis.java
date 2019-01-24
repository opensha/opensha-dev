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
import org.opensha.sha.earthquake.calc.recurInterval.LognormalDistCalc;

import com.google.common.io.Files;

import scratch.UCERF3.inversion.CommandLineInversionRunner;

public class Analysis {
	
	/* Data files came from:
	*  http://opensha.usc.edu/ftp/kmilner/ucerf3/eal_calcs/2018_11_20-ucerf3-ngaw2-cea-consolidate/
	*/
	
	final static String ROOT_DIR = "src/scratch//ned/U3_TreeValuation/";
	
	String inputFileName, outDirName;

	String infoString = "";
	double[] branchWt, branchEAL, branchNormEAL;
	HashMap<String, String[]> allBranchValuesMap;
	HashMap<String, ArrayList<String>> optionsForBranchHashMap;
	
	int totNumBranches;
	double totMeanEAL, totMedianNormEAL, totModalNormEAL;
	double totStdDevEAL, totCOV_EAL;
	double totOrigWeight;
	HashMap<String, String> closestValuesForBranchMap; // this stores the branch value that has a mean closest to the totalMeanEAL
	
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
	HashMap<String, Double> covIfBrRemovedMap;
	HashMap<String, Double> covDiffIfBrRemovedMap;
	HashMap<String, Double> fracWithIn10percIfBrRemovedMap;
	HashMap<String, Double> fracWithIn10percDiffIfBrRemovedMap;
	HashMap<String, Double> factor95percIfBrRemovedMap;
	HashMap<String, Double> mfactor95percRatioDiffIfBrRemovedMap;


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

	
	public Analysis(String fileInputName, String outDirName) {
		
		boolean popUpPlots = true;
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
		
//		readBranchLevelSummaryDataFromFile();
		readAllBranchDataFromFile();
		
		generateBranchValueResults(false, true);
//		
//		makeFractionChangePlot(true, true, meanDiffForBrValMap, "meanDiffForBrValMap");
		
//		makeEAL_Historgram(popUpPlots, savePlots);
		
		System.out.println(infoString);
		
		//		makeLog10_EAL_Historgram();

		
		// hard coded assignment of branch type:
		erf_branches = new ArrayList<String>();
		gmm_branches = new ArrayList<String>();

		erf_branches.add("Fault Model");
		erf_branches.add("Deformation Model");
		erf_branches.add("Scaling Relationship");
		erf_branches.add("Slip Along Rupture Model (Dsr)");
//		erf_branches.add("Inversion Model");
		erf_branches.add("Total Mag 5 Rate");
		erf_branches.add("MMax Off Fault");
//		erf_branches.add("Moment Rate Fixes");
		erf_branches.add("Spatial Seismicity PDF");
		erf_branches.add("ERF Probability Model");

		gmm_branches.add("Ground Motion Model");
		gmm_branches.add("GMM Additional Epistemic Uncertainty");
		gmm_branches.add("Vs30 Model");

//		make_ERF_vs_GMM_UncertHists(true, true);
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
		HashMap<String, ArrayList<String>> optionsForBranchHashMapNew;
		
		ArrayList<String> optionsList = optionsForBranchHashMap.get(branchName);
		
		int oldNumOptions = optionsList.size();
		int totNumBranchesNew = totNumBranches*(oldNumOptions-1)/oldNumOptions;
		double[] branchWtNew = new double[totNumBranchesNew];
		double[] branchNormEAL_New = new double[totNumBranchesNew];
		
		optionsList.remove(branchValue);
		int newNumOptions = optionsForBranchHashMap.get(branchName).size();
		System.out.println("oldNumOptions="+oldNumOptions+"; newNumOptions="+newNumOptions);
		
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
		
		double newMean = this.computeWeightedAverage(branchWtNew, branchNormEAL_New);
		for(int i=0; i<branchNormEAL_New.length;i++) branchNormEAL_New[i] /= newMean;

		totMeanEAL = newMean;
		totStdDevEAL = computeWeightedStdDev(branchWtNew, branchNormEAL_New, totMeanEAL);
		totCOV_EAL = totStdDevEAL/totMeanEAL;


	}
	
	
	/**
	 * This reads the branch level data
	 * @return
	 */
	private void readAllBranchDataFromFile() {
		
		File file = new File(ROOT_DIR+inputFileName);
		List<String> fileLines;
		double[] branchEAL;
		
		infoString += "\nBranch Names:\n\n";
		
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			int numData = fileLines.size()-1;
			numData = 4*numData/5;	// removing "IDRISS_2014" GMM data
			branchWt = new double[numData];
			branchEAL = new double[numData];
			branchNormEAL = new double[numData];
			allBranchValuesMap = new HashMap<String, String[]>();
			String str = fileLines.get(0);
			String[] colNamesArray = str.split(",");
			for(int c=5;c<18; c++) {
				String keyName = colNamesArray[c];
				String[] strArray = new String[numData];
				allBranchValuesMap.put(keyName, strArray);
				infoString += "\t"+keyName+"\n";
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
				// filter out "IDRISS_2014" GMM
				if(split[15].equals("IDRISS_2014")) {
					continue;
				}
//for(int j=0;j<split.length;j++)
//	System.out.println(j+"\t"+split[j]);
//System.exit(-1);
				branchWt[arrayIndex] = Double.parseDouble(split[1]);
				branchEAL[arrayIndex] = Double.parseDouble(split[2])*1e-6;	// convert to Billions
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
					valArray[arrayIndex] = split[c];
				}
				arrayIndex +=1;
			}
			
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
			for(int i=0;i<branchNormEAL.length;i++) {
				branchNormEAL[i] = branchEAL[i]/totMeanEAL;
				branchWt[i] /= totOrigWeight;
			}
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
			
			optionsForBranchHashMap = new HashMap<String, ArrayList<String>>();
			for(String brName:allBranchValuesMap.keySet()) {
				String[] allValues = allBranchValuesMap.get(brName);
				ArrayList<String> optionsList = new ArrayList<String>();
				for(String opt:allValues)
					if(!optionsList.contains(opt))
						optionsList.add(opt);
				optionsForBranchHashMap.put(brName,optionsList);
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
	

	public void generateBranchValueResults(boolean popUpWindows, boolean saveResults) {

		double min = 0;
		double max = 6;
		int num = 120;

		ArrayList<String> tableLinesList = new ArrayList<String>();
		String tableHeaderLine = "branchName\tbranchValue\tweight\tmean\tmeanDiff\twtedMeanDiff\tcov\tfactor95perc\twithin10perc\tmeanIfRemoved\t"+
				"meanIfRemovedDiff\tcovIfRemoved\tfactor95percIfRemoved\twithin10percIfRemoved";	
		tableLinesList.add(tableHeaderLine);
		
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
		covIfBrRemovedMap = new HashMap<String, Double>();
		covDiffIfBrRemovedMap = new HashMap<String, Double>();
		fracWithIn10percIfBrRemovedMap = new HashMap<String, Double>();
		fracWithIn10percDiffIfBrRemovedMap = new HashMap<String, Double>();
		factor95percIfBrRemovedMap = new HashMap<String, Double>();
		mfactor95percRatioDiffIfBrRemovedMap = new HashMap<String, Double>();

		infoString += "\n\nExpected COV etc if Branches Removed:\n\n" ;
		infoString += "name\tcov\tcovRatio\tfact95perc\twithin10perc\n";

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
						(float)meanForBr + "\t" +
						(float)(meanForBr-1.0) + "\t" +
						(float)((meanForBr-1.0)*wtForBr) + "\t"+
						(float)covForBr + "\t" +
						(float)factor95percForBr + "\t" +
						(float)within10percForBr + "\t" +
						
						(float)meanForOtherBrs + "\t" +
						(float)(meanForOtherBrs-1.0) + "\t" +
						(float)covForOtherBrs + "\t" +
						(float)factor95percForOtherBrs + "\t" +
						(float)within10percForOtherBrs + "\t";
				tableLinesList.add(tableLine);
				
				String combinedName = branchName+" = "+brOpt;
				meanForBrValMap.put(combinedName, meanForBr);
				meanDiffForBrValMap.put(combinedName, meanForBr-1.0);

			}
			
			infoString += branchName+"\t"+(float)expCOV_ifBranchesRemoved+"\t"+(float)(expCOV_ifBranchesRemoved/totCOV_EAL)+
					"\t"+(float)expFactor95percIfBranchesRemoved+"\t"+(float)expwithin10percIfBranchesRemoved+"\n";
			
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

				if(saveResults) {
					try{
						FileWriter fw = new FileWriter(ROOT_DIR+outDirName+"/branchValueStats.txt");
						for(String line:	 tableLinesList)		
							fw.write(line+"\n");
						fw.close();
					}catch(Exception e) {
						e.printStackTrace();
					}			
				}

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
			plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colorMap.get(branchOption)));
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
		totMeanLine.set(totMeanEAL_here,1.0);
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

		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		plotChars2.addAll(plotCharsLines);
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars2.addAll(plotCharsLines);

		PlotSpec specCombined2 = new PlotSpec(funcs2, plotChars2, branchName, "Normalized Loss", "Density");
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
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(26);
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
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
	}


	
	
	public void old_makeStackedHistograms(boolean popUpWindows, String dirName) {
		
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdirs();
				
		double min = 0;
		double max = 6;
		int num = 120;
		
		double xMinPlot=0;
		double xMaxPlot=3;
		
		closestValuesForBranchMap = new HashMap<String, String> ();
		
		for(String branchName:allBranchValuesMap.keySet()) {
//System.out.println("working on" + branchName);
			String[] branchValuesArray = allBranchValuesMap.get(branchName);
			HashMap<String, HistogramFunction> histogramsHashMap = new HashMap<String, HistogramFunction>();
			HashMap<String, Double> meanHashMap = new HashMap<String, Double>();		
			HashMap<String, Double> weightHashMap = new HashMap<String, Double>();	
			for(int i=0; i<branchNormEAL.length; i++) {
				String branchValue = branchValuesArray[i];
				HistogramFunction hist;
				if(!histogramsHashMap.keySet().contains(branchValue)) {
					hist = new HistogramFunction(min+(max-min)/(num*2),max-(max-min)/(num*2), num);
					hist.setName(branchValue);
					histogramsHashMap.put(branchValue, hist);
					meanHashMap.put(branchValue, branchNormEAL[i]*branchWt[i]);
					weightHashMap.put(branchValue, branchWt[i]);
				}
				else {
					hist = histogramsHashMap.get(branchValue);
					double newMean = meanHashMap.get(branchValue)+branchNormEAL[i]*branchWt[i];
					meanHashMap.put(branchValue, newMean);
					double newWt = weightHashMap.get(branchValue)+branchWt[i];
					weightHashMap.put(branchValue,newWt);
				}
				hist.add(branchNormEAL[i], branchWt[i]);
			}
			if(histogramsHashMap.keySet().size() == 1) {
				continue;
			}
			
			// finish weights and means and correct histograms
			double maxHistValue=0;
			double totWt =0;
			double minMeanDiff = Double.MAX_VALUE;
			String closestValue=null;
			for(Double wt:weightHashMap.values())
				totWt += wt;
			double meanTest = 0;
			for(String branchValue: histogramsHashMap.keySet()) {
				double mean = meanHashMap.get(branchValue)/weightHashMap.get(branchValue);
				meanHashMap.put(branchValue,mean);
				double wt = weightHashMap.get(branchValue)/totWt;
				meanTest += wt*mean;
				weightHashMap.put(branchValue, wt);
				HistogramFunction hist = histogramsHashMap.get(branchValue);
				
				hist.normalizeBySumOfY_Vals();
				hist.scale(wt/hist.getDelta());

//				double test = hist.calcSumOfY_Vals()*hist.getDelta();
//System.out.println("test: "+(float)(wt/test));
				if(maxHistValue <hist.getMaxY())
					maxHistValue = hist.getMaxY();
				double meanDiff = Math.abs(1.0-mean);
				if(minMeanDiff>meanDiff) {
					minMeanDiff=meanDiff;
					closestValue = branchValue;
				}
			}
			
			if(Math.abs(meanTest-1.0)>0.00001)
				throw new RuntimeException("something wrong with mean value");
			
			closestValuesForBranchMap.put(branchName, closestValue);
			
			
			ArrayList<HistogramFunction> histArrayList = new ArrayList<HistogramFunction>();
//			for(HistogramFunction hist:histogramsHashMap.values())
//				histArrayList.add(hist);
			for(String branchValue:meanHashMap.keySet())
				histArrayList.add(histogramsHashMap.get(branchValue));
//System.out.println(histArrayList.size()+"; "+ histogramsHashMap.keySet().size());
			List<HistogramFunction> stackedHistograms = HistogramFunction.getStackedHists(histArrayList, false);
			
			
			
			HashMap<String, Color> colorMap = getColorMap(meanHashMap);
			ArrayList<PlotCurveCharacterstics> plotCharsHists = new ArrayList<PlotCurveCharacterstics>();
			ArrayList<PlotCurveCharacterstics> plotCharsLines = new ArrayList<PlotCurveCharacterstics>();
			ArrayList<PlotCurveCharacterstics> plotCharsValuePDFs = new ArrayList<PlotCurveCharacterstics>();

			for(String branchOption: meanHashMap.keySet()) {
//System.out.println(branchName+"\t"+branchOption+"\t"+colorMap.get(branchOption)+"\t"+meanHashMap.get(branchOption));
				plotCharsHists.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colorMap.get(branchOption)));
				plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colorMap.get(branchOption)));
				plotCharsValuePDFs.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colorMap.get(branchOption)));
			}

			
			ArrayList<DefaultXY_DataSet> meanLinesList = new ArrayList<DefaultXY_DataSet>();
			double totMeanEAL_here=0;
			// the following are to plot the mean lines offset from the graph
			double maxYhist = stackedHistograms.get(0).getMaxY();
//System.out.println("maxYhist before = "+maxYhist);
			maxYhist = Math.ceil(maxYhist*100)/100.0;
//System.out.println("maxYhist after = "+maxYhist);
			double yPlotMax = 1.3*maxYhist;
//System.out.println("yPlotMax = "+yPlotMax);
			double yOffset=maxYhist*1.03;
//System.out.println("yOffset = "+yOffset);
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
//				System.out.println(branchName+"\t"+branchValue+"\t"+(float)wt+"\t"+(float)mean+"\t"+meanFromHist+"\t"+stdevFromHist);
			}
			
			DefaultXY_DataSet totMeanLine = new DefaultXY_DataSet();
			totMeanLine.set(totMeanEAL_here,0.0);
			totMeanLine.set(totMeanEAL_here,1.0);
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
			funcs.addAll(histArrayList);
//			funcs.add(meanLinesPlatform);

			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.addAll(plotCharsHists);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
			plotChars.addAll(plotCharsLines);
			plotChars.addAll(plotCharsValuePDFs);
//			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));

			PlotSpec specCombined = new PlotSpec(funcs, plotChars, branchName, "Normalized Loss", "Density");
			specCombined.setLegendVisible(false);
			specCombined.setLegendLocation(RectangleEdge.RIGHT);
			
			if(dirName != null) {
				String fname = branchName.replace(" ", "_");
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setUserBounds(xMinPlot, xMaxPlot, 0d, yPlotMax);
				gp.setTickLabelFontSize(23);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(26);
				gp.setBackgroundColor(Color.WHITE);
				gp.drawGraphPanel(specCombined);
				try {
//					File file = new File(dirName, fname);
//					gp.getChartPanel().setSize(1000, 800);
//					gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//					gp.saveAsPNG(file.getAbsolutePath() + ".png");
//					gp.saveAsTXT(file.getAbsolutePath() + ".txt");
					File file = new File(dirName, fname);
					gp.getChartPanel().setSize(400, 300);
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//					gp.saveAsPNG(file.getAbsolutePath()+".png");				

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
//				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
				gw.setAxisRange(xMinPlot, xMaxPlot, 0d, yPlotMax);
				gw.setTickLabelFontSize(19);
				gw.setAxisLabelFontSize(20);
				gw.setPlotLabelFontSize(21);
//				gw.setBackgroundColor(Color.WHITE);	
			}

			if(dirName != null) {
				String fname = branchName.replace(" ", "_")+"_wLegend";
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setUserBounds(xMinPlot, xMaxPlot, 0d, yPlotMax);
				gp.setTickLabelFontSize(23);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(26);
				gp.setBackgroundColor(Color.WHITE);
				gp.drawGraphPanel(specCombinedLegend);
				try {
//					File file = new File(dirName, fname);
//					gp.getChartPanel().setSize(1000, 800);
//					gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//					gp.saveAsPNG(file.getAbsolutePath() + ".png");
//					gp.saveAsTXT(file.getAbsolutePath() + ".txt");
					File file = new File(dirName, fname);
					gp.getChartPanel().setSize(500, 300);
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//					gp.saveAsPNG(file.getAbsolutePath()+".png");				

				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		System.out.println("\nClosest branches to totalMean EAL:\n");
		for(String branchName:closestValuesForBranchMap.keySet()) {
			System.out.println(branchName+"\t"+closestValuesForBranchMap.get(branchName));
		}
	}
	
	
	
	
	

	
	public void old_writeBranchValueStats(String dirName) {
		
		ArrayList<String> outputLinesList = new ArrayList<String>();
		String headerLine = "branchName\tbranchValue\tweight\tmean\tmeanDiff\twtedMeanDiff\tmeanIfRemoved\tmeanIfRemovedDiff";	
		outputLinesList.add(headerLine);
		
		// for branch-value averages
		weightForBrValMap = new HashMap<String, Double>();
		meanForBrValMap = new HashMap<String, Double>();
		meanDiffForBrValMap = new HashMap<String, Double>();
		meanDiffWtedForBrValMap = new HashMap<String, Double>();
		meanIfBrRemovedMap = new HashMap<String, Double>();
		meanDiffIfBrRemovedMap = new HashMap<String, Double>();
		
		for(String branchName:allBranchValuesMap.keySet()) {
			String[] branchValuesArray = allBranchValuesMap.get(branchName);
			HashMap<String, Double> meanHashMap = new HashMap<String, Double>();
			HashMap<String, Double> meanDiffHashMap = new HashMap<String, Double>();
			HashMap<String, Double> wtedMeanDiffHashMap = new HashMap<String, Double>();
			HashMap<String, Double> weightHashMap = new HashMap<String, Double>();
			for(int i=0; i<branchNormEAL.length; i++) {
				String branchValue = branchValuesArray[i];
				HistogramFunction hist;
				if(!meanHashMap.keySet().contains(branchValue)) {
					meanHashMap.put(branchValue, branchNormEAL[i]*branchWt[i]);
					weightHashMap.put(branchValue, branchWt[i]);
				}
				else {
					double newMean = meanHashMap.get(branchValue)+branchNormEAL[i]*branchWt[i];
					meanHashMap.put(branchValue, newMean);
					double newWt = weightHashMap.get(branchValue)+branchWt[i];
					weightHashMap.put(branchValue,newWt);
				}
			}
			if(meanHashMap.keySet().size() == 1) {
				continue;
			}
			
			// finish weights and means
			double totWt =0;
			for(Double wt:weightHashMap.values())
				totWt += wt;
			double meanTest=0d;
			for(String branchValue: meanHashMap.keySet()) {
				double mean = meanHashMap.get(branchValue)/weightHashMap.get(branchValue);
				meanHashMap.put(branchValue,mean);
				double wt = weightHashMap.get(branchValue)/totWt;
				meanTest+=mean*wt;
				weightHashMap.put(branchValue, wt);
				meanDiffHashMap.put(branchValue, mean-1.0);
				wtedMeanDiffHashMap.put(branchValue, (mean-1.0)*wt);
			}
			
			if(Math.abs(meanTest-1.0) >0.0001)
				throw new RuntimeException("something wrong with mean value");
			
			HashMap<String, Double> meanIfRemovedHashMap = new HashMap<String, Double>();
			HashMap<String, Double> meanIfRemovedDiffHashMap = new HashMap<String, Double>();
			for(String branchValue: meanHashMap.keySet()) {
				double wtTotal=0;
				double meanWithout=0;
				for(String branchValue2: meanHashMap.keySet()) {
					if(!branchValue2.equals(branchValue)) {
						wtTotal += weightHashMap.get(branchValue2);
						meanWithout += meanHashMap.get(branchValue2)*weightHashMap.get(branchValue2);
					}
				}
				meanWithout /= wtTotal;
				meanIfRemovedHashMap.put(branchValue, meanWithout);
				meanIfRemovedDiffHashMap.put(branchValue, meanWithout-1.0);
			}


			
			for(String branchValue: meanHashMap.keySet()) {
				String line = branchName+"\t"+branchValue+"\t"+weightHashMap.get(branchValue).floatValue()+
						"\t"+meanHashMap.get(branchValue).floatValue()+"\t"+meanDiffHashMap.get(branchValue).floatValue()+"\t"+
						wtedMeanDiffHashMap.get(branchValue).floatValue()+"\t"+meanIfRemovedHashMap.get(branchValue).floatValue()+"\t"+
						meanIfRemovedDiffHashMap.get(branchValue).floatValue();	
				outputLinesList.add(line);
			}	
			
			// save to total combined hashmap with combined name keys
			for(String branchValue: meanHashMap.keySet()) {
				String combinedName = branchName+" = "+branchValue;
				weightForBrValMap.put(combinedName, weightHashMap.get(branchValue));
				meanForBrValMap.put(combinedName, meanHashMap.get(branchValue));
				meanDiffForBrValMap.put(combinedName, meanDiffHashMap.get(branchValue));
				meanDiffWtedForBrValMap.put(combinedName, wtedMeanDiffHashMap.get(branchValue));
				meanIfBrRemovedMap.put(combinedName, meanIfRemovedHashMap.get(branchValue));
				meanDiffIfBrRemovedMap.put(combinedName, meanIfRemovedDiffHashMap.get(branchValue));
			}
		}
		
//		for(String line:	 outputLinesList)		
//			System.out.println(line);
		
		if(dirName != null) {
			try{
				FileWriter fw = new FileWriter(this.ROOT_DIR+dirName+"/branchValueStats.txt");
				for(String line:	 outputLinesList)		
					fw.write(line+"\n");
				fw.close();
			}catch(Exception e) {
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
	public ArrayList<String> getLogicTreeBranchValueOrder() {
		ArrayList<String> list = new ArrayList<String>();
		list.add("Fault Model = FM3_1");
		list.add("Fault Model = FM3_2");
		list.add("Deformation Model = GEOL");
		list.add("Deformation Model = ABM");
		list.add("Deformation Model = NEOK");
		list.add("Deformation Model = ZENGBB");
		list.add("Scaling Relationship = Shaw09Mod");
		list.add("Scaling Relationship = EllB");
		list.add("Scaling Relationship = HB08");
		list.add("Scaling Relationship = EllBsqrtLen");
		list.add("Scaling Relationship = ShConStrDrp");
		list.add("Slip Along Rupture Model (Dsr) = Tap");
		list.add("Slip Along Rupture Model (Dsr) = Uni");
		list.add("Total Mag 5 Rate = 6.5");
		list.add("Total Mag 5 Rate = 7.9");
		list.add("Total Mag 5 Rate = 9.6");
		list.add("MMax Off Fault = 7.3");
		list.add("MMax Off Fault = 7.6");
		list.add("MMax Off Fault = 7.9");
		list.add("Spatial Seismicity PDF = U2");
		list.add("Spatial Seismicity PDF = U3");
		list.add("ERF Probability Model = LOW_VALUES");
		list.add("ERF Probability Model = MID_VALUES");
		list.add("ERF Probability Model = HIGH_VALUES");
		list.add("ERF Probability Model = POISSON");
		list.add("Vs30 Model = Wills2015");
		list.add("Vs30 Model = WaldAllen");
		list.add("Ground Motion Model = ASK_2014");
		list.add("Ground Motion Model = BSSA_2014");
		list.add("Ground Motion Model = CB_2014");
		list.add("Ground Motion Model = CY_2014");
		list.add("GMM Additional Epistemic Uncertainty = LOWER");
		list.add("GMM Additional Epistemic Uncertainty = NONE");
		list.add("GMM Additional Epistemic Uncertainty = UPPER");
		return list;
	}
	
	public void makeFractionChangePlot(boolean popUpWindows, boolean savePlots, HashMap<String, Double> fractChangeMap, String fileName) {
		
		HashMap<String, Color> colorMap = new HashMap<String, Color>();
		colorMap.put("Fault Model", Color.BLACK);
		colorMap.put("Deformation Model", Color.GREEN);
		colorMap.put("Scaling Relationship", new Color(0, 0, 255));
		colorMap.put("Slip Along Rupture Model (Dsr)", new Color(80, 80, 255));
		colorMap.put("Total Mag 5 Rate", new Color(160, 160, 255));
		colorMap.put("MMax Off Fault", new Color(128, 0, 255));
		colorMap.put("Spatial Seismicity PDF", new Color(0, 128, 255));
		colorMap.put("ERF Probability Model", new Color(200, 150, 0));
		colorMap.put("Ground Motion Model", Color.RED);
		colorMap.put("GMM Additional Epistemic Uncertainty", Color.MAGENTA);
		colorMap.put("Vs30 Model", Color.ORANGE);

		
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
		
		double maxDiff = Math.max(-sortedVauesList.get(0), sortedVauesList.get(sortedVauesList.size()-1));
		double xMax = Math.ceil(maxDiff*10.0)/10.0;
		
		//override to the order in logic-tree figure
		sortedNamesList = getLogicTreeBranchValueOrder();
		Collections.reverse(sortedNamesList);
		sortedVauesList = new ArrayList<Double>();
		for(String nameVal:sortedNamesList)
			sortedVauesList.add(fractChangeMap.get(nameVal));
		
		
		
//		// override to TORNADO TYPE PLOT:
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
			xyFunc.set(0.0, yVal);
			xyFunc.setName(sortedNamesList.get(i));
			funcList.add(xyFunc);
			String[] branchName = sortedNamesList.get(i).split(" = ");
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 10f, colorMap.get(branchName[0])));
		}
		
		String stringForInfo="";
		for(int i=sortedNamesList.size()-1;i>=0;i--) {
			stringForInfo += sortedNamesList.get(i)+"\n";
		}
		funcList.get(0).setInfo(stringForInfo);

		
		double xMinPlot = -xMax;
		double xMaxPlot = xMax;
		double yMinPlot = 0;
		double yMaxPlot = sortedNamesList.size();

		//add 10% change lines
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

		

		PlotSpec specCombined = new PlotSpec(funcList, plotChars, fileName, "Fractional Change", "Branch");
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
//				File file = new File(dirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(400, 500);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
		// now with legend
		PlotSpec specCombinedLegend = new PlotSpec(funcList, plotChars, fileName, "Fractional Change", "Branch");
		specCombinedLegend.setLegendVisible(true);
		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);

		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(specCombinedLegend);
//			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			gw.setAxisRange(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
//			gw.setBackgroundColor(Color.WHITE);	
		}

		if(savePlots) {
			String fname = fileName+"_wLegend";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xMinPlot, xMaxPlot, yMinPlot, yMaxPlot);
			gp.setTickLabelFontSize(23);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(22);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(specCombinedLegend);
			try {
//				File file = new File(dirName, fname);
//				gp.getChartPanel().setSize(1000, 800);
//				gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//				gp.saveAsPNG(file.getAbsolutePath() + ".png");
//				gp.saveAsTXT(file.getAbsolutePath() + ".txt");
				File file = new File(ROOT_DIR+outDirName, fname);
				gp.getChartPanel().setSize(800, 500);
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		
	}
	
	

	public void makeEAL_Historgram(boolean popUpWindows, boolean savePlots) {
		double min = 0;
		double max =6;
		int num = 120;
		
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
		
		double factor95percent = get95percConfFactorForValue(1.0, histCumulative);
		infoString += "\n95% Conf Bounds Factor About Mean: " + (float)factor95percent+"\n";
		double within10percent = histCumulative.getInterpolatedY(1.0*1.1) - histCumulative.getInterpolatedY(1.0*0.9);
		infoString += "Fract within 10% of Mean: " + (float)within10percent+" (chance it's outside: "+(float)(1.0-within10percent)+")\n";

		
		// Median info
		totMedianNormEAL = histCumulative.getFirstInterpolatedX(0.5);
		infoString += "\nMedian (normalized) from Cumulative Dist: " + (float)+totMedianNormEAL+"\n";
		int minDiffWtedIndex = getBranchClosestToValue(totMedianNormEAL, true);
		factor95percent = get95percConfFactorForValue(totMedianNormEAL, histCumulative);
		infoString += "95% Conf Bounds Factor About Median: " + (float)factor95percent+"\n";
		within10percent = histCumulative.getInterpolatedY(totMedianNormEAL*1.1) - histCumulative.getInterpolatedY(totMedianNormEAL*0.9);
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
		
		PlotSpec spec = new PlotSpec(funcsArray, plotChars, "Normalized EAL Distribution", "EAL (Billion $)", "Density");
//		specCombinedLegend.setLegendVisible(false);
//		specCombinedLegend.setLegendLocation(RectangleEdge.RIGHT);

		if(popUpWindows) {
			GraphWindow gw = new GraphWindow(spec);
			gw.setAxisRange(0d, 3d, 0d, 1.5);
			gw.setTickLabelFontSize(19);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
		}

		if(savePlots) {
			String fname = "EAL_Histogram";
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(0d, 3d, 0d, 1.5);
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
//				gp.saveAsPNG(file.getAbsolutePath()+".png");				

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

		PlotSpec specCum = new PlotSpec(funcsArray2, plotChars2, "Normalized Cumulative EAL Distribution", "EAL (Billion $)", "Density");
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
//				gp2.saveAsTXT(file.getAbsolutePath() + ".txt");
				gp2.getChartPanel().setSize(500, 400);
				gp2.saveAsPDF(file.getAbsolutePath()+".pdf");
//				gp2.saveAsPNG(file.getAbsolutePath()+".png");				

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
	private double get95percConfFactorForValue(double normEAL_value, HistogramFunction cumDist) {
		
		double maxX = cumDist.getMaxX()/normEAL_value;
		EvenlyDiscretizedFunc fractWithinScaleFactFunc = new EvenlyDiscretizedFunc(1.1,maxX,200);
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
		


	public static void main(String[] args) {
		
		Analysis analysis = new Analysis("all_branch_results.csv", "fullRegionResults");
		
	}

}
