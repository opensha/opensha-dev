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
	HashMap<String, String[]> branchValuesHashMap;
	
	double totMeanEAL;
	double totStdDevEAL;
	double totWeight;
	HashMap<String, String> closestValuesForBranchMap; // this stores the branch value that has a mean closest to the totalMeanEAL
	
	// for branch-value averages
	HashMap<String, Double> weightHashMap;
	HashMap<String, Double> meanHashMap;
	HashMap<String, Double> meanDiffHashMap;
	HashMap<String, Double> wtedMeanDiffHashMap;
	HashMap<String, Double> meanIfRemovedHashMap;
	HashMap<String, Double> meanIfRemovedDiffHashMap;

	
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
		
		// make directory if it doesn't exist
		File dir = new File(ROOT_DIR+outDirName);
		if(!dir.exists())
			dir.mkdirs();

		infoString += "Input File Name: "+inputFileName+"\n";
		infoString += "Output Dir Name: "+outDirName+"\n";
		
//		readBranchLevelSummaryDataFromFile();
		readAllBranchDataFromFile();
		
		writeBranchValueStats(outDirName);

		
		makeEAL_Historgram(popUpPlots, savePlots);
		
		System.out.println(infoString);
		
		
//		makeEAL_Historgram(true, ROOT_DIR+"stackedHistograms");
//		makeLog10_EAL_Historgram();

		
//		makeStackedHistograms(true, ROOT_DIR+"stackedHistograms");
		
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

//		make_ERF_vs_GMM_UncertHists(true, ROOT_DIR+"stackedHistograms");
	}
	
	/**
	 * This reads the branch level data
	 * @return
	 */
	private void readBranchLevelSummaryDataFromFile() {
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
			branchWt = new double[fileLines.size()-1];
			branchEAL = new double[fileLines.size()-1];
			branchNormEAL = new double[fileLines.size()-1];
			branchValuesHashMap = new HashMap<String, String[]>();
			String str = fileLines.get(0);
			String[] colNamesArray = str.split(",");
			for(int c=5;c<18; c++) {
				String keyName = colNamesArray[c];
				String[] strArray = new String[fileLines.size()-1];
				branchValuesHashMap.put(keyName, strArray);
				infoString += "\t"+keyName+"\n";
//				System.out.println(keyName);
			}
			
			double minEAL = Double.MAX_VALUE;
			double maxEAL = 0;
			int maxLineNum=-1;
			int minLineNum=-1;
			
			totWeight = 0;
			for(int i=1; i<fileLines.size();i++ ) {
				str = fileLines.get(i);
				String[] split = str.split(",");
				branchWt[i-1] = Double.parseDouble(split[1]);
				branchEAL[i-1] = Double.parseDouble(split[2])*1e-6;	// convert to Billions
				totWeight+=branchWt[i-1];
				if(minEAL > branchEAL[i-1]) {
					minEAL = branchEAL[i-1];
					minLineNum = i;
				}
				if(maxEAL < branchEAL[i-1]) {
					maxEAL = branchEAL[i-1];
					maxLineNum = i;
				}
				// set branch values
				for(int c=5;c<18; c++) {
					String keyName = colNamesArray[c];
					String[] valArray = branchValuesHashMap.get(keyName);
					valArray[i-1] = split[c];
				}
			}
			
			totMeanEAL = computeWeightedAverage(branchWt, branchEAL);
//			totMeanEAL = Precision.round(totMeanEAL, 3);
			totStdDevEAL = computeWeightedStdDev(branchWt, branchEAL, totMeanEAL);
//			totStdDevEAL = Precision.round(totStdDevEAL, 3);
			
			infoString += "\ntotMeanEAL = "+ Precision.round(totMeanEAL, 3);
			infoString += "\ntotStdDevEAL = "+ Precision.round(totStdDevEAL, 3);
			infoString += "\ntotWeight = "+(float)totWeight;
			infoString += "\nsize = "+branchWt.length;
			infoString += "\nMinEAL:\t"+(float)minEAL+"; fileline:";
			infoString += "\n\t"+fileLines.get(minLineNum);
			infoString += "\nMaxEAL:\t"+(float)maxEAL+"; fileline:";
			infoString += "\n\t"+fileLines.get(maxLineNum);

			
			for(int i=0;i<branchNormEAL.length;i++)
				branchNormEAL[i] = branchEAL[i]/totMeanEAL;
			
//			System.out.println("totMeanEAL:\t"+(float)totMeanEAL);
//			System.out.println("totStdDevEAL:\t"+(float)totStdDevEAL);
//			System.out.println("totWeight:\t"+(float)totWeight);
//			System.out.println("size:\t"+branchWt.length);
//			System.out.println("MinEAL:\t"+(float)minEAL);
//			System.out.println("\t"+fileLines.get(minLineNum));
//			System.out.println("MaxEAL:\t"+(float)maxEAL);
//			System.out.println("\t"+fileLines.get(maxLineNum));

			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
		
	
	public void makeStackedHistograms(boolean popUpWindows, String dirName) {
				
		double min = 0;
		double max = 20;
		int num = 80;
		
		double xMinPlot=0;
		double xMaxPlot=10;
		
		closestValuesForBranchMap = new HashMap<String, String> ();
		
		for(String branchName:branchValuesHashMap.keySet()) {
//System.out.println("working on" + branchName);
			String[] branchValuesArray = branchValuesHashMap.get(branchName);
			HashMap<String, HistogramFunction> histogramsHashMap = new HashMap<String, HistogramFunction>();
			HashMap<String, Double> meanHashMap = new HashMap<String, Double>();
			HashMap<String, Double> weightHashMap = new HashMap<String, Double>();
			for(int i=0; i<branchEAL.length; i++) {
				String branchValue = branchValuesArray[i];
				HistogramFunction hist;
				if(!histogramsHashMap.keySet().contains(branchValue)) {
					hist = new HistogramFunction(min+(max-min)/(num*2),max-(max-min)/(num*2), num);
					hist.setName(branchValue);
					histogramsHashMap.put(branchValue, hist);
					meanHashMap.put(branchValue, branchEAL[i]*branchWt[i]);
					weightHashMap.put(branchValue, branchWt[i]);
				}
				else {
					hist = histogramsHashMap.get(branchValue);
					double newMean = meanHashMap.get(branchValue)+branchEAL[i]*branchWt[i];
					meanHashMap.put(branchValue, newMean);
					double newWt = weightHashMap.get(branchValue)+branchWt[i];
					weightHashMap.put(branchValue,newWt);
				}
				hist.add(branchEAL[i], branchWt[i]);
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
			for(String branchValue: histogramsHashMap.keySet()) {
				double mean = meanHashMap.get(branchValue)/weightHashMap.get(branchValue);
				meanHashMap.put(branchValue,mean);
				double wt = weightHashMap.get(branchValue)/totWt;
				weightHashMap.put(branchValue, wt);
				HistogramFunction hist = histogramsHashMap.get(branchValue);
				double scaleCorr = wt / (hist.calcSumOfY_Vals()*hist.getDelta());
				hist.scale(scaleCorr);
//				double test = hist.calcSumOfY_Vals()*hist.getDelta();
//System.out.println("test: "+(float)(wt/test));
				if(maxHistValue <hist.getMaxY())
					maxHistValue = hist.getMaxY();
				double meanDiff = Math.abs(totMeanEAL-mean);
				if(minMeanDiff>meanDiff) {
					minMeanDiff=meanDiff;
					closestValue = branchValue;
				}
			}
			
			closestValuesForBranchMap.put(branchName, closestValue);
			
			
			ArrayList<HistogramFunction> histArrayList = new ArrayList<HistogramFunction>();
			for(HistogramFunction hist:histogramsHashMap.values())
				histArrayList.add(hist);
//System.out.println(histArrayList.size()+"; "+ histogramsHashMap.keySet().size());
			List<HistogramFunction> stackedHistograms = HistogramFunction.getStackedHists(histArrayList, true);
			
			
			
			HashMap<String, Color> colorMap = this.getColorMap(meanHashMap);
			ArrayList<PlotCurveCharacterstics> plotCharsHists = new ArrayList<PlotCurveCharacterstics>();
			ArrayList<PlotCurveCharacterstics> plotCharsLines = new ArrayList<PlotCurveCharacterstics>();

			for(String branchOption: meanHashMap.keySet()) {
//System.out.println(branchName+"\t"+branchOption+"\t"+colorMap.get(branchOption)+"\t"+meanHashMap.get(branchOption));
				plotCharsHists.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colorMap.get(branchOption)));
				plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colorMap.get(branchOption)));
			}

			
			ArrayList<DefaultXY_DataSet> meanLinesList = new ArrayList<DefaultXY_DataSet>();
			double totMeanEAL_here=0;
			// the following are to plot the mean lines offset from the graph
			double yPlotMax = 0.16;
			double yOffset=0.11;
			double yScale=(yPlotMax-yOffset)/0.8;	// 0.8 is the max branch weight (MMax)
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
//				System.out.println(branchName+"\t"+branchValue+"\t"+(float)wt+"\t"+(float)mean+"\t"+meanFromHist+"\t"+stdevFromHist);
			}
			
			DefaultXY_DataSet totMeanLine = new DefaultXY_DataSet();
			totMeanLine.set(totMeanEAL_here,0.0);
			totMeanLine.set(totMeanEAL_here,1.0);
			totMeanLine.setName("Mean");
			totMeanLine.setInfo("totMean="+(float)totMeanEAL_here);
			meanLinesList.add(totMeanLine);
			plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
			
			// this is for combined stack hists and means
			ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
			funcs.addAll(stackedHistograms);
			funcs.addAll(meanLinesList);
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.addAll(plotCharsHists);
			plotChars.addAll(plotCharsLines);

			PlotSpec specCombined = new PlotSpec(funcs, plotChars, branchName, "EAL (Billion $)", "Density");
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
			PlotSpec specCombinedLegend = new PlotSpec(funcs, plotChars, branchName, "EAL (Billion $)", "Density");
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
	
	
	
	
	
	public void oldMakeStackedHistograms(boolean popUpWindows, String dirName) {
		
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdirs();
		
		double min = 0;
		double max = 20;
		int num = 80;
		
		double xMinPlot=0;
		double xMaxPlot=12;
		
		for(String branchName:branchValuesHashMap.keySet()) {
//System.out.println("working on" + branchName);
			String[] branchValuesArray = branchValuesHashMap.get(branchName);
			HashMap<String, HistogramFunction> histogramsHashMap = new HashMap<String, HistogramFunction>();
			HashMap<String, Double> meanHashMap = new HashMap<String, Double>();
			HashMap<String, Double> weightHashMap = new HashMap<String, Double>();
			for(int i=0; i<branchEAL.length; i++) {
				String branchValue = branchValuesArray[i];
				HistogramFunction hist;
				if(!histogramsHashMap.keySet().contains(branchValue)) {
					hist = new HistogramFunction(min+(max-min)/(num*2),max-(max-min)/(num*2), num);
					hist.setName(branchValue);
					histogramsHashMap.put(branchValue, hist);
					meanHashMap.put(branchValue, branchEAL[i]*branchWt[i]);
					weightHashMap.put(branchValue, branchWt[i]);
				}
				else {
					hist = histogramsHashMap.get(branchValue);
					double newMean = meanHashMap.get(branchValue)+branchEAL[i]*branchWt[i];
					meanHashMap.put(branchValue, newMean);
					double newWt = weightHashMap.get(branchValue)+branchWt[i];
					weightHashMap.put(branchValue,newWt);
				}
				hist.add(branchEAL[i], branchWt[i]);
			}
			if(histogramsHashMap.keySet().size() == 1) {
				continue;
			}
			
			// finish weights and means and correct histograms
			double maxHistValue=0;
			double totWt =0;
			for(Double wt:weightHashMap.values())
				totWt += wt;
			for(String branchValue: histogramsHashMap.keySet()) {
				double mean = meanHashMap.get(branchValue)/weightHashMap.get(branchValue);
				meanHashMap.put(branchValue,mean);
				double wt = weightHashMap.get(branchValue)/totWt;
				weightHashMap.put(branchValue, wt);
				HistogramFunction hist = histogramsHashMap.get(branchValue);
				double scaleCorr = wt / (hist.calcSumOfY_Vals()*hist.getDelta());
				hist.scale(scaleCorr);
				double test = hist.calcSumOfY_Vals()*hist.getDelta();
//System.out.println("test: "+(float)(wt/test));
				if(maxHistValue <hist.getMaxY())
					maxHistValue = hist.getMaxY();
			}
			
			
			ArrayList<HistogramFunction> histArrayList = new ArrayList<HistogramFunction>();
			for(HistogramFunction hist:histogramsHashMap.values())
				histArrayList.add(hist);
//System.out.println(histArrayList.size()+"; "+ histogramsHashMap.keySet().size());
			List<HistogramFunction> stackedHistograms = HistogramFunction.getStackedHists(histArrayList, true);
			
			
			
			HashMap<String, Color> colorMap = this.getColorMap(meanHashMap);
			ArrayList<PlotCurveCharacterstics> plotCharsHists = new ArrayList<PlotCurveCharacterstics>();
			ArrayList<PlotCurveCharacterstics> plotCharsLines = new ArrayList<PlotCurveCharacterstics>();

			for(String branchOption: meanHashMap.keySet()) {
//System.out.println(branchName+"\t"+branchOption+"\t"+colorMap.get(branchOption)+"\t"+meanHashMap.get(branchOption));
				plotCharsHists.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colorMap.get(branchOption)));
				plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, colorMap.get(branchOption)));
			}

			PlotSpec spec = new PlotSpec(stackedHistograms, plotCharsHists, branchName, "EAL (Billion $)", "Density");
			spec.setLegendVisible(false);
			spec.setLegendLocation(RectangleEdge.RIGHT);
			
			if(popUpWindows) {
				GraphWindow gw = new GraphWindow(spec);
//				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
				gw.setAxisRange(xMinPlot, xMaxPlot, 0d, 0.1);
				gw.setTickLabelFontSize(18);
				gw.setAxisLabelFontSize(20);
				gw.setPlotLabelFontSize(21);
//				gw.setBackgroundColor(Color.WHITE);	
			}

			if(dirName != null) {
				String fname = branchName.replace(" ", "_");
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setUserBounds(xMinPlot, xMaxPlot, 0d, 0.12);
				gp.setTickLabelFontSize(20);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.setBackgroundColor(Color.WHITE);
				gp.drawGraphPanel(spec);
				try {
//					File file = new File(dirName, fname);
//					gp.getChartPanel().setSize(1000, 800);
//					gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//					gp.saveAsPNG(file.getAbsolutePath() + ".png");
//					gp.saveAsTXT(file.getAbsolutePath() + ".txt");
					File file = new File(dirName, fname);
					gp.getChartPanel().setSize(500, 200);
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//					gp.saveAsPNG(file.getAbsolutePath()+".png");				

				} catch (IOException e) {
					e.printStackTrace();
				}
				
				// now make plots with legend
				PlotSpec spec2 = new PlotSpec(stackedHistograms, plotCharsHists, branchName, "EAL (Billion $)", "Density");
				spec2.setLegendVisible(true);
				spec2.setLegendLocation(RectangleEdge.RIGHT);
				HeadlessGraphPanel gp2 = new HeadlessGraphPanel();
				gp2.setUserBounds(xMinPlot, xMaxPlot, 0d, 0.12);
				gp2.setTickLabelFontSize(18);
				gp2.setAxisLabelFontSize(20);
				gp2.setPlotLabelFontSize(21);
				gp2.setBackgroundColor(Color.WHITE);
				gp2.drawGraphPanel(spec2);
				try {
//					File file = new File(dirName, fname);
//					gp2.getChartPanel().setSize(1000, 800);
//					gp2.saveAsPDF(file.getAbsolutePath() + ".pdf");
//					gp2.saveAsPNG(file.getAbsolutePath() + ".png");
//					gp2.saveAsTXT(file.getAbsolutePath() + ".txt");
					File file = new File(dirName, fname+"_wLegend");
					gp2.getChartPanel().setSize(500, 200);
					gp2.saveAsPDF(file.getAbsolutePath()+".pdf");
//					gp2.saveAsPNG(file.getAbsolutePath()+".png");				

				} catch (IOException e) {
					e.printStackTrace();
				}

			}
			
			
			ArrayList<DefaultXY_DataSet> meanLinesList = new ArrayList<DefaultXY_DataSet>();
			double totMean=0;
			for(String branchValue: histogramsHashMap.keySet()) {
				HistogramFunction hist = histogramsHashMap.get(branchValue);
				DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
				double mean = meanHashMap.get(branchValue);
				double wt = weightHashMap.get(branchValue);
				meanLine.set(mean,0.0);
				meanLine.set(mean,wt);
				meanLine.setName(hist.getName());
				double meanFromHist = hist.computeMean();
				double stdevFromHist = hist.computeStdDev();
				totMean += mean*wt;

				meanLine.setInfo("mean="+(float)mean+"\nweight="+(float)wt+"\tmeanFromHist="+(float)meanFromHist+
						"\nstdevFromHist="+(float)stdevFromHist);
				meanLinesList.add(meanLine);
//				System.out.println(branchName+"\t"+branchValue+"\t"+(float)wt+"\t"+(float)mean+"\t"+meanFromHist+"\t"+stdevFromHist);
			}
			
			DefaultXY_DataSet totMeanLine = new DefaultXY_DataSet();
			totMeanLine.set(totMean,0.0);
			totMeanLine.set(totMean,1.0);
			totMeanLine.setName("Mean");
			totMeanLine.setInfo("totMean="+(float)totMean);
			meanLinesList.add(totMeanLine);
			plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
			
			PlotSpec specMeans = new PlotSpec(meanLinesList, plotCharsLines, branchName, "EAL (Billion $)", "Weight");
			specMeans.setLegendVisible(false);
			specMeans.setLegendLocation(RectangleEdge.RIGHT);
			
			if(popUpWindows) {
				GraphWindow gw = new GraphWindow(specMeans);
//				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
				gw.setAxisRange(xMinPlot, xMaxPlot, 0d, 1.0);
				gw.setTickLabelFontSize(18);
				gw.setAxisLabelFontSize(20);
				gw.setPlotLabelFontSize(21);
			}

			if(dirName != null) {
				String fname = branchName.replace(" ", "_")+"_Means";
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setUserBounds(xMinPlot, xMaxPlot, 0d, 1.0);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.setBackgroundColor(Color.WHITE);
				gp.drawGraphPanel(specMeans);
				try {
//					File file = new File(dirName, fname);
//					gp.getChartPanel().setSize(1000, 800);
//					gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//					gp.saveAsPNG(file.getAbsolutePath() + ".png");
//					gp.saveAsTXT(file.getAbsolutePath() + ".txt");
					File file = new File(dirName, fname);
					gp.getChartPanel().setSize(500, 200);
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//					gp.saveAsPNG(file.getAbsolutePath()+".png");				

				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			
			// this is for combined stack hists and means
			ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
			funcs.addAll(stackedHistograms);
			
			
			double wtScale = 0.2;

			for(DefaultXY_DataSet xy_data: meanLinesList) {
				xy_data.scale(wtScale);
				funcs.add(xy_data);
				funcs.add(xy_data);
			}
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.addAll(plotCharsHists);
			for(String branchOption: meanHashMap.keySet()) {
				plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY));
				plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, colorMap.get(branchOption)));
			}
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));


			PlotSpec specCombinedd = new PlotSpec(funcs, plotChars, branchName, "EAL (Billion $)", "Density");
			specCombinedd.setLegendVisible(false);
			specCombinedd.setLegendLocation(RectangleEdge.RIGHT);
			
			if(popUpWindows) {
				GraphWindow gw = new GraphWindow(specCombinedd);
//				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
				gw.setAxisRange(0d, 12d, 0d, 0.15);
				gw.setTickLabelFontSize(18);
				gw.setAxisLabelFontSize(20);
				gw.setPlotLabelFontSize(21);
//				gw.setBackgroundColor(Color.WHITE);	
			}

			if(dirName != null) {
				String fname = branchName.replace(" ", "_")+"_Conbined";
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setUserBounds(0d, 12d, 0d, 0.2);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.setBackgroundColor(Color.WHITE);
				gp.drawGraphPanel(specCombinedd);
				try {
//					File file = new File(dirName, fname);
//					gp.getChartPanel().setSize(1000, 800);
//					gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
//					gp.saveAsPNG(file.getAbsolutePath() + ".png");
//					gp.saveAsTXT(file.getAbsolutePath() + ".txt");
					File file = new File(dirName, fname);
					gp.getChartPanel().setSize(500, 400);
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
//					gp.saveAsPNG(file.getAbsolutePath()+".png");				

				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			
			
			
			
//			// now may non-stacked, normalized histograms with means
//			// normalize histograms
//			for(HistogramFunction hist:histArrayList)
//				hist.scale(1/(hist.calcSumOfY_Vals()*hist.getDelta()));
//			for(DefaultXY_DataSet xy_data: meanLinesList)
//				xy_data.scale(1.0);
//			ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
//			funcs.addAll(histArrayList);
//			HistogramFunction totalHist = stackedHistograms.get(0);
//			totalHist.scale(1/(totalHist.calcSumOfY_Vals()*totalHist.getDelta()));
//			totalHist.setName("totalHist");
//			totalHist.setInfo("");
//			funcs.add(totalHist);
//			funcs.addAll(meanLinesList);
//			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
//			plotChars.addAll(plotCharsLines);
//			plotChars.addAll(plotCharsLines);
//
//			PlotSpec specNonStacked = new PlotSpec(funcs, plotChars, branchName, "EAL (Billion $)", "Density");
//			specNonStacked.setLegendVisible(true);
//			specNonStacked.setLegendLocation(RectangleEdge.RIGHT);
//			
//			if(popUpWindows) {
//				GraphWindow gw = new GraphWindow(specNonStacked);
////				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
//				gw.setAxisRange(0d, 12d, 0d, 1.0);
//				gw.setTickLabelFontSize(18);
//				gw.setAxisLabelFontSize(20);
//				gw.setPlotLabelFontSize(21);
////				gw.setBackgroundColor(Color.WHITE);	
//			}

				


			
		}
	}

	
	public void writeBranchValueStats(String dirName) {
		
		ArrayList<String> outputLinesList = new ArrayList<String>();
		String headerLine = "branchName\tbranchValue\tweight\tmean\tmeanDiff\twtedMeanDiff\tmeanIfRemoved\tmeanIfRemovedDiff";	
		outputLinesList.add(headerLine);

		
		for(String branchName:branchValuesHashMap.keySet()) {
			String[] branchValuesArray = branchValuesHashMap.get(branchName);
			meanHashMap = new HashMap<String, Double>();
			meanDiffHashMap = new HashMap<String, Double>();
			wtedMeanDiffHashMap = new HashMap<String, Double>();
			weightHashMap = new HashMap<String, Double>();
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
				throw new RuntimeException("somethig wrong with mean value");
			
			meanIfRemovedHashMap = new HashMap<String, Double>();
			meanIfRemovedDiffHashMap = new HashMap<String, Double>();
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
		}
		
		for(String line:	 outputLinesList)		
			System.out.println(line);
		
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

	
	public void make_ERF_vs_GMM_UncertHists(boolean popUpWindows, String dirName) {
		
		double min = 0;
		double max = 20;
		int num = 40;
		
		// make ERF histogram
		double mean_erf =0;
		double totWt_erf=0;
		HistogramFunction hist_erf = new HistogramFunction(0.0+max/(num*2),max-max/(num*2), num);
		for(int i=0; i<branchEAL.length; i++) {
				// use only gmm branches that are the closest one to the mean
				boolean useThis = true;
				for(String gmmBranchName:gmm_branches) {
					String branchValue = branchValuesHashMap.get(gmmBranchName)[i];
					String closestValue = closestValuesForBranchMap.get(gmmBranchName);
					if(!branchValue.equals(closestValue))
						useThis = false;
				}
				if(useThis) {
					hist_erf.add(branchEAL[i], branchWt[i]);
					mean_erf += branchEAL[i]*branchWt[i];
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
		for(int i=0; i<branchEAL.length; i++) {
			// use only erf branches that are the closest one to the mean
				boolean useThis = true;
				for(String erfBranchName:erf_branches) {
					String branchValue = branchValuesHashMap.get(erfBranchName)[i];
					String closestValue = this.closestValuesForBranchMap.get(erfBranchName);
					if(!branchValue.equals(closestValue))
						useThis = false;
				}
				if(useThis) {
					hist_gmm.add(branchEAL[i], branchWt[i]);
					mean_gmm += branchEAL[i]*branchWt[i];
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
		
		HistogramFunction hist_tot = new HistogramFunction(0.0+max/(num*2),max-max/(num*2), num);
		for(int i=0;i<hist_tot.size();i++) {
			hist_tot.set(i, hist_erf.getY(i)*totWt_erf+hist_gmm.getY(i)*totWt_gmm);
		}
		hist_tot.scale(1.0/(hist_tot.getDelta()*hist_tot.calcSumOfY_Vals()));
		hist_tot.setName("hist_tot");


		ArrayList<HistogramFunction> histArrayList = new ArrayList<HistogramFunction>();
		histArrayList.add(hist_erf);
		histArrayList.add(hist_gmm);
		histArrayList.add(hist_tot);
		
		ArrayList<PlotCurveCharacterstics> plotCharsLines = new ArrayList<PlotCurveCharacterstics>();
		plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotCharsLines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));

		
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

		if(dirName != null) {
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
				File file = new File(dirName, fname);
				gp.getChartPanel().setSize(500, 300);
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
		int num = 60;
		
		HistogramFunction hist = new HistogramFunction(0.0+(max-min)/(num*2),max-(max-min)/(num*2), num);
		for(int i=0; i<branchNormEAL.length; i++) {
			hist.add(branchNormEAL[i], branchWt[i]);
		}
		
		hist.normalizeBySumOfY_Vals();
		hist.scale(1.0/hist.getDelta());
		
		String histName = "Normalized PDF; orig mean = "+this.totMeanEAL+" and orig stdDev = "+this.totStdDevEAL;
		hist.setName(histName);
		
		double meanFromPDF = hist.computeMean();
		String info = "\nMean from Norm PDF = "+(float)meanFromPDF+"\n";
		double stdDevFromPDF = hist.computeStdDev();
		info += "StdDev from Norm PDF: " + (float)stdDevFromPDF+"\n";
		info += "COV from data Norm PDF = " + (float)(stdDevFromPDF/meanFromPDF)+"\n";
		info += "COV from orig (non-norm) data = " + (float)(totStdDevEAL/totMeanEAL)+"\n";
		
		HistogramFunction histCumulative = hist.getCumulativeDistFunctionWithHalfBinOffset();
		histCumulative.scale(hist.getDelta());
		histCumulative.setName("Cumulative Dist");
		
		info += "Median from Cumulative Dist: " + (float)histCumulative.getFirstInterpolatedX(0.5)+"\n";
		double low95 = (float)histCumulative.getFirstInterpolatedX(0.025);
		double hi95 = (float)histCumulative.getFirstInterpolatedX(0.975);
		info += "95% Conf Bounds from Cumulative Dist: " + (float)low95+", "+(float)hi95+"\n";
		double within10percent = histCumulative.getInterpolatedY(1.1) - histCumulative.getInterpolatedY(0.9);
		info += "Fract within 10% of mean: " + (float)within10percent+"\n";
		info += "Likelihood it's more than 10% from mean: " + (float)(1.0-within10percent)+"\n";
		hist.setInfo(info);
		histCumulative.setInfo(info);
		infoString += "\n"+info;
		
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
