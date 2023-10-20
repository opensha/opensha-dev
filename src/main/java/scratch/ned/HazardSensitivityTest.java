package scratch.ned;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;

import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagAreaRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.rupForecastImpl.FloatingPoissonFaultERF;
import org.opensha.sha.earthquake.rupForecastImpl.PoissonFaultERF;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SingleMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.param.SimpleFaultParameter;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.alessandro.FaultSystemRuptureRateInversion;

public class HazardSensitivityTest {
	
	static double hazCurveLnMin = Math.log(0.001);
	static double hazCurveLnMax = Math.log(10);
	static int hazCurveNum = 50;
	static double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);
	
	
	/**
	 */
	private static EvenlyDiscretizedFunc getHazardCurveLnX(ERF erf, Location location, double saPeriod) {
		
		boolean D= false;
		
		// write out parameter values
		if(D) {
			ParameterList paramList = erf.getAdjustableParameterList();
			for(int i=0;i<paramList.size(); i++) {
				Parameter<?> param = paramList.getByIndex(i);
				System.out.println(param.getName()+"\t"+param.getValue());
			}			
		}

		
		// chose attenuation relationship (GMPE)
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		imr.setParamDefaults();
		
		EvenlyDiscretizedFunc curveLogXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		
		// make the site object and set values
		Site site = new Site(location);
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
//			System.out.println(param.getName()+"\t"+param.getValue());
		}
		
		// set the IMT
		ArbitrarilyDiscretizedFunc curveLinearXvalues;
		if(saPeriod == 0) {
			imr.setIntensityMeasure(PGA_Param.NAME);
		}
		else {
			SA_Param saParam = (SA_Param)imr.getParameter(SA_Param.NAME);
			saParam.getPeriodParam().setValue(saPeriod);
			imr.setIntensityMeasure(saParam);
		}
	
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
		
		return curveLogXvalues;
	}

	

	public static String compareHazardCurves1(ERF[] erfArray, Location location, double saPeriod, String dirName, 
			boolean popupWindow, String plotTitle, double duration) {
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		

		EvenlyDiscretizedFunc[] curveLogXvaluesArray = new EvenlyDiscretizedFunc[erfArray.length];
		for(int i=0;i<erfArray.length;i++) {
			curveLogXvaluesArray[i] = getHazardCurveLnX(erfArray[i], location, saPeriod);
		}


		// convert to linear x valules
		double[] twoIn50valueArray = new double[erfArray.length];
		double[] tenIn50valueArray = new double[erfArray.length];
		ArbitrarilyDiscretizedFunc[] curveLinearXvaluesArray = new ArbitrarilyDiscretizedFunc[erfArray.length];
		for(int i=0;i<erfArray.length;i++) {
			curveLinearXvaluesArray[i] = new ArbitrarilyDiscretizedFunc();
			for(int j=0;j<curveLogXvaluesArray[i].size();j++) {
				curveLinearXvaluesArray[i].set(Math.exp(curveLogXvaluesArray[i].getX(j)), curveLogXvaluesArray[i].getY(j));
			}
			twoIn50valueArray[i] = Math.exp(curveLogXvaluesArray[i].getFirstInterpolatedX(0.02));
			tenIn50valueArray[i] = Math.exp(curveLogXvaluesArray[i].getFirstInterpolatedX(0.1));
			curveLinearXvaluesArray[i].setInfo("2in50 value: "+twoIn50valueArray[i]+"\n"+"10in50 value: "+tenIn50valueArray[i]+
					"\nLocation: "+location.getLatitude()+", "+location.getLongitude());

		}

		ArbitrarilyDiscretizedFunc twoPercentLine = new ArbitrarilyDiscretizedFunc();
		twoPercentLine.set(curveLinearXvaluesArray[0].getX(0),0.02);
		twoPercentLine.set(curveLinearXvaluesArray[0].getX(curveLinearXvaluesArray[0].size()-1),0.02);
		ArbitrarilyDiscretizedFunc tenPercentLine = new ArbitrarilyDiscretizedFunc();
		tenPercentLine.set(curveLinearXvaluesArray[0].getX(0),0.1);
		tenPercentLine.set(curveLinearXvaluesArray[0].getX(curveLinearXvaluesArray[0].size()-1),0.1);
		
		double minTwoIn50 = Double.MAX_VALUE;
		double maxTwoIn50 = 0d;
		double minTenIn50 = Double.MAX_VALUE;
		double maxTenIn50 = 0d;
		for(int i=0;i<twoIn50valueArray.length;i++) {
			if(minTwoIn50>twoIn50valueArray[i]) minTwoIn50=twoIn50valueArray[i];
			if(maxTwoIn50<twoIn50valueArray[i]) maxTwoIn50=twoIn50valueArray[i];
			if(minTenIn50>tenIn50valueArray[i]) minTenIn50=tenIn50valueArray[i];
			if(maxTenIn50<tenIn50valueArray[i]) maxTenIn50=tenIn50valueArray[i];
		}
		
		
			
		// ratio between the max and min mag
//		double twoIn50valueRatio = twoIn50valueArray[twoIn50valueArray.length-1]/twoIn50valueArray[0];
//		double tenIn50valueRatio = tenIn50valueArray[tenIn50valueArray.length-1]/tenIn50valueArray[0];
		double twoIn50valueRatio = maxTwoIn50/minTwoIn50;
		double tenIn50valueRatio = maxTenIn50/minTenIn50;
		
		double covTwoIn50 = ((maxTwoIn50-minTwoIn50)/Math.sqrt(12)) / (0.5*(maxTwoIn50+minTwoIn50));	//stdDev/mean
		double covTenIn50 = ((maxTenIn50-minTenIn50)/Math.sqrt(12)) / (0.5*(maxTenIn50+minTenIn50));	//stdDev/mean

	
		String returnString = (float)twoIn50valueRatio + "\t"+(float)tenIn50valueRatio+ "\t"+(float)covTwoIn50+ "\t"+(float)covTenIn50;
		
		// make the plot
		plottingFuncsArray.add(twoPercentLine);
		plottingFuncsArray.add(tenPercentLine);
		for(int i=0;i<erfArray.length;i++)
			plottingFuncsArray.add(curveLinearXvaluesArray[i]);
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.ORANGE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));

		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/hazardCurve";
		String xAxisLabel = imtString;
		String yAxisLabel = "Probability (in "+duration+" yr)";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = true;
		boolean logY = true;

		writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		return returnString;

	}
	
	
	
	public static String compareHazardCurves2(ERF[] erfArray, ERF grERF, Location location, double saPeriod, String dirName, 
			boolean popupWindow, String plotTitle, double duration) {
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		EvenlyDiscretizedFunc grCurveLogXvalues = getHazardCurveLnX(grERF, location, saPeriod);
		ArbitrarilyDiscretizedFunc grCurveLinearXvalues = new ArbitrarilyDiscretizedFunc();
		for(int j=0;j<grCurveLogXvalues.size();j++) {
			grCurveLinearXvalues.set(Math.exp(grCurveLogXvalues.getX(j)), grCurveLogXvalues.getY(j));
		}
		double grTwoIn50value = Math.exp(grCurveLogXvalues.getFirstInterpolatedX(0.02));
		double grTenIn50value = Math.exp(grCurveLogXvalues.getFirstInterpolatedX(0.1));
		grCurveLinearXvalues.setInfo("grTwoIn50value: "+grTwoIn50value+"\n"+"grTenIn50value: "+grTenIn50value);
		grCurveLinearXvalues.setName("grCurveLinearXvalues");


		EvenlyDiscretizedFunc[] curveLogXvaluesArray = new EvenlyDiscretizedFunc[erfArray.length];
		for(int i=0;i<erfArray.length;i++) {
			curveLogXvaluesArray[i] = getHazardCurveLnX(erfArray[i], location, saPeriod);
		}

		// convert to linear x valules
		double[] twoIn50valueArray = new double[erfArray.length];
		double[] tenIn50valueArray = new double[erfArray.length];
		ArbitrarilyDiscretizedFunc[] curveLinearXvaluesArray = new ArbitrarilyDiscretizedFunc[erfArray.length];
		for(int i=0;i<erfArray.length;i++) {
			curveLinearXvaluesArray[i] = new ArbitrarilyDiscretizedFunc();
			for(int j=0;j<curveLogXvaluesArray[i].size();j++) {
				curveLinearXvaluesArray[i].set(Math.exp(curveLogXvaluesArray[i].getX(j)), curveLogXvaluesArray[i].getY(j));
			}
			twoIn50valueArray[i] = Math.exp(curveLogXvaluesArray[i].getFirstInterpolatedX(0.02));
			tenIn50valueArray[i] = Math.exp(curveLogXvaluesArray[i].getFirstInterpolatedX(0.1));
			curveLinearXvaluesArray[i].setInfo("2in50 value: "+twoIn50valueArray[i]+"\n"+"10in50 value: "+tenIn50valueArray[i]+
					"\nLocation: "+location.getLatitude()+", "+location.getLongitude());

		}

		ArbitrarilyDiscretizedFunc twoPercentLine = new ArbitrarilyDiscretizedFunc();
		twoPercentLine.set(curveLinearXvaluesArray[0].getX(0),0.02);
		twoPercentLine.set(curveLinearXvaluesArray[0].getX(curveLinearXvaluesArray[0].size()-1),0.02);
		ArbitrarilyDiscretizedFunc tenPercentLine = new ArbitrarilyDiscretizedFunc();
		tenPercentLine.set(curveLinearXvaluesArray[0].getX(0),0.1);
		tenPercentLine.set(curveLinearXvaluesArray[0].getX(curveLinearXvaluesArray[0].size()-1),0.1);
		
		double minTwoIn50 = Double.MAX_VALUE;
		double maxTwoIn50 = 0d;
		double minTenIn50 = Double.MAX_VALUE;
		double maxTenIn50 = 0d;
		for(int i=0;i<twoIn50valueArray.length;i++) {
			if(minTwoIn50>twoIn50valueArray[i]) minTwoIn50=twoIn50valueArray[i];
			if(maxTwoIn50<twoIn50valueArray[i]) maxTwoIn50=twoIn50valueArray[i];
			if(minTenIn50>tenIn50valueArray[i]) minTenIn50=tenIn50valueArray[i];
			if(maxTenIn50<tenIn50valueArray[i]) maxTenIn50=tenIn50valueArray[i];
		}
		
		
			
		// ratio between the max and min mag
//		double twoIn50valueRatio = twoIn50valueArray[twoIn50valueArray.length-1]/twoIn50valueArray[0];
//		double tenIn50valueRatio = tenIn50valueArray[tenIn50valueArray.length-1]/tenIn50valueArray[0];
		double twoIn50valueRatio = maxTwoIn50/minTwoIn50;
		double tenIn50valueRatio = maxTenIn50/minTenIn50;

			
		// frat diff between the max and min mag
		double twoIn50valueMinMag = (twoIn50valueArray[0]/grTwoIn50value);
		double tenIn50valueMinMag = (tenIn50valueArray[0]/grTenIn50value);
		double twoIn50valueMaxMag = (twoIn50valueArray[twoIn50valueArray.length-1]/grTwoIn50value);
		double tenIn50valueMaxMag = (tenIn50valueArray[tenIn50valueArray.length-1]/grTenIn50value);
	
		String returnString = (float)twoIn50valueMinMag + "\t"+(float)twoIn50valueMaxMag + 
				"\t"+(float)tenIn50valueMinMag + "\t"+(float)tenIn50valueMaxMag +
				"\t"+(float)twoIn50valueRatio + "\t"+(float)tenIn50valueRatio;
		
		// make the plot
		plottingFuncsArray.add(twoPercentLine);
		plottingFuncsArray.add(tenPercentLine);
		plottingFuncsArray.add(grCurveLinearXvalues);
		for(int i=0;i<erfArray.length;i++)
			plottingFuncsArray.add(curveLinearXvaluesArray[i]);
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.ORANGE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));

		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/hazardCurve";
		String xAxisLabel = imtString;
		String yAxisLabel = "Probability (in "+duration+" yr)";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = true;
		boolean logY = true;

		writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		return returnString;

	}
	
	
	public static String compareHazardCurves3(ERF erf1, ERF erf2, Location location, double saPeriod, String dirName, 
			boolean popupWindow, String plotTitle, double duration) {
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		EvenlyDiscretizedFunc curve1_LogXvalues = getHazardCurveLnX(erf1, location, saPeriod);
		ArbitrarilyDiscretizedFunc curve1_LinearXvalues = new ArbitrarilyDiscretizedFunc();
		for(int j=0;j<curve1_LogXvalues.size();j++) {
			curve1_LinearXvalues.set(Math.exp(curve1_LogXvalues.getX(j)), curve1_LogXvalues.getY(j));
		}
//		double grTwoIn50value = Math.exp(curve1_LogXvalues.getFirstInterpolatedX(0.02));
//		double grTenIn50value = Math.exp(curve1_LogXvalues.getFirstInterpolatedX(0.1));
//		curve1_LinearXvalues.setInfo("grTwoIn50value: "+grTwoIn50value+"\n"+"grTenIn50value: "+grTenIn50value);
		curve1_LinearXvalues.setName("erf1 curve");
		
		
		EvenlyDiscretizedFunc curve2_LogXvalues = getHazardCurveLnX(erf2, location, saPeriod);
		ArbitrarilyDiscretizedFunc curve2_LinearXvalues = new ArbitrarilyDiscretizedFunc();
		for(int j=0;j<curve2_LogXvalues.size();j++) {
			curve2_LinearXvalues.set(Math.exp(curve2_LogXvalues.getX(j)), curve2_LogXvalues.getY(j));
		}
//		double grTwoIn50value = Math.exp(curve2_LogXvalues.getFirstInterpolatedX(0.02));
//		double grTenIn50value = Math.exp(curve2_LogXvalues.getFirstInterpolatedX(0.1));
//		curve2_LinearXvalues.setInfo("grTwoIn50value: "+grTwoIn50value+"\n"+"grTenIn50value: "+grTenIn50value);
		curve2_LinearXvalues.setName("erf1 curve");


		ArbitrarilyDiscretizedFunc twoPercentLine = new ArbitrarilyDiscretizedFunc();
		twoPercentLine.set(curve2_LinearXvalues.getX(0),0.02);
		twoPercentLine.set(curve2_LinearXvalues.getX(curve2_LinearXvalues.size()-1),0.02);
		ArbitrarilyDiscretizedFunc tenPercentLine = new ArbitrarilyDiscretizedFunc();
		tenPercentLine.set(curve2_LinearXvalues.getX(0),0.1);
		tenPercentLine.set(curve2_LinearXvalues.getX(curve2_LinearXvalues.size()-1),0.1);
		
		
//		String returnString = (float)twoIn50valueMinMag + "\t"+(float)twoIn50valueMaxMag + 
//				"\t"+(float)tenIn50valueMinMag + "\t"+(float)tenIn50valueMaxMag +
//				"\t"+(float)twoIn50valueRatio + "\t"+(float)tenIn50valueRatio;
		
		String returnString = "";
				
		// make the plot
		plottingFuncsArray.add(twoPercentLine);
		plottingFuncsArray.add(tenPercentLine);
		plottingFuncsArray.add(curve1_LinearXvalues);
		plottingFuncsArray.add(curve2_LinearXvalues);
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.ORANGE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));

		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/hazardCurve";
		String xAxisLabel = imtString;
		String yAxisLabel = "Probability (in "+duration+" yr)";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = true;
		boolean logY = true;

		writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		return returnString;

	}



	
	/**
	 * The general x-y plotting method
	 * @param funcs
	 * @param plotChars
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	protected static void writeAndOrPlotFuncs(
			ArrayList<XY_DataSet> funcs, 
			ArrayList<PlotCurveCharacterstics> plotChars, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			String fileNamePrefix, 
			boolean popupWindow) {
		
		
		if(popupWindow) {
			GraphWindow graph = new GraphWindow(funcs, plotName);
			if(xAxisRange != null)
				graph.setX_AxisRange(xAxisRange.getLowerBound(),xAxisRange.getUpperBound());
			if(yAxisRange != null)
				graph.setY_AxisRange(yAxisRange.getLowerBound(),yAxisRange.getUpperBound());
			graph.setXLog(logX);
			graph.setYLog(logY);
			graph.setPlotChars(plotChars);
			graph.setX_AxisLabel(xAxisLabel);
			graph.setY_AxisLabel(yAxisLabel);
			graph.setTickLabelFontSize(22);
			graph.setAxisLabelFontSize(24);
			
			if(fileNamePrefix != null) {
				try {
					graph.saveAsPDF(fileNamePrefix+".pdf");
					graph.saveAsPNG(fileNamePrefix+".png");
					graph.saveAsTXT(fileNamePrefix+".txt");
				} catch (IOException e) {
					e.printStackTrace();
				}
			}		
		}
		else if (fileNamePrefix != null){	// no pop-up window; just save plot
			PlotSpec spec = new PlotSpec(funcs, plotChars, plotName, xAxisLabel, yAxisLabel);
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(24);
			gp.setAxisLabelFontSize(26);
			gp.drawGraphPanel(spec, logX, logY);
			gp.getChartPanel().setSize(1000, 800);
			
			try {
				gp.saveAsPNG(fileNamePrefix+".png");
				gp.saveAsPDF(fileNamePrefix+".pdf");
				gp.saveAsTXT(fileNamePrefix+".txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}
	
	
	/**
	 * put all moment into one magnitude comparison
	 */
	public static void doIt1() {
		
		double[] siteDistKmArray = {0, 10};
		double[] saPeriod = {0.0, 1.0};
		double mMax = 8.3;
//		double[] magArray = {7.8, mMax};
		double[] magArray = {6.3, 6.7, 7.1, 7.6, 8.0, mMax};

		
		double fltWidth = 12e3;
		double slipRate = 30*1e-3;	// 30 mm/yr = 0.03 m/yr
		double rake = 0;
		// set the forecast duration
		double duration = 50;

		String resultsString = "";
		
//		Location centerLoc = new Location(0.0-120,0.0);
		
		// Fault attributes
		Location locCenter = new Location(0.0,-120.0,0.0);
		double gridSpacing = 1.0;
		ArrayList<Double> lats = new ArrayList<Double>();
		ArrayList<Double> lons; 		// these will be solved for
		ArrayList<Double> dips = new ArrayList<Double>();
		ArrayList<Double> depths = new ArrayList<Double>();
		lats.add(locCenter.getLatitude());
		lats.add(locCenter.getLatitude());
		dips.add(90.);
		depths.add(0.);
		depths.add(12.);		
		// lons are solved for below
		
		String dirName = null;
		boolean popupWindow = true;

		ArrayList<ScalingRelationships> scalingRelList = new ArrayList<ScalingRelationships>();
		scalingRelList.add(ScalingRelationships.ELLSWORTH_B);
		scalingRelList.add(ScalingRelationships.HANKS_BAKUN_08);
		scalingRelList.add(ScalingRelationships.SHAW_2009_MOD);
		scalingRelList.add(ScalingRelationships.ELLB_SQRT_LENGTH);
		scalingRelList.add(ScalingRelationships.SHAW_CONST_STRESS_DROP);
		

		int ithPlot = 1;
		for(int i=0;i<scalingRelList.size();i++) {

			System.out.println(scalingRelList.get(i).getName());
			double maxLengthKm = 1e-3*scalingRelList.get(i).getArea(mMax, fltWidth)/fltWidth;	
			
			double moRate = FaultMomentCalc.getMoment(maxLengthKm*1e3*fltWidth, slipRate);
			System.out.println("\n\tmoRate = "+moRate);
			
			double maxLengthDegrees = (maxLengthKm/111.2);	
			lons = new ArrayList<Double>();
			lons.add(locCenter.getLongitude()-maxLengthDegrees/2.0);
			lons.add(locCenter.getLongitude()+maxLengthDegrees/2.0);

			// test
			Location loc1 = new Location(lats.get(0),lons.get(0),0.0);
			Location loc2 = new Location(lats.get(1),lons.get(1),0.0);
			double testDist = LocationUtils.horzDistance(loc1, loc2);
			//			System.out.println(loc1);
			//			System.out.println(loc2);
			System.out.println("\tLenth for M="+mMax+": "+(float)maxLengthKm);
			System.out.println("\tTestLength (ratio): "+(float)testDist+"  ("+(float)(maxLengthKm/testDist)+")");
			//			System.exit(-1);

			// make list of ERFs
			double[] tempRate = new double[magArray.length];
			ERF[] erfArray = new ERF[magArray.length];
			for(int m=0;m<magArray.length;m++) {
				double mag=magArray[m];
				double area = scalingRelList.get(i).getArea(mag, fltWidth);
				double length = area/fltWidth;
				double slip = scalingRelList.get(i).getAveSlip(area, length, fltWidth, fltWidth, rake);
				length *= 1e-3; // convert to km
				double eventRate = (slipRate/slip)*(maxLengthKm/length);
				tempRate[m] = eventRate;
				System.out.println("\n\tmag="+mag+"\n\t\tlength="+(float)length+
						"\n\t\tslip="+(float)slip+"\n\t\trate="+(float)eventRate+"\n");

				// create ERF
				//						PoissonFaultERF erf = new PoissonFaultERF();
				FloatingPoissonFaultERF erf = new FloatingPoissonFaultERF();
				erf.getParameter(FloatingPoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
				SingleMagFreqDist mfd = new SingleMagFreqDist(mag,2,0.1, mag,1.0);
				mfd.scaleToIncrRate(mag, eventRate);
				erf.getParameter(FloatingPoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(mfd);
				erf.getParameter(FloatingPoissonFaultERF.OFFSET_PARAM_NAME).setValue(5.0);
				erf.getParameter(FloatingPoissonFaultERF.FLOATER_TYPE_PARAM_NAME).setValue(FloatingPoissonFaultERF.FLOATER_TYPE_FULL_DDW);
				SimpleFaultParameter faultParam = (SimpleFaultParameter)erf.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
				faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
				faultParam.setEvenlyGriddedSurfaceFromParams();
				erf.getTimeSpan().setDuration(duration);
				erf.updateForecast();
				System.out.println("\t\tmoRate = "+(float)ERF_Calculator.getTotalMomentRateInRegion(erf, null));
				erfArray[m] = erf;
				System.out.println("\t\tnumSrc = "+erf.getNumSources()+"\n\t\tnumRupForFirstSrc = "+erf.getSource(0).getNumRuptures());
			}
			
			for(int p=0;p<saPeriod.length;p++) {
				for(int d=0;d<siteDistKmArray.length;d++) {
					Location location = new Location(siteDistKmArray[d]/111.2, -120,0);	// 
					String title = scalingRelList.get(i).getShortName()+"; Per="+saPeriod[p]+"; Dist="+siteDistKmArray[d];
					String ratioString = compareHazardCurves1(erfArray, location, saPeriod[p], dirName, popupWindow, title, duration);

					if(magArray.length == 2) {
						System.out.println("\tRateRatio="+(float)(tempRate[0]/tempRate[1]));
						System.out.println("\t"+ratioString+"\n");
					}
					resultsString += ithPlot+"\t"+scalingRelList.get(i).getShortName()+"\t"+saPeriod[p]+"\t"+siteDistKmArray[d]+"\t"+ratioString+"\n";		
					ithPlot += 1;
				}
			}
		}
		System.out.println(resultsString);

	}

	
	/**
	 * Adjust b-value of U3 on-fault MFD and preserve moment rates
	 */
	public static void doit2() {
		
		double[] siteDistKmArray = {0, 10};
		double[] saPeriod = {0.0, 1.0};
//		double[] bValueChangeArray = {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
		double[] bValueChangeArray = {-2.0, -1.0, 0.0, 1.0, 2.0};

		
		double fltWidth = 12e3;
		double rake = 0;
		// set the forecast duration
		double duration = 50;

		String resultsString = "";
		
		// Fault attributes
		Location locCenter = new Location(0.0,-120.0,0.0);
		double gridSpacing = 1.0;
		ArrayList<Double> lats = new ArrayList<Double>();
		ArrayList<Double> lons; 		// these will be solved for
		ArrayList<Double> dips = new ArrayList<Double>();
		ArrayList<Double> depths = new ArrayList<Double>();
		lats.add(locCenter.getLatitude());
		lats.add(locCenter.getLatitude());
		dips.add(90.);
		depths.add(0.);
		depths.add(12.);		
		// lons are solved for below
		
		String dirName = null;
		boolean popupWindow = true;
		
		SummedMagFreqDist origMFD = getU3_LongTermOnFaultIncrementalMFD(false);
		double mMax = origMFD.getMaxX();

		int ithPlot = 1;

		WC1994_MagAreaRelationship magAreaRel = new WC1994_MagAreaRelationship();
		double areaKm = magAreaRel.getMedianArea(mMax);
		
		double maxLengthKm = areaKm/(fltWidth*1e-3);	
			
		double moRate = origMFD.getTotalMomentRate();
		System.out.println("\n\tmoRate = "+moRate+";\tmaxLengthKm="+maxLengthKm);
			
		double maxLengthDegrees = (maxLengthKm/111.2);	
		lons = new ArrayList<Double>();
		lons.add(locCenter.getLongitude()-maxLengthDegrees/2.0);
		lons.add(locCenter.getLongitude()+maxLengthDegrees/2.0);

		// test
		Location loc1 = new Location(lats.get(0),lons.get(0),0.0);
		Location loc2 = new Location(lats.get(1),lons.get(1),0.0);
		double testDist = LocationUtils.horzDistance(loc1, loc2);
		System.out.println("\tTestLength (ratio): "+(float)testDist+"  ("+(float)(maxLengthKm/testDist)+")");

		FloatingPoissonFaultERF origERF = new FloatingPoissonFaultERF();
		origERF.getParameter(FloatingPoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
		origERF.getParameter(FloatingPoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(origMFD);
		origERF.getParameter(FloatingPoissonFaultERF.OFFSET_PARAM_NAME).setValue(5.0);
		origERF.getParameter(FloatingPoissonFaultERF.FLOATER_TYPE_PARAM_NAME).setValue(FloatingPoissonFaultERF.FLOATER_TYPE_FULL_DDW);
		SimpleFaultParameter faultParam = (SimpleFaultParameter)origERF.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
		faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
		faultParam.setEvenlyGriddedSurfaceFromParams();
		origERF.getTimeSpan().setDuration(duration);
		origERF.updateForecast();
		// make list of ERFs
		ERF[] erfArray = new ERF[bValueChangeArray.length];
		ArrayList<SummedMagFreqDist> mfdList = new ArrayList<SummedMagFreqDist>();
		int erfIndex=0;
		for(double bValChange : bValueChangeArray) {
			SummedMagFreqDist mfd = modify_bValueOfMFD(origMFD, bValChange);
			mfd.setInfo("bValua change ="+bValChange);
			mfdList.add(mfd);
			// create ERF
			FloatingPoissonFaultERF erf = new FloatingPoissonFaultERF();
			erf.getParameter(FloatingPoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
			erf.getParameter(FloatingPoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(mfd);
			erf.getParameter(FloatingPoissonFaultERF.OFFSET_PARAM_NAME).setValue(5.0);
			erf.getParameter(FloatingPoissonFaultERF.FLOATER_TYPE_PARAM_NAME).setValue(FloatingPoissonFaultERF.FLOATER_TYPE_FULL_DDW);
			faultParam = (SimpleFaultParameter)erf.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
			faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
			faultParam.setEvenlyGriddedSurfaceFromParams();
			erf.getTimeSpan().setDuration(duration);
			erf.updateForecast();
			erfArray[erfIndex] = erf;
			erfIndex+=1;
//			System.out.println("\t\tnumSrc = "+erf.getNumSources()+"\n\t\tnumRupForFirstSrc = "+erf.getSource(0).getNumRuptures());
		}
		
			
			for(int p=0;p<saPeriod.length;p++) {
				for(int d=0;d<siteDistKmArray.length;d++) {
					Location location = new Location(siteDistKmArray[d]/111.2, -120,0);	// 
					String title = "Per="+saPeriod[p]+"; Dist="+siteDistKmArray[d];
			//		String ratioString = compareHazardCurves1(erfArray, location, saPeriod[p], dirName, popupWindow, title, duration);
					String ratioString = compareHazardCurves2(erfArray, origERF, location, saPeriod[p], dirName, 
							popupWindow, title, duration);
					resultsString += ithPlot+"\t"+saPeriod[p]+"\t"+siteDistKmArray[d]+"\t"+ratioString+"\n";		
					ithPlot += 1;
				}
			}
			
			GraphWindow graph = new GraphWindow(mfdList, "MFDs"); 
			graph.setX_AxisLabel("Magnitude");
			graph.setY_AxisLabel("Incremental Rate (per yr)");
			graph.setYLog(true);


		System.out.println(resultsString);

	}

	
	/**
	 * compare different MFDs (with same moment rate) for a single-fault-based ERF
	 */
	public static void doit3() {
		
		double[] siteDistKmArray = {0, 10};
		double[] saPeriod = {0.0, 1.0};

		
		double fltWidth = 12e3;
		double rake = 0;
		// set the forecast duration
		double duration = 50;

		String resultsString = "";
		
		// Fault attributes
		Location locCenter = new Location(0.0,-120.0,0.0);
		double gridSpacing = 1.0;
		ArrayList<Double> lats = new ArrayList<Double>();
		ArrayList<Double> lons; 		// these will be solved for
		ArrayList<Double> dips = new ArrayList<Double>();
		ArrayList<Double> depths = new ArrayList<Double>();
		lats.add(locCenter.getLatitude());
		lats.add(locCenter.getLatitude());
		dips.add(90.);
		depths.add(0.);
		depths.add(12.);		
		// lons are solved for below
		
		String dirName = null;
		boolean popupWindow = true;
		
		SummedMagFreqDist mfd1 = getU3_SanJacintoNuclIncrMFD();
		double mMax = mfd1.getMaxX();
		SummedMagFreqDist mfd2 = getU2_SanJacintoNuclIncrMFD();
		
		double moRate = mfd1.getTotalMomentRate();
		mfd1.scaleToTotalMomentRate(moRate*10);
		mfd2.scaleToTotalMomentRate(moRate*10);

		double moRateRatio = mfd1.getTotalMomentRate()/mfd2.getTotalMomentRate();
		System.out.println("moment rate ratio: "+(float)moRateRatio);
		
		ArrayList<EvenlyDiscretizedFunc> mfdList = new ArrayList<EvenlyDiscretizedFunc>();
		mfdList.add(mfd1);
		mfdList.add(mfd2);
		mfdList.add(mfd1.getCumRateDistWithOffset());
		mfdList.add(mfd2.getCumRateDistWithOffset());
		GraphWindow graph = new GraphWindow(mfdList, "MFDs"); 
		graph.setX_AxisLabel("Magnitude");
		graph.setY_AxisLabel("Incremental Rate (per yr)");
		graph.setYLog(true);


		int ithPlot = 1;

		WC1994_MagAreaRelationship magAreaRel = new WC1994_MagAreaRelationship();
		double areaKm = magAreaRel.getMedianArea(mMax);
		
		double maxLengthKm = areaKm/(fltWidth*1e-3);	
			
		System.out.println("\n\tmoRate = "+moRate+";\tmaxLengthKm="+maxLengthKm);

		double maxLengthDegrees = (maxLengthKm/111.2);	
		lons = new ArrayList<Double>();
		lons.add(locCenter.getLongitude()-maxLengthDegrees/2.0);
		lons.add(locCenter.getLongitude()+maxLengthDegrees/2.0);

		// test
		Location loc1 = new Location(lats.get(0),lons.get(0),0.0);
		Location loc2 = new Location(lats.get(1),lons.get(1),0.0);
		double testDist = LocationUtils.horzDistance(loc1, loc2);
		System.out.println("\tTestLength (ratio): "+(float)testDist+"  ("+(float)(maxLengthKm/testDist)+")");

		PoissonFaultERF origERF = new PoissonFaultERF();
		origERF.getParameter(PoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
		origERF.getParameter(PoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(mfd1);
		SimpleFaultParameter faultParam = (SimpleFaultParameter)origERF.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
		faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
		faultParam.setEvenlyGriddedSurfaceFromParams();
		origERF.getTimeSpan().setDuration(duration);
		origERF.updateForecast();
		System.out.println("\t\tnumSrc = "+origERF.getNumSources()+"\n\t\tnumRupForFirstSrc = "+origERF.getSource(0).getNumRuptures());

		
		PoissonFaultERF altERF = new PoissonFaultERF();
		altERF.getParameter(PoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
		altERF.getParameter(PoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(mfd2);
		SimpleFaultParameter faultParam2 = (SimpleFaultParameter)altERF.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
		faultParam2.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
		faultParam2.setEvenlyGriddedSurfaceFromParams();
		altERF.getTimeSpan().setDuration(duration);
		altERF.updateForecast();
		System.out.println("\t\tnumSrc = "+altERF.getNumSources()+"\n\t\tnumRupForFirstSrc = "+altERF.getSource(0).getNumRuptures());

		for(int p=0;p<saPeriod.length;p++) {
				for(int d=0;d<siteDistKmArray.length;d++) {
					Location location = new Location(siteDistKmArray[d]/111.2, -120,0);	// 
					String title = "Per="+saPeriod[p]+"; Dist="+siteDistKmArray[d];
			//		String ratioString = compareHazardCurves1(erfArray, location, saPeriod[p], dirName, popupWindow, title, duration);
					String ratioString = compareHazardCurves3(altERF, origERF, location, saPeriod[p], dirName, 
							popupWindow, title, duration);
					resultsString += ithPlot+"\t"+saPeriod[p]+"\t"+siteDistKmArray[d]+"\t"+ratioString+"\n";		
					ithPlot += 1;
				}
			}
			


		System.out.println(resultsString);

	}

	
	
	public static void old_doIt2() {
		
		double[] siteDistKmArray = {0, 10};
		double[] saPeriod = {0.0, 1.0};
		double[] magArray = {6.3, 8.0};
//		double[] magArray = {6.3, 6.7, 7.1, 7.6, 8.0};

		
		double fltWidth = 12e3;
		double slipRate = 30*1e-3;	// 30 mm/yr
		double rake = 0;
		// set the forecast duration
		double duration = 50;

		String resultsString = "";
		
		// Fault attributes
		double gridSpacing = 1.0;
		ArrayList lats = new ArrayList();
		ArrayList lons = new ArrayList();
		ArrayList dips = new ArrayList();
		ArrayList depths = new ArrayList();
		lats.add(36.);
		lats.add(36.);
		lons.add(-119.);
		lons.add(-121.);
		dips.add(90.);
		depths.add(0.);
		depths.add(12.);		
		
		String dirName = null;
		boolean popupWindow = true;

		ArrayList<ScalingRelationships> scalingRelList = new ArrayList<ScalingRelationships>();
		scalingRelList.add(ScalingRelationships.ELLSWORTH_B);
		scalingRelList.add(ScalingRelationships.HANKS_BAKUN_08);
		scalingRelList.add(ScalingRelationships.SHAW_2009_MOD);
		scalingRelList.add(ScalingRelationships.ELLB_SQRT_LENGTH);
		scalingRelList.add(ScalingRelationships.SHAW_CONST_STRESS_DROP);
		
//		for(int i=0;i<scalingRelList.size();i++) {
//			System.out.println(scalingRelList.get(i).getShortName()+"\t"+scalingRelList.get(i).getMag(12*12*1e6, 12*1e3));
//		}
//		System.exit(-1);


		for(int i=0;i<scalingRelList.size();i++) {
			ScalingRelationships scaleRel = scalingRelList.get(i);
			magArray[0] = scaleRel.getMag(12*12*1e6, 12*1e3, 12*1e3, 12*1e3, rake);
			for(int p=0;p<saPeriod.length;p++) {
				for(int d=0;d<siteDistKmArray.length;d++) {
					Location location = new Location(36+siteDistKmArray[d]/111, -120,0);	// 
					System.out.println(scaleRel.getName());
					double[] tempRate = new double[magArray.length];
					ERF[] erfArray = new ERF[magArray.length];
					
					
					// MAKE GR Dist
					double maxMag = magArray[magArray.length-1];
					double minMag = magArray[0];
					int num = (int)Math.ceil((maxMag-minMag)/0.1);
					double delta = (maxMag-minMag)/num;
					GutenbergRichterMagFreqDist grDist = new GutenbergRichterMagFreqDist(minMag, num+1, delta, minMag, maxMag, 1.0,0.0);
					double testSlipRate =0;
					for(int m=0;m<grDist.size();m++) {
						double mag=grDist.getX(m);
						double area = scaleRel.getArea(mag, fltWidth);
						double length = area/fltWidth;
						double slip = scaleRel.getAveSlip(area, length, fltWidth, fltWidth, rake);
						testSlipRate+=slip*grDist.getY(m);
					}
					grDist.scale(slipRate/testSlipRate);
					// test slip rates:
					testSlipRate =0;
					for(int m=0;m<grDist.size();m++) {
						double mag=grDist.getX(m);
						double area = scaleRel.getArea(mag, fltWidth);
						double length = area/fltWidth;
						double slip = scaleRel.getAveSlip(area, length, fltWidth, fltWidth, rake);
						testSlipRate+=slip*grDist.getY(m);
					}
					float ratio = (float)(slipRate/testSlipRate);
					if(ratio != 1)
						throw new RuntimeException(ratio+" should be 1.0");
					
					PoissonFaultERF grERF = new PoissonFaultERF();
					grERF.getParameter(PoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
					grERF.getParameter(PoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(grDist);
					SimpleFaultParameter faultParam = (SimpleFaultParameter)grERF.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
					faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
					faultParam.setEvenlyGriddedSurfaceFromParams();
					grERF.getTimeSpan().setDuration(duration);
					grERF.updateForecast();
//					System.out.println(grDist);
//					System.exit(-1);
					
					// Now make each ERF
					for(int m=0;m<magArray.length;m++) {
						double mag=magArray[m];
						double area = scaleRel.getArea(mag, fltWidth);
						double length = area/fltWidth;
						double slip = scaleRel.getAveSlip(area, length, fltWidth, fltWidth, rake);
						length *= 1e-3;
						double eventRate = slipRate/slip;
						tempRate[m] = eventRate;
						System.out.println("\n\tmag="+mag+"\n\t\tlength="+(float)length+
								"\n\t\tslip="+(float)slip+"\n\t\trate="+(float)eventRate+"\n");
						
						// create ERF
						PoissonFaultERF erf = new PoissonFaultERF();
						erf.getParameter(PoissonFaultERF.RAKE_PARAM_NAME).setValue(rake);
						SingleMagFreqDist mfd = new SingleMagFreqDist(mag,2,0.1, mag,1.0);
						mfd.scaleToIncrRate(mag, eventRate);
						erf.getParameter(PoissonFaultERF.MAG_DIST_PARAM_NAME).setValue(mfd);
						faultParam = (SimpleFaultParameter)erf.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
						faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
						faultParam.setEvenlyGriddedSurfaceFromParams();
						erf.getTimeSpan().setDuration(duration);
						erf.updateForecast();
						erfArray[m] = erf;
					}
					String title = scaleRel.getShortName()+"_"+siteDistKmArray[d]+"kmDist";
					String ratioString = compareHazardCurves2(erfArray, grERF, location, saPeriod[p], dirName, popupWindow, title, duration);

					if(magArray.length == 2) {
						System.out.println("\tRateRatio="+(float)(tempRate[0]/tempRate[1]));
						System.out.println("\t"+ratioString+"\n");
					}
					resultsString += title+"\t"+saPeriod[p]+"\t"+siteDistKmArray[d]+"\t"+ratioString+"\n";				
				}
			}
		}
		System.out.println(resultsString);

	}
	
	public static GutenbergRichterMagFreqDist getApproxGRdist(double[] magArray, double slipRate, ScalingRelationships scalingRel, double fltWidth) {
		// make the GR dist ERF
		double maxMag = magArray[magArray.length-1];
		double minMag = magArray[0];
		int num = (int)Math.ceil((maxMag-minMag)/0.1);
		double delta = (maxMag-minMag)/num;
		GutenbergRichterMagFreqDist grDist = new GutenbergRichterMagFreqDist(minMag, num+1, delta, minMag, maxMag, 1.0,0.0);
		double testSlipRate =0;
		for(int m=0;m<grDist.size();m++) {
			double mag=grDist.getX(m);
			double area = scalingRel.getArea(mag, fltWidth);
			double length = area/fltWidth;
			double slip = scalingRel.getAveSlip(area, length, fltWidth, fltWidth, Double.NaN);
			testSlipRate+=slip*grDist.getY(m);
		}
		grDist.scale(slipRate/testSlipRate);
		
		// test slip rates:
		testSlipRate =0;
		for(int m=0;m<grDist.size();m++) {
			double mag=grDist.getX(m);
			double area = scalingRel.getArea(mag, fltWidth);
			double length = area/fltWidth;
			double slip = scalingRel.getAveSlip(area, length, fltWidth, fltWidth, Double.NaN);
			testSlipRate+=slip*grDist.getY(m);
		}
		
		float ratio = (float)(slipRate/testSlipRate);
		if(ratio != 1)
			throw new RuntimeException(ratio+" should be 1.0");
		
		System.out.println(ratio);
		
		return grDist;

	}

	
	public static double[] epsilon_test(double mag, double rate, String dirName, boolean popupWindow) {
		
		double distance = 0;
		
		String infoString = "";
		
		Location loc = new Location(34.,-118.,0.);
		double rake = 0;
		double prob = 1;	// this is ingnored
		PointSurface pointSurface = new PointSurface(loc);
		pointSurface.setAveDip(90);
		ProbEqkRupture eqkRup = new ProbEqkRupture(mag, rake, prob, pointSurface, loc);
		
		ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
		imr.setParamDefaults();
		
		double hazCurveLnMin = Math.log(0.001);
		double hazCurveLnMax = Math.log(20);
		int hazCurveNum = 1000;
		double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);
		EvenlyDiscretizedFunc curveLnXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		
		// make the site object and set values
		Location siteLoc = new Location(34.+distance/111.,-118.,0.);

		Site site = new Site(siteLoc);
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
			infoString += param.getName()+"\t"+param.getValue()+"\n";
		}
		
		// set the IMT
		imr.setIntensityMeasure(PGA_Param.NAME);
		
		// set site and eqkRup
		imr.setSite(site);
		imr.setEqkRupture(eqkRup);
		
		for (Parameter<?> param : imr.getPropagationEffectParams()) {
			infoString += param.getName()+"\t"+param.getValue()+"\n";
		}
		for (Parameter<?> param : imr.getEqkRuptureParams()) {
			infoString += param.getName()+"\t"+param.getValue()+"\n";
		}
		for (Parameter<?> param : imr.getOtherParams()) {
			infoString += param.getName()+"\t"+param.getValue()+"\n";
		}
		infoString += "\n";
		
		double meanInLnSpace = imr.getMean();
		double stdev = imr.getStdDev();
		double lnPlusOneSigmaPGA = meanInLnSpace+stdev;
		imr.getExceedProbabilities(curveLnXvalues);
		curveLnXvalues.scale(rate);
		
		double twoIn50_prob = 4e-4;
		double pluOneSigmaProb = 0.16;
		double lnTwoIn50_PGA = curveLnXvalues.getFirstInterpolatedX(twoIn50_prob);
		double twoIn50_epsilon = (lnTwoIn50_PGA-meanInLnSpace)/stdev;
		infoString += "mag="+mag+"\nrate="+rate+
				"\nmeanInLnSpace="+(float)meanInLnSpace+
				"\nstdevInLnSpace="+(float)stdev+
				"\nmedianPGA="+(float)Math.exp(meanInLnSpace)+
				"\n84thPercPGA="+(float)Math.exp(meanInLnSpace+stdev)+
				"\nln(twoIn50_PGA)="+(float)lnTwoIn50_PGA+
				"\ntwoIn50_PGA="+(float)Math.exp(lnTwoIn50_PGA)+
				"\ntwoIn50_epsilon="+(float)twoIn50_epsilon+"\n";
		
		if(rate == 1) {
			double lnPlusOneSigmaPGA_test = curveLnXvalues.getFirstInterpolatedX(pluOneSigmaProb);
			double test_epsilon = (lnPlusOneSigmaPGA_test-meanInLnSpace)/stdev;
			infoString += "ln(pluOneSigmaPGA)="+(float)lnPlusOneSigmaPGA+
					"\nln(pluOneSigmaPGA_test)="+(float)lnPlusOneSigmaPGA_test+
					"\ntest_epsilon="+(float)test_epsilon;
		}
		infoString +="\n";

		
		ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
		for(int i=0;i<curveLnXvalues.size();i++) {
			curveLinearXvalues.set(Math.exp(curveLnXvalues.getX(i)),curveLnXvalues.getY(i));
		}
		curveLinearXvalues.setInfo(infoString);
//		System.out.println(infoString);
		double[] returnDoubleArray = {rate,twoIn50_epsilon,mag,meanInLnSpace,stdev,lnTwoIn50_PGA};
		System.out.println((float)rate+"\t"+(float)twoIn50_epsilon+"\t"+mag+"\t"+(float)meanInLnSpace+"\t"+(float)stdev+"\t"+(float)lnTwoIn50_PGA);
		
		double xAxisMin = 1e-2;
		double xAxisMax = 20;
		double yAxisMin = 1e-8;
		double yAxisMax = 2;
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		plottingFuncsArray.add(curveLinearXvalues);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
//		ArbitrarilyDiscretizedFunc pointsOfInteret = new ArbitrarilyDiscretizedFunc();
//		pointsOfInteret.set(Math.exp(meanInLnSpace),0.5);
//		pointsOfInteret.set(Math.exp(lnPlusOneSigmaPGA),pluOneSigmaProb);
//		pointsOfInteret.set(Math.exp(lnTwoIn50_PGA),twoIn50_prob);
//		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 12f, Color.RED));
		
		if(rate==1) {
			DefaultXY_DataSet medianHorzLine = new DefaultXY_DataSet();
			medianHorzLine.set(xAxisMin,0.5);
			medianHorzLine.set(xAxisMax,0.5);	
			plottingFuncsArray.add(medianHorzLine);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		}
		DefaultXY_DataSet medianVertLine = new DefaultXY_DataSet();
		medianVertLine.set(Math.exp(meanInLnSpace),yAxisMin);
		medianVertLine.set(Math.exp(meanInLnSpace),yAxisMax);
		plottingFuncsArray.add(medianVertLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));

		if(rate==1) {
			DefaultXY_DataSet plusOneSigmaHorzLine = new DefaultXY_DataSet();
			plusOneSigmaHorzLine.set(xAxisMin,pluOneSigmaProb);
			plusOneSigmaHorzLine.set(xAxisMax,pluOneSigmaProb);
			plottingFuncsArray.add(plusOneSigmaHorzLine);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, new Color(155,0,155)));
		}
		DefaultXY_DataSet plusOneSigmaVertLine = new DefaultXY_DataSet();
		plusOneSigmaVertLine.set(Math.exp(lnPlusOneSigmaPGA),yAxisMin);
		plusOneSigmaVertLine.set(Math.exp(lnPlusOneSigmaPGA),yAxisMax);
		plottingFuncsArray.add(plusOneSigmaVertLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, new Color(155,0,155)));

		DefaultXY_DataSet twoIn50HorzLine = new DefaultXY_DataSet();
		twoIn50HorzLine.set(xAxisMin,twoIn50_prob);
		twoIn50HorzLine.set(xAxisMax,twoIn50_prob);
		DefaultXY_DataSet twoIn50VertLine = new DefaultXY_DataSet();
		twoIn50VertLine.set(Math.exp(lnTwoIn50_PGA),yAxisMin);
		twoIn50VertLine.set(Math.exp(lnTwoIn50_PGA),yAxisMax);
		plottingFuncsArray.add(twoIn50HorzLine);
		plottingFuncsArray.add(twoIn50VertLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));

		String plotTitle = "Mag = "+mag+" & Event Rate = "+(float) rate;
		String fileNamePrefix = null;
		if(dirName != null) {
			String temp = dirName+"/hazardCurve_Mag"+mag+"_Rate"+rate;
//			System.out.println(temp);
			fileNamePrefix = temp.replaceAll("\\.", "pt");
//			System.out.println(fileNamePrefix);
		}
		String xAxisLabel = "PGA (g)";
		String yAxisLabel = "Probability of Exceedance (per year)";
		Range xAxisRange = new Range(xAxisMin,xAxisMax);
		Range yAxisRange = new Range(yAxisMin,yAxisMax);
		boolean logX = true;
		boolean logY = true;

		writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		return returnDoubleArray;

	}
	
	
	
public static void epsilon_test_mfd(String dirName, boolean popupWindow) {
		
		double distance = 10;
				
		Location loc = new Location(34.,-118.,0.);
		double rake = 0;
		double prob = 1;	// this is ingnored
		PointSurface pointSurface = new PointSurface(loc);
		pointSurface.setAveDip(90);
		double bValue=0.0;
		double totCumRate = 0.01;
		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(bValue, totCumRate, 6.0, 8.0, 11);
		double[] meanArray = new double[mfd.size()];
		double[] stdevArray = new double[mfd.size()];
		
		ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
		imr.setParamDefaults();
		
		double hazCurveLnMin = Math.log(0.001);
		double hazCurveLnMax = Math.log(20);
		int hazCurveNum = 100;
		double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);
		EvenlyDiscretizedFunc totalCurveLnXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		ArrayList<EvenlyDiscretizedFunc> curvesLnXvaluesList = new ArrayList<EvenlyDiscretizedFunc>();
		ArbitrarilyDiscretizedFunc totalCurveLinearXvalues = new ArbitrarilyDiscretizedFunc();
		for(int i=0;i<totalCurveLnXvalues.size();i++) {
			totalCurveLinearXvalues.set(Math.exp(totalCurveLnXvalues.getX(i)),0.0);
		}
		ArrayList<ArbitrarilyDiscretizedFunc> allCurves = new ArrayList<ArbitrarilyDiscretizedFunc>();
		totalCurveLinearXvalues.setName("totalCurveLinearXvalues");
		allCurves.add(totalCurveLinearXvalues);
		
		// make the site object and set values
		Location siteLoc = new Location(34.+distance/111.,-118.,0.);

		Site site = new Site(siteLoc);
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
		}
		
		// set the IMT
		imr.setIntensityMeasure(PGA_Param.NAME);
		
		// set site and eqkRup
		imr.setSite(site);
				
		double twoIn50_prob = 4e-4;
		double pluOneSigmaProb = 0.16;
		
		for(int m=0;m<mfd.size();m++) {
			double mag = mfd.getX(m);
			double rate = mfd.getY(m);
			ProbEqkRupture eqkRup = new ProbEqkRupture(mag, rake, prob, pointSurface, loc);
			imr.setEqkRupture(eqkRup);
			double meanInLnSpace = imr.getMean();
			double stdev = imr.getStdDev();
			double lnPlusOneSigmaPGA = meanInLnSpace+stdev;
			
			meanArray[m] = meanInLnSpace;
			stdevArray[m] = stdev;
			
			EvenlyDiscretizedFunc curveLnXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
			imr.getExceedProbabilities(curveLnXvalues);
			curveLnXvalues.scale(rate);
			curvesLnXvaluesList.add(curveLnXvalues);
			
			double lnTwoIn50_PGA = curveLnXvalues.getFirstInterpolatedX(twoIn50_prob);
			double twoIn50_epsilon = (lnTwoIn50_PGA-meanInLnSpace)/stdev;
			String infoString = "mag="+mag+"\nrate="+rate+
					"\nmeanInLnSpace="+(float)meanInLnSpace+
					"\nstdevInLnSpace="+(float)stdev+
					"\nmedianPGA="+(float)Math.exp(meanInLnSpace)+
					"\n84thPercPGA="+(float)Math.exp(meanInLnSpace+stdev)+
					"\nln(twoIn50_PGA)="+(float)lnTwoIn50_PGA+
					"\ntwoIn50_PGA="+(float)Math.exp(lnTwoIn50_PGA)+
					"\ntwoIn50_epsilon="+(float)twoIn50_epsilon+"\n";

			
			ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
			for(int i=0;i<curveLnXvalues.size();i++) {
				curveLinearXvalues.set(Math.exp(curveLnXvalues.getX(i)),curveLnXvalues.getY(i));
				totalCurveLinearXvalues.set(i,totalCurveLinearXvalues.getY(i)+curveLnXvalues.getY(i));
				totalCurveLnXvalues.set(i,totalCurveLnXvalues.getY(i)+curveLnXvalues.getY(i));
			}
			curveLinearXvalues.setInfo(infoString);
			
			allCurves.add(curveLinearXvalues);

		}
		
		double lnTwoIn50_PGA = totalCurveLnXvalues.getFirstInterpolatedX(twoIn50_prob);
		
		String infoString = "lnTwoIn50_PGA="+lnTwoIn50_PGA+"\nTwoIn50_PGA="+Math.exp(lnTwoIn50_PGA)+"\n";
		
		double ave_mag = 0;
		double ave_epsilon = 0;
		double totWt = 0;
		double maxFract =0;
		
		infoString += "mag\tfract\tepsilon\tmean\tstdev\n";
		for(int m=0;m<mfd.size();m++) {
			double epsilon = (lnTwoIn50_PGA-meanArray[m])/stdevArray[m];
			double fract = curvesLnXvaluesList.get(m).getInterpolatedY(lnTwoIn50_PGA)/twoIn50_prob;
			infoString += (float)mfd.getX(m)+"\t"+(float)fract+"\t"+(float)epsilon+"\t"+
					(float)meanArray[m]+"\t"+(float)stdevArray[m]+"\n";
			ave_mag += fract*mfd.getX(m);
			ave_epsilon+= fract*epsilon;
			totWt+=fract;
			if(fract>maxFract)
				maxFract = fract;
			// mag, epsilon, fract, mag
		}
		infoString +="\nave_mag = "+ (float)ave_mag+"\nave_epsilon = "+(float)ave_epsilon+"\ntotFract (test) = "+(float)totWt;

//		//test
//		double altRate=0;
//		for(int m=0;m<mfd.size();m++) {
//			double fract = curvesLnXvaluesList.get(m).getInterpolatedY(lnTwoIn50_PGA)/twoIn50_prob;
//			altRate += (fract/maxFract)*curvesLnXvaluesList.get(m).getY(0);
//			System.out.println(fract+"\t"+maxFract+"\t"+curvesLnXvaluesList.get(m).getY(0));
//		}
		
		// get epsilon for ave_mag
		ProbEqkRupture eqkRup = new ProbEqkRupture(ave_mag, rake, prob, pointSurface, loc);
		imr.setEqkRupture(eqkRup);
		double meanInLnSpace = imr.getMean();
		double stdev = imr.getStdDev();
		double epsilon = (lnTwoIn50_PGA-meanInLnSpace)/stdev;
	
		infoString += "\nepsilon for ave_mag: "+(float)epsilon;
		
		double meanPlusOneSigmaPGA = Math.exp(meanInLnSpace+stdev);
		infoString += "\nmeanPlusOneSigmaPGA for ave_mag: "+(float)meanPlusOneSigmaPGA+"\n";

		totalCurveLinearXvalues.setInfo(infoString);

		// make curve for ave_mag
		EvenlyDiscretizedFunc curveLnXvaluesForMeanMag = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		imr.getExceedProbabilities(curveLnXvaluesForMeanMag);
		double scaleFactor = twoIn50_prob/curveLnXvaluesForMeanMag.getInterpolatedY(lnTwoIn50_PGA);
		curveLnXvaluesForMeanMag.scale(scaleFactor);
		ArbitrarilyDiscretizedFunc curveLinearXvaluesForMeanMag = new ArbitrarilyDiscretizedFunc();
		for(int i=0;i<curveLnXvaluesForMeanMag.size();i++) {
			curveLinearXvaluesForMeanMag.set(Math.exp(curveLnXvaluesForMeanMag.getX(i)),curveLnXvaluesForMeanMag.getY(i));
		}
		curveLinearXvaluesForMeanMag.setName("Scaled curve for ave_mag");
		double rate = curveLinearXvaluesForMeanMag.getY(0);
		double relativeRate = rate/totalCurveLinearXvalues.getY(0);
		curveLinearXvaluesForMeanMag.setInfo("ave_mag = "+(float)ave_mag+"\nrate = "+(float)rate+"\nrate/tot_curve_rate = "+relativeRate);
		allCurves.add(curveLinearXvaluesForMeanMag);
		
		// test
		double scaleRate = (1./400.)/rate;
		curveLnXvaluesForMeanMag.scale(scaleRate);
		double testLnPGA = curveLnXvaluesForMeanMag.getFirstInterpolatedX(twoIn50_prob);
		double testEpsilon = (testLnPGA-meanInLnSpace)/stdev;
		System.out.println("testEpsilon = "+(float)testEpsilon+" (should be 1.0)");

		
		double xAxisMin = 1e-2;
		double xAxisMax = 20;
		double yAxisMin = 1e-8;
		double yAxisMax = 2;
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		plottingFuncsArray.addAll(allCurves);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		for(int m=0;m<mfd.size();m++) {
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		}
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.ORANGE));

		
		
		DefaultXY_DataSet twoIn50HorzLine = new DefaultXY_DataSet();
		twoIn50HorzLine.set(xAxisMin,twoIn50_prob);
		twoIn50HorzLine.set(xAxisMax,twoIn50_prob);
		DefaultXY_DataSet twoIn50VertLine = new DefaultXY_DataSet();
		twoIn50VertLine.set(Math.exp(lnTwoIn50_PGA),yAxisMin);
		twoIn50VertLine.set(Math.exp(lnTwoIn50_PGA),yAxisMax);
		plottingFuncsArray.add(twoIn50HorzLine);
		plottingFuncsArray.add(twoIn50VertLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));

		String plotTitle = "MFD; b-value="+bValue;
		String fileNamePrefix = null;
		if(dirName != null) {
			String temp = dirName+"/hazardCurve_MFD";
		}
		String xAxisLabel = "PGA (g)";
		String yAxisLabel = "Probability of Exceedance (per year)";
		Range xAxisRange = new Range(xAxisMin,xAxisMax);
		Range yAxisRange = new Range(yAxisMin,yAxisMax);
		boolean logX = true;
		boolean logY = true;

		writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
	}




public static void epsilon_test_oneVsManySrc(String dirName, boolean popupWindow) {
	
	double distance = 10;
			
	Location loc = new Location(34.,-118.,0.);
	double rake = 0;
	double prob = 1;	// this is ingnored
	PointSurface pointSurface = new PointSurface(loc);
	pointSurface.setAveDip(90);
	double bValue=0.0;
	double totCumRate = 0.01;
	int numSrc = 20;
	double srcRate = totCumRate/(double)numSrc;
	double mag = 7;
	double[] meanArray = new double[numSrc];
	double[] stdevArray = new double[numSrc];
	
	ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
	imr.setParamDefaults();
	
	double hazCurveLnMin = Math.log(0.001);
	double hazCurveLnMax = Math.log(20);
	int hazCurveNum = 100;
	double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);
	EvenlyDiscretizedFunc totalCurveLnXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
	ArrayList<EvenlyDiscretizedFunc> curvesLnXvaluesList = new ArrayList<EvenlyDiscretizedFunc>();
	ArbitrarilyDiscretizedFunc totalCurveLinearXvalues = new ArbitrarilyDiscretizedFunc();
	for(int i=0;i<totalCurveLnXvalues.size();i++) {
		totalCurveLinearXvalues.set(Math.exp(totalCurveLnXvalues.getX(i)),0.0);
	}
	ArrayList<ArbitrarilyDiscretizedFunc> allCurves = new ArrayList<ArbitrarilyDiscretizedFunc>();
	totalCurveLinearXvalues.setName("totalCurveLinearXvalues");
	allCurves.add(totalCurveLinearXvalues);
	
	// make the site object and set values
	Location siteLoc = new Location(34.+distance/111.,-118.,0.);

	Site site = new Site(siteLoc);
	for (Parameter<?> param : imr.getSiteParams()) {
		site.addParameter(param);
	}
	
	// set the IMT
	imr.setIntensityMeasure(PGA_Param.NAME);
	
	// set site and eqkRup
	imr.setSite(site);
			
	double twoIn50_prob = 4e-4;
	double pluOneSigmaProb = 0.16;
	
	for(int m=0;m<numSrc;m++) {
		double rate = srcRate;
		ProbEqkRupture eqkRup = new ProbEqkRupture(mag, rake, prob, pointSurface, loc);
		imr.setEqkRupture(eqkRup);
		double meanInLnSpace = imr.getMean();
		double stdev = imr.getStdDev();
		double lnPlusOneSigmaPGA = meanInLnSpace+stdev;
		
		meanArray[m] = meanInLnSpace;
		stdevArray[m] = stdev;
		
		EvenlyDiscretizedFunc curveLnXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		imr.getExceedProbabilities(curveLnXvalues);
		curveLnXvalues.scale(rate);
		curvesLnXvaluesList.add(curveLnXvalues);
		
		double lnTwoIn50_PGA = curveLnXvalues.getFirstInterpolatedX(twoIn50_prob);
		double twoIn50_epsilon = (lnTwoIn50_PGA-meanInLnSpace)/stdev;
		String infoString = "mag="+mag+"\nrate="+rate+
				"\nmeanInLnSpace="+(float)meanInLnSpace+
				"\nstdevInLnSpace="+(float)stdev+
				"\nmedianPGA="+(float)Math.exp(meanInLnSpace)+
				"\n84thPercPGA="+(float)Math.exp(meanInLnSpace+stdev)+
				"\nln(twoIn50_PGA)="+(float)lnTwoIn50_PGA+
				"\ntwoIn50_PGA="+(float)Math.exp(lnTwoIn50_PGA)+
				"\ntwoIn50_epsilon="+(float)twoIn50_epsilon+"\n";

		
		ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
		for(int i=0;i<curveLnXvalues.size();i++) {
			curveLinearXvalues.set(Math.exp(curveLnXvalues.getX(i)),curveLnXvalues.getY(i));
			totalCurveLinearXvalues.set(i,totalCurveLinearXvalues.getY(i)+curveLnXvalues.getY(i));
			totalCurveLnXvalues.set(i,totalCurveLnXvalues.getY(i)+curveLnXvalues.getY(i));
		}
		curveLinearXvalues.setInfo(infoString);
		
		allCurves.add(curveLinearXvalues);

	}
	
	double lnTwoIn50_PGA = totalCurveLnXvalues.getFirstInterpolatedX(twoIn50_prob);
	
	String infoString = "lnTwoIn50_PGA="+lnTwoIn50_PGA+"\nTwoIn50_PGA="+Math.exp(lnTwoIn50_PGA)+"\n";
	
	double ave_mag = 0;
	double ave_epsilon = 0;
	double totWt = 0;
	double maxFract =0;
	
	infoString += "mag\tfract\tepsilon\tmean\tstdev\n";
	for(int m=0;m<numSrc;m++) {
		double epsilon = (lnTwoIn50_PGA-meanArray[m])/stdevArray[m];
		double fract = curvesLnXvaluesList.get(m).getInterpolatedY(lnTwoIn50_PGA)/twoIn50_prob;
		infoString += mag+"\t"+(float)fract+"\t"+(float)epsilon+"\t"+
				(float)meanArray[m]+"\t"+(float)stdevArray[m]+"\n";
		ave_mag += fract*mag;
		ave_epsilon+= fract*epsilon;
		totWt+=fract;
		if(fract>maxFract)
			maxFract = fract;
		// mag, epsilon, fract, mag
	}
	infoString +="\nave_mag = "+ (float)ave_mag+"\nave_epsilon = "+(float)ave_epsilon+"\ntotFract (test) = "+(float)totWt;

//	//test
//	double altRate=0;
//	for(int m=0;m<mfd.size();m++) {
//		double fract = curvesLnXvaluesList.get(m).getInterpolatedY(lnTwoIn50_PGA)/twoIn50_prob;
//		altRate += (fract/maxFract)*curvesLnXvaluesList.get(m).getY(0);
//		System.out.println(fract+"\t"+maxFract+"\t"+curvesLnXvaluesList.get(m).getY(0));
//	}
	
	// get epsilon for ave_mag
	ProbEqkRupture eqkRup = new ProbEqkRupture(ave_mag, rake, prob, pointSurface, loc);
	imr.setEqkRupture(eqkRup);
	double meanInLnSpace = imr.getMean();
	double stdev = imr.getStdDev();
	double epsilon = (lnTwoIn50_PGA-meanInLnSpace)/stdev;

	infoString += "\nepsilon for ave_mag: "+(float)epsilon;
	
	double meanPlusOneSigmaPGA = Math.exp(meanInLnSpace+stdev);
	infoString += "\nmeanPlusOneSigmaPGA for ave_mag: "+(float)meanPlusOneSigmaPGA+"\n";

	totalCurveLinearXvalues.setInfo(infoString);

	// make curve for ave_mag
	EvenlyDiscretizedFunc curveLnXvaluesForMeanMag = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
	imr.getExceedProbabilities(curveLnXvaluesForMeanMag);
	double scaleFactor = twoIn50_prob/curveLnXvaluesForMeanMag.getInterpolatedY(lnTwoIn50_PGA);
	curveLnXvaluesForMeanMag.scale(scaleFactor);
	ArbitrarilyDiscretizedFunc curveLinearXvaluesForMeanMag = new ArbitrarilyDiscretizedFunc();
	for(int i=0;i<curveLnXvaluesForMeanMag.size();i++) {
		curveLinearXvaluesForMeanMag.set(Math.exp(curveLnXvaluesForMeanMag.getX(i)),curveLnXvaluesForMeanMag.getY(i));
	}
	curveLinearXvaluesForMeanMag.setName("Scaled curve for ave_mag");
	double rate = curveLinearXvaluesForMeanMag.getY(0);
	double relativeRate = rate/totalCurveLinearXvalues.getY(0);
	curveLinearXvaluesForMeanMag.setInfo("ave_mag = "+(float)ave_mag+"\nrate = "+(float)rate+"\nrate/tot_curve_rate = "+relativeRate);
	allCurves.add(curveLinearXvaluesForMeanMag);
	
	// test
	double scaleRate = (1./400.)/rate;
	curveLnXvaluesForMeanMag.scale(scaleRate);
	double testLnPGA = curveLnXvaluesForMeanMag.getFirstInterpolatedX(twoIn50_prob);
	double testEpsilon = (testLnPGA-meanInLnSpace)/stdev;
	System.out.println("testEpsilon = "+(float)testEpsilon+" (should be 1.0)");

	
	double xAxisMin = 1e-2;
	double xAxisMax = 20;
	double yAxisMin = 1e-8;
	double yAxisMax = 2;
	
	ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
	
	plottingFuncsArray.addAll(allCurves);
	plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
	for(int m=0;m<numSrc;m++) {
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
	}
	plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.ORANGE));

	
	
	DefaultXY_DataSet twoIn50HorzLine = new DefaultXY_DataSet();
	twoIn50HorzLine.set(xAxisMin,twoIn50_prob);
	twoIn50HorzLine.set(xAxisMax,twoIn50_prob);
	DefaultXY_DataSet twoIn50VertLine = new DefaultXY_DataSet();
	twoIn50VertLine.set(Math.exp(lnTwoIn50_PGA),yAxisMin);
	twoIn50VertLine.set(Math.exp(lnTwoIn50_PGA),yAxisMax);
	plottingFuncsArray.add(twoIn50HorzLine);
	plottingFuncsArray.add(twoIn50VertLine);
	plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
	plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));

	String plotTitle = "MFD; b-value="+bValue;
	String fileNamePrefix = null;
	if(dirName != null) {
		String temp = dirName+"/hazardCurve_MFD";
	}
	String xAxisLabel = "PGA (g)";
	String yAxisLabel = "Probability of Exceedance (per year)";
	Range xAxisRange = new Range(xAxisMin,xAxisMax);
	Range yAxisRange = new Range(yAxisMin,yAxisMax);
	boolean logX = true;
	boolean logY = true;

	writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
			xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
	
}


	
	public static void epsilonVersusRatePlot(String dirName) {
		double mag = 8;
		double minLogRate = Math.log10(5e-4);
		double maxLogRate = Math.log10(1.0);
		int numPoints = 50;
		EvenlyDiscretizedFunc logRateFunc = new EvenlyDiscretizedFunc(minLogRate, maxLogRate,numPoints);
//		System.out.println(logRateFunc);
		
		ArbitrarilyDiscretizedFunc epsilonVersusRateFunc = new ArbitrarilyDiscretizedFunc();
		
		for(int i=0;i< logRateFunc.size(); i++) {
			double rate = Math.pow(10, logRateFunc.getX(i));
			double[] resultArray = epsilon_test(mag,rate,null,false);
			epsilonVersusRateFunc.set(resultArray[0],resultArray[1]);
		}
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		plottingFuncsArray.add(epsilonVersusRateFunc);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		// the following verifies that epsilon=1 corresponds to a rate of 1/400 (0.0025) and that epsilon=1 is for a rate of 8e-4
//		double rateForEpsilon0 = epsilonVersusRateFunc.getFirstInterpolatedX(0.0);
//		double rateForEpsilon1 = epsilonVersusRateFunc.getFirstInterpolatedX(1.0);
//		System.out.println("rateForEpsilon0="+(float)rateForEpsilon0);
//		System.out.println("rateForEpsilon1="+(float)rateForEpsilon1);
//		epsilon_test(mag,rateForEpsilon0,dirName,false);
//		epsilon_test(mag,8e-4,dirName,false);
//		epsilon_test(mag,rateForEpsilon1,dirName,false);
//		epsilon_test(mag,0.0025,dirName,false);
//		epsilon_test(mag,0.00251,dirName,false);
//		epsilon_test(mag,0.00252,dirName,false);
		
		double xAxisMin = 4e-4;
		double xAxisMax = 1.5;
		double yAxisMin = -1;
		double yAxisMax = 3.5;
		
		DefaultXY_DataSet epsilon1_vertLine = new DefaultXY_DataSet();
		epsilon1_vertLine.set(0.0025, yAxisMin);
		epsilon1_vertLine.set(0.0025, yAxisMax);	
		plottingFuncsArray.add(epsilon1_vertLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		DefaultXY_DataSet epsilon1_horzLine = new DefaultXY_DataSet();
		epsilon1_horzLine.set(xAxisMin,1);
		epsilon1_horzLine.set(xAxisMax,1);	
		plottingFuncsArray.add(epsilon1_horzLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));

		DefaultXY_DataSet epsilon0_vertLine = new DefaultXY_DataSet();
		epsilon0_vertLine.set(8e-4, yAxisMin);
		epsilon0_vertLine.set(8e-4, yAxisMax);	
		plottingFuncsArray.add(epsilon0_vertLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		DefaultXY_DataSet epsilon0_horzLine = new DefaultXY_DataSet();
		epsilon0_horzLine.set(xAxisMin,0);
		epsilon0_horzLine.set(xAxisMax,0);	
		plottingFuncsArray.add(epsilon0_horzLine);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));


		GraphWindow graph = new GraphWindow(plottingFuncsArray, "2% in 50-year Epsilon");
		graph.setX_AxisRange(xAxisMin,xAxisMax);
		graph.setY_AxisRange(yAxisMin,yAxisMax);
		graph.setXLog(true);
//		graph.setYLog(true);
		graph.setPlotChars(plotChars);
		graph.setX_AxisLabel("Event Rate (per year)");
		graph.setY_AxisLabel("Epsilon");
		graph.setTickLabelFontSize(22);
		graph.setAxisLabelFontSize(24);
		
		String fileNamePrefix=null;
		if(dirName != null)
			fileNamePrefix = dirName+"/epsilonVersusRatePlot";
		if(fileNamePrefix != null) {
			try {
				graph.saveAsPDF(fileNamePrefix+".pdf");
				graph.saveAsPNG(fileNamePrefix+".png");
				graph.saveAsTXT(fileNamePrefix+".txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}		

	}
	
	
	/**
	 * The two EllB lines overlap and the two Shaw lines overlap
	 */
	public static void partMFD_ForGR_Fault() {
		
		double mMax = 8.25;
		double mMin = 6.35;
		double ddw = 12e3;
		
		String dirName = null;
		boolean popupWindow = true;

		ArrayList<ScalingRelationships> scalingRelList = new ArrayList<ScalingRelationships>();
		scalingRelList.add(ScalingRelationships.ELLSWORTH_B);
		scalingRelList.add(ScalingRelationships.HANKS_BAKUN_08);
		scalingRelList.add(ScalingRelationships.SHAW_2009_MOD);
		scalingRelList.add(ScalingRelationships.ELLB_SQRT_LENGTH);
		scalingRelList.add(ScalingRelationships.SHAW_CONST_STRESS_DROP);
		
//		for(int i=0;i<scalingRelList.size();i++) {
//			System.out.println(scalingRelList.get(i).getShortName()+"\t"+scalingRelList.get(i).getMag(12*12*1e6, 12*1e3));
//		}
//		System.exit(-1);
		
		GutenbergRichterMagFreqDist grDist = new GutenbergRichterMagFreqDist(mMin, 20, 0.1, mMin, mMax, 1.0,1.0);
		grDist.scaleToCumRate(mMin, 1.0);

		ArrayList<XY_DataSet> mfd_List = new ArrayList<XY_DataSet>();
		grDist.setName("GR Dist");
		mfd_List.add(grDist);
		
		double minY = grDist.getMinY();
		double maxY = grDist.getMaxY();

		IncrementalMagFreqDist meanPartMFD = new IncrementalMagFreqDist(mMin, 20, 0.1);
		meanPartMFD.setName("Mean Part. MFD");

		for(ScalingRelationships scaleRel:scalingRelList) {
			double maxLength = scaleRel.getArea(mMax, ddw)/ddw;	// meters
			IncrementalMagFreqDist partMFD = new IncrementalMagFreqDist(mMin, 20, 0.1);
			for(int i=0;i<partMFD.size();i++) {
				double mag = partMFD.getX(i);
				double area = scaleRel.getArea(mag, ddw);
				double rupLengh = area/ddw;
				double aveSlip = scaleRel.getAveSlip(area, rupLengh, ddw, ddw, Double.NaN);
				if(i==0)
					System.out.println(scaleRel.getName()+" length for M "+mMin+": "+rupLengh/1e3);
				System.out.println((float)mag+"\t"+(float)partMFD.getY(i)+"\t"+(float)aveSlip);
				double partRate = grDist.getY(i)*rupLengh/maxLength;
				partMFD.set(i,partRate);
				meanPartMFD.set(i,meanPartMFD.getY(i)+partRate/scalingRelList.size());
			}
			partMFD.setName(scaleRel.getName()+" Part MFD");
			mfd_List.add(partMFD);
			
			if(minY>partMFD.getMinY())
				minY = partMFD.getMinY();
			if(maxY<partMFD.getMaxY())
				maxY = partMFD.getMaxY();

		}
		System.out.println(grDist.toString());
		
		mfd_List.add(meanPartMFD);

		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.ORANGE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.CYAN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
		
		Range xAxisRange = new Range(6.3,8.3);
//		Range yAxisRange = new Range(minY,maxY);
		Range yAxisRange = new Range(1e-3,0.3);

		writeAndOrPlotFuncs(mfd_List, plotChars, "Part MFDs for GR Fault","Magnitude","Rate",
				xAxisRange,yAxisRange,false,true,dirName, popupWindow);
		
	}
	
	/**
	 * 
	 */
	public static void slipPDFatPoint_ForGR_Fault() {
		
		double mMax = 8.25;
		double mMin = 6.35;
		double ddw = 12e3;
		
		String dirName = null;
		boolean popupWindow = true;

		ArrayList<ScalingRelationships> scalingRelList = new ArrayList<ScalingRelationships>();
		scalingRelList.add(ScalingRelationships.ELLSWORTH_B);
		scalingRelList.add(ScalingRelationships.HANKS_BAKUN_08);
		scalingRelList.add(ScalingRelationships.SHAW_2009_MOD);
		scalingRelList.add(ScalingRelationships.ELLB_SQRT_LENGTH);
		scalingRelList.add(ScalingRelationships.SHAW_CONST_STRESS_DROP);
		
//		for(int i=0;i<scalingRelList.size();i++) {
//			System.out.println(scalingRelList.get(i).getShortName()+"\t"+scalingRelList.get(i).getMag(12*12*1e6, 12*1e3));
//		}
//		System.exit(-1);
		
		GutenbergRichterMagFreqDist grDist = new GutenbergRichterMagFreqDist(mMin, 200, 0.01, mMin, mMax, 1.0,1.0);
		grDist.scaleToCumRate(mMin, 1.0);

		ArrayList<XY_DataSet> slip_pdf_List = new ArrayList<XY_DataSet>();
		grDist.setName("GR Dist");
		
		
		
		ArrayList<DefaultXY_DataSet> slipVsMagList = new ArrayList<DefaultXY_DataSet>();
		
		double minY = grDist.getMinY();
		double maxY = grDist.getMaxY();

		HistogramFunction meanSlipPDF = new HistogramFunction(0.5, 14, 1.0);
		HistogramFunction tempNumHist = new HistogramFunction(0.5, 14, 1.0);
		meanSlipPDF.setName("meanSlipPDF");

		for(ScalingRelationships scaleRel:scalingRelList) {
			double maxLength = scaleRel.getArea(mMax, ddw)/ddw;	// meters
			HistogramFunction pdf = new HistogramFunction(0.5, 14, 1.0);
			DefaultXY_DataSet slipVsMag = new DefaultXY_DataSet();
			for(int i=0;i<grDist.size();i++) {
				double mag = grDist.getX(i);
				double area = scaleRel.getArea(mag, ddw);
				double rupLengh = area/ddw;
				double partRate = grDist.getY(i)*rupLengh/maxLength;
				double aveSlip = scaleRel.getAveSlip(area, rupLengh, ddw, ddw, Double.NaN);
				pdf.add(aveSlip, partRate);
//				System.out.println((float)mag+"\t"+(float)partMFD.getY(i)+"\t"+(float)aveSlip);
				tempNumHist.add(aveSlip, 1.0);
				slipVsMag.set(mag,aveSlip);
			}
			pdf.setName(scaleRel.getName()+" slip PDF");
			pdf.normalizeBySumOfY_Vals();
			slip_pdf_List.add(pdf);
			
			slipVsMagList.add(slipVsMag);
			

		}
		
		for(int i=0;i<meanSlipPDF.size();i++) {
			double num = 0;
			double sum = 0;
			for(XY_DataSet pdf:slip_pdf_List) {
				if(pdf.getY(i) > 0.0) {
					sum +=  pdf.getY(i);
					num += 1;;
				}
			}
			if(num>0)
				meanSlipPDF.set(i,sum/num);
		}
		meanSlipPDF.normalizeBySumOfY_Vals();
		
		slip_pdf_List.add(meanSlipPDF);

		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.ORANGE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.CYAN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
		
		Range xAxisRange = new Range(0,12);
		Range yAxisRange = new Range(1e-5,1.0);
		
		GraphWindow graph = new GraphWindow(slipVsMagList, "slip vs mag", plotChars);

		writeAndOrPlotFuncs(slip_pdf_List, plotChars, "Slip PDF for GR Fault","Slip","Density",
				xAxisRange,yAxisRange,false,true,dirName, popupWindow);
		
	}

	
	public static SummedMagFreqDist getU3_LongTermOnFaultIncrementalMFD(boolean mkPlot) {
		File file = new File("src/scratch/ned/U3_LongTermOnFaultMFD_Incremental.txt");
		
		SummedMagFreqDist mfd  = new SummedMagFreqDist(6.35, 8.35, 21);
		
		boolean first = true;
		try {
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				if(first) {
					first = false;
					continue;
				}
				//System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				double xVal = Double.valueOf(split[0]);
//				System.out.println(Double.valueOf(split[0])+"\t"+Double.valueOf(split[1]));
				if(mfd.hasX(xVal) && xVal>6.3) {
					mfd.add(xVal,Double.valueOf(split[1]));
				}
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		System.out.println(mfd);
		
		if(mkPlot) {
			ArrayList<IncrementalMagFreqDist> mfdList = new ArrayList<IncrementalMagFreqDist>();
			mfdList.add(mfd);
			mfdList.add(modify_bValueOfMFD(mfd,-0.5));
			mfdList.add(modify_bValueOfMFD(mfd,-1.0));
			mfdList.add(modify_bValueOfMFD(mfd,-1.5));
			mfdList.add(modify_bValueOfMFD(mfd,-2.0));
			mfdList.add(modify_bValueOfMFD(mfd,+0.5));
			mfdList.add(modify_bValueOfMFD(mfd,+1.0));
			mfdList.add(modify_bValueOfMFD(mfd,+1.5));
			mfdList.add(modify_bValueOfMFD(mfd,+2.0));
			
			GraphWindow graph = new GraphWindow(mfdList, "MFD"); 
			graph.setX_AxisLabel("Magnitude");
			graph.setY_AxisLabel("Incremental Rate (per yr)");
			graph.setYLog(true);
		}

		return mfd;
	}
	
	
	public static SummedMagFreqDist getU3_SanJacintoNuclIncrMFD() {
		File file = new File("src/scratch/ned/U3_San_Jacinto_Stepovers_Combined__nucleation.txt");
		
		SummedMagFreqDist mfd  = new SummedMagFreqDist(6.05, 8.55, 26);
		
		boolean first = true;
		try {
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				if(first) {
					first = false;
					continue;
				}
				//System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				double xVal = Double.valueOf(split[0]);
//				System.out.println(Double.valueOf(split[0])+"\t"+Double.valueOf(split[1]));
				if(mfd.hasX(xVal) && xVal>6.3) {
					mfd.add(xVal,Double.valueOf(split[1]));
				}
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		System.out.println(mfd);
		
		return mfd;
	}
	
	public static SummedMagFreqDist getU2_SanJacintoNuclIncrMFD() {
		File file = new File("src/scratch/ned/U2_San_Jacinto_Stepovers_Combined__nucleation.txt");
		
		SummedMagFreqDist mfd  = new SummedMagFreqDist(6.05, 8.55, 26);
		
		boolean first = true;
		try {
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				if(first) {
					first = false;
					continue;
				}
				//System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				double xVal = Double.valueOf(split[0]);
//				System.out.println(Double.valueOf(split[0])+"\t"+Double.valueOf(split[1]));
				if(mfd.hasX(xVal) && xVal>6.3) {
					mfd.add(xVal,Double.valueOf(split[1]));
				}
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		System.out.println(mfd);
		
		return mfd;
	}


	
	public static SummedMagFreqDist modify_bValueOfMFD(IncrementalMagFreqDist mfd, double deltaBvalue) {
		SummedMagFreqDist modMFD = new SummedMagFreqDist(mfd.getMinX(), mfd.getMaxX(), mfd.size());
		double moRate = mfd.getTotalMomentRate();
		GutenbergRichterMagFreqDist grDist = new GutenbergRichterMagFreqDist(deltaBvalue, 1.0,
				mfd.getMinX(), mfd.getMaxX(), mfd.size());
		for(int i=0;i<mfd.size();i++)
			modMFD.add(i, mfd.getY(i)*grDist.getY(i));
		modMFD.scaleToTotalMomentRate(moRate);
		
		return modMFD;
	}
	
	
	/**
	 * This computes the COV of 2in50 values for the LA site from Peter's U3-TI analyses.
	 */
	public static void computeLA_Haz_COV() {
		// Read section rate constraints
		ArrayList<Double> wtList = new  ArrayList<Double>();
		ArrayList<Double> twoIn50Val = new  ArrayList<Double>();
		ArbDiscrEmpiricalDistFunc distFunc = new ArbDiscrEmpiricalDistFunc();
		HistogramFunction histDist = new HistogramFunction(0.0,1.95, 200);
//		HistogramFunction histDist = new HistogramFunction(0.45,0.95, 50);

//		File file = new File("/Users/field/Desktop/U3_HazDataFromPeter/LOS_ANGELES/summary.txt");
		File file = new File("/Users/field/Desktop/U3_HazDataFromPeter/LOS_ANGELES+GMM+EPI/distros_PE2IN50.txt");
		
		boolean first = true;
		try {
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				if(first) {
					first = false;
					continue;
				}
				//System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				wtList.add(Double.valueOf(split[0]));
				twoIn50Val.add(Double.valueOf(split[1]));
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println(wtList.size());
		
		for(int i=0;i<wtList.size();i++) {
			distFunc.set(twoIn50Val.get(i), wtList.get(i));
			double yVal = histDist.getY(twoIn50Val.get(i));
			histDist.set(twoIn50Val.get(i),wtList.get(i)+yVal);
		}
		
		System.out.println("SumTest="+distFunc.getSumOfAllY_Values());
		System.out.println("Mean="+distFunc.getMean());
		System.out.println("StdDev="+distFunc.getStdDev());
		System.out.println("COV="+distFunc.getCOV());
		
		System.out.println("Mean="+histDist.computeMean());
		System.out.println("StdDev="+histDist.computeStdDev());
		System.out.println("COV="+histDist.computeCOV());
//		System.out.println(histDist);
		
		GraphWindow graph = new GraphWindow(histDist.getCumulativeDistFunctionWithHalfBinOffset(), "Test"); 
		graph.setX_AxisLabel("IMT");
		graph.setY_AxisLabel("Cum Wt");
//		graph.setX_AxisRange(0.4, 1200);
//		graph.setY_AxisRange(1e-6, 1);
//		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.RED));
//		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.GREEN));
//		graph.setPlotChars(plotChars);
//		graph.setYLog(true);
//		graph.setXLog(true);


	}
	
	
	public static void writeAveLengthForMag(double mag, double ddwKm) {
		
		ArrayList<ScalingRelationships> scalingRelList = new ArrayList<ScalingRelationships>();
		scalingRelList.add(ScalingRelationships.ELLSWORTH_B);
		scalingRelList.add(ScalingRelationships.HANKS_BAKUN_08);
		scalingRelList.add(ScalingRelationships.SHAW_2009_MOD);
		scalingRelList.add(ScalingRelationships.ELLB_SQRT_LENGTH);
		scalingRelList.add(ScalingRelationships.SHAW_CONST_STRESS_DROP);
		
		System.out.println("M="+mag+"; ddw="+ddwKm+":");
		double ddwMeters = ddwKm*1e3;
		double aveLength=0;
		for(ScalingRelationships scRel:scalingRelList) {
			double lengthKm = (scRel.getArea(mag, ddwMeters)/ddwMeters)*1e-3;
			System.out.println("\t"+(float)lengthKm+" for "+scRel.getName());
			aveLength+=lengthKm/scalingRelList.size();		
		}
		System.out.println("aveLength="+(float)aveLength);

	}

	public static void main(String[] args) {
		
		doit3();
		
		
//		slipPDFatPoint_ForGR_Fault();
		
//		writeAveLengthForMag(8.0, 11.0);

		//doit2();
		
//		computeLA_Haz_COV();
		
//		partMFD_ForGR_Fault();
		
//		getU3_LongTermOnFaultIncrementalMFD(false);
		
//		doIt1();
		
//		// This is for the epsilon analysis
//
//		String dirName = "src/scratch/ned/epsilon_analysis";
//	    File file = new File(dirName);
//	    file.mkdirs();
//	    
////	    epsilon_test_oneVsManySrc(dirName, true);
////	    epsilon_test_mfd(dirName, true);
//
////		double mag = 5;
//////		double rate = 1./(400.*3.125);
////		double rate = 1e-2;
////		epsilon_test(5,1.0,dirName,true);
////		epsilon_test(8,1.0,dirName,true);
////		epsilon_test(5,1e-2,dirName,true);
////		epsilon_test(8,1e-2,dirName,true);
//		
//		epsilonVersusRatePlot(dirName);
		
	}

}
