package scratch.ned;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
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
import org.opensha.sha.earthquake.rupForecastImpl.PoissonFaultERF;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.SingleMagFreqDist;
import org.opensha.sha.param.SimpleFaultParameter;

import scratch.UCERF3.FaultSystemSolution;
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
			
		// ratio between the max and min mag
		double twoIn50valueRatio = twoIn50valueArray[twoIn50valueArray.length-1]/twoIn50valueArray[0];
		double tenIn50valueRatio = tenIn50valueArray[tenIn50valueArray.length-1]/tenIn50valueArray[0];
	
		String returnString = (float)twoIn50valueRatio + "\t"+(float)tenIn50valueRatio;
		
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
			
		// frat diff between the max and min mag
		double twoIn50valueMinMag = (twoIn50valueArray[0]/grTwoIn50value);
		double tenIn50valueMinMag = (tenIn50valueArray[0]/grTenIn50value);
		double twoIn50valueMaxMag = (twoIn50valueArray[twoIn50valueArray.length-1]/grTwoIn50value);
		double tenIn50valueMaxMag = (tenIn50valueArray[tenIn50valueArray.length-1]/grTenIn50value);
	
		String returnString = (float)twoIn50valueMinMag + "\t"+(float)twoIn50valueMaxMag + "\t"+(float)tenIn50valueMinMag + "\t"+(float)tenIn50valueMaxMag;
		
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
	
	public static void doIt1() {
		
		double[] siteDistKmArray = {0, 10};
		double[] saPeriod = {0.0, 1.0};
//		double[] magArray = {6.3, 8.0};
		double[] magArray = {6.3, 6.7, 7.1, 7.6, 8.0};

		
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
			for(int p=0;p<saPeriod.length;p++) {
				for(int d=0;d<siteDistKmArray.length;d++) {
					Location location = new Location(36+siteDistKmArray[d]/111, -120,0);	// 
					System.out.println(scalingRelList.get(i).getName());
					double[] tempRate = new double[magArray.length];
					ERF[] erfArray = new ERF[magArray.length];
					for(int m=0;m<magArray.length;m++) {
						double mag=magArray[m];
						double area = scalingRelList.get(i).getArea(mag, fltWidth);
						double length = area/fltWidth;
						double slip = scalingRelList.get(i).getAveSlip(area, length, fltWidth);
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
						SimpleFaultParameter faultParam = (SimpleFaultParameter)erf.getParameter(PoissonFaultERF.FAULT_PARAM_NAME);
						faultParam.setAll(gridSpacing, lats, lons, dips, depths, SimpleFaultParameter.STIRLING);
						faultParam.setEvenlyGriddedSurfaceFromParams();
						erf.getTimeSpan().setDuration(duration);
						erf.updateForecast();
						erfArray[m] = erf;
					}
					String title = scalingRelList.get(i).getName();
					String ratioString = compareHazardCurves1(erfArray, location, saPeriod[p], dirName, popupWindow, title, duration);

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

	
	
	public static void doIt2() {
		
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
			magArray[0] = scaleRel.getMag(12*12*1e6, 12*1e3);
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
						double slip = scaleRel.getAveSlip(area, length, fltWidth);
						testSlipRate+=slip*grDist.getY(m);
					}
					grDist.scale(slipRate/testSlipRate);
					// test slip rates:
					testSlipRate =0;
					for(int m=0;m<grDist.size();m++) {
						double mag=grDist.getX(m);
						double area = scaleRel.getArea(mag, fltWidth);
						double length = area/fltWidth;
						double slip = scaleRel.getAveSlip(area, length, fltWidth);
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
						double slip = scaleRel.getAveSlip(area, length, fltWidth);
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
			double slip = scalingRel.getAveSlip(area, length, fltWidth);
			testSlipRate+=slip*grDist.getY(m);
		}
		grDist.scale(slipRate/testSlipRate);
		
		// test slip rates:
		testSlipRate =0;
		for(int m=0;m<grDist.size();m++) {
			double mag=grDist.getX(m);
			double area = scalingRel.getArea(mag, fltWidth);
			double length = area/fltWidth;
			double slip = scalingRel.getAveSlip(area, length, fltWidth);
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

	public static void main(String[] args) {
		
//		doIt2();
		
		// This is for the epsilon analysis

		String dirName = "src/scratch/ned/epsilon_analysis";
	    File file = new File(dirName);
	    file.mkdirs();
	    
//	    epsilon_test_oneVsManySrc(dirName, true);
//	    epsilon_test_mfd(dirName, true);

//		double mag = 5;
////		double rate = 1./(400.*3.125);
//		double rate = 1e-2;
//		epsilon_test(5,1.0,dirName,true);
//		epsilon_test(8,1.0,dirName,true);
//		epsilon_test(5,1e-2,dirName,true);
//		epsilon_test(8,1e-2,dirName,true);
		
		epsilonVersusRatePlot(dirName);
		
	}

}
