package scratch.ned;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.rupForecastImpl.PoissonFaultERF;
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
			graph.setTickLabelFontSize(18);
			graph.setAxisLabelFontSize(20);
			
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


	public static void main(String[] args) {
		
		doIt2();
	}

}
