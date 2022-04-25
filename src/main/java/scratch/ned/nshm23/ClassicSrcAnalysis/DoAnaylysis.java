package scratch.ned.nshm23.ClassicSrcAnalysis;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.CoolingScheduleType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.rupForecastImpl.FaultRuptureSource;
import org.opensha.sha.earthquake.rupForecastImpl.FloatingPoissonFaultSource;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GaussianMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SparseGutenbergRichterSolver;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sra.rtgm.RTGM;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import scratch.UCERF3.U3FaultSystemRupSet;
import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.ned.FSS_Inversion2019.FaultSystemRuptureRateInversion;
import scratch.ned.FSS_Inversion2019.PlottingUtils;
import scratch.ned.FSS_Inversion2019.SectionRateConstraint;
import scratch.ned.FSS_Inversion2019.SegmentationConstraint;
import scratch.ned.FSS_Inversion2019.SlipRateSegmentationConstraint;
import scratch.ned.FSS_Inversion2019.SimpleFaultInversion.InversionSolutionType;
import scratch.ned.FSS_Inversion2019.SimpleFaultInversion.MFD_TargetType;
import scratch.ned.FSS_Inversion2019.SimpleFaultInversion.SlipRateProfileType;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.ScalingRelationshipEnum;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.SlipAlongRuptureModelEnum;

public class DoAnaylysis {
		
	final static boolean D = true;
	
	public final static String ROOT_PATH = "src/main/java/scratch/ned/nshm23/ClassicSrcAnalysis/Results/";

	
	final static double hazGridSpacing = 0.1;
	
	final static double[] hazardProbArray = {0.02, 0.10};
	final static String[] hazardProbNameArray = {"2in50", "10in50"};
	final static double hazardDurationYrs = 50;

	final static double hazCurveLnMin = Math.log(0.001);
	final static double hazCurveLnMax = Math.log(10);
	final static int hazCurveNum = 40;
	final static double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);

	/**
	 * 
	 */
	private static GriddedGeoDataSet readHazardMapDataFromFile(String fileName, GriddedRegion region) {
		GriddedGeoDataSet griddedDataSet = new GriddedGeoDataSet(region, true);
		try {
			File file = new File(fileName);
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			Preconditions.checkState(fileLines.size()-1 == griddedDataSet.size(), "File and gridded region have different numbers of points (%s vs %s, respectively)", fileLines.size(), griddedDataSet.size());
			for(int i=0;i<griddedDataSet.size();i++) {
				String str = fileLines.get(i+1);	 // skip header
				String[] split = str.split("\t");
				int index = Integer.parseInt(split[0]);
				Preconditions.checkState(index == i, "Incompatible indices (%s vs %s, respectively)", index, i);
				griddedDataSet.set(i, Double.parseDouble(split[1]));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return griddedDataSet;	
	}

	
	/**
	 * 
	 * @param fileName1 - numerator data
	 * @param fileName2 - denominator data
	 * @param label - label for the plot and file name prefix
	 * @param dirName - directory (full path); this cannot be null
	 * @param popupWindow - whether to show results in a pop up window
	 */
	public static void makeHazardMapRatio(String fileName1, String fileName2, String label, String dirName, 
			boolean popupWindow, GriddedRegion region, FaultTrace faultTrace) {
		
		GriddedGeoDataSet data1 = readHazardMapDataFromFile(fileName1, region);
		GriddedGeoDataSet data2 = readHazardMapDataFromFile(fileName2, region);
		HistogramFunction ratioHistogram = new HistogramFunction(0.01,300,0.02);

		
		double[] ratioArray = new double[data1.size()];
		double fractMoreThan10percASway = 0;
		double aveAbsValDiff = 0;
		double min = Double.MAX_VALUE;
		double max = 0;
		double distFrom1 = Double.MAX_VALUE;
		Location maxLoc=null, minLoc=null, equalLoc=null;
		for(int i=0;i<data1.size();i++) {
			double ratio = data1.get(i)/data2.get(i);
//			if(ratio < ratioHistogram.getMaxX())
			ratioHistogram.add(ratio, 1.0);
			// put log10 ratio in first data set:
			data1.set(i, Math.log10(ratio));
			ratioArray[i] = ratio;
			aveAbsValDiff += Math.abs(ratio-1.0)/ratioArray.length;
			if(ratio>1.1 || ratio<0.9) {
				fractMoreThan10percASway += 1;
			}
			if(min>ratio) {
				min=ratio;
				minLoc = data1.getLocation(i);
			}
			if(max<ratio) {
				max=ratio;
				maxLoc = data1.getLocation(i);
			}
			if(Math.abs(ratio-1.0)<distFrom1) {
				distFrom1 = Math.abs(ratio-1);
				equalLoc = data1.getLocation(i);
			}
		}
		fractMoreThan10percASway /= (double)data1.size();
		double mean = StatUtils.mean(ratioArray);
		double cov = Math.sqrt(StatUtils.populationVariance(ratioArray)) / mean;
		aveAbsValDiff /= mean; // normalize so it's fractional diff
				
		try {
//			CPT cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance().rescale(-0.48, 0.48); // factor of 3
			CPT cpt = new CPT();
			Color purple = new Color(128,0,128);
			Color orange = new Color(255,165,0); // this looks better than the default
			cpt.add(new CPTVal((float)Math.log10(1.0/5.0), purple, (float)Math.log10(1.0/3.0), Color.BLUE));
			cpt.add(new CPTVal((float)Math.log10(1.0/3.0), Color.BLUE, (float)Math.log10(1.0/2.0), Color.CYAN));
			cpt.add(new CPTVal((float)Math.log10(1.0/2.0), Color.CYAN, (float)Math.log10(1.0/1.1), Color.GREEN));
			cpt.add(new CPTVal((float)Math.log10(1.0/1.1), Color.GREEN, (float)Math.log10(1/1.05), Color.LIGHT_GRAY));
			cpt.add(new CPTVal((float)Math.log10(1/1.05), Color.LIGHT_GRAY, (float)Math.log10(1.05), Color.LIGHT_GRAY));
//			cpt.add(new CPTVal((float)Math.log10(1.0/1.1), Color.GREEN, (float)Math.log10(1.0), Color.LIGHT_GRAY));
//			cpt.add(new CPTVal((float)Math.log10(1.0), Color.LIGHT_GRAY, (float)Math.log10(1.1), Color.YELLOW));
			cpt.add(new CPTVal((float)Math.log10(1.05), Color.LIGHT_GRAY, (float)Math.log10(1.1), Color.YELLOW));
			cpt.add(new CPTVal((float)Math.log10(1.1), Color.YELLOW, (float)Math.log10(2.0), orange));
			cpt.add(new CPTVal((float)Math.log10(2.0), orange, (float)Math.log10(3.0), Color.RED));
			cpt.add(new CPTVal((float)Math.log10(3.0), Color.RED, (float)Math.log10(5.0), Color.MAGENTA));
			cpt.setBelowMinColor(purple);
			cpt.setAboveMaxColor(Color.MAGENTA);
			
			ArrayList<LocationList> faults = new ArrayList<LocationList>();
			faults.add(faultTrace);
			double[] values = new double[faults.size()];
			for(int i=0;i<values.length;i++)
				values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;

			GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, data1, hazGridSpacing, region, true, "Log10 "+label);

			// override default trace width
			for(PSXYPolygon fltTracePoly: map.getPolys())
				fltTracePoly.setPenWidth(1.1);
			
			// remove coast and political boundaries
			map.setCoast(null);
			
			// I din't see any effect for the following
//			map.setInterpSettings(GMT_InterpolationSettings.getDefaultSettings());

			try {
				FaultBasedMapGen.SAVE_ZIPS=true;
				FaultBasedMapGen.plotMap(new File(dirName), label, popupWindow, map);	
			} catch (GMT_MapException e) {
				e.printStackTrace();
			}

			// write out values to a text file
			FileWriter fw = new FileWriter(dirName+"/"+label+".txt");
			fw.write("index\tvalue\tlatitude\tlongitude\n");
			for(int i=0;i<data1.size(); i++)	 {
				Location loc = data1.getLocation(i);
				fw.write(i+"\t"+data1.get(i)+"\t"+loc.getLatitude()+"\t"+loc.getLongitude()+"\n");
			}
			fw.close();
			
			ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
			ratioHistogram.setName("ratioHistogram");
			ratioHistogram.setInfo("Num Data = "+(float)ratioHistogram.calcSumOfY_Vals()+"\nMean = "+mean+
					"\nFraction More Than 10% Away From 1.0: "+
					fractMoreThan10percASway+"\nCOV = "+(float)cov+
					"\nAveAbsValDiff (mean normalized) = "+(float)aveAbsValDiff+
					"\nminRatio = "+(float)min+" at lat/lon: "+minLoc.getLatitude()+", "+minLoc.getLongitude()+
					"\nmaxRatio = "+(float)max+" at lat/lon: "+maxLoc.getLatitude()+", "+maxLoc.getLongitude()+
					"\nLoc with ratio closest to 1.0: "+equalLoc.getLatitude()+", "+equalLoc.getLongitude());
			ratioHistogram.scale(1.0/ratioHistogram.calcSumOfY_Vals()/ratioHistogram.getDelta());;
			funcs.add(ratioHistogram);
			
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
			
			DefaultXY_DataSet tenPercentDiffLine1 = new DefaultXY_DataSet();
			DefaultXY_DataSet tenPercentDiffLine2 = new DefaultXY_DataSet();
			tenPercentDiffLine1.set(0.9,0.0);
			tenPercentDiffLine1.set(0.9,ratioHistogram.getMaxY());
			tenPercentDiffLine2.set(1.1,0.0);
			tenPercentDiffLine2.set(1.1,ratioHistogram.getMaxY());
			PlotCurveCharacterstics tenPercentDiffLinesPlotChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK);
			funcs.add(tenPercentDiffLine1);
			funcs.add(tenPercentDiffLine2);
			plotChars.add(tenPercentDiffLinesPlotChar);
			plotChars.add(tenPercentDiffLinesPlotChar);

			
			double minX=0.5;
			double maxX=1.5;

			String fileNamePrefix = null;
			if(dirName != null)
				fileNamePrefix = dirName+"/"+label+"_Histogram";
			String plotName =label+"_Histogram";
			String xAxisLabel = "Ratio";
			String yAxisLabel = "PDF";
			Range xAxisRange = new Range(minX,maxX);
//			Range yAxisRange = null;
			Range yAxisRange = new Range(0,10);
			boolean logX = false;
			boolean logY = false;

			PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, 
					xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

			yAxisRange=null;		

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}

	
	private static void writeHazardMapGriddedData(GriddedGeoDataSet hazData, String fileName) {
		FileWriter fw;
		try {
			fw = new FileWriter(fileName);
			fw.write("index\tvalue\tlatitude\tlongitude\n");
			for(int i=0;i<hazData.size(); i++)	 {
				Location loc = hazData.getLocation(i);
				fw.write(i+"\t"+hazData.get(i)+"\t"+loc.getLatitude()+"\t"+loc.getLongitude()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	
	/**
	 * This is makes to hazard maps, one for 10% in 50-year and one for 2% in 50-year ground motions.  
	 * For saPeriod of 1.0 or 0.2, this also make RTGM maps.
	 * @param fltSysRupInversion
	 * @param saPeriod - set as 0 for PGA; this will crash if SA period is not supported
	 * @param dirName
	 * @param popupWindow
	 */
	public static void makeHazardMaps(AbstractERF erf, double saPeriod, String hazDirName, boolean popupWindow,
			GriddedRegion region, ScalarIMR imr, FaultTrace faultTrace) {
		
		if(D) System.out.println("Making hazard map");
		
		File hazDirFile = null;
		if(hazDirName != null) {
			hazDirFile = new File(hazDirName);
			hazDirFile.mkdirs();
		}
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		// the following makes RTGM if sa period is 1.0 or 0.2
		ArrayList<GriddedGeoDataSet> griddedDataList = makeHazardMapGriddedData(erf, saPeriod, true, region, imr);	
		
		GriddedGeoDataSet rtgmGriddedDataSetArray=null;
		if(griddedDataList.size()>2)
			rtgmGriddedDataSetArray = griddedDataList.get(2);
		
		try {
			for(int p=0;p<hazardProbNameArray.length;p++) {
				CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(griddedDataList.get(p).getMinZ(), griddedDataList.get(p).getMaxZ());
//				CPT cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance().rescale(gridData.getMinZ(), gridData.getMaxZ());
				ArrayList<LocationList> faults = new ArrayList<LocationList>();
				faults.add(faultTrace);
				double[] values = new double[faults.size()];
				for(int i=0;i<values.length;i++)
					values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;
				
				String label = imtString+"_"+hazardProbNameArray[p];
								
				GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, griddedDataList.get(p), hazGridSpacing, region, true, label);
//				map.setGenerateKML(true); // this tells it to generate a KML file that can be loaded into Google Earthe
				
				// override default trace width
				for(PSXYPolygon fltTracePoly: map.getPolys())
					fltTracePoly.setPenWidth(1);
				
				// remove coast and poitical boundaries
				map.setCoast(null);
				
				try {
					FaultBasedMapGen.SAVE_ZIPS=true;
					FaultBasedMapGen.plotMap(hazDirFile, label, popupWindow, map);	
				} catch (GMT_MapException e) {
					e.printStackTrace();
				}
				
				// write out values to a text file
				if(hazDirName != null) {
					String fileName = hazDirName +"/"+label+".txt";
					writeHazardMapGriddedData(griddedDataList.get(p), fileName);
				}
			}
			
			if(rtgmGriddedDataSetArray != null) {
				CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(rtgmGriddedDataSetArray.getMinZ(), rtgmGriddedDataSetArray.getMaxZ());
				ArrayList<LocationList> faults = new ArrayList<LocationList>();
				faults.add(faultTrace);
				double[] values = new double[faults.size()];
				for(int i=0;i<values.length;i++)
					values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;
				
				String label = imtString+"_RTGM";
				
				GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, rtgmGriddedDataSetArray, hazGridSpacing, region, true, label);
//				map.setGenerateKML(true); // this tells it to generate a KML file that can be loaded into Google Earthe
				
				// override default trace width
				for(PSXYPolygon fltTracePoly: map.getPolys())
					fltTracePoly.setPenWidth(1);
				
				// remove coast and poitical boundaries
				map.setCoast(null);

				try {
					FaultBasedMapGen.SAVE_ZIPS=true;
					FaultBasedMapGen.plotMap(hazDirFile, label, popupWindow, map);	
				} catch (GMT_MapException e) {
					e.printStackTrace();
				}
				
				// write out values to a text file
				if(hazDirName != null) {
					String fileName = hazDirName +"/"+label+".txt";
					writeHazardMapGriddedData(rtgmGriddedDataSetArray, fileName);
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	// this is the region for hazard maps
	private static GriddedRegion getGriddedRegion(FaultTrace faultTrace, double hazGridSpacing) {
//		Location loc1 = parentFaultData.getFaultTrace().get(0);
//		Location loc2 = parentFaultData.getFaultTrace().get(parentFaultData.getFaultTrace().size()-1);
//		System.out.println(loc1+"\n"+loc2);
//		System.exit(-1);
		
		// spatial buffer is 1 degree
		Location firstLoc = faultTrace.get(0);
		Location lastLoc = faultTrace.get(faultTrace.size()-1);
		double minLat = firstLoc.getLatitude()-1-hazGridSpacing/2.0;
		double maxLat = lastLoc.getLatitude()+1-hazGridSpacing/2.0;
		double minLon = firstLoc.getLongitude()-1;
		double maxLon = firstLoc.getLongitude()+1;
		Location regionCorner1 = new Location(maxLat, minLon);
		Location regionCorner2 = new Location(minLat, maxLon);
		GriddedRegion region = new GriddedRegion(regionCorner1, regionCorner2, hazGridSpacing, hazGridSpacing, null);
		return region;
	}

	
	/**
	 * This is makes to hazard maps, one for 10% in 50-year and one for 2% in 50-year ground motions.  For saPeriod of 1.0 or 0.2,
	 * this also make RTGM maps.
	 * @param fltSysRupInversion
	 * @param saPeriod - set as 0 for PGA; this will crash if SA period is not supported
	 * @param includeRTGRM
	 * @return - array with 2in50 and then 10in50 data, plus RTGM if saPeriod is 1.0 or 0.2
	 */
	public static ArrayList<GriddedGeoDataSet> makeHazardMapGriddedData(AbstractERF erf, double saPeriod, boolean includeRTGRM, 
			GriddedRegion region, ScalarIMR imr) {
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		if(D) System.out.println("Making hazard map data for "+imtString);

		
		
		// make gridded data sets
		GriddedGeoDataSet[] griddedDataSetArray = new GriddedGeoDataSet[hazardProbArray.length];
		for(int i=0;i<hazardProbArray.length;i++) {
			griddedDataSetArray[i] = new GriddedGeoDataSet(region, true);
		}
		
		GriddedGeoDataSet rtgmGriddedDataSetArray=null;
		if((saPeriod == 1.0 || saPeriod == 0.2) && includeRTGRM) // 1 or 5 Hz SA
			rtgmGriddedDataSetArray = new GriddedGeoDataSet(region, true);
		
				
		// set the attenuation relationship (GMPE)
		imr.setParamDefaults();
		
		// set the IMT & curve x-axis values
		ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
		
		// adding this to be consistent with computeHazardCurveLnX (and I commented out setting curveLinearXvalues below)
		EvenlyDiscretizedFunc tempCurveLogXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		for(int i=0;i<tempCurveLogXvalues.size();i++)
			curveLinearXvalues.set(Math.exp(tempCurveLogXvalues.getX(i)),1.0);

		
		if(saPeriod == 0) {
			imr.setIntensityMeasure(PGA_Param.NAME);
//			curveLinearXvalues = IMT_Info.getUSGS_PGA_Function();
		}
		else {
			SA_Param saParam = (SA_Param)imr.getParameter(SA_Param.NAME);
			saParam.getPeriodParam().setValue(saPeriod);
			imr.setIntensityMeasure(saParam);
//			curveLinearXvalues = IMT_Info.getUSGS_SA_01_AND_02_Function();
		}
		
		// make the site object and set values
		Site site = new Site(new Location(33d,-117d));	// Location will get over written
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
			System.out.println(param.getName()+"\t"+param.getValue());
		}
		
		ArbitrarilyDiscretizedFunc curveLogXvalues = HazardCurveSetCalculator.getLogFunction(curveLinearXvalues); // this is what the calculator expects
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
//		// test
//		calc.getHazardCurve(curveLogXvalues, site, imr, erf);
//		System.out.println(curveLogXvalues);
//		System.exit(0);
		
		// loop over sites
		int counter = -1;
		for(int i=0;i<griddedDataSetArray[0].size();i++) {
			counter+=1;
			if(counter == 500) {
				System.out.println(i+" of "+griddedDataSetArray[0].size());
				counter = 0;
			}
			site.setLocation(griddedDataSetArray[0].getLocation(i));
			calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
			// get the points on the hazard curve
			for(int p=0;p<hazardProbArray.length;p++) {
				if(hazardProbArray[p] >= curveLogXvalues.getMinY() && hazardProbArray[p] <= curveLogXvalues.getMaxY()) {
//					double logVal = curveLogXvalues.getFirstInterpolatedX(hazardProbArray[p]);
					// this is more accurate:
					double logVal = curveLogXvalues.getFirstInterpolatedX_inLogYDomain(hazardProbArray[p]);
					griddedDataSetArray[p].set(i, Math.exp(logVal));					
				}
				else {
					griddedDataSetArray[p].set(i,0);					
					System.out.println("prob outside range at i="+i+":\t"+griddedDataSetArray[0].getLocation(i));
//					System.out.println(curveLogXvalues);
//					System.exit(-1);
				}
			}
			// do RTGM is SA is 1 or 5 Hz
			if(rtgmGriddedDataSetArray != null) {
				// get annual rate curve with linear x-axis values
				boolean hasInfValue=false;
				for(int j=0;j<curveLogXvalues.size();j++) {
					// convert probability to annual rate
					double rate = -Math.log(1.0-curveLogXvalues.getY(j))/hazardDurationYrs;
					curveLinearXvalues.set(j, rate);
					if(Double.isInfinite(rate)) {
						hasInfValue = true;
//						System.out.println("Infinite Rate from: "+curveLogXvalues.getY(j));
					}
				}
				// remove Inf values that can occur if above curveLogXvalues.getY(j) = 1.0;
				ArbitrarilyDiscretizedFunc curveLinearXvaluesCleaned;
				if(hasInfValue) {
					curveLinearXvaluesCleaned = new ArbitrarilyDiscretizedFunc();
					for(int j=0;j<curveLinearXvalues.size();j++)
						if(!Double.isInfinite(curveLinearXvalues.getY(j)))
								curveLinearXvaluesCleaned.set(curveLinearXvalues.getX(j),curveLinearXvalues.getY(j));
				}
				else
					curveLinearXvaluesCleaned = curveLinearXvalues;
				// now get RTGM
//				Frequency freq;
//				if(saPeriod == 1.0)
//					freq = Frequency.SA_1P00;
//				else
//					freq = Frequency.SA_0P20;
				RTGM rtgm;
				
				// this is to write out test value to compare with the USGS on-line calculator
				try {
//					rtgm = RTGM.create(curveLinearXvaluesCleaned, freq, 0.8).call();
					rtgm = RTGM.create(curveLinearXvaluesCleaned, null, null).call();
					rtgmGriddedDataSetArray.set(i, rtgm.get());
					// this is to write one out to test against USGS web calculator; it checked out
//					if(i==0) {
//						System.out.println("TEST VALUES\n\n");
//						for(int j=0;j<curveLinearXvalues.size();j++)
//							System.out.print(curveLinearXvalues.getX(j)+",");
//						System.out.println("\n\n");
//						for(int j=0;j<curveLinearXvalues.size();j++)
//							System.out.print((float)curveLinearXvalues.getY(j)+",");
//						System.out.println("\n\nTest RTGM = "+rtgm.get());
//					}
				} catch (Exception e) {
					System.out.println("RTGM Error; Hazard curve is:\n"+curveLinearXvaluesCleaned);
					System.out.println("Location is:\n"+griddedDataSetArray[0].getLocation(i));
					e.printStackTrace();
					System.exit(0);
				}
			}
		}
		
		ArrayList<GriddedGeoDataSet> returnList = new ArrayList<GriddedGeoDataSet>();
		for(GriddedGeoDataSet dataSet : griddedDataSetArray)
			returnList.add(dataSet);
		if(rtgmGriddedDataSetArray != null) {
			returnList.add(rtgmGriddedDataSetArray);
		}
		return returnList;
	}

	
	
	/**
	 * This plots the classic MFDs for the paper
	 */
	public static void plotMFDs(ArrayList<IncrementalMagFreqDist> mfdList, double yMin, double yMax, String fileNamePrefix) {
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(mfdList.get(0));
		funcs.add(mfdList.get(1));
		funcs.add(mfdList.get(2));
		funcs.add(mfdList.get(3));
		funcs.add(mfdList.get(4));
		funcs.add(mfdList.get(0).getCumRateDistWithOffset());
		funcs.add(mfdList.get(1).getCumRateDistWithOffset());
		funcs.add(mfdList.get(2).getCumRateDistWithOffset());
		funcs.add(mfdList.get(3).getCumRateDistWithOffset());
		funcs.add(mfdList.get(4).getCumRateDistWithOffset());
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.MAGENTA));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.MAGENTA));
		
		Range xRange = new Range(6.5,8.2);
		Range yRange = new Range(yMin, yMax);
		
		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Magnitude", "Rate (per yr)", xRange, yRange, 
				false, true, 3.5, 3.0, fileNamePrefix, true);	
	}
	
	
	public static void equivBvalVsMmax() {
		
		double grMagMin = 6.55;
		double toMoRate = 1e18;
		double grBval = 1d;
		for(double maxMag=6.75; maxMag<8; maxMag+=0.1) {
			GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(grMagMin,20,0.1); // up to M ~8.5
			grMFD.setAllButTotCumRate(grMagMin, maxMag, toMoRate/3.0, grBval);
			grMFD.setName("GR");

			GaussianMagFreqDist charMFD = new GaussianMagFreqDist(grMagMin, 20, 0.1, 
					maxMag, 0.12, 2.0*toMoRate/3.0, 2.0, 2);
			
			SummedMagFreqDist summedMFD = new SummedMagFreqDist(grMagMin,20,0.1);
			summedMFD.setName("summedMFD");
			summedMFD.addIncrementalMagFreqDist(grMFD);
			summedMFD.addIncrementalMagFreqDist(charMFD);
			double bValFit = summedMFD.compute_bValueAlt(grMagMin,maxMag);
			System.out.println((float)maxMag+"\t"+(float)bValFit);
		}

	}
	
	public static U3FaultSystemSolution getClassicFSS(FaultSectionPrefData faultData, double ddw,
			GutenbergRichterMagFreqDist grMFD, GaussianMagFreqDist charMFD, Ellsworth_B_WG02_MagAreaRel ellB_magArea) {
		
		ArrayList<FaultSectionPrefData> faultSectionDataList = faultData.getSubSectionsList(ddw/2.0+0.01); // add a bit to get the correct number of subsections
		
//		System.out.println(faultSectionDataList.get(0).toString());
//		System.exit(0);
		
		// set all the section ids and compute length
		double lengthSum = 0;
		for(int i=0; i<faultSectionDataList.size() ;i++) {
			faultSectionDataList.get(i).setSectionId(i);
			lengthSum += faultSectionDataList.get(i).getTraceLength();
		}
		
		System.out.println("num subsect: "+faultSectionDataList.size());
		System.out.println("subsection trace length: "+faultSectionDataList.get(0).getTraceLength());
		
//		double totLengthFromSections = faultSectionDataList.get(0).getTraceLength()*faultSectionDataList.size();
//		System.out.println("Length Tests (parent, subSec sums): "+faultData.getTraceLength()+"\t"+totLengthFromSections+"\t"+lengthSum);

		if(Math.abs(faultData.getTraceLength()-lengthSum)>1)
			throw new RuntimeException("Trace length problem in subsections");
		
		int NUM_SUBSECT_PER_RUP = 2;
		int numSections = faultSectionDataList.size();
		int numRupturesGR=0;
		int curNumRup = numSections - (NUM_SUBSECT_PER_RUP-1);
		while(curNumRup>0) {
			numRupturesGR += curNumRup;
			curNumRup -= 1;
		}
		int rupIndex=0;
		int[][] rupSectionMatrixGR = new int[numSections][numRupturesGR];
		for(int curSectPerRup = NUM_SUBSECT_PER_RUP; curSectPerRup<numSections+1; curSectPerRup++) {
			for(int s=0;s<numSections-curSectPerRup+1;s++) {
				int firstSect = s;
				int lastSect = s+curSectPerRup-1;
//				System.out.println(rupIndex+":\t"+firstSect+"\t"+lastSect);
				for(int col=firstSect; col <= lastSect; col++)
					rupSectionMatrixGR[col][rupIndex] = 1;
				rupIndex += 1;
			}
		}		


		List<List<Integer>> sectionsForRups = Lists.newArrayList();;
		for(int r=0;r<numRupturesGR;r++) {
			ArrayList<Integer> sectList = new ArrayList<Integer>();
			for(int s=0;s<numSections;s++) {
				if(rupSectionMatrixGR[s][r] == 1)
					sectList.add(s);
			}
			sectionsForRups.add(sectList);
		}
		
		double[] sectSlipRate = new double[numSections];
		double[] sectSlipRateStdDev = new double[numSections];
		double[] sectArea = new double[numSections];
		for(int s=0;s<faultSectionDataList.size();s++) {
			FaultSectionPrefData fltData = faultSectionDataList.get(s);
			sectSlipRate[s] = fltData.getOrigAveSlipRate();
			sectSlipRateStdDev[s] = fltData.getOrigAveSlipRate()/10.;
			sectArea[s] = fltData.getArea(false);
		}
		
		ArrayList<Integer> charNonZeroIndices = new ArrayList<Integer>();
		for(int i=0;i<charMFD.size();i++)
			if(charMFD.getY(i)>0)
				charNonZeroIndices.add(i);
		
		int totNumRups = numRupturesGR + charNonZeroIndices.size();
		double[] rupMeanMag = new double[totNumRups];
		double[] rupAveRake = new double[totNumRups];
		double[] rupArea = new double[totNumRups];
		double[] rupLength = new double[totNumRups];
		
		double[] rupRateSolution = new double[totNumRups];
		
		// fill in GR rupture attributes
		HistogramFunction numRupForMagFunc = new HistogramFunction(grMFD.getMinX(), grMFD.size(), grMFD.getDelta());
		for(int r=0;r<sectionsForRups.size();r++) {
			int numSectForRup = sectionsForRups.get(r).size();
			rupLength[r] = numSectForRup*faultSectionDataList.get(0).getTraceLength();
			rupArea[r] = rupLength[r]*ddw;
			rupAveRake[r]= faultSectionDataList.get(0).getAveRake();
			double mag = rupMeanMag[r] = ellB_magArea.getMedianMag(rupArea[r]);
			if(mag<6.5) mag = 6.51; // this converts the 2-section ruptures from 6.49 to 6.51, which also eliminates the zero bin
			rupMeanMag[r] = mag;
//			if(r< 120) System.out.println(rupMeanMag[r]);
			numRupForMagFunc.add(rupMeanMag[r], 1.);
		}
		System.out.println("numRupForMagFunc:\n"+numRupForMagFunc);
		
		
//		ArrayList<Double> mags = new ArrayList<Double>();
//		for(int i=0;i<numRupForMagFunc.size();i++)
//			if(numRupForMagFunc.getY(i) > 0.0)
//				mags.add(numRupForMagFunc.getX(i));		
//		System.out.println("mags:\n"+mags);
//		IncrementalMagFreqDist grMFD_wZeros = SparseGutenbergRichterSolver.getEquivGR(grMFD, mags, grMFD.getTotalMomentRate(), grMFD.get_bValue());

		for(int r=0;r<sectionsForRups.size();r++) {
			double mag = rupMeanMag[r];
			int xIndex = grMFD.getClosestXIndex(mag);
			rupRateSolution[r] = grMFD.getY(xIndex)/numRupForMagFunc.getY(xIndex);
		}
		
		// fill in char rupture attributes
		for(int r=numRupturesGR;r<numRupturesGR+charNonZeroIndices.size();r++) {
			rupLength[r] = rupLength[numRupturesGR-1];
			rupArea[r] = rupArea[numRupturesGR-1];
			rupAveRake[r]= rupAveRake[numRupturesGR-1];
			sectionsForRups.add(sectionsForRups.get(numRupturesGR-1));
			rupMeanMag[r] = charMFD.getX(charNonZeroIndices.get(r-numRupturesGR));
			rupRateSolution[r]= charMFD.getY(charNonZeroIndices.get(r-numRupturesGR));
		}

		U3FaultSystemRupSet rupSet = new U3FaultSystemRupSet(
				faultSectionDataList,
				sectSlipRate,
				sectSlipRateStdDev,
				sectArea,
				sectionsForRups,
				rupMeanMag,
				rupAveRake,
				rupArea,
				rupLength,
				"ClassicModelFRS");

		return new U3FaultSystemSolution(rupSet, rupRateSolution);
		
	}
	
	
	
	/**
	 * This generates a solution and diagnostic info for the current parameter settings
	 */
	public static FaultSystemRuptureRateInversion getSA_Solution(IncrementalMagFreqDist targetMFD, ScalingRelationshipEnum scalingRel, 
		String dirName, boolean popUpPlots, FaultSectionPrefData faultData, double ddw) {	

		int NUM_SUBSECT_PER_RUP = 2;
		
		long randomSeed = System.currentTimeMillis();
		
		CoolingScheduleType saCooling = CoolingScheduleType.FAST_SA;
		GenerationFunctionType perturbationFunc = GenerationFunctionType.UNIFORM_0p001;
		CompletionCriteria completionCriteria = new IterationCompletionCriteria((long) 1e7);
				
		IncrementalMagFreqDist mfdSigma = null;
		if(targetMFD != null && targetMFD.getMinY()>0) {
			mfdSigma = targetMFD.deepClone();
			mfdSigma.scale(0.1); // uncertainty is 10%
		}
		
		ArrayList<FaultSectionPrefData> faultSectionDataList = faultData.getSubSectionsList(ddw/NUM_SUBSECT_PER_RUP+0.01); // add a bit to get the correct number of subsections
		
		// set all the section ids and compute length
		double lengthSum = 0;
		for(int i=0; i<faultSectionDataList.size() ;i++) {
			faultSectionDataList.get(i).setSectionId(i);
			lengthSum += faultSectionDataList.get(i).getTraceLength();
		}
		System.out.println("num subsect: "+faultSectionDataList.size());
		System.out.println("subsection trace length: "+faultSectionDataList.get(0).getTraceLength());
//		double totLengthFromSections = faultSectionDataList.get(0).getTraceLength()*faultSectionDataList.size();
//		System.out.println("Length Tests (parent, subSec sums): "+faultData.getTraceLength()+"\t"+totLengthFromSections+"\t"+lengthSum);
		if(Math.abs(faultData.getTraceLength()-lengthSum)>1)
			throw new RuntimeException("Trace length problem in subsections");
		
		int numSections = faultSectionDataList.size();
		int numRuptures=0;
		int curNumRup = numSections - (NUM_SUBSECT_PER_RUP-1);
		while(curNumRup>0) {
			numRuptures += curNumRup;
			curNumRup -= 1;
		}
		int rupIndex=0;
		int[][] rupSectionMatrix = new int[numSections][numRuptures];
		for(int curSectPerRup = NUM_SUBSECT_PER_RUP; curSectPerRup<numSections+1; curSectPerRup++) {
			for(int s=0;s<numSections-curSectPerRup+1;s++) {
				int firstSect = s;
				int lastSect = s+curSectPerRup-1;
//				System.out.println(rupIndex+":\t"+firstSect+"\t"+lastSect);
				for(int col=firstSect; col <= lastSect; col++)
					rupSectionMatrix[col][rupIndex] = 1;
				rupIndex += 1;
			}
		}		


		List<List<Integer>> sectionsForRups = Lists.newArrayList();;
		for(int r=0;r<numRuptures;r++) {
			ArrayList<Integer> sectList = new ArrayList<Integer>();
			for(int s=0;s<numSections;s++) {
				if(rupSectionMatrix[s][r] == 1)
					sectList.add(s);
			}
			sectionsForRups.add(sectList);
		}
		
		String solutionName = "SimulatedAnnealing";
		SlipAlongRuptureModelEnum slipModelType = SlipAlongRuptureModelEnum.UNIFORM; 
		ArrayList<SlipRateSegmentationConstraint> slipRateSegmentationConstraintList = new ArrayList<SlipRateSegmentationConstraint> ();
		ArrayList<SectionRateConstraint> sectionRateConstraintList = new ArrayList<SectionRateConstraint>();
		double relativeSectRateWt = 0; 
		double relative_aPrioriRupWt = 0; 
		String aPrioriRupRateFilename = null;
		boolean wtedInversion = true;
		double minRupRate = 0.0;
		boolean applyProbVisible = false;
		double moRateReduction = 0.0; 
		double relativeMFD_constraintWt = 1.0;
		ArrayList<SegmentationConstraint> segmentationConstraintList = new ArrayList<SegmentationConstraint> ();
		double relative_segConstraintWt = 0.0;
		double totalRateConstraint = 0.0;
		double totalRateSigma = 0.1*totalRateConstraint;
		double relativeTotalRateConstraintWt = 0;
		ArrayList<int[]> smoothnessConstraintList = null;
		double relativeSmoothnessConstraintWt = 0;
		double magAareaAleatoryVariability = 0;
		
		
		// this will be used to keep track of runtimes
		long startTimeMillis = System.currentTimeMillis();
		if(D)
			System.out.println("Starting Inversion");
		
		// create an instance of the inversion class with the above settings
		FaultSystemRuptureRateInversion fltSysRupInversion = new  FaultSystemRuptureRateInversion(
				solutionName,
				"UniformSlipRate",
				faultSectionDataList, 
				rupSectionMatrix, 
				slipModelType, 
				scalingRel, 
				slipRateSegmentationConstraintList,
				sectionRateConstraintList, 
				relativeSectRateWt, 
				relative_aPrioriRupWt, 
				aPrioriRupRateFilename,
				wtedInversion, 
				minRupRate, 
				applyProbVisible, 
				moRateReduction,
				targetMFD,
				mfdSigma,
				relativeMFD_constraintWt,
				segmentationConstraintList,
				relative_segConstraintWt,
				totalRateConstraint,
				totalRateSigma,
				relativeTotalRateConstraintWt,
				smoothnessConstraintList,
				relativeSmoothnessConstraintWt,
				magAareaAleatoryVariability);

		
		// make the directory for storing results
		if(dirName != null) {
		    File file = new File(dirName);
		    file.mkdirs();			
		}
	    
	    // set a-prior rates from MFD so these can be applied as initial model
		boolean setAprioriRatesFromMFD_Constraint = false;
	    if(setAprioriRatesFromMFD_Constraint)
	    	fltSysRupInversion.setAprioriRupRatesFromMFD_Constrint();
	    
	    // write the setup info to a file
		if(dirName != null)
			fltSysRupInversion.writeInversionSetUpInfoToFile(dirName);

		// do this before overriding solutionType for rePlotOnly case (otherwise rupSamplerMFD won't be plotted)
		IntegerPDF_FunctionSampler rupSampler = null;
		boolean applyRuptureSampler = true;
		// set the rupture sampler; default is GR if targetMFD == null
		if(applyRuptureSampler) {
			double[] rupSampleProbArray;
			rupSampleProbArray = fltSysRupInversion.getRupRatesForTargetMFD(targetMFD, false);
			rupSampler = new IntegerPDF_FunctionSampler(rupSampleProbArray.length);
			for(int r=0; r<rupSampleProbArray.length; r++)
				rupSampler.set(r,rupSampleProbArray[r]);	
			// do this here so re-plot includes this MFD
			// System.out.println("targetMFD:\n"+targetMFD);
			fltSysRupInversion.computeRupSamplerMFD(rupSampler);
		}
	    
		// set the initial state & rupSampler if SA
	    double[] initialState = null;
	    //  ALL_ZEROS:
	    	initialState = new double[fltSysRupInversion.getNumRuptures()];
	    //  A_PRIORI_RATES:
//	    	double tempArray[] = fltSysRupInversion.getAprioriRuptureRates();
//	    	if(tempArray.length != fltSysRupInversion.getNumRuptures())
//	    		throw new RuntimeException("aPrior rates are not set for each rupture");
//	    	initialState = tempArray;
	    //	FROM_MFD_CONSTRAINT:
//	    		initialState = fltSysRupInversion.getRupRatesForTargetMFD(targetMFD, false);
	    
    	fltSysRupInversion.doInversionSA(completionCriteria, initialState, randomSeed, saCooling, perturbationFunc, rupSampler);
    	
		double runTimeSec = ((double)(System.currentTimeMillis()-startTimeMillis))/1000.0;
		
		String runTimeString = "Done with Inversion after "+(float)runTimeSec+" seconds.";
		if(D) System.out.println(runTimeString);
		
		// Write out info & make plots
		if(dirName != null) { 
			fltSysRupInversion.addToModelRunInfoString("\n"+runTimeString+"\n");
			// write results to file
			fltSysRupInversion.writeInversionRunInfoToFile(dirName);
			fltSysRupInversion.writeRuptureRatesToFile(dirName);
		
			// plot various model related histograms
			fltSysRupInversion.writeAndOrPlotMagHistograms(dirName, popUpPlots, 3.5, 3.0);

			// plot solution MFDs
			fltSysRupInversion.writeAndOrPlotMFDs(dirName, popUpPlots, null, new Range(1e-5,1), 3.5, 3.0);
		
// for this one need to cut and paste mothod from SimpleFaultInversion:
			// plot various things as a function of section index
//			makeResultsVsSectionPlots(fltSysRupInversion);
		
			// plot normalized residual versus row index
			fltSysRupInversion.writeAndOrPlotNormalizedResiduals(dirName, popUpPlots, 3.5, 3.0);
		
			// plot section boundary rates
			fltSysRupInversion.writeAndOrPlotSectionBoundaryRates(dirName, popUpPlots, null, 6.5, 1.7);
		
			// plot rupture rate versus index (big and not very informative plots)
//			fltSysRupInversion.writeAndOrPlotRupRateVsIndex(dirName, popUpPlots, 6.5, 4.0);
		
			// plot rupture event and slip rate as a function of section for all non-zero ruptures (big and not very informative, except to confirm that tapered slip working)
//			fltSysRupInversion.writeAndOrPlotNonZeroRateRups(dirName, popUpPlots, 9.0, 6.5);

			// plot the rate of each rupture (triangle plot with rate of first section)
			Range zRange = new Range(-9, -1);
			fltSysRupInversion.writeAndOrPlotRupRateOfFirstSection(dirName, popUpPlots, 3.5, 4.0, zRange);

		}
		
		if(D) System.out.println("Done getting solution");

		return fltSysRupInversion;
		
	}

	
	public static void mkMFD_PlotForSSA_Talk() {
		double grMagMin = 6.55;
		double grMagMax = 7.75;
		double grBval = 1d;
		double toMoRate = 1e18;  // 1e18 has problems with 10% in 50 yr
		
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(grMagMin,200,0.01); // up to M ~8.5
		grMFD.setAllButTotCumRate(grMagMin, grMagMax, toMoRate/3.0, grBval);
		grMFD.setName("GR");

		GaussianMagFreqDist charMFD = new GaussianMagFreqDist(grMagMin, 200, 0.01, 
				grMagMax, 0.12, 2.0*toMoRate/3.0, 2.0, 2);
		
		SummedMagFreqDist summedMFD = new SummedMagFreqDist(grMagMin,200,0.01);
		summedMFD.setName("summedMFD");
		summedMFD.addIncrementalMagFreqDist(grMFD);
		summedMFD.addIncrementalMagFreqDist(charMFD);
		double bValFit = summedMFD.compute_bValueAlt(grMagMin,grMagMax);
		summedMFD.setInfo("Equiv. b-value = "+(float)bValFit);
		
		// Plot MFDs
		ArrayList<XY_DataSet> mfdList = new ArrayList<XY_DataSet>();
		mfdList.add(grMFD);
		mfdList.add(charMFD);
		mfdList.add(summedMFD);
		
	    String fileNameMFD = "MFD_ForSSA2022";

		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		
		Range xRange = new Range(6.5,8.0);
		Range yRange = new Range(5e-6, 5e-4);
		
		PlottingUtils.writeAndOrPlotFuncs(mfdList, plotChars, null, "Magnitude", "Rate (per yr)", xRange, yRange, 
				false, true, 3.5, 3.0, fileNameMFD, true);	
	}
	


	public static void main(String[] args) {
		
		mkMFD_PlotForSSA_Talk();
		System.exit(0);
		
//		equivBvalVsMmax();
//		System.exit(0);
		
		boolean makeMapData=true;
		double grMagMin = 6.55;
//		double grTargetMagMax = 6.95;
//		double grTargetMagMax = 7.45;
		double grTargetMagMax = 7.95;
		double grBval = 1d;
		double toMoRate = 1e19;  // 1e18 has problems with 10% in 50 yr
		double floatOffset = 1;
		double ddw = 14;
		double rake = 0;
		double duration = hazardDurationYrs;
		
		Ellsworth_B_WG02_MagAreaRel ellB_magArea = new Ellsworth_B_WG02_MagAreaRel();
		
		double faultLength = ellB_magArea.getMedianArea(grTargetMagMax)/ddw;
		GutenbergRichterMagFreqDist junkMFD = new GutenbergRichterMagFreqDist(grMagMin,20,0.1); // up to M ~8.5

		// get distance that is an integer multiple of ddw/2
		System.out.println("faultLength="+faultLength+"\ntargetMag="+grTargetMagMax);
		int num = (int)(faultLength/(ddw/2.0));
		double cutLength = num*(ddw/2.0);
		double cutMag = ellB_magArea.getMedianMag(cutLength*ddw);
		System.out.println("orig num="+num+"\tcutLength="+cutLength+"\tcutMag="+cutMag);
		num+=1;
		double finalLength = num*(ddw/2.0);
		double finalMag = ellB_magArea.getMedianMag(finalLength*ddw);
		System.out.println("final num="+num+"\tfinalLength="+finalLength+"\tfinalMag="+(float)finalMag);
		
		double finalRoundedMag = junkMFD.getX(junkMFD.getClosestXIndex(finalMag));
		System.out.println("finalRoundedMag="+(float)finalRoundedMag);
		
		if(Math.abs(finalRoundedMag-grTargetMagMax)>0.01)
			throw new RuntimeException("Bad rounding");
	
		double grMagMax = finalRoundedMag;
		faultLength = finalLength;
		
//		   area: the fault area (in square Meters)
//		   moment:(in Newton-Meters) or moment rate 
//		   * @return Slip (in meters) or slip rate if moment-rate is given
		double slipRate = FaultMomentCalc.getSlip(faultLength*ddw*1e6, toMoRate)*1000; // mm/yr
		System.out.println("slip rate (MM/yr) ="+slipRate);

		FaultTrace faultTrace = new FaultTrace("FaultTrace");
		faultTrace.add(new Location(33d,-117d));
		faultTrace.add(new Location(33d+faultLength/111.195052,-117));
		System.out.println("Fault trace length: "+faultTrace.getTraceLength());
		FaultSectionPrefData faultData = new FaultSectionPrefData();
		faultData.setFaultTrace(faultTrace);
		faultData.setAveDip(90);
		faultData.setAveRake(rake);
		faultData.setAveUpperDepth(0.0);
		faultData.setAveLowerDepth(ddw);
		faultData.setAveSlipRate(slipRate);
		faultData.setSlipRateStdDev(slipRate*0.1);
		StirlingGriddedSurface faultSurf = faultData.getFaultSurface(1.0);
		
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(grMagMin,20,0.1); // up to M ~8.5
		grMFD.setAllButTotCumRate(grMagMin, grMagMax, toMoRate/3.0, grBval);
		grMFD.setName("GR");

		GaussianMagFreqDist charMFD = new GaussianMagFreqDist(grMagMin, 20, 0.1, 
				grMagMax, 0.12, 2.0*toMoRate/3.0, 2.0, 2);
		
		SummedMagFreqDist summedMFD = new SummedMagFreqDist(grMagMin,20,0.1);
		summedMFD.setName("summedMFD");
		summedMFD.addIncrementalMagFreqDist(grMFD);
		summedMFD.addIncrementalMagFreqDist(charMFD);
		double bValFit = summedMFD.compute_bValueAlt(grMagMin,grMagMax);
		summedMFD.setInfo("Equiv. b-value = "+(float)bValFit);
		
		GutenbergRichterMagFreqDist summedEquivGR_MFD = new GutenbergRichterMagFreqDist(grMagMin,20,0.1); // up to M ~8.5
		summedEquivGR_MFD.setAllButBvalue(grMagMin, grMagMax, toMoRate, summedMFD.getTotalIncrRate());
		summedEquivGR_MFD.setName("summedEquivGR_MFD");
		summedEquivGR_MFD.setInfo("b-value = "+summedEquivGR_MFD.get_bValue());
		
		
		// Classic model as FSS:
		U3FaultSystemSolution classicFSS = getClassicFSS(faultData, ddw, grMFD, charMFD, ellB_magArea);
		IncrementalMagFreqDist classicFSS_MFD = classicFSS.calcTotalNucleationMFD(summedMFD.getMinX(), summedMFD.getMaxX(), summedMFD.getDelta());
		classicFSS_MFD.setName("classicFSS_MFD");
		// Create the ERF
		FaultSystemSolutionERF classicFSS_ERF = new FaultSystemSolutionERF(classicFSS);
		classicFSS_ERF.setName("classicFSS ERF");
		classicFSS_ERF.getTimeSpan().setDuration(duration);
		classicFSS_ERF.updateForecast();
		
		// Simulated Annealing Inversion
		int numPts = (int)Math.round((grMagMax-grMagMin)/0.1) + 1;
		GutenbergRichterMagFreqDist scratchTargetMFD = new GutenbergRichterMagFreqDist(grMagMin-0.1,numPts+1,0.1); // up to M ~8.5
		ArrayList<Double> mags = new ArrayList<Double>();
		mags.add(scratchTargetMFD.getX(0));
		for(int i=2;i<scratchTargetMFD.size();i++) // skip the second element
			mags.add(scratchTargetMFD.getX(i));
		IncrementalMagFreqDist inversionMFD = SparseGutenbergRichterSolver.getEquivGR(scratchTargetMFD, mags, toMoRate, summedEquivGR_MFD.get_bValue());
//		System.out.println(summedEquivGR_MFD+"\n\n"+inversionMFD);
		ScalingRelationshipEnum scalingRel = ScalingRelationshipEnum.ELLSWORTH_B;
		boolean popUpPlots = true;
		String dirNameSA = ROOT_PATH+"/SA_Run_Mmax"+(float)grMagMax;
		FaultSystemRuptureRateInversion fltInversion = getSA_Solution(inversionMFD, scalingRel, 
				dirNameSA, popUpPlots, faultData, ddw);
		FaultSystemSolutionERF erf_sa = new FaultSystemSolutionERF(fltInversion.getFaultSystemSolution());
		erf_sa.setName("SA ERF");
		erf_sa.getTimeSpan().setDuration(duration);
		erf_sa.updateForecast();
		
		// Plot MFDs
		ArrayList<IncrementalMagFreqDist> mfdList = new ArrayList<IncrementalMagFreqDist>();
		mfdList.add(grMFD);
		mfdList.add(charMFD);
		mfdList.add(summedMFD);
		mfdList.add(summedEquivGR_MFD);
		mfdList.add(classicFSS_MFD);

		double minPlotY = Math.pow(10, Math.floor(Math.log10(grMFD.getY(grMagMax))));
		double maxPlotY = Math.pow(10, Math.ceil(Math.log10(summedMFD.getTotalIncrRate())));
		
		
	    String fileNameMFD = ROOT_PATH+"/MFD_Mmax"+(float)grMagMax;

		plotMFDs(mfdList, minPlotY, maxPlotY, fileNameMFD);
		
//		System.exit(0);

		
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);

//		System.out.println(grMFD);
//		System.out.println(charMFD);
//		System.out.println(summedMFD);
		
		ArrayList<ProbEqkSource> sourceList = new ArrayList<ProbEqkSource>();
		
		FaultRuptureSource charSrc = new FaultRuptureSource(charMFD, faultSurf, rake, duration);
		
		FloatingPoissonFaultSource grSrc = new FloatingPoissonFaultSource(grMFD, faultSurf,
				ellB_magArea, 0d, 1d, floatOffset, rake, duration, grMagMin, 0, grMagMax);
		
		
		sourceList.add(charSrc);
		sourceList.add(grSrc);
		ArbSrcListERF erf = new ArbSrcListERF(sourceList, duration);
		erf.getTimeSpan().setDuration(50);
		System.out.println("erf.getNumSources = "+erf.getNumSources());
//		for(int s=0;s<erf.getNumSources();s++) {
//			ProbEqkSource src = erf.getSource(s);
//			for(int r=0;r<src.getNumRuptures();r++) {
//				ProbEqkRupture rup = src.getRupture(r);
//				System.out.println((float)rup.getMag()+"\t"+(float)rup.getProbability()+"\t"+(float)rup.getRuptureSurface().getFirstLocOnUpperEdge().getLatitude());
//			}
//		}
		System.out.println("ERF calc MoRate = "+(float)ERF_Calculator.getTotalMomentRateInRegion(erf, null));
//		System.exit(0);
		
		
		double[] saPeriodForHazArray = {0.0, 1.0};
//		double[] saPeriodForHazArray = {0.0};
		boolean popupWindow = true;
		GriddedRegion region = getGriddedRegion(faultTrace,hazGridSpacing);
		
		String solutionName1 = "ClassicNSHMP_Mmax"+(float)grMagMax; 
		String dirName = ROOT_PATH+solutionName1;

		for(double saPeriodForHaz : saPeriodForHazArray) {
			String hazDirName = null;
			if(dirName != null) {
				hazDirName = dirName+"/hazardMaps";
			}
			// make mean map
			
			if(makeMapData)
				makeHazardMaps(erf, saPeriodForHaz, hazDirName, popupWindow, region, imr, faultTrace);
		}
		
		
		
		// NOW do GR Equiv Case
		FloatingPoissonFaultSource grEquivSrc = new FloatingPoissonFaultSource(summedEquivGR_MFD, faultSurf,
				ellB_magArea, 0d, 1d, floatOffset, rake, duration, grMagMin, 0, grMagMax);
		sourceList = new ArrayList<ProbEqkSource>();
		sourceList.add(grEquivSrc);
		erf = new ArbSrcListERF(sourceList, duration);
		erf.getTimeSpan().setDuration(50);

		String solutionName2 = "GR_EquivToClassic_Mmax"+(float)grMagMax; 
		dirName = ROOT_PATH+solutionName2;
		for(double saPeriodForHaz : saPeriodForHazArray) {
			String hazDirName = null;
			if(dirName != null) {
				hazDirName = dirName+"/hazardMaps";
			}
			// make mean map
			
			if(makeMapData)
				makeHazardMaps(erf, saPeriodForHaz, hazDirName, popupWindow, region, imr, faultTrace);
		}
		
		
		String solutionName3 = "ClassicFSS_Mmax"+(float)grMagMax; 
		dirName = ROOT_PATH+solutionName3;
		for(double saPeriodForHaz : saPeriodForHazArray) {
			String hazDirName = null;
			if(dirName != null) {
				hazDirName = dirName+"/hazardMaps";
			}
			if(makeMapData)
				makeHazardMaps(classicFSS_ERF, saPeriodForHaz, hazDirName, popupWindow, region, imr, faultTrace);
		}
		
		String solutionName4 = "SA_InversionMmax"+(float)grMagMax; 
		dirName = ROOT_PATH+solutionName4;
		for(double saPeriodForHaz : saPeriodForHazArray) {
			String hazDirName = null;
			if(dirName != null) {
				hazDirName = dirName+"/hazardMaps";
			}
			if(makeMapData)
				makeHazardMaps(erf_sa, saPeriodForHaz, hazDirName, popupWindow, region, imr, faultTrace);
		}


		
		
		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
		for(String name: nameArray) {
		    String fileName1 = ROOT_PATH+solutionName2+"/hazardMaps/"+name+".txt";
		    String fileName2 = ROOT_PATH+solutionName1+"/hazardMaps/"+name+".txt";
		    dirName = ROOT_PATH+solutionName2+"/hazardMaps";
		    String label = name+"_RatioToClassic";
			makeHazardMapRatio(fileName1, fileName2, label, dirName, true, region, faultTrace);		
			
		    fileName1 = ROOT_PATH+solutionName3+"/hazardMaps/"+name+".txt";
		    dirName = ROOT_PATH+solutionName3+"/hazardMaps";
		    label = name+"_RatioToClassic";
			makeHazardMapRatio(fileName1, fileName2, label, dirName, true, region, faultTrace);		
		
		    fileName1 = ROOT_PATH+solutionName4+"/hazardMaps/"+name+".txt";
		    dirName = ROOT_PATH+solutionName4+"/hazardMaps";
		    label = name+"_RatioToClassic";
			makeHazardMapRatio(fileName1, fileName2, label, dirName, true, region, faultTrace);		

		}
	}
}
