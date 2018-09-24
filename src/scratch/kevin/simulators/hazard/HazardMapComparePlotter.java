package scratch.kevin.simulators.hazard;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPOutputStream;

import javax.imageio.ImageIO;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculator;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.DisaggregationPlotViewerWindow;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.NEHRP_TestCity;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

public class HazardMapComparePlotter {
	
	private static final boolean map_parallel = true;
	private static final boolean disagg_parallel = true;

	@SuppressWarnings("unused")
	public static void main(String[] args) throws Exception {
		
		boolean u3SupraMinMag = false;
		double minMag = 6.5d;
		
		File hazardJobDir, catalogsBaseDir, mainOutputDir;
		
		if (args.length == 0) {
			hazardJobDir = new File("/home/kevin/Simulators/hazard/");
			catalogsBaseDir = new File("/data/kevin/simulators/catalogs");
			mainOutputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		} else {
			Preconditions.checkArgument(args.length >= 3, "USAGE: <root-haz-dir> <catalog-base-dir> <git-catalogs-dir> "
					+ "[<catalog name> <dir-name-1> ... <dir-name-N>]");
			hazardJobDir = new File(args[0]);
			catalogsBaseDir = new File(args[1]);
			mainOutputDir = new File(args[2]);
		}
		
		List<File> jobDirs = Lists.newArrayList();
		
		String catalogName = "RSQSim";
		String catalogFileName = "rsqsim";
		RSQSimCatalog catalog;
		
		if (args.length > 3) {
			catalog = Catalogs.valueOf(args[3]).instance(catalogsBaseDir);
			
			for (int i=4; i<args.length; i++)
				jobDirs.add(new File(hazardJobDir, args[i]));
		} else {
			FaultBasedMapGen.LOCAL_MAPGEN = true;
			catalog = Catalogs.BRUCE_2585.instance(catalogsBaseDir);
			
			jobDirs.add(new File(hazardJobDir, "2018_02_16-bruce2585-m6.5-sectArea0.2-skip5000yr-pga-8xPoints-maxDist1000"));
			jobDirs.add(new File(hazardJobDir, "2018_02_16-bruce2585-m6.5-sectArea0.2-skip5000yr-sa-0.2s-8xPoints-maxDist1000"));
			jobDirs.add(new File(hazardJobDir, "2018_02_16-bruce2585-m6.5-sectArea0.2-skip5000yr-sa-1.0s-8xPoints-maxDist1000"));
			jobDirs.add(new File(hazardJobDir, "2018_02_16-bruce2585-m6.5-sectArea0.2-skip5000yr-sa-2.0s-8xPoints-maxDist1000"));
			jobDirs.add(new File(hazardJobDir, "2018_02_16-bruce2585-m6.5-sectArea0.2-skip5000yr-sa-5.0s-8xPoints-maxDist1000"));
			jobDirs.add(new File(hazardJobDir, "2018_02_16-bruce2585-m6.5-sectArea0.2-skip5000yr-sa-10.0s-8xPoints-maxDist1000"));
		}
		
		File catOutDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catOutDir.exists() || catOutDir.mkdir());
		
		
//		String catalogName = "UCERF2";
//		String catalogFileName = "ucerf2";
////		jobDirs.add(new File(mainDir, "2017_09_05-ucerf2-faults-m6.5-pga-8xPoints"));
////		jobDirs.add(new File(mainDir, "2017_09_05-ucerf2-faults-m6.5-sa-0.2s-8xPoints"));
//		jobDirs.add(new File(mainJobDir, "2017_09_05-ucerf2-faults-m6.5-sa-1.0s-8xPoints"));
////		jobDirs.add(new File(mainDir, "2017_09_05-ucerf2-faults-m6.5-sa-2.0s-8xPoints"));
////		jobDirs.add(new File(mainDir, "2017_09_05-ucerf2-faults-m6.5-sa-5.0s-8xPoints"));
////		jobDirs.add(new File(mainDir, "2017_09_05-ucerf2-faults-m6.5-sa-10.0s-8xPoints"));
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		Map<Double, File> u3FullFiles = new HashMap<>();
		u3FullFiles.put(Double.NaN, new File("/home/kevin/OpenSHA/UCERF3/maps/"
				+ "2017_07_14-ucerf3-gridded-tests/full/curves/imrs1.bin"));
		u3FullFiles.put(1d, new File("/home/kevin/OpenSHA/UCERF3/maps/"
				+ "2017_09_06-ucerf3-geol-gridded-tests-sa-1.0s/full/curves/imrs1.bin"));
		
		boolean plotAll = true;
		
		boolean plotMaps = true || plotAll;
		boolean plotHist = true || plotAll;
		boolean plotCurves = true || plotAll;
//		boolean plotDisagg = false || (plotAll && !catalogName.contains("UCERF2"));
		boolean plotDisagg = false;
		boolean plotCurveProfile = false; // only if specified
		
		int[] mapRPs = { 1000, 2500, 10000};
		int[] histRPs = { 1000, 2500, 5000, 10000 };
		int[] nehrpRPs = { 2500 };
		int[] disaggRPs = { 2500, 10000 };
		int histHighlightIndex = 1;
		
		Map<String, Location> curveLocs = new HashMap<>();
		curveLocs.put("Pasadena", new Location(34.148426, -118.17119));
		curveLocs.put("San Bernardino", new Location(34.064986, -117.29201));
		curveLocs.put("USC", new Location(34.0192, -118.286));
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
			curveLocs.put(city.toString(), city.location());
		}
		
		File compareDir = new File(hazardJobDir, "ucerf-comparisons");
		File u3Dir;
		if (u3SupraMinMag)
			u3Dir = new File(compareDir, "ucerf3-supra");
		else
			u3Dir = new File(compareDir, "ucerf3-m"+(float)minMag);
		Map<String, Map<Location, DiscretizedFunc>> u3CurvesMap = new HashMap<>();
		Map<String, Map<Location, DiscretizedFunc>> u3FullCurvesMap = new HashMap<>();
		
		AbstractERF u3ERF = null;
		File u3SolFile = new File(u3Dir, "ucerf3_sol_filtered.zip");
		FaultSystemSolution u3Sol = null;
		
		for (int d=0; d<jobDirs.size(); d++) {
			File jobDir = jobDirs.get(d);
			System.out.println("Processing "+(d+1)+"/"+jobDirs.size()+": "+jobDir);
			
			String imtLabel = "PGA (g)";
			String imtRatioLabel = "PGA";
			String imt = PGA_Param.NAME;
			String imtFileLabel;
			double period = Double.NaN;
			// try to detect other IMT
			if (jobDir.getName().toLowerCase().contains("-sa-")) {
				imt = SA_Param.NAME;
				String jobName = jobDir.getName().toLowerCase();
				String periodStr = jobName.substring(jobName.indexOf("-sa-")+4);
				periodStr = periodStr.substring(0, periodStr.indexOf("s"));
				period = Double.parseDouble(periodStr);
				imtLabel = (float)period+"s SA (g)";
				imtRatioLabel = (float)period+"s SA";
				System.out.println("Detected IMT: "+imtLabel);
				imtFileLabel = "sa_"+(float)period+"s";
			} else {
				System.out.println("Assuming IMT: "+imtLabel);
				imtFileLabel = "pga";
			}
			
			double stdDev = -1;
			if (jobDir.getName().contains("-stdDev")) {
				String jobName = jobDir.getName();
				String stdDevStr = jobName.substring(jobName.indexOf("-stdDev")+7);
				if (stdDevStr.contains("-"))
					stdDevStr = stdDevStr.substring(0, stdDevStr.indexOf("-"));
				stdDev = Double.parseDouble(stdDevStr);
				
				imtLabel += ", σ="+(float)stdDev;
				imtFileLabel += "_sigma"+(float)stdDev;
			}
			
			if (jobDir.getName().contains("-gmpe")) {
				String jobName = jobDir.getName();
				String gmpeStr = jobName.substring(jobName.indexOf("-gmpe")+5);
				if (gmpeStr.contains("-"))
					gmpeStr = gmpeStr.substring(0, gmpeStr.indexOf("-"));
				String gmpeName = gmpeStr;
				
				imtLabel += ", "+gmpeName;
				imtFileLabel += "_gmpe"+gmpeName;
			}
			
			List<String> lines = new ArrayList<>();
			
			// header
			lines.add("# Hazard Comparisons");
			lines.add("");
			lines.add("*IMT: "+imtLabel+"*");
			lines.add("");
			
			String outputName = "hazard_"+imtFileLabel;
			if (jobDir.getName().contains("sectArea")) {
				String str = jobDir.getName();
				str = str.substring(str.indexOf("sectArea")+"sectArea".length());
				str = str.substring(0, str.indexOf("-"));
				double minFractForInclusion = Double.parseDouble(str);
				if ((float)minFractForInclusion != (float)RSQSimBBP_Config.MIN_SUB_SECT_FRACT)
					outputName += "_sectArea"+str;
				lines.add("*Subsections participate in a rupture if at least "+(float)(minFractForInclusion*100d)+" % of its area ruptures*");
				lines.add("");
			}
			lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
			lines.add("");
			
			File catHazardOutDir = new File(catOutDir, outputName);
			Preconditions.checkState(catHazardOutDir.exists() || catHazardOutDir.mkdir());
			File resourcesDir = new File(catHazardOutDir, "resources");
			Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
			
			
			
			int tocIndex = lines.size();
			String topLink = "*[(top)](#table-of-contents)*";
			
			Map<Location, DiscretizedFunc> u3Curves;
			Map<Location, DiscretizedFunc> u3FullCurves;
			if (u3CurvesMap.containsKey(imtLabel)) {
				u3Curves = u3CurvesMap.get(imtLabel);
				u3FullCurves = u3FullCurvesMap.get(imtLabel);
			} else {
				File myU3Dir;
				if (stdDev >= 0)
					myU3Dir = new File(u3Dir.getAbsolutePath()+"-stdDev"+(float)stdDev);
				else
					myU3Dir = u3Dir;
				if (jobDir.getName().contains("xPoints")) {
					String jobName = jobDir.getName();
					jobName = jobName.substring(0, jobName.indexOf("xPoints"));
					jobName = jobName.substring(jobName.lastIndexOf("-")+1);
					myU3Dir = new File(myU3Dir.getAbsolutePath()+"-"+jobName+"xPoints");
				}
				if (jobDir.getName().contains("-maxDist")) {
					String jobName = jobDir.getName();
					jobName = jobName.substring(jobName.indexOf("-maxDist")+1);
					if (jobName.contains("-"))
						jobName = jobName.substring(0, jobName.indexOf("-"));
					myU3Dir = new File(myU3Dir.getAbsolutePath()+"-"+jobName);
				}
				File u3File = new File(myU3Dir, HazardMapCompareScriptGen.getCurveDirName(imt, period)+"/imrs1.bin");
				System.out.print("Loading UCERF3 from: "+u3File.getAbsolutePath()+" ...");
				BinaryHazardCurveReader u3Reader = new BinaryHazardCurveReader(u3File.getAbsolutePath());
				System.out.println("DONE");
				u3Curves = asLightFixedXMap(u3Reader.getCurveMap());
				u3FullCurves = null;
				if (plotCurves && u3FullFiles.containsKey(period) && u3FullFiles.get(period).exists()) {
					BinaryHazardCurveReader u3FullReader = new BinaryHazardCurveReader(u3FullFiles.get(period).getAbsolutePath());
					u3FullCurves = asLightFixedXMap(u3FullReader.getCurveMap());
				}
				u3CurvesMap.put(imtLabel, u3Curves);
				u3FullCurvesMap.put(imtLabel, u3FullCurves);
			}
			
			File curveFile = new File(jobDir, HazardMapCompareScriptGen.getCurveDirName(imt, period)+"/imrs1.bin");
			BinaryHazardCurveReader curveReader = new BinaryHazardCurveReader(curveFile.getAbsolutePath());
			File rsSolFile = new File(jobDir, "rsqsim_solution.zip");
			FaultSystemSolution rsSol = null;
			
			Map<Location, DiscretizedFunc> curves = asLightFixedXMap(curveReader.getCurveMap());
			
			CPT hazardCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
			hazardCPT = hazardCPT.rescale(0d, 1.2d);
			if (imt.equalsIgnoreCase(SA_Param.NAME)) {
				if (period < 1d)
					hazardCPT = hazardCPT.rescale(0d, 2d);
				else if (period < 2d)
					hazardCPT = hazardCPT.rescale(0d, 1.2);
				else if (period < 5d)
					hazardCPT = hazardCPT.rescale(0d, 1d);
				else
					hazardCPT = hazardCPT.rescale(0d, 0.6);
			}
			hazardCPT.setNanColor(Color.WHITE);
			
			if (plotMaps) {
				System.out.println("Plotting maps");
				for (int mapRP : mapRPs) {
					boolean isProbAtIML = false;
					double level = 1d/(double)mapRP;
					
					String durationLabel = mapRP+"yr";
					
					GriddedGeoDataSet rsqsimData = loadFromBinary(gridReg, curves, isProbAtIML, level);
					GriddedGeoDataSet u3Data = loadFromBinary(gridReg, u3Curves, isProbAtIML, level);
					
					File csvFile = new File(resourcesDir, "map_"+durationLabel+".csv.gz");
					CSVFile<String> csv = new CSVFile<>(true);
					csv.addLine("Index", "Longitude", "Latitude", catalogName, "UCERF3");
					for (int i=0; i<gridReg.getNodeCount(); i++) {
						Location loc = rsqsimData.getLocation(i);
						csv.addLine(i+"", (float)loc.getLongitude()+"", (float)loc.getLatitude()+"",
								(float)rsqsimData.get(i)+"", (float)u3Data.get(i)+"");
					}
					System.out.println("Writing CSV to "+csvFile.getAbsolutePath());
					OutputStream gzFileOut = new GZIPOutputStream(new FileOutputStream(csvFile));
					csv.writeToStream(gzFileOut);
					
					plotMaps(resourcesDir, "map_"+durationLabel+"_"+catalogFileName, rsqsimData, region,
							(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), catalogName+", "+durationLabel+", "+imtLabel, hazardCPT, false);
					plotMaps(resourcesDir, "map_"+durationLabel+"_u3", u3Data, region,
							(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), "UCERF3, "+durationLabel+", "+imtLabel, hazardCPT, false);
					
					GriddedGeoDataSet ratioData = new GriddedGeoDataSet(gridReg, false);
					for (int i=0; i<gridReg.getNodeCount(); i++)
						ratioData.set(i, rsqsimData.get(i) / u3Data.get(i));
					ratioData.log();
					
					CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance();
					ratioCPT.setNanColor(Color.WHITE);
					ratioCPT = ratioCPT.rescale(-0.2d, 0.2d);
					plotMaps(resourcesDir, "map_"+durationLabel+"_ratio_log_tight", ratioData, region, (double)ratioCPT.getMinValue(),
							(double)ratioCPT.getMaxValue(), "Ln("+catalogName+" / UCERF3), "+durationLabel+", "+imtRatioLabel, ratioCPT, false);
					ratioCPT = ratioCPT.rescale(-0.5d, 0.5d);
					plotMaps(resourcesDir, "map_"+durationLabel+"_ratio_log", ratioData, region, (double)ratioCPT.getMinValue(),
							(double)ratioCPT.getMaxValue(), "Ln("+catalogName+" / UCERF3), "+durationLabel+", "+imtRatioLabel, ratioCPT, false);
				}
			}
			
			if (plotHist) {
				System.out.println("Plotting hists");
				plotHists(u3Curves, curves, catalogName, gridReg, histRPs, histHighlightIndex, resourcesDir, true, imtLabel);
				plotHists(u3Curves, curves, catalogName, gridReg, histRPs, histHighlightIndex, resourcesDir, false, imtLabel);
				plotNEHRP_Hists(u3Curves, curves, catalogName, gridReg, nehrpRPs, resourcesDir, imtLabel);
//				plotMeanStdDevTrend(500, 30000, u3Curves, rsqsimCurves, gridReg, resourcesDir, imtLabel);
				plotMeanStdDevTrend(1e-6, u3Curves, curves, catalogName, gridReg, resourcesDir, imtLabel);
			}
			
			if (plotCurves) {
				System.out.println("Plotting curves");
				File curveDir = new File(resourcesDir, "curves");
				Preconditions.checkState(curveDir.exists() || curveDir.mkdir());
				
				for (String siteName : curveLocs.keySet()) {
					Location siteLoc = curveLocs.get(siteName);
//					Location gridLoc = rsqsimData.getLocation(rsqsimData.indexOf(siteLoc));
					Location gridLoc = gridReg.locationForIndex(gridReg.indexForLocation(siteLoc));
					double dist = LocationUtils.horzDistance(siteLoc, gridLoc);
					System.out.println(siteName+" is "+(float)dist+" km away from nearest grid point");
					
					DiscretizedFunc u3Curve = u3Curves.get(gridLoc);
					DiscretizedFunc rsqsimCurve = curves.get(gridLoc);
					DiscretizedFunc u3FullCurve = null;
					if (u3FullCurves != null)
						u3FullCurve = u3FullCurves.get(gridLoc);
					
					plotCurves(u3Curve, u3FullCurve, rsqsimCurve, catalogName, curveDir, siteName, imtLabel, mapRPs);
				}
			}
			
			if (plotDisagg) {
				System.out.println("Plotting disagg");
				System.out.println("Building RSQSim ERF");
				rsSol = FaultSystemIO.loadSol(rsSolFile);
				FaultSystemSolutionERF rsERF = buildERF(rsSol);
				if (u3ERF == null) {
					System.out.println("Building UCERF3 ERF");
					u3Sol = FaultSystemIO.loadSol(u3SolFile);
					u3ERF = buildERF(u3Sol);
				}
				
				File disaggDir = new File(resourcesDir, "disagg");
				Preconditions.checkState(disaggDir.exists() || disaggDir.mkdir());
				
				double disaggMinMag;
				if (u3SupraMinMag)
					disaggMinMag = 6;
				else
					disaggMinMag = minMag;
				
				calcDisagg(u3ERF, rsERF, curveLocs, catalogName, disaggDir, imt, period, disaggRPs, disaggMinMag);
			}
			
			if (plotCurveProfile) {
				Location startLoc = new Location(40, -120);
				plotCurveProfile(u3Curves, curves, gridReg, imtLabel, resourcesDir, "curve_profile", startLoc, 15, 0d, 0.02);
			}
			
			waitOnFutures(false);
			
			if (plotMaps || new File(resourcesDir, "map_"+mapRPs[0]+"yr_"+catalogFileName+".png").exists()) {
				lines.add("## Hazard Maps");
				lines.add(topLink); lines.add("");
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine("Return Period", catalogName, "UCERF3", "Ratio", "Tight Ratio");
				for (int mapRP : mapRPs) {
					String durationLabel = mapRP+"yr";
					table.initNewLine();
					table.addColumn("**"+mapRP+" yr**");
					File rsPlot = new File(resourcesDir, "map_"+mapRP+"yr_"+catalogFileName+".png");
					Preconditions.checkState(rsPlot.exists());
					table.addColumn("![Catalog Map]("+resourcesDir.getName()+"/"+rsPlot.getName()+")");
					File u3Plot = new File(resourcesDir, "map_"+mapRP+"yr_u3.png");
					Preconditions.checkState(u3Plot.exists());
					table.addColumn("![UCERF3 Map]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
					File ratioPlot = new File(resourcesDir, "map_"+mapRP+"yr_ratio_log.png");
					Preconditions.checkState(ratioPlot.exists());
					table.addColumn("![Ratio Map]("+resourcesDir.getName()+"/"+ratioPlot.getName()+")");
					File tightRatioPlot = new File(resourcesDir, "map_"+mapRP+"yr_ratio_log_tight.png");
					Preconditions.checkState(tightRatioPlot.exists());
					table.addColumn("![Tight Ratio Map]("+resourcesDir.getName()+"/"+tightRatioPlot.getName()+")");
					
					table.finalizeLine();
				}
				lines.add("");
				lines.addAll(table.build());
				lines.add("");
			}
			
			if (plotHist || new File(resourcesDir, "hist_0d.png").exists()) {
				lines.add("## Histograms");
				lines.add("");
				lines.add("## 0-D Histogram");
				lines.add(topLink); lines.add("");
				lines.add("![0-D Hist]("+resourcesDir.getName()+"/hist_0d.png)");
				lines.add("");
				lines.add("## 1-D Histogram");
				lines.add(topLink); lines.add("");
				lines.add("![1-D Hist]("+resourcesDir.getName()+"/hist_1d.png)");
				lines.add("");
				lines.add("## 2-D Histograms");
				lines.add(topLink); lines.add("");
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine("Return Period", "2-D Histogram");
				for (int rp : histRPs) {
					table.initNewLine();
					table.addColumn("**"+rp+" yr**");
					table.addColumn("![2-D Hist]("+resourcesDir.getName()+"/hist_2d_"+rp+"yr.png)");
					table.finalizeLine();
				}
				lines.addAll(table.build());
				lines.add("");
			}
			
			if (plotCurves || new File(resourcesDir, "curves").exists()) {
				lines.add("## Hazard Curves");
				lines.add(topLink); lines.add("");
				
				List<String> siteNames = new ArrayList<>(curveLocs.keySet());
				Collections.sort(siteNames);
				
				File curveDir = new File(resourcesDir, "curves");
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				for (String siteName : siteNames) {
					File plot = new File(curveDir, "curves_"+siteName.replaceAll(" ", "_")+".png");
					Preconditions.checkState(plot.exists());
					table.addColumn("!["+siteName+"]("+resourcesDir.getName()+"/curves/"+plot.getName()+")");
				}
				table.finalizeLine();
				table.wrap(4, 0);
				
				lines.addAll(table.build());
				lines.add("");
			}
			
			if (plotDisagg || new File(resourcesDir, "disagg").exists()) {
				lines.add("## Hazard Disaggregations");
				lines.add(topLink); lines.add("");
				
				List<String> siteNames = new ArrayList<>(curveLocs.keySet());
				Collections.sort(siteNames);
				
				File disaggDir = new File(resourcesDir, "disagg");
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				for (String siteName : siteNames) {
					for (int rp : disaggRPs) {
						File plot = new File(disaggDir, "disagg_"+siteName.replaceAll(" ", "_")+"_"+rp+"yr_combined.png");
//						Preconditions.checkState(plot.exists());
						table.addLine("**"+siteName+", "+rp+" yr**", "!["+siteName+"]("+resourcesDir.getName()+"/disagg/"+plot.getName()+")");
					}
				}
				
				lines.addAll(table.build());
				lines.add("");
			}
			
			// add TOC
			lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
			lines.add(tocIndex, "## Table Of Contents");
			
			// write markdown
			MarkdownUtils.writeReadmeAndHTML(lines, catHazardOutDir);
		}
		
		catalog.writeMarkdownSummary(catOutDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
		
		waitOnFutures(true);
	}
	
	private static ArrayDeque<ScalarIMR> gmpeDeque;
	
	private synchronized static ScalarIMR checkOutGMPE(String imt, double period) {
		if (gmpeDeque == null)
			gmpeDeque = new ArrayDeque<>();
		ScalarIMR gmpe;
		if (gmpeDeque.isEmpty()) {
			// build a new one
			gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		} else {
			gmpe = gmpeDeque.pop();
		}
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(imt);
		if (imt.equals(SA_Param.NAME))
			SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
		return gmpe;
	}
	
	private synchronized void checkInGMPE(ScalarIMR gmpe) {
		gmpeDeque.push(gmpe);
	}
	
	private static FaultSystemSolutionERF buildERF(FaultSystemSolution sol) {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		return erf;
	}
	
	private static Map<Location, DiscretizedFunc> asLightFixedXMap(Map<Location, ? extends DiscretizedFunc> orig) {
		Map<Location, DiscretizedFunc> ret = new HashMap<>();
		double[] xVals = null;
		
		for (Location loc : orig.keySet()) {
			DiscretizedFunc curve = orig.get(loc);
			
			if (xVals == null) {
				xVals = new double[curve.size()];
				for (int i=0; i<xVals.length; i++)
					xVals[i] = curve.getX(i);
			}
			
			double[] yVals = new double[curve.size()];
			for (int i=0; i<yVals.length; i++)
				yVals[i] = curve.getY(i);
			
			ret.put(loc, new LightFixedXFunc(xVals, yVals));
		}
		
		return ret;
	}
	
	static GriddedGeoDataSet loadFromBinary(GriddedRegion gridReg, Map<Location, ? extends DiscretizedFunc> curves,
			boolean isProbAtIML, double level) {
		GriddedGeoDataSet data = new GriddedGeoDataSet(gridReg, false);
		
		for (Location loc : curves.keySet()) {
			DiscretizedFunc curve = curves.get(loc);
			double val = HazardDataSetLoader.getCurveVal(curve, isProbAtIML, level);
			data.set(loc, val);
		}
		
		return data;
	}
	
	static void plotMaps(File outputDir, String prefix, GriddedGeoDataSet data, Region region,
			Double customMin, Double customMax, String label, CPT cpt, boolean rescaleCPT)
					throws GMT_MapException, IOException {
		
		System.out.println("Creating map instance...");
		GMT_Map map = new GMT_Map(region, data, data.getRegion().getLatSpacing(), cpt);
		
		map.setCustomLabel(label);
		map.setLabelTickSize(20);
		map.setLabelSize(24);
//		map.setTopoResolution(TopographicSlopeFile.US_SIX);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setCoast(new CoastAttributes(Color.GRAY, 1));
		map.setLogPlot(false);
		map.setDpi(300);
		map.setCustomScaleMin(customMin);
		map.setCustomScaleMax(customMax);
		map.setRescaleCPT(rescaleCPT);
		map.setBlackBackground(false);
//		map.setPDFFileName(null);
		map.setJPGFileName(null);
		
		System.out.println("Making map...");
		
		Runnable run = new Runnable() {
			
			@Override
			public void run() {
				try {
					FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		};
		
		if (map_parallel)
			submitRunnable(run, 1000);
		else
			run.run();
	}
	
	private static synchronized void submitRunnable(Runnable run, long sleepMillis) {
		if (exec == null) {
			exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
			futures = new LinkedList<>();
		}
		futures.add(exec.submit(run));
		if (sleepMillis > 0) {
			try {
				Thread.sleep(sleepMillis);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	static void waitOnFutures(boolean shutdown) {
		if (futures == null)
			return;
		System.out.println("Waiting on "+futures.size()+" maps...");
		int count = 0;
		while (!futures.isEmpty()) {
			Future<?> future = futures.removeFirst();
			try {
				future.get();
			} catch (Exception e) {
				e.printStackTrace();
			}
			count++;
			if (count % 10 == 0)
				System.out.println("\t"+count+" Futures done, "+futures.size()+" left");
		}
		System.out.println("DONE");
		if (shutdown)
			exec.shutdown();
	}
	
	private static ExecutorService exec;
	private static LinkedList<Future<?>> futures;
	
	public static void plotHists(Map<Location, ? extends DiscretizedFunc> ucerf3Curves,
			Map<Location, ? extends DiscretizedFunc> rsqsimCurves, String catalogName, GriddedRegion gridReg,
			int[] returnPeriods, int hightlightIndex, File outputDir, boolean log, String imtLabel) throws IOException {
		
		PlotPreferences plotPrefs = buildPlotPrefs();
		
		List<List<Double>> logRatioValsList = new ArrayList<>();
		
		for (int rp : returnPeriods) {
			GriddedGeoDataSet ucerf3 = loadFromBinary(gridReg, ucerf3Curves, false, 1d/(double)rp);
			GriddedGeoDataSet rsqsim = loadFromBinary(gridReg, rsqsimCurves, false, 1d/(double)rp);
			
			// min non zero
			double minZ = Double.POSITIVE_INFINITY;
			for (int i=0; i<ucerf3.size(); i++) {
				if (ucerf3.get(i) > 0)
					minZ = Math.min(minZ, ucerf3.get(i));
				if (rsqsim.get(i) > 0)
					minZ = Math.min(minZ, rsqsim.get(i));
			}
			double maxZ = Math.max(ucerf3.getMaxZ(), rsqsim.getMaxZ());
			
			double delta = 0.02;
			if (log) {
				rsqsim = rsqsim.copy();
				ucerf3 = ucerf3.copy();
				rsqsim.log();
				ucerf3.log();
				minZ = Math.log(minZ);
				maxZ = Math.log(maxZ);
				delta = 0.05;
			}
			
			HistogramFunction sampleHist = HistogramFunction.getEncompassingHistogram(minZ, maxZ, delta);
			while (sampleHist.size() < 30) {
				delta /= 2;
				sampleHist = HistogramFunction.getEncompassingHistogram(minZ, maxZ, delta);
			}
			
			EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(sampleHist.size(), sampleHist.size(), minZ, minZ, delta);
			
			List<Double> logRatioVals = new ArrayList<>();
			logRatioValsList.add(logRatioVals);
			
			for (int i=0; i<ucerf3.size(); i++) {
				double u3 = ucerf3.get(i);
				double rs = rsqsim.get(i);
				if (log && (Double.isNaN(u3) || Double.isNaN(rs)))
					continue;
				if (!log && (u3 == 0 || rs == 0))
					continue;
				int xInd = sampleHist.getClosestXIndex(u3);
				int yInd = sampleHist.getClosestXIndex(rs);
				
				xyz.set(xInd, yInd, xyz.get(xInd, yInd) + 1d);
				
				if (log)
					// already in log space, subract
					logRatioVals.add(rs - u3);
			}
			
			xyz.log10();
			
			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
			cpt = cpt.rescale(Math.log10(1d), xyz.getMaxZ());
			cpt.setNanColor(Color.WHITE);
			cpt.setBelowMinColor(Color.WHITE);
			
			String xAxisLabel = "UCERF3 1/"+rp+"yr, "+imtLabel;
			String yAxisLabel = catalogName+" 1/"+rp+"yr, "+imtLabel;
			if (log) {
				xAxisLabel = "Ln "+xAxisLabel;
				yAxisLabel = "Ln "+yAxisLabel;
			}
			
			XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, "Hazard Histogram", xAxisLabel, yAxisLabel, "Log10(Number)");
			
			Range plotRange = new Range(minZ - 0.5*delta, maxZ + 0.5*delta);
			
			List<XY_DataSet> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
			oneToOne.set(plotRange.getLowerBound(), plotRange.getLowerBound());
			oneToOne.set(plotRange.getUpperBound(), plotRange.getUpperBound());
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			xyzSpec.setXYElems(funcs);
			xyzSpec.setXYChars(chars);
			
			XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
			xyzGP.drawPlot(xyzSpec, false, false, plotRange, plotRange);
			if (!log) {
				TickUnits tus = new TickUnits();
				TickUnit tu = new NumberTickUnit(0.25);
				tus.add(tu);
				xyzGP.getXAxis().setStandardTickUnits(tus);
				xyzGP.getYAxis().setStandardTickUnits(tus);
			}
			// write plot
			xyzGP.getChartPanel().setSize(850, 900);
			String prefix = "hist_2d_"+rp+"yr";
			xyzGP.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			xyzGP.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		}
		
		CPT rpCPT = getRPlogCPT(returnPeriods);
		
		if (log) {
			boolean logY = true;
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			double highlightMean = 0;
			double highlightStdDev = 0;
			
			double histMinNonZero = Double.POSITIVE_INFINITY;
			double histMax = 0;
			
			for (int i=0; i<returnPeriods.length; i++) {
				int rp = returnPeriods[i];
				Color color = rpCPT.getColor((float)Math.log10(rp));
				
				double[] logRatioArray = Doubles.toArray(logRatioValsList.get(i));
				HistogramFunction oneDHist = HistogramFunction.getEncompassingHistogram(
						StatUtils.min(logRatioArray), StatUtils.max(logRatioArray), 0.05);
				for (double val : logRatioArray)
					oneDHist.add(val, 1d);
				double mean = StatUtils.mean(logRatioArray);
				double stdDev = Math.sqrt(StatUtils.variance(logRatioArray));
				
				oneDHist.normalizeBySumOfY_Vals();
				oneDHist.setName(rp+"yr");
				
				for (Point2D pt : oneDHist)
					if (pt.getY() > 0)
						histMinNonZero = Math.min(histMinNonZero, pt.getY());
				histMax = Math.max(histMax, oneDHist.getMaxY());
				
				if (i == hightlightIndex) {
					funcs.add(0, oneDHist);
					chars.add(0, new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
					highlightMean = mean;
					highlightStdDev = stdDev;
				} else {
					funcs.add(oneDHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
				}
			}
			
			Range yRange = null;
			if (logY) {
				yRange = new Range(histMinNonZero*0.7, histMax*1.5);
			}
			
			if (hightlightIndex >= 0) {
				DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
				meanLine.set(highlightMean, 0d);
				if (logY) {
					meanLine.set(highlightMean, yRange.getLowerBound());
					meanLine.set(highlightMean, yRange.getUpperBound());
				}
				meanLine.set(highlightMean, histMax);
				meanLine.setName(returnPeriods[hightlightIndex]+"yr Mean: "+fourDigits.format(highlightMean)
					+", Std Dev: "+fourDigits.format(highlightStdDev));
				
				funcs.add(meanLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Hazard 1D Histogram", "Ln("+catalogName+"/UCERF3), "+imtLabel, "Fraction");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			gp.drawGraphPanel(spec, false, logY, null, yRange);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(new File(outputDir, "hist_1d.png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, "hist_1d.pdf").getAbsolutePath());
		}
	}
	
	private static final DecimalFormat fourDigits = new DecimalFormat("0.0000");
	
	static void plotMeanStdDevTrend(double minProb, Map<Location, ? extends DiscretizedFunc> ucerf3Curves,
			Map<Location, ? extends DiscretizedFunc> rsqsimCurves, String catalogName,
			GriddedRegion gridReg, File outputDir, String imtLabel) throws IOException {
		double logMinProb = Math.log10(minProb);
		double logMaxProb = 0;
		int numProbs = (int)(logMaxProb - logMinProb)*5 + 1;
		
		EvenlyDiscretizedFunc meanFunc = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
		EvenlyDiscretizedFunc stdDevFunc = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
		EvenlyDiscretizedFunc medianFunc = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
		EvenlyDiscretizedFunc meanAbsDiffFunc = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
		
//		double[] pga_thresholds = {0.0, 0.1, 0.2, 0.4, 0.8, 1.2};
		double[] pga_thresholds = null;
		CPT pgaWeightCPT = null;
		EvenlyDiscretizedFunc[] pgaWeightFuncs = null;
		EvenlyDiscretizedFunc[] pgaWeightStdDevFuncs = null;
		if (pga_thresholds != null && pga_thresholds.length > 0) {
			pgaWeightCPT = new CPT(pga_thresholds[0], pga_thresholds[pga_thresholds.length-1], Color.GRAY, Color.BLACK);
			pgaWeightFuncs = new EvenlyDiscretizedFunc[pga_thresholds.length];
			for (int p=0; p<pga_thresholds.length; p++)
				pgaWeightFuncs[p] = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
			pgaWeightStdDevFuncs = new EvenlyDiscretizedFunc[pga_thresholds.length];
			for (int p=0; p<pga_thresholds.length; p++)
				pgaWeightStdDevFuncs[p] = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
		}
		
		EvenlyDiscretizedFunc nehrpMeanFunc = new EvenlyDiscretizedFunc(logMinProb, logMaxProb, numProbs);
		
		HashSet<Integer> nehrpGridIndexes = new HashSet<>();
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
			Location loc = city.getSite().getLocation();
			nehrpGridIndexes.add(gridReg.indexForLocation(loc));
		}
		
		for (int i=0; i<numProbs; i++) {
			double prob = Math.pow(10, meanFunc.getX(i));
			GriddedGeoDataSet ucerf3 = loadFromBinary(gridReg, ucerf3Curves, false, prob);
			GriddedGeoDataSet rsqsim = loadFromBinary(gridReg, rsqsimCurves, false, prob);
			
			List<Double> ratioVals = new ArrayList<>();
			
			List<List<Double>> pgaWeights = null;
			if (pga_thresholds != null && pga_thresholds.length > 0) {
				pgaWeights = new ArrayList<>();
				for (int p=0; p<pga_thresholds.length; p++)
					pgaWeights.add(new ArrayList<>());
			}
			
			List<Double> nehrpRatios = new ArrayList<>();
			
			for (int j=0; j<ucerf3.size(); j++) {
				double u3 = ucerf3.get(j);
				double rs = rsqsim.get(j);
				if (u3 == 0 || rs == 0 || Double.isNaN(u3) || Double.isNaN(rs))
					continue;
				
				double ratio = Math.log(rs/u3);
				Preconditions.checkState(Double.isFinite(ratio), "Bad ratio: Ln(%s/%s) = %s", rs, u3, ratio);
				
				if (nehrpGridIndexes.contains(j))
					nehrpRatios.add(ratio);
				
				ratioVals.add(Math.log(rs/u3));
				if (pga_thresholds != null && pga_thresholds.length > 0) {
					for (int p=0; p<pga_thresholds.length; p++) {
						double weight = u3 - pga_thresholds[p];
						if (weight < 0)
							weight = 0;
						pgaWeights.get(p).add(weight);
					}
				}
			}
			
			double[] ratioArray = Doubles.toArray(ratioVals);
			
			// calc MAD
			double meanAbsDiff, median;
			if (ratioVals.isEmpty()) {
				meanAbsDiff = Double.NaN;
				median = Double.NaN;
			} else {
				median = DataUtils.median(ratioArray);
				meanAbsDiff = 0d;
				for (int j=0; j<ratioArray.length; j++)
					meanAbsDiff += Math.abs(ratioArray[j]);
				meanAbsDiff /= ratioArray.length;
//				int onePercent = ratioArray.length / 100;
//				for (int j=0; j<ratioArray.length; j++) {
//					if (j % onePercent == 0)
//						System.out.println("j="+j+"/"+ratioArray.length);
//					for (int k=0; k<ratioArray.length; k++)
//						meanAbsDiff += Math.abs(ratioArray[j] - ratioArray[k]);
//				}
//				meanAbsDiff /= ratioArray.length*ratioArray.length;
			}
			
			double mean = StatUtils.mean(ratioArray);
			double stdDev = Math.sqrt(StatUtils.variance(ratioArray));
			
			meanFunc.set(i, mean);
			stdDevFunc.set(i, stdDev);
			medianFunc.set(i, median);
			meanAbsDiffFunc.set(i, meanAbsDiff);
			
			double nehrpMean = StatUtils.mean(Doubles.toArray(nehrpRatios));
			nehrpMeanFunc.set(i, nehrpMean);
			
			if (pga_thresholds != null && pga_thresholds.length > 0) {
				for (int p=0; p<pga_thresholds.length; p++) {
					ArbDiscrEmpiricalDistFunc func = new ArbDiscrEmpiricalDistFunc();
					for (int j=0; j<ratioVals.size(); j++)
						func.set(ratioVals.get(j), pgaWeights.get(p).get(j));
					pgaWeightFuncs[p].set(i, func.getMean());
					pgaWeightStdDevFuncs[p].set(i, func.getStdDev());
				}
			}
			
//			System.out.println(rp+" "+mean+" "+stdDev+" "+ratioVals.size());
		}
		
		List<DiscretizedFunc> logFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (pga_thresholds != null && pga_thresholds.length > 0) {
			for (int p=0; p<pga_thresholds.length; p++) {
				logFuncs.add(pgaWeightFuncs[p]);
				String name;
				if (p == 0 || p == pga_thresholds.length-1)
					name = "PSA_0="+(float)pga_thresholds[p];
				else
					name = (float)pga_thresholds[p]+"";
				pgaWeightFuncs[p].setName(name);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f+2f*(float)pga_thresholds[p],
						pgaWeightCPT.getColor((float)pga_thresholds[p])));
			}
			
			for (int p=0; p<pga_thresholds.length; p++) {
				logFuncs.add(pgaWeightStdDevFuncs[p]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f+2f*(float)pga_thresholds[p],
						pgaWeightCPT.getColor((float)pga_thresholds[p])));
			}
		}
		
		logFuncs.add(nehrpMeanFunc);
		nehrpMeanFunc.setName("NEHRP City Mean");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
		logFuncs.add(meanFunc);
		meanFunc.setName("Mean");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		logFuncs.add(stdDevFunc);
		stdDevFunc.setName("StdDev");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		
		logFuncs.add(medianFunc);
		medianFunc.setName("Median");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.MAGENTA));
		
		logFuncs.add(meanAbsDiffFunc);
		meanAbsDiffFunc.setName("MeanAbsDiff");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.ORANGE.darker()));
		
		List<DiscretizedFunc> linearFuncs = new ArrayList<>();
		for (DiscretizedFunc logFunc : logFuncs)
			linearFuncs.add(toLinearX(logFunc));
		
		for (int i=linearFuncs.size(); --i>=0;) {
			if (linearFuncs.get(i).size() == 0) {
				linearFuncs.remove(i);
				chars.remove(i);
			}
		}
		
		PlotSpec spec = new PlotSpec(linearFuncs, chars, "Mean/StdDev Trend", "Annual Probability", "Ln("+catalogName+"/UCERF3), "+imtLabel);
		spec.setLegendVisible(true);
		
		PlotPreferences plotPrefs = buildPlotPrefs();
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		gp.drawGraphPanel(spec, true, false, new Range(minProb, 1d), new Range(-0.6, 0.6));
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsTXT(new File(outputDir, "hist_0d.txt").getAbsolutePath());
		gp.saveAsPNG(new File(outputDir, "hist_0d.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "hist_0d.pdf").getAbsolutePath());
	}
	
	private static DiscretizedFunc toLinearX(DiscretizedFunc logXFunc) {
		ArbitrarilyDiscretizedFunc linearFunc = new ArbitrarilyDiscretizedFunc();
		linearFunc.setName(logXFunc.getName());
		
		for (int i=0; i<logXFunc.size(); i++)
			if (Double.isFinite(logXFunc.getY(i)))
				linearFunc.set(Math.pow(10d, logXFunc.getX(i)), logXFunc.getY(i));
		
		return linearFunc;
	}
	
	public static CPT getRPlogCPT(int[] rps) {
		int minRP = Integer.MAX_VALUE;
		int maxRP = 0;
		
		for (int rp : rps) {
			if (rp < minRP)
				minRP = rp;
			if (rp > maxRP)
				maxRP = rp;
		}
		return new CPT(Math.log10(minRP), Math.log10(maxRP), Color.LIGHT_GRAY, Color.BLACK);
	}
	
	static void plotNEHRP_Hists(Map<Location, ? extends DiscretizedFunc> ucerf3Curves,
			Map<Location, ? extends DiscretizedFunc> rsqsimCurves, String catalogName, GriddedRegion gridReg,
			int[] returnPeriods, File outputDir, String imtLabel) throws IOException {
		HashSet<Integer> nehrpGridIndexes = new HashSet<>();
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
			Location loc = city.getSite().getLocation();
			nehrpGridIndexes.add(gridReg.indexForLocation(loc));
		}
		
		for (int rp : returnPeriods) {
			List<Double> ratioVals = new ArrayList<>();
			
			GriddedGeoDataSet ucerf3 = loadFromBinary(gridReg, ucerf3Curves, false, 1d/(double)rp);
			GriddedGeoDataSet rsqsim = loadFromBinary(gridReg, rsqsimCurves, false, 1d/(double)rp);
			
			for (int index : nehrpGridIndexes) {
				double u3 = ucerf3.get(index);
				double rs = rsqsim.get(index);
//				if (u3 == 0 || rs == 0 || Double.isNaN(u3) || Double.isNaN(rs))
//					continue;
				
				double ratio = Math.log(rs/u3);
				Preconditions.checkState(Double.isFinite(ratio), "Bad ratio: Ln(%s/%s) = %s", rs, u3, ratio);
				
				ratioVals.add(ratio);
			}
			
			double[] ratioArray = Doubles.toArray(ratioVals);
			
			double mean = StatUtils.mean(ratioArray);
			double stdDev = Math.sqrt(StatUtils.variance(ratioArray));
			
			double minX = -0.5;
			double maxX = 0.5;
			
			HistogramFunction hist = HistogramFunction.getEncompassingHistogram(minX, maxX, 0.05);
			for (double val : ratioArray)
				hist.add(hist.getClosestXIndex(val), 1d);
			hist.normalizeBySumOfY_Vals();
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.DARK_GRAY));
			
			Range xRange = new Range(minX, maxX);
			Range yRange = new Range(hist.getMinY()*0.9, hist.getMaxY()*1.1);
			
			PlotSpec spec = new PlotSpec(funcs, chars, "NEHRP City Residuals, "+rp+"yr",
					"Ln("+catalogName+"/UCERF3), "+imtLabel, "Fraction");
			
			double x = minX + 0.1;
			double y1 = yRange.getUpperBound()*0.9;
			double y2 = yRange.getUpperBound()*0.8;
			List<XYTextAnnotation> anns = new ArrayList<>();
			XYTextAnnotation meanAnn = new XYTextAnnotation("Mean: "+fourDigits.format(mean), x, y1);
			meanAnn.setTextAnchor(TextAnchor.TOP_LEFT);
			meanAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			anns.add(meanAnn);
			XYTextAnnotation stdDevAnn = new XYTextAnnotation("StdDev: "+fourDigits.format(stdDev), x, y2);
			stdDevAnn.setTextAnchor(TextAnchor.TOP_LEFT);
			stdDevAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			anns.add(stdDevAnn);
			spec.setPlotAnnotations(anns);
			
			PlotPreferences plotPrefs = buildPlotPrefs();
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			gp.getChartPanel().setSize(800, 600);
			String prefix = "nehrp_hist1d_"+rp;
			gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
		}
	}
	
	static void plotCurves(DiscretizedFunc ucerf3Compare, DiscretizedFunc ucerf3Full,
			DiscretizedFunc rsqsim, String catalogName, File outputDir, String siteName,
			String imt, int[] rps) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = new Range(1e-3, 1e1);
		Range yRange = new Range(1e-8, 1e0);
		
		if (ucerf3Full != null) {
			ucerf3Full.setName("UCERF3 Full");
			funcs.add(ucerf3Full);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		}
		
		ucerf3Compare.setName("UCERF3");
		funcs.add(ucerf3Compare);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		rsqsim.setName(catalogName);
		funcs.add(rsqsim);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		
		if (rps != null && rps.length > 0) {
			CPT rpCPT = getRPlogCPT(rps);
			for (int rp : rps) {
				Color color = rpCPT.getColor((float)Math.log10(rp));
				double probLevel = 1d/(double)rp;
				DiscretizedFunc probLine = new ArbitrarilyDiscretizedFunc();
				probLine.set(xRange.getLowerBound(), probLevel);
				probLine.set(xRange.getUpperBound(), probLevel);
				probLine.setName(rp+"yr");
				funcs.add(probLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, color));
			}
		}

		PlotSpec spec = new PlotSpec(funcs, chars, siteName+" Hazard Curves", imt, "Annual Probabilitiy");
		spec.setLegendVisible(true);
		
		PlotPreferences plotPrefs = buildPlotPrefs();
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		String prefix = "curves_"+siteName.replaceAll(" ", "_");
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(new File(outputDir, prefix).getAbsolutePath()+".png");
		gp.saveAsPDF(new File(outputDir, prefix).getAbsolutePath()+".pdf");
		gp.saveAsTXT(new File(outputDir, prefix).getAbsolutePath()+".txt");
	}
	
	static void calcDisagg(AbstractERF u3ERF, AbstractERF rsERF,
			Map<String, Location> curveLocs, String catalogName, File outputDir,
			String imt, double period, int[] rps, double minMag) throws IOException {
		DisaggCalculator calc = new DisaggCalculator(u3ERF, rsERF, catalogName, outputDir, imt, period, rps, minMag);
		
		for (String siteName : curveLocs.keySet())
			calc.calcForSite(siteName, curveLocs.get(siteName));
	}
	
	private static class DisaggCalculator {
		
		private final AbstractERF u3ERF;
		private final AbstractERF rsERF;
		private final String catalogName;
		private final File outputDir;
		private final String imt;
		private final double period;
		private final int[] rps;
		
		private final DiscretizedFunc xVals;
		private final DiscretizedFunc logXVals;
		
		private final ParameterList disaggParams;
		
		// disagg plot settings
		private final double minMag;
		private final int numMags;
		private final double deltaMag;

		private static final int numSourcesForDisag = 100;

		private static final boolean showSourceDistances = true;

		private static final double maxZAxis = Double.NaN;

		public DisaggCalculator(AbstractERF u3ERF, AbstractERF rsERF, String catalogName, File outputDir,
				String imt, double period, int[] rps, double minMag) {
			this.u3ERF = u3ERF;
			this.rsERF = rsERF;
			this.catalogName = catalogName;
			this.outputDir = outputDir;
			this.imt = imt;
			this.period = period;
			this.rps = rps;
			
			this.minMag = minMag;
			this.deltaMag = 0.2;
			this.numMags = (int)((8.6d - minMag)/deltaMag + 0.5);

			xVals = new IMT_Info().getDefaultHazardCurve(imt);
			logXVals = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : xVals)
				logXVals.set(Math.log(pt.getX()), 1d);
			
//			DisaggregationCalculator disaggCalc = new DisaggregationCalculator();
			disaggParams = DisaggregationCalculator.getDefaultParams();
		}
		
		public void calcForSite(String siteName, Location loc) {
			Runnable run = new Runnable() {
				
				@Override
				public void run() {
					try {
						HazardCurveCalculator curveCalc = new HazardCurveCalculator();
						DisaggregationCalculator disaggCalc = new DisaggregationCalculator();
						
						System.out.println("Disaggregating for "+siteName);
						
						ScalarIMR gmpe = checkOutGMPE(imt, period);
						
						Site site = new Site(loc);
						site.addParameterList(gmpe.getSiteParams());
						
						List<File> rsPNGs = new ArrayList<>();
						List<File> u3PNGs = new ArrayList<>();
						
						String sitePrefix = "disagg_"+siteName.replaceAll(" ", "_");
						
						for (boolean rs : new boolean[] {true, false}) {
							AbstractERF erf;
							if (rs) {
								erf = rsERF;
							} else {
								erf = u3ERF;
							}
							DiscretizedFunc logCurve = logXVals.deepClone();
							curveCalc.getHazardCurve(logCurve, site, gmpe, erf);
							DiscretizedFunc linearCurve = new ArbitrarilyDiscretizedFunc();
							for (int i=0; i<logCurve.size(); i++)
								linearCurve.set(xVals.getX(i), logCurve.getY(i));
							
							for (int rp : rps) {
								String prefix;
								if (rs) {
									prefix = sitePrefix+"_"+rp+"yr_"+catalogName.replaceAll(" ", "_").toLowerCase();
								} else {
									prefix = sitePrefix+"_"+rp+"yr_ucerf3";
								}
								
								double prob = 1d/(double)rp;
								double iml = HazardDataSetLoader.getCurveVal(linearCurve, false, prob); // iml at prob
								if (!Double.isFinite(iml)) {
									System.out.println("Couldn't get IML for "+siteName+", "+rp+"yr. Skipping disagg!");
									return;
								}

								System.out.println("Disaggregating for prob="+prob+", iml="+iml);
								disaggCalc.setMagRange(minMag, numMags, deltaMag);
								disaggCalc.setNumSourcestoShow(numSourcesForDisag);
								disaggCalc.setShowDistances(showSourceDistances);
								boolean success = disaggCalc.disaggregate(Math.log(iml), site, gmpe, erf, disaggParams);
								if (!success)
									throw new RuntimeException("Disagg calc failed (see errors above, if any).");
								disaggCalc.setMaxZAxisForPlot(maxZAxis);
								System.out.println("Done Disaggregating");
								String metadata = "temp metadata";

								System.out.println("Fetching plot...");
								String address = disaggCalc.getDisaggregationPlotUsingServlet(metadata);

								String meanModeText = disaggCalc.getMeanAndModeInfo();
								String binDataText = disaggCalc.getBinData();
								String sourceDataText = disaggCalc.getDisaggregationSourceInfo();

								File outputFile = new File(outputDir, prefix);

								String metadataText = "Custom disagg";

								File pdfFile = new File(outputFile.getAbsolutePath()+".pdf");
								File pngFile = new File(outputFile.getAbsolutePath()+".png");
								DisaggregationPlotViewerWindow.saveAsPDF(
										address+DisaggregationCalculator.DISAGGREGATION_PLOT_PDF_NAME,
										pdfFile.getAbsolutePath(), meanModeText, metadataText, binDataText, sourceDataText);
								FileUtils.downloadURL(address+DisaggregationCalculator.DISAGGREGATION_PLOT_PNG_NAME,
										pngFile);
								DisaggregationPlotViewerWindow.saveAsTXT(outputFile.getAbsolutePath()+".txt", meanModeText, metadataText,
										binDataText, sourceDataText);
								
								if (rs)
									rsPNGs.add(pngFile);
								else
									u3PNGs.add(pngFile);
							}
						}
						
						// now write combined PNG
						for (int i=0; i<rps.length; i++) {
							String prefix = sitePrefix+"_"+rps[i]+"yr_combined";
							BufferedImage rsImg = ImageIO.read(rsPNGs.get(i));
							BufferedImage u3Img = ImageIO.read(u3PNGs.get(i));
							int height = rsImg.getHeight();
							if (u3Img.getHeight() > rsImg.getHeight())
								height = u3Img.getHeight();
							int width = rsImg.getWidth()+u3Img.getWidth();
							
							BufferedImage comb = new BufferedImage(width, height, rsImg.getType());
							
							int xOffset = 0;
							for (int x=0; x<rsImg.getWidth(); x++)
								for (int y=0; y<rsImg.getHeight(); y++)
									comb.setRGB(x+xOffset, y, rsImg.getRGB(x, y));
							xOffset = rsImg.getWidth();
							for (int x=0; x<u3Img.getWidth(); x++)
								for (int y=0; y<u3Img.getHeight(); y++)
									comb.setRGB(x+xOffset, y, u3Img.getRGB(x, y));
							
							Graphics g = comb.getGraphics();
							g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 30));
							g.setColor(Color.BLACK);
							g.drawString(catalogName, 20, 40);
							g.drawString("UCERF3", rsImg.getWidth()+20, 40);
							
							ImageIO.write(comb, "png", new File(outputDir, prefix+".png"));
						}
					} catch (IOException e) {
						ExceptionUtils.throwAsRuntimeException(e);
					}
				}
			};
			if (disagg_parallel)
				submitRunnable(run, 100);
			else
				run.run();
		}
		
	}
	
	private static void plotCurveProfile(Map<Location, ? extends DiscretizedFunc> u3Curves,
			Map<Location, ? extends DiscretizedFunc> rsCurves, GriddedRegion gridReg, String imt,
			File outputDir, String prefix, Location startLoc, int num, double deltaLat, double deltaLon)
					throws IOException {
		CPT u3CPT = new CPT(0d, num, Color.BLUE, new Color(125, 125, 255));
		CPT rsCPT = new CPT(0d, num, Color.RED, new Color(255, 125, 125));
		
		// snap to grid
		startLoc = gridReg.getLocation(gridReg.indexForLocation(startLoc));
		double lat = startLoc.getLatitude();
		double lon = startLoc.getLongitude();
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<num; i++) {
			Location loc = new Location(lat+i*deltaLat, lon+i*deltaLon);
			// snap to grid
			loc = gridReg.getLocation(gridReg.indexForLocation(loc));
			
			Preconditions.checkNotNull(u3Curves.get(loc));
			
			funcs.add(u3Curves.get(loc));
			funcs.add(rsCurves.get(loc));
			if (i == 0) {
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, PlotSymbol.FILLED_CIRCLE, 2f, u3CPT.getColor((float)i)));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, PlotSymbol.FILLED_CIRCLE, 2f, rsCPT.getColor((float)i)));
			} else {
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, u3CPT.getColor((float)i)));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, rsCPT.getColor((float)i)));
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Hazard Curve Profile from "+(float)lat+", "+(float)lon, imt,
				"Annual Probabilitiy");
//		spec.setLegendVisible(true);
		
		PlotPreferences plotPrefs = buildPlotPrefs();
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		Range xRange = new Range(1e-2, 1e0);
		Range yRange = new Range(1e-6, 1e-2);
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(new File(outputDir, prefix).getAbsolutePath()+".png");
		gp.saveAsPDF(new File(outputDir, prefix).getAbsolutePath()+".pdf");
		gp.saveAsTXT(new File(outputDir, prefix).getAbsolutePath()+".txt");
	}
	
	private static PlotPreferences buildPlotPrefs() {
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(24);
		plotPrefs.setAxisLabelFontSize(26);
		plotPrefs.setPlotLabelFontSize(28);
		plotPrefs.setLegendFontSize(22);
		plotPrefs.setBackgroundColor(Color.WHITE);
		return plotPrefs;
	}

}
