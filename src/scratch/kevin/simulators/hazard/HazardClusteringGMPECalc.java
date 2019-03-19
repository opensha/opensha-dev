package scratch.kevin.simulators.hazard;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.utils.SimulatorUtils;
import org.opensha.sha.util.NEHRP_TestCity;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import oracle.net.aso.e;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.kevin.simCompare.GroundMotionScatterPlot;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import scratch.kevin.util.ReturnPeriodUtils;

public class HazardClusteringGMPECalc {

	public static void main(String[] args) throws IOException, DocumentException {
		File catalogsBaseDir = new File("/data/kevin/simulators/catalogs");
		File mainOutputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2381.instance(catalogsBaseDir);
		
		List<Double> rps = new ArrayList<>();
		
		double hazMinMag = 6.5;
		double skipYears = 5000;
		double minFractForInclusion = 0.2;
		double cutoffDist = 200d;
		
		double referenceDuration = 50;
		double[] referenceProbs = { 0.02, 0.10 };
		double[] tiOneYearComps = new double[referenceProbs.length];
		for (int i=0; i<referenceProbs.length; i++)
			tiOneYearComps[i] = ReturnPeriodUtils.calcExceedanceProb(referenceProbs[i], referenceDuration, 1d);
		
		rps.add(7d/365.25);
		rps.add(1d/12d);
		rps.add(1d);
		rps.add(50d);
		for (double referenceProb : referenceProbs)
			rps.add(ReturnPeriodUtils.calcDurationWithExceedanceProb(0.5, referenceProb, referenceDuration));
		rps.add(2500d);
		
		// for the plot as a function of RP
		double minFuncRP = 1d;
		double maxFuncRP = 2500;
		int numRPs = 50;
		
		double[] countMinMags = { 6d, 6.5, 7d };
		
		ScalarIMR gmpe = AttenRelRef.BSSA_2014.instance(null);
		gmpe.setParamDefaults();
		
		String imt = PGA_Param.NAME;
		double period = Double.NaN;
		
		String imtLabel = imt;
		String imtFileLabel = imt.toLowerCase();
		
		gmpe.setIntensityMeasure(imt);
		if (imt.equals(SA_Param.NAME)) {
			SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
			imtLabel = (float)period+"s "+imt;
			imtFileLabel += "_"+(float)period+"s";
		}
		
		List<String> lines = new ArrayList<>();
		
		File catOutDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catOutDir.exists() || catOutDir.mkdir());
		File catHazardOutDir = new File(catOutDir, "hazard_clustering_"+imtFileLabel);
		Preconditions.checkState(catHazardOutDir.exists() || catHazardOutDir.mkdir());
		File resourcesDir = new File(catHazardOutDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		// header
		lines.add("# Hazard Clustering Comparisons");
		lines.add("");
		lines.add("* GMPE: "+gmpe.getName());
		lines.add("* IMT: "+imtLabel);
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		Map<String, Location> curveLocs = new HashMap<>();
		curveLocs.put("Pasadena", new Location(34.148426, -118.17119));
		curveLocs.put("San Bernardino", new Location(34.064986, -117.29201));
		curveLocs.put("USC", new Location(34.0192, -118.286));
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA())
			curveLocs.put(city.toString(), city.location());
		
		List<String> siteNamesSorted = new ArrayList<>(curveLocs.keySet());
		Collections.sort(siteNamesSorted);
		
		File seedFile = new File(resourcesDir, "seed.txt");
		long seed;
		if (seedFile.exists()) {
			seed = Long.parseLong(Files.readLines(seedFile, Charset.defaultCharset()).get(0).trim());
			System.out.println("Loaded previous random seed: "+seed);
		} else {
			seed = System.currentTimeMillis();
			System.out.println("Writing random seed: "+seed);
			FileWriter fw = new FileWriter(seedFile);
			fw.write(seed+"\n");
			fw.close();
		}
		RandomGenerator rand = new Well19937c(seed);
		
		System.out.println("Loading catalog...");
		List<RSQSimEvent> events = catalog.loader().minMag(Math.min(hazMinMag, StatUtils.min(countMinMags)))
				.skipYears(skipYears).withinCutoffDist(cutoffDist, curveLocs.values()).load();
		double startYears = events.get(0).getTimeInYears();
		double durationYears = SimulatorUtils.getSimulationDurationYears(events);
		System.out.println("Loaded "+events.size()+" events ("+(float)durationYears+" years)");
		
		System.out.print("Building solution...");
		FaultSystemSolution sol = catalog.buildSolution(events, hazMinMag, minFractForInclusion);
		System.out.println("DONE");
		
		System.out.println("Calculating TI curves");
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(imt);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		Map<String, DiscretizedFunc> siteCurves = new HashMap<>();
		Map<String, Site> sitesMap = new HashMap<>();
		Map<String, double[]> siteEventIMLs = new HashMap<>();
		Map<String, RegionIden> siteIdens = new HashMap<>();
		for (String name : curveLocs.keySet()) {
			Location loc = curveLocs.get(name);
			Site site = new Site(loc);
			site.addParameterList(gmpe.getSiteParams());
			sitesMap.put(name, site);
			
			System.out.print("Calculating "+name+"...");
			calc.getHazardCurve(logXVals, site, gmpe, erf);
			System.out.println("DONE");
			DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<xVals.size(); i++)
				curve.set(xVals.getX(i), logXVals.getY(i));
			siteCurves.put(name, curve);
			double[] imls = new double[events.size()];
			for (int i=0; i<imls.length; i++)
				imls[i] = Double.NaN;
			siteEventIMLs.put(name, imls);
			siteIdens.put(name, new RegionIden(new Region(loc, cutoffDist)));
		}
		
		Color[] fractColors = { Color.BLUE.darker(), Color.GREEN.darker(), Color.RED.darker() };
		
		System.out.println("Calculating ground motions for all sites/events");
		int mod = (int)Math.round(events.size()/20d);
		for (int i=0; i<events.size(); i++) {
			if (i % mod == 0) {
				int percent = (int)Math.round(100d*i/events.size());
				System.out.println("Processing event "+i+"/"+events.size()+" ("+percent+" %)");
			}
			RSQSimEvent event = events.get(i);
			if (event.getMagnitude() < hazMinMag)
				continue;
			EqkRupture rup = catalog.getMappedSubSectRupture(event, minFractForInclusion);
			gmpe.setEqkRupture(rup);
			
			for (String name : siteNamesSorted) {
				if (!siteIdens.get(name).isMatch(event))
					continue;
				Site site = sitesMap.get(name);
				gmpe.setSite(site);
				
				double logMean = gmpe.getMean();
				double stdDev = gmpe.getStdDev();
				
				NormalDistribution normDist = new NormalDistribution(rand, logMean, stdDev);
				double logSample = normDist.sample();
				Preconditions.checkState(Double.isFinite(logSample), "bad log motion: %s", logSample);
				double motion = Math.exp(logSample);
				Preconditions.checkState(Double.isFinite(motion), "bad motion. e^%s = %s", logSample, motion);
				siteEventIMLs.get(name)[i] = motion;
			}
		}
		
		for (double rp : rps) {
			double[] exceedProbs = new double[referenceProbs.length];
			for (int i=0; i<referenceProbs.length; i++)
				exceedProbs[i] = ReturnPeriodUtils.calcExceedanceProb(referenceProbs[i], referenceDuration, rp);
			
			int numWindows = (int)(durationYears/rp);
			
			int[][] counts = new int[countMinMags.length][];
			int[][] poissonCounts = new int[countMinMags.length][];
			
			String lenStr;
			if (rp < 1d) {
				if (rp == 1d/12d)
					lenStr = "1mo";
				else if (rp == 7d/365.25)
					lenStr = "1wk";
				else
					lenStr = rpDF.format(rp)+"d";
			} else {
				lenStr = rpDF.format(rp)+"yr";
			}
			
			lines.add("## "+lenStr+" Sub-Catalogs");
			lines.add(topLink); lines.add("");
			lines.add("**"+numWindows+" sub-catalogs**");
			lines.add("");
			TableBuilder table = MarkdownUtils.tableBuilder();
			table.addLine("Fractile", "Time-Independent Probability");
			for (int i=0; i<exceedProbs.length; i++)
				table.addLine("**"+formatProb(exceedProbs[i])+"**", "**"+formatProb(tiOneYearComps[i])+"**");
			lines.addAll(table.build());
			
			System.out.println("Calculating for rp="+rp+" with "+numWindows+" windows");
			if (numWindows < 10) {
				System.out.println("Skipping, not enough windows");
				continue;
			}

			Map<String, double[]> sitePeakMotions = calcSitePeakMotions(events, siteNamesSorted,
					siteEventIMLs, rp, false, rand, hazMinMag, countMinMags, counts);
			Map<String, double[]> sitePeakPoissonMotions = calcSitePeakMotions(events, siteNamesSorted,
					siteEventIMLs, rp, true, rand, hazMinMag, countMinMags, poissonCounts);
			
			DefaultXY_DataSet[] scatters = new DefaultXY_DataSet[exceedProbs.length];
			for (int j=0; j<scatters.length; j++)
				scatters[j] = new DefaultXY_DataSet();
			
			Map<String, String> sitePlotNames = new HashMap<>();
			for (String name : curveLocs.keySet()) {
				ArbitrarilyDiscretizedFunc fractFunc = calcSiteFractileFunc(sitePeakMotions.get(name));
				ArbitrarilyDiscretizedFunc fractFuncPoisson = calcSiteFractileFunc(sitePeakPoissonMotions.get(name));
				
				double maxIML = Math.max(fractFunc.getMaxX(), fractFuncPoisson.getMaxX());
				double minIML = Double.POSITIVE_INFINITY;
				for (Point2D pt : fractFunc)
					if (pt.getX() > 0)
						minIML = Math.min(minIML, pt.getX());
				for (Point2D pt : fractFuncPoisson)
					if (pt.getX() > 0)
						minIML = Math.min(minIML, pt.getX());
				
				boolean logX = true;
				double minX = minIML;
				double maxX = maxIML;
				
				boolean logY = true;
				double minY = 1d/(double)numWindows;
				double maxY = 1d;
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				funcs.add(fractFuncPoisson);
				fractFuncPoisson.setName("Poisson");
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY));
				
				funcs.add(fractFunc);
				fractFunc.setName("Simulation");
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
				
				DiscretizedFunc tiCurve = siteCurves.get(name);
				
				for (int i=0; i<exceedProbs.length; i++) {
					double exceedProb = exceedProbs[i];
					double tiComp = tiOneYearComps[i];
					
					double catVal = fractFunc.getFirstInterpolatedX(exceedProb);
					double tiVal = tiCurve.getFirstInterpolatedX_inLogXLogYDomain(tiComp);
					
					scatters[i].set(tiVal, catVal);
					
					PlotCurveCharacterstics fractChar = new PlotCurveCharacterstics(
							PlotLineType.SOLID, 2f, fractColors[i % fractColors.length]);
					PlotCurveCharacterstics tiChar = new PlotCurveCharacterstics(
							PlotLineType.DASHED, 2f, fractColors[i % fractColors.length].brighter());
					
					funcs.add(line("p"+formatProb(exceedProb)+": "+imlDF.format(catVal), minX, exceedProb, catVal, exceedProb));
					chars.add(fractChar);
					funcs.add(line(null, catVal, minY, catVal, exceedProb));
					chars.add(fractChar);
					
					funcs.add(line("TI Comp = "+imlDF.format(tiVal), tiVal, minY, tiVal, maxY));
					chars.add(tiChar);
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, name+", "+lenStr+" Comparisons", imtLabel, lenStr+" Exceed. Prob");
				spec.setLegendVisible(true);
				
				String prefix = name.replaceAll(" ", "_")+"_"+lenStr;
				
				PlotPreferences plotPrefs = PlotPreferences.getDefault();
				plotPrefs.setTickLabelFontSize(18);
				plotPrefs.setAxisLabelFontSize(20);
				plotPrefs.setPlotLabelFontSize(21);
				plotPrefs.setBackgroundColor(Color.WHITE);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
				
				gp.drawGraphPanel(spec, logX, logY, new Range(minX, maxX), new Range(minY, maxY));
				gp.getChartPanel().setSize(800, 600);
				File pngFile = new File(resourcesDir, prefix+".png");
				gp.saveAsPNG(pngFile.getAbsolutePath());
				sitePlotNames.put(name, prefix+".png");
			}
			
			lines.add("### "+lenStr+" Scatters");
			lines.add(topLink); lines.add("");
			for (int i=0; i<exceedProbs.length; i++) {
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				funcs.add(scatters[i]);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, Color.BLACK));
				
				GroundMotionScatterPlot.PLOT_WIDTH = 600;
				GroundMotionScatterPlot.WRITE_PDF = false;
				GroundMotionScatterPlot.YELLOW_REGION = false;
				GroundMotionScatterPlot.SCATTER_QUANTITY_NAME = "Sites";
				List<String> binDescriptions = new ArrayList<>();
				binDescriptions.add("Cat Exceed. Prob = "+formatProb(exceedProbs[i]));
				binDescriptions.add("TI Prob = "+formatProb(tiOneYearComps[i]));
				String prefix = lenStr+"_scatter_p"+formatProb(exceedProbs[i]);
				GroundMotionScatterPlot.plot(scatters[i], "Time-Independent "+imtLabel, "Sub-Catalogs "+imtLabel,
						binDescriptions, lenStr+" Scatter, p="+formatProb(exceedProbs[i]), resourcesDir, prefix);
				
				lines.add("!["+formatProb(exceedProbs[i])+" Fract Scatter]("+resourcesDir.getName()+"/"+prefix+".png)");
				lines.add("");
			}
			
			lines.add("### "+lenStr+" Event Count Hists");
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("Min Mag", "Simulation", "Poisson");
			
			boolean[] poissons =  { false, true };
			
			for (int m=0; m<countMinMags.length; m++) {
				table.initNewLine();
				table.addColumn("M≥"+(float)countMinMags[m]);
				HistogramFunction[] hists = new HistogramFunction[poissons.length];
				double[][] countsDoubles = new double[poissons.length][];
				
				int min = Integer.MAX_VALUE;
				int max = 0;
				for (int[] myCounts : new int[][] { counts[m], poissonCounts[m] }) {
					for (int count : myCounts) {
						if (count < min)
							min = count;
						if (count > max)
							max = count;
					}
				}
				int span = max - min;
				double delta;
				if (span > 1000)
					delta = 100;
				else if (span > 500)
					delta = 50;
				else if (span > 100)
					delta = 10;
				else if (span > 50)
					delta = 5;
				else if (span > 20)
					delta = 2;
				else
					delta = 1;
				
				double minX = Double.POSITIVE_INFINITY;
				double maxX = 0;
				double maxY = 0;
				for (int p=0; p<poissons.length; p++) {
					boolean poisson = poissons[p];
					int[] myCounts;
					if (poisson)
						myCounts = poissonCounts[m];
					else
						myCounts = counts[m];
					
					hists[p] = HistogramFunction.getEncompassingHistogram((double)min, (double)max, delta);
					hists[p] = new HistogramFunction(hists[p].getMinX()-0.5, hists[p].getMaxX()+0.5, hists[p].size()+1);
					countsDoubles[p] = new double[myCounts.length];
					for (int i=0; i<myCounts.length; i++) {
						int count = myCounts[i];
						hists[p].add((double)count, 1d);
						countsDoubles[p][i] = (double)count;
					}
					
					minX = Math.min(minX, hists[p].getMinX()-0.75*delta);
					maxX = Math.max(maxX, hists[p].getMaxX()+0.75*delta);
					maxY = Math.max(maxY, hists[p].getMaxY()*1.1);
				}
				for (int p=0; p<poissons.length; p++) {
					
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					funcs.add(hists[p]);
					hists[p].setName("Histogram");
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
					
					double mean = StatUtils.mean(countsDoubles[p]);
					double median = DataUtils.median(countsDoubles[p]);
					
					funcs.add(line("Mean = "+imlDF.format(mean), mean, 0, mean, maxY));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
					
					funcs.add(line("Median = "+imlDF.format(median), median, 0, median, maxY));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
					
					String title, prefix;
					String minMagStr;
					if (countMinMags[m] == Math.round(countMinMags[m]))
						minMagStr = (int)countMinMags[m]+"";
					else
						minMagStr = (float)countMinMags[m]+"";
					if (poissons[p]) {
						title = "Poisson M≥"+minMagStr+" Count";
						prefix = lenStr+"_m"+minMagStr+"_hist_poisson";
					} else {
						title = "Simulation M≥"+minMagStr+" Count";
						prefix = lenStr+"_m"+minMagStr+"_hist";
					}
					PlotSpec spec = new PlotSpec(funcs, chars, title, "Count", "Num Sub-Catalogs");
					spec.setLegendVisible(true);
					
					PlotPreferences plotPrefs = PlotPreferences.getDefault();
					plotPrefs.setTickLabelFontSize(18);
					plotPrefs.setAxisLabelFontSize(20);
					plotPrefs.setPlotLabelFontSize(21);
					plotPrefs.setBackgroundColor(Color.WHITE);
					
					HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
					
					gp.drawGraphPanel(spec, false, false, new Range(minX, maxX), new Range(0, maxY));
					gp.getChartPanel().setSize(800, 600);
					File pngFile = new File(resourcesDir, prefix+".png");
					gp.saveAsPNG(pngFile.getAbsolutePath());
					table.addColumn("!["+title+"]("+resourcesDir.getName()+"/"+prefix+".png)");
				}
			}
			
			table.finalizeLine();
			lines.addAll(table.build());
			lines.add("");
			
			lines.add("### "+lenStr+" Site Fractile Curves");
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("Site", "Fractile Curve Plot");
			
			for (String name : siteNamesSorted)
				table.addLine("**"+name+"**", "!["+name+"]("+resourcesDir.getName()+"/"+sitePlotNames.get(name)+")");
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// now as a function of RP
		EvenlyDiscretizedFunc logRPfunc = new EvenlyDiscretizedFunc(Math.log10(minFuncRP), Math.log10(maxFuncRP), numRPs);
		double[] funcRPs = new double[logRPfunc.size()];
		for (int i=0; i<logRPfunc.size(); i++)
			funcRPs[i] = Math.pow(10, logRPfunc.getX(i));
		
		Map<String, ArbitrarilyDiscretizedFunc> siteFuncRPs = new HashMap<>();
		Map<String, ArbitrarilyDiscretizedFunc> siteFuncPoissonRPs = new HashMap<>();
		for (String name : siteNamesSorted) {
			siteFuncRPs.put(name, new ArbitrarilyDiscretizedFunc());
			siteFuncPoissonRPs.put(name, new ArbitrarilyDiscretizedFunc());
		}
		
		double r1 = 0.02;
		double t1 = 50;
		String tiRPfuncLabel = "2% in 50 years";
		
		for (double rp : funcRPs) {
			Map<String, double[]> sitePeakMotions =
					calcSitePeakMotions(events, siteNamesSorted, siteEventIMLs, rp, false, rand, hazMinMag, null, null);
			Map<String, double[]> sitePoissonPeakMotions =
					calcSitePeakMotions(events, siteNamesSorted, siteEventIMLs, rp, true, rand, hazMinMag, null, null);
			
//			double t2 = rp;
//			double r2Star = r1Star*t2/t1;
//			
//			double r2 = Math.sqrt(1 + 2*r2Star) - 1;
			double r2 = ReturnPeriodUtils.calcExceedanceProb(r1, t1, rp);
			
			for (String name : siteNamesSorted) {
				ArbitrarilyDiscretizedFunc fractFunc = calcSiteFractileFunc(sitePeakMotions.get(name));
				ArbitrarilyDiscretizedFunc fractPoissonFunc = calcSiteFractileFunc(sitePoissonPeakMotions.get(name));
				
				siteFuncRPs.get(name).set(rp, fractFunc.getFirstInterpolatedX(r2));
				siteFuncPoissonRPs.get(name).set(rp, fractPoissonFunc.getFirstInterpolatedX(r2));
			}
		}
		
		lines.add("## Site "+tiRPfuncLabel+" RP Functions");
		lines.add(topLink); lines.add("");
//		lines.add("**Time Independent Probability Level: "+(float)tiRPfuncComp+"**");
//		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Site", "Fractile Curve Plot");
		
		double tiRPfuncComp = ReturnPeriodUtils.calcExceedanceProb(r1, t1, 1d);
		
		for (String name : siteNamesSorted) {
			double tiVal = siteCurves.get(name).getFirstInterpolatedX_inLogXLogYDomain(tiRPfuncComp);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			ArbitrarilyDiscretizedFunc siteRPfunc = siteFuncRPs.get(name);
			ArbitrarilyDiscretizedFunc siteRPfuncPoisson = siteFuncPoissonRPs.get(name);
			
			funcs.add(line("TI "+tiRPfuncLabel, minFuncRP, tiVal, maxFuncRP, tiVal));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, fractColors[0]));
			
			funcs.add(siteRPfuncPoisson);
			siteRPfuncPoisson.setName("Poisson");
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY));
			
			funcs.add(siteRPfunc);
			siteRPfunc.setName("Simulation");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, name+" RP Comparisons", "Sub-Catalog Length", tiRPfuncLabel+", "+imtLabel);
			spec.setLegendVisible(true);
			
			String prefix = name.replaceAll(" ", "_")+"_rp_func";
			
			PlotPreferences plotPrefs = PlotPreferences.getDefault();
			plotPrefs.setTickLabelFontSize(18);
			plotPrefs.setAxisLabelFontSize(20);
			plotPrefs.setPlotLabelFontSize(21);
			plotPrefs.setBackgroundColor(Color.WHITE);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			gp.drawGraphPanel(spec, true, false, new Range(minFuncRP, maxFuncRP), null);
			gp.getChartPanel().setSize(800, 600);
			File pngFile = new File(resourcesDir, prefix+".png");
			gp.saveAsPNG(pngFile.getAbsolutePath());
			table.addLine("**"+name+"**", "!["+name+"]("+resourcesDir.getName()+"/"+prefix+".png)");
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, catHazardOutDir);
		
		catalog.writeMarkdownSummary(catOutDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
	}
	
	private static DecimalFormat rpDF = new DecimalFormat("0.#");
	
	private static DecimalFormat imlDF = new DecimalFormat("0.###");
	
	private static DecimalFormat probDecimalDF = new DecimalFormat("0.##");
	private static DecimalFormat probScientificDF = new DecimalFormat("0.#E0");
	private static String formatProb(double prob) {
		if (prob < 0.005)
			return probScientificDF.format(prob);
		return probDecimalDF.format(prob);
	}
	
	private static XY_DataSet line(String name, double x0, double y0, double x1, double y1) {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		xy.set(x0, y0);
		xy.set(x1, y1);
		xy.setName(name);
		
		return xy;
	}
	
	private static Map<String, double[]> calcSitePeakMotions(List<RSQSimEvent> events, List<String> siteNamesSorted,
			Map<String, double[]> siteEventIMLs, double rp, boolean poisson, RandomGenerator rand, double hazardMin,
			double[] countMags, int[][] counts) {
		double durationYears = SimulatorUtils.getSimulationDurationYears(events);
		double startYears = events.get(0).getTimeInYears();
		
		Preconditions.checkState(countMags == null || countMags.length == counts.length);
		
		int numWindows = (int)(durationYears/rp);
		if (countMags != null) {
			Preconditions.checkState(countMags.length == counts.length);
			for (int i=0; i<counts.length; i++)
				counts[i] = new int[numWindows];
		}
		
		Map<String, double[]> sitePeakMotions = new HashMap<>();
		
		for (String name : siteNamesSorted)
			sitePeakMotions.put(name, new double[numWindows]);
		
		for (int e=0; e<events.size(); e++) {
			RSQSimEvent event = events.get(e);
			double time; // relative to the start of the catalog
			if (poisson)
				time = rand.nextDouble()*durationYears;
			else
				time = event.getTimeInYears() - startYears;
			int w = (int)(time/(double)rp);
			Preconditions.checkState(w >= 0 && w <= numWindows, "Bad w=%s with t=%s");
			if (w == numWindows)
				// partial window at the end
				continue;
			
			if (countMags != null) {
				for (int i=0; i<countMags.length; i++)
					if (event.getMagnitude() >= countMags[i])
						counts[i][w]++;
			}
			
			if (event.getMagnitude() >= hazardMin) {
				for (String name : siteNamesSorted) {
					double motion = siteEventIMLs.get(name)[e];
					if (Double.isNaN(motion))
						// more than 200km
						continue;
					
					double[] sitePeaks = sitePeakMotions.get(name);
					sitePeaks[w] = Math.max(sitePeaks[w], motion);
				}
			}
		}
		
		return sitePeakMotions;
	}
	
	private static ArbitrarilyDiscretizedFunc calcSiteFractileFunc(double[] sitePeaks) {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		
		Arrays.sort(sitePeaks);
		
		for (int p=0; p<sitePeaks.length; p++) {
			int numExceeds = sitePeaks.length - p;
			double fractExceeds = (double)numExceeds/(double)sitePeaks.length;
			if (!func.hasX(sitePeaks[p]))
				func.set(sitePeaks[p], fractExceeds);
		}
		
		return func;
	}

}
