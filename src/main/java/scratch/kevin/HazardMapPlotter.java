package scratch.kevin;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.ReturnPeriodUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculator;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.gui.infoTools.DisaggregationPlotViewerWindow;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class HazardMapPlotter {

	public static void main(String[] args) throws Exception {
		File jobDir = new File("/home/kevin/OpenSHA/maps/2018_01_23-ucerf3-ba-la-zoomed-site-effects-pga");
//		File jobDir = new File("/home/kevin/OpenSHA/maps/2018_01_23-ucerf3-ba-la-zoomed-site-effects-1s");
//		File jobDir = new File("/home/kevin/OpenSHA/maps/2018_01_23-ucerf3-ba-la-zoo	med-pga-rock");
//		Region region = new CaliforniaRegions.LA_BOX();
		Region region = new Region(new Location(34.5, -117.2), new Location(33.6, -119));
		double spacing = 0.005;
		File inputFile = new File(jobDir, "inputs.xml");
		File curvesFile = new File(new File(jobDir, "curves"), "imrs1.bin");
		File localFSSDir = new File("/home/kevin/.opensha/ucerf3_erf");
		
		List<SiteData<Double>> dataProvs = new ArrayList<>();
		if (!jobDir.getName().contains("rock")) {
			dataProvs.add(new WillsMap2015());
			dataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			dataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		}
		
		boolean plotMaps = true;
		boolean plotCurves = false;
		boolean plotDisagg = false;
		
		File mapDir = new File(jobDir, "maps");
		Preconditions.checkState(mapDir.exists() || mapDir.mkdir());
		
		File curveDir = new File(jobDir, "curves");
		Preconditions.checkState(curveDir.exists() || curveDir.mkdir());
		
		File disaggDir = new File(jobDir, "disagg");
		Preconditions.checkState(disaggDir.exists() || disaggDir.mkdir());
		
		double referenceDuration = 50d;
		double[] refProbs = { 0.02, 0.1 };
//		double[] mapMaxes = { 2d, 1.5d }
		double[] mapMaxes = null;
		
		Map<String, Location> curveSites = new HashMap<>();
		curveSites.put("USC", new Location(34.0192, -118.286));
		curveSites.put("Caltech", new Location(34.137652, -118.125129));
		
		Document doc = XMLUtils.loadDocument(inputFile);
		CalculationInputsXMLFile.UDPATE_FORECAST = false;
		CalculationInputsXMLFile inputs = CalculationInputsXMLFile.loadXML(doc);
		
		ScalarIMR gmpe = inputs.getIMRMaps().iterator().next().values().iterator().next();
		String imtLabel = gmpe.getIntensityMeasure().getName()+" ("+gmpe.getIntensityMeasure().getUnits()+")";
		if (imtLabel.equals(SA_Param.NAME))
			imtLabel = (float)SA_Param.getPeriodInSA_Param(gmpe.getIntensityMeasure())+"s "+imtLabel;
		
		ERF erf = inputs.getERF();
		double calcDuration = erf.getTimeSpan().getDuration();
		
		// load curves
		System.out.println("Loading curves");
		BinaryHazardCurveReader read = new BinaryHazardCurveReader(curvesFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> curvesMap = read.getCurveMap();
		
		CPT hazardCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		hazardCPT.setNanColor(Color.GRAY);
		
		String[] probLabels = new String[refProbs.length];
		String[] probFileLabels = new String[refProbs.length];
		double[] calcProbs = new double[refProbs.length];
		for (int i=0; i<refProbs.length; i++) {
			probLabels[i] = optionalDigitDF.format(refProbs[i]*100d)+"% in "+optionalDigitDF.format(referenceDuration)+"yr";
			calcProbs[i] = ReturnPeriodUtils.calcExceedanceProb(refProbs[i], referenceDuration, calcDuration);
			probFileLabels[i] = probLabels[i].replace("%", "p").replaceAll(" ", "_");
		}
		
		if (plotMaps) {
			for (int i=0; i<refProbs.length; i++) {
				String label = probLabels[i]+", "+imtLabel;
				double calcProb = calcProbs[i];
				String fileLabel = probFileLabels[i];
				System.out.println("Generating map for "+label);
				System.out.println(label+" = "+(float)calcProb+" in "+optionalDigitDF.format(calcDuration));
				GeoDataSet baseMap = HazardDataSetLoader.extractPointFromCurves(curvesMap, false, calcProb);
				
				double mapMax = baseMap.getMaxZ();
				Double customMin = 0d;
				Double customMax;
				if (mapMaxes != null) {
					customMax = mapMaxes[i];
				} else {
					if (mapMax < 0.75)
						customMax = 0.6;
					else if (mapMax < 1.15)
						customMax = 1d;
					else if (mapMax < 1.65)
						customMax = 1.5;
					else
						customMax = 2d;
				}
				
				String prefix = fileLabel;
				
				plotMap(mapDir, prefix, baseMap, region, spacing, customMin, customMax, label, hazardCPT, true);
			}
			
			for (SiteData<Double> prov : dataProvs) {
				GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
				ArrayList<Double> datas = prov.getValues(gridReg.getNodeList());
				
				GeoDataSet map = new GriddedGeoDataSet(gridReg, false);
				for (int i=0; i<datas.size(); i++)
					map.set(i, datas.get(i));
				
				String prefix = "site_data_"+prov.getShortName();
				String label;
				Double customMin, customMax;
				if (prov.getDataType().equals(SiteData.TYPE_VS30)) {
					customMin = map.getMinZ();
					customMax = map.getMaxZ();
					prefix += "_vs30";
					label = prov.getName()+", Vs30 (m/s)";
				} else if (prov.getDataType().equals(SiteData.TYPE_DEPTH_TO_1_0)) {
					customMin = 0d;
					customMax = map.getMaxZ();
					prefix += "_z10";
					label = prov.getShortName()+", Z1.0 (km)";
				} else if (prov.getDataType().equals(SiteData.TYPE_DEPTH_TO_2_5)) {
					customMin = 0d;
					customMax = map.getMaxZ();
					prefix += "_z25";
					label = prov.getShortName()+", Z2.5 (km)";
				} else {
					customMin = null;
					customMax = null;
					label = prov.getName()+", "+prov.getDataType();
				}
				
				plotMap(mapDir, prefix, map, region, spacing, customMin, customMax, label, hazardCPT, true);
			}
		}
		
		if (plotCurves || plotDisagg) {
			if (localFSSDir != null && erf instanceof FaultSystemSolutionERF) {
				File curFSS = (File) erf.getAdjustableParameterList().getParameter(FaultSystemSolutionERF.FILE_PARAM_NAME).getValue();
				File fssFile = new File(localFSSDir, curFSS.getName());
				Preconditions.checkState(fssFile.exists(), "FSS file doesn't exist: %s", fssFile.getAbsolutePath());
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, fssFile);
			}
			erf.updateForecast();
			
			SiteTranslator trans = new SiteTranslator();
			HazardCurveCalculator calc = new HazardCurveCalculator();
			DisaggregationCalculator disaggCalc = new DisaggregationCalculator();
			
			DiscretizedFunc xVals = inputs.getCalcSettings().getXValues(gmpe.getIntensityMeasure().getName());
			DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : xVals)
				logXVals.set(Math.log(pt.getX()), 1d);
			
			for (String siteName : curveSites.keySet()) {
				System.out.println("Calculating curve for "+siteName);
				Location loc = curveSites.get(siteName);
				Site site = new Site(loc);
				
				List<SiteDataValue<?>> siteDatas = new ArrayList<>();
				for (SiteData<?> prov : dataProvs)
					siteDatas.add(prov.getAnnotatedValue(loc));
				for (Parameter<?> param : gmpe.getSiteParams()) {
					param = (Parameter<?>) param.clone();
					trans.setParameterValue(param, siteDatas);
					site.addParameter(param);
				}
				
				DiscretizedFunc logCurve = logXVals.deepClone();
				calc.getHazardCurve(logCurve, site, gmpe, erf);
				DiscretizedFunc linearCurve = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<logCurve.size(); i++)
					linearCurve.set(xVals.getX(i), logCurve.getY(i));
				
				if (plotCurves) {
					List<DiscretizedFunc> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					Range xRange = new Range(1e-3, 1e1);
					Range yRange = new Range(1e-6, 1e0);
					
					linearCurve.setName("Hazard Curve");
					funcs.add(linearCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					
					if (calcProbs != null && calcProbs.length > 0) {
						CPT probCPT = getProbLogCPT(calcProbs);
						for (int i=0; i<calcProbs.length; i++) {
							double prob = calcProbs[i];
							Color color = probCPT.getColor((float)Math.log10(prob));
							DiscretizedFunc probLine = new ArbitrarilyDiscretizedFunc();
							probLine.set(xRange.getLowerBound(), prob);
							probLine.set(xRange.getUpperBound(), prob);
							probLine.setName(probLabels[i]);
							funcs.add(probLine);
							chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, color));
						}
					}

					PlotSpec spec = new PlotSpec(funcs, chars, siteName+" Hazard Curve",
							imtLabel, optionalDigitDF.format(calcDuration)+"yr Probability");
					spec.setLegendVisible(true);
					
					PlotPreferences plotPrefs = PlotPreferences.getDefault();
					plotPrefs.setTickLabelFontSize(18);
					plotPrefs.setAxisLabelFontSize(20);
					plotPrefs.setPlotLabelFontSize(21);
					plotPrefs.setBackgroundColor(Color.WHITE);
					
					HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
					
					String prefix = "curve_"+siteName.replaceAll(" ", "_");
					
					gp.drawGraphPanel(spec, true, true, xRange, yRange);
					gp.getChartPanel().setSize(800, 600);
					gp.saveAsPNG(new File(curveDir, prefix).getAbsolutePath()+".png");
					gp.saveAsPDF(new File(curveDir, prefix).getAbsolutePath()+".pdf");
					gp.saveAsTXT(new File(curveDir, prefix).getAbsolutePath()+".txt");
				}
				
				if (plotDisagg) {
					double minMag = 6d;
					double deltaMag = 0.2;
					int numMags = (int)((8.6d - minMag)/deltaMag + 0.5);
					
					int numSourcesForDisag = 100;
					
					boolean showSourceDistances = true;
					
					double maxZAxis = Double.NaN;
					
					ParameterList disaggParams = DisaggregationCalculator.getDefaultParams();
					
					for (int i=0; i<calcProbs.length; i++) {
						double prob = calcProbs[i];
						double iml = HazardDataSetLoader.getCurveVal(linearCurve, false, prob); // iml at prob
						if (!Double.isFinite(iml)) {
							System.out.println("Couldn't get IML for "+siteName+", prob="+prob+"yr. Skipping disagg!");
							return;
						}

						System.out.println("Disaggregating for prob="+prob+", iml="+iml);
						disaggCalc.setMagRange(minMag, numMags, deltaMag);
						disaggCalc.setNumSourcestoShow(numSourcesForDisag);
						disaggCalc.setShowDistances(showSourceDistances);
						boolean success = disaggCalc.disaggregate(Math.log(iml), site, gmpe, (AbstractERF)erf, disaggParams);
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

						String prefix = "disagg_"+siteName.replaceAll(" ", "_")+"_"+probFileLabels[i]
								+"_"+threeDigitDF.format(iml)+gmpe.getIntensityMeasure().getUnits();
						File outputFile = new File(disaggDir, prefix);

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
					}
				}
			}
		}
	}
	
	private static void plotMap(File outputDir, String prefix, GeoDataSet baseMap, Region region, double spacing,
			Double customMin, Double customMax, String label, CPT cpt, boolean rescaleCPT)
					throws IOException, ClassNotFoundException, GMT_MapException {
		GMT_Map map = new GMT_Map(region, baseMap, spacing, cpt);
		map.setCustomLabel(label);
		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
		map.setUseGMTSmoothing(true);
		map.setLogPlot(false);
//		map.setDpi(300);
//		map.setImageWidth(10d);
		map.setXyzFileName("base_map.xyz");
		map.setCustomScaleMin(customMin);
		map.setCustomScaleMax(customMax);
		map.setBlackBackground(false);
		map.setRescaleCPT(rescaleCPT);
		map.setJPGFileName(null);
		
		if (!rescaleCPT) {
			map.setCPTEqualSpacing(false);
			map.setCPTCustomInterval(0.1);
		}
		
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
	}
	
	private static final DecimalFormat optionalDigitDF = new DecimalFormat("0.#");
	
	private static final DecimalFormat threeDigitDF = new DecimalFormat("0.###");
	
	public static CPT getProbLogCPT(double[] probs) {
		double minProb = Double.POSITIVE_INFINITY;
		double maxProb = 0;
		
		for (double prob : probs) {
			if (prob < minProb)
				minProb = prob;
			if (prob > maxProb)
				maxProb = prob;
		}
		return new CPT(Math.log10(minProb), Math.log10(maxProb), Color.LIGHT_GRAY, Color.DARK_GRAY);
	}

}
