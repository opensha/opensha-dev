package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.COMBINED_DIR;
import static scratch.kevin.prvi25.figures.PRVI_Paths.COMBINED_SOL;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.MarkdownUtils.TableTextAlignment;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculator.EpsilonCategories;
import org.opensha.sha.calc.disaggregation.DisaggregationPlotData;
import org.opensha.sha.calc.disaggregation.chart3d.PureJavaDisaggPlotter;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.faultSysSolution.util.SolSiteHazardCalc;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import net.mahdilamb.colormap.Colors;

public class SiteHazardInvestigationsPageGen {

	public static void main(String[] args) throws IOException {
		List<Site> sites = PRVI25_RegionLoader.loadHazardSites();
		File mainOoutputDir = new File(COMBINED_DIR, "site_hazard_investigations");
		Preconditions.checkState(mainOoutputDir.exists() || mainOoutputDir.mkdir());
		
		FaultSystemSolution combSol = FaultSystemSolution.load(COMBINED_SOL);
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(combSol);
		erf.getTimeSpan().setDuration(1d);
		erf.setCacheGridSources(true);
		erf.setGriddedSeismicitySettings(erf.getGriddedSeismicitySettings().forSupersamplingSettings(GridCellSupersamplingSettings.DEFAULT));
		erf.updateForecast();
		
		BaseFaultSystemSolutionERF disaggIndvERF = new BaseFaultSystemSolutionERF();
		disaggIndvERF.setSolution(combSol);
		disaggIndvERF.getTimeSpan().setDuration(1d);
		disaggIndvERF.setCacheGridSources(true);
		disaggIndvERF.setGriddedSeismicitySettings(erf.getGriddedSeismicitySettings().forSupersamplingSettings(null));
		disaggIndvERF.updateForecast();
		
		HazardModel prevModel = HazardModel.load(Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-prvi-2003-main"));
		// needs to be include for slab and interface
		NshmErf prevERF = new NshmErf(prevModel, EnumSet.allOf(TectonicRegionType.class),
				IncludeBackgroundOption.INCLUDE);
		prevERF.getTimeSpan().setDuration(1d);
		prevERF.updateForecast();
		
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		List<Future<?>> futures = new ArrayList<>();
		
		for (Site site : sites) {
			futures.add(exec.submit(new Runnable() {
				
				@Override
				public void run() {
					try {
						writeSitePage(mainOoutputDir, site, erf, disaggIndvERF, prevERF);
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
				}
			}));
		}
		
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		exec.shutdown();
	}
	
	private static void writeSitePage(File mainOutputDir, Site site, BaseFaultSystemSolutionERF erf,
			BaseFaultSystemSolutionERF disaggIndvERF, NshmErf prevERF) throws IOException {
		String siteName = site.getName();
		System.out.println("Starting "+siteName);
		String sitePrefix = FileNameUtils.simplify(siteName);
		File outputDir = new File(mainOutputDir, sitePrefix);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());

		DecimalFormat oDF = new DecimalFormat("0.##");
		DecimalFormat pDF = new DecimalFormat("0.##%");
		
		Range mfdXRange = new Range(5d, 9.2d);
		Range mfdYRange = new Range(1e-6, 1e-1);
		
		double[] radii = { 20d, 50d, 100d, 200d };
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 9.2);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# "+siteName+" Hazard Investigations");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Nearby Source MFDs");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder().textAlign(TableTextAlignment.CENTER);
		
		Map<TectonicRegionType, Color> trtColors = Map.of(TectonicRegionType.ACTIVE_SHALLOW, Colors.tab_blue,
				TectonicRegionType.SUBDUCTION_INTERFACE, Colors.tab_orange,
				TectonicRegionType.SUBDUCTION_SLAB, Colors.tab_green);
		
		table.addLine("__Incremental__", "__Cumulative__");
		for (double radius : radii) {
			table.addLine("__Within "+oDF.format(radius)+" km (in 3D)__", "");
			
			Map<TectonicRegionType, IncrementalMagFreqDist> mfds = calcTRT_MFDs(erf, site, radius, refMFD);
			Map<TectonicRegionType, IncrementalMagFreqDist> prevMFDs = calcTRT_MFDs(prevERF, site, radius, refMFD);
			
			List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			IncrementalMagFreqDist activeMFD = mfds.get(TectonicRegionType.ACTIVE_SHALLOW);
			activeMFD.setName("Crustal (2025)");
			incrFuncs.add(activeMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
			
			IncrementalMagFreqDist interfaceMFD = mfds.get(TectonicRegionType.SUBDUCTION_INTERFACE);
			interfaceMFD.setName("Interface (2025)");
			incrFuncs.add(interfaceMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
			
			IncrementalMagFreqDist slabMFD = mfds.get(TectonicRegionType.SUBDUCTION_SLAB);
			slabMFD.setName("Intraslab (2025)");
			incrFuncs.add(slabMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
			
			IncrementalMagFreqDist prevActiveMFD = prevMFDs.get(TectonicRegionType.ACTIVE_SHALLOW);
			prevActiveMFD.setName("Crustal (2003)");
			incrFuncs.add(prevActiveMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_blue));
			
			IncrementalMagFreqDist prevInterfaceMFD = prevMFDs.get(TectonicRegionType.SUBDUCTION_INTERFACE);
			prevInterfaceMFD.setName("Interface (2003)");
			incrFuncs.add(prevInterfaceMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_orange));
			
			IncrementalMagFreqDist prevSlabMFD = prevMFDs.get(TectonicRegionType.SUBDUCTION_SLAB);
			prevSlabMFD.setName("Intraslab (2003)");
			incrFuncs.add(prevSlabMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_green));
			
			List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
			for (IncrementalMagFreqDist mfd : incrFuncs)
				cmlFuncs.add(mfd.getCumRateDistWithOffset());
			
			PlotSpec incrPlot = new PlotSpec(incrFuncs, chars, "Within "+oDF.format(radius)+" km of "+siteName,
					"Magnitude", "Incremental Rate (/year)");
			incrPlot.setLegendInset(true);
			
			PlotSpec cmlPlot = new PlotSpec(cmlFuncs, chars, "Within "+oDF.format(radius)+" km of "+siteName,
					"Magnitude", "Cumulative Rate (/year)");
			cmlPlot.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			
			gp.drawGraphPanel(incrPlot, false, true, mfdXRange, mfdYRange);
			
			String prefix = "sources_within_"+oDF.format(radius)+"km";
			
			PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 800, true, true, false);
			
			gp.drawGraphPanel(cmlPlot, false, true, mfdXRange, mfdYRange);
			
			PlotUtils.writePlots(resourcesDir, prefix+"_cml", gp, 800, 800, true, true, false);
			
			table.addLine("![Incremental]("+resourcesDir.getName()+"/"+prefix+".png)",
					"![Cumulative]("+resourcesDir.getName()+"/"+prefix+"_cml.png)");
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		
		double[] periods = {0d, 0.2d, 1d, 5d};
		
		ReturnPeriods[] rps = {ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50};
		
		ScalarIMR crustalGMM = AttenRelRef.USGS_PRVI_ACTIVE.get();
		ScalarIMR interfaceGMM = AttenRelRef.USGS_PRVI_INTERFACE.get();
		ScalarIMR slabGMM = AttenRelRef.USGS_PRVI_SLAB.get();
		
		site.addParameterList(crustalGMM.getSiteParams());
		
		Map<TectonicRegionType, ScalarIMR> gmmsMap = Map.of(
				TectonicRegionType.ACTIVE_SHALLOW, crustalGMM,
				TectonicRegionType.SUBDUCTION_INTERFACE, interfaceGMM,
				TectonicRegionType.SUBDUCTION_SLAB, slabGMM);
		
		for (ScalarIMR gmm : gmmsMap.values()) {
			gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
			gmm.getParameter(SigmaTruncLevelParam.NAME).setValue(3d);
		}
		
		IMT_Info imtInfo = new IMT_Info();
		
		Range curveXRange = new Range(1e-2, 1e1);
		Range curveYRange = new Range(1e-6, 1e0);
		
		EvenlyDiscretizedFunc disaggMagRange = SolSiteHazardCalc.disaggRange(5d, 9.1d, 0.5, false);
		EvenlyDiscretizedFunc disaggDistRange = SolSiteHazardCalc.disaggRange(10, 300, 20, true);
		
		for (int p=0; p<periods.length; p++) {
			String perLabel, perPrefix;
			
			DiscretizedFunc xVals;
			if (periods[p] == 0d) {
				xVals = imtInfo.getDefaultHazardCurve(PGA_Param.NAME);
				perLabel = "PGA";
				perPrefix = "pga";
				for (ScalarIMR gmm : gmmsMap.values())
					gmm.setIntensityMeasure(PGA_Param.NAME);
			} else {
				xVals = imtInfo.getDefaultHazardCurve(SA_Param.NAME);
				perLabel = oDF.format(periods[p])+"s SA";
				perPrefix = oDF.format(periods[p])+"s";
				for (ScalarIMR gmm : gmmsMap.values()) {
					gmm.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), periods[p]);
				}
			}
			DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : xVals)
				logXVals.set(Math.log(pt.getX()), pt.getY());
			
			lines.add("## "+perLabel+" Hazard Contributions");
			lines.add(topLink); lines.add("");
			
			TectonicRegionType[] plotTRTs = {
					null,
					TectonicRegionType.ACTIVE_SHALLOW,
					TectonicRegionType.SUBDUCTION_INTERFACE,
					TectonicRegionType.SUBDUCTION_SLAB
			};
			
			List<XY_DataSet> curveFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> curveChars = new ArrayList<>();
			
			DiscretizedFunc newTotalCurve = null;
			DiscretizedFunc prevTotalCurve = null;
			
			DisaggregationCalculator disagg = new DisaggregationCalculator();
			disagg.setMagRange(disaggMagRange.getMinX(), disaggMagRange.size(), disaggMagRange.getDelta());
			disagg.setDistanceRange(disaggDistRange.getMinX(), disaggDistRange.size(), disaggDistRange.getDelta());
			
			List<DisaggregationPlotData[]> disaggPlotDatas = new ArrayList<>();
			
			for (TectonicRegionType plotTRT : plotTRTs) {
				AbstractERF calcERF, calcPrevERF;
				Color color;
				float thickness;
				String type;
				if (plotTRT == null) {
					calcERF = erf;
					calcPrevERF = prevERF;
					color = Color.BLACK;
					thickness = 4f;
					type = "Full";
				} else {
					calcERF = new TRTFilteredERF(erf, plotTRT);
					calcPrevERF = new TRTFilteredERF(prevERF, plotTRT);
					color = trtColors.get(plotTRT);
					thickness = 2f;
					switch (plotTRT) {
					case ACTIVE_SHALLOW:
						type = "Crustal";
						break;
					case SUBDUCTION_INTERFACE:
						type = "Interface";
						break;
					case SUBDUCTION_SLAB:
						type = "Intraslab";
						break;

					default:
						throw new IllegalStateException();
					}
				}
				calc.getHazardCurve(logXVals, site, gmmsMap, calcPrevERF);
				DiscretizedFunc prevCurve = xVals.deepClone();
				for (int i=0; i<xVals.size(); i++)
					prevCurve.set(i, logXVals.getY(i));
				
				calc.getHazardCurve(logXVals, site, gmmsMap, calcERF);
				DiscretizedFunc newCurve = xVals.deepClone();
				for (int i=0; i<xVals.size(); i++)
					newCurve.set(i, logXVals.getY(i));
				
				newCurve.setName(type+" (2025)");
				curveFuncs.add(newCurve);
				curveChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, color));
				
				prevCurve.setName(type+" (2003)");
				curveFuncs.add(prevCurve);
				curveChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, thickness, color));
				
				if (plotTRT == null) {
					newTotalCurve = newCurve;
					prevTotalCurve = prevCurve;
				}
				
				DisaggregationPlotData[] disaggDatas = new DisaggregationPlotData[rps.length];
				for (int r=0; r<rps.length; r++) {
					ReturnPeriods rp = rps[r];
					double iml = newTotalCurve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
					disagg.disaggregate(Math.log(iml), site, gmmsMap, calcERF, calc.getSourceFilters(), calc.getAdjustableParams());
					
					disaggDatas[r] = disagg.getDisaggPlotData();
					
					if (plotTRT != null) {
						// scale to my contribution
						double myProb = newCurve.getInterpolatedY_inLogXLogYDomain(iml);
						double scale = myProb/rp.oneYearProb;
						double[][][] pdf3D = disaggDatas[r].getPdf3D();
						for (int d=0; d<pdf3D.length; d++)
							for (int m=0; m<pdf3D[d].length; m++)
								for (int e=0; e<pdf3D[d][m].length; e++)
									pdf3D[d][m][e] *= scale;
					}
				}
				
				disaggPlotDatas.add(disaggDatas);
			}
			
			List<XYAnnotation> anns = SolSiteHazardCalc.addRPAnnotations(curveFuncs, curveChars, curveXRange, curveYRange, rps, false);
			
			PlotSpec plot = new PlotSpec(curveFuncs, curveChars, siteName,
					perLabel+" (g)", "Annual Probability of Exceedance");
			plot.setLegendInset(true);
			plot.setPlotAnnotations(anns);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			
			gp.drawGraphPanel(plot, true, true, curveXRange, curveYRange);
			
			String prefix = "hazard_curve_"+perPrefix;
			
			PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 800, true, true, false);
			
			lines.add("![Hazard Curves]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			for (int r=0; r<rps.length; r++) {
				ReturnPeriods rp = rps[r];
				lines.add("### "+perLabel+" "+rp.label+" disaggregations");
				lines.add(topLink); lines.add("");
				
				double iml = newTotalCurve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
				double prevIML = prevTotalCurve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
				
				table = MarkdownUtils.tableBuilder();
				table.addLine("", "2025", "2003");
				table.addLine(rp.label, (float)iml, (float)prevIML);
				
				lines.addAll(table.build());
				lines.add("");
				
				for (int t=0; t<plotTRTs.length; t++) {
					TectonicRegionType trt = plotTRTs[t];
					String disaggPrefix = prefix+"_"+rp.name()+"_disagg";
					if (trt != null)
						disaggPrefix += "_"+trt.name();
					
					DisaggregationPlotData plotData = disaggPlotDatas.get(t)[r];
						
					PureJavaDisaggPlotter.writeChartPlot(resourcesDir, disaggPrefix,
							plotData, 800, 800, true, false);
					
					String type;
					if (trt == null) {
						type = "Total";
					} else {
						switch (trt) {
						case ACTIVE_SHALLOW:
							type = "Crustal";
							break;
						case SUBDUCTION_INTERFACE:
							type = "Interface";
							break;
						case SUBDUCTION_SLAB:
							type = "Intraslab";
							break;

						default:
							throw new IllegalStateException();
						}
					}
					
					if (trt != null) {
						lines.add("#### "+perLabel+" "+rp.label+" "+type+" Disaggregation");
						lines.add(topLink); lines.add("");
					}
					
					lines.add("![Disagg plot]("+resourcesDir.getName()+"/"+disaggPrefix+".png)");
					lines.add("");
					
					double[][][] pdf = plotData.getPdf3D();
					
					double trtSumContrib = 0d;
					double[] trtEpsilonContribs = new double[plotData.getNUM_E()];
					for (int d=0; d<disaggDistRange.size(); d++) {
						for (int m=0; m<disaggMagRange.size(); m++) {
							for (int e=0; e<trtEpsilonContribs.length; e++) {
								trtSumContrib += pdf[d][m][e];
								trtEpsilonContribs[e] += pdf[d][m][e];
							}
						}
					}
					
					table = MarkdownUtils.tableBuilder().initNewLine();
					if (trt != null)
						table.addColumn(type+" Contribution");
					Preconditions.checkState(plotData.getNUM_E() == EpsilonCategories.values().length);
					for (EpsilonCategories cat : EpsilonCategories.values())
						table.addColumn(cat.label);
					table.finalizeLine().initNewLine();
					if (trt != null)
						table.addColumn(pDF.format(trtSumContrib/100d));
					for (double trtEpsilonContrib : trtEpsilonContribs)
						table.addColumn(pDF.format(trtEpsilonContrib/100d));
					table.finalizeLine();
					lines.add("__Epsilon breakdown__");
					lines.add("");
					lines.addAll(table.build());
					lines.add("");
					
					if (trt != null ) {
						ScalarIMR gmm = gmmsMap.get(trt);
						
						Map<int[], Double> indexToContribMap = new HashMap<>();
						for (int d=0; d<disaggDistRange.size(); d++) {
							for (int m=0; m<disaggMagRange.size(); m++) {
								double contrib = 0d;
								for (double val : pdf[d][m])
									contrib += val;
								if (contrib > 0d)
									indexToContribMap.put(new int[] {d,m}, contrib);
							}
						}
						List<int[]> sorted = ComparablePairing.getSortedData(indexToContribMap);
						Collections.reverse(sorted);
						
						table = MarkdownUtils.tableBuilder();
						table.addLine("Magnitude", "Distance (km)", "Contribution %", "Max Contributing Rupture",
								"Rup Median GM", "Rup Epsilon", "Rup P(GM>"+oDF.format(iml)+")");
						
						List<List<ProbEqkRupture>> distBinnedRups = new ArrayList<>();
						for (int d=0; d<disaggDistRange.size(); d++)
							distBinnedRups.add(new ArrayList<>());
						double maxDist = disaggDistRange.getMaxX()+0.5*disaggDistRange.getDelta();
						for (ProbEqkSource source : disaggIndvERF) {
							if (source.getTectonicRegionType() == trt) {
								double sourceDist = source.getMinDistance(site);
								if (sourceDist < maxDist) {
									for (ProbEqkRupture rup : source) {
										double rupDist = rup.getRuptureSurface().getDistanceRup(site.getLocation());
										if (rupDist < maxDist) {
											distBinnedRups.get(disaggDistRange.getClosestXIndex(rupDist)).add(rup);
										}
									}
								}
							}
						}
						
						for (int i=0; i<sorted.size() && i<10; i++) {
							int[] indexes = sorted.get(i);
							int d = indexes[0];
							int m = indexes[1];
							
							double sumContrib = 0d;
							for (double val : pdf[d][m])
								sumContrib += val;
							
							if (sumContrib == 0d)
								break;
							
							table.initNewLine();
							
							double middleMag = disaggMagRange.getX(m);
							double lowerMag = middleMag - 0.5*disaggMagRange.getDelta();
							double upperMag = middleMag + 0.5*disaggMagRange.getDelta();
							
							table.addColumn("M"+oDF.format(lowerMag)+"-"+oDF.format(upperMag));
							
							double middleDist = disaggDistRange.getX(d);
							double lowerDist = middleDist - 0.5*disaggDistRange.getDelta();
							double upperDist = middleDist + 0.5*disaggDistRange.getDelta();
							
							table.addColumn(oDF.format(lowerDist)+"-"+oDF.format(upperDist)+" km");
							table.addColumn(pDF.format(sumContrib/100d));
							
							double maxContrib = 0d;
							ProbEqkRupture maxContribRup = null;
							double maxContribMedian = Double.NaN;
							double maxContribEpsilon = Double.NaN;
							double maxContribPOE = Double.NaN;
							for (ProbEqkRupture rup : distBinnedRups.get(d)) {
								if ((float)rup.getMag() >= (float)lowerMag && (float)rup.getMag() < (float)upperMag) {
									gmm.setEqkRupture(rup);
									double poe = gmm.getExceedProbability(Math.log(iml));
									double contrib = poe*rup.getMeanAnnualRate(1d);
									if (contrib > maxContrib) {
										maxContrib = contrib;
										maxContribMedian = Math.exp(gmm.getMean());
										maxContribEpsilon = gmm.getEpsilon(Math.log(iml));
										maxContribPOE = poe;
										maxContribRup = rup;
									}
								}
							}
							
							table.addColumn("M"+oDF.format(maxContribRup.getMag())+" @ "
									+(int)(maxContribRup.getRuptureSurface().getDistanceRup(site.getLocation())+0.5)+" km");
							table.addColumn(oDF.format(maxContribMedian));
							table.addColumn(oDF.format(maxContribEpsilon));
							table.addColumn(pDF.format(maxContribPOE));
							
							table.finalizeLine();
						}
						
						lines.add("__Top "+type+" contributors__");
						lines.add("");
						lines.addAll(table.build());
						lines.add("");
					}
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
		System.out.println("DONE with "+siteName);
	}
	
	private static Map<TectonicRegionType, IncrementalMagFreqDist> calcTRT_MFDs(
			AbstractERF erf, Site site, double radius, EvenlyDiscretizedFunc refMFD) {
		IncrementalMagFreqDist activeMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		IncrementalMagFreqDist interfaceMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		IncrementalMagFreqDist slabMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		
		for (ProbEqkSource source : erf) {
			double sourceDist = source.getMinDistance(site);
			if (sourceDist < radius+200) {
				TectonicRegionType trt = source.getTectonicRegionType();
				for (ProbEqkRupture rup : source) {
					double dist = rup.getRuptureSurface().getDistanceRup(site.getLocation());
					if (dist < radius) {
						double mag = rup.getMag();
						double rate = rup.getMeanAnnualRate(1d);
						int magIndex = refMFD.getClosestXIndex(mag);
						switch (trt) {
						case ACTIVE_SHALLOW:
							activeMFD.add(magIndex, rate);
							break;
						case SUBDUCTION_INTERFACE:
							interfaceMFD.add(magIndex, rate);
							break;
						case SUBDUCTION_SLAB:
							slabMFD.add(magIndex, rate);
							break;

						default:
							throw new IllegalStateException("Unexpected tectonic region type: "+trt);
						}
					}
				}
			}
		}
		
		return Map.of(TectonicRegionType.ACTIVE_SHALLOW, activeMFD,
				TectonicRegionType.SUBDUCTION_INTERFACE, interfaceMFD,
				TectonicRegionType.SUBDUCTION_SLAB, slabMFD);
	}
	
	private static class TRTFilteredERF extends AbstractERF {
		
		private List<ProbEqkSource> sources;
		
		public TRTFilteredERF(AbstractERF erf, TectonicRegionType trt) {
			sources = new ArrayList<>();
			for (ProbEqkSource source : erf)
				if (source.getTectonicRegionType() == trt)
					sources.add(source);
		}

		@Override
		public int getNumSources() {
			return sources.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return sources.get(idx);
		}

		@Override
		public void updateForecast() {
			
		}

		@Override
		public String getName() {
			return "Filtered";
		}
		
	}
	

}
