package scratch.kevin.nshm23.hazardValidation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.Precision;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueListList;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.FixedDistanceCutoffFilter;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupCartoonGenerator;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.util.SiteTranslator;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import gov.usgs.earthquake.nshmp.gmm.Gmm;
import gov.usgs.earthquake.nshmp.gmm.GmmInput;
import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import net.mahdilamb.colormap.Colors;
import scratch.kevin.nshm23.hazardValidation.WrapperRupSetRupMapper.WrapperMatch;

public class ComparisonPageGen {
	
	public static void main(String[] args) throws IOException {
//		String imtName = "PGA";
//		String imtDir = "PGA";
//		double period = 0d;
		
		String imtName = "1s SA";
		String imtDir = "SA1P0";
		double period = 1d;
		
		double closestDist = Double.POSITIVE_INFINITY;
		for (Location loc : NSHM23_RegionLoader.GridSystemRegions.CEUS_STABLE.load().getBorder())
			if (loc.lat < 49d)
				closestDist = Math.min(closestDist, LocationUtils.horzDistanceFast(loc, new Location(loc.lat, -115)));
		System.out.println("Closest CEUS grid node to the -115 boundary: "+(int)closestDist+" km");
		
//		ReturnPeriods[] rps = ReturnPeriods.values();
		ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50 };
		
		EnumSet<TectonicRegionType> trts = EnumSet.of(
				TectonicRegionType.ACTIVE_SHALLOW
				, TectonicRegionType.STABLE_SHALLOW
//				, TectonicRegionType.SUBDUCTION_INTERFACE
//				, TectonicRegionType.SUBDUCTION_SLAB
		);
		
		int vs30 = 760;
		
		String nameMine = "OpenSHA";
		String nameTheirs = "NSHMP-Haz";
		String suffixMine = "opensha";
		String suffixTheirs = "ext";
		
		File fssDir = new File("/home/kevin/OpenSHA/fss_inversions");
		File nshm23Dir = new File(fssDir, "2024_02_02-nshm23_branches-WUS_FM_v3/");
		
		EnumSet<IncludeBackgroundOption> bgOps = EnumSet.allOf(IncludeBackgroundOption.class);
//		EnumSet<IncludeBackgroundOption> bgOps = EnumSet.of(IncludeBackgroundOption.EXCLUDE);
//		EnumSet<IncludeBackgroundOption> bgOps = EnumSet.of(IncludeBackgroundOption.ONLY);
//		File inputDir = new File(fssDir, "2026_07_30-nshm23-hazard_validation-WUS");
//		File inputDir = new File(fssDir, "2026_07_30-nshm23-hazard_validation-WUS-maxDist300");
		File inputDir = new File(fssDir, "2026_07_30-nshm23-hazard_validation-WUS-maxDist300-nshmpIMLs");
		File inputSolFile = new File(nshm23Dir, "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip");
		boolean rupMFDs = false;

		File compDirActiveSub = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
				+ "conus_2023.R2-no_gmm_region-by_source-active_subduction-vs760-0p1-20260515-95ae7e82fbb85d");
//				+ "conus_2023.R2-no_gmm_region-by_source-all-vs760-0p1-20260515-ad25787f499c13"); // hack for all, not just active
		File compDirStable = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
				+ "conus_2023.R2-no_gmm_region-by_source-stable-vs760-0p1-20260515-ea28700aef6c6f");
		File outputDir = new File(inputDir, "nshmp_haz_comparisons_"+imtDir);
		
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
//		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0");
		
		boolean doWrapperCalc = false;
		
		GriddedRegion mapReg = GriddedRegion.fromFeature(Feature.read(new File(inputDir, "gridded_region.geojson")));
		boolean plotTraces = false;
		
		Region zoomReg = null;

		Range curveXRange = new Range(1e-3, 1e1);
		Range curveYRange = new Range(1e-6, 1e-0);
//		Range curveXRange = new Range(1e-2, 1e1);
//		Range curveYRange = new Range(1e-6, 1e-2);
		
		// the bay area paleo issue
//		zoomReg = new Region(new Location(37.5, -122.5), new Location(38.5, -121.5));
		
		// the rake-bin issue in NE CA 
//		curveXRange = new Range(1e-1, 1e1);
//		curveYRange = new Range(1e-5, 2e-3);
//		zoomReg = new Region(new Location(40, -122), new Location(41, -120));
		
		// the dip issue near Santa Barbara
//		curveXRange = new Range(1e-1, 1e1);
//		curveYRange = new Range(1e-5, 2e-3);
//		zoomReg = new Region(new Location(34, -121), new Location(35, -118));
		
		// gridded degugging in Montana
//		zoomReg = new Region(new Location(45, -110), new Location(49, -105));
		
		// wider gridded degugging in NE corner
//		zoomReg = new Region(new Location(40, -115), new Location(49, -105));
		
		// Nevada-ish
		zoomReg = new Region(new Location(35, -121), new Location(43, -113));
		
		SourceFilterManager sourceFilters = new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS);
//		SourceFilterManager sourceFilters = new SourceFilterManager(SourceFilters.FIXED_DIST_CUTOFF);
//		sourceFilters.getFilterInstance(FixedDistanceCutoffFilter.class).setMaxDistance(300d);
		
		Map<TectonicRegionType, AttenRelRef> gmmRefs = Map.of(
				TectonicRegionType.ACTIVE_SHALLOW, AttenRelRef.USGS_NSHM23_ACTIVE,
				TectonicRegionType.STABLE_SHALLOW, AttenRelRef.USGS_NSHM23_STABLE);
//				TectonicRegionType.ACTIVE_SHALLOW, AttenRelRef.NGAWest_2014_AVG_NOIDRISS);
		Map<TectonicRegionType, ScalarIMR> gmms = new HashMap<>();
		Map<TectonicRegionType, ScalarIMR> wrapperGMMs = gmms;
		for (TectonicRegionType trt : gmmRefs.keySet()) {
			ScalarIMR gmm = gmmRefs.get(trt).get();
			gmm.setParamDefaults();
			gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
			gmm.getParameter(SigmaTruncLevelParam.NAME).setValue(3d);
			if (period == 0d) {
				gmm.setIntensityMeasure(PGA_Param.NAME);
			} else {
				gmm.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
			}
			gmms.put(trt, gmm);
		}
		
//		NSHMP_GMM_Wrapper mixedForWrapperGMM = new NSHMP_GMM_Wrapper.WeightedCombination(
//				Map.of(Gmm.TOTAL_TREE_CONUS_ACTIVE_CRUST_2023, 2d/3d,
//						Gmm.TOTAL_TREE_CONUS_STABLE_CRUST_2023, 1d/3d), "Hack Mixed GMM", "HackMixed", false);
//		mixedForWrapperGMM.setParamDefaults();
//		mixedForWrapperGMM.getParameter(SigmaTruncLevelParam.NAME).setValue(3d);
//		if (period == 0d) {
//			mixedForWrapperGMM.setIntensityMeasure(PGA_Param.NAME);
//		} else {
//			mixedForWrapperGMM.setIntensityMeasure(SA_Param.NAME);
//			SA_Param.setPeriodInSA_Param(mixedForWrapperGMM.getIntensityMeasure(), period);
//		}
//		wrapperGMMs = Map.of(TectonicRegionType.ACTIVE_SHALLOW, gmms.get(TectonicRegionType.ACTIVE_SHALLOW),
//				TectonicRegionType.STABLE_SHALLOW, mixedForWrapperGMM);
		
		File sourcesDirActiveSub = new File(compDirActiveSub, "vs30-"+vs30+"/"+imtDir+"/source");
		File sourcesDirStable = new File(compDirStable, "vs30-"+vs30+"/"+imtDir+"/source");
		boolean convertToProb = true;

		System.out.println("Output dir: "+outputDir.getAbsolutePath());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		List<Site> sites = SolHazardMapCalc.loadSites(mapReg, gmmRefs);
		
		Map<IncludeBackgroundOption, DiscretizedFunc[]> myCurvesMap = new EnumMap<>(IncludeBackgroundOption.class);
		for (IncludeBackgroundOption bgOp : bgOps)
			myCurvesMap.put(bgOp, loadRegularCurves(inputDir, bgOp, mapReg, period));
		
		Map<IncludeBackgroundOption, DiscretizedFunc[]> extCurvesMap = new EnumMap<>(IncludeBackgroundOption.class);
		boolean loadFault = bgOps.contains(IncludeBackgroundOption.INCLUDE) || bgOps.contains(IncludeBackgroundOption.EXCLUDE);
		boolean loadGrid = bgOps.contains(IncludeBackgroundOption.INCLUDE) || bgOps.contains(IncludeBackgroundOption.ONLY);
		
		for (TectonicRegionType trt : trts) {
			File sourceDir = trt == TectonicRegionType.STABLE_SHALLOW ? sourcesDirStable : sourcesDirActiveSub;
			if (trt == TectonicRegionType.ACTIVE_SHALLOW || trt == TectonicRegionType.STABLE_SHALLOW) {
				if (loadFault) {
					DiscretizedFunc[] curves = loadExtCurves(new File(sourceDir, "FAULT/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.EXCLUDE);
					if (trt == TectonicRegionType.STABLE_SHALLOW) {
						// "fault cluster" is interface in the active_subduction dir, but crustal for stable
						curves = loadExtCurves(new File(sourceDir, "FAULT_CLUSTER/curves.csv"), mapReg);
						addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.EXCLUDE);
					}
					curves = loadExtCurves(new File(sourceDir, "FAULT_SYSTEM/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.EXCLUDE);
				}
				if (loadGrid) {
					DiscretizedFunc[] curves = loadExtCurves(new File(sourceDir, "GRID/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.ONLY);
					// TODO: is zone grid or fault?
					curves = loadExtCurves(new File(sourceDir, "ZONE/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.ONLY);
				}
			} else if (trt == TectonicRegionType.SUBDUCTION_INTERFACE) {
				if (loadFault) {
					DiscretizedFunc[] curves = loadExtCurves(new File(sourceDir, "INTERFACE/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.EXCLUDE);
					// "fault cluster" is interface in the active_subduction dir, but crustal for stable
					curves = loadExtCurves(new File(sourceDir, "FAULT_CLUSTER/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.EXCLUDE);
				}
			} else if (trt == TectonicRegionType.SUBDUCTION_SLAB) {
				if (loadGrid) {
					DiscretizedFunc[] curves = loadExtCurves(new File(sourceDir, "SLAB/curves.csv"), mapReg);
					addTo(extCurvesMap, curves, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.ONLY);
				}
			}
		}
		
		// convert to probabilities
		if (convertToProb) {
			for (IncludeBackgroundOption bgOp : List.copyOf(extCurvesMap.keySet()))
				extCurvesMap.put(bgOp, ratesToProbs(extCurvesMap.get(bgOp)));
		}
		
		FaultSystemSolution sol = null;
		if (plotTraces || doWrapperCalc || zoomReg != null)
			sol = FaultSystemSolution.load(inputSolFile);
		GeographicMapMaker fullMapMaker = new GeographicMapMaker(mapReg);
		GeographicMapMaker zoomMapMaker = null;
		if (zoomReg != null) {
			zoomMapMaker = new GeographicMapMaker(zoomReg);
			zoomMapMaker.setFaultSections(sol.getRupSet().getFaultSectionDataList());
			zoomMapMaker.setSectOutlineChar(null);
		}
		if (plotTraces) {
			fullMapMaker.setFaultSections(sol.getRupSet().getFaultSectionDataList());
			fullMapMaker.setSectOutlineChar(null);
		}
//		mapMaker.setDefaultPlotWidth(1000);
		
		BaseFaultSystemSolutionERF solERF = new BaseFaultSystemSolutionERF();
		if (doWrapperCalc) {
			solERF.setSolution(sol);
			solERF.setGriddedSeismicitySettings(solERF.getGriddedSeismicitySettings().forSupersamplingSettings(GridCellSupersamplingSettings.QUICK));
			if (sol.hasModule(RupMFDsModule.class))
				solERF.setParameter(UseRupMFDsParam.NAME, rupMFDs);
			solERF.getTimeSpan().setDuration(1d);
		}
		
		WrapperMatch[] wrapperMappings = null;
		
		Color transparent = new Color(255, 255, 255, 0);
		
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		hazCPT.setLog10(true);
		hazCPT.setNanColor(transparent);
		
//		CPT pDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance(), 10d, 50d);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-10d, 10d);
		pDiffCPT.setNanColor(transparent);
		
		HazardModel model = null;
		if (doWrapperCalc)
			model = HazardModel.load(modelDir.toPath());
		
		double diffScale;
		if (period == 0d)
			diffScale = 0.1;
		else
			diffScale = 0.05;
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-diffScale, diffScale);
		diffCPT.setNanColor(transparent);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Hazard Comparisons, "+nameMine+" vs "+nameTheirs);
		lines.add("");
		
		lines.add("This page compares "+nameMine+" and "+nameTheirs+" hazard maps.");
		lines.add("");
		
		IncludeBackgroundOption[] bgOrder = {
				IncludeBackgroundOption.INCLUDE,
				IncludeBackgroundOption.EXCLUDE,
				IncludeBackgroundOption.ONLY
		};
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		for (IncludeBackgroundOption bgOp : bgOrder) {
			if (!bgOps.contains(bgOp))
				continue;
			DiscretizedFunc[] myCurves = myCurvesMap.get(bgOp);
			DiscretizedFunc[] extCurves = extCurvesMap.get(bgOp);
			
			NshmErf wrapperERF = null;
			if (doWrapperCalc) {
				wrapperERF = new NshmErf(model, trts, bgOp);
				wrapperERF.getTimeSpan().setDuration(1d);
				wrapperERF.updateForecast();
				
				solERF.setParameter(IncludeBackgroundParam.NAME, bgOp);
				solERF.updateForecast();
			}
			
			String mapLabelAdd;
			switch (bgOp) {
			case INCLUDE:
				lines.add("## Total hazard (fault+gridded)");
				mapLabelAdd = "";
				break;
			case EXCLUDE:
				lines.add("## On-fault hazard");
				mapLabelAdd = "On-fault ";
				break;
			case ONLY:
				lines.add("## Gridded hazard");
				mapLabelAdd = "Gridded ";
				break;
			default:
				throw new IllegalArgumentException("Unexpected value: " + bgOp);
			}
			lines.add(topLink); lines.add("");
			
			for (ReturnPeriods rp : rps) {
				GriddedGeoDataSet myMap = curvestoMap(myCurves, mapReg, rp);
				GriddedGeoDataSet extMap = curvestoMap(extCurves, mapReg, rp);
				
				System.out.println("Doing "+bgOp+", "+rp);
				
				lines.add("### "+rp.label);
				lines.add(topLink); lines.add("");
				
				boolean[] zooms = zoomReg == null ? new boolean[] {false} : new boolean[] {false,true};
				
				for (boolean zoom : zooms) {
					String hazLabel = imtName+", "+rp.label;
					String prefix = bgOp.name()+"_"+rp.name();
					
					boolean curves = zoom || zooms.length == 1;
					
					GeographicMapMaker mapMaker = zoom ? zoomMapMaker : fullMapMaker;
					Region plotReg = zoom ? zoomReg : mapReg;
					if (zoom) {
						prefix += "_zoom";
						
						lines.add("#### Zoomed region");
						lines.add(topLink); lines.add("");
					}
					
					TableBuilder table = MarkdownUtils.tableBuilder();
					
					table.addLine(nameMine, nameTheirs);
					
					table.initNewLine();
					
					mapMaker.plotXYZData(myMap, hazCPT, mapLabelAdd+nameMine+", "+hazLabel+" (g)");
					mapMaker.plot(resourcesDir, prefix+"_"+suffixMine, " ");
					table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_"+suffixMine+".png)");
					mapMaker.plotXYZData(extMap, hazCPT, mapLabelAdd+nameTheirs+", "+hazLabel+" (g)");
					mapMaker.plot(resourcesDir, prefix+"_"+suffixTheirs, " ");
					table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_"+suffixTheirs+".png)");
					
					table.finalizeLine();
					
					table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
					
					GriddedGeoDataSet pDiff = mapPDiff(myMap, extMap);
					GriddedGeoDataSet diff = mapDiff(myMap, extMap);
					
					table.initNewLine();
					
					mapMaker.plotXYZData(pDiff, pDiffCPT, mapLabelAdd+nameMine+" vs "+nameTheirs+", % Change, "+hazLabel);
					mapMaker.plot(resourcesDir, prefix+"_pDiff", " ");
					table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_pDiff.png)");
					mapMaker.plotXYZData(diff, diffCPT, mapLabelAdd+nameMine+" - "+nameTheirs+", "+hazLabel+" (g)");
					mapMaker.plot(resourcesDir, prefix+"_diff", " ");
					table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_diff.png)");
					
					table.finalizeLine();
					table.addLine(diffStr(pDiff, true, plotReg), diffStr(diff, false, plotReg));
					
					lines.addAll(table.build());
					lines.add("");
					
					if (curves) {
						double maxDiff = Double.NEGATIVE_INFINITY;
						double minDiff = Double.POSITIVE_INFINITY;
						int maxDiffIndex = -1;
						int minDiffIndex = -1;
						for (int i=0; i<diff.size(); i++) {
							if (!plotReg.contains(diff.getLocation(i)))
								continue;
							double v = diff.get(i);
							if (Double.isFinite(v)) {
								if (v > maxDiff) {
									maxDiff = v;
									maxDiffIndex = i;
								}
								if (v < minDiff) {
									minDiff = v;
									minDiffIndex = i;
								}
							}
						}
						
						table = MarkdownUtils.tableBuilder();
						
						table.addLine("Min difference", "Max difference");
						
						table.initNewLine();
						Site minSite = sites.get(minDiffIndex);
						Site maxSite = sites.get(maxDiffIndex);
						HazardCurveCalculator calc = new HazardCurveCalculator(sourceFilters);
						for (boolean min : new boolean[] {true,false}) {
							int index = min ? minDiffIndex : maxDiffIndex;
							Preconditions.checkState(index >= 0);
							DiscretizedFunc myCurve = myCurves[index];
							DiscretizedFunc extCurve = extCurves[index];
							Location loc = mapReg.getLocation(index);
							
							List<XY_DataSet> funcs = new ArrayList<>();
							List<PlotCurveCharacterstics> chars = new ArrayList<>();
							
							myCurve.setName(nameMine);
							extCurve.setName(nameTheirs);
							
							funcs.add(extCurve);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
							
							Site site = min ? minSite : maxSite;
							
							if (doWrapperCalc) {
								DiscretizedFunc wrappedCurve = extCurve.deepClone();
								DiscretizedFunc logCurve = new ArbitrarilyDiscretizedFunc();
								for (int i=0; i<wrappedCurve.size(); i++)
									logCurve.set(Math.log(wrappedCurve.getX(i)), 0d);
								
								
								calc.getHazardCurve(logCurve, site, wrapperGMMs, wrapperERF);
								
								for (int i=0; i<wrappedCurve.size(); i++)
									wrappedCurve.set(i, logCurve.getY(i));
								
								wrappedCurve.setName("Wrapper");
								funcs.add(wrappedCurve);
								chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
								
							}
							
							funcs.add(myCurve);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
							
							if (doWrapperCalc) {
								DiscretizedFunc fssCurve = myCurve.deepClone();
								DiscretizedFunc logCurve = new ArbitrarilyDiscretizedFunc();
								for (int i=0; i<fssCurve.size(); i++)
									logCurve.set(Math.log(fssCurve.getX(i)), 0d);
								
								
								calc.getHazardCurve(logCurve, site, gmms, solERF);
								
								for (int i=0; i<fssCurve.size(); i++)
									fssCurve.set(i, logCurve.getY(i));
								
								fssCurve.setName(nameMine+" (recalc)");
								funcs.add(fssCurve);
								chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_lightblue));
								
							}
							
							PlotSpec plot = new PlotSpec(funcs, chars, (float)loc.lat+", "+(float)loc.lon, hazLabel, "Annual Probability of Exceedance");
							plot.setLegendInset(true);
							
							HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
							
							gp.drawGraphPanel(plot, true, true, curveXRange, curveYRange);
							
							String curvePrefix = min ? prefix+"_curve_min" : prefix+"_curve_max";
							
							PlotUtils.writePlots(resourcesDir, curvePrefix, gp, 800, 800, true, true, false);
							
							table.addColumn("![Curve]("+resourcesDir.getName()+"/"+curvePrefix+".png)");
						}
						table.finalizeLine().initNewLine();
						for (boolean min : new boolean[] {true,false}) {
							int index = min ? minDiffIndex : maxDiffIndex;
							Location loc = mapReg.getLocation(index);
							String str = "Site: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude();
							Site site = min ? minSite : maxSite;
							for (Parameter<?> param : site)
								str += "<p>"+param.getName()+": "+param.getValue();
							table.addColumn(str);
						}
						table.finalizeLine();
						
						if (doWrapperCalc && bgOp == IncludeBackgroundOption.EXCLUDE) {
							// add distance hists
							double maxDist = 50d;
							double maxQuickDist = 60d;
							int bins = 50;
							
							if (wrapperMappings == null)
								wrapperMappings = WrapperRupSetRupMapper.map(sol, wrapperERF);
							WrapperMatch[] myWrapperMappings = wrapperMappings;
							
							List<Integer> minFSSRups = getRupIndexesWithinCutoff(sol, mapReg.getLocation(minDiffIndex), maxQuickDist);
							List<ProbEqkSource> minWrapperSources = minFSSRups.stream().map(
									rupIndex->myWrapperMappings[rupIndex].wrapperSource()).collect(Collectors.toList());
							List<ProbEqkRupture> minWrapperRups = minFSSRups.stream().map(
									rupIndex->myWrapperMappings[rupIndex].wrapperRup()).collect(Collectors.toList());
							List<Integer> maxFSSRups = getRupIndexesWithinCutoff(sol, mapReg.getLocation(maxDiffIndex), maxQuickDist);
							List<ProbEqkSource> maxWrapperSources = maxFSSRups.stream().map(
									rupIndex->myWrapperMappings[rupIndex].wrapperSource()).collect(Collectors.toList());
							List<ProbEqkRupture> maxWrapperRups = maxFSSRups.stream().map(
									rupIndex->myWrapperMappings[rupIndex].wrapperRup()).collect(Collectors.toList());
							
							int minRupIndex = -1;
							int maxRupIndex = -1;
							
							System.out.println("Max rate source within "+(float)maxDist+" km");
							table.initNewLine();
							for (boolean min : new boolean[] {true,false}) {
								Site site = min ? minSite : maxSite;
								
								List<XY_DataSet> funcs = new ArrayList<>();
								List<PlotCurveCharacterstics> chars = new ArrayList<>();
								
								List<Integer> fssIndexes = min ? minFSSRups : maxFSSRups;
								List<ProbEqkRupture> wrapperRups = min ? minWrapperRups : maxWrapperRups;
								Preconditions.checkState(fssIndexes.size() == wrapperRups.size());
								
								if (fssIndexes.isEmpty()) {
									table.addColumn("");
									continue;
								}
								
								// now find the rupture with the biggest change at our gm
								int siteIndex = min ? minDiffIndex : maxDiffIndex;
								double targetIM = curveVal(myCurves[siteIndex], rp);

								double maxDiffForRateDiff = 0d;
								double maxRateDiff = 0d;
								int maxDiffRupIndex = -1;
								DiscretizedFunc maxDiffMyIMs = null;
								DiscretizedFunc maxDiffTheirIMs = null;
								ProbEqkRupture fssRup = null;
								
								for (int r=0; r<fssIndexes.size(); r++) {
									int rupIndex = fssIndexes.get(r);
									ProbEqkSource mySource = solERF.getSource(solERF.getSrcIndexForFltSysRup(rupIndex));
									Preconditions.checkState(mySource.getNumRuptures() == 1);
									ProbEqkRupture myRup = mySource.getRupture(0);
									ProbEqkRupture theirRup = wrapperRups.get(r);
									
									// mine first
									DiscretizedFunc myIMs = myCurves[siteIndex].deepClone();
									DiscretizedFunc logIMs = new ArbitrarilyDiscretizedFunc();
									for (int i=0; i<myIMs.size(); i++)
										logIMs.set(Math.log(myIMs.getX(i)), 0d);
									
									ScalarIMR gmm = gmms.size() == 1 ? gmms.values().iterator().next() : gmms.get(mySource.getTectonicRegionType());
									
									gmm.setSite(site);
									
									gmm.getExceedProbabilities(myRup, logIMs);
									for (int i=0; i<myIMs.size(); i++)
										myIMs.set(i, logIMs.getY(i));
									
									double myIM = myIMs.getInterpolatedY_inLogXLogYDomain(targetIM);
									
									// theirs first
									DiscretizedFunc theirIMs = extCurves[siteIndex].deepClone();
									logIMs = new ArbitrarilyDiscretizedFunc();
									for (int i=0; i<theirIMs.size(); i++)
										logIMs.set(Math.log(theirIMs.getX(i)), 0d);
									
									gmm.getExceedProbabilities(theirRup, logIMs);
									for (int i=0; i<theirIMs.size(); i++)
										theirIMs.set(i, logIMs.getY(i));
									
									double theirIM = theirIMs.getInterpolatedY_inLogXLogYDomain(targetIM);
									
									double myDiff = min ? theirIM - myIM : myIM - theirIM;
									
									double rateDiff = myDiff * sol.getRateForRup(rupIndex);
									if (rateDiff > maxRateDiff) {
										maxRateDiff = rateDiff;
										maxDiffMyIMs = myIMs;
										maxDiffTheirIMs = theirIMs;
										maxDiffRupIndex = rupIndex;
										fssRup = myRup;
										maxDiffForRateDiff = myDiff;
									}
								}
								
								if (min)
									minRupIndex = maxDiffRupIndex;
								else
									maxRupIndex = maxDiffRupIndex;
								
								ProbEqkRupture wrapperRup = wrapperMappings[maxDiffRupIndex].wrapperRup();
								TectonicRegionType trt = wrapperMappings[maxDiffRupIndex].wrapperSource().getTectonicRegionType();
								ScalarIMR gmm = gmms.size() == 1 ? gmms.values().iterator().next() : gmms.get(trt);
								String fssRupName = "M"+(float)fssRup.getMag()+", rake="+(float)fssRup.getAveRake()
										+", P="+(float)fssRup.getProbability();
								String wrapperRupName = "M"+(float)wrapperRup.getMag()+", rake="+(float)wrapperRup.getAveRake()
										+", P="+(float)wrapperRup.getProbability();
								
								NSHMP_GMM_Wrapper wrapGMM = gmm instanceof NSHMP_GMM_Wrapper ? (NSHMP_GMM_Wrapper)gmm : null;
								
								if (min)
									System.out.println("Min change rup for IM="+(float)targetIM+": "+maxDiffRupIndex);
								else
									System.out.println("Max change rup for IM="+(float)targetIM+": "+maxDiffRupIndex);
								System.out.println("\tDiff: "+maxDiffForRateDiff);
								System.out.println("\tRate: "+sol.getRateForRup(maxDiffRupIndex));
								System.out.println("\tRate-Diff: "+maxRateDiff);
								System.out.println("\tFSS:\t"+fssRupName);
								System.out.println("\t\tPOE:\t"+maxDiffMyIMs.getInterpolatedY_inLogXLogYDomain(targetIM));
								if (wrapGMM != null) {
									wrapGMM.setEqkRupture(fssRup);
									System.out.println("\t\tGmmInput:\t"+wrapGMM.getCurrentGmmInput());
								}
								System.out.println("\tWrapper:\t"+wrapperRupName);
								System.out.println("\t\tPOE:\t"+maxDiffTheirIMs.getInterpolatedY_inLogXLogYDomain(targetIM));
								if (wrapGMM != null) {
									wrapGMM.setEqkRupture(wrapperRup);
									System.out.println("\t\tGmmInput:\t"+wrapGMM.getCurrentGmmInput());
								}
								
								maxDiffMyIMs.setName(nameMine+" "+fssRupName);
								funcs.add(maxDiffMyIMs);
								chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
								
								maxDiffTheirIMs.setName(nameTheirs+" "+wrapperRupName);
								funcs.add(maxDiffTheirIMs);
								chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
								
								if (wrapGMM != null) {
									// extra tests
									wrapGMM.setEqkRupture(fssRup);
									GmmInput fssInput = wrapGMM.getCurrentGmmInput();
									wrapGMM.setEqkRupture(wrapperRup);
									GmmInput wrapInput = wrapGMM.getCurrentGmmInput();
									
									List<GmmInput> updates = new ArrayList<>();
									List<PlotCurveCharacterstics> updateChars = new ArrayList<>();
									List<String> updateNames = new ArrayList<>();
									
									updates.add(GmmInput.builder().fromCopy(fssInput)
											.rake(wrapInput.rake).build());
									updateNames.add("wrapper rake");
									updateChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_orange));
									
									updates.add(GmmInput.builder().fromCopy(fssInput)
											.rake(wrapInput.rake).dip(wrapInput.dip).build());
									updateNames.add("wrapper rake/dip");
									updateChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_purple));
									
									updates.add(GmmInput.builder().fromCopy(fssInput)
											.rake(wrapInput.rake).dip(wrapInput.dip).zHyp(wrapInput.zHyp).build());
									updateNames.add("wrapper rake/dip/zHyp");
									updateChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Colors.tab_olive));
									
									for (int u=0; u<updates.size(); u++) {
										wrapGMM.setCurrentGmmInput(updates.get(u));
										DiscretizedFunc modIMs = maxDiffMyIMs.deepClone();
										ArbitrarilyDiscretizedFunc logIMs = new ArbitrarilyDiscretizedFunc();
										for (int i=0; i<modIMs.size(); i++)
											logIMs.set(Math.log(modIMs.getX(i)), 0d);
										wrapGMM.getExceedProbabilities(logIMs);
										for (int i=0; i<modIMs.size(); i++)
											modIMs.set(i, logIMs.getY(i));
										modIMs.setName(nameMine+" w/ "+updateNames.get(u));
										funcs.add(modIMs);
										chars.add(updateChars.get(u));
									}
								}
								
								Range yRange = new Range(1e-5, 1e0);
								
								DefaultXY_DataSet imFunc = new DefaultXY_DataSet(targetIM, yRange.getLowerBound(), targetIM, yRange.getUpperBound());
								imFunc.setName(rp.label+": "+(float)targetIM);
								funcs.add(imFunc);
								chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY));
								
								Location loc = site.getLocation();
								PlotSpec plot = new PlotSpec(funcs, chars, (float)loc.lat+", "+(float)loc.lon, hazLabel, "Conditional Probability of Exceedance");
								plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
								
								HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
								
								gp.getPlotPrefs().setLegendFontSize(14);
								
//								gp.drawGraphPanel(plot, true, false, curveXRange, new Range(0d, 1d));
								gp.drawGraphPanel(plot, true, true, curveXRange, yRange);
								
								String curvePrefix = min ? prefix+"_cond_prob_min" : prefix+"_cond_prob_max";
								
								PlotUtils.writePlots(resourcesDir, curvePrefix, gp, 800, 800, true, true, false);
								
								table.addColumn("![Curve]("+resourcesDir.getName()+"/"+curvePrefix+".png)");
							}
							table.finalizeLine();
							
							// rupture cartoons
							if (minRupIndex >= 0 && maxRupIndex >= 0) {
								FaultSystemRupSet rupSet = sol.getRupSet();
								ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
								RupCartoonGenerator.plotRupture(resourcesDir, prefix+"_rup_min", cRups.get(minRupIndex),
										"Rup "+minRupIndex+" M"+(float)rupSet.getMagForRup(minRupIndex)
										+" rake="+(float)rupSet.getAveRakeForRup(minRupIndex)
										+" dip="+(float)rupSet.getSurfaceForRupture(minRupIndex, 1d).getAveDip(), false, true);
								RupCartoonGenerator.plotRupture(resourcesDir, prefix+"_rup_max", cRups.get(maxRupIndex),
										"Rup "+maxRupIndex+" M"+(float)rupSet.getMagForRup(maxRupIndex)
										+" rake="+(float)rupSet.getAveRakeForRup(maxRupIndex)
										+" dip="+(float)rupSet.getSurfaceForRupture(maxRupIndex, 1d).getAveDip(), false, true);
								
								table.addLine("![Min rup]("+resourcesDir.getName()+"/"+prefix+"_rup_min.png)",
										"![Max rup]("+resourcesDir.getName()+"/"+prefix+"_rup_max.png)");
							}
							
							// rake hist
							table.initNewLine();
							for (boolean min : new boolean[] {true,false}) {
								int index = min ? minDiffIndex : maxDiffIndex;
								Preconditions.checkState(index >= 0);
								Location loc = mapReg.getLocation(index);
								List<ProbEqkRupture> wrapperRups = min ? minWrapperRups : maxWrapperRups;
								List<Integer> fssIndexes = min ? minFSSRups : maxFSSRups;
								List<ProbEqkSource> fssSources = fssIndexes.stream().map(
										I->solERF.getSource(solERF.getSrcIndexForFltSysRup(I))).collect(Collectors.toList());
								
								EvenlyDiscretizedFunc rakeHistFSS = new EvenlyDiscretizedFunc(-180d, 180d, 181);
								EvenlyDiscretizedFunc rakeHistWrapper = new EvenlyDiscretizedFunc(-180d, 180d, 181);
								
								for (ProbEqkSource source : fssSources)
									for (ProbEqkRupture rup : source)
										rakeHistFSS.add(rakeHistFSS.getClosestXIndex(rup.getAveRake()), rup.getMeanAnnualRate(1d));
								for (ProbEqkRupture rup : wrapperRups)
									rakeHistWrapper.add(rakeHistWrapper.getClosestXIndex(rup.getAveRake()), rup.getMeanAnnualRate(1d));
								
								List<DiscretizedFunc> funcs = new ArrayList<>();
								List<PlotCurveCharacterstics> chars = new ArrayList<>();
								
								rakeHistWrapper.setName("Wrapper");
								funcs.add(rakeHistWrapper);
								chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, trans(Colors.tab_green, 127)));
								
								rakeHistFSS.setName(nameMine);
								funcs.add(rakeHistFSS);
								chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, trans(Colors.tab_blue, 127)));
								
								PlotSpec plot = new PlotSpec(funcs, chars, " ", "Nearby rupture rake (degrees)", "Rate");
								plot.setLegendInset(RectangleAnchor.TOP_LEFT);
								
								HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
								
								gp.drawGraphPanel(plot, false, false, new Range(-180, 180d), null);
								
								String histPrefix = prefix+"_hist_rake";
								if (min)
									histPrefix += "_min";
								else
									histPrefix += "_max";
								
								PlotUtils.writePlots(resourcesDir, histPrefix, gp, 800, 800, true, true, false);
								
								table.addColumn("![Rake hist]("+resourcesDir.getName()+"/"+histPrefix+".png)");
							}
							table.finalizeLine();
							
							for (int d=0; d<3; d++) {
								String distName;
								HistogramFunction histFSS;
								BiFunction<RuptureSurface, Location, Double> distFunc;
								if (d == 0) {
									distName = "Rrup";
									histFSS = new HistogramFunction(0d, maxDist, bins);
									distFunc = (S,L) -> S.getDistanceRup(L);
								} else if (d == 1) {
									distName = "Rjb";
									histFSS = new HistogramFunction(0d, maxDist, bins);
									distFunc = (S,L) -> S.getDistanceJB(L);
								} else {
									distName = "RX";
									histFSS = new HistogramFunction(-maxDist, maxDist, bins);
									distFunc = (S,L) -> S.getDistanceX(L);
								}
								
								EvenlyDiscretizedFunc histWrapper = histFSS.deepClone();
								
								table.initNewLine();
								for (boolean min : new boolean[] {true,false}) {
									int index = min ? minDiffIndex : maxDiffIndex;
									Preconditions.checkState(index >= 0);
									Location loc = mapReg.getLocation(index);
									List<ProbEqkSource> wrapperSources = min ? minWrapperSources : maxWrapperSources;
									List<Integer> fssIndexes = min ? minFSSRups : maxFSSRups;
									List<ProbEqkSource> fssSources = fssIndexes.stream().map(
											I->solERF.getSource(solERF.getSrcIndexForFltSysRup(I))).collect(Collectors.toList());
									
									histFSS.scale(0d);
									histWrapper.scale(0d);
									
									fillDistanceHist(fssSources, distFunc, histFSS, loc);
									fillDistanceHist(wrapperSources, distFunc, histWrapper, loc);
									
									List<DiscretizedFunc> funcs = new ArrayList<>();
									List<PlotCurveCharacterstics> chars = new ArrayList<>();
									
									histWrapper.setName(null);
									funcs.add(histWrapper);
									chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, trans(Colors.tab_green, 127)));
									
									histFSS.setName(null);
									funcs.add(histFSS);
									chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, trans(Colors.tab_blue, 127)));
									
									EvenlyDiscretizedFunc cmlFSS = new EvenlyDiscretizedFunc(histFSS.getMinX()-0.5*histFSS.getDelta(), histFSS.size(), histFSS.getDelta());
									double sum = 0d;
									for (int i=0; i<cmlFSS.size(); i++) {
										sum += histFSS.getY(i);
										cmlFSS.set(i, sum);
									}
									EvenlyDiscretizedFunc cmlWrapper = cmlFSS.deepClone();
									sum = 0d;
									for (int i=0; i<cmlWrapper.size(); i++) {
										sum += histWrapper.getY(i);
										cmlWrapper.set(i, sum);
									}
									
									cmlFSS.setName("Wrapper");
									funcs.add(cmlFSS);
									chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_green));
									
									cmlWrapper.setName(nameMine);
									funcs.add(cmlWrapper);
									chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_blue));
									
									PlotSpec plot = new PlotSpec(funcs, chars, " ", distName+" (km)", "Rate");
									plot.setLegendInset(RectangleAnchor.TOP_LEFT);
									
									HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
									
									gp.drawGraphPanel(plot, false, false, new Range(cmlFSS.getMinX(), cmlFSS.getMaxX()+cmlFSS.getDelta()), null);
									
									String histPrefix = prefix+"_hist_"+distName;
									if (min)
										histPrefix += "_min";
									else
										histPrefix += "_max";
									
									PlotUtils.writePlots(resourcesDir, histPrefix, gp, 800, 800, true, true, false);
									
									table.addColumn("!["+distName+" hist]("+resourcesDir.getName()+"/"+histPrefix+".png)");
								}
								table.finalizeLine();
							}
						}
						
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
	}
	
	private static File getHazardDir(File runDir, IncludeBackgroundOption bgOp,
			GriddedRegion gridReg) {
		String hazardPrefix = "hazard_"+(float)gridReg.getSpacing()+"deg";
		hazardPrefix += "_grid_seis_";
		hazardPrefix += bgOp.name();
		File resultsDir = new File(runDir, "results");
		return new File(resultsDir, hazardPrefix);
	}
	
	private static DiscretizedFunc[] loadRegularCurves(File runDir, IncludeBackgroundOption bgOp,
			GriddedRegion gridReg, double period) throws IOException {
		File subDir = getHazardDir(runDir, bgOp, gridReg);
		Preconditions.checkState(subDir.exists(), "Doesn't exist: %s", subDir.getAbsolutePath());
		
		File curvesFile = new File(subDir, SolHazardMapCalc.getCSV_FileName("curves", period));
		if (!curvesFile.exists())
			curvesFile = new File(curvesFile.getAbsolutePath()+".gz");
		Preconditions.checkState(curvesFile.exists(), "Doesn't exist: %s", curvesFile.getAbsolutePath());
		
		return SolHazardMapCalc.loadCurvesCSV(CSVFile.readFile(curvesFile, true), gridReg);
	}
	
	private static DiscretizedFunc[] loadExtCurves(File csvFile, GriddedRegion gridReg) throws IOException {
		System.out.println("Loading external curves from "+csvFile.getAbsolutePath());
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		double[] xVals = new double[csv.getNumCols()-2];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = csv.getDouble(0, i+2);
		
		DiscretizedFunc[] curves = new DiscretizedFunc[gridReg.getNodeCount()];
		
		int numSkipped = 0;
		int numFilled = 0;
		for (int row=1; row<csv.getNumRows(); row++) {
			double lon = csv.getDouble(row, 0);
			double lat = csv.getDouble(row, 1);
			Location loc = new Location(lat, lon);
			int index = gridReg.indexForLocation(loc);
			if (index < 0) {
				numSkipped++;
			} else {
				numFilled++;
				double[] yVals = new double[xVals.length];
				for (int i=0; i<yVals.length; i++)
					yVals[i] = csv.getDouble(row, i+2);
				curves[index] = new LightFixedXFunc(xVals, yVals);
			}
		}
		System.out.println("\tFilled in "+numFilled+"/"+gridReg.getNodeCount()+" curves (skipped "+numSkipped+")");
		
		return curves;
	}
	
	private static DiscretizedFunc[] ratesToProbs(DiscretizedFunc[] curves) {
		DiscretizedFunc[] ret = new DiscretizedFunc[curves.length];
		DiscretizedFunc curve0 = null;
		for (int i=0; i<curves.length; i++) {
			if (curves[i] != null) {
				curve0 = curves[i];
				break;
			}
		}
		Preconditions.checkNotNull(curve0, "All curves are null");
		double[] xVals = new double[curve0.size()];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = curve0.getX(i);
		for (int n=0; n<curves.length; n++) {
			DiscretizedFunc curve = curves[n];
			if (curve == null)
				continue;
			double[] probs = new double[xVals.length];
			for (int i=0; i<probs.length; i++)
				probs[i] = 1d - Math.exp(-curve.getY(i));
			ret[n] = new LightFixedXFunc(xVals, probs);
		}
		return ret;
	}
	
	private static DiscretizedFunc[] interpolateToMatch(DiscretizedFunc[] curves, DiscretizedFunc targetSpacing) {
		DiscretizedFunc[] ret = new DiscretizedFunc[curves.length];
		double[] xVals = new double[targetSpacing.size()];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = targetSpacing.getX(i);
		for (int n=0; n<curves.length; n++) {
			DiscretizedFunc curve = curves[n];
			double[] interpY = new double[xVals.length];
			for (int i=0; i<interpY.length; i++) {
				double x = xVals[i];
				if (x < curve.getMinX())
					interpY[i] = curve.getY(0);
				else if (x > curve.getMaxX())
					interpY[i] = 0d;
				else
					interpY[i] = curve.getInterpolatedY_inLogXDomain(x);
			}
			ret[n] = new LightFixedXFunc(xVals, interpY);
		}
		return ret;
	}
	
	private static void addTo(Map<IncludeBackgroundOption, DiscretizedFunc[]> map, DiscretizedFunc[] curves,
			IncludeBackgroundOption... bgOps) {
		for (IncludeBackgroundOption bgOp : bgOps)
			map.put(bgOp, addTo(map.get(bgOp), curves));
	}
	
	private static DiscretizedFunc[] addTo(DiscretizedFunc[] current, DiscretizedFunc[]... allCurves) {
//		System.out.println("addTo; current: "+(current == null ? "null" : "non-null"));
//		System.out.println("allCurves.length: "+allCurves.length);
//		for (int i=0; i<allCurves.length; i++) {
//			System.out.println("allCurves["+i+"].length: "+allCurves[i].length);
//			int numNull = 0;
//			for (DiscretizedFunc curve : allCurves[i])
//				if (curve == null)
//					numNull++;
//			System.out.println("\t"+numNull+"/"+allCurves[i].length+" are null");
//		}
		if (current == null) {
			if (allCurves.length == 0)
				return allCurves[0];
			return add(allCurves);
		} else if (allCurves.length > 1) {
			DiscretizedFunc[] comb = add(allCurves);
			return add(current, comb);
		} else {
			Preconditions.checkState(allCurves.length == 1);
			return add(current, allCurves[0]);
		}
	}
	
	private static DiscretizedFunc[] add(DiscretizedFunc[]... allCurves) {
		DiscretizedFunc[] ret = new DiscretizedFunc[allCurves[0].length];
		
		double[] xVals = new double[allCurves[0][0].size()];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = allCurves[0][0].getX(i);
		
		for (int i=0; i<ret.length; i++) {
			boolean anyNonNull = false;
			double[] yVals = new double[xVals.length];
			for (DiscretizedFunc[] curves : allCurves) {
				DiscretizedFunc curve = curves[i];
				if (curve == null)
					continue;
				else
					anyNonNull = true;
				Preconditions.checkState(curve.size() == xVals.length,
						"X-value count mismatch: %s vs %s; minX: %s vs %s; maxX: %s vs %s",
						curve.size(), xVals.length, (float)curve.getMinX(), (float)xVals[0],
						(float)curve.getMaxX(), (float)xVals[xVals.length-1]);
				for (int j=0; j<xVals.length; j++) {
					Preconditions.checkState((float)xVals[j] == (float)curve.getX(j));
					yVals[j] += curve.getY(j);
				}
			}
			if (anyNonNull)
				ret[i] = new LightFixedXFunc(xVals, yVals);
		}
		
		return ret;
	}
	
	private static GriddedGeoDataSet curvestoMap(DiscretizedFunc[] curves, GriddedRegion gridRegion, ReturnPeriods rp) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(gridRegion);
		for (int i=0; i<ret.size(); i++)
			ret.set(i, Double.NaN);
		
		for (int i=0; i<curves.length; i++)
			ret.set(i, curveVal(curves[i], rp));
		
		return ret;
	}
	
	private static GriddedGeoDataSet zerosToNaNs(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		for (int i=0; i<xyz.size(); i++)
			if (xyz.get(i) == 0d)
				xyz.set(i, Double.NaN);
		return xyz;
	}
	
	public static GriddedGeoDataSet mapPDiff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++) {
			double v1 = map1.get(i);
			double v2 = map2.get(i);
			if (Double.isNaN(v1) || Double.isNaN(v2))
				ret.set(i, Double.NaN);
			else
				ret.set(i, 100d*(v1 - v2)/v2);
		}
		return ret;
	}
	
	public static GriddedGeoDataSet mapDiff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++) {
			double v1 = map1.get(i);
			double v2 = map2.get(i);
			if (Double.isNaN(v1) || Double.isNaN(v2))
				ret.set(i, Double.NaN);
			else
				ret.set(i, v1 - v2);
		}
		return ret;
	}
	
	private static String diffStr(GriddedGeoDataSet diff, boolean isPDiff, Region plotReg) {
		int numNan = 0;
		int numInf = 0;
		double max = Double.NEGATIVE_INFINITY;
		double min = Double.POSITIVE_INFINITY;
		double maxAbs = 0d;
		double sumDiff = 0d;
		double sumAbsDiff = 0d;
		
		if (isPDiff) {
			diff = diff.copy();
			diff.scale(0.01);
		}

		List<Double> allConsidered = new ArrayList<>();
		List<Double> allConsideredAbs = new ArrayList<>();
		
		for (int i=0; i<diff.size(); i++) {
			if (!plotReg.contains(diff.getLocation(i)))
				continue;
			double v = diff.get(i);
			if (Double.isNaN(v)) {
				numNan++;
			} else if (Double.isInfinite(v)) {
				numInf++;
			} else {
				Preconditions.checkState(Double.isFinite(v), "Non-finite (but non-nan) value: %s", v);
				max = Math.max(max, v);
				min = Math.min(min, v);
				double abs = Math.abs(v);
				maxAbs = Math.max(maxAbs, abs);
				sumDiff += v;
				sumAbsDiff += abs;
				allConsidered.add(v);
				allConsideredAbs.add(abs);
			}
		}
		
		int numConsidered = allConsidered.size();
		double avg = sumDiff/(double)numConsidered;
		double avgAbs = sumAbsDiff/(double)numConsidered;

		double median = StatUtils.mean(Doubles.toArray(allConsidered));
		double medianAbs = StatUtils.mean(Doubles.toArray(allConsideredAbs));
		
		DecimalFormat df;
		if (isPDiff)
			df = new DecimalFormat("0.00%");
		else
			df = new DecimalFormat("0.000");
		
		return "Range=["+df.format(min)+", "+df.format(max)+"]; maxAbs="+df.format(maxAbs)
				+"<p>avg="+df.format(avg)+"; avgAbs="+df.format(avgAbs)
				+"<p>median="+df.format(median)+"; medianAbs="+df.format(medianAbs)
				+"<p>"+numNan+" NaN; "+numInf+" inf";
	}
	
	private static List<ProbEqkSource> getSourcesWithinCutoff(AbstractERF erf, Location loc, double cutoff) {
		Site site = new Site(loc);
		return erf.getSourceList().parallelStream().filter(S->(float)S.getMinDistance(site) <= (float)cutoff).collect(Collectors.toList());
	}
	
	private static List<Integer> getRupIndexesWithinCutoff(FaultSystemSolution sol, Location loc, double cutoff) {
		List<Integer> ret = new ArrayList<>();
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			if (sol.getRateForRup(rupIndex) == 0d)
				continue;
			RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 1d);
			if (surf.getQuickDistance(loc) < cutoff)
				ret.add(rupIndex);
		}
		return ret;
	}
	
	private static void fillDistanceHist(List<ProbEqkSource> sources, BiFunction<RuptureSurface, Location, Double> distFunc,
			EvenlyDiscretizedFunc hist, Location loc) {
		float min = (float)(hist.getMinX() < 0 ? hist.getMinX() - 0.5*hist.getDelta() : 0d);
		float max = (float)(hist.getMaxX() + 0.5*hist.getDelta());
		sources.parallelStream().forEach((source)->{
			for (ProbEqkRupture rup : source) {
				double dist = distFunc.apply(rup.getRuptureSurface(), loc);
				if ((float)dist >= min && (float)dist <= max) {
					int index = hist.getClosestXIndex(dist);
					synchronized (hist) {
						hist.add(index, rup.getMeanAnnualRate(1d));
					}
				}
			}
		});
	}
	
	private static Color trans(Color c, int a) {
		return new Color(c.getRed(), c.getGreen(), c.getBlue(), a);
	}
	
	private static double curveVal(DiscretizedFunc curve, ReturnPeriods rp) {
		if (curve == null)
			return Double.NaN;
		if (rp.oneYearProb > curve.getMaxY())
			return 0d;
		if (rp.oneYearProb < curve.getMinY())
			// saturated
			return curve.getMaxX();
		return curve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
//		return curve.getFirstInterpolatedX_inLogYDomain(rp.oneYearProb);
	}

}
