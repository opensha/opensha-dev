package scratch.kevin.prvi25;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Supplier;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.FixedDistanceCutoffFilter;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.hazard.QuickGriddedHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysHazardCalcSettings;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.utils.PointSourceDistanceCorrections;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenRelSupplier;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.imr.param.OtherParams.TectonicRegionTypeParam;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class HazardCurveDifferenceDebug {

	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException {
//		Site site = new Site(new Location(18, -68));
//		Site site = new Site(new Location(18.35, -67.25));
		Site site = new Site(new Location(18.325, -67.275));
		Location siteLoc = site.getLocation();
		GriddedRegion gridReg = new GriddedRegion(
				new Region(siteLoc, 1d),
				1d, siteLoc);
		Location grid0 = gridReg.getLocation(0);
		Preconditions.checkState(gridReg.getNodeCount() == 1, "Expected 1, have %s", gridReg.getNodeCount());
		Preconditions.checkState(LocationUtils.areSimilar(siteLoc, grid0), "Location mismatch: %s != %s", siteLoc, grid0);
		
		TectonicRegionType calcTRT = null;
//		TectonicRegionType calcTRT = TectonicRegionType.ACTIVE_SHALLOW;
//		TectonicRegionType calcTRT = TectonicRegionType.SUBDUCTION_INTERFACE;
//		TectonicRegionType calcTRT = TectonicRegionType.SUBDUCTION_SLAB;
		
		GriddedSeismicitySettings fullSettings = BaseFaultSystemSolutionERF.GRID_SETTINGS_DEFAULT;
//		fullSettings = fullSettings.forDistanceCorrections(PointSourceDistanceCorrections.NONE).forPointSourceMagCutoff(5d);
		GriddedSeismicitySettings quickSettings = fullSettings;
		
		fullSettings = fullSettings.forSupersamplingSettings(GridCellSupersamplingSettings.DEFAULT);
		quickSettings = quickSettings.forSupersamplingSettings(GridCellSupersamplingSettings.QUICK);
		
//		fullSettings = fullSettings.forSupersamplingSettings(null);
//		quickSettings = quickSettings.forSupersamplingSettings(null);
		
		boolean mfdsFull = true;
		boolean mfdsQuick = true;
		Double truncFull = null;
		Double truncQuick = null;
		boolean quickDoTree = false;
		boolean quickUseSLT = true;
		
		Map<TectonicRegionType, ScalarIMR> gmms = new HashMap<>();
		gmms.put(TectonicRegionType.ACTIVE_SHALLOW, AttenRelRef.USGS_PRVI_ACTIVE.get());
		gmms.put(TectonicRegionType.SUBDUCTION_INTERFACE, AttenRelRef.USGS_PRVI_INTERFACE.get());
		gmms.put(TectonicRegionType.SUBDUCTION_SLAB, AttenRelRef.USGS_PRVI_SLAB.get());
//		gmms.put(TectonicRegionType.ACTIVE_SHALLOW, AttenRelRef.ASK_2014.get());
//		gmms.put(TectonicRegionType.SUBDUCTION_INTERFACE, AttenRelRef.ASK_2014.get());
//		gmms.put(TectonicRegionType.SUBDUCTION_SLAB, AttenRelRef.ASK_2014.get());
		if (calcTRT != null) {
			List<TectonicRegionType> origTRTs = List.copyOf(gmms.keySet());
			for (TectonicRegionType trt : origTRTs)
				if (trt != calcTRT)
					gmms.remove(trt);
		}
		
		GridSourceList overrideGridList = null;
//		TectonicRegionType tempTRT = calcTRT == null ? TectonicRegionType.ACTIVE_SHALLOW : calcTRT;
//		overrideGridList = new GridSourceList.Precomputed(gridReg, calcTRT, List.of(List.of(
//				new GriddedRupture(0, grid0,
//						new GriddedRuptureProperties(8d, 90d, 50d, Double.NaN, null, 5d, 15d, 100d, 10d, Double.NaN, tempTRT)
//						, 0.1d))));
		
		Map<TectonicRegionType, Supplier<ScalarIMR>> gmmSuppliers = new HashMap<>();
		for (TectonicRegionType trt : gmms.keySet()) {
			ScalarIMR gmm = gmms.get(trt);
			gmmSuppliers.put(trt, new Supplier<ScalarIMR>() {
				@Override
				public ScalarIMR get() {
					return gmm;
				}
			});
			System.out.println(trt+": "+gmm.getName());
			TectonicRegionTypeParam trtParam = (TectonicRegionTypeParam)gmm.getParameter(TectonicRegionTypeParam.NAME);
			Preconditions.checkState(trtParam != null, "Multiple GMPEs supplied, but GMPE "+gmm.getShortName()+" doesn't have a TRT");
			System.out.println("\tGMM reports TRT: "+trtParam.getValueAsTRT());
			for (Parameter<?> param : gmm.getSiteParams())
				if (!site.containsParameter(param.getName()))
					site.addParameter((Parameter<?>)param.clone());
			
			gmm.setIntensityMeasure(PGA_Param.NAME);
		}
		
		List<Map<TectonicRegionType, ? extends Supplier<ScalarIMR>>> quickGMMSuppliers;
		List<Double> quickGMMSupplierWeights;
		if (quickDoTree) {
			List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels;
			if (calcTRT == TectonicRegionType.ACTIVE_SHALLOW)
				gmmLevels = PRVI25_LogicTree.levelsCrustalGMM;
			else if (calcTRT == TectonicRegionType.SUBDUCTION_INTERFACE)
				gmmLevels = PRVI25_LogicTree.levelsInterfaceGMM;
			else if (calcTRT == TectonicRegionType.SUBDUCTION_SLAB)
				gmmLevels = PRVI25_LogicTree.levelsSlabGMM;
			else
				gmmLevels = PRVI25_LogicTree.levelsCombinedGMM;
			LogicTree<?> gmmTree = LogicTree.buildExhaustive(gmmLevels, true);
			quickGMMSuppliers = new ArrayList<>();
			quickGMMSupplierWeights = new ArrayList<>();
			for (LogicTreeBranch<?> branch : gmmTree) {
				Map<TectonicRegionType, AttenRelSupplier> suppliers = FaultSysHazardCalcSettings.getGMM_Suppliers(branch, null, false);
				quickGMMSuppliers.add(suppliers);
				quickGMMSupplierWeights.add(gmmTree.getBranchWeight(branch));
			}
		} else {
			quickGMMSuppliers = List.of(gmmSuppliers);
			quickGMMSupplierWeights = List.of(1d);
		}
		
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		
		File invsDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/");
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File(invsDir,"2025_08_01-prvi25_crustal_subduction_combined_branches/"
						+ "combined_branch_averaged_solution.zip"));
		
		if (overrideGridList != null)
			sol.setGridSourceProvider(overrideGridList);
		erf.setSolution(sol);
		erf.getTimeSpan().setDuration(1d);
		
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		
		erf.setGriddedSeismicitySettings(fullSettings);
		erf.setParameter(UseRupMFDsParam.NAME, mfdsFull);
		
		System.out.println("Full ERF params:");
		for (Parameter<?> param : erf.getAdjustableParameterList())
			System.out.println(param.getName()+":\t"+param.getValue());
		System.out.println();
		
		erf.updateForecast();
		
		AbstractERF calcERF;
		if (calcTRT == null)
			calcERF = erf;
		else
			calcERF = new GMMTreeCalcDebug.TRTWrappedERF(erf, calcTRT, null);
		
		Map<Class<?>, Integer> sourceClassCounts = new HashMap<>();
		for (ProbEqkSource source : calcERF) {
			Class<? extends ProbEqkSource> clazz = source.getClass();
			int prevCount = sourceClassCounts.containsKey(clazz) ? sourceClassCounts.get(clazz) : 0;
			sourceClassCounts.put(clazz, prevCount+1);
		}
		System.out.println("\nCalc ERF has "+calcERF.getNumSources()+" sources and "+sourceClassCounts.size()+" unique types:");
		for (Class<?> clazz : sourceClassCounts.keySet())
			System.out.println("\t"+ClassUtils.getClassNameWithoutPackage(clazz)+":\t"+sourceClassCounts.get(clazz));
		System.out.println();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		ArbitrarilyDiscretizedFunc newXVals = new ArbitrarilyDiscretizedFunc();
		newXVals.set(1e-9, 1d);
		newXVals.set(1e-8, 1d);
		newXVals.set(1e-7, 1d);
		newXVals.set(1e-6, 1d);
		newXVals.set(1e-5, 1d);
		for (Point2D pt : xVals)
			newXVals.set(pt);
		xVals = newXVals;
		
		for (ScalarIMR gmm : gmms.values()) {
			if (truncFull == null) {
				gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_NONE);
			} else {
				gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
				gmm.getParameter(SigmaTruncLevelParam.NAME).setValue(truncFull);
			}
		}
		
		DiscretizedFunc fullCurve = new CurveCalc(site, calcERF, gmms, xVals).get();
		
		if (quickUseSLT) {
			SolutionLogicTree slt = SolutionLogicTree.load(new File(invsDir,
					"2025_08_01-prvi25_crustal_subduction_combined_branches-ba_only-gmTreeCalcs-vs760/fake_erf_slt.zip"));
			sol = slt.forBranch(slt.getLogicTree().getBranch(0));
			
			sol.write(new File("/tmp/fss_from_slt.zip"));
			
			if (overrideGridList != null)
				sol.setGridSourceProvider(overrideGridList);
			erf.setSolution(sol);
		}
		
		erf.setGriddedSeismicitySettings(quickSettings);
		erf.setParameter(UseRupMFDsParam.NAME, mfdsQuick);
		
		System.out.println("Quick ERF params:");
		for (Parameter<?> param : erf.getAdjustableParameterList())
			System.out.println(param.getName()+":\t"+param.getValue());
		System.out.println();
		
		erf.updateForecast();
		
		if (calcTRT == null)
			calcERF = erf;
		else
			calcERF = new GMMTreeCalcDebug.TRTWrappedERF(erf, calcTRT, null);
		
		for (ScalarIMR gmm : gmms.values()) {
			if (truncQuick == null) {
				gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_NONE);
			} else {
				gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
				gmm.getParameter(SigmaTruncLevelParam.NAME).setValue(truncQuick);
			}
		}
		
		ExecutorService exec = Executors.newSingleThreadExecutor();
		DiscretizedFunc avgQuickCurve = xVals.deepClone();
		double sumQuickWeight = 0d;
		avgQuickCurve.scale(0d);
		for (int b=0; b<quickGMMSuppliers.size(); b++) {
			Map<TectonicRegionType, ? extends Supplier<ScalarIMR>> suppliers = quickGMMSuppliers.get(b);
			double weight = quickGMMSupplierWeights.get(b);
			sumQuickWeight += weight;
			
			QuickGriddedHazardMapCalc quickCalc = new QuickGriddedHazardMapCalc(suppliers, 0, xVals, sourceFilters, quickSettings);
			
			GridSourceProvider gridProv;
			if (calcTRT == null) {
				gridProv = erf.getGridSourceProvider();
			} else {
				GridSourceList fullGridProv = (GridSourceList) erf.getGridSourceProvider();
				List<List<GriddedRupture>> ruptureLists = new ArrayList<>();
				for (int l=0; l<fullGridProv.getNumLocations(); l++) {
					List<GriddedRupture> ruptures = new ArrayList<>();
					ruptures.addAll(fullGridProv.getRuptures(calcTRT, l));
					ruptureLists.add(ruptures);
				}
				Preconditions.checkState(ruptureLists.size() == fullGridProv.getNumLocations());
				gridProv = new GridSourceList.Precomputed(fullGridProv.getGriddedRegion(), calcTRT, ruptureLists);
			}
			DiscretizedFunc quickCurve = quickCalc.calc(gridProv, gridReg, exec, 1)[0];
			if (calcTRT != TectonicRegionType.SUBDUCTION_SLAB) {
				// need to add on-fault
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
				erf.updateForecast();
				
				if (calcTRT == null)
					calcERF = erf;
				else
					calcERF = new GMMTreeCalcDebug.TRTWrappedERF(erf, calcTRT, null);
				
				Map<TectonicRegionType, ScalarIMR> quickGMMs = new HashMap<>();
				for (TectonicRegionType trt : suppliers.keySet())
					quickGMMs.put(trt, suppliers.get(trt).get());
				
				DiscretizedFunc excludeCurve = new CurveCalc(site, calcERF, quickGMMs, xVals).get();
				ArbitrarilyDiscretizedFunc combCurve = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<quickCurve.size(); i++)
					combCurve.set(quickCurve.getX(i), 1d - (1d - quickCurve.getY(i))*(1d - excludeCurve.getY(i)));
				quickCurve = combCurve;
			}
			for (int i=0; i<quickCurve.size(); i++)
				avgQuickCurve.set(i, avgQuickCurve.getY(i) + weight*quickCurve.getY(i));
		}
		
		exec.shutdown();
		
		if ((float)sumQuickWeight != 1f)
			avgQuickCurve.scale(1d/sumQuickWeight);
		
		for (int i=0; i<fullCurve.size(); i++) {
			double x = fullCurve.getX(i);
			double y1 = fullCurve.getY(i);
			double y2 = avgQuickCurve.getY(i);
			double diff = y2 - y1;
			double pDiff = 100d * diff / y1;
			System.out.println("X: "+(float)x+"\tFull: "+(float)y1+"\tQuick: "+(float)y2+"\tDiff: "+(float)diff+" ("+(float)pDiff+" %)");
		}
		
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		double rpProb = rp.oneYearProb;
		System.out.println("Return period: "+rp+" at y="+rpProb);
		double rp1 = fullCurve.getFirstInterpolatedX_inLogXLogYDomain(rpProb);
		double rp2 = avgQuickCurve.getFirstInterpolatedX_inLogXLogYDomain(rpProb);
		double diff = rp2 - rp1;
		double pDiff = 100d * diff / rp1;
		System.out.println("Full 2in50: "+(float)rp1+"\tQuick 2in50: "+(float)rp2+"\tDiff: "+(float)diff+" ("+(float)pDiff+" %)");
	}
	
//	private static SourceFilterManager sourceFilters = new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS);
	private static SourceFilterManager sourceFilters = new SourceFilterManager(SourceFilters.FIXED_DIST_CUTOFF);
	static {
		sourceFilters.getFilterInstance(FixedDistanceCutoffFilter.class).setMaxDistance(10000d);
	}
	
	private static class CurveCalc implements Supplier<DiscretizedFunc> {
		
		private Site site;
		private AbstractERF erf;
		private Map<TectonicRegionType, ScalarIMR> gmms;
		private DiscretizedFunc xVals;

		public CurveCalc(Site site, AbstractERF erf, Map<TectonicRegionType, ScalarIMR> gmms, DiscretizedFunc xVals) {
			this.site = site;
			this.erf = erf;
			this.gmms = gmms;
			this.xVals = xVals;
		}

		@Override
		public DiscretizedFunc get() {
			double[] logXValsArray = new double[xVals.size()];
			for (int i=0; i<xVals.size(); i++)
				logXValsArray[i] = Math.log(xVals.getX(i));
			LightFixedXFunc logXVals = new LightFixedXFunc(logXValsArray, new double[xVals.size()]);
			HazardCurveCalculator calc = new HazardCurveCalculator(sourceFilters);
			
			for (ScalarIMR gmm : gmms.values())
				gmm.setIntensityMeasure(PGA_Param.NAME);
			
			calc.getHazardCurve(logXVals, site, gmms, erf);
			
			DiscretizedFunc ret = xVals.deepClone();
			for (int i=0; i<xVals.size(); i++)
				ret.set(i, logXVals.getY(i));
			
			return ret;
		}
		
	}

}
