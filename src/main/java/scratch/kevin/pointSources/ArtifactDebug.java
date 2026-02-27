package scratch.kevin.pointSources;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.function.Supplier;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Stopwatch;
import com.google.common.collect.ImmutableList;

public class ArtifactDebug {

	public static void main(String[] args) throws IOException {
		Region subReg = new Region(new Location(36.5, -119.75), new Location(38, -118.75));
		GriddedRegion gridReg = new GriddedRegion(subReg, 0.02, GriddedRegion.ANCHOR_0_0);
		
		File cacheDir = new File("/tmp");
		String cachePrefix = "artifact_debug_cache";
		
		FaultSystemSolution fullSol = FaultSystemSolution.load(ORIG_SOL_FILE);
		
		double period = 0d;
		int threads = FaultSysTools.defaultNumThreads();
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;

		CPT cpt = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-5d, 5d);
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-0.8, -0.3);
		hazCPT.setLog10(true);

//		Models testModel = Models.SPINNING_AVG_CENTERED_M5;
//		Models refModel = Models.SPINNING_AVG_CENTERED_M6;
//		Models testModel = Models.SPINNING_DIST_5X_UNCENTERED;
//		Models refModel = Models.SPINNING_DIST_5X_CENTERED;
		Models testModel = Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5;
		Models refModel = Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN;

		GridSourceList origSources = fullSol.requireModule(GridSourceList.class);
		GridSourceList gridSources = origSources;
		
//		Models testModel = Models.AS_PUBLISHED;
//		Models refModel = Models.TRUE_POINT;
//		cpt = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-12d, 12d);
		
//		Supplier<ScalarIMR> gmmSupplier = AttenRelRef.WRAPPED_ASK_2014;
		Supplier<ScalarIMR> gmmSupplier = AttenRelRef.USGS_NSHM23_ACTIVE;
		
//		String suffix = "";
//		Consumer<SolHazardMapCalc> parameterizer = null;
		
		String suffix = "_len_shift";
		Consumer<SolHazardMapCalc> parameterizer = C -> {
			fullSol.setGridSourceProvider(new LengthShifter(fullSol.requireModule(GridSourceList.class), 1.5));
		};
		
//		String suffix = "_no_ss";
//		Consumer<SolHazardMapCalc> parameterizer = C -> {
//			C.setSupersamplingSettings(null);
//		};
		
//		String suffix = "_far_ss";
//		Consumer<SolHazardMapCalc> parameterizer =C -> {
//			C.setSupersamplingSettings(new GridCellSupersamplingSettings(1d, 60, 100, 0, true));
//		};
		
//		String suffix = "_no_snap_rup_props";
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setSnapRupProps(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setSnapRupProps(false);
//		Consumer<SolHazardMapCalc> parameterizer = null;
		
//		String suffix = "_no_cache_no_snap_rup_props";
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setSnapRupProps(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setSnapRupProps(false);
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setCacheForRups(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setCacheForRups(false);
//		Consumer<SolHazardMapCalc> parameterizer = null;
		
		// implies no snap
//		String suffix = "_no_cache_dists";
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setCacheForRups(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setCacheForRups(false);
//		Consumer<SolHazardMapCalc> parameterizer = null;
		
//		String suffix = "_no_grid_optimize";
//		Consumer<SolHazardMapCalc> parameterizer = C -> {
//			C.setPointSourceOptimizations(false);
//		};
		
//		String suffix = "_no_cache_dists_no_grid_optimize";
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setCacheForRups(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setCacheForRups(false);
//		Consumer<SolHazardMapCalc> parameterizer = C -> {
//			C.setPointSourceOptimizations(false);
//		};
		
//		String suffix = "_no_cache_no_snap_rup_props_no_grid_optimize";
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setSnapRupProps(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setSnapRupProps(false);
//		((DistanceDistributionCorrection)testModel.distCorr.get()).setCacheForRups(false);
//		((DistanceDistributionCorrection)refModel.distCorr.get()).setCacheForRups(false);
//		Consumer<SolHazardMapCalc> parameterizer = C -> {
//			C.setPointSourceOptimizations(false);
//		};
		
		File excludeCurveFile = new File(cacheDir, SolHazardMapCalc.getCSV_FileName(cachePrefix, period));
		SolHazardMapCalc excludeCalc;
		if (excludeCurveFile.exists()) {
			excludeCalc = SolHazardMapCalc.loadCurves(fullSol, gridReg, new double[] {period}, cacheDir, cachePrefix);
		} else {
			excludeCalc = new SolHazardMapCalc(fullSol, gmmSupplier, gridReg, period);
			
			System.out.println("Calculating exluding");
			excludeCalc.setBackSeisOption(IncludeBackgroundOption.EXCLUDE);
			excludeCalc.calcHazardCurves(threads);
			
			excludeCalc.writeCurvesCSV(excludeCurveFile, period);
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(subReg);
		mapMaker.setFaultSections(fullSol.getRupSet().getFaultSectionDataList());
		mapMaker.setSectOutlineChar(null);
		mapMaker.plotScatters(gridSources.getGriddedRegion().getNodeList(), Color.BLACK);
		mapMaker.setScatterSymbol(PlotSymbol.CIRCLE, 3f);
		
		GriddedGeoDataSet excludeMap = excludeCalc.buildMap(period, rp);
		mapMaker.plotXYZData(excludeMap, hazCPT, "PGA 2 in 50");
		mapMaker.plot(cacheDir, "artifact_debug_exclude", " ");
		
		System.out.println("Calculating "+testModel);
		if (testModel.getGridModFunction() != null)
			fullSol.setGridSourceProvider(testModel.getGridModFunction().apply(gridSources));
		SolHazardMapCalc calc = new SolHazardMapCalc(fullSol, gmmSupplier, gridReg, period);
		calc.setBackSeisOption(IncludeBackgroundOption.ONLY);
		calc.setGriddedSeismicitySettings(testModel.getGridProps());
		if (parameterizer != null)
			parameterizer.accept(calc);
		Stopwatch watch = Stopwatch.createStarted();
		calc.calcHazardCurves(threads, excludeCalc);
		watch.stop();
		System.out.println("Took "+(float)(watch.elapsed(TimeUnit.MILLISECONDS)/1000d)+" s");
		fullSol.setGridSourceProvider(gridSources);
		
		GriddedGeoDataSet testMap = calc.buildMap(period, rp);
		
		System.out.println("Calculating "+refModel);
		if (refModel.getGridModFunction() != null)
			fullSol.setGridSourceProvider(refModel.getGridModFunction().apply(gridSources));
		calc = new SolHazardMapCalc(fullSol, gmmSupplier, gridReg, period);
		calc.setBackSeisOption(IncludeBackgroundOption.ONLY);
		calc.setGriddedSeismicitySettings(refModel.getGridProps());
		if (parameterizer != null)
			parameterizer.accept(calc);
		watch = Stopwatch.createStarted();
		calc.calcHazardCurves(threads, excludeCalc);
		watch.stop();
		System.out.println("Took "+(float)(watch.elapsed(TimeUnit.MILLISECONDS)/1000d)+" s");
		
		GriddedGeoDataSet refMap = calc.buildMap(period, rp);
		
		GriddedGeoDataSet pDiff = new GriddedGeoDataSet(gridReg);
		for (int i=0; i<pDiff.size(); i++)
			pDiff.set(i, 100d*(testMap.get(i) - refMap.get(i))/refMap.get(i));
		
		System.out.println("Plotting map");
		mapMaker.plotXYZData(pDiff, cpt, "% change");
		mapMaker.plot(cacheDir, "artifact_debug_"+testModel.name()+"_"+refModel.name()+suffix, " ");
		mapMaker.plotXYZData(testMap, hazCPT, "PGA 2 in 50");
		mapMaker.plot(cacheDir, "artifact_debug_raw_"+testModel.name()+suffix, " ");
		mapMaker.plotXYZData(refMap, hazCPT, "PGA 2 in 50");
		mapMaker.plot(cacheDir, "artifact_debug_raw_"+refModel.name()+suffix, " ");
		
		System.out.println("DONE");
	}
	
	private static class LengthShifter extends GridSourceList.DynamicallyBuilt {
		
		private GridSourceList origSources;
		private double lenScalar;

		public LengthShifter(GridSourceList origSources, double lenScalar) {
			super(origSources.getTectonicRegionTypes(), origSources.getGriddedRegion(), origSources.getRefMFD());
			this.origSources = origSources;
			this.lenScalar = lenScalar;
			
		}

		@Override
		public int getNumSources() {
			return origSources.getNumSources();
		}
		
		@Override
		public int getLocationIndexForSource(int sourceIndex) {
			return origSources.getLocationIndexForSource(sourceIndex);
		}
		
		@Override
		public TectonicRegionType tectonicRegionTypeForSourceIndex(int sourceIndex) {
			return origSources.tectonicRegionTypeForSourceIndex(sourceIndex);
		}
		
		@Override
		public Set<Integer> getAssociatedGridIndexes(int sectionIndex) {
			return origSources.getAssociatedGridIndexes(sectionIndex);
		}
		
		@Override
		protected List<GriddedRupture> buildRuptures(TectonicRegionType tectonicRegionType, int gridIndex) {
			ImmutableList<GriddedRupture> rups = origSources.getRuptures(tectonicRegionType, gridIndex);
			List<GriddedRupture> modRups = new ArrayList<>(rups.size());
			for (GriddedRupture rup : rups) {
				GriddedRuptureProperties modProps = new GriddedRuptureProperties(rup.properties.magnitude,
						rup.properties.rake, rup.properties.dip, rup.properties.strike, rup.properties.strikeRange,
						rup.properties.upperDepth, rup.properties.lowerDepth, rup.properties.length*lenScalar,
						rup.properties.hypocentralDepth, rup.properties.hypocentralDAS, rup.properties.tectonicRegionType);
				modRups.add(new GriddedRupture(rup.gridIndex, rup.location, modProps, rup.rate));
			}
			return modRups;
		}
		
	}

}
