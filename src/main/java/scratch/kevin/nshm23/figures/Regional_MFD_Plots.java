package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.LocalRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.NSHM23_BaseRegion;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.StitchedRegions;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;

class Regional_MFD_Plots {
	
	private static final File NSHM18 = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.3.0");
	
	private static final File NSHM23_WRAPPED = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.b.1");
	
	private static final File NSHM23_DIR = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
			+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
	private static final File NSHM23_PLOTS_DIR = new File(NSHM23_DIR, "misc_plots");
	private static final File NSHM23_SOL = new File(NSHM23_DIR, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
	
	private static final File U3_SOL = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/branch_avgs_combined.zip");
	
	// TODO update
	private static final File METHODS_SOL = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/branch_avgs_combined.zip");

	public static void main(String[] args) throws IOException {
		doCompCascadia();
//		doCompU3();
		doCompNSHM18();
		doCompEast();
//		doMethodsCompU3();
	}
	
	private static void doCompEast() throws IOException {
		File outputDir = new File(NSHM23_PLOTS_DIR, "reg_mfds_comp_east");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_BaseRegion[] analysis = {
				AnalysisRegions.CONUS_EAST
		};
		
		HazardModel nshm23 = HazardModel.load(NSHM23_WRAPPED.toPath());
		HazardModel nshm18 = HazardModel.load(NSHM18.toPath());
		boolean subduction = false;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,
				TectonicRegionType.STABLE_SHALLOW);
		if (subduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);

		System.out.println("Calculating NSHM23 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> modelMFDs = calcModelMFDs(nshm23, trts, analysis, refMFD);
		System.out.println("Calculating NSHM18 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> compMFDs = calcModelMFDs(nshm18, trts, analysis, refMFD);
		
		String modelName = "NSHM23";
		String compName = "NSHM18";
		
		Range xRange = new Range(5d, 8.5d);
		Range yRange = new Range(1e-7, 1e0);
		
		boolean plotAllTypes = true;
		boolean plotTarget = false;
//		DIST_TRANS_COLOR = TRANS_LIGHT_GRAY;
		
		for (NSHM23_BaseRegion reg : analysis) {
			writePlot(reg, refMFD, null, false, plotAllTypes, null, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, compMFDs.row(reg));
		}
	}
	
	private static void doCompCascadia() throws IOException {
		File outputDir = new File(NSHM23_PLOTS_DIR, "reg_mfds_comp_cascadia");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_BaseRegion[] analysis = {
				NSHM23_RegionLoader.CATCH_ALL_REGION
		};
		
		HazardModel nshm23 = HazardModel.load(NSHM23_WRAPPED.toPath());
		HazardModel nshm18 = HazardModel.load(NSHM18.toPath());
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(9.45);

		System.out.println("Calculating NSHM23 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> modelMFDs = calcCascadiaModelMFDs(nshm23, analysis, refMFD);
		System.out.println("Calculating NSHM18 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> compMFDs = calcCascadiaModelMFDs(nshm18, analysis, refMFD);
		
		String modelName = "NSHM23";
		String compName = "NSHM18";
		
		Range xRange = new Range(5d, 9.5d);
		Range yRange = new Range(1e-7, 1e0);
		
		boolean plotAllTypes = true;
		boolean plotTarget = false;
//		DIST_TRANS_COLOR = TRANS_LIGHT_GRAY;
		
		for (NSHM23_BaseRegion reg : analysis) {
			writePlot(reg, refMFD, null, true, plotAllTypes, null, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, compMFDs.row(reg));
		}
	}
	
	private static void doCompNSHM18() throws IOException {
		File outputDir = new File(NSHM23_PLOTS_DIR, "reg_mfds_comp_nshm18");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_BaseRegion[] analysis = {
				AnalysisRegions.CONUS_U3_RELM,
				AnalysisRegions.CONUS_IMW,
				AnalysisRegions.CONUS_PNW,
				LocalRegions.CONUS_SF_BAY,
				LocalRegions.CONUS_LA_BASIN,
				LocalRegions.CONUS_NEW_MADRID,
				LocalRegions.CONUS_PUGET,
				LocalRegions.CONUS_WASATCH,
				StitchedRegions.CONUS_WEST
		};
		
		FaultSystemSolution modelSol = FaultSystemSolution.load(NSHM23_SOL);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);
		
		System.out.println("Calculating Model MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> modelMFDs = calcSolMFDs(modelSol, analysis, refMFD);
		
		HazardModel nshm18 = HazardModel.load(NSHM18.toPath());
		boolean subduction = false;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,
				TectonicRegionType.STABLE_SHALLOW);
		if (subduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}
		
		System.out.println("Calculating NSHM18 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> compMFDs = calcModelMFDs(nshm18, trts, analysis, refMFD);
		
		String modelName = "NSHM23";
		String compName = "NSHM18";
		
		Range xRange = new Range(5d, 8.5d);
		Range yRange = new Range(1e-6, 1e1);
		
		boolean plotAllTypes = true;
		boolean plotTarget = false;
//		DIST_TRANS_COLOR = TRANS_LIGHT_GRAY;
		
		for (NSHM23_BaseRegion reg : analysis) {
			writePlot(reg, refMFD, MFDType.SUPRA_ONLY, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, compMFDs.row(reg));
			writePlot(reg, refMFD, MFDType.GRID_ONLY, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, compMFDs.row(reg));
			writePlot(reg, refMFD, MFDType.SUM, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, compMFDs.row(reg));
			writePlot(reg, refMFD, null, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, compMFDs.row(reg));
		}
	}
	
	private static void doCompU3() throws IOException {
		File outputDir = new File(NSHM23_PLOTS_DIR, "reg_mfds_comp_u3");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_BaseRegion[] analysis = {
				AnalysisRegions.CONUS_U3_RELM,
				LocalRegions.CONUS_SF_BAY,
				LocalRegions.CONUS_LA_BASIN
		};
		
		FaultSystemSolution u3Sol = FaultSystemSolution.load(U3_SOL);
		
		FaultSystemSolution modelSol = FaultSystemSolution.load(NSHM23_SOL);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);
		
		System.out.println("Calculating U3 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> u3MFDs = calcSolMFDs(u3Sol, analysis, refMFD);
		System.out.println("Calculating Model MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> modelMFDs = calcSolMFDs(modelSol, analysis, refMFD);
		
		String modelName = "NSHM23";
		String compName = "UCERF3";
		
		Range xRange = new Range(5d, 8.5d);
		Range yRange = new Range(1e-6, 1e1);
		
		boolean plotAllTypes = true;
		boolean plotTarget = false;
//		DIST_TRANS_COLOR = TRANS_LIGHT_GRAY;
		
		for (NSHM23_BaseRegion reg : analysis) {
			writePlot(reg, refMFD, MFDType.SUPRA_ONLY, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, u3MFDs.row(reg));
			writePlot(reg, refMFD, MFDType.GRID_ONLY, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, u3MFDs.row(reg));
			writePlot(reg, refMFD, MFDType.SUM, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, u3MFDs.row(reg));
			writePlot(reg, refMFD, null, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					null, compName, u3MFDs.row(reg));
		}
	}
	
	private static void doMethodsCompU3() throws IOException {
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/reg_mfds_comp_u3");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_BaseRegion[] analysis = {
				AnalysisRegions.CONUS_U3_RELM,
				LocalRegions.CONUS_SF_BAY,
				LocalRegions.CONUS_LA_BASIN
		};
		
		FaultSystemSolution u3Sol = FaultSystemSolution.load(U3_SOL);
		
		FaultSystemSolution methodsSol = FaultSystemSolution.load(METHODS_SOL);
		
		FaultSystemSolution modelSol = FaultSystemSolution.load(NSHM23_SOL);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);
		
		System.out.println("Calculating U3 MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> u3MFDs = calcSolMFDs(u3Sol, analysis, refMFD);
		System.out.println("Calculating Methods MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> methodsMFDs = calcSolMFDs(methodsSol, analysis, refMFD);
		System.out.println("Calculating Model MFDs");
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> modelMFDs = calcSolMFDs(modelSol, analysis, refMFD);
		
		String modelName = "NSHM23";
		String compName = "UCERF3";
		
		Range xRange = new Range(6d, 8.5d);
		Range yRange = new Range(1e-6, 2e0);
		
		boolean plotAllTypes = false;
		boolean plotTarget = true;
//		DIST_TRANS_COLOR = TRANS_LIGHT_RED;
		
		for (NSHM23_BaseRegion reg : analysis) {
			writePlot(reg, refMFD, MFDType.SUPRA_ONLY, false, plotAllTypes, modelSol, plotTarget, outputDir, xRange, yRange, " ", modelName, modelMFDs.row(reg),
					methodsMFDs.row(reg), compName, u3MFDs.row(reg));
		}
	}
	
	private static enum ModelType {
		COMP,
		METHODS,
		MODEL
	}
	
	private static void writePlot(NSHM23_BaseRegion aReg, EvenlyDiscretizedFunc refMFD, MFDType type, boolean subduction,
			boolean plotAllTypes, FaultSystemSolution sol, boolean plotTarget, File outputDir, Range xRange, Range yRange, String title,
			String modelName, Map<MFDType, IncrementalMagFreqDist> modelMFDs,
			Map<MFDType, IncrementalMagFreqDist> methodsMFDs,
			String compName, Map<MFDType, IncrementalMagFreqDist> compMFDs)
					throws IOException {
		System.out.println("Plotting for "+aReg+", "+type);
		Region region = aReg.load();
		
		List<XY_DataSet> incrFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
		List<XY_DataSet> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> cmlChars = new ArrayList<>();
		
		List<String> csvLabels = new ArrayList<>();
		List<IncrementalMagFreqDist> csvFuncs = new ArrayList<>();
		
		if (aReg != NSHM23_RegionLoader.CATCH_ALL_REGION) {
			UncertainBoundedIncrMagFreqDist dataBounds = NSHM23_RegionalSeismicity.getRemapped(region,
					NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
			dataBounds.setName("Observed");
			dataBounds.setBoundName("95% Bounds");
			UncertainBoundedIncrMagFreqDist dataForBounds = dataBounds.deepClone();
			dataForBounds.setName(dataBounds.getBoundName());
			
			// add observed bounds
			Color obsColor = new Color(125, 80, 145); // "indigo"
			incrFuncs.add(dataBounds);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, obsColor));
			EvenlyDiscretizedFunc dataCumulative = getCmlAsFakeIncr(dataBounds);
			dataCumulative.setName(dataBounds.getName());
			cmlFuncs.add(dataCumulative);
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, obsColor));
			dataCumulative = dataCumulative.deepClone();
			dataCumulative.setName(dataBounds.getName());
			
			csvLabels.add("Observed");
			csvFuncs.add(dataBounds);
			csvLabels.add("Observed, p2.5");
			csvFuncs.add(dataBounds.getLower());
			csvLabels.add("Observed, p97.5");
			csvFuncs.add(dataBounds.getUpper());
			
			incrFuncs.add(dataForBounds);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60)));
			
			EvenlyDiscretizedFunc upperCumulative = getCmlAsFakeIncr(dataBounds.getUpper());
			EvenlyDiscretizedFunc lowerCumulative = getCmlAsFakeIncr(dataBounds.getLower());
			Preconditions.checkState(dataCumulative.size() == upperCumulative.size());
			for (int i=0; i<dataCumulative.size(); i++) {
				upperCumulative.set(i, Math.max(dataCumulative.getY(i), upperCumulative.getY(i)));
				lowerCumulative.set(i, Math.max(0, Math.min(dataCumulative.getY(i), lowerCumulative.getY(i))));
			}
			
			UncertainArbDiscFunc cmlBounded = new UncertainArbDiscFunc(dataCumulative, lowerCumulative, upperCumulative);
			cmlBounded.setName(dataBounds.getBoundName());
			cmlFuncs.add(cmlBounded);
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60)));
		}
		
		Preconditions.checkState(incrFuncs.size() == cmlFuncs.size());
		Preconditions.checkState(incrFuncs.size() == incrChars.size());
		Preconditions.checkState(incrFuncs.size() == cmlChars.size());
		List<Integer> onTopIndexes = new ArrayList<>();

		String plotModelName = modelName;
		if (methodsMFDs != null && !methodsMFDs.isEmpty() && (type == null || methodsMFDs.containsKey(type)))
			plotModelName += " Model";
		
		PlotLineType sumLine = PlotLineType.SOLID;
		PlotLineType supraLine = PlotLineType.DASHED;
		PlotLineType gridLine = PlotLineType.DOTTED;
		
		if (type == null)
			plotAllTypes = true;
		
//		boolean highlightsSolid = true;
//		float highlightThickness = plotAllTypes && type != null ? 5f : 4f;
//		float regularThickness = 2f;
		
		boolean highlightsSolid = false;
		float highlightThickness = 3f;
		float regularThickness = 3f;
		
		MFDType[] sortedTypes = {
				MFDType.SUM,
				MFDType.SUPRA_ONLY,
				MFDType.GRID_ONLY
		};
		
		boolean addTypeToLabel = true;
		
		if (type == null) {
			// make sure we have nonzero gridded
			boolean haveGridded = false;
			if (modelMFDs != null && modelMFDs.containsKey(MFDType.GRID_ONLY) && modelMFDs.get(MFDType.GRID_ONLY).calcSumOfY_Vals() > 0d)
				haveGridded = true;
			if (compMFDs != null && compMFDs.containsKey(MFDType.GRID_ONLY) && compMFDs.get(MFDType.GRID_ONLY).calcSumOfY_Vals() > 0d)
				haveGridded = true;
			if (methodsMFDs != null && methodsMFDs.containsKey(MFDType.GRID_ONLY) && methodsMFDs.get(MFDType.GRID_ONLY).calcSumOfY_Vals() > 0d)
				haveGridded = true;
			
			if (!haveGridded) {
				sortedTypes = new MFDType[] { MFDType.SUPRA_ONLY };
				supraLine = PlotLineType.SOLID;
				addTypeToLabel = false;
			}
		}
		
		for (ModelType modelType : ModelType.values()) {
			String name;
			Color color;
			Map<MFDType, IncrementalMagFreqDist> map;
			if (modelType == ModelType.COMP) {
				map = compMFDs;
				name = compName;
				color = Color.BLUE;
			} else if (modelType == ModelType.METHODS) {
				map = methodsMFDs;
				name = modelName+" Methods";
				color = Color.GREEN.darker();
			} else if (modelType == ModelType.MODEL) {
				map = modelMFDs;
				name = plotModelName;
				color = Color.RED;
			} else {
				throw new IllegalStateException();
			}
			if (map == null)
				continue;
			
			for (MFDType myType : sortedTypes) {
				if (map.get(myType) == null)
					continue;
				if (!plotAllTypes && type != myType)
					continue;
				boolean highlight = type == myType;
				
				float thickness = highlight ? highlightThickness : regularThickness;
				PlotLineType line;
				String label = name;
				switch (myType) {
				case SUM:
					line = sumLine;
					if (addTypeToLabel)
						label += ", Total";
					break;
				case SUPRA_ONLY:
					line = supraLine;
					if (addTypeToLabel)
						label += subduction ? ", Interface" : ", Faults";
					break;
				case GRID_ONLY:
					line = gridLine;
					if (addTypeToLabel)
						label += subduction ? ", Slab" : ", Gridded";
					break;
				default:
					throw new IllegalStateException();
				}
				if (!plotAllTypes || highlight && highlightsSolid)
					// always solid if we're just plotting one thing, or if we're highlighting and chose to always
					// make the highlight solid
					line = PlotLineType.SOLID;
				
				onTopIndexes.add(incrFuncs.size());
				addMFDs(map.get(myType), new PlotCurveCharacterstics(line, thickness, color),
						label, incrFuncs, incrChars, cmlFuncs, cmlChars, csvLabels, csvFuncs);
				
				if (modelType == ModelType.MODEL && type == myType && sol != null) {
					addSolRegionalBounds(type, sol, xRange, region, incrFuncs, incrChars,
							cmlFuncs, cmlChars, csvLabels, csvFuncs, refMFD);
					
					if (type == MFDType.SUPRA_ONLY && plotTarget) {
						// add model target
						FaultSystemRupSet modelRupSet = sol.getRupSet();
						List<? extends IncrementalMagFreqDist> targetMFDs = modelRupSet.getModule(InversionTargetMFDs.class).getOnFaultSupraSeisNucleationMFDs();
						if (targetMFDs != null) {
							IncrementalMagFreqDist target = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
							double[] modelSectFracts = modelRupSet.getFractSectsInsideRegion(region, false);
							for (int s=0; s<modelRupSet.getNumSections(); s++) {
								if (modelSectFracts[s] > 0) {
									IncrementalMagFreqDist mfd = targetMFDs.get(s);
									for (int i=0; i<mfd.size(); i++)
										target.set(i, target.getY(i)+mfd.getY(i)*modelSectFracts[s]);
								}
							}
							target.setName("Target");
							incrFuncs.add(target);
							incrChars.add(new PlotCurveCharacterstics(supraLine, 3f, Color.GRAY));
							cmlFuncs.add(target.getCumRateDistWithOffset());
							cmlChars.add(new PlotCurveCharacterstics(supraLine, 3f, Color.GRAY));
						}
					}
				}
			}
		}
		
		Preconditions.checkState(incrFuncs.size() == cmlFuncs.size());
		Preconditions.checkState(incrFuncs.size() == incrChars.size());
		Preconditions.checkState(incrFuncs.size() == cmlChars.size());
		
		// now add again on top without names
		for (int i : onTopIndexes) {
			XY_DataSet incr = incrFuncs.get(i);
			incr = incr.deepClone();
			incr.setName(null);
			incrFuncs.add(incr);
			incrChars.add(incrChars.get(i));
			XY_DataSet cml = cmlFuncs.get(i);
			cml = cml.deepClone();
			cml.setName(null);
			cmlFuncs.add(cml);
			cmlChars.add(cmlChars.get(i));
		}
		
		Preconditions.checkState(incrFuncs.size() == cmlFuncs.size());
		Preconditions.checkState(incrFuncs.size() == incrChars.size());
		Preconditions.checkState(incrFuncs.size() == cmlChars.size());
		
		PlotSpec incrSpec = new PlotSpec(incrFuncs, incrChars, title, "Magnitude", "Incremental Rate (/yr)");
		incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		PlotSpec cmlSpec = new PlotSpec(cmlFuncs, cmlChars, title, "Magnitude", "Cumulative Rate (/yr)");
		cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setPlotLabelFontSize(26);
		gp.setAxisLabelFontSize(28);
		gp.setTickLabelFontSize(24);
		
		gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
		
		String prefix = aReg.name()+"_mfds";
		if (type != null) {
			prefix += "_"+type.name();
			if (!plotAllTypes)
				prefix += "_only";
		}
		
		int width = xRange.getLowerBound() < 5.5d || xRange.getUpperBound() > 8.6 ? 950 : 800;
		
		PlotUtils.writePlots(outputDir, prefix+"_incr", gp, width, 750, true, true, false);
		
		gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_cml", gp, width, 750, true, true, false);
		
		System.out.println(aReg.name());
		for (XY_DataSet func : cmlFuncs) {
			if (func.getName() != null) {
				for (Point2D pt : func) {
					if ((float)pt.getX() == (float)xRange.getLowerBound()) {
						System.out.println(func.getName()+": rate M>"+(float)pt.getX()+"="+(float)pt.getY());
					}
				}
			}
		}
		
		// csv files
		CSVFile<String> incrCSV = new CSVFile<>(true);
		CSVFile<String> cmlCSV = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("Magnitude");
		header.addAll(csvLabels);
		
		incrCSV.addLine(header);
		cmlCSV.addLine(header);
		
		List<EvenlyDiscretizedFunc> csvCmlFuncs = new ArrayList<>();
		EvenlyDiscretizedFunc cmlRef = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta()).getCumRateDistWithOffset();
		for (int i=0; i<csvFuncs.size(); i++) {
			String label = csvLabels.get(i);
			IncrementalMagFreqDist mfd = csvFuncs.get(i);
			csvCmlFuncs.add(mfd.getCumRateDistWithOffset());
			Preconditions.checkState(mfd.size() <= refMFD.size(),
					"Bad size=%s (expected %s) for MFD %s", mfd.size(), refMFD.size(), label);
			Preconditions.checkState(mfd.getMinX() == refMFD.getMinX(),
					"Bad minX=%s (expected %s) for MFD %s", mfd.getMinMagWithNonZeroRate(), refMFD.getMinX(), label);
		}
		
		for (int i=0; i<refMFD.size(); i++) {
			List<String> incrLine = new ArrayList<>(header.size());
			List<String> cmlLine = new ArrayList<>(header.size());
			
			incrLine.add((float)refMFD.getX(i)+"");
			cmlLine.add((float)cmlRef.getX(i)+"");
			
			for (IncrementalMagFreqDist mfd : csvFuncs) {
				if (i >= mfd.size())
					incrLine.add("0.0");
				else
					incrLine.add((float)mfd.getY(i)+"");
			}
			for (EvenlyDiscretizedFunc mfd : csvCmlFuncs) {
				if (i >= mfd.size())
					cmlLine.add("0.0");
				else
					cmlLine.add((float)mfd.getY(i)+"");
			}
			
			incrCSV.addLine(incrLine);
			cmlCSV.addLine(cmlLine);
		}
		
		incrCSV.writeToFile(new File(outputDir, prefix+"_incr.csv"));
		cmlCSV.writeToFile(new File(outputDir, prefix+"_cml.csv"));
	}
	
	static EvenlyDiscretizedFunc getCmlAsFakeIncr(IncrementalMagFreqDist incrMFD) {
		EvenlyDiscretizedFunc rawCML = incrMFD.getCumRateDistWithOffset();
		IncrementalMagFreqDist scaledIncr = incrMFD.deepClone();
		int pinToIndex = incrMFD.getClosestXIndex(5.05); 
		scaledIncr.scaleToIncrRate(pinToIndex, rawCML.getY(pinToIndex));
		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(rawCML.getMinX(), rawCML.size(), rawCML.getDelta());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, scaledIncr.getY(i));
		return ret;
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	private static final Color TRANS_LIGHT_RED = new Color(255, 0, 0, 30);
	private static final Color TRANS_LIGHT_GRAY = new Color(0, 0, 0, 50);
	private static Color DIST_TRANS_COLOR = TRANS_LIGHT_RED;

	private static void addSolRegionalBounds(MFDType type, FaultSystemSolution sol, Range xRange, Region region,
			List<XY_DataSet> incrFuncs, List<PlotCurveCharacterstics> incrChars, List<XY_DataSet> cmlFuncs,
			List<PlotCurveCharacterstics> cmlChars, List<String> csvLabels, List<IncrementalMagFreqDist> csvFuncs,
			EvenlyDiscretizedFunc refMFD) {
		// add full model bounds
		BranchRegionalMFDs branchMFDs = sol.getModule(BranchRegionalMFDs.class);
		List<Region> rois = sol.getRupSet().getModule(RegionsOfInterest.class).getRegions();
		if (branchMFDs != null && branchMFDs.hasRegionalBranchMFDs(type) && rois != null) {
			for (int i=0; i<rois.size(); i++) {
				Region testReg = rois.get(i);
				if (region.equalsRegion(testReg)) {
					System.out.println("Found region match!");

					IncrementalMagFreqDist[] incrPercentiles = branchMFDs.calcRegionalIncrementalFractiles(
							type, i, SolMFDPlot.standardFractiles);
					for (int p=0; p<incrPercentiles.length; p++) {
						csvLabels.add("p"+oDF.format(SolMFDPlot.standardFractiles[p]*100d));
						// need to remap
						IncrementalMagFreqDist rawMFD = incrPercentiles[p];
						IncrementalMagFreqDist remappedMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
						double min = refMFD.getMinX()-0.5*refMFD.getDelta();
						double max = refMFD.getMaxX()+0.5*refMFD.getDelta();
						for (int j=0; j<rawMFD.size(); j++) {
							double mag = rawMFD.getX(j);
							if ((float)mag >= (float)min && (float)mag <= (float)max)
								remappedMFD.set(refMFD.getClosestXIndex(mag), rawMFD.getY(j));
						}
						csvFuncs.add(remappedMFD);
					}

					PlotCurveCharacterstics minMaxChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, DIST_TRANS_COLOR);

					for (UncertainIncrMagFreqDist bounds : SolMFDPlot.processIncrFractiles(incrPercentiles)) {
						incrFuncs.add(bounds);
						incrChars.add(minMaxChar);
					}

					EvenlyDiscretizedFunc[] cmlPercentiles = branchMFDs.calcRegionalCumulativeFractiles(
							type, i, SolMFDPlot.standardFractiles);
					for (UncertainArbDiscFunc cmlBounds : SolMFDPlot.processCmlFractiles(cmlPercentiles, xRange.getLowerBound())) {
						cmlFuncs.add(cmlBounds);
						cmlChars.add(minMaxChar);
					}
				}
			}
		}
	}
	
	private static void addMFDs(IncrementalMagFreqDist incrMFD, PlotCurveCharacterstics pChar,
			String name, List<XY_DataSet> incrFuncs, List<PlotCurveCharacterstics> incrChars,
			List<XY_DataSet> cmlFuncs, List<PlotCurveCharacterstics> cmlChars,
			List<String> csvLabels, List<IncrementalMagFreqDist> csvFuncs) {
		csvLabels.add(name);
		csvFuncs.add(incrMFD);
		
		incrMFD.setName(name);
		incrFuncs.add(incrMFD);
		incrChars.add(pChar);
		
		EvenlyDiscretizedFunc cmlMFD = incrMFD.getCumRateDistWithOffset();
		cmlMFD.setName(name);
		cmlFuncs.add(cmlMFD);
		cmlChars.add(pChar);
	}

	private static Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> calcSolMFDs(
			FaultSystemSolution sol, NSHM23_BaseRegion[] regions, EvenlyDiscretizedFunc refMFD) throws IOException {
		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> ret = HashBasedTable.create();
		
		for (NSHM23_BaseRegion aReg : regions) {
			Region region = aReg.load();
			
			IncrementalMagFreqDist incrMFD = sol.calcNucleationMFD_forRegion(region,
					refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
			
			ret.put(aReg, MFDType.SUPRA_ONLY, incrMFD);
		}
		
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		if (gridProv != null) {
			float mfdMin = (float)(refMFD.getMinX()-0.5*refMFD.getDelta());
			float mfdMax = (float)(refMFD.getMaxX()+0.5*refMFD.getDelta());
			for (NSHM23_BaseRegion aReg : regions) {
				Region region = aReg.load();
				
				IncrementalMagFreqDist gridMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				GriddedRegion gridReg = gridProv.getGriddedRegion();
				for (int i=0; i<gridReg.getNodeCount(); i++) {
					IncrementalMagFreqDist nodeMFD = gridProv.getMFD(i);
					if (nodeMFD == null)
						continue;
					if (region.contains(gridReg.getLocation(i))) {
						for (int m=0; m<nodeMFD.size(); m++) {
							double mag = nodeMFD.getX(m);
							if ((float)mag <= mfdMax && (float)mag >= mfdMin)
								gridMFD.add(gridMFD.getClosestXIndex(mag), nodeMFD.getY(m));
						}
					}
				}
//				System.out.println("Total gridded MFD for "+aReg+": "+gridMFD);
				if (gridMFD.calcSumOfY_Vals() > 0d) {
					ret.put(aReg, MFDType.GRID_ONLY, gridMFD);
					
					IncrementalMagFreqDist supraMFD = ret.get(aReg, MFDType.SUPRA_ONLY);
					
					ret.put(aReg, MFDType.SUM, sum(supraMFD, gridMFD));
				}
			}
		}
		
		return ret;
	}
	
	private static Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> calcModelMFDs(
			HazardModel model, Set<TectonicRegionType> trts, NSHM23_BaseRegion[] regions, EvenlyDiscretizedFunc refMFD) throws IOException {
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		NshmErf faultERF = new NshmErf(model, trts, IncludeBackgroundOption.EXCLUDE);
		faultERF.getTimeSpan().setDuration(1d);
		faultERF.updateForecast();
		
		IncrementalMagFreqDist[] faultMFDs = erfMFDCalc(regions, faultERF, refMFD, exec);
		
		NshmErf gridERF = new NshmErf(model, trts, IncludeBackgroundOption.ONLY);
		gridERF.getTimeSpan().setDuration(1d);
		gridERF.updateForecast();
		
		IncrementalMagFreqDist[] gridMFDs = erfMFDCalc(regions, gridERF, refMFD, exec);
		
		exec.shutdown();

		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> ret = HashBasedTable.create();
		for (int r=0; r<regions.length; r++) {
			ret.put(regions[r], MFDType.SUPRA_ONLY, faultMFDs[r]);
			ret.put(regions[r], MFDType.GRID_ONLY, gridMFDs[r]);
			ret.put(regions[r], MFDType.SUM, sum(faultMFDs[r], gridMFDs[r]));
		}
		
		return ret;
	}
	
	private static Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> calcCascadiaModelMFDs(
			HazardModel model, NSHM23_BaseRegion[] regions, EvenlyDiscretizedFunc refMFD) throws IOException {
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		NshmErf faultERF = new NshmErf(model, Set.of(TectonicRegionType.SUBDUCTION_INTERFACE), IncludeBackgroundOption.EXCLUDE);
		faultERF.getTimeSpan().setDuration(1d);
		faultERF.updateForecast();
		
		IncrementalMagFreqDist[] faultMFDs = erfMFDCalc(regions, faultERF, refMFD, exec);
		
		NshmErf gridERF = new NshmErf(model, Set.of(TectonicRegionType.SUBDUCTION_SLAB), IncludeBackgroundOption.EXCLUDE);
		gridERF.getTimeSpan().setDuration(1d);
		gridERF.updateForecast();
		
		IncrementalMagFreqDist[] gridMFDs = erfMFDCalc(regions, gridERF, refMFD, exec);
		
		exec.shutdown();

		Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> ret = HashBasedTable.create();
		for (int r=0; r<regions.length; r++) {
			ret.put(regions[r], MFDType.SUPRA_ONLY, faultMFDs[r]);
			ret.put(regions[r], MFDType.GRID_ONLY, gridMFDs[r]);
			ret.put(regions[r], MFDType.SUM, sum(faultMFDs[r], gridMFDs[r]));
		}
		
		return ret;
	}
	
	private static IncrementalMagFreqDist[] erfMFDCalc(NSHM23_BaseRegion[] regions, AbstractERF erf,
			EvenlyDiscretizedFunc refMFD, ExecutorService exec) throws IOException {
		IncrementalMagFreqDist[] ret = new IncrementalMagFreqDist[regions.length];
		Region[] regionInstances = new Region[regions.length];
		for (int i=0; i<ret.length; i++) {
			ret[i] = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			regionInstances[i] = regions[i].load();
		}
		List<Future<List<double[]>>> fractInRegCalls = new ArrayList<>();
		
		ConcurrentMap<RuptureSurface, SurfInRegionsResult> subSurfInsideCache = new ConcurrentHashMap<>();
		Deque<Region[]> regDeque = new ArrayDeque<>();
		
		for (ProbEqkSource source : erf)
			fractInRegCalls.add(exec.submit(new SourceFractInRegCall(source, regionInstances,
					subSurfInsideCache, regDeque)));
		
		for (int sourceID=0; sourceID<fractInRegCalls.size(); sourceID++) {
			List<double[]> fractsInReg;
			try {
				fractsInReg = fractInRegCalls.get(sourceID).get();
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				throw ExceptionUtils.asRuntimeException(e);
			}
			
			for (int rupID=0; rupID<fractsInReg.size(); rupID++) {
				ProbEqkRupture rup = erf.getRupture(sourceID, rupID);
				double mag = rup.getMag();
				int magIndex = refMFD.getClosestXIndex(mag);
				double rate = rup.getMeanAnnualRate(1d);
				
				double[] fracts = fractsInReg.get(rupID);
				
				for (int r=0; r<regions.length; r++) {
					if (fracts[r] > 0d)
						ret[r].add(magIndex, rate*fracts[r]);
				}
			}
		}
		
		return ret;
	}
	
	private static class SurfInRegionsResult {
		final int numLocs;
		final int[] inside;
		public SurfInRegionsResult(int numLocs, int[] inside) {
			super();
			this.numLocs = numLocs;
			this.inside = inside;
		}
	}
	
	private static class SourceFractInRegCall implements Callable<List<double[]>> {
		
		private ProbEqkSource source;
		private Region[] regions;
		private ConcurrentMap<RuptureSurface, SurfInRegionsResult> subSurfInsideCache;
		private Deque<Region[]> regDeque;

		public SourceFractInRegCall(ProbEqkSource source, Region[] regions,
				ConcurrentMap<RuptureSurface, SurfInRegionsResult> subSurfInsideCache,
				Deque<Region[]> regDeque) {
			this.source = source;
			this.regions = regions;
			this.subSurfInsideCache = subSurfInsideCache;
			this.regDeque = regDeque;
		}

		@Override
		public List<double[]> call() throws Exception {
			// get a thread-specific version of the regions
			// regions use area, which uses synchronized vectors, and region.clone() won't get around this
			Region[] regions = null;
			synchronized (regDeque) {
				if (!regDeque.isEmpty())
					regions = regDeque.pop();
			}
			if (regions == null) {
				// build one
				regions = new Region[this.regions.length];
				for (int r=0; r<regions.length; r++)
					regions[r] = new Region(this.regions[r].getBorder(), null);
			}
			
			List<double[]> ret = new ArrayList<>();
			
			for (ProbEqkRupture rup : source) {
				RuptureSurface surf = rup.getRuptureSurface();
				int[] countsInside = new int[regions.length];
				int count = 0;
				if (surf instanceof CompoundSurface) {
					for (RuptureSurface subSurf : ((CompoundSurface)surf).getSurfaceList()) {
						SurfInRegionsResult cached = subSurfInsideCache.get(subSurf);
						if (cached == null) {
							int subCount = 0;
							int[] subInsides = new int[regions.length];
							for (Location loc : subSurf.getEvenlyDiscritizedListOfLocsOnSurface()) {
								subCount++;
								for (int r=0; r<regions.length; r++)
									if (quickContains(regions[r], loc))
										subInsides[r]++;
							}
							cached = new SurfInRegionsResult(subCount, subInsides);
							subSurfInsideCache.putIfAbsent(subSurf, cached);
						}
						count += cached.numLocs;
						for (int r=0; r<regions.length; r++)
							countsInside[r] += cached.inside[r];
					}
				} else {
					for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
						count++;
						for (int r=0; r<regions.length; r++)
							if (quickContains(regions[r], loc))
								countsInside[r]++;
					}
				}
				
				double[] fracts = new double[regions.length];
				for (int r=0; r<fracts.length; r++) {
					Preconditions.checkState(countsInside[r] <= count);
					fracts[r] = (double)countsInside[r]/(double)count;
				}
				ret.add(fracts);
			}
			
			synchronized (regDeque) {
				regDeque.push(regions);
			}
			
			return ret;
		}
		
		private boolean quickContains(Region region, Location loc) {
			if (loc.lat < region.getMinLat() || loc.lat > region.getMaxLat()
					|| loc.lon < region.getMinLon() || loc.lon > region.getMaxLon())
				return false;
			return region.contains(loc);
		}
		
	}
	
	private static IncrementalMagFreqDist sum(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		Preconditions.checkState(mfd1.size() == mfd2.size());
		Preconditions.checkState(mfd1.getMinX() == mfd2.getMinX());
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(mfd1.getMinX(), mfd1.size(), mfd1.getDelta());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, mfd1.getY(i) + mfd2.getY(i));
		return ret;
	}

}
