package scratch.kevin.prvi25;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.util.Precision;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;

public class SlabTestModGridProvWriter {

	public static void main(String[] args) throws IOException {
		/*
		 * to figure out what is driving slab hazard changes, write out the following:
		 * 
		 * Rates	| PDFs		| Depths	| Mmax
		 * 2003		| 2003		| 2003		| 7.2		| 2003
		 * 2003		| 2003		| 2025		| 7.2		| 2003_at_slab2
		 * 2003		| 2025		| 2025		| 8.0		| 2025_scaled_to_2003_rate
		 * 2025		| 2025		| 2025		| 7.2		| 2025_mmax_7p2
		 * 2003		| 2025		| 2025		| 7.2		| 2025_scaled_to_2003_rate_mmax_7p2
		 */
		
		FaultSystemSolution subSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		File outputDir = new File(INV_DIR, SUBDUCTION_DIR.getName()+"-slab_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Region calcReg = PRVI25_RegionLoader.loadPRVI_MapExtents();
//		double hazResolution = 0.2;
		double hazResolution = 0.05;
//		double hazResolution = 0.025;
		GriddedRegion hazardGrid = new GriddedRegion(calcReg, hazResolution, GriddedRegion.ANCHOR_0_0);
		double[] periods = {0d, 0.2, 1d, 5d};
//		ReturnPeriods[] rps = SolHazardMapCalc.MAP_RPS;
		ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50 };
		AttenRelRef gmmRef = AttenRelRef.USGS_PRVI_SLAB;
		
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		DiscretizedFunc[][] addCurves = calcWithoutSlab(hazardGrid, exec, periods);
		
		GriddedRegion fullGrid = new GriddedRegion(PRVI25_SeismicityRegions.CRUSTAL.load(), 0.1, GriddedRegion.ANCHOR_0_0);
		
		// write out 2003 as is, and keep track of the overall rate and PDF while we're at it
		HazardModel prevModel = HazardModel.load(Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-prvi-2003-main"));
		// needs to be include for slab and interface
		NshmErf subERF = new NshmErf(prevModel, Set.of(TectonicRegionType.SUBDUCTION_SLAB),
				IncludeBackgroundOption.INCLUDE);
		subERF.getTimeSpan().setDuration(1d);
		subERF.updateForecast();
		
		List<List<GriddedRupture>> gridRups2003 = new ArrayList<>(fullGrid.getNodeCount());
		for (int i=0; i<fullGrid.getNodeCount(); i++)
			gridRups2003.add(new ArrayList<>());
		
		GriddedGeoDataSet pdf2003 = new GriddedGeoDataSet(fullGrid);
		double totRate2003 = 0d;
		double rateSkipped = 0;
		for (ProbEqkSource source : subERF) {
			Preconditions.checkState(source.getTectonicRegionType() == TectonicRegionType.SUBDUCTION_SLAB);
			for (ProbEqkRupture rup : source) {
				if (rup.getMag() < 5d)
					continue;
				RuptureSurface sourceSurf = rup.getRuptureSurface();
				LocationList sourceLocs = sourceSurf.getEvenlyDiscritizedListOfLocsOnSurface();
				Preconditions.checkState(sourceLocs.size() == 1);
				double rupRate = rup.getMeanAnnualRate(1d);
				
				totRate2003 += rupRate;
				double rupRateEach = rupRate / sourceLocs.size();
				for (Location loc : sourceLocs) {
					int index = fullGrid.indexForLocation(loc);
					if (index >= 0) {
						pdf2003.set(index, rupRateEach);
						GriddedRuptureProperties props = new GriddedRuptureProperties(rup.getMag(), 0d, 90d,
								Double.NaN, null, loc.depth, loc.depth, 0d, Double.NaN, Double.NaN, TectonicRegionType.SUBDUCTION_SLAB);
						gridRups2003.get(index).add(new GriddedRupture(index, fullGrid.getLocation(index), props, rupRateEach));
					} else {
						rateSkipped += rupRateEach;
					}
				}
			}
		}
		System.out.println("2003 rate M>5: "+(float)totRate2003);
		if (rateSkipped > 0d)
			System.out.println("Skipped "+rateSkipped+" M5/year");
		GridSourceList gridSources2003 = new GridSourceList.Precomputed(fullGrid, TectonicRegionType.SUBDUCTION_SLAB, gridRups2003);
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		ArrayDeque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();
		
		System.out.println("Calculating 2003 maps");
		GriddedGeoDataSet[][] maps2003 = calcMaps(subSol, gridSources2003, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		
		writeMap(outputDir, "slab_2003", "2003", maps2003, periods, rps);
		
		// first just build 2025 provs
		LogicTreeBranch<LogicTreeNode> slabBranch = new LogicTreeBranch<>(PRVI25_LogicTreeBranch.levelsSubductionGridded);
		slabBranch.setValue(PRVI25_DeclusteringAlgorithms.AVERAGE);
		slabBranch.setValue(PRVI25_SeisSmoothingAlgorithms.AVERAGE);
		slabBranch.setValue(PRVI25_SubductionCaribbeanSeismicityRate.PREFFERRED);
		slabBranch.setValue(PRVI25_SubductionMuertosSeismicityRate.PREFFERRED);
		GridSourceList carSources2025 = PRVI25_GridSourceBuilder.buildSlabGridSourceList(slabBranch, PRVI25_SeismicityRegions.CAR_INTRASLAB);
		GridSourceList mueSources2025 = PRVI25_GridSourceBuilder.buildSlabGridSourceList(slabBranch, PRVI25_SeismicityRegions.MUE_INTRASLAB);
		GridSourceList gridSources2025 = GridSourceList.combine(carSources2025, mueSources2025);
		
		System.out.println("Calculating 2025 maps");
		GriddedGeoDataSet[][] maps2025 = calcMaps(subSol, gridSources2025, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		
		writeMap(outputDir, "slab_2025", "2025", maps2025, periods, rps);
		
		writeMapComparisons(outputDir, "slab_2025_vs_2003", "[2025] vs [2003]", maps2025, maps2003, periods, rps);
		
		// now rescale to use the 2003 PDF
		double sumMappedPDF = 0d;
		double totRate2025 = 0d;
		for (int l=0; l<gridSources2025.getNumLocations(); l++) {
			Location loc = gridSources2025.getLocation(l);
			int mappedIndex = fullGrid.indexForLocation(loc);
			if (mappedIndex >= 0)
				sumMappedPDF += pdf2003.get(mappedIndex);
			for (GriddedRupture rup : gridSources2025.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, l))
				if (rup.properties.magnitude >= 5d)
					totRate2025 += rup.rate;
		}
		System.out.println(sumMappedPDF+" of PDF is mapped");
		
		List<List<GriddedRupture>> gridSources2025pdf2003 = new ArrayList<>();
		double modPDFrate = 0d;
		for (int l=0; l<gridSources2025.getNumLocations(); l++) {
			Location loc = gridSources2025.getLocation(l);
			ImmutableList<GriddedRupture> slabRups = gridSources2025.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, l);
			int mappedIndex = fullGrid.indexForLocation(loc);
			if (mappedIndex < 0) {
				gridSources2025pdf2003.add(new ArrayList<>());
				continue;
			}
			double curRate = 0d;
			for (GriddedRupture rup : slabRups)
				if (rup.properties.magnitude >= 5d)
					curRate += rup.rate;
			double curFract = curRate/totRate2025;
			double targetFract = pdf2003.get(mappedIndex)/sumMappedPDF;
			double scalar = targetFract/curFract;
			List<GriddedRupture> rups = new ArrayList<>();
			for (GriddedRupture rup : slabRups) {
				if (rup.properties.magnitude >= 5d) {
					rups.add(new GriddedRupture(l, loc, rup.properties, rup.rate*scalar));
					modPDFrate += rup.rate*scalar;
				}
			}
			gridSources2025pdf2003.add(rups);
		}
		Preconditions.checkState(Precision.equals(totRate2025, modPDFrate, 0.0001), "%s != %s", (float)totRate2003, (float)modPDFrate);
		GridSourceList gridSources2025pdf03 = new GridSourceList.Precomputed(gridSources2025.getGriddedRegion(),
				TectonicRegionType.SUBDUCTION_SLAB, gridSources2025pdf2003);
		GriddedGeoDataSet[][] maps2025pdf03 = calcMaps(subSol, gridSources2025pdf03, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		
		writeMap(outputDir, "slab_2025_pdf03", "2025 w/ 2003 PDF", maps2025pdf03, periods, rps);
		
		writeMapComparisons(outputDir, "slab_2025_vs_2025_pdf03", "[2025] vs [2025 w/ 2003 PDF]", maps2025, maps2025pdf03, periods, rps);
		
		// now 2003 at Slab2 depths
		List<List<GriddedRupture>> gridRups2003atSlab2 = new ArrayList<>(fullGrid.getNodeCount());
		
		for (int l=0; l<gridRups2003.size(); l++) {
			Location loc = fullGrid.getLocation(l);
			int indexInCar = carSources2025.getGriddedRegion().indexForLocation(loc);
			int indexInMue = mueSources2025.getGriddedRegion().indexForLocation(loc);
			List<GriddedRupture> inRups = gridRups2003.get(l);
			List<GriddedRupture> outRups = new ArrayList<>();
			if (indexInCar >= 0) {
				double depth = carSources2025.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, indexInCar).get(0).properties.upperDepth;
				double rateScalar = indexInMue >= 0 ? 0.5 : 1d;
				for (GriddedRupture inRup : inRups) {
					GriddedRuptureProperties inProps = inRup.properties;
					GriddedRuptureProperties outProps = new GriddedRuptureProperties(inProps.magnitude, 0d, 90d,
							Double.NaN, null, depth, depth, 0d, Double.NaN, Double.NaN, TectonicRegionType.SUBDUCTION_SLAB);
					GriddedRupture outRup = new GriddedRupture(l, loc, outProps, inRup.rate*rateScalar);
					outRups.add(outRup);
				}
			}
			if (indexInMue >= 0) {
				double depth = mueSources2025.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, indexInMue).get(0).properties.upperDepth;
				double rateScalar = indexInCar >= 0 ? 0.5 : 1d;
				for (GriddedRupture inRup : inRups) {
					GriddedRuptureProperties inProps = inRup.properties;
					GriddedRuptureProperties outProps = new GriddedRuptureProperties(inProps.magnitude, 0d, 90d,
							Double.NaN, null, depth, depth, 0d, Double.NaN, Double.NaN, TectonicRegionType.SUBDUCTION_SLAB);
					GriddedRupture outRup = new GriddedRupture(l, loc, outProps, inRup.rate*rateScalar);
					outRups.add(outRup);
				}
			}
			gridRups2003atSlab2.add(outRups);
		}
		GridSourceList gridSources2003atSlab2 = new GridSourceList.Precomputed(fullGrid, TectonicRegionType.SUBDUCTION_SLAB, gridRups2003atSlab2);
		System.out.println("Calculating 2003@Slab2 maps");
		GriddedGeoDataSet[][] maps2003atSlab2 = calcMaps(subSol, gridSources2003atSlab2, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		
		writeMap(outputDir, "slab_2003_at_slab2", "2003@Slab2", maps2003atSlab2, periods, rps);
		
		writeMapComparisons(outputDir, "slab_2025_vs_2003_at_slab2", "[2025] vs [2003@Slab2]", maps2025, maps2003atSlab2, periods, rps);
		writeMapComparisons(outputDir, "slab_2003_at_slab2_vs_2003", "[2003@Slab2] vs [2003]", maps2003atSlab2, maps2003, periods, rps);
		
		PRVI25_GridSourceBuilder.SLAB_MMAX = 7.15;
		GridSourceList combSources2025mmax7p2 = PRVI25_GridSourceBuilder.buildSlabGridSourceList(slabBranch);
		
		System.out.println("Calculating 2025 Mmax=7.2 maps");
		GriddedGeoDataSet[][] maps2025mmax7p2 = calcMaps(subSol, combSources2025mmax7p2, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		
		writeMap(outputDir, "slab_2025_7p2", "2025 Mmax=7.2", maps2025mmax7p2, periods, rps);
		
		writeMapComparisons(outputDir, "slab_2025_7p2_vs_2003", "[2025 Mmax=7.2] vs [2003]", maps2025mmax7p2, maps2003, periods, rps);
		
		writeMapComparisons(outputDir, "slab_2025_vs_2025_7p2", "[2025] vs [2025 Mmax=7.2]", maps2025, maps2025mmax7p2, periods, rps);
		
		double scaleTo2003 = totRate2003/totRate2025;
		System.out.println("Scalar to 2003: "+(float)totRate2003+" / "+(float)totRate2025+" = "+(float)scaleTo2003);
		
		GridSourceList combSources2025scaledTo2003 = getScaled(gridSources2025, scaleTo2003);
		
		GridSourceList combSources2025mmax7p2ScaledTo2003 = getScaled(combSources2025mmax7p2, scaleTo2003);
		
		System.out.println("Calculating 2025 Mmax=7.2 Scale03 maps");
		GriddedGeoDataSet[][] maps2025mmax7p2ScaledTo2003 = calcMaps(subSol, combSources2025mmax7p2ScaledTo2003, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		writeMap(outputDir, "slab_2025_7p2_scale03", "2025 Mmax=7.2", maps2025mmax7p2ScaledTo2003, periods, rps);
		writeMapComparisons(outputDir, "slab_2025_7p2_rate03_vs_2003", "[2025 w/ 2003 Rate & Mmax=7.2] vs [2003]", maps2025mmax7p2ScaledTo2003, maps2003, periods, rps);
		
		System.out.println("Calculating 2025 Scale03 maps");
		GriddedGeoDataSet[][] maps2025ScaledTo2003 = calcMaps(subSol, combSources2025scaledTo2003, gmmRef, gmmDeque, calcDeque, periods, rps, hazardGrid, exec, addCurves);
		writeMap(outputDir, "slab_2025_rate03", "2025 w/ 2003 Rate", maps2025ScaledTo2003, periods, rps);
		writeMapComparisons(outputDir, "slab_2025_rate03_vs_2003", "[2025 w/ 2003 Rate] vs [2003]", maps2025ScaledTo2003, maps2003, periods, rps);
		writeMapComparisons(outputDir, "slab_2025_vs_2025_rate03", "[2025] vs [2025 w/ 2003 Rate]", maps2025, maps2025ScaledTo2003, periods, rps);
		
		exec.shutdown();
	}
	
	private static GridSourceList getScaled(GridSourceList gridSources, double scalar) {
		List<List<GriddedRupture>> scaledRups = new ArrayList<>();
		for (int l=0; l<gridSources.getNumLocations(); l++) {
			List<GriddedRupture> rups = new ArrayList<>();
			for (GriddedRupture rup : gridSources.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, l))
				rups.add(new GriddedRupture(l, rup.location, rup.properties, rup.rate*scalar));
			scaledRups.add(rups);
		}
		return new GridSourceList.Precomputed(gridSources.getGriddedRegion(), TectonicRegionType.SUBDUCTION_SLAB, scaledRups);
	}
	
	private static GriddedGeoDataSet[][] calcMaps(FaultSystemSolution sol, GridSourceList gridSources,
			AttenRelRef gmmRef, ArrayDeque<ScalarIMR> gmmDeque, ArrayDeque<HazardCurveCalculator> calcDeque,
			double[] periods, ReturnPeriods[] rps, GriddedRegion gridReg, ExecutorService exec, DiscretizedFunc[][] addCurves) {
		
		sol.setGridSourceProvider(gridSources);
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.getTimeSpan().setDuration(1d);
		erf.setCacheGridSources(true);
		if (gridReg.getSpacing() < 0.1)
			erf.setGriddedSeismicitySettings(erf.getGriddedSeismicitySettings()
					.forSupersamplingSettings(GridCellSupersamplingSettings.QUICK));
		erf.updateForecast();
		
		DiscretizedFunc[][] curves = calcCurves(erf, gmmRef, gmmDeque, calcDeque, periods, gridReg, exec);
		
		GriddedGeoDataSet[][] ret = new GriddedGeoDataSet[periods.length][rps.length];
		for (int p=0; p<periods.length; p++) {
			for (int r=0; r<rps.length; r++)
				ret[p][r] = new GriddedGeoDataSet(gridReg);
			for (int s=0; s<gridReg.getNodeCount(); s++) {
				DiscretizedFunc curve = curves[p][s];
				if (addCurves != null)
					curve = addCurves(curve, addCurves[p][s]);
				for (int r=0; r<rps.length; r++) {
					double prob = rps[r].oneYearProb;
					double val;
					if (prob > curve.getMaxY())
						val = 0d;
					else if (prob < curve.getMinY())
						val = Math.exp(curve.getMaxX());
					else
						val = Math.exp(curve.getFirstInterpolatedX(prob));
					ret[p][r].set(s, val);
				}
			}
		}
		
		return ret;
	}
	
	private static DiscretizedFunc[][] calcCurves(BaseFaultSystemSolutionERF erf,
			AttenRelRef gmmRef, ArrayDeque<ScalarIMR> gmmDeque, ArrayDeque<HazardCurveCalculator> calcDeque,
			double[] periods, GriddedRegion gridReg, ExecutorService exec) {
		IMT_Info imtInto = new IMT_Info();
		DiscretizedFunc[][] ret = new DiscretizedFunc[periods.length][];
		for (int p=0; p<periods.length; p++) {
			System.out.println("Calculating for "+gridReg.getNodeCount()+" sites, T="+(float)periods[p]);
			List<Future<DiscretizedFunc>> siteFutures = new ArrayList<>(gridReg.getNodeCount());
			
			DiscretizedFunc xVals;
			if (periods[p] == 0d)
				xVals = imtInto.getDefaultHazardCurve(PGA_Param.NAME);
			else
				xVals = imtInto.getDefaultHazardCurve(SA_Param.NAME);
			double[] logXValsArray = new double[xVals.size()];
			for (int i = 0; i < xVals.size(); i++)
				logXValsArray[i] = Math.log(xVals.getX(i));
			DiscretizedFunc logXVals = new LightFixedXFunc(logXValsArray, new double[logXValsArray.length]);
			
			for (int s=0; s<gridReg.getNodeCount(); s++)
				siteFutures.add(exec.submit(new SiteCalcCallable(erf, gmmRef, gmmDeque, calcDeque, periods[p],
						logXVals, gridReg.getLocation(s))));
			
			ret[p] = new DiscretizedFunc[gridReg.getNodeCount()];
			for (int s=0; s<siteFutures.size(); s++) {
				try {
					ret[p][s] = siteFutures.get(s).get();
					System.out.print(".");
					if ((s+1) % 100 == 0 && s < siteFutures.size()-1)
						System.out.println(" "+(s+1)+"/"+siteFutures.size());
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			System.out.println(" DONE");
		}
		return ret;
	}
	
	private static DiscretizedFunc addCurves(DiscretizedFunc... curves) {
		DiscretizedFunc ret = curves[0].deepClone();
		for (int i=1; i<curves.length; i++)
			Preconditions.checkState(ret.size() == curves[0].size());
		for (int i=0; i<ret.size(); i++) {
			double prod = 1d;
			for (DiscretizedFunc curve : curves)
				prod *= 1d - curve.getY(i);
			ret.set(i, 1d - prod);
		}
		return ret;
	}
	
	private static class SiteCalcCallable implements Callable<DiscretizedFunc> {
		
		private BaseFaultSystemSolutionERF erf;
		private ArrayDeque<ScalarIMR> gmmDeque;
		private ArrayDeque<HazardCurveCalculator> calcDeque;
		private double period;
		private Location loc;
		private AttenRelRef gmmRef;
		private DiscretizedFunc logXVals;

		public SiteCalcCallable(BaseFaultSystemSolutionERF erf, AttenRelRef gmmRef, ArrayDeque<ScalarIMR> gmmDeque,
				ArrayDeque<HazardCurveCalculator> calcDeque, double period, DiscretizedFunc logXVals, Location loc) {
			this.erf = erf;
			this.gmmRef = gmmRef;
			this.gmmDeque = gmmDeque;
			this.calcDeque = calcDeque;
			this.period = period;
			this.logXVals = logXVals;
			this.loc = loc;
		}

		@Override
		public DiscretizedFunc call() throws Exception {
			ScalarIMR gmm = null;
			synchronized (gmmDeque) {
				if (!gmmDeque.isEmpty())
					gmm = gmmDeque.pop();
			}
			if (gmm == null)
				gmm = gmmRef.get();
			HazardCurveCalculator calc = null;
			synchronized (calcDeque) {
				if (!calcDeque.isEmpty())
					calc = calcDeque.pop();
			}
			if (calc == null) {
				calc = new HazardCurveCalculator();
				calc.setMaxSourceDistance(TectonicRegionType.SUBDUCTION_SLAB.defaultCutoffDist());
			}
			
			if (period == 0d) {
				gmm.setIntensityMeasure(PGA_Param.NAME);
			} else {
				gmm.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
			}
			
			Site site = new Site(loc);
			site.addParameterList(gmm.getSiteParams());
			
			DiscretizedFunc curve = logXVals.deepClone();
			calc.getHazardCurve(curve, site, gmm, erf);
			
//			double[] ret = new double[rps.length];
//			for (int r=0; r<rps.length; r++) {
//				double prob = rps[r].oneYearProb;
//				if (prob > curve.getMaxY())
//					ret[r] = 0d;
//				else if (prob < curve.getMinY())
//					ret[r] = Math.exp(curve.getMaxX());
//				else
//					ret[r] = Math.exp(curve.getFirstInterpolatedX(prob));
//			}
//			
//			if (FIRST_CURVE_DEBUG) {
//				synchronized (SiteCalcCallable.class) {
//					if (FIRST_CURVE_DEBUG) {
//						System.out.println("\nDEBUG");
//						System.out.println("Log curve: "+curve);
//						for (int r = 0; r < rps.length; r++)
//							System.out.println(rps[r].name() + ": " + (float) ret[r]);
//						System.out.flush();
//						FIRST_CURVE_DEBUG = false;
//					}
//				}
//			}
			
			synchronized (gmmDeque) {
				gmmDeque.push(gmm);
			}
			synchronized (calcDeque) {
				calcDeque.push(calc);
			}
			return curve;
		}
		
	}
	
	private static boolean FIRST_CURVE_DEBUG = true;
	
	private static void writeMap(File outputDir, String prefix, String label,
			GriddedGeoDataSet[][] maps, double[] periods, ReturnPeriods[] rps) throws IOException {
		GeographicMapMaker mapMaker = new GeographicMapMaker(maps[0][0].getRegion());
		
		CPT logCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		
		for (int r=0; r<rps.length; r++) {
			String rpPrefix = prefix+"_"+rps[r].name();
			String rpLabel = rps[r].label;
			
			List<PlotSpec> plots = new ArrayList<>();
			
			List<Range> xRanges = new ArrayList<>();
			List<Range> yRanges = new ArrayList<>();
			List<String> periodNames = new ArrayList<>();
			
			for (int p=0; p<periods.length; p++) {
				GriddedGeoDataSet map = maps[p][r].copy();
				map.log10();
				
				mapMaker.plotXYZData(map, logCPT, rpLabel);
				PlotSpec plot = mapMaker.buildPlot(" ");
				plots.add(plot);
				if (xRanges.isEmpty())
					xRanges.add(mapMaker.getXRange());
				yRanges.add(mapMaker.getYRange());
				
				if (periods[p] == 0d)
					periodNames.add("PGA");
				else
					periodNames.add((float)periods[p]+"s SA");
				
				if (p == 0)
					plot.setTitle(label);
			}
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plots, false, false, xRanges, yRanges);
			
			PlotUtils.addSubplotTitles(gp, periodNames, new Font(Font.SANS_SERIF, Font.PLAIN, 22));
			
			PlotUtils.writePlots(outputDir, rpPrefix, gp, 800, true, true, false, false);
		}
	}
	
	private static void writeMapComparisons(File outputDir, String prefix, String label,
			GriddedGeoDataSet[][] testMaps, GriddedGeoDataSet[][] refMaps, double[] periods, ReturnPeriods[] rps) throws IOException {
		GeographicMapMaker mapMaker = new GeographicMapMaker(refMaps[0][0].getRegion());
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		
		for (int r=0; r<rps.length; r++) {
			String rpPrefix = prefix+"_"+rps[r].name();
			String rpLabel = "% change due to intraslab, "+rps[r].label;
			
			List<PlotSpec> plots = new ArrayList<>();
			
			List<Range> xRanges = new ArrayList<>();
			List<Range> yRanges = new ArrayList<>();
			List<String> periodNames = new ArrayList<>();
			
			for (int p=0; p<periods.length; p++) {
				GriddedGeoDataSet testMap = testMaps[p][r];
				GriddedGeoDataSet refMap = refMaps[p][r];
				
				GriddedGeoDataSet pDiff = new GriddedGeoDataSet(refMap.getRegion());
				for (int i=0; i<pDiff.size(); i++) {
					double ref = refMap.get(i);
					double test = testMap.get(i);
					pDiff.set(i, 100d*(test-ref)/ref);
				}
				
				mapMaker.plotXYZData(pDiff, pDiffCPT, rpLabel);
				PlotSpec plot = mapMaker.buildPlot(" ");
				plots.add(plot);
				if (xRanges.isEmpty())
					xRanges.add(mapMaker.getXRange());
				yRanges.add(mapMaker.getYRange());
				
				if (periods[p] == 0d)
					periodNames.add("PGA");
				else
					periodNames.add((float)periods[p]+"s SA");
				
				if (p == 0)
					plot.setTitle(label);
			}
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plots, false, false, xRanges, yRanges);
			
			PlotUtils.addSubplotTitles(gp, periodNames, new Font(Font.SANS_SERIF, Font.PLAIN, 22));
			
			PlotUtils.writePlots(outputDir, rpPrefix, gp, 800, true, true, false, false);
		}
	}
	
	private static DiscretizedFunc[][] calcWithoutSlab(GriddedRegion gridReg, ExecutorService exec, double[] periods) throws IOException {
		FaultSystemSolution crustalSol = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED);
		
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(crustalSol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.setCacheGridSources(true);
		if (gridReg.getSpacing() < 0.1)
			erf.setGriddedSeismicitySettings(erf.getGriddedSeismicitySettings().forSupersamplingSettings(GridCellSupersamplingSettings.QUICK));
		erf.updateForecast();
		
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		ArrayDeque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();
		
		System.out.println("Calculating crustal curves");
		DiscretizedFunc[][] crustalCurves = calcCurves(erf, AttenRelRef.USGS_PRVI_ACTIVE, gmmDeque, calcDeque, periods, gridReg, exec);
		
		gmmDeque.clear();
		FaultSystemSolution interfaceSol = FaultSystemSolution.load(SUBDUCTION_SOLS_COMBINED);
		GridSourceList gridList = interfaceSol.requireModule(GridSourceList.class);
		List<List<GriddedRupture>> interfaceRups = new ArrayList<>();
		for (int l=0; l<gridList.getNumLocations(); l++)
			interfaceRups.add(gridList.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, l));
		interfaceSol.setGridSourceProvider(new GridSourceList.Precomputed(
				gridList.getGriddedRegion(), TectonicRegionType.SUBDUCTION_INTERFACE, interfaceRups));
		erf.setSolution(interfaceSol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.setCacheGridSources(true);
		if (gridReg.getSpacing() < 0.1)
			erf.setGriddedSeismicitySettings(erf.getGriddedSeismicitySettings().forSupersamplingSettings(GridCellSupersamplingSettings.QUICK));
		erf.updateForecast();
		System.out.println("Calculating interface curves");
		DiscretizedFunc[][] interfaceCurves = calcCurves(erf, AttenRelRef.USGS_PRVI_ACTIVE, gmmDeque, calcDeque, periods, gridReg, exec);
		
		DiscretizedFunc[][] ret = new DiscretizedFunc[periods.length][gridReg.getNodeCount()];
		
		for (int p=0; p<periods.length; p++)
			for (int s=0; s<gridReg.getNodeCount(); s++)
				ret[p][s] = addCurves(crustalCurves[p][s], interfaceCurves[p][s]);
		
		return ret;
	}

}
