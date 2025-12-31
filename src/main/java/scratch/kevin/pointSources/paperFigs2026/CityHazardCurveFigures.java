package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.PointSourceOptimizedExceedProbCalc;
import org.opensha.sha.calc.RuptureExceedProbCalculator;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.util.SolSiteHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenRelSupplier;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.util.NEHRP_TestCity;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class CityHazardCurveFigures {
	
	public static void main(String[] args) throws IOException {
		double[] periods = {0d, 1d};
		
		Region calcReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		List<NEHRP_TestCity> cities = new ArrayList<>();
		for (NEHRP_TestCity city : NEHRP_TestCity.values())
			if (calcReg.contains(city.location()))
				cities.add(city);
		
		File outputDir = new File(FIGURES_DIR, "site_hazard");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution origSol = FaultSystemSolution.load(ORIG_SOL_FILE);
		
		AttenRelRef rafGmmRef = AttenRelRef.USGS_NSHM23_ACTIVE;
		AttenRelSupplier gmmSupplier = new AttenRelSupplier() {
			
			@Override
			public String getName() {
				return rafGmmRef.getName();
			}
			
			@Override
			public String getShortName() {
				return rafGmmRef.getShortName();
			}
			
			@Override
			public ScalarIMR get() {
				ScalarIMR gmm = rafGmmRef.get();
				gmm.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
				gmm.getParameter(SigmaTruncLevelParam.NAME).setValue(3d);
				return gmm;
			}
		};
		
		ReturnPeriods[] plotRPs = { ReturnPeriods.TEN_IN_50, ReturnPeriods.TWO_IN_50 };
		ReturnPeriods sortRP = ReturnPeriods.TWO_IN_50;
		
		List<BaseFaultSystemSolutionERF> erfs = new ArrayList<>();
		List<String> erfNames = new ArrayList<>();
		List<Color> erfColors = new ArrayList<>();
		
		erfs.add(initERF(origSol, Models.AS_PUBLISHED.getGridProps()));
		erfNames.add("NSHM23 As Published");
		erfColors.add(Colors.tab_brown);
		
		erfs.add(initERF(origSol, PROPOSED_DIST_CORR_MODEL.getGridProps()));
		erfNames.add("Distribution Distance Correction");
		erfColors.add(Colors.tab_blue);
		
		FaultSystemSolution proposedModSol = FaultSystemSolution.load(ORIG_SOL_FILE);
		proposedModSol.setGridSourceProvider(PROPOSED_FULL_MODEL.getGridModFunction().apply(origSol.requireModule(GridSourceList.class)));
		erfs.add(initERF(proposedModSol, PROPOSED_FULL_MODEL.getGridProps()));
		erfNames.add("Proposed Model");
		erfColors.add(Colors.tab_orange);
		
		Supplier<HazardCurveCalculator> calcSupplier = () -> {
			return new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		};
		
		LightFixedXFunc[] xValFuncs = new LightFixedXFunc[periods.length];
		LightFixedXFunc[] logXValFuncs = new LightFixedXFunc[periods.length];
		for (int p=0; p<periods.length; p++) {
			DiscretizedFunc xValsFunc;
			if (periods[p] == 0d)
				xValsFunc = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
			else
				xValsFunc = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
			double[] xVals = new double[xValsFunc.size()];
			double[] logXVals = new double[xVals.length];
			for (int i=0; i<logXVals.length; i++) {
				xVals[i] = xValsFunc.getX(i);
				logXVals[i] = Math.log(xVals[i]);
			}
			xValFuncs[p] = new LightFixedXFunc(xVals, new double[logXVals.length]);
			logXValFuncs[p] = new LightFixedXFunc(logXVals, new double[logXVals.length]);
		}
		
		RuptureExceedProbCalculator exceedCalc = new PointSourceOptimizedExceedProbCalc();
		
		List<CompletableFuture<DiscretizedFunc[][]>> cityFutures = new ArrayList<>();
		
		System.out.println("Calculating for "+cities.size()+" cities");
		for (NEHRP_TestCity city : cities) {
			ScalarIMR gmm = gmmSupplier.get();
			cityFutures.add(CompletableFuture.supplyAsync(() -> {
				System.out.println("Calculating for "+city);
				HazardCurveCalculator calc = calcSupplier.get();
				
				Site site = new Site(city.location());
				site.addParameterList(gmm.getSiteParams());
				
				DiscretizedFunc[][] ret = new DiscretizedFunc[periods.length][erfs.size()];
				
				for (int e=0; e<erfs.size(); e++) {
					for (int p=0; p<periods.length; p++) {
						if (periods[p] == 0d) {
							gmm.setIntensityMeasure(PGA_Param.NAME);
						} else {
							gmm.setIntensityMeasure(SA_Param.NAME);
							SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), periods[p]);
						}
						
						LightFixedXFunc logCurve = new LightFixedXFunc(logXValFuncs[p].getXVals(), new double[logXValFuncs[p].size()]);
						calc.getHazardCurve(logCurve, site, Map.of(TectonicRegionType.ACTIVE_SHALLOW, gmm), erfs.get(e), exceedCalc);
						
						ret[p][e] = new LightFixedXFunc(xValFuncs[p].getXVals(), logCurve.getYVals());
					}
				}
				System.out.println("Done with "+city);
				
				return ret;
			}));
		}
		
		List<DiscretizedFunc[][]> citiesCurves = new ArrayList<>(cities.size());
		
		for (CompletableFuture<DiscretizedFunc[][]> future : cityFutures)
			citiesCurves.add(future.join());
		
		System.out.println("Done calculating");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		Range xRange = new Range(1e-3, 3e0);
		Range yRange = new Range(1e-6, 1e0);

		DecimalFormat oDF = new DecimalFormat("0.#");
		DecimalFormat twoDF = new DecimalFormat("0.00");
		
		for (int p=0; p<periods.length; p++) {
			String perName, perPrefix;
			if (periods[p] == 0d) {
				perName = "PGA";
				perPrefix = "pga";
			} else {
				perName = oDF.format(periods[p])+"s SA";
				perPrefix = "sa_"+oDF.format(periods[p])+"s";
			}
			
			List<CityResult> results = new ArrayList<>(cities.size());
			for (int c=0; c<cities.size(); c++) {
				NEHRP_TestCity city = cities.get(c);
				DiscretizedFunc[] cityCurves = citiesCurves.get(c)[p];
				results.add(new CityResult(city, cityCurves, sortRP));
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				for (int e=0; e<erfs.size(); e++) {
					cityCurves[e].setName(erfNames.get(e));
					funcs.add(cityCurves[e]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, erfColors.get(e)));
				}
				
				List<XYAnnotation> anns = SolSiteHazardCalc.addRPAnnotations(funcs, chars, xRange, yRange, plotRPs, true);
				
				PlotSpec plot = new PlotSpec(funcs, chars, city.toString(), perName+" (g)", "Annual Probability of Exceedance");
				plot.setPlotAnnotations(anns);
				plot.setLegendInset(true);
				
				gp.drawGraphPanel(plot, true, true, xRange, yRange);
				
				String prefix = FileNameUtils.simplify(city.toString())+"_"+perPrefix;
				
				PlotUtils.writePlots(outputDir, prefix, gp, 800, 800, true, true, false);
			}
			
			Collections.sort(results);
			
			// smallest to largest
			System.out.println("Results for "+perName);
			for (CityResult result : results) {
				System.out.println("\t"+result.city+": "+(float)result.change+" (g), "+twoDF.format(result.percentChange)+"%");
			}
		}
	}
	
	private static class CityResult implements Comparable<CityResult> {
		private final NEHRP_TestCity city;
		private final DiscretizedFunc[] curves;
		private final double[] values;
		private final double change;
		private final double percentChange;
		private final double absChange;
		
		private CityResult(NEHRP_TestCity city, DiscretizedFunc[] curves, ReturnPeriods sortRP) {
			super();
			this.city = city;
			this.curves = curves;
			this.values = new double[curves.length];
			for (int i=0; i<curves.length; i++)
				values[i] = getCurveVal(curves[i], sortRP);
			change = values[values.length-1] - values[0];
			absChange = Math.abs(change);
			percentChange = (change*100d)/values[0];
		}

		@Override
		public int compareTo(CityResult o) {
			return Double.compare(absChange, o.absChange);
		}
	}
	
	private static BaseFaultSystemSolutionERF initERF(FaultSystemSolution sol, GriddedSeismicitySettings gridSettings) {
		System.out.println("Initializing ERF for "+gridSettings);
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.setGriddedSeismicitySettings(gridSettings);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		return erf;
	}
	
	private static double getCurveVal(DiscretizedFunc curve, ReturnPeriods rp) {
		if (rp.oneYearProb > curve.getY(0))
			return curve.getX(0);
		if (rp.oneYearProb < curve.getMinY())
			return curve.getMaxX();
		return curve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
	}

}
