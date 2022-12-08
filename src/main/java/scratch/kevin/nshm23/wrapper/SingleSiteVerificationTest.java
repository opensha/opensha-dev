package scratch.kevin.nshm23.wrapper;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSurface;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class SingleSiteVerificationTest {

	public static void main(String[] args) throws IOException {
		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.a.3");
		File solFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_11_10-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File outputDir = new File("/tmp/wrapper_tests");
		Location testLoc = new Location(39, -122);
//		Location testLoc = new Location(37.5, -110);
		
		boolean gridded = true;
		boolean subduction = false;
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NshmErf erf = new NshmErf(erfPath, subduction, gridded);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(PGA_Param.NAME);
		
		Site site = new Site(testLoc);
		site.addParameterList(gmpe.getSiteParams());
		
		double[] distTests = { 1d, 10d, 50d, 100d };
		Site[] distSites = new Site[distTests.length];
		for (int i=0; i<distTests.length; i++) {
			distSites[i] = new Site(LocationUtils.location(testLoc, 0d, distTests[i]));
			distSites[i].addParameterList(gmpe.getSiteParams());
		}
		
		FaultSystemSolution baSol = FaultSystemSolution.load(solFile);
		baSol.removeModuleInstances(RupMFDsModule.class);
		GridSourceProvider gridProv = baSol.getGridSourceProvider();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 0d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		calc.setMaxSourceDistance(500);
		
		GriddedRegion gridReg = gridProv.getGriddedRegion();
		
		int testIndex = gridReg.indexForLocation(testLoc);
		ProbEqkSource wrapperTestSrc = null;
		
		for (ProbEqkSource source : erf) {
			for (ProbEqkRupture rup : source) {
				if (rup.getMag() >= 5d) {
					RuptureSurface surf = rup.getRuptureSurface();
					Location centroid;
					if (surf instanceof NshmSurface)
						centroid = ((NshmSurface)rup.getRuptureSurface()).centroid();
					else
						continue;
					int index = gridReg.indexForLocation(centroid);
					if (index == testIndex) {
						Preconditions.checkState(wrapperTestSrc == source || wrapperTestSrc == null);
						double dist = LocationUtils.linearDistanceFast(centroid, testLoc);
						Preconditions.checkState(dist < 0.01,
								"Wrapper test source distance: %s, centroid=%s, testLoc=%s", dist, centroid, testLoc);
						wrapperTestSrc = source;
					}
				}
			}
		}
		
		calc.getHazardCurve(logXVals, site, gmpe, erf);
		DiscretizedFunc wrapperFullCurve = toLinear(logXVals, xVals);
		calc.getHazardCurve(logXVals, site, gmpe, new SingleSourceERF(wrapperTestSrc));
		DiscretizedFunc wrapperSourceCurve = toLinear(logXVals, xVals);
		DiscretizedFunc[] wrapperSourceDistCurves = new DiscretizedFunc[distSites.length];
		for (int i=0; i<wrapperSourceDistCurves.length; i++) {
			calc.getHazardCurve(logXVals, distSites[i], gmpe, new SingleSourceERF(wrapperTestSrc));
			wrapperSourceDistCurves[i] = toLinear(logXVals, xVals);
		}
		
		gridded = false;
		erf = new NshmErf(erfPath, subduction, gridded);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		calc.getHazardCurve(logXVals, site, gmpe, erf);
		DiscretizedFunc wrapperFaultCurve = toLinear(logXVals, xVals);
		
		ProbEqkSource modelTestSrc = gridProv.getSource(testIndex, 1d, false, BackgroundRupType.POINT);
		FaultSystemSolutionERF modelERF = new FaultSystemSolutionERF(baSol);
		modelERF.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		modelERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		modelERF.getTimeSpan().setDuration(1d);
		modelERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		modelERF.updateForecast();
		
		calc.getHazardCurve(logXVals, site, gmpe, modelERF);
		DiscretizedFunc modelFullCurve = toLinear(logXVals, xVals);
		calc.getHazardCurve(logXVals, site, gmpe, new SingleSourceERF(modelTestSrc));
		DiscretizedFunc modelSourceCurve = toLinear(logXVals, xVals);
		DiscretizedFunc[] modelSourceDistCurves = new DiscretizedFunc[distSites.length];
		for (int i=0; i<modelSourceDistCurves.length; i++) {
			calc.getHazardCurve(logXVals, distSites[i], gmpe, new SingleSourceERF(modelTestSrc));
			modelSourceDistCurves[i] = toLinear(logXVals, xVals);
		}
		
		modelERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		modelERF.updateForecast();
		calc.getHazardCurve(logXVals, site, gmpe, modelERF);
		DiscretizedFunc modelFaultCurve = toLinear(logXVals, xVals);
		
		curveCompPlot(outputDir, "full_erfs", "Full ERF Test", wrapperFullCurve, modelFullCurve);
		curveCompPlot(outputDir, "fault_only", "Fault-Only Test", wrapperFaultCurve, modelFaultCurve);
		curveCompPlot(outputDir, "single_source", "Single Source Test", wrapperSourceCurve, modelSourceCurve);
		for (int i=0; i<distTests.length; i++)
			curveCompPlot(outputDir, "single_source_"+(int)distTests[i]+"km",
					"Single Source Test, "+(int)distTests[i]+" km away", wrapperSourceDistCurves[i], modelSourceDistCurves[i]);
		
		multiDistPlot(outputDir, "dists_compare", "Source Dists Compare", distTests, wrapperSourceDistCurves, modelSourceDistCurves);
	}
	
	private static DiscretizedFunc toLinear(DiscretizedFunc logCurve, DiscretizedFunc xVals) {
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<logCurve.size(); i++)
			ret.set(xVals.getX(i), logCurve.getY(i));
		return ret;
	}
	
	private static class SingleSourceERF extends AbstractERF {
		private ProbEqkSource source;

		public SingleSourceERF(ProbEqkSource source) {
			this.source = source;
			
		}

		@Override
		public int getNumSources() {
			return 1;
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			Preconditions.checkState(idx == 0);
			return source;
		}

		@Override
		public void updateForecast() {}

		@Override
		public String getName() {
			return "Single source";
		}
	}
	
	private static void curveCompPlot(File outputDir, String prefix, String title,
			DiscretizedFunc wrapperCurve, DiscretizedFunc modelCurve) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		modelCurve.setName("Model Curve");
		funcs.add(modelCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		wrapperCurve.setName("Wrapper Curve");
		funcs.add(wrapperCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		
		double minX = Double.POSITIVE_INFINITY;
		double maxX = 0d;
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0d;
		
		for (DiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 0) {
					minX = Math.min(minX, pt.getX());
					maxX = Math.max(maxX, pt.getX());
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		minX = Math.pow(10, Math.floor(Math.log10(minX)));
		maxX = Math.pow(10, Math.ceil(Math.log10(maxX)));
		minY = Math.max(1e-10, Math.pow(10, Math.floor(Math.log10(minY))));
		maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
		
		Range xRange = new Range(minX, maxX);
		Range yRange = new Range(minY, maxY);
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "PGA (g)", "Annual Probability of Exceedance");
		spec.setLegendInset(true);
		
		// now % difference
		DiscretizedFunc pDiff = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<wrapperCurve.size(); i++) {
			double x = wrapperCurve.getX(i);
			double y1 = wrapperCurve.getY(i);
			double y2 = modelCurve.getY(i);
			
			if (y1 > 0 && y2 > 0)
				pDiff.set(x, 100d*(y1-y2)/y2);
		}
		
		funcs = new ArrayList<>();
		chars = new ArrayList<>();
		funcs.add(pDiff);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		Range diffRange = new Range(-30, 30);
		
		PlotSpec rangeSpec = new PlotSpec(funcs, chars, title, spec.getXAxisLabel(), "Wrapper / Model, % Difference");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(List.of(spec, rangeSpec), List.of(true), List.of(true, false),
				List.of(xRange), List.of(yRange, diffRange));
		
		PlotUtils.setSubPlotWeights(gp, 5, 2);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 1000, 1200, true, false, false);
	}
	
	private static void multiDistPlot(File outputDir, String prefix, String title,
			double[] dists, DiscretizedFunc[] wrapperCurves, DiscretizedFunc[] modelCurves) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		PlotLineType[] lineTypes = { PlotLineType.SOLID, PlotLineType.DASHED, PlotLineType.DOTTED, PlotLineType.DOTTED_AND_DASHED };
		
		for (int i=0; i<dists.length; i++) {
			DiscretizedFunc curve = modelCurves[i];
			curve.setName((i == 0 ? "Model Curve, " : "")+(int)dists[i]+"km");
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(lineTypes[i % lineTypes.length], 3f, Color.BLACK));
		}
		for (int i=0; i<dists.length; i++) {
			DiscretizedFunc curve = wrapperCurves[i];
			curve.setName((i == 0 ? "Wrapper Curve, " : "")+(int)dists[i]+"km");
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(lineTypes[i % lineTypes.length], 3f, Color.RED));
		}
		
		double minX = Double.POSITIVE_INFINITY;
		double maxX = 0d;
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0d;
		
		for (DiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 0) {
					minX = Math.min(minX, pt.getX());
					maxX = Math.max(maxX, pt.getX());
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		minX = Math.pow(10, Math.floor(Math.log10(minX)));
		maxX = Math.pow(10, Math.ceil(Math.log10(maxX)));
		minY = Math.max(1e-10, Math.pow(10, Math.floor(Math.log10(minY))));
		maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
		
		Range xRange = new Range(minX, maxX);
		Range yRange = new Range(minY, maxY);
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "PGA (g)", "Annual Probability of Exceedance");
		spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 1000, 850, true, false, false);
	}

}
