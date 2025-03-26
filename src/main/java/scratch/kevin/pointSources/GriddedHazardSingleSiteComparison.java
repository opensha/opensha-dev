package scratch.kevin.pointSources;

import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_WUS_BranchAveragedERF;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedFiniteRuptureSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.utils.PointSourceDistanceCorrections;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.gmm.Gmm;
import net.mahdilamb.colormap.Colors;

public class GriddedHazardSingleSiteComparison {

	public static void main(String[] args) throws IOException {
		BaseFaultSystemSolutionERF erf = new NSHM23_WUS_BranchAveragedERF();
		erf.updateForecast();
		
//		PointSourceDistanceCorrections corrType = PointSourceDistanceCorrections.ANALYTICAL_TWENTY_POINT;
//		PointSourceDistanceCorrections corrType = PointSourceDistanceCorrections.ANALYTICAL_FIVE_POINT;
		PointSourceDistanceCorrections corrType = PointSourceDistanceCorrections.SUPERSAMPLING_0p1_FIVE_POINT_RJB_DIST;
//		PointSourceDistanceCorrections corrType = PointSourceDistanceCorrections.ANALYTICAL_MEDIAN;
//		GridCellSupersamplingSettings ssSettings = GridCellSupersamplingSettings.DEFAULT;
		GridCellSupersamplingSettings ssSettings = null;
		GridCellSupersamplingSettings compSSSettings = GridCellSupersamplingSettings.DEFAULT.forApplyToFinite(true);
//		GridCellSupersamplingSettings compSSSettings = null;
		int numFinite = 12;
		String title = corrType.toString();
		
		// switch to single-mechanism
////		FocalMech mech = FocalMech.STRIKE_SLIP;
//		FocalMech mech = FocalMech.REVERSE;
//		title += ", "+mech.toString()+" Only";
//		FaultSystemSolution sol = erf.getSolution();
//		sol.setGridSourceProvider(getAsOnlyMech(sol.requireModule(GridSourceList.class), mech));
//		erf = new BaseFaultSystemSolutionERF();
//		erf.setSolution(sol);
		
		GriddedSeismicitySettings settings = erf.getGriddedSeismicitySettings();
		
		settings = settings.forSurfaceType(BackgroundRupType.POINT)
				.forDistanceCorrections(corrType)
				.forMinimumMagnitude(5d).forPointSourceMagCutoff(5d)
				.forSupersamplingSettings(ssSettings);
		
		GriddedSeismicitySettings compSettings = settings.forSurfaceType(BackgroundRupType.FINITE)
				.forDistanceCorrections(PointSourceDistanceCorrections.NONE)
				.forFiniteRuptureSettings(GriddedFiniteRuptureSettings.DEFAULT_CROSSHAIR.forNumSurfaces(numFinite))
				.forSupersamplingSettings(compSSSettings);
		
		ScalarIMR gmm = new NSHMP_GMM_Wrapper(Gmm.COMBINED_ACTIVE_CRUST_2023);
		gmm.setIntensityMeasure(PGA_Param.NAME);
		
//		Site site = new Site(new Location(37.6, -118.9));
		Site site = new Site(new Location(37.7, -119.2));
		site.addParameterList(gmm.getSiteParams());
		
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmm.getIntensityMeasure());
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		
//		ScalarIMR gmm = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.get();
		System.out.println("Calculating exclude curve");
		calc.getHazardCurve(logXVals, site, gmm, erf);
		
		DiscretizedFunc excludeCurve = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			excludeCurve.set(xVals.getX(i), logXVals.getY(i));
		
		System.out.println("DONE");
		
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.setGriddedSeismicitySettings(settings);
		erf.setCacheGridSources(settings.surfaceType == BackgroundRupType.POINT || numFinite <= 4
				|| (numFinite < 10 && settings.supersamplingSettings == null));
		erf.updateForecast();
		
		System.out.println("Calculating for gridded settings: "+settings);
		System.out.println("	We have "+erf.getTotNumRups()+" ruptures");
		
		calc.getHazardCurve(logXVals, site, gmm, erf);
		
		DiscretizedFunc primaryCurve = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			primaryCurve.set(xVals.getX(i), 1d - (1d-logXVals.getY(i))*(1d-excludeCurve.getY(i)));
		
		System.out.println("DONE");
		
		double primary2In50 = primaryCurve.getFirstInterpolatedX_inLogXLogYDomain(ReturnPeriods.TWO_IN_50.oneYearProb); 
		
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.setGriddedSeismicitySettings(compSettings);
		erf.setCacheGridSources(compSettings.surfaceType == BackgroundRupType.POINT || numFinite <= 4
				|| (numFinite < 10 && compSettings.supersamplingSettings == null));
		erf.updateForecast();
		
		System.out.println("Calculating for comp gridded settings: "+compSettings);
		System.out.println("	We have "+erf.getTotNumRups()+" ruptures");
		
		calc.getHazardCurve(logXVals, site, gmm, erf);
		
		DiscretizedFunc compCurve = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			compCurve.set(xVals.getX(i), 1d - (1d-logXVals.getY(i))*(1d-excludeCurve.getY(i)));
		
		System.out.println("DONE");
		
		double comp2In50 = compCurve.getFirstInterpolatedX_inLogXLogYDomain(ReturnPeriods.TWO_IN_50.oneYearProb);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		primaryCurve.setName(settings.surfaceType == BackgroundRupType.POINT ? "Point Surfaces" : "Finite Surfaces");
		compCurve.setName(compSettings.surfaceType == BackgroundRupType.POINT ? "Comparision Point Surfaces" : "Comparison Finite Surfaces");
		
		funcs.add(primaryCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
		
		funcs.add(compCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
		
		excludeCurve.setName("Fault-only");
		funcs.add(excludeCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Colors.tab_grey));
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, gmm.getIntensityMeasure().getName()+" (g)", "Annual Probability of Exceedance");
		plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		DecimalFormat twoDigitsDF = new DecimalFormat("0.00");
		DecimalFormat pDF = new DecimalFormat("0.00%");
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 25);
		XYTextAnnotation ann = new XYTextAnnotation("Primary: "+twoDigitsDF.format(primary2In50)+" g", 2, 0.7);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(annFont);
		plot.addPlotAnnotation(ann);
		ann = new XYTextAnnotation("Comparison: "+twoDigitsDF.format(comp2In50)+" g", 2, 0.4);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(annFont);
		plot.addPlotAnnotation(ann);
		ann = new XYTextAnnotation("Diff: "+twoDigitsDF.format(primary2In50-comp2In50)+" ("+pDF.format((primary2In50-comp2In50)/comp2In50)+")", 2, 0.2325);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		ann.setFont(annFont);
		plot.addPlotAnnotation(ann);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(plot, true, true, new Range(1e-2, 3e0), new Range(1e-5, 1e0));
		
		PlotUtils.writePlots(new File("/tmp"), "gridded_curve_comparison", gp, 1000, 900, true, false, false);
	}
	
	private static GridSourceList getAsOnlyMech(GridSourceList orig, FocalMech mech) {
		List<List<GriddedRupture>> rupsList = new ArrayList<>();
		
		for (int l=0; l<orig.getNumLocations(); l++) {
			List<GriddedRupture> rups = new ArrayList<>();
			rupsList.add(rups);
			
			for (GriddedRupture rup : orig.getRuptures(TectonicRegionType.ACTIVE_SHALLOW, l)) {
				GriddedRuptureProperties props = rup.properties;
				GriddedRuptureProperties modProps = new GriddedRuptureProperties(props.magnitude, mech.rake(), mech.dip(),
						props.strike, props.strikeRange, props.upperDepth, props.lowerDepth, props.length,
						props.hypocentralDepth, props.hypocentralDAS, props.tectonicRegionType);
				rups.add(new GriddedRupture(l, orig.getLocation(l), modProps, rup.rate));
			}
		}
		
		return new GridSourceList.Precomputed(orig.getGriddedRegion(), TectonicRegionType.ACTIVE_SHALLOW, rupsList);
	}

}
