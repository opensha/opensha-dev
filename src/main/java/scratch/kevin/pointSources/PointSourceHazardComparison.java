package scratch.kevin.pointSources;

import static java.lang.Math.sin;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;

import org.apache.commons.math3.util.Precision;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.PointEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.PointSource13b;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm.PointSurfaceNshm;
import org.opensha.sha.earthquake.rupForecastImpl.PointToFiniteSource;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr.Type;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.FocalMech;

import com.google.common.base.Preconditions;

public class PointSourceHazardComparison {
	
	enum PointSourceType {
		TRUE_POINT_SOURCE("True Point Source", false, Type.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample) {
				return new PointEqkSource(centerLoc, mfd, 1d, aveRake, aveDip);
			}
		},
		POINT_SOURCE_13b_NO_CORR("PointSource13b (no distance correction)", false, Type.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample) {
				return new PointSource13b(centerLoc, mfd, 1d, new double[] {5.0, 1.0}, mechWtMapForRake(aveRake));
			}
		},
		POINT_SOURCE_13b_NSHMP_CORR("PointSource13b (NSHMP08 distance correction)", false, Type.NSHMP08) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample) {
				return new CorrOverrideProbEqkSource(POINT_SOURCE_13b_NO_CORR.buildSource(
						centerLoc, mfd, aveRake, aveDip, isSupersample), Type.NSHMP08);
			}
		},
		POINT_SOURCE_NSHM("PointSourceNshm", false, Type.NONE) { // none is fine here, it's baked in
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample) {
				return new PointSourceNshm(centerLoc, mfd, 1d, mechWtMapForRake(aveRake));
			}
		},
		APROX_SUPERSAMPLE_POINT_SOURCE_NSHM("Approx. Supersampled PointSourceNshm", false, Type.NONE) { // none is fine here, it's baked in
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample) {
				ProbEqkSource source = POINT_SOURCE_NSHM.buildSource(centerLoc, mfd, aveRake, aveDip, isSupersample);
				if (!isSupersample) {
					// not externally supersampled, approximate it here
					List<ProbEqkRupture> modRups = new ArrayList<>();
					for (ProbEqkRupture rup : source) {
						PointSurfaceNshm ptSurf = (PointSurfaceNshm)rup.getRuptureSurface();
						rup.setRuptureSurface(new SupersampledApproxPointSurfaceNshm(ptSurf, rup.getMag(), 0.1, SUPERSAMPLE_SPACING));
						modRups.add(rup);
					}
					source = new RupListSource(modRups, source);
				}
				return source;
			}
		},
		POINT_TO_FINITE("PointToFiniteSource", true, Type.NONE) {
			private WC1994_MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample) {
				FocalMechanism focalMech = new FocalMechanism(Double.NaN, aveDip, aveRake);
				double dipRad = Math.toRadians(aveDip);
				List<ProbEqkRupture> rups = new ArrayList<>();
				ProbEqkSource source0 = null;
				for (int i=0; i<mfd.size(); i++) {
					double mag = mfd.getX(i);
					double depth = (float)mag <= 6.5f ? 5d : 1d;
					double length = WC94.getMedianLength(mag);
				    double aspectWidth = length / 1.5;
				    double ddWidth = (14.0 - depth) / Math.sin(dipRad);
				    ddWidth = Math.min(aspectWidth, ddWidth);
				    double lower = depth + ddWidth * Math.sin(dipRad);
				    
				    IncrementalMagFreqDist singleMFD = new IncrementalMagFreqDist(mag, mag, 1);
				    singleMFD.set(0, mfd.getY(i));
				    
				    PointToFiniteSource source = new PointToFiniteSource(singleMFD, centerLoc, focalMech, depth, WC94, lower, 1d, singleMFD.getMinX(), false);
				    if (source0 == null)
				    	source0 = source;
				    rups.addAll(source.getRuptureList());
				}
				return new RupListSource(rups, source0);
			}
		};
		
		private String label;
		private boolean stochastic;
		private Type distCorrType;

		private PointSourceType(String label, boolean stochastic, Type distCorrType) {
			this.label = label;
			this.stochastic = stochastic;
			Preconditions.checkNotNull(distCorrType);
			this.distCorrType = distCorrType;
		}
		
		public abstract ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
				double aveRake, double aveDip, boolean isSupersample);
	}
	
	private static Map<FocalMech, Double> mechWtMapForRake(double rake) {
		Map<FocalMech, Double> mechWtMap = new HashMap<>();
		mechWtMap.put(mechForRake(rake), 1d);
		for (FocalMech mech : FocalMech.values())
			if (!mechWtMap.containsKey(mech))
				mechWtMap.put(mech, 0d);
		return mechWtMap;
	}
	
	private static FocalMech mechForRake(double rake) {
		double minDiff = Double.POSITIVE_INFINITY;
		FocalMech match = null;
		for (FocalMech mech : FocalMech.values()) {
			double diff = FaultUtils.getAbsAngleDiff(rake, mech.rake());
			if (rake == 180d || rake == -180d)
				// try 0
				diff = Math.min(diff, FaultUtils.getAbsAngleDiff(0d, mech.rake()));
			if (diff < minDiff) {
				minDiff = diff;
				match = mech;
			}
		}
		Preconditions.checkState(minDiff < 45d);
		return match;
	}
	
	private static class RupListSource extends ProbEqkSource {
		
		private ProbEqkSource source0;
		private List<ProbEqkRupture> rups;
		
		public RupListSource(ProbEqkSource... sources) {
			this.source0 = sources[0];
			rups = new ArrayList<>();
			for (ProbEqkSource source : sources)
				rups.addAll(source.getRuptureList());
		}
		
		public RupListSource(List<ProbEqkRupture> rups, ProbEqkSource source0) {
			this.source0 = source0;
			this.rups = rups;
		}

		@Override
		public LocationList getAllSourceLocs() {
			return source0.getAllSourceLocs();
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return source0.getSourceSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return source0.getMinDistance(site);
		}

		@Override
		public int getNumRuptures() {
			return rups.size();
		}

		@Override
		public ProbEqkRupture getRupture(int nRupture) {
			return rups.get(nRupture);
		}
		
	}
	
	private static class CorrOverrideProbEqkSource extends ProbEqkSource {
		
		private ProbEqkSource source;
		private Type corrType;

		public CorrOverrideProbEqkSource(ProbEqkSource source, Type corrType) {
			this.source = source;
			this.corrType = corrType;
		}

		@Override
		public LocationList getAllSourceLocs() {
			return source.getAllSourceLocs();
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return source.getSourceSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return source.getMinDistance(site);
		}

		@Override
		public int getNumRuptures() {
			return source.getNumRuptures();
		}

		@Override
		public ProbEqkRupture getRupture(int nRupture) {
			ProbEqkRupture rup = source.getRupture(nRupture);
			Preconditions.checkState(rup.getRuptureSurface() instanceof PointSurface);
			((PointSurface)rup.getRuptureSurface()).setDistCorrMagAndType(rup.getMag(), corrType);
			return rup;
		}
		
	}
	
	private static class PointSourceCalcERF extends AbstractERF {
		
		private Type distCorrType;
		private List<ProbEqkSource> sources;
		
		public PointSourceCalcERF(PointSourceType type, Location loc, IncrementalMagFreqDist mfd,
				double rake, double dip, int numPerStochastic) {
			this.distCorrType = type.distCorrType;
			if (type.stochastic) {
				sources = new ArrayList<>(numPerStochastic);
				IncrementalMagFreqDist mfdEach = mfd.deepClone();
				mfdEach.scale(1d/numPerStochastic);
				for (int i=0; i<numPerStochastic; i++)
					sources.add(type.buildSource(loc, mfdEach, rake, dip, false));
			} else {
				sources = List.of(type.buildSource(loc, mfd, rake, dip, false));
			}
		}
		
		public PointSourceCalcERF(PointSourceType type, GriddedRegion locs, IncrementalMagFreqDist mfd,
				double rake, double dip, int numPerStochastic) {
			this.distCorrType = type.distCorrType;
			IncrementalMagFreqDist mfdEach = mfd.deepClone();
			sources = new ArrayList<>(type.stochastic ? numPerStochastic*locs.getNodeCount() : locs.getNodeCount());
			mfdEach.scale(1d/locs.getNodeCount());
			for (Location loc : locs.getNodeList()) {
				if (type.stochastic) {
					IncrementalMagFreqDist subMFDEach = mfdEach.deepClone();
					subMFDEach.scale(1d/numPerStochastic);
					for (int i=0; i<numPerStochastic; i++)
						sources.add(type.buildSource(loc, subMFDEach, rake, dip, true));
				} else {
					sources.add(type.buildSource(loc, mfdEach, rake, dip, true));
				}
			}
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
		public void updateForecast() {}

		@Override
		public String getName() {
			return null;
		}
		
	}
	
	enum GridType {
		COLOCATED("Colocated", "colocated"),
		OFFSET("Offset", "offset"),
		HIGH_RES("High Resolution", "high_res");
		
		private String label;
		private String prefix;

		private GridType(String label, String prefix) {
			this.label = label;
			this.prefix = prefix;
		}
	}
	
	private static final double SUPERSAMPLE_SPACING = 0.01;

	public static void main(String[] args) throws IOException {
//		PointSourceType mainType = PointSourceType.POINT_SOURCE_13b_NSHMP_CORR;
		PointSourceType mainType = PointSourceType.POINT_SOURCE_NSHM;
//		PointSourceType mainType = PointSourceType.APROX_SUPERSAMPLE_POINT_SOURCE_NSHM;
//		PointSourceType mainType = PointSourceType.POINT_TO_FINITE;
		PointSourceType compType = null;
//		PointSourceType compType = PointSourceType.POINT_SOURCE_NSHM;
//		PointSourceType compType = PointSourceType.POINT_SOURCE_13b_NSHMP_CORR;
//		PointSourceType compType = PointSourceType.POINT_SOURCE_13b_NO_CORR;
//		PointSourceType compType = PointSourceType.POINT_TO_FINITE;
		
		AttenRelRef gmmRef = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
		double[] periods = {0d, 1d};
		double[] rakes = {0d, 90d};
		
		ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50 };
		
		Map<Double, Double> rakeToDipMap = new HashMap<>();
		rakeToDipMap.put(0d, 90d);
		rakeToDipMap.put(90d, 50d);
		
		File mainOutputDir = new File("/home/kevin/markdown/nshm23-misc/point_source_corr");
		Preconditions.checkState(mainOutputDir.exists() || mainOutputDir.mkdir());
		
		File outputDir;
		if (compType == null)
			outputDir = new File(mainOutputDir, mainType.name());
		else
			outputDir = new File(mainOutputDir, mainType.name()+"_vs_"+compType.name());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		Region calcRegion = new Region(new Location(0d, 0d), new Location(1d, 1d));
		double spacing = 0.1;
		double halfSpacing = 0.5*spacing;
		GriddedRegion colocatedGridded = new GriddedRegion(calcRegion, spacing, GriddedRegion.ANCHOR_0_0);
		Location gridCenter = new Location(0.5d, 0.5d);
		int gridCenterIndex = colocatedGridded.indexForLocation(gridCenter);
		Preconditions.checkState(LocationUtils.areSimilar(gridCenter, colocatedGridded.getLocation(gridCenterIndex)));
		
		EnumMap<GridType, GriddedRegion> siteGrids = new EnumMap<>(GridType.class);
		siteGrids.put(GridType.COLOCATED, colocatedGridded);
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(0d, 100d, 500);
		List<Location> distFuncLocs = new ArrayList<>();
		for (int i=0; i<distFunc.size(); i++) {
			double dist = distFunc.getX(i);
			if (dist == 0d)
				distFuncLocs.add(gridCenter);
			else
				distFuncLocs.add(LocationUtils.location(gridCenter, 0d, dist));
		}
		
		GriddedRegion offsetCalcGridded = new GriddedRegion(calcRegion, spacing, new Location(halfSpacing, halfSpacing));
		siteGrids.put(GridType.OFFSET, offsetCalcGridded);
		
		siteGrids.put(GridType.HIGH_RES, new GriddedRegion(calcRegion, spacing*0.2, GriddedRegion.ANCHOR_0_0));
		
		int numPerStochasticCentered = 1000;
		int numPerStochasticSupersampled = 50;
		
		Region cellReg = new Region(new Location(gridCenter.lat-halfSpacing, gridCenter.lon-halfSpacing),
				new Location(gridCenter.lat+halfSpacing, gridCenter.lon+halfSpacing));
		GriddedRegion superSampledCell = new GriddedRegion(cellReg, SUPERSAMPLE_SPACING,
				new Location(0.5*SUPERSAMPLE_SPACING, 0.5*SUPERSAMPLE_SPACING));
		
		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(5.05, 30, 0.1);
		mfd.setAllButTotMoRate(mfd.getMinX(), mfd.getMaxX(), 1d, 1d);
//		GriddedRegion calcRegion = new GriddedR
		
		double[] distPlotMags = { 5.05, 6.05, 7.05, mfd.getMaxX() };
		
		List<String> lines = new ArrayList<>();
		
		if (compType == null)
			lines.add("# Point Source Hazard Comparisons, "+mainType.label);
		else
			lines.add("# Point Source Hazard Comparisons, "+mainType.label+" vs "+compType.label);
		lines.add("");
		
		lines.add("Point source hazard calculations using the _"+gmmRef.getName()+"_ GMM.");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Calculation Setup");
		lines.add(topLink); lines.add("");
		
		PlotCurveCharacterstics siteChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, Color.BLACK);
		PlotCurveCharacterstics centerChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 10f, Color.RED);
		PlotCurveCharacterstics superSampledChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, Color.BLUE);
		
		for (boolean offset : new boolean[] {false,true}) {
			GeographicMapMaker mapMaker = new GeographicMapMaker(calcRegion);
			mapMaker.setWriteGeoJSON(false);
			
			List<Location> scatterLocs = new ArrayList<>();
			List<PlotCurveCharacterstics> scatterChars = new ArrayList<>();
			List<String> legendItems = new ArrayList<>();
			List<PlotCurveCharacterstics> legendChars = new ArrayList<>();
			GriddedRegion siteReg = offset ? offsetCalcGridded : colocatedGridded;
			for (Location loc : siteReg.getNodeList()) {
				scatterLocs.add(loc);
				scatterChars.add(siteChar);
				if (legendItems.isEmpty()) {
					if (offset)
						legendItems.add("Offset Site Locations");
					else
						legendItems.add("Colocated Site Locations");
					legendChars.add(scatterChars.get(0));
				}
			}
			for (Location loc : superSampledCell.getNodeList()) {
				scatterLocs.add(loc);
				scatterChars.add(superSampledChar);
			}
			legendItems.add("Supersampled Sources");
			legendChars.add(scatterChars.get(scatterChars.size()-1));
			scatterLocs.add(gridCenter);
			scatterChars.add(centerChar);
			legendItems.add("Source Center");
			legendChars.add(scatterChars.get(scatterChars.size()-1));
			
			legendItems.add("Source Cell");
			legendChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 14f, Color.BLACK));
			mapMaker.plotInsetRegions(List.of(cellReg), new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK), null, 0f);
			
			mapMaker.plotScatters(scatterLocs, scatterChars, null);
			mapMaker.setCustomLegendItems(legendItems, legendChars);
			if (offset)
				mapMaker.plot(resourcesDir, "offset_site_grid", "Offset Site Grid");
			else
				mapMaker.plot(resourcesDir, "colocated_site_grid", "Colocated Site Grid");
		}
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("![colocated]("+resourcesDir.getName()+"/colocated_site_grid.png)",
				"![offset]("+resourcesDir.getName()+"/offset_site_grid.png)");
		lines.addAll(table.build());
		lines.add("");
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		Deque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();
		Deque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		CPT batlowCategorical = GMT_CPT_Files.CATEGORICAL_BATLOW_UNIFORM.instance();
		Color mainColor = batlowCategorical.get(0).minColor;
		Color compColor = batlowCategorical.get(1).minColor;
		
		Range curveYRange = new Range(1e-5, 1e0);
		GeographicMapMaker mapMaker = new GeographicMapMaker(calcRegion);
		mapMaker.setWriteGeoJSON(false);
		
		List<Location> scatterLocs = new ArrayList<>();
		List<PlotCurveCharacterstics> scatterChars = new ArrayList<>();
		scatterLocs.add(gridCenter);
		centerChar = new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 10f, Color.BLACK);
		scatterChars.add(centerChar);
		mapMaker.plotScatters(scatterLocs, scatterChars, null);
		
		mapMaker.plotInsetRegions(List.of(cellReg), new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK), null, 0f);
		
//		List<String> legendItems = new ArrayList<>();
//		List<PlotCurveCharacterstics> legendChars = new ArrayList<>();
//		legendItems.add("Source Center");
//		legendChars.add(scatterChars.get(scatterChars.size()-1));
//		legendItems.add("Source Cell");
//		legendChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 14f, Color.BLACK));
//		
//		mapMaker.setCustomLegendItems(legendItems, legendChars);
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		
		boolean[] falseTrue = {false,true};
		boolean[] falseOnly = {false};
		
		for (double rake : rakes) {
			FocalMech mech = mechForRake(rake);
			String mechLabel = mech.toString();
			String mechPrefix = mech.name();
			
			double dip = rakeToDipMap.get(rake);
			PointSourceCalcERF centerERF = new PointSourceCalcERF(mainType, gridCenter, mfd, rake, dip, numPerStochasticCentered);
			PointSourceCalcERF supersampledERF = new PointSourceCalcERF(mainType, superSampledCell, mfd, rake, dip, numPerStochasticSupersampled);
			
			PointSourceCalcERF compCenterERF = null;
			PointSourceCalcERF compSupersampledERF = null;
			boolean[] compBools;
			if (compType != null) {
				compBools = falseTrue;
				compCenterERF = new PointSourceCalcERF(compType, gridCenter, mfd, rake, dip, numPerStochasticCentered);
				compSupersampledERF = new PointSourceCalcERF(compType, superSampledCell, mfd, rake, dip, numPerStochasticSupersampled);
			} else {
				compBools = falseOnly;
			}
			
			lines.add("## "+mechLabel);
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("DistanceJB", "DistanceRup");
			
			for (double mag : distPlotMags) {
				Range distRange = new Range(0d, distFunc.getMaxX());
				Range ratioRange = new Range(0, 2);
				Range diffRange = new Range(-30d, 10d);
				
				EvenlyDiscretizedFunc centeredRJBFunc = distFunc.deepClone();
				EvenlyDiscretizedFunc centeredRRupFunc = distFunc.deepClone();
				System.out.println("Calculating distances for M"+(float)mag+", "+mechLabel);
				calcDistFuncs(centerERF, distFuncLocs, centeredRJBFunc, centeredRRupFunc, mag);
				EvenlyDiscretizedFunc supersampledRJBFunc = distFunc.deepClone();
				EvenlyDiscretizedFunc supersampledRRupFunc = distFunc.deepClone();
				calcDistFuncs(supersampledERF, distFuncLocs, supersampledRJBFunc, supersampledRRupFunc, mag);
				EvenlyDiscretizedFunc compCenteredRJBFunc = compType == null ? null : distFunc.deepClone();
				EvenlyDiscretizedFunc compCenteredRRupFunc = compType == null ? null : distFunc.deepClone();
				EvenlyDiscretizedFunc compSupersampledRJBFunc = compType == null ? null : distFunc.deepClone();
				EvenlyDiscretizedFunc compSupersampledRRupFunc = compType == null ? null : distFunc.deepClone();
				if (compType != null) {
					calcDistFuncs(compCenterERF, distFuncLocs, compCenteredRJBFunc, compCenteredRRupFunc, mag);
					calcDistFuncs(compSupersampledERF, distFuncLocs, compSupersampledRJBFunc, compSupersampledRRupFunc, mag);
				}
				
				table.initNewLine();
				for (boolean rRup : falseTrue) {
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					EvenlyDiscretizedFunc centeredFunc = rRup ? centeredRRupFunc : centeredRJBFunc;
					EvenlyDiscretizedFunc supersampledFunc = rRup ? supersampledRRupFunc : supersampledRJBFunc;
					EvenlyDiscretizedFunc compCenteredFunc = null;
					EvenlyDiscretizedFunc compSupersampledFunc = null;
					if (compType != null) {
						compCenteredFunc = rRup ? compCenteredRRupFunc : compCenteredRJBFunc;
						compSupersampledFunc = rRup ? compSupersampledRRupFunc : compSupersampledRJBFunc;
					}
					
					centeredFunc.setName("Centered Source");
					funcs.add(asDistRatio(centeredFunc));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mainColor));
					
					if (compType != null) {
						compCenteredFunc.setName("Comparison Cenetered Source");
						funcs.add(asDistRatio(compCenteredFunc));
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
					}
					
					supersampledFunc.setName("Supersampled Source");
					funcs.add(asDistRatio(supersampledFunc));
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
					
					if (compType != null) {
						compSupersampledFunc.setName("Comparison Supersampled Source");
						funcs.add(asDistRatio(compSupersampledFunc));
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
					}
					
					String title = mechLabel+" M"+oDF.format(mag)+" Distance"+(rRup ? "Rup" : "JB")+" Comparison";
					String xAxisLabel = "Horizontal Distance to Source Center (km)";
					PlotSpec ratioSpec = new PlotSpec(funcs, chars, title,
							xAxisLabel, rRup ? "DistanceRup Ratio" : "DistanceJB Ratio");
					ratioSpec.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					funcs = new ArrayList<>();
					chars = new ArrayList<>();
					
					funcs.add(asDistDiff(centeredFunc));
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mainColor));
					
					if (compType != null) {
						funcs.add(asDistDiff(compCenteredFunc));
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
					}
					
					funcs.add(asDistDiff(supersampledFunc));
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
					
					if (compType != null) {
						funcs.add(asDistDiff(compSupersampledFunc));
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
					}
					
					PlotSpec diffSpec = new PlotSpec(funcs, chars, title,
							xAxisLabel, rRup ? "DistanceRup Difference" : "DistanceJB Difference");
					diffSpec.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					HeadlessGraphPanel gp = PlotUtils.initHeadless();
					
					gp.drawGraphPanel(List.of(ratioSpec, diffSpec), false, false, List.of(distRange), List.of(ratioRange, diffRange));
					
					gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
					
					String prefix = mechPrefix+"_m"+oDF.format(mag)+"_"+(rRup ? "rRup" : "rJB");
					PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 900, true, true, false);
					
					table.addColumn("![Curves]("+resourcesDir.getName()+"/"+prefix+".png)");
				}
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
			
			for (double period : periods) {
				String perLabel, perPrefix, perUnits;
				Range curveXRange;
				DiscretizedFunc xVals;
				if (period == -1d) {
					perLabel = "PGV";
					perPrefix = "pgv";
					perUnits = "cm/s";
					xVals = new IMT_Info().getDefaultHazardCurve(PGV_Param.NAME);
					curveXRange = new Range(1e-1, 1e2);
				} else if (period == 0d) {
					perLabel = "PGA";
					perPrefix = "pga";
					perUnits = "g";
					xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
					curveXRange = new Range(1e-2, 1e1);
				} else {
					Preconditions.checkState(period > 0d);
					perLabel = oDF.format(period)+"s SA";
					perPrefix = oDF.format(period)+"s";
					perUnits = "g";
					xVals = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
					curveXRange = new Range(1e-2, 1e1);
				}
				
				CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(
						Math.log10(curveXRange.getLowerBound()), Math.log10(curveXRange.getUpperBound()));
				DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
				for (Point2D pt : xVals)
					logXVals.set(Math.log(pt.getX()), 0d);
				if (periods.length > 1) {
					lines.add("### "+mechLabel+", "+perLabel);
					lines.add(topLink); lines.add("");
				}
				
				String mechPerPrefix = mechPrefix+"_"+perPrefix;
				
				Map<GridType, GriddedGeoDataSet[]> centeredMaps = new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> supersampledMaps = new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> compCenteredMaps = compType == null ? null : new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> compSupersampledMaps = compType == null ? null : new EnumMap<>(GridType.class);
				
				
				for (GridType grid : siteGrids.keySet()) {
					GriddedRegion siteGrid = siteGrids.get(grid);
					
					for (boolean comp : compBools) {
						PointSourceCalcERF myCenterERF = comp ? compCenterERF : centerERF;
						PointSourceCalcERF mySupersampledERF = comp ? compSupersampledERF : supersampledERF;
						System.out.println("Calculating "+grid.label+", centered, "+mechLabel+", "+perLabel+", comp="+comp);
						List<DiscretizedFunc> centeredCurves = calcCurves(siteGrid.getNodeList(), myCenterERF, gmmRef,
								gmmDeque, calcDeque, period, xVals, logXVals, exec);
						
						System.out.println("Calculating "+grid.label+", supersampled, "+mechLabel+", "+perLabel+", comp="+comp);
						List<DiscretizedFunc> supersampledCurves = calcCurves(siteGrid.getNodeList(), mySupersampledERF, gmmRef,
								gmmDeque, calcDeque, period, xVals, logXVals, exec);
						
						GriddedGeoDataSet[] myCenteredMaps = new GriddedGeoDataSet[rps.length];
						GriddedGeoDataSet[] mySupersampledMaps = new GriddedGeoDataSet[rps.length];
						
						for (int r=0; r<rps.length; r++) {
							myCenteredMaps[r] = curvesToMap(siteGrid, centeredCurves, rps[r]);
							mySupersampledMaps[r] = curvesToMap(siteGrid, supersampledCurves, rps[r]);
						}
						
						if (comp) {
							compCenteredMaps.put(grid, myCenteredMaps);
							compSupersampledMaps.put(grid, mySupersampledMaps);
						} else {
							centeredMaps.put(grid, myCenteredMaps);
							supersampledMaps.put(grid, mySupersampledMaps);
						}
					}
				}
				DiscretizedFunc colocatedCurve = null;
				DiscretizedFunc colocatedSupersampledCurve = null;
				List<DiscretizedFunc> distFuncCenteredCurves = null;
				List<DiscretizedFunc> distFuncSupersampledCurves = null;
				
				DiscretizedFunc compColocatedCurve = null;
				DiscretizedFunc compColocatedSupersampledCurve = null;
				List<DiscretizedFunc> compDistFuncCenteredCurves = null;
				List<DiscretizedFunc> compDistFuncSupersampledCurves = null;
				
				for (boolean comp : compBools) {
					PointSourceCalcERF myCenterERF = comp ? compCenterERF : centerERF;
					PointSourceCalcERF mySupersampledERF = comp ? compSupersampledERF : supersampledERF;
					
					System.out.println("Calculating dist func, centered, "+mechLabel+", "+perLabel+", comp="+comp);
					List<DiscretizedFunc> calcDistFuncCenteredCurves = calcCurves(distFuncLocs, myCenterERF, gmmRef,
							gmmDeque, calcDeque, period, xVals, logXVals, exec);
					
					System.out.println("Calculating dist func, supersampled, "+mechLabel+", "+perLabel+", comp="+comp);
					List<DiscretizedFunc> calcDistFuncSupersampledCurves = calcCurves(distFuncLocs, mySupersampledERF, gmmRef,
							gmmDeque, calcDeque, period, xVals, logXVals, exec);
					
					if (comp) {
						compColocatedCurve = calcDistFuncCenteredCurves.get(0);
						compColocatedSupersampledCurve = calcDistFuncSupersampledCurves.get(0);
						compDistFuncCenteredCurves = calcDistFuncCenteredCurves;
						compDistFuncSupersampledCurves = calcDistFuncSupersampledCurves;
					} else {
						colocatedCurve = calcDistFuncCenteredCurves.get(0);
						colocatedSupersampledCurve = calcDistFuncSupersampledCurves.get(0);
						distFuncCenteredCurves = calcDistFuncCenteredCurves;
						distFuncSupersampledCurves = calcDistFuncSupersampledCurves;
					}
				}
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				colocatedCurve.setName("Centered Source");
				funcs.add(colocatedCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mainColor));
				
				if (compType != null) {
					compColocatedCurve.setName("Comparison Cenetered Source");
					funcs.add(compColocatedCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
				}
				
				colocatedSupersampledCurve.setName("Supersampled Source");
				funcs.add(colocatedSupersampledCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
				
				if (compType != null) {
					compColocatedSupersampledCurve.setName("Comparison Supersampled Source");
					funcs.add(compColocatedSupersampledCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, "Zero-Distance Hazard Curves",
						perLabel+" ("+perUnits+")", "Annual Probability of Exceedance");
				spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(spec, true, true, curveXRange, curveYRange);
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				PlotUtils.writePlots(resourcesDir, mechPerPrefix+"_curves", gp, 800, 800, true, true, false);
				
				lines.add("![Curves]("+resourcesDir.getName()+"/"+mechPerPrefix+"_curves.png)");
				lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				for (int r=0; r<rps.length; r++) {
					funcs = new ArrayList<>();
					chars = new ArrayList<>();
					
					EvenlyDiscretizedFunc xy = distValFunc(distFunc, distFuncCenteredCurves, rps[r]);
					xy.setName("Centered Source");
					funcs.add(xy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mainColor));
					
					if (compType != null) {
						EvenlyDiscretizedFunc compXY = distValFunc(distFunc, compDistFuncCenteredCurves, rps[r]);
						compXY.setName("Comparison Centered Source");
						funcs.add(compXY);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
					}
					
					xy = distValFunc(distFunc, distFuncSupersampledCurves, rps[r]);
					xy.setName("Supersampled Source");
					funcs.add(xy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
					
					if (compType != null) {
						EvenlyDiscretizedFunc compXY = distValFunc(distFunc, compDistFuncSupersampledCurves, rps[r]);
						compXY.setName("Comparison Supersampled Source");
						funcs.add(compXY);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
					}
					
					spec = new PlotSpec(funcs, chars, rps[r].label+" vs Distance",
							"Distance (km)", rps[r].label+", "+perLabel+" ("+perUnits+")");
					spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					gp.drawGraphPanel(spec, false, true, new Range(0d, xy.getMaxX()), logBufferedYRange(funcs));
					
					gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
					
					String myPrefix = mechPerPrefix+"_"+rps[r].name()+"_vs_dist";
					PlotUtils.writePlots(resourcesDir, myPrefix, gp, 800, 650, true, true, false);
					
					table.addLine("![Scatter]("+resourcesDir.getName()+"/"+myPrefix+".png)");
				}
				
				lines.addAll(table.build());
				lines.add("");
				
				for (int r=0; r<rps.length; r++) {
					if (rps.length > 1) {
						if (periods.length > 1)
							lines.add("#### "+mechLabel+", "+perLabel+", "+rps[r].label);
						else
							lines.add("### "+mechLabel+", "+rps[r].label);
						lines.add(topLink); lines.add("");
					}
					String mapPrefix = mechPerPrefix+"_"+rps[r].name();
					
					String hazLabel = rps[r].label+", "+perLabel+" ("+perUnits+")";
					String hazChangeLabel = rps[r].label+", "+perLabel+", % Change";
					
					for (GridType grid : siteGrids.keySet()) {
						
						lines.add("__"+grid.label+" Grid__");
						lines.add("");
						
						table = MarkdownUtils.tableBuilder();
						if (compType != null)
							table.addLine(mainType.label, compType.label, "% Change");
						else
							table.initNewLine();
						
						for (boolean supersampled : falseTrue) {
							String mainPrefix = mapPrefix;
							String title;
							GriddedGeoDataSet mainMap;
							GriddedGeoDataSet compMap = null;
							
							mainPrefix += "_"+grid.prefix;
							title = grid.label+" Sites";
							mainMap = supersampled ? supersampledMaps.get(grid)[r] : centeredMaps.get(grid)[r];
							if (compType != null)
								compMap = supersampled ? compSupersampledMaps.get(grid)[r] : compCenteredMaps.get(grid)[r];
							if (supersampled) {
								mainPrefix += "_supersampled";
								title += ", Supersampled Source";
							} else {
								mainPrefix += "_centered";
								title += ", Centered Source";
							}
							
							mapMaker.plotXYZData(log10(mainMap), hazCPT, hazLabel);
							mapMaker.plot(resourcesDir, mainPrefix, title);
							
							if (compType != null)
								table.initNewLine();
							table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+".png)");
							if (compType != null) {
								mapMaker.plotXYZData(log10(compMap), hazCPT, hazLabel);
								mapMaker.plot(resourcesDir, mainPrefix+"_comp", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+"_comp.png)");
								
								GriddedGeoDataSet pDiff = mapPDiff(mainMap, compMap);
								mapMaker.plotXYZData(pDiff, pDiffCPT, hazChangeLabel);
								mapMaker.plot(resourcesDir, mainPrefix+"_comp_pDiff", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+"_comp_pDiff.png)");
								table.finalizeLine();
							}
						}
						
						// now add centered vs colocated
						String title = "Centered vs Supersampled Source";
						String diffPrefix = mapPrefix+"_centered_supersampled_pDiff";
						
						GriddedGeoDataSet pDiff = mapPDiff(centeredMaps.get(grid)[r], supersampledMaps.get(grid)[r]);
						GriddedGeoDataSet compPDiff = null;
						if (compType != null)
							compPDiff = mapPDiff(compCenteredMaps.get(grid)[r], compSupersampledMaps.get(grid)[r]);
						title += ", "+grid.label+" Sites";
						diffPrefix += "_"+grid.prefix;
						
						if (compType != null)
							table.initNewLine();
						mapMaker.plotXYZData(pDiff, pDiffCPT, hazChangeLabel);
						mapMaker.plot(resourcesDir, diffPrefix, title);
						table.addColumn("![map]("+resourcesDir.getName()+"/"+diffPrefix+".png)");
						
						if (compType != null) {
							mapMaker.plotXYZData(compPDiff, pDiffCPT, hazChangeLabel);
							mapMaker.plot(resourcesDir, diffPrefix+"_comp", title);
							table.addColumn("![map]("+resourcesDir.getName()+"/"+diffPrefix+"_comp.png)");
							
							table.addColumn("");
						}
						table.finalizeLine();
						
						lines.addAll(table.build());
						lines.add("");
					}
				}
			}
		}
		
		exec.shutdown();
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static List<DiscretizedFunc> calcCurves(List<Location> siteLocs, PointSourceCalcERF erf,
			AttenRelRef gmmRef, Deque<ScalarIMR> gmmDeque, Deque<HazardCurveCalculator> calcDeque,
			double period, DiscretizedFunc xVals,
			DiscretizedFunc logXVals, ExecutorService exec) {
		List<Future<DiscretizedFunc>> futures = new ArrayList<>(siteLocs.size());
		
		for (Location siteLoc : siteLocs) {
			futures.add(exec.submit(new Callable<DiscretizedFunc>() {

				@Override
				public DiscretizedFunc call() throws Exception {
					HazardCurveCalculator calc = null;
					synchronized (calcDeque) {
						if (!calcDeque.isEmpty())
							calc = calcDeque.pop();
					}
					if (calc == null)
						calc = new HazardCurveCalculator();
					ScalarIMR gmm = null;
					synchronized (gmmDeque) {
						if (!gmmDeque.isEmpty())
							gmm = gmmDeque.pop();
					}
					if (gmm == null)
						gmm = gmmRef.get();
					if (period == -1d) {
						gmm.setIntensityMeasure(PGV_Param.NAME);
					} else if (period == 0d) {
						gmm.setIntensityMeasure(PGA_Param.NAME);
					} else {
						Preconditions.checkState(period > 0d);
						gmm.setIntensityMeasure(SA_Param.NAME);
						SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
					}
					DiscretizedFunc logCurve = logXVals.deepClone();
					
					Site site = new Site(siteLoc);
					site.addParameterList(gmm.getSiteParams());
					
					calc.setPtSrcDistCorrType(erf.distCorrType);
					calc.getHazardCurve(logCurve, site, gmm, erf);
					
					synchronized (calcDeque) {
						calcDeque.push(calc);
					}
					synchronized (gmmDeque) {
						gmmDeque.push(gmm);
					}
					
					DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
					for (int i=0; i<xVals.size(); i++)
						ret.set(xVals.getX(i), logCurve.getY(i));
					return ret;
				}
			}));
		}
		
		List<DiscretizedFunc> ret = new ArrayList<>();
		for (Future<DiscretizedFunc> future : futures) {
			try {
				ret.add(future.get());
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		return ret;
	}
	
	private static GriddedGeoDataSet curvesToMap(GriddedRegion reg, List<DiscretizedFunc> curves, ReturnPeriods rp) {
		Preconditions.checkState(curves.size() == reg.getNodeCount());
		GriddedGeoDataSet map = new GriddedGeoDataSet(reg);
		for (int i=0; i<curves.size(); i++) {
			DiscretizedFunc curve = curves.get(i);
			map.set(i, curveVal(curve, rp));
		}
		return map;
	}
	
	private static double curveVal(DiscretizedFunc curve, ReturnPeriods rp) {
		double val;
		if (rp.oneYearProb > curve.getMaxY())
			val = 0d;
		else if (rp.oneYearProb < curve.getMinY())
			// saturated
			val = curve.getMaxX();
		else
			val = curve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
		return val;
	}
	
	private static EvenlyDiscretizedFunc distValFunc(EvenlyDiscretizedFunc distFunc, List<DiscretizedFunc> curves, ReturnPeriods rp) {
		EvenlyDiscretizedFunc ret = distFunc.deepClone();
		
		for (int i=0; i<ret.size(); i++)
			ret.set(i, curveVal(curves.get(i), rp));
		
		return ret;
	}
	
	private static Range logBufferedYRange(List<? extends XY_DataSet> funcs) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (XY_DataSet func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 0d) {
					min = Math.min(min, pt.getY());
					max = Math.max(max, pt.getY());
				}
			}
		}
		
		double logMax = Math.log10(max);
		if (max < 2.5*Math.pow(10, Math.floor(logMax)))
			max = 3*Math.pow(10, Math.floor(logMax));
		else
			max = Math.pow(10, Math.ceil(logMax));
		double logMin = Math.log10(min);
		if (min > 8.5*Math.pow(10, Math.ceil(logMin)-1))
			min = 8*Math.pow(10, Math.ceil(logMin)-1);
		else
			min = Math.pow(10, Math.floor(logMin));
		
		return new Range(min, max);
	}
	
	private static GriddedGeoDataSet log10(GriddedGeoDataSet map) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map.getRegion());
		
		for (int i=0; i<ret.size(); i++) {
			ret.set(i, Math.log10(map.get(i)));
		}
		return ret;
	}
	
	private static GriddedGeoDataSet mapPDiff(GriddedGeoDataSet numerator, GriddedGeoDataSet denominator) {
		Preconditions.checkState(numerator.size() == denominator.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(numerator.getRegion());
		
		for (int i=0; i<ret.size(); i++) {
			double num = numerator.get(i);
			double den = denominator.get(i);
			ret.set(i, 100d*(num - den)/den);
		}
		return ret;
	}
	
	private static void calcDistFuncs(PointSourceCalcERF erf, List<Location> distLocs,
			EvenlyDiscretizedFunc rJB, EvenlyDiscretizedFunc rRup, double mag) {
		List<ProbEqkRupture> rupsAtMag = new ArrayList<>();
		for (ProbEqkSource source : erf)
			for (ProbEqkRupture rup : source)
				if (Precision.equals(rup.getMag(), mag))
					rupsAtMag.add(rup);
		Preconditions.checkState(!rupsAtMag.isEmpty());
		
		IntStream.range(0, distLocs.size()).parallel().forEach(new IntConsumer() {
			
			@Override
			public void accept(int i) {
				List<Double> rJBs = new ArrayList<>(rupsAtMag.size());
				List<Double> rRups = new ArrayList<>(rupsAtMag.size());
				Location loc = distLocs.get(i);
				for (ProbEqkRupture rup : rupsAtMag) {
					RuptureSurface surf = rup.getRuptureSurface();
					rJBs.add(surf.getDistanceJB(loc));
					rRups.add(surf.getDistanceRup(loc));
				}
				rJB.set(i, rJBs.stream().mapToDouble(D->D).average().getAsDouble());
				rRup.set(i, rRups.stream().mapToDouble(D->D).average().getAsDouble());
			}
		});
	}
	
	private static EvenlyDiscretizedFunc asDistRatio(EvenlyDiscretizedFunc distFunc) {
		EvenlyDiscretizedFunc ret = distFunc.deepClone();
		
		for (int i=0; i<ret.size(); i++)
			ret.set(i, distFunc.getY(i)/distFunc.getX(i));
		
		return ret;
	}
	
	private static EvenlyDiscretizedFunc asDistDiff(EvenlyDiscretizedFunc distFunc) {
		EvenlyDiscretizedFunc ret = distFunc.deepClone();
		
		for (int i=0; i<ret.size(); i++)
			ret.set(i, distFunc.getY(i) - distFunc.getX(i));
		
		return ret;
	}

}