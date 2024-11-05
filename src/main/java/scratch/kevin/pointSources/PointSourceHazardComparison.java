package scratch.kevin.pointSources;

import java.awt.Color;
import java.awt.Font;
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
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.WeightedValue;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
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
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.earthquake.rupForecastImpl.PointToFiniteSource;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_AbstractGridSourceProvider;
import org.opensha.sha.faultSurface.FiniteApproxPointSurface;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PointSourceDistanceCorrection;
import org.opensha.sha.faultSurface.utils.PointSourceDistanceCorrections;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
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
import com.google.common.base.Stopwatch;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Table;

import net.mahdilamb.colormap.Colors;

import com.google.common.collect.ImmutableMap.Builder;

import scratch.kevin.pointSources.InvCDF_RJBCorrPointSurface.RJBCorrInvPDFs;
import scratch.kevin.pointSources.TableBackedDistCorrPointSurface.DistanceCorrTables;

public class PointSourceHazardComparison {
	
	enum PointSourceType {
		TRUE_POINT_SOURCE("True Point Source", false, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample, Random r) {
//				return new PointEqkSource(centerLoc, mfd, 1d, aveRake, aveDip);
				return PointSource.poissonBuilder(centerLoc)
						.truePointSources()
						.forMFDAndFocalMech(mfd, new FocalMechanism(Double.NaN, aveDip, aveRake))
						.duration(1d)
						.build();
			}
		},
		POINT_SOURCE_NSHM("PointSourceNshm", false, PointSourceDistanceCorrections.NSHM_2013) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample, Random r) {
				return new PointSourceNshm(centerLoc, mfd, 1d, mechWtMapForRake(aveRake), getDistCorrs());
			}
		},
		APROX_SUPERSAMPLE_POINT_SOURCE_NSHM("Approx. Supersampled PointSourceNshm", false, PointSourceDistanceCorrections.NONE) { // none is fine here, it's baked in
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
					double aveRake, double aveDip, boolean isSupersample, Random r) {
				ProbEqkSource source = POINT_SOURCE_NSHM.buildSource(centerLoc, mfd, aveRake, aveDip, isSupersample, r);
				if (!isSupersample) {
					// not externally supersampled, approximate it here
					List<ProbEqkRupture> modRups = new ArrayList<>();
					for (ProbEqkRupture rup : source) {
						FiniteApproxPointSurface ptSurf = (FiniteApproxPointSurface)rup.getRuptureSurface();
						rup.setRuptureSurface(new SupersampledApproxPointSurfaceNshm(ptSurf, rup.getMag(), 0.1, SUPERSAMPLE_SPACING));
						modRups.add(rup);
					}
					source = new RupListSource(modRups, source);
				}
				return source;
			}
		},
		POINT_TO_FINITE("PointToFiniteSource", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
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
		},
		SIMPLE_APPROX_FINITE_POINT_NO_CORR("Simple Approx. Finite, No Dist. Corr.", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildPointSource(centerLoc, null, mfd, aveRake, aveDip, r, getDistCorrs());
			}
		},
		SIMPLE_APPROX_FINITE_POINT_NSHMP_CORR("Simple Approx. Finite, NSHMP08 Dist. Corr.", true, PointSourceDistanceCorrections.NSHM_2008) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildPointSource(centerLoc, null, mfd, aveRake, aveDip, r, getDistCorrs());
			}
		},
		SINGLE_QUAD_TRACE_CENTERED("Single Quad Surface, Trace Centered", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.fractionalHypocentralDepth(0d);
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, 1);
			}
		},
		SINGLE_QUAD_SURF_CENTER("Single Quad Surface", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.fractionalHypocentralDepth(0.5);
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, 1);
			}
		},
		CROSSHAIR_QUAD("Crosshair Quad Surface", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.fractionalHypocentralDepth(0.5d);
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, 2);
			}
		},
		QUAD_QUAD("4x Quad Surface", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.fractionalHypocentralDepth(0.5d);
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, isSupersample ? 2 : 4);
			}
		},
		OCT_QUAD("8x Quad Surface", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.fractionalHypocentralDepth(0.5d);
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, isSupersample ? 2 : 8);
			}
		},
		OCT_QUAD_RAND_DAS_DD("8x Quad Surface, Random DAS & DD", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.sampleDASs().sampleHypocentralDepths();
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, isSupersample ? 2 : 8);
			}
		},
		OCT_QUAD_RAND_CELL("8x Quad Surface, Random Cell Locations", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.fractionalHypocentralDepth(0.5d);
				if (!isSupersample)
					builder.sampleFromCell(cellRegion(centerLoc));
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, isSupersample ? 2 : 8);
			}
		},
		OCT_QUAD_RAND_DAS_DD_CELL("8x Quad Surface, Random DAS, DD, & Cell Locations", true, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
				builder.dip(aveDip);
				builder.sampleDASs().sampleHypocentralDepths();
				if (!isSupersample)
					builder.sampleFromCell(cellRegion(centerLoc));
				builder.random(r);
				return buildQuadSource(builder, mfd, aveRake, aveDip, isSupersample ? 2 : 8);
			}
		},
		TABLE_FINITE_APPROX_POINT_SOURCE("Table-Based Finite Approx", false, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildTableBackedSource(centerLoc, mfd, aveRake, aveDip, false, r);
			}
		},
		TABLE_FINITE_APPROX_POINT_SOURCE_SUPERSAMPLE("Table-Based Finite Supersample Approx", false, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildTableBackedSource(centerLoc, mfd, aveRake, aveDip, isSupersample ? false : true, r);
			}
		},
		HAZ_EQUIV_APPROX_POINT_SOURCE("Hazard Equivalent Finite Approx", false, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildHazEquivSource(centerLoc, mfd, aveRake, aveDip, this);
			}
		},
		CDF_BASED_RJB_CORR_APPROX_POINT_SOURCE("RJB Distribution Approx Finite Point Source", false, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildRJBCorrSource(centerLoc, mfd, aveRake, aveDip, OCT_QUAD, 5);
			}
		},
		CDF_BASED_RJB_CORR_APPROX_SUPERSAMPLED_POINT_SOURCE("RJB Distribution Approx Supersampled Finite Point Source", false, PointSourceDistanceCorrections.NONE) {
			@Override
			public ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
					double aveDip, boolean isSupersample, Random r) {
				return buildRJBCorrSource(centerLoc, mfd, aveRake, aveDip, isSupersample ? OCT_QUAD : OCT_QUAD_RAND_CELL, 5);
			}
		};
		
		private String label;
		private boolean stochastic;
		private WeightedList<PointSourceDistanceCorrection> distCorrs;

		private PointSourceType(String label, boolean stochastic, PointSourceDistanceCorrections distCorrType) {
			this(label, stochastic, distCorrType.get());
		}

		private PointSourceType(String label, boolean stochastic, WeightedList<PointSourceDistanceCorrection> distCorrs) {
			this.label = label;
			this.stochastic = stochastic;
			Preconditions.checkNotNull(distCorrs);
			this.distCorrs = distCorrs;
		}
		
		public WeightedList<PointSourceDistanceCorrection> getDistCorrs() {
			return distCorrs;
		}
		
		public abstract ProbEqkSource buildSource(Location centerLoc, IncrementalMagFreqDist mfd,
				double aveRake, double aveDip, boolean isSupersample, Random r);
		
		@Override
		public String toString() {
			return label;
		}
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
	
	private static WC1994_MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();
	
	private static Region cellRegion(Location loc) {
		double half = 0.05;
		return new Region(new Location(loc.lat-half, loc.lon-half), new Location(loc.lat+half, loc.lon+half));
	}
	
	private static ProbEqkSource buildQuadSource(PointSurfaceBuilder builder, IncrementalMagFreqDist mfd,
			double aveRake, double aveDip, int numEach) {
		double dipRad = Math.toRadians(aveDip);
		List<ProbEqkRupture> rups = new ArrayList<>(mfd.size()*numEach);
		for (int i=0; i<mfd.size(); i++) {
			double mag = mfd.getX(i);
			double depth = (float)mag < 6.5f ? 5d : 1d;
			double length = WC94.getMedianLength(mag);
			double aspectWidth = length / 1.5;
			double ddWidth = (14.0 - depth) / Math.sin(dipRad);
			ddWidth = Math.min(aspectWidth, ddWidth);
			double lower = depth + ddWidth * Math.sin(dipRad);

			builder.length(length);
			builder.upperDepth(depth).lowerDepth(lower);
			
			double rateEach = mfd.getY(i)/(double)numEach;
			double probEach = 1-Math.exp(-rateEach);
			for (QuadSurface surf : builder.buildRandQuadSurfaces(numEach))
				rups.add(new ProbEqkRupture(mag, aveRake, probEach, surf, null));
		}
		return new RupListSource(rups, null);
	}
	
	private static ProbEqkSource buildPointSource(Location centerLoc, Region cell, IncrementalMagFreqDist mfd,
			double aveRake, double aveDip, Random r, WeightedList<PointSourceDistanceCorrection> distCorrs) {
		PointSurfaceBuilder builder = new PointSurfaceBuilder(centerLoc);
		builder.dip(aveDip);
		builder.random(r);
		builder.sampleFromCell(cell);
		double dipRad = Math.toRadians(aveDip);
		int expectedSize = mfd.size();
		if (aveDip != 90d)
			expectedSize *= 2;
		if (distCorrs != null)
			expectedSize *= distCorrs.size();
		List<ProbEqkRupture> rups = new ArrayList<>(expectedSize);
		for (int i=0; i<mfd.size(); i++) {
			double mag = mfd.getX(i);
			double depth = (float)mag < 6.5f ? 5d : 1d;
			double length = WC94.getMedianLength(mag);
			double aspectWidth = length / 1.5;
			double ddWidth = (14.0 - depth) / Math.sin(dipRad);
			ddWidth = Math.min(aspectWidth, ddWidth);
			double lower = depth + ddWidth * Math.sin(dipRad);

			builder.length(length);
			builder.upperDepth(depth).lowerDepth(lower);
			
			double totRate = mfd.getY(i);
			WeightedList<FiniteApproxPointSurface> surfs = builder.buildFiniteApproxPointSurfaces(distCorrs);
			
			for (WeightedValue<FiniteApproxPointSurface> surfWeight : surfs) {
				double myRate = totRate*surfWeight.weight;
				double myProb = 1-Math.exp(-myRate);
				ProbEqkRupture rup = new ProbEqkRupture(mag, aveRake, myProb, surfWeight.value, null);
				rups.add(rup);
			}
		}
		return new RupListSource(rups, null);
	}
	
	private static File REF_MAIN_DIR = new File("/home/kevin/markdown/nshm23-misc/point_source_corr");
	
	private static DistanceCorrTables vertCenteredTables = null;
	private static DistanceCorrTables vertSupersampleTables = null;
	private static DistanceCorrTables dippingCenteredTables = null;
	private static DistanceCorrTables dippingSupersampleTables = null;
	private static File refTablesDir = new File(REF_MAIN_DIR, "OCT_QUAD/resources/");
	
	public static ProbEqkSource buildTableBackedSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
			double aveDip, boolean isSupersample, Random r) {
		DistanceCorrTables distCorr = getLoadTables(aveDip, isSupersample);
		
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

			TableBackedDistCorrPointSurface surf = new TableBackedDistCorrPointSurface(
					distCorr, centerLoc, aveDip, depth, lower, length, mag, r);
			double prob = 1-Math.exp(-mfd.getY(i));
			rups.add(new ProbEqkRupture(mag, aveRake, prob, surf, null));
		}
		return new RupListSource(rups, source0);
	}
	
	private static synchronized DistanceCorrTables getLoadTables(double aveDip, boolean isSupersample) {
		try {
			if (aveDip == 90d) {
				if (isSupersample) {
					if (vertSupersampleTables == null)
						vertSupersampleTables = TableBackedDistCorrPointSurface.loadCorrTables(refTablesDir, "STRIKE_SLIP_supersampled");
					return vertSupersampleTables;
				} else {
					if (vertCenteredTables == null)
						vertCenteredTables = TableBackedDistCorrPointSurface.loadCorrTables(refTablesDir, "STRIKE_SLIP_centered");
					return vertCenteredTables;
				}
			} else {
				if (isSupersample) {
					if (dippingSupersampleTables == null)
						dippingSupersampleTables = TableBackedDistCorrPointSurface.loadCorrTables(refTablesDir, "REVERSE_supersampled");
					return dippingSupersampleTables;
				} else {
					if (dippingCenteredTables == null)
						dippingCenteredTables = TableBackedDistCorrPointSurface.loadCorrTables(refTablesDir, "REVERSE_centered");
					return dippingCenteredTables;
				}
			}
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}

	private static Table<PointSourceType, Double, EvenlyDiscrXYZ_DataSet> refHazEquivRJBs = HashBasedTable.create();
	private static Map<PointSourceType, File> refHazEquivDirs = new EnumMap<>(PointSourceType.class);
	static {
		refHazEquivDirs.put(PointSourceType.HAZ_EQUIV_APPROX_POINT_SOURCE,
				new File(REF_MAIN_DIR, "SIMPLE_APPROX_FINITE_POINT_NO_CORR_vs_OCT_QUAD/hazard_equiv_rJBs"));
	}
	
	public static ProbEqkSource buildHazEquivSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
			double aveDip, PointSourceType type) {
		EvenlyDiscrXYZ_DataSet distCorr = getHazEquivRJBTable(type, aveDip);
		
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
			
			if (aveDip < 90d) {
				// one for each footwall boolean
				for (boolean footwall : new boolean[] {false,true}) {
					HazardEquivalentRJBCorrPointSurface surf = new HazardEquivalentRJBCorrPointSurface(
							distCorr, centerLoc, aveDip, depth, lower, length, mag, footwall);
					double prob = 1-Math.exp(-mfd.getY(i)*0.5);
					rups.add(new ProbEqkRupture(mag, aveRake, prob, surf, null));
				}
			} else {
				HazardEquivalentRJBCorrPointSurface surf = new HazardEquivalentRJBCorrPointSurface(
						distCorr, centerLoc, aveDip, depth, lower, length, mag, true);
				double prob = 1-Math.exp(-mfd.getY(i));
				rups.add(new ProbEqkRupture(mag, aveRake, prob, surf, null));
			}
		}
		return new RupListSource(rups, source0);
	}
	
	private static synchronized EvenlyDiscrXYZ_DataSet getHazEquivRJBTable(PointSourceType type, double aveDip) {
		File refDir = refHazEquivDirs.get(type);
		Preconditions.checkNotNull(refDir, "no refDir for type %s", type);
		Preconditions.checkState(refDir.exists(), "Doesn't exist: %s", refDir.getAbsolutePath());
		EvenlyDiscrXYZ_DataSet ret = refHazEquivRJBs.get(type, aveDip);
		if (ret != null)
			return ret;
		try {
			if (aveDip == 90d)
				ret = TableBackedDistCorrPointSurface.loadXYZ_CSV(new File(refDir, "STRIKE_SLIP_equiv_rJBs.csv"));
			else
				ret = TableBackedDistCorrPointSurface.loadXYZ_CSV(new File(refDir, "REVERSE_equiv_rJBs.csv"));
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		refHazEquivRJBs.put(type, aveDip, ret);
		return ret;
	}
	
	private static Table<PointSourceType, Double, RJBCorrInvPDFs> refRJBCorrInvPDFs = HashBasedTable.create();
	
	public static ProbEqkSource buildRJBCorrSource(Location centerLoc, IncrementalMagFreqDist mfd, double aveRake,
			double aveDip, PointSourceType type, int numProbs) {
		RJBCorrInvPDFs corrTable = getRJBCorrTable(type, aveDip);
		
		DiscretizedFunc cdfProbs = InvCDF_RJBCorrPointSurface.getEvenlySpacedProbs(numProbs);
//		DiscretizedFunc cdfProbs = InvCDF_RJBCorrPointSurface.getEvenlySpacedProbs(501);
//		DiscretizedFunc cdfProbs = InvCDF_RJBCorrPointSurface.getSigmaSpacedProbs(5);
		
		boolean[] footwalls = aveDip == 90d ? new boolean[] {true} : new boolean[] {false,true};
		
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
			
			double totRate = mfd.getY(i);
			for (Point2D cdfPt : cdfProbs) {
				double cdfProb = cdfPt.getX();
				double weight =  cdfPt.getY();
				double rateEach = weight*totRate/(double)footwalls.length;
				double probEach = 1-Math.exp(-rateEach);
				for (boolean footwall : footwalls) {
					InvCDF_RJBCorrPointSurface surf = new InvCDF_RJBCorrPointSurface(corrTable, centerLoc, aveDip,
							depth, lower, footwall, length, mag, cdfProb);
					rups.add(new ProbEqkRupture(mag, aveRake, probEach, surf, null));
				}
			}
			
		}
		return new RupListSource(rups, source0);
	}
	
	private static synchronized RJBCorrInvPDFs getRJBCorrTable(PointSourceType type, double aveDip) {
		File refDir = new File(REF_MAIN_DIR, type.name()+"/rjb_dists");
		Preconditions.checkNotNull(refDir, "no refDir for type %s", type);
		Preconditions.checkState(refDir.exists(), "Doesn't exist: %s", refDir.getAbsolutePath());
		RJBCorrInvPDFs ret = refRJBCorrInvPDFs.get(type, aveDip);
		if (ret != null)
			return ret;
		try {
			File csvFile;
			if (aveDip == 90d)
				csvFile = new File(refDir, "STRIKE_SLIP_rJB_cumDists.csv");
			else
				csvFile = new File(refDir, "REVERSE_rJB_cumDists.csv");
			Preconditions.checkState(csvFile.exists());
			CSVFile<String> csv = CSVFile.readFile(csvFile, true);
			ret = new RJBCorrInvPDFs(csv);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		refRJBCorrInvPDFs.put(type, aveDip, ret);
		return ret;
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
			if (source0 == null)
				return rups.get(0).getRuptureSurface().getQuickDistance(site.getLocation());
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
		private PointSourceDistanceCorrection distCorr;

		public CorrOverrideProbEqkSource(ProbEqkSource source, WeightedList<PointSourceDistanceCorrection> distCorrs) {
			this.source = source;
			Preconditions.checkState(distCorrs.size() == 1);
			this.distCorr = distCorrs.getValue(0);
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
			((PointSurface)rup.getRuptureSurface()).setDistanceCorrection(distCorr, rup);
			return rup;
		}
		
	}
	
	static class PointSourceCalcERF extends AbstractERF {
		
		private List<ProbEqkSource> sources;
		private List<Callable<ProbEqkSource>> sourceCalls;
		
		public PointSourceCalcERF(PointSourceType type, Location loc, IncrementalMagFreqDist mfd,
				double rake, double dip, int numPerStochastic, Random r) {
			if (type.stochastic) {
				sources = new ArrayList<>(numPerStochastic);
				IncrementalMagFreqDist mfdEach = mfd.deepClone();
				mfdEach.scale(1d/numPerStochastic);
				for (int i=0; i<numPerStochastic; i++)
					sources.add(type.buildSource(loc, mfdEach, rake, dip, false, r));
			} else {
				sources = List.of(type.buildSource(loc, mfd, rake, dip, false, r));
			}
		}
		
		public PointSourceCalcERF(PointSourceType type, GriddedRegion locs, IncrementalMagFreqDist mfd,
				double rake, double dip, int numPerStochastic, Random r) {
			IncrementalMagFreqDist mfdEach = mfd.deepClone();
			sources = new ArrayList<>(type.stochastic ? numPerStochastic*locs.getNodeCount() : locs.getNodeCount());
			mfdEach.scale(1d/locs.getNodeCount());
			for (Location loc : locs.getNodeList()) {
				if (type.stochastic) {
					IncrementalMagFreqDist subMFDEach = mfdEach.deepClone();
					subMFDEach.scale(1d/numPerStochastic);
					for (int i=0; i<numPerStochastic; i++)
						sources.add(type.buildSource(loc, subMFDEach, rake, dip, true, r));
				} else {
					sources.add(type.buildSource(loc, mfdEach, rake, dip, true, r));
				}
			}
		}
		
		public PointSourceCalcERF(PointSourceType type, GridSourceProvider gridProv, double supersampleDisc,
				int numPerStochastic, Random r, Region calcRegion, double maxDist) {
			sourceCalls = new ArrayList<>();
			GriddedRegion gridReg = gridProv.getGriddedRegion();
			for (int i=0; i<gridProv.getNumLocations(); i++) {
				Location centerLoc = gridReg.getLocation(i);
				if (calcRegion.distanceToLocation(centerLoc) > maxDist)
					continue;
				IncrementalMagFreqDist mfd = gridProv.getMFD(i, 5.05d);
				if (mfd == null || mfd.calcSumOfY_Vals() == 0d)
					continue;
				List<Location> locs;
				IncrementalMagFreqDist fullMFDEach;
				if (supersampleDisc > 0d) {
					double halfSpacing = 0.5*gridReg.getSpacing();
					Region cell = new Region(new Location(centerLoc.lat-halfSpacing, centerLoc.lon-halfSpacing),
							new Location(centerLoc.lat+halfSpacing, centerLoc.lon+halfSpacing));
					GriddedRegion superSampledCell = new GriddedRegion(cell, supersampleDisc,
							new Location(0.5*supersampleDisc, 0.5*supersampleDisc));
					locs = superSampledCell.getNodeList();
					fullMFDEach = mfd.deepClone();
					fullMFDEach.scale(1d/locs.size());
				} else {
					locs = List.of(centerLoc);
					fullMFDEach = mfd;
				}
				for (FocalMech mech : FocalMech.values()) {
					double fract;
					switch (mech) {
					case NORMAL:
						fract = gridProv.getFracNormal(i);
						break;
					case REVERSE:
						fract = gridProv.getFracReverse(i);
						break;
					case STRIKE_SLIP:
						fract = gridProv.getFracStrikeSlip(i);
						break;

					default:
						throw new IllegalStateException();
					}
					if (fract == 0d)
						continue;
					IncrementalMagFreqDist mfdEach = fullMFDEach.deepClone();
					mfdEach.scale(fract);
					for (Location loc : locs) {
						if (type.stochastic && numPerStochastic > 1) {
							IncrementalMagFreqDist subMFDEach = mfdEach.deepClone();
							subMFDEach.scale(1d/numPerStochastic);
							for (int j=0; j<numPerStochastic; j++)
								sourceCalls.add(new Callable<ProbEqkSource>() {
									
									@Override
									public ProbEqkSource call() throws Exception {
										return type.buildSource(loc, subMFDEach, mech.rake(), mech.dip(), locs.size()>1, r);
									}
								});
						} else {
							sourceCalls.add(new Callable<ProbEqkSource>() {
								
								@Override
								public ProbEqkSource call() throws Exception {
									return type.buildSource(loc, mfdEach, mech.rake(), mech.dip(), locs.size()>1, r);
								}
							});
						}
					}
				}
			}
		}
		
		public double calcTotalRate() {
			Preconditions.checkNotNull(sources);
			double tot = 0d;
			for (ProbEqkSource source : sources) {
				for (ProbEqkRupture rup : source) {
					tot += rup.getMeanAnnualRate(1d);
				}
			}
			return tot;
		}

		@Override
		public int getNumSources() {
			return sources != null ? sources.size() : sourceCalls.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			try {
				return sources != null ? sources.get(idx) : sourceCalls.get(idx).call();
			} catch (Exception e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
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
	
	enum RupSurfProps {
		ZTOP("Top Depth (km)") {
			@Override
			public double calcForRupture(ProbEqkRupture rup) {
				return rup.getRuptureSurface().getAveRupTopDepth();
			}
		},
		ZBOT("Bottom Depth (km)") {
			@Override
			public double calcForRupture(ProbEqkRupture rup) {
				RuptureSurface surf = rup.getRuptureSurface();
				return surf.getAveRupTopDepth() + surf.getAveWidth() * Math.sin(Math.toRadians(surf.getAveDip()));
			}
		},
		WIDTH("Width (km)") {
			@Override
			public double calcForRupture(ProbEqkRupture rup) {
				return rup.getRuptureSurface().getAveWidth();
			}
		},
		LENGTH("Length (km)") {
			@Override
			public double calcForRupture(ProbEqkRupture rup) {
				return rup.getRuptureSurface().getAveLength();
			}
		},
		DIP("Dip") {
			@Override
			public double calcForRupture(ProbEqkRupture rup) {
				return rup.getRuptureSurface().getAveDip();
			}
		},
		RAKE("Rake") {
			@Override
			public double calcForRupture(ProbEqkRupture rup) {
				return rup.getRuptureSurface().getAveDip();
			}
		};
		
		private String label;

		private RupSurfProps(String label) {
			this.label = label;
		}
		
		public abstract double calcForRupture(ProbEqkRupture rup);
	}
	
	enum RupSurfMapProps {
		FRACT_FOOTWALL("Fraction On Hanging Wall (rX>=0)") {
			@Override
			public double calcForRuptureAndLocation(ProbEqkRupture rup, Location siteLoc) {
				if (rup.getRuptureSurface().getDistanceX(siteLoc) >= 0d)
					return 1d;
				return 0;
			}

			@Override
			public CPT getCPT() throws IOException {
				return GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
			}
		};
		
		private String label;

		private RupSurfMapProps(String label) {
			this.label = label;
		}
		
		public abstract double calcForRuptureAndLocation(ProbEqkRupture rup, Location siteLoc);

		public abstract CPT getCPT() throws IOException;
	}
	
	private static final double SUPERSAMPLE_SPACING = 0.01;

	// TODO: this is just a marker to make it easy to find the top of main
	public static void main(String[] args) throws IOException {
//		PointSourceType mainType = PointSourceType.POINT_SOURCE_13b_NSHMP_CORR;
		PointSourceType mainType = PointSourceType.POINT_SOURCE_NSHM;
//		PointSourceType mainType = PointSourceType.QUAD_QUAD;
//		PointSourceType mainType = PointSourceType.QUAD_QUAD_RAND_CELL;
//		PointSourceType mainType = PointSourceType.OCT_QUAD_RAND_CELL;
//		PointSourceType mainType = PointSourceType.OCT_QUAD;
//		PointSourceType mainType = PointSourceType.APROX_SUPERSAMPLE_POINT_SOURCE_NSHM;
//		PointSourceType mainType = PointSourceType.POINT_TO_FINITE;
//		PointSourceType mainType = PointSourceType.CROSSHAIR_QUAD;
//		PointSourceType mainType = PointSourceType.QUAD_QUAD;
//		PointSourceType compType = null;
//		PointSourceType compType = PointSourceType.TABLE_FINITE_APPROX_POINT_SOURCE;
//		PointSourceType compType = PointSourceType.HAZ_EQUIV_APPROX_POINT_SOURCE;
//		PointSourceType compType = PointSourceType.CDF_BASED_RJB_CORR_APPROX_POINT_SOURCE;
//		PointSourceType compType = PointSourceType.CDF_BASED_RJB_CORR_APPROX_SUPERSAMPLED_POINT_SOURCE;
//		PointSourceType compType = PointSourceType.SIMPLE_APPROX_FINITE_POINT_NSHMP_CORR;
//		PointSourceType compType = PointSourceType.SIMPLE_APPROX_FINITE_POINT_NO_CORR;
//		PointSourceType compType = PointSourceType.TABLE_FINITE_APPROX_POINT_SOURCE_SUPERSAMPLE;
		PointSourceType compType = PointSourceType.OCT_QUAD;
//		PointSourceType compType = PointSourceType.OCT_QUAD_RAND_CELL;
//		PointSourceType compType = PointSourceType.OCT_QUAD_RAND_DAS_DD;
//		PointSourceType compType = PointSourceType.POINT_SOURCE_NSHM;
//		PointSourceType compType = PointSourceType.POINT_SOURCE_13b_NSHMP_CORR;
//		PointSourceType compType = PointSourceType.POINT_SOURCE_13b_NO_CORR;
//		PointSourceType compType = PointSourceType.POINT_TO_FINITE;
//		PointSourceType compType = PointSourceType.POINT_TO_SIMPLE_QUAD;
//		PointSourceType compType = PointSourceType.QUAD_QUAD_RAND_CELL;
		
//		int sleepMins = 60*6;
//		System.out.println("Sleeping for "+sleepMins+" mins");
//		try {
//			Thread.sleep(sleepMins*60l*1000l);
//		} catch (InterruptedException e) {}
		
		boolean doSingleCellHazard = true;
		boolean doNSHMModelHazard = true;
		boolean doHighRes = false;
		boolean doSupersample = true;
		boolean forceWriteIntermediate = true;
		boolean writeTables = false;
		boolean strikeSlipOnly = true;
		
		AttenRelRef gmmRef = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
		double[] periods = {0d, 1d};
		double[] rakes = strikeSlipOnly ? new double[] {0d} : new double[] {0d, 90d};
		
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
		if (!doSupersample)
			outputDir = new File(mainOutputDir, outputDir.getName()+"_no_supersample");
		if (strikeSlipOnly)
			outputDir = new File(mainOutputDir, outputDir.getName()+"_ss_only");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		boolean writeIntermediate = forceWriteIntermediate || !new File(outputDir, "index.html").exists();
		
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
		
		EvenlyDiscretizedFunc tableDistFunc = writeTables ?
				new EvenlyDiscretizedFunc(0d, 300d, 301) : new EvenlyDiscretizedFunc(0d, 100d, 101);
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(0d, 100d, 101);
		Range distPlotRange = new Range(0d, 100d);
		List<List<Location>> tableDistFuncLocs = new ArrayList<>();
		int numEachDist = 10;
		double maxActualDist = 0d;
		for (int i=0; i<tableDistFunc.size(); i++) {
			double dist = tableDistFunc.getX(i);
			if (dist == 0d) {
				tableDistFuncLocs.add(List.of(gridCenter));
			} else {
				List<Location> locs = new ArrayList<>(numEachDist);
				double deltaEach = 2d*Math.PI/numEachDist;
				for (int j=0; j<numEachDist; j++)
					locs.add(LocationUtils.location(gridCenter, deltaEach*j, dist));
				for (Location loc : locs) {
					double myDist = LocationUtils.horzDistanceFast(loc, gridCenter);
					maxActualDist = Math.max(maxActualDist, myDist);
				}
				tableDistFuncLocs.add(locs);
			}
		}
		List<List<Location>> distFuncLocs = tableDistFuncLocs.subList(0, distFunc.size());
		System.out.println("Maximum relocated distance was "+(float)maxActualDist+", max intended was "+(float)tableDistFunc.getMaxX());
		
		GriddedRegion offsetCalcGridded = new GriddedRegion(calcRegion, spacing, new Location(halfSpacing, halfSpacing));
		siteGrids.put(GridType.OFFSET, offsetCalcGridded);
		
		if (doHighRes)
			siteGrids.put(GridType.HIGH_RES, new GriddedRegion(calcRegion, spacing*0.2, GriddedRegion.ANCHOR_0_0));
		
		Region cellReg = new Region(new Location(gridCenter.lat-halfSpacing, gridCenter.lon-halfSpacing),
				new Location(gridCenter.lat+halfSpacing, gridCenter.lon+halfSpacing));
		GriddedRegion superSampledCell = new GriddedRegion(cellReg, SUPERSAMPLE_SPACING,
				new Location(0.5*SUPERSAMPLE_SPACING, 0.5*SUPERSAMPLE_SPACING));
		
//		int numPerStochasticCentered = 1000;
//		int numPerStochasticSupersampled = 50;
//		int numPerStochasticCentered = 1000;
//		int numPerStochasticSupersampled = 20;
		int numPerStochasticCentered = 500;
		int numPerStochasticSupersampled = Integer.max(10, numPerStochasticCentered/superSampledCell.getNodeCount());
//		int numPerStochasticCentered = 100;
//		int numPerStochasticSupersampled = 5;
		if (mainType.stochastic || (compType != null && compType.stochastic))
			System.out.println("Will do "+numPerStochasticCentered+" stocastic relizations ("+numPerStochasticSupersampled+" when supersampling)");
		
		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(5.05, 30, 0.1);
		mfd.setAllButTotMoRate(mfd.getMinX(), mfd.getMaxX(), 1d, 1d);
//		GriddedRegion calcRegion = new GriddedR
		
		double[] distPlotMags = { 5.05, 6.05, 7.05, mfd.getMaxX() };
		double[] distPlotDists = { 0d, 10d, 50d, 100d, 200d };
		
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
			if (doSupersample) {
				for (Location loc : superSampledCell.getNodeList()) {
					scatterLocs.add(loc);
					scatterChars.add(superSampledChar);
				}
				legendItems.add("Supersampled Sources");
				legendChars.add(scatterChars.get(scatterChars.size()-1));
			}
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
		Deque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();
		Deque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
//		CPT batlowCategorical = GMT_CPT_Files.CATEGORICAL_BATLOW_UNIFORM.instance();
//		Color mainColor = batlowCategorical.get(0).minColor;
//		Color compColor = batlowCategorical.get(1).minColor;
		
//		CPT tab10 = GMT_CPT_Files.CATEGORICAL_TAB10.instance();
//		Color mainColor = tab10.getColor(0);
//		Color compColor = tab10.getColor(1);
		
		Color mainColor = Colors.tab_blue;
		Color compColor = Colors.tab_orange;
		
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
		
		String[] perLabels = new String[periods.length];
		String[] perPrefixes = new String[periods.length];
		String[] perUnits = new String[periods.length];
		Range[] curveXRanges = new Range[periods.length];
		DiscretizedFunc[] xVals = new DiscretizedFunc[periods.length];
		DiscretizedFunc[] logXVals = new DiscretizedFunc[periods.length];
		CPT[] hazCPTs = new CPT[periods.length];
		for (int p=0; p<periods.length; p++) {
			double period = periods[p];
			if (period == -1d) {
				perLabels[p] = "PGV";
				perPrefixes[p] = "pgv";
				perUnits[p] = "cm/s";
				xVals[p] = new IMT_Info().getDefaultHazardCurve(PGV_Param.NAME);
				curveXRanges[p] = new Range(1e-1, 1e2);
			} else if (period == 0d) {
				perLabels[p] = "PGA";
				perPrefixes[p] = "pga";
				perUnits[p] = "g";
				xVals[p] = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
				curveXRanges[p] = new Range(1e-2, 1e1);
			} else {
				Preconditions.checkState(period > 0d);
				perLabels[p] = oDF.format(period)+"s SA";
				perPrefixes[p] = oDF.format(period)+"s";
				perUnits[p] = "g";
				xVals[p] = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
				curveXRanges[p] = new Range(1e-2, 1e1);
			}
			
			hazCPTs[p] = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(
					Math.log10(curveXRanges[p].getLowerBound()), Math.log10(curveXRanges[p].getUpperBound()));
			logXVals[p] = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : xVals[p])
				logXVals[p].set(Math.log(pt.getX()), 0d);
		}
		
		
		Random rand = new Random(123456789l);
		
		boolean[] falseTrue = {false,true};
		boolean[] falseOnly = {false};
		
		for (double rake : rakes) {
			FocalMech mech = mechForRake(rake);
			String mechLabel = mech.toString();
			String mechPrefix = mech.name();
			
			double dip = rakeToDipMap.get(rake);
			PointSourceCalcERF centerERF = new PointSourceCalcERF(mainType, gridCenter, mfd, rake, dip, numPerStochasticCentered, rand);
			PointSourceCalcERF supersampledERF = doSupersample ?
					new PointSourceCalcERF(mainType, superSampledCell, mfd, rake, dip, numPerStochasticSupersampled, rand) : null;
			
			System.out.println("Centered ERF, total rate: "+centerERF.calcTotalRate());
			if (doSupersample)
				System.out.println("Supersampled ERF, total rate: "+supersampledERF.calcTotalRate());
			
			PointSourceCalcERF compCenterERF = null;
			PointSourceCalcERF compSupersampledERF = null;
			boolean[] compBools;
			if (compType != null) {
				compBools = falseTrue;
				compCenterERF = new PointSourceCalcERF(compType, gridCenter, mfd, rake, dip, numPerStochasticCentered, rand);
				if (doSupersample)
					compSupersampledERF = new PointSourceCalcERF(compType, superSampledCell, mfd, rake, dip, numPerStochasticSupersampled, rand);
				System.out.println("Comparison centered ERF, total rate: "+compCenterERF.calcTotalRate());
				if (doSupersample)
					System.out.println("Comparison supersampled ERF, total rate: "+compSupersampledERF.calcTotalRate());
			} else {
				compBools = falseOnly;
			}
			
			lines.add("## "+mechLabel);
			lines.add(topLink); lines.add("");
			
			lines.add("### "+mechLabel+" Distances");
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("DistanceJB", "DistanceRup");
			
			for (double mag : distPlotMags) {
				Range ratioRange = new Range(0, 2);
				Range diffRange = new Range(-30d, 10d);
				
				EvenlyDiscretizedFunc centeredRJBFunc = distFunc.deepClone();
				EvenlyDiscretizedFunc centeredRRupFunc = distFunc.deepClone();
				System.out.println("Calculating distances for M"+(float)mag+", "+mechLabel);
				calcDistFuncsAtMag(centerERF, distFuncLocs, centeredRJBFunc, centeredRRupFunc, mag);
				EvenlyDiscretizedFunc supersampledRJBFunc = distFunc.deepClone();
				EvenlyDiscretizedFunc supersampledRRupFunc = distFunc.deepClone();
				if (doSupersample)
					calcDistFuncsAtMag(supersampledERF, distFuncLocs, supersampledRJBFunc, supersampledRRupFunc, mag);
				EvenlyDiscretizedFunc compCenteredRJBFunc = compType == null ? null : distFunc.deepClone();
				EvenlyDiscretizedFunc compCenteredRRupFunc = compType == null ? null : distFunc.deepClone();
				EvenlyDiscretizedFunc compSupersampledRJBFunc = compType == null ? null : distFunc.deepClone();
				EvenlyDiscretizedFunc compSupersampledRRupFunc = compType == null ? null : distFunc.deepClone();
				if (compType != null) {
					calcDistFuncsAtMag(compCenterERF, distFuncLocs, compCenteredRJBFunc, compCenteredRRupFunc, mag);
					if (doSupersample)
						calcDistFuncsAtMag(compSupersampledERF, distFuncLocs, compSupersampledRJBFunc, compSupersampledRRupFunc, mag);
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
					
					if (doSupersample) {
						supersampledFunc.setName("Supersampled Source");
						funcs.add(asDistRatio(supersampledFunc));
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
						
						if (compType != null) {
							compSupersampledFunc.setName("Comparison Supersampled Source");
							funcs.add(asDistRatio(compSupersampledFunc));
							chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
						}
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
					
					if (doSupersample) {
						funcs.add(asDistDiff(supersampledFunc));
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
						
						if (compType != null) {
							funcs.add(asDistDiff(compSupersampledFunc));
							chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
						}
					}
					
					PlotSpec diffSpec = new PlotSpec(funcs, chars, title,
							xAxisLabel, rRup ? "DistanceRup Difference" : "DistanceJB Difference");
					diffSpec.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
					
					HeadlessGraphPanel gp = PlotUtils.initHeadless();
					
					gp.drawGraphPanel(List.of(ratioSpec, diffSpec), false, false, List.of(distPlotRange), List.of(ratioRange, diffRange));
					
					gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
					
					String prefix = mechPrefix+"_m"+oDF.format(mag)+"_"+(rRup ? "rRup" : "rJB");
					PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 900, true, true, false);
					
					table.addColumn("![Curves]("+resourcesDir.getName()+"/"+prefix+".png)");
				}
				table.finalizeLine();
			}
			
			// now at fixed distances as a function of mag
			Range magRange = new Range(mfd.getMinX()-0.5*mfd.getDelta(), mfd.getMaxX()+0.5*mfd.getDelta());
			List<PlotSpec> magSpecsRRup = new ArrayList<>(distPlotDists.length);
			List<PlotSpec> magSpecsRJB = new ArrayList<>(distPlotDists.length);
			List<Range> magPlotYRanges = new ArrayList<>(distPlotDists.length);
			for (double dist : distPlotDists) {
				List<Location> distLocs = new ArrayList<>(dist == 0d ? 1 : numEachDist);
				if (dist == 0d) {
					distLocs.add(gridCenter);
				} else {
					for (int i=0; i<numEachDist; i++) {
						double deltaEach = 2d*Math.PI/numEachDist;
						for (int j=0; j<numEachDist; j++)
							distLocs.add(LocationUtils.location(gridCenter, deltaEach*j, dist));
					}
				}
				EvenlyDiscretizedFunc centeredRJBFunc = mfd.deepClone();
				EvenlyDiscretizedFunc centeredRRupFunc = mfd.deepClone();
				System.out.println("Calculating distances for at fixed dist="+(float)dist+", "+mechLabel);
				calcMagFuncsAtDist(centerERF, distLocs, centeredRJBFunc, centeredRRupFunc, dist);
				EvenlyDiscretizedFunc supersampledRJBFunc = mfd.deepClone();
				EvenlyDiscretizedFunc supersampledRRupFunc = mfd.deepClone();
				if (doSupersample)
					calcMagFuncsAtDist(supersampledERF, distLocs, supersampledRJBFunc, supersampledRRupFunc, dist);
				EvenlyDiscretizedFunc compCenteredRJBFunc = compType == null ? null : mfd.deepClone();
				EvenlyDiscretizedFunc compCenteredRRupFunc = compType == null ? null : mfd.deepClone();
				EvenlyDiscretizedFunc compSupersampledRJBFunc = compType == null ? null : mfd.deepClone();
				EvenlyDiscretizedFunc compSupersampledRRupFunc = compType == null ? null : mfd.deepClone();
				if (compType != null) {
					calcMagFuncsAtDist(compCenterERF, distLocs, compCenteredRJBFunc, compCenteredRRupFunc, dist);
					if (doSupersample)
						calcMagFuncsAtDist(compSupersampledERF, distLocs, compSupersampledRJBFunc, compSupersampledRRupFunc, dist);
				}
				
				double minY = Double.POSITIVE_INFINITY;
				double maxY = 0d;
				
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
					
					DefaultXY_DataSet fixedDist = new DefaultXY_DataSet();
					fixedDist.set(magRange.getLowerBound(), dist);
					fixedDist.set(magRange.getUpperBound(), dist);
					fixedDist.setName("Bin Center Distance");
					funcs.add(fixedDist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.GRAY));
					
					centeredFunc.setName("Centered Source");
					funcs.add(centeredFunc);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mainColor));
					
					if (compType != null) {
						compCenteredFunc.setName("Comparison Cenetered Source");
						funcs.add(compCenteredFunc);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
					}
					
					if (doSupersample) {
						supersampledFunc.setName("Supersampled Source");
						funcs.add(supersampledFunc);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
						
						if (compType != null) {
							compSupersampledFunc.setName("Comparison Supersampled Source");
							funcs.add(compSupersampledFunc);
							chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
						}
					}
					
					for (XY_DataSet func : funcs) {
						minY = Math.min(minY, func.getMinY());
						maxY = Math.max(maxY, func.getMaxY());
					}
					
					String title = mechLabel+" Distance"+(rRup ? "Rup" : "JB")+" Comparison";
					
					String xAxisLabel = "Magnitude";
					PlotSpec spec = new PlotSpec(funcs, chars, title,
							xAxisLabel, (rRup ? "DistanceRup" : "DistanceJB")+" at "+oDF.format(dist)+" km");
					
					if (dist == distPlotDists[distPlotDists.length-1])
						spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
					
					if (rRup)
						magSpecsRRup.add(spec);
					else
						magSpecsRJB.add(spec);
				}
				double scalar = Math.max(1d, dist/20d);
				maxY = Math.ceil(1.01*maxY/scalar)*scalar;
				minY = Math.floor(0.99*minY/scalar)*scalar;
				Range range = new Range(Math.max(0, minY), maxY);
				magPlotYRanges.add(range);
				
				XYTextAnnotation ann = new XYTextAnnotation(oDF.format(dist)+" km",
						magRange.getLowerBound() + 0.95*magRange.getLength(),
						range.getLowerBound() + 0.975*range.getLength());
				ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));
				ann.setTextAnchor(TextAnchor.TOP_RIGHT);
				magSpecsRJB.get(magSpecsRJB.size()-1).addPlotAnnotation(ann);
				magSpecsRRup.get(magSpecsRRup.size()-1).addPlotAnnotation(ann);
			}
			
			table.initNewLine();
			for (boolean rRup : falseTrue) {
				List<PlotSpec> specs = rRup ? magSpecsRRup : magSpecsRJB;
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(specs, false, false, List.of(magRange), magPlotYRanges);
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				String prefix = mechPrefix+"_fixed_"+(rRup ? "rRup" : "rJB");
				PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 100+specs.size()*300, true, true, false);
				
				table.addColumn("![Curves]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			if (writeIntermediate)
				writeMarkdown(outputDir, lines, tocIndex);
			
			if (writeTables) {
				System.out.println("Calculating mag/dist tables");
				writeMagDistFuncs(centerERF, tableDistFunc, tableDistFuncLocs, mfd, resourcesDir, mechPrefix+"_centered");
				if (doSupersample)
					writeMagDistFuncs(supersampledERF, tableDistFunc, tableDistFuncLocs, mfd, resourcesDir, mechPrefix+"_supersampled");
			}
			
			lines.add("### "+mechLabel+" Surface Properties");
			lines.add(topLink); lines.add("");
			
			// now rupture properties as a function of magnitude
			List<List<PlotSpec>> magPropSpecs = new ArrayList<>();
			List<List<Range>> magPropYRanges = new ArrayList<>();
			int numPerBundle = 3;
			for (RupSurfProps prop : RupSurfProps.values()) {
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				EvenlyDiscretizedFunc centeredFunc = calcMagProps(centerERF, mfd, prop);
				
				centeredFunc.setName("Centered Source");
				funcs.add(centeredFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mainColor));
				
				if (compType != null) {
					EvenlyDiscretizedFunc compCenteredFunc = calcMagProps(compCenterERF, mfd, prop);
					compCenteredFunc.setName("Comparison Cenetered Source");
					funcs.add(compCenteredFunc);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
				}
				
				double minY = Double.POSITIVE_INFINITY;
				double maxY = Double.NEGATIVE_INFINITY;
				for (XY_DataSet func : funcs) {
					minY = Math.min(minY, func.getMinY());
					maxY = Math.max(maxY, func.getMaxY());
				}
				
				String title = mechLabel+" Rupture Properties";
				
				List<PlotSpec> curSpecBundle;
				List<Range> curRangeBundle;
				if (magPropSpecs.isEmpty() || magPropSpecs.get(magPropSpecs.size()-1).size() == numPerBundle) {
					curSpecBundle = new ArrayList<>();
					magPropSpecs.add(curSpecBundle);
					curRangeBundle = new ArrayList<>();
					magPropYRanges.add(curRangeBundle);
				} else {
					curSpecBundle = magPropSpecs.get(magPropSpecs.size()-1);
					curRangeBundle = magPropYRanges.get(magPropYRanges.size()-1);
				}
				
				String xAxisLabel = "Magnitude";
				PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, prop.label);
				if (curSpecBundle.isEmpty())
					spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				curSpecBundle.add(spec);
				curRangeBundle.add(null);
			}
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (int i=0; i<magPropSpecs.size(); i++) {
				List<PlotSpec> specs = magPropSpecs.get(i);
				List<Range> yRanges = magPropYRanges.get(i);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(specs, false, false, List.of(magRange), yRanges);
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				String prefix = mechPrefix+"_surf_props_"+i;
				PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 100+specs.size()*300, true, true, false);
				
				table.addColumn("![Properties]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
			
			if (writeIntermediate)
				writeMarkdown(outputDir, lines, tocIndex);
			
			for (RupSurfMapProps props : RupSurfMapProps.values()) {
				System.out.println("Doing "+mechLabel+", "+props.label);
				lines.add("### "+mechLabel+" "+props.label);
				lines.add(topLink); lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				if (compType != null)
					table.addLine(mainType.label, compType.label);
				else
					table.initNewLine();
				
				CPT cpt = props.getCPT();
				
				String propPrefix = mechPrefix+"_"+props.name();
				
				for (GridType grid : siteGrids.keySet()) {
					GriddedRegion siteGrid = siteGrids.get(grid);
					
					String gridPrefix = propPrefix+"_"+grid.prefix;
					
					GriddedGeoDataSet centeredXYZ = calcMapProps(centerERF, siteGrid, props);
					
					mapMaker.plotXYZData(centeredXYZ, cpt, props.label);
					addRangeAnns(mapMaker, centeredXYZ, grid == GridType.COLOCATED ? gridCenterIndex : -1, true);
					mapMaker.plot(resourcesDir, gridPrefix+"_centered", grid.label+", Centered Sources");
					mapMaker.clearAnnotations();
					
					if (doSupersample) {
						GriddedGeoDataSet supersampledXYZ = calcMapProps(supersampledERF, siteGrid, props);
						addRangeAnns(mapMaker, supersampledXYZ, grid == GridType.COLOCATED ? gridCenterIndex : -1, true);
						mapMaker.plotXYZData(supersampledXYZ, cpt, props.label);
						mapMaker.plot(resourcesDir, gridPrefix+"_supersampled", grid.label+", Supersampled Sources");
						mapMaker.clearAnnotations();
					}
					
					if (compType != null) {
						GriddedGeoDataSet compCenteredXYZ = calcMapProps(compCenterERF, siteGrid, props);
						
						mapMaker.plotXYZData(compCenteredXYZ, cpt, props.label);
						addRangeAnns(mapMaker, compCenteredXYZ, grid == GridType.COLOCATED ? gridCenterIndex : -1, true);
						mapMaker.plot(resourcesDir, gridPrefix+"_centered_comp", grid.label+", Centered Sources");
						mapMaker.clearAnnotations();
						
						if (doSupersample) {
							GriddedGeoDataSet compSupersampledXYZ = calcMapProps(compSupersampledERF, siteGrid, props);
							addRangeAnns(mapMaker, compSupersampledXYZ, grid == GridType.COLOCATED ? gridCenterIndex : -1, true);
							mapMaker.plotXYZData(compSupersampledXYZ, cpt, props.label);
							mapMaker.plot(resourcesDir, gridPrefix+"_supersampled_comp", grid.label+", Supersampled Sources");
							mapMaker.clearAnnotations();
						}
					}
					
					table.addColumn("![Map]("+resourcesDir.getName()+"/"+gridPrefix+"_centered.png)");
					
					if (compType != null) {
						table.addColumn("![Map]("+resourcesDir.getName()+"/"+gridPrefix+"_centered_comp.png)");
						table.finalizeLine();
						if (doSupersample)
							table.initNewLine();
					}
					
					if (doSupersample) {
						table.addColumn("![Map]("+resourcesDir.getName()+"/"+gridPrefix+"_supersampled.png)");
						
						if (compType != null)
							table.addColumn("![Map]("+resourcesDir.getName()+"/"+gridPrefix+"_supersampled_comp.png)");
						table.finalizeLine();
					} else if (compType == null) {
						table.finalizeLine();
					}
				}
				
//				if (compType == null && !doSupersample)
//					table.finalizeLine();
				lines.addAll(table.build());
				lines.add("");
				
				if (writeIntermediate)
					writeMarkdown(outputDir, lines, tocIndex);
			}
			
			if (!doSingleCellHazard)
				continue;
			
			for (int p=0; p<periods.length; p++) {
				double period = periods[p];
				if (periods.length > 1) {
					lines.add("### "+mechLabel+", "+perLabels[p]);
					lines.add(topLink); lines.add("");
				}
				
				String mechPerPrefix = mechPrefix+"_"+perPrefixes[p];
				
				Map<GridType, GriddedGeoDataSet[]> centeredMaps = new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> supersampledMaps = new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> compCenteredMaps = compType == null ? null : new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> compSupersampledMaps = compType == null ? null : new EnumMap<>(GridType.class);
				
				for (GridType grid : siteGrids.keySet()) {
					GriddedRegion siteGrid = siteGrids.get(grid);
					
					for (boolean comp : compBools) {
						PointSourceCalcERF myCenterERF = comp ? compCenterERF : centerERF;
						System.out.println("Calculating "+grid.label+", centered, "+mechLabel+", "+perLabels[p]+", comp="+comp);
						List<DiscretizedFunc> centeredCurves = calcCurves(siteGrid.getNodeList(), myCenterERF, gmmRef,
								gmmDeque, calcDeque, period, xVals[p], logXVals[p], exec);
						
						List<DiscretizedFunc> supersampledCurves = null;
						if (doSupersample) {
							PointSourceCalcERF mySupersampledERF = comp ? compSupersampledERF : supersampledERF;
							System.out.println("Calculating "+grid.label+", supersampled, "+mechLabel+", "+perLabels[p]+", comp="+comp);
							supersampledCurves = calcCurves(siteGrid.getNodeList(), mySupersampledERF, gmmRef,
									gmmDeque, calcDeque, period, xVals[p], logXVals[p], exec);
						}
						
						GriddedGeoDataSet[] myCenteredMaps = new GriddedGeoDataSet[rps.length];
						GriddedGeoDataSet[] mySupersampledMaps = doSupersample ? new GriddedGeoDataSet[rps.length] : null;
						
						for (int r=0; r<rps.length; r++) {
							myCenteredMaps[r] = curvesToMap(siteGrid, centeredCurves, rps[r]);
							if (doSupersample)
								mySupersampledMaps[r] = curvesToMap(siteGrid, supersampledCurves, rps[r]);
						}
						
						if (comp) {
							compCenteredMaps.put(grid, myCenteredMaps);
							if (doSupersample)
								compSupersampledMaps.put(grid, mySupersampledMaps);
						} else {
							centeredMaps.put(grid, myCenteredMaps);
							if (doSupersample)
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
					
					System.out.println("Calculating dist func, centered, "+mechLabel+", "+perLabels[p]+", comp="+comp);
					List<DiscretizedFunc> calcDistFuncCenteredCurves = calcAverageCurves(distFuncLocs, myCenterERF, gmmRef,
							gmmDeque, calcDeque, period, xVals[p], logXVals[p], exec);
					
					List<DiscretizedFunc> calcDistFuncSupersampledCurves = null;
					if (doSupersample) {
						System.out.println("Calculating dist func, supersampled, "+mechLabel+", "+perLabels[p]+", comp="+comp);
						calcDistFuncSupersampledCurves = calcAverageCurves(distFuncLocs, mySupersampledERF, gmmRef,
								gmmDeque, calcDeque, period, xVals[p], logXVals[p], exec);
					}
					
					if (comp) {
						compColocatedCurve = calcDistFuncCenteredCurves.get(0);
						compDistFuncCenteredCurves = calcDistFuncCenteredCurves;
						if (doSupersample) {
							compDistFuncSupersampledCurves = calcDistFuncSupersampledCurves;
							compColocatedSupersampledCurve = calcDistFuncSupersampledCurves.get(0);
						}
					} else {
						colocatedCurve = calcDistFuncCenteredCurves.get(0);
						distFuncCenteredCurves = calcDistFuncCenteredCurves;
						if (doSupersample) {
							distFuncSupersampledCurves = calcDistFuncSupersampledCurves;
							colocatedSupersampledCurve = calcDistFuncSupersampledCurves.get(0);
						}
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
				
				if (doSupersample) {
					colocatedSupersampledCurve.setName("Supersampled Source");
					funcs.add(colocatedSupersampledCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, mainColor));
					
					if (compType != null) {
						compColocatedSupersampledCurve.setName("Comparison Supersampled Source");
						funcs.add(compColocatedSupersampledCurve);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, compColor));
					}
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, "Zero-Distance Hazard Curves",
						perLabels[p]+" ("+perUnits[p]+")", "Annual Probability of Exceedance");
				spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(spec, true, true, curveXRanges[p], curveYRange);
				
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
					
					if (doSupersample) {
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
					}
					
					spec = new PlotSpec(funcs, chars, rps[r].label+" vs Distance",
							"Distance (km)", rps[r].label+", "+perLabels[p]+" ("+perUnits[p]+")");
					spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					gp.drawGraphPanel(spec, false, true, distPlotRange, logBufferedYRange(funcs));
					
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
							lines.add("#### "+mechLabel+", "+perLabels[p]+", "+rps[r].label);
						else
							lines.add("### "+mechLabel+", "+rps[r].label);
						lines.add(topLink); lines.add("");
					}
					String mapPrefix = mechPerPrefix+"_"+rps[r].name();
					
					String hazLabel = rps[r].label+", "+perLabels[p]+" ("+perUnits[p]+")";
					String hazChangeLabel = rps[r].label+", "+perLabels[p]+", % Change";
					
					for (GridType grid : siteGrids.keySet()) {
						
						lines.add("__"+grid.label+" Grid__");
						lines.add("");
						
						table = MarkdownUtils.tableBuilder();
						if (compType != null)
							table.addLine(mainType.label, compType.label, "% Change");
						else
							table.initNewLine();
						
						boolean[] superBools = doSupersample ? falseTrue : falseOnly;
						for (boolean supersampled : superBools) {
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
							
							mapMaker.plotXYZData(log10(mainMap), hazCPTs[p], hazLabel);
							mapMaker.plot(resourcesDir, mainPrefix, title);
							
							if (compType != null)
								table.initNewLine();
							table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+".png)");
							if (compType != null) {
								mapMaker.plotXYZData(log10(compMap), hazCPTs[p], hazLabel);
								mapMaker.plot(resourcesDir, mainPrefix+"_comp", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+"_comp.png)");
								
								GriddedGeoDataSet pDiff = mapPDiff(mainMap, compMap);
								mapMaker.plotXYZData(pDiff, pDiffCPT, hazChangeLabel);
								addPDiffAnns(mapMaker, pDiff, grid == GridType.COLOCATED ? gridCenterIndex : -1, false);
								mapMaker.plot(resourcesDir, mainPrefix+"_comp_pDiff", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+"_comp_pDiff.png)");
								table.finalizeLine();
								
								mapMaker.clearAnnotations();
							}
						}
						
						// now add centered vs supersampled
						if (doSupersample) {
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
							addPDiffAnns(mapMaker, pDiff, grid == GridType.COLOCATED ? gridCenterIndex : -1, false);
							mapMaker.plot(resourcesDir, diffPrefix, title);
							table.addColumn("![map]("+resourcesDir.getName()+"/"+diffPrefix+".png)");
							mapMaker.clearAnnotations();
							
							if (compType != null) {
								mapMaker.plotXYZData(compPDiff, pDiffCPT, hazChangeLabel);
								addPDiffAnns(mapMaker, compPDiff, grid == GridType.COLOCATED ? gridCenterIndex : -1, false);
								mapMaker.plot(resourcesDir, diffPrefix+"_comp", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+diffPrefix+"_comp.png)");
								mapMaker.clearAnnotations();
								
								table.addColumn("");
							}
							table.finalizeLine();
						} else if (compType == null) {
							table.finalizeLine();
						}
						
						lines.addAll(table.build());
						lines.add("");
					}
				}
				
				if (writeIntermediate)
					writeMarkdown(outputDir, lines, tocIndex);
			}
		}
		
		if (doNSHMModelHazard) {
			String modelLabel = "NSHM23";
			System.out.println("Doing "+modelLabel);
			GridSourceProvider gridProv = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
					+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip")).getGridSourceProvider();
			
			if (strikeSlipOnly) {
				GriddedRegion region = gridProv.getGriddedRegion();
				int nodeCount = region.getNodeCount();
				Builder<Integer, IncrementalMagFreqDist> subSeisBuilder = ImmutableMap.builder();
				Builder<Integer, IncrementalMagFreqDist> unassociatedBuilder = ImmutableMap.builder();
				double[] fracStrikeSlip = new double[nodeCount];
				double[] fracNormal = new double[nodeCount];
				double[] fracReverse = new double[nodeCount];
				for (int i=0; i<nodeCount; i++) {
					IncrementalMagFreqDist subSeis = gridProv.getMFD_SubSeisOnFault(i);
					if (subSeis != null)
						subSeisBuilder.put(i, subSeis);
					IncrementalMagFreqDist unassociated = gridProv.getMFD_Unassociated(i);
					if (unassociated != null)
						unassociatedBuilder.put(i, unassociated);
					fracStrikeSlip[i] = 1d;
					fracNormal[i] = 0d;
					fracReverse[i] = 0d;
				}
				ImmutableMap<Integer, IncrementalMagFreqDist> nodeSubSeisMFDs = subSeisBuilder.build();
				ImmutableMap<Integer, IncrementalMagFreqDist> nodeUnassociatedMFDs = unassociatedBuilder.build();
				gridProv = new NSHM23_AbstractGridSourceProvider.Precomputed(region, nodeSubSeisMFDs, nodeUnassociatedMFDs,
						fracStrikeSlip, fracNormal, fracReverse, ((MFDGridSourceProvider)gridProv).getTectonicRegionTypeArray());
			}
			
			lines.add("## "+modelLabel+" Gridded Seismicity Hazard");
			lines.add(topLink); lines.add("");
			
//			calcRegion = NSHM23_SingleStates.UT.loadRegion();
			calcRegion = new Region(new Location(37, -120), new Location(39, -118)); // SS and Reverse, CA-NV Border
//			calcRegion = new Region(new Location(42, -109), new Location(44, -107)); // SS-only region, Wyoming
			
			colocatedGridded = new GriddedRegion(calcRegion, spacing, GriddedRegion.ANCHOR_0_0);
			
			siteGrids = new EnumMap<>(GridType.class);
			siteGrids.put(GridType.COLOCATED, colocatedGridded);
			
			offsetCalcGridded = new GriddedRegion(calcRegion, spacing, new Location(halfSpacing, halfSpacing));
			siteGrids.put(GridType.OFFSET, offsetCalcGridded);
			
//			if (doHighRes)
//				siteGrids.put(GridType.HIGH_RES, new GriddedRegion(calcRegion, spacing*0.2, GriddedRegion.ANCHOR_0_0));
			
			String compName = compType == null ? null : compType.label;
			if (compType == null && mainType.stochastic) {
				// compare with another realization
				compType = mainType;
				compName = "Second Random Realization";
			}
			
			numPerStochasticCentered = 1;
			numPerStochasticSupersampled = 1;
			double maxDistOutside = 100d;
			double modelSuperSampleSpacing = Math.max(0.02, SUPERSAMPLE_SPACING);
//			double modelSuperSampleSpacing = SUPERSAMPLE_SPACING;
			
			System.out.println("Building "+modelLabel+" ERF, "+mainType.label);
			PointSourceCalcERF centerERF = new PointSourceCalcERF(mainType, gridProv, 0d,
					numPerStochasticCentered, rand, calcRegion, maxDistOutside);
			PointSourceCalcERF supersampledERF = null;
			if (doSupersample) {
				System.out.println("Building "+modelLabel+" ERF, "+mainType.label+", Supersampled");
				supersampledERF = new PointSourceCalcERF(mainType, gridProv, modelSuperSampleSpacing,
						numPerStochasticCentered, rand, calcRegion, maxDistOutside);
			}
			
			PointSourceCalcERF compCenterERF = null;
			PointSourceCalcERF compSupersampledERF = null;
			boolean[] compBools;
			if (compType != null) {
				compBools = falseTrue;
				System.out.println("Building "+modelLabel+" ERF, "+compName);
				compCenterERF = new PointSourceCalcERF(compType, gridProv, 0d,
						numPerStochasticCentered, rand, calcRegion, maxDistOutside);
				if (doSupersample) {
					System.out.println("Building "+modelLabel+" ERF, "+compName+", Supersampled");
					compSupersampledERF = new PointSourceCalcERF(compType, gridProv, modelSuperSampleSpacing,
							numPerStochasticCentered, rand, calcRegion, maxDistOutside);
				}
			} else {
				compBools = falseOnly;
			}
			
			System.out.println("Done building ERFs");
			
			for (int p=0; p<periods.length; p++) {
				double period = periods[p];
				if (periods.length > 1) {
					lines.add("### NSHM23 "+perLabels[p]);
					lines.add(topLink); lines.add("");
				}
				
				String modelPerPrefix = "model_"+perPrefixes[p];
				
				Map<GridType, GriddedGeoDataSet[]> centeredMaps = new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> supersampledMaps = new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> compCenteredMaps = compType == null ? null : new EnumMap<>(GridType.class);
				Map<GridType, GriddedGeoDataSet[]> compSupersampledMaps = compType == null ? null : new EnumMap<>(GridType.class);
				
				table = MarkdownUtils.tableBuilder();
				table.initNewLine().addColumn("");
				for (GridType type : siteGrids.keySet())
					table.addColumn(type.label);
				table.finalizeLine();
				
				for (boolean comp : compBools) {
					List<String> centeredTimes = new ArrayList<>();
					List<String> supersampledTimes = new ArrayList<>();
					for (GridType grid : siteGrids.keySet()) {
						GriddedRegion siteGrid = siteGrids.get(grid);
						PointSourceCalcERF myCenterERF = comp ? compCenterERF : centerERF;
						PointSourceCalcERF mySupersampledERF = comp ? compSupersampledERF : supersampledERF;
						System.out.println("Calculating "+grid.label+", centered, "+modelLabel+", "+perLabels[p]+", comp="+comp);
						Stopwatch watch = Stopwatch.createStarted();
						List<DiscretizedFunc> centeredCurves = calcCurves(siteGrid.getNodeList(), myCenterERF, gmmRef,
								gmmDeque, calcDeque, period, xVals[p], logXVals[p], exec);
						watch.stop();
						double time = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
						System.out.println("Took "+timeStr(time));
						centeredTimes.add(timeStr(time)+" ("+timeStr(time/siteGrid.getNodeCount())+"/curve)");
						
						List<DiscretizedFunc> supersampledCurves = null;
						if (doSupersample) {
							System.out.println("Calculating "+grid.label+", supersampled, "+modelLabel+", "+perLabels[p]+", comp="+comp);
							watch = Stopwatch.createStarted();
							supersampledCurves = calcCurves(siteGrid.getNodeList(), mySupersampledERF, gmmRef,
									gmmDeque, calcDeque, period, xVals[p], logXVals[p], exec);
							watch.stop();
							time = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
							System.out.println("Took "+timeStr(time));
							supersampledTimes.add(timeStr(time)+" ("+timeStr(time/siteGrid.getNodeCount())+"/curve)");
						}
						
						GriddedGeoDataSet[] myCenteredMaps = new GriddedGeoDataSet[rps.length];
						GriddedGeoDataSet[] mySupersampledMaps = doSupersample ? new GriddedGeoDataSet[rps.length] : null;
						
						for (int r=0; r<rps.length; r++) {
							myCenteredMaps[r] = curvesToMap(siteGrid, centeredCurves, rps[r]);
							if (doSupersample)
								mySupersampledMaps[r] = curvesToMap(siteGrid, supersampledCurves, rps[r]);
						}
						
						if (comp) {
							compCenteredMaps.put(grid, myCenteredMaps);
							if (doSupersample)
								compSupersampledMaps.put(grid, mySupersampledMaps);
						} else {
							centeredMaps.put(grid, myCenteredMaps);
							if (doSupersample)
								supersampledMaps.put(grid, mySupersampledMaps);
						}
					}
					
					table.initNewLine();
					if (comp) {
						table.addColumn(compName+", Centered");
					} else {
						table.addColumn(mainType.label+", Centered");
					}
					table.addColumns(centeredTimes);
					table.finalizeLine();
					
					if (doSupersample) {
						table.initNewLine();
						if (comp) {
							table.addColumn(compName+", Supersampled");
						} else {
							table.addColumn(mainType.label+", Supersampled");
						}
						table.addColumns(supersampledTimes);
						table.finalizeLine();
					}
				}
				
				lines.add("__Calculation Times:__");
				lines.add("");
				lines.addAll(table.build());
				lines.add("");
				
				mapMaker.setRegion(calcRegion);
				
				for (int r=0; r<rps.length; r++) {
					if (rps.length > 1) {
						if (periods.length > 1)
							lines.add("#### "+modelLabel+", "+perLabels[p]+", "+rps[r].label);
						else
							lines.add("### "+modelLabel+", "+rps[r].label);
						lines.add(topLink); lines.add("");
					}
					String mapPrefix = modelPerPrefix+"_"+rps[r].name();
					
					String hazLabel = rps[r].label+", "+perLabels[p]+" ("+perUnits[p]+")";
					String hazChangeLabel = rps[r].label+", "+perLabels[p]+", % Change";
					
					for (GridType grid : siteGrids.keySet()) {
						
						lines.add("__"+grid.label+" Grid__");
						lines.add("");
						
						table = MarkdownUtils.tableBuilder();
						if (compType != null)
							table.addLine(mainType.label, compName, "% Change");
						else
							table.initNewLine();
						
						boolean[] superBools = doSupersample ? falseTrue : falseOnly;
						for (boolean supersampled : superBools) {
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
							
							mapMaker.plotXYZData(log10(mainMap), hazCPTs[p], hazLabel);
							mapMaker.plot(resourcesDir, mainPrefix, title);
							
							if (compType != null)
								table.initNewLine();
							table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+".png)");
							if (compType != null) {
								mapMaker.plotXYZData(log10(compMap), hazCPTs[p], hazLabel);
								mapMaker.plot(resourcesDir, mainPrefix+"_comp", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+"_comp.png)");
								
								GriddedGeoDataSet pDiff = mapPDiff(mainMap, compMap);
								mapMaker.plotXYZData(pDiff, pDiffCPT, hazChangeLabel);
								addPDiffAnns(mapMaker, pDiff, -1, true);
								mapMaker.plot(resourcesDir, mainPrefix+"_comp_pDiff", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+mainPrefix+"_comp_pDiff.png)");
								table.finalizeLine();
								
								mapMaker.clearAnnotations();
							}
						}
						
						// now add centered vs supersampled
						if (doSupersample) {
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
							addPDiffAnns(mapMaker, pDiff, -1, true);
							mapMaker.plot(resourcesDir, diffPrefix, title);
							table.addColumn("![map]("+resourcesDir.getName()+"/"+diffPrefix+".png)");
							mapMaker.clearAnnotations();
							
							if (compType != null) {
								mapMaker.plotXYZData(compPDiff, pDiffCPT, hazChangeLabel);
								addPDiffAnns(mapMaker, compPDiff, -1, true);
								mapMaker.plot(resourcesDir, diffPrefix+"_comp", title);
								table.addColumn("![map]("+resourcesDir.getName()+"/"+diffPrefix+"_comp.png)");
								mapMaker.clearAnnotations();
								
								table.addColumn("");
							}
							table.finalizeLine();
						} else if (compType == null) {
							table.finalizeLine();
						}
						
						lines.addAll(table.build());
						lines.add("");
					}
				}
				
				if (writeIntermediate)
					writeMarkdown(outputDir, lines, tocIndex);
			}
		}
		
		exec.shutdown();
		
		writeMarkdown(outputDir, lines, tocIndex);
	}
	
	private static String timeStr(double secs) {
		if (secs < 60d)
			return oDF.format(secs)+" s";
		double mins = secs/60d;
		if (mins < 60d)
			return oDF.format(mins)+" m";
		double hours = mins / 60d;
		return oDF.format(hours)+" h";
	}
	
	private static void writeMarkdown(File outputDir, List<String> lines, int tocIndex) throws IOException {
		lines = new ArrayList<>(lines);
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		System.out.println("Writing markdown and HTML to "+outputDir.getAbsolutePath());
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}

	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	private static final Font ANN_FONT = new Font(Font.SANS_SERIF, Font.BOLD, 40);
	private static final DecimalFormat TWO_DIGITS = new DecimalFormat("0.00");
	private static void addPDiffAnns(GeographicMapMaker mapMaker, GriddedGeoDataSet pDiff,
			int centerIndex, boolean doAvg) {
		addRangeAnns(mapMaker, pDiff, centerIndex, doAvg, "%");
	}
	
	private static void addRangeAnns(GeographicMapMaker mapMaker, GriddedGeoDataSet xyz,
			int centerIndex, boolean doAvg) {
		addRangeAnns(mapMaker, xyz, centerIndex, doAvg, "");
	}

	private static void addRangeAnns(GeographicMapMaker mapMaker, GriddedGeoDataSet xyz,
			int centerIndex, boolean doAvg, String suffix) {
		Range yRange = mapMaker.getYRange();
		Range xRange = mapMaker.getXRange();
		double x = xRange.getLowerBound() + 0.975*xRange.getLength();
		double y = yRange.getLowerBound() + 0.975*yRange.getLength();
		XYTextAnnotation ann = new XYTextAnnotation("["+TWO_DIGITS.format(xyz.getMinZ())
			+suffix+", "+TWO_DIGITS.format(xyz.getMaxZ())+suffix+"]", x, y);
		ann.setFont(ANN_FONT);
		ann.setTextAnchor(TextAnchor.TOP_RIGHT);
		
		mapMaker.addAnnotation(ann);
		
		if (centerIndex >= 0) {
			y -= 0.075*yRange.getLength();
			ann = new XYTextAnnotation("Center: "+TWO_DIGITS.format(xyz.get(centerIndex))+suffix, x, y);
			ann.setFont(ANN_FONT);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			
			mapMaker.addAnnotation(ann);
		}
		
		if (doAvg) {
			y -= 0.075*yRange.getLength();
			double avg = 0d;
			for (int i=0; i<xyz.size(); i++)
				avg += xyz.get(i);
			avg /= (double)xyz.size();
			ann = new XYTextAnnotation("Mean: "+TWO_DIGITS.format(avg)+suffix, x, y);
			ann.setFont(ANN_FONT);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			
			mapMaker.addAnnotation(ann);
		}
	}
	
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static List<DiscretizedFunc> calcCurves(List<Location> siteLocs, PointSourceCalcERF erf,
			AttenRelRef gmmRef, Deque<ScalarIMR> gmmDeque, Deque<HazardCurveCalculator> calcDeque,
			double period, DiscretizedFunc xVals,
			DiscretizedFunc logXVals, ExecutorService exec) {
		List<Future<DiscretizedFunc>> futures = new ArrayList<>(siteLocs.size());
		
		for (Location siteLoc : siteLocs)
			futures.add(exec.submit(new CurveCalcCall(siteLoc, erf, gmmRef, gmmDeque, calcDeque, period, xVals, logXVals)));
		
		System.out.println("Calculating "+siteLocs.size()+" curves");;
		List<DiscretizedFunc> ret = new ArrayList<>();
		for (Future<DiscretizedFunc> future : futures) {
			try {
				ret.add(future.get());
				System.out.print(".");
				if (ret.size() % 100 == 0)
					System.out.println(" "+ret.size()+" ("+pDF.format((double)ret.size()/(double)siteLocs.size())+")");
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		if (ret.size() % 100 != 0)
			System.out.println(" "+ret.size()+" ("+pDF.format(1d)+")");
		return ret;
	}
	
	private static List<DiscretizedFunc> calcAverageCurves(List<List<Location>> siteLocs, PointSourceCalcERF erf,
			AttenRelRef gmmRef, Deque<ScalarIMR> gmmDeque, Deque<HazardCurveCalculator> calcDeque,
			double period, DiscretizedFunc xVals, DiscretizedFunc logXVals, ExecutorService exec) {
		List<List<Future<DiscretizedFunc>>> futures = new ArrayList<>(siteLocs.size());
		
		for (List<Location> siteLocBundle : siteLocs) {
			List<Future<DiscretizedFunc>> bundle = new ArrayList<>(siteLocBundle.size());
			futures.add(bundle);
			for (Location siteLoc : siteLocBundle)
				bundle.add(exec.submit(new CurveCalcCall(siteLoc, erf, gmmRef, gmmDeque, calcDeque, period, xVals, logXVals)));
		}
		
		List<DiscretizedFunc> ret = new ArrayList<>();
		for (List<Future<DiscretizedFunc>> bundle : futures) {
			List<DiscretizedFunc> curves = new ArrayList<>(bundle.size());
			for (Future<DiscretizedFunc> future : bundle) {
				try {
					curves.add(future.get());
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			if (curves.size() == 1) {
				ret.add(curves.get(0));
			} else {
				double scalar = 1d/curves.size();
				DiscretizedFunc curve = xVals.deepClone();
				for (int x=0; x<curve.size(); x++) {
					double sumY = 0d;
					for (DiscretizedFunc subCurve : curves)
						sumY += subCurve.getY(x);
					curve.set(x, sumY*scalar);
				}
				ret.add(curve);
			}
			System.out.print(".");
			if (ret.size() % 100 == 0)
				System.out.println(" "+ret.size()+" ("+pDF.format((double)ret.size()/(double)siteLocs.size())+")");
		}
		if (ret.size() % 100 != 0)
			System.out.println(" "+ret.size()+" ("+pDF.format(1d)+")");
		return ret;
	}
	
	private static class CurveCalcCall implements Callable<DiscretizedFunc> {
		
		private Location loc;
		private PointSourceCalcERF erf;
		private AttenRelRef gmmRef;
		private Deque<ScalarIMR> gmmDeque;
		private Deque<HazardCurveCalculator> calcDeque;
		private double period;
		private DiscretizedFunc xVals;
		private DiscretizedFunc logXVals;

		public CurveCalcCall(Location loc, PointSourceCalcERF erf,
				AttenRelRef gmmRef, Deque<ScalarIMR> gmmDeque, Deque<HazardCurveCalculator> calcDeque,
				double period, DiscretizedFunc xVals,
				DiscretizedFunc logXVals) {
					this.loc = loc;
					this.erf = erf;
					this.gmmRef = gmmRef;
					this.gmmDeque = gmmDeque;
					this.calcDeque = calcDeque;
					this.period = period;
					this.xVals = xVals;
					this.logXVals = logXVals;
			
		}

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
			
			Site site = new Site(loc);
			site.addParameterList(gmm.getSiteParams());
			
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
	
	private static List<ProbEqkRupture> getRupsAtMag(PointSourceCalcERF erf, double mag) {
		List<ProbEqkRupture> rupsAtMag = new ArrayList<>();
		for (ProbEqkSource source : erf)
			for (ProbEqkRupture rup : source)
				if (Precision.equals(rup.getMag(), mag))
					rupsAtMag.add(rup);
		Preconditions.checkState(!rupsAtMag.isEmpty());
		return rupsAtMag;
	}
	
	private static void calcDistFuncsAtMag(PointSourceCalcERF erf, List<List<Location>> distLocs,
			EvenlyDiscretizedFunc rJB, EvenlyDiscretizedFunc rRup, double mag) {
		List<ProbEqkRupture> rupsAtMag = getRupsAtMag(erf, mag);
		IntStream.range(0, distLocs.size()).parallel().forEach(i -> {
			double sumWeight = 0d;
			double sumWeightedRJB = 0d;
			double sumWeightedRRup = 0d;
			for (Location loc : distLocs.get(i)) {
				for (ProbEqkRupture rup : rupsAtMag) {
					double weight = rup.getMeanAnnualRate(1);
					RuptureSurface surf = rup.getRuptureSurface();
					sumWeight += weight;
					sumWeightedRJB += weight*surf.getDistanceJB(loc);
					sumWeightedRRup += weight*surf.getDistanceRup(loc);
				}
			}
			rJB.set(i, sumWeightedRJB/sumWeight);
			rRup.set(i, sumWeightedRRup/sumWeight);
		});
	}
	
	private static void calcMagFuncsAtDist(PointSourceCalcERF erf, List<Location> distLocs,
			EvenlyDiscretizedFunc rJB, EvenlyDiscretizedFunc rRup, double dist) {
		IntStream.range(0, rJB.size()).parallel().forEach(i -> {
			List<ProbEqkRupture> rupsAtMag = getRupsAtMag(erf, rJB.getX(i));
			double sumWeight = 0d;
			double sumWeightedRJB = 0d;
			double sumWeightedRRup = 0d;
			for (Location loc : distLocs) {
				for (ProbEqkRupture rup : rupsAtMag) {
					double weight = rup.getMeanAnnualRate(1);
					RuptureSurface surf = rup.getRuptureSurface();
					sumWeight += weight;
					sumWeightedRJB += weight*surf.getDistanceJB(loc);
					sumWeightedRRup += weight*surf.getDistanceRup(loc);
				}
			}
			rJB.set(i, sumWeightedRJB/sumWeight);
			rRup.set(i, sumWeightedRRup/sumWeight);
		});
	}
	
	private static EvenlyDiscretizedFunc calcMagProps(PointSourceCalcERF erf, IncrementalMagFreqDist mfd,
			RupSurfProps prop) {
		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(mfd.getMinX(), mfd.size(), mfd.getDelta());
		IntStream.range(0, mfd.size()).parallel().forEach(i -> {
			List<ProbEqkRupture> rupsAtMag = getRupsAtMag(erf, mfd.getX(i));
			double sumWeight = 0d;
			double sumWeightedProps = 0d;
			for (ProbEqkRupture rup : rupsAtMag) {
				double weight = rup.getMeanAnnualRate(1);
				sumWeight += weight;
				sumWeightedProps += weight*prop.calcForRupture(rup);
			}
			ret.set(i, sumWeightedProps/sumWeight);
		});
		return ret;
	}
	
	private static GriddedGeoDataSet calcMapProps(PointSourceCalcERF erf, GriddedRegion gridReg,
			RupSurfMapProps prop) {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg);
		IntStream.range(0, gridReg.getNodeCount()).parallel().forEach(i -> {
			Location loc = gridReg.getLocation(i);
			double sumWeight = 0d;
			double sumWeightedProps = 0d;
			for (ProbEqkSource source : erf) {
				for (ProbEqkRupture rup : source) {
					double weight = rup.getMeanAnnualRate(1);
					sumWeight += weight;
					sumWeightedProps += weight*prop.calcForRuptureAndLocation(rup, loc);
				}
			}
			xyz.set(i, sumWeightedProps/sumWeight);
		});
		return xyz;
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
	
	private static void writeMagDistFuncs(PointSourceCalcERF erf, EvenlyDiscretizedFunc distFunc,
			List<List<Location>> distLocs, IncrementalMagFreqDist mfd, File outputDir, String prefix) throws IOException {
		System.out.println("Calculating mag/dist funcs, will write with prefix "+prefix);
		Preconditions.checkState(distLocs.size() == distFunc.size());
		double[][] rJBs = new double[distFunc.size()][mfd.size()];
		double[][] rRups = new double[distFunc.size()][mfd.size()];
		double[][] rXs = new double[distFunc.size()][mfd.size()];
		double[][] rXFractPos = new double[distFunc.size()][mfd.size()];
		
		IntStream.range(0, distFunc.size()).parallel().forEach(d -> {
			int[] numEntries = new int[mfd.size()];
			
			for (Location loc : distLocs.get(d)) {
				for (ProbEqkSource source : erf) {
					for (ProbEqkRupture rup : source) {
						double mag = rup.getMag();
						int m = mfd.getClosestXIndex(mag);
						RuptureSurface surf = rup.getRuptureSurface();
						rJBs[d][m] += surf.getDistanceJB(loc);
						rRups[d][m] += surf.getDistanceRup(loc);
						double rX = surf.getDistanceX(loc);
						rXs[d][m] += Math.abs(rX);
						if (rX >= 0d)
							rXFractPos[d][m] += 1d;
						numEntries[m]++;
					}
				}
			}
			
			for (int m=0; m<mfd.size(); m++) {
				rJBs[d][m] /= (double)numEntries[m];
				rRups[d][m] /= (double)numEntries[m];
				rXs[d][m] /= (double)numEntries[m];
				rXFractPos[d][m] /= (double)numEntries[m];
			}
		});
		
		System.out.println("Writing mag/dist funcs");
		writeMagDistCSV(new File(outputDir, prefix+"_rJB.csv"), rJBs, distFunc, mfd);
		writeMagDistCSV(new File(outputDir, prefix+"_rRup.csv"), rRups, distFunc, mfd);
		writeMagDistCSV(new File(outputDir, prefix+"_rX.csv"), rXs, distFunc, mfd);
		writeMagDistCSV(new File(outputDir, prefix+"_rX_positive.csv"), rXFractPos, distFunc, mfd);
	}
	
	private static void writeMagDistCSV(File outputFile, double[][] vals,
			EvenlyDiscretizedFunc distFunc, IncrementalMagFreqDist mfd) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>(mfd.size()+1);
		header.add("Horizontal Distance");
		for (int m=0; m<mfd.size(); m++)
			header.add((float)mfd.getX(m)+"");
		csv.addLine(header);
		for (int d=0; d<vals.length; d++) {
			List<String> line = new ArrayList<>(header.size());
			line.add((float)distFunc.getX(d)+"");
			for (int m=0; m<mfd.size(); m++)
				line.add((float)vals[d][m]+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(outputFile);
	}

}