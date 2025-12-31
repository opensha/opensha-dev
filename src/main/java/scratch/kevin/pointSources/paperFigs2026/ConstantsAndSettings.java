package scratch.kevin.pointSources.paperFigs2026;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Leonard2010_MagLengthRelationship;
import org.opensha.commons.data.Named;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedFiniteRuptureSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class ConstantsAndSettings {
	
	public static final File PAPER_DIR =new File("/home/kevin/Documents/papers/2026_nshm_grid_seis_dist_corr/papers-2026-nshm-grid-seis-dist-corrs");
	public static final File FIGURES_DIR =new File(PAPER_DIR, "Figures");
	
	public static final File INV_DIR = new File("/data/kevin/nshm23/batch_inversions/");
	
	public static final File ORIG_SOL_DIR = new File(INV_DIR, "2024_02_02-nshm23_branches-WUS_FM_v3");
	public static final File ORIG_SOL_FILE = new File(ORIG_SOL_DIR, "results_WUS_FM_v3_branch_averaged_gridded.zip");
	
	public static final Models REF_FINITE_MODEL = Models.FINITE_20X_UNCENTERED;
	public static final Models PROPOSED_DIST_CORR_MODEL = Models.SPINNING_DIST_5X_UNCENTERED;
	public static final Models PROPOSED_FULL_MODEL = Models.SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5;
	
	public static final String HAZARD_MODEL_PREFIX = "2025_12_29-nshm23_pt_src_tests";
	public static final String HAZARD_MODEL_ZOOM_PREFIX = HAZARD_MODEL_PREFIX+"-zoom";
	
	public static final GriddedRegion FULL_GRID_REG;
	static {
		try {
			FULL_GRID_REG = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.1, GriddedRegion.ANCHOR_0_0);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	public static final GriddedRegion ZOOM_GRID_REG = new GriddedRegion(
			new Region(new Location(36, -120), new Location(39, -117)), 0.02, GriddedRegion.ANCHOR_0_0);
	
	public enum Models implements Named {
		/*
		 * Point source correction tests
		 */
		AS_PUBLISHED("NSHM23 as-published", "AsPublished",
				5d, PointSourceDistanceCorrections.NSHM_2013, 6d, true, null,
				"as_published"),
		SPINNING_AVG_CENTERED_M6("Spinning average, centered, Correct M>6", "AvgCenteredMFiveCorrMSix",
				5d, PointSourceDistanceCorrections.AVERAGE_SPINNING_CENTERED, 6d, true, null,
				"spinning_average-centered-m5-corrM6"),
		SPINNING_AVG_CENTERED_M5("Spinning average, centered, M>5", "AvgCenteredMFive",
				5d, PointSourceDistanceCorrections.AVERAGE_SPINNING_CENTERED, 5d, true, null,
				"spinning_average-centered-m5"),
		FINITE_20X_CENTERED("Finite ruptures (20x), centered", "FiniteCenteredTwenty",
				5d, 20, false, 5d, true, null,
				"finite-20x-centered-m5"),
		// Reference to compare point source calculations
		FINITE_20X_UNCENTERED("Finite ruptures (20x), uncentered", "FiniteUncenteredTwenty",
				5d, 20, true, 5d, true, null,
				"finite-20x-uncentered-m5"),
		FINITE_50X_UNCENTERED("Finite ruptures (50x), uncentered", "FiniteUncenteredFifty",
				5d, 50, true, 5d, true, null,
				"finite-50x-uncentered-m5"),
		FINITE_100X_UNCENTERED("Finite ruptures (100x), uncentered", "FiniteUncenteredHundred",
				5d, 100, true, 5d, true, null,
				"finite-100x-uncentered-m5"),
		FINITE_1X_UNCENTERED("Finite ruptures (1x), uncentered", "FiniteSingleUncentered",
				5d, 1, true, 5d, true, null,
				"finite-1x-uncentered-m5"),
		FINITE_1X_UNCENTERED_ALT_RAND("Finite ruptures (1x), uncentered, alt random", "FiniteSingleUncenteredAltRand",
				5d, 1, true, 5d, true, 111234511l,
				"finite-1x-uncentered-m5-alt_rand"),
		FINITE_2X_UNCENTERED("Finite ruptures (2x), uncentered", "FiniteDoubleUncentered",
				5d, 2, true, 5d, true, null,
				"finite-2x-uncentered-m5"),
		FINITE_2X_UNCENTERED_ALT_RAND("Finite ruptures (2x), uncentered, alt random", "FiniteDoubleUncenteredAltRand",
				5d, 2, true, 5d, true, 221234522l,
				"finite-2x-uncentered-m5-alt_rand"),
		SPINNING_DIST_5X_CENTERED("Spinning distribution (5x), centered", "SpinningDistCentered",
				5d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_CENTERED, 5d, true, null,
				"spinning_dist-centered-m5"),
		/*
		 * Proposed point source correction
		 */
		SPINNING_DIST_5X_UNCENTERED("Spinning distribution (5x), uncentered", "SpinningDistUncentered",
				5d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST, 5d, true, null,
				"spinning_dist-uncentered-m5"),
		/*
		 * Further recommendations/tests, all using proposed point source calculation
		 */
		SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR("Spinning distribution (5x), uncentered, updated Ztor", "SpinningDistUncenteredModZtor",
				5d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST, 5d, true, null,
				"spinning_dist-uncentered-m5-mod_ztor") {

					@Override
					public Function<GridSourceList, GridSourceList> getGridModFunction() {
						return original -> { return new GridPropertyUpdate(original, null, true);};
					}
			
		},
		SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN("Spinning distribution (5x), uncentered, updated Ztor and lengths", "SpinningDistUncenteredModZtorLen",
				5d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST, 5d, true, null,
				"spinning_dist-uncentered-m5-mod_ztor_len") {

					@Override
					public Function<GridSourceList, GridSourceList> getGridModFunction() {
						return original -> { return new GridPropertyUpdate(original, LEONARD, true);};
					}
			
		},
		SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M4("Spinning distribution (5x), uncentered, updated Ztor and lengths, M>4", "SpinningDistUncenteredModZtorLenMFour",
				4d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST, 5d, true, null,
				"spinning_dist-uncentered-m4-mod_ztor_len") {

					@Override
					public Function<GridSourceList, GridSourceList> getGridModFunction() {
						return original -> { return new GridPropertyUpdate(original, LEONARD, true);};
					}
			
		},
		SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3p5("Spinning distribution (5x), uncentered, updated Ztor and lengths, M>3.5", "SpinningDistUncenteredModZtorLenMThreeFive",
				3.5d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST, 5d, true, null,
				"spinning_dist-uncentered-m3p5-mod_ztor_len") {

					@Override
					public Function<GridSourceList, GridSourceList> getGridModFunction() {
						return original -> { return new GridPropertyUpdate(original, LEONARD, true);};
					}
			
		},
		SPINNING_DIST_5X_UNCENTERED_MOD_ZTOR_LEN_M3("Spinning distribution (5x), uncentered, updated Ztor and lengths, M>3", "SpinningDistUncenteredModZtorLenMThree",
				3d, PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST, 5d, true, null,
				"spinning_dist-uncentered-m3-mod_ztor_len") {

					@Override
					public Function<GridSourceList, GridSourceList> getGridModFunction() {
						return original -> { return new GridPropertyUpdate(original, LEONARD, true);};
					}
			
		},
		;
		
		public final String name;
		public final String texName;
		public final BackgroundRupType bgRupType;
		public final double pointMinMag;
		public final PointSourceDistanceCorrections distCorr;
		public final boolean finiteUncentered;
		public final int finiteNum;
		public final double pointFiniteMinMag;
		public final boolean supersample;
		public final Long customRandSeed;
		
		public final File mapDir;
		public final File zoomMapDir;
		
		// point surface constructor
		private Models(String name, String texName,
				double pointMinMag,
				PointSourceDistanceCorrections distCorr,
				double pointFiniteMinMag,
				boolean supersample,
				Long customRandSeed,
				String mapDirSuffix) {
			this(name, texName, BackgroundRupType.POINT, pointMinMag, distCorr, -1, false,
					pointFiniteMinMag, supersample, customRandSeed, mapDirSuffix);
		}
		
		// finite surface constructor
		private Models(String name, String texName,
				double pointMinMag,
				int finiteNum,
				boolean finiteUncentered,
				double pointFiniteMinMag,
				boolean supersample,
				Long customRandSeed,
				String mapDirSuffix) {
			this(name, texName, BackgroundRupType.FINITE, pointMinMag, null, finiteNum, finiteUncentered,
					pointFiniteMinMag, supersample, customRandSeed, mapDirSuffix);
		}

		private Models(String name, String texName,
				BackgroundRupType bgRupType,
				double pointMinMag,
				PointSourceDistanceCorrections distCorr,
				int finiteNum,
				boolean finiteUncentered,
				double pointFiniteMinMag,
				boolean supersample,
				Long customRandSeed,
				String mapDirSuffix) {
			this.name = name;
			this.texName = texName;
			this.bgRupType = bgRupType;
			this.pointMinMag = pointMinMag;
			this.distCorr = distCorr;
			this.finiteNum = finiteNum;
			this.finiteUncentered = finiteUncentered;
			this.pointFiniteMinMag = pointFiniteMinMag;
			this.supersample = supersample;
			this.customRandSeed = customRandSeed;
			this.mapDir = new File(INV_DIR, HAZARD_MODEL_PREFIX+"-"+mapDirSuffix);
			this.zoomMapDir = new File(INV_DIR, HAZARD_MODEL_ZOOM_PREFIX+"-"+mapDirSuffix);
		}

		@Override
		public String getName() {
			return name;
		}
		
		public String getGridArgs() {
			String argz = "--point-source-type "+bgRupType.name()+" --point-min-mag "+(float)pointMinMag;
			if (bgRupType == BackgroundRupType.FINITE) {
				argz += " --point-finite-num-rand-surfaces "+finiteNum;
				if (finiteUncentered)
					argz += " --point-finite-sample-along-strike --point-finite-sample-down-dip";
			} else {
				argz += " --dist-corr "+distCorr.name();
			}
			argz += " --point-finite-min-mag "+(float)pointFiniteMinMag;
			if (supersample)
				argz += " --supersample-finite";
			if (customRandSeed != null)
				argz += " -D"+PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME+"="+customRandSeed;
			
			return argz;
		}
		
		public GriddedSeismicitySettings getGridProps() {
			GriddedSeismicitySettings settings = GriddedSeismicitySettings.DEFAULT
					.forSurfaceType(bgRupType)
					.forMinimumMagnitude(pointMinMag);
			if (bgRupType == BackgroundRupType.FINITE) {
				GriddedFiniteRuptureSettings finiteSettings = GriddedFiniteRuptureSettings.DEFAULT.forNumSurfaces(finiteNum);
				
				if (finiteUncentered)
					finiteSettings = finiteSettings.forSampleAlongStrike(true).forSampleDownDip(true);
				
				settings = settings.forFiniteRuptureSettings(finiteSettings);
			} else {
				settings = settings.forDistanceCorrection(distCorr.get());
			}
			settings = settings.forPointSourceMagCutoff(pointFiniteMinMag);
			GridCellSupersamplingSettings ssSettings = null;
			if (supersample) {
				ssSettings = GridCellSupersamplingSettings.DEFAULT.forApplyToFinite(true);
			}
			if (customRandSeed == null)
				System.clearProperty(PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME);
			else
				System.setProperty(PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME, customRandSeed+"");
			return settings.forSupersamplingSettings(ssSettings);
		}
		
		public Function<GridSourceList, GridSourceList> getGridModFunction() {
			return null;
		}
	}
	
	static {
		// sanity checks
		HashSet<String> names = new HashSet<>(Models.values().length);
		HashSet<String> texNames = new HashSet<>(Models.values().length);
		HashSet<String> dirNames = new HashSet<>(Models.values().length);
		for (Models model : Models.values()) {
			Preconditions.checkState(!names.contains(model.name), "Duplicate model name: %s", model.name);
			Preconditions.checkState(!texNames.contains(model.texName), "Duplicate model tex name: %s", model.texName);
			String dirName = model.mapDir.getName();
			Preconditions.checkState(!dirNames.contains(dirName), "Duplicate model directory name: %s", dirName);
			names.add(model.name);
			texNames.add(model.texName);
			dirNames.add(dirName);
		}
	}
	
	public static GriddedGeoDataSet buildLandMask(GriddedRegion gridReg) {
		GriddedGeoDataSet mask;
		mask = new GriddedGeoDataSet(gridReg);
		for (XY_DataSet polBound : PoliticalBoundariesData.loadDefaultOutlines(gridReg)) {
			LocationList list = new LocationList();
			for (Point2D pt : polBound)
				list.add(new Location(pt.getY(), pt.getX()));
			Region reg = new Region(list, BorderType.MERCATOR_LINEAR);
			double halfLatSpacing = gridReg.getLatSpacing()*0.5;
			double halfLonSpacing = gridReg.getLonSpacing()*0.5;
			for (int i=0; i<gridReg.getNodeCount(); i++) {
				if (mask.get(i) == 1d)
					continue;
				Location center = gridReg.getLocation(i);
				Location[] testLocs = {
						center,
						new Location(center.lat+halfLatSpacing, center.lon+halfLonSpacing),
						new Location(center.lat+halfLatSpacing, center.lon-halfLonSpacing),
						new Location(center.lat-halfLatSpacing, center.lon+halfLonSpacing),
						new Location(center.lat-halfLatSpacing, center.lon-halfLonSpacing),
				};
				for (Location loc : testLocs) {
					if (reg.contains(loc)) {
						mask.set(i, 1d);
						break;
					}
				}
			}
		}
		System.out.println("Mask kept "+(float)(100d*mask.getSumZ()/gridReg.getNodeCount())+" % of locations");
		for (int i=0; i<mask.size(); i++)
			if (mask.get(i) < 1d)
				mask.set(i, Double.NaN);
		return mask;
	}
	
	public static class GridPropertyUpdate extends GridSourceList.DynamicallyBuilt {

		private GridSourceList original;
		private Function<GriddedRuptureProperties, Double> ml;
		private boolean updateZtor;

		public GridPropertyUpdate(GridSourceList original, Function<GriddedRuptureProperties, Double> ml, boolean updateZtor) {
			super(original.getTectonicRegionTypes(), original.getGriddedRegion(), original.getRefMFD());
			this.original = original;
			this.ml = ml;
			this.updateZtor = updateZtor;
		}

		@Override
		public int getNumSources() {
			return original.getNumSources();
		}
		
		@Override
		public int getLocationIndexForSource(int sourceIndex) {
			return original.getLocationIndexForSource(sourceIndex);
		}
		
		@Override
		public TectonicRegionType tectonicRegionTypeForSourceIndex(int sourceIndex) {
			return original.tectonicRegionTypeForSourceIndex(sourceIndex);
		}
		
		@Override
		public Set<Integer> getAssociatedGridIndexes(int sectionIndex) {
			return original.getAssociatedGridIndexes(sectionIndex);
		}
		
		@Override
		protected List<GriddedRupture> buildRuptures(TectonicRegionType tectonicRegionType, int gridIndex) {
			List<GriddedRupture> orig = original.getRuptures(tectonicRegionType, gridIndex);
			if (orig == null || orig.isEmpty())
				return orig;
			List<GriddedRupture> ret = new ArrayList<>(orig.size());
			for (GriddedRupture rup : orig) {
				GriddedRuptureProperties props = rup.properties;
				double mag = props.magnitude;
				double dipRad = Math.toRadians(props.dip);
				double length;
				if (ml == null)
					length = props.length;
				else
					length = ml.apply(props);
				double upper, lower;
				if (updateZtor) {
					double[] depths = calcDepths(mag, length, dipRad);
					upper = depths[0];
					lower = depths[1];
				} else {
					upper = props.upperDepth;
					lower = props.lowerDepth;
				}
				props = new GriddedRuptureProperties(mag, props.rake, props.dip,
						props.strike, props.strikeRange, upper, lower, length,
						Double.NaN, Double.NaN, tectonicRegionType);
				rup = new GriddedRupture(rup.gridIndex, rup.location, props, rup.rate,
						rup.associatedSections, rup.associatedSectionFracts);
				ret.add(rup);
			}
			return ret;
		}
		
	}
	
	private static final double zCenter = 5d;	// km
	private static final double zTopMin = 1d;	// km
	private static final double zBotMax = 14d;	// km
	private static final double maxVertThickness = zBotMax - zTopMin;	// 13 km
	
	public static double[] calcDepths(double mag, double length, double dipRad) {
		// Same aspect control as before
		double aspectWidthDD = length / 1.5;	// down-dip width (km)

		// Cap by maximum possible down-dip width given seismogenic thickness
		double sinDip = Math.sin(dipRad);
		double maxWidthDD = maxVertThickness / sinDip;
		double ddWidth = Math.min(aspectWidthDD, maxWidthDD);

		// Convert to vertical thickness (km)
		double vertThickness = ddWidth * sinDip;

		// Center vertically at 5 km
		double upper = zCenter - 0.5 * vertThickness;
		double lower = zCenter + 0.5 * vertThickness;

		// Enforce saturation bounds by shifting the interval if needed.
		// (Keep thickness fixed unless it exceeds the max thickness.)
		if (vertThickness >= maxVertThickness) {
			upper = zTopMin;
			lower = zBotMax;
		} else if (upper < zTopMin) {
			upper = zTopMin;
			lower = upper + vertThickness;
		} else if (lower > zBotMax) {
			lower = zBotMax;
			upper = lower - vertThickness;
		}
		return new double[] {upper, lower};
	}
	
	public static final Function<GriddedRuptureProperties, Double> LEONARD =
			P -> P.dip < 89.99 ? Leonard2010_MagLengthRelationship.DIP_SLIP.getMedianLength(P.magnitude)
					: Leonard2010_MagLengthRelationship.STRIKE_SLIP.getMedianLength(P.magnitude);
	
	/**
	 * WC 94 RLD (subsurface rupture length) formula
	 */
	public static final Function<GriddedRuptureProperties, Double> WC94_RLD =
			P -> Math.pow(10.0,-2.44+0.59*P.magnitude);

}
