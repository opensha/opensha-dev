package scratch.kevin.pointSources.paperFigs2026;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.data.Named;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;

public class CalcPaths {
	
	public static final File PAPER_DIR =new File("/home/kevin/Documents/papers/2026_nshm_grid_seis_dist_corr/papers-2026-nshm-grid-seis-dist-corrs");
	public static final File FIGURES_DIR =new File(PAPER_DIR, "Figures");
	
	public static final File INV_DIR = new File("/data/kevin/nshm23/batch_inversions/");
	
	public static final Models REF_MODEL = Models.FINITE_20X_UNCENTERED;
	public static final Models PROPOSED_MODEL = Models.SPINNING_DIST_5X_UNCENTERED;
	
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
	
	enum Models implements Named {
		AS_PUBLISHED("NSHM23 as-published", "AsPublished",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-as_published-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-as_published-supersample-ba_only")),
		SPINNING_AVG("Spinning average, centered", "AvgCentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-average_spinning-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-average_spinning-supersample-ba_only")),
		SPINNING_AVG_M6("Spinning average, centered, M>6", "AvgCenteredMSix",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-average_spinning_m6-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-average_spinning_m6-supersample-ba_only")),
		SPINNING_AVG_UNCENTERED("Spinning average, uncentered", "AvgUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-average_spinning_along-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-average_spinning_along-supersample-ba_only")),
		FINITE_20X_CENTERED("Finite ruptures (20x), centered", "FiniteCentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_20x-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-finite_20x-supersample-ba_only")),
		FINITE_1X_UNCENTERED("Finite ruptures (1x), uncentered", "FiniteSingleUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_1x_along-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-finite_1x_along-supersample-ba_only")),
		FINITE_1X_UNCENTERED_ALT("Finite ruptures (1x), uncentered, alt random", "FiniteSingleUncenteredAlt",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_1x_along_alt_rand-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-finite_1x_along_alt_rand-supersample-ba_only")),
		FINITE_2X_UNCENTERED("Finite ruptures (2x), uncentered", "FiniteDoubleUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_2x_along-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-finite_2x_along-supersample-ba_only")),
		FINITE_2X_UNCENTERED_ALT("Finite ruptures (2x), uncentered, alt random", "FiniteDoubleUncenteredAlt",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_2x_along_alt_rand-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-finite_2x_along_alt_rand-supersample-ba_only")),
		FINITE_20X_UNCENTERED("Finite ruptures (20x), uncentered", "FiniteUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_20x_along-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-finite_20x_along-supersample-ba_only")),
		FINITE_50X_UNCENTERED("Finite ruptures (50x), uncentered", "FiniteFiftyUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_50x_along-supersample-ba_only"),
				null),
		FINITE_100X_UNCENTERED("Finite ruptures (100x), uncentered", "FiniteHundredUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-finite_100x_along-supersample-ba_only"),
				null),
		SPINNING_DIST_5X_CENTERED("Spinning distribution (5x), centered", "SpinningDistCentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-five_pt_spinning-supersample-ba_only")),
		SPINNING_DIST_5X_UNCENTERED("Spinning distribution (5x), uncentered", "SpinningDistUncentered",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning_along-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-five_pt_spinning_along-supersample-ba_only")),
		SPINNING_DIST_5X_UNCENTERED_M3("Spinning distribution (5x), uncentered, M>3", "SpinningDistUncenteredMThree",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning_along-m3-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-five_pt_spinning_along-m3-supersample-ba_only")),
		SPINNING_DIST_5X_UNCENTERED_M4("Spinning distribution (5x), uncentered, M>4", "SpinningDistUncenteredMFour",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning_along-m4-supersample-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-five_pt_spinning_along-m4-supersample-ba_only")),
		SPINNING_DIST_5X_UNCENTERED_ALT_DEPTH("Spinning distribution (5x), uncentered, alt grid depths", "SpinningDistUncenteredAltDepths",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning_along-alt_grid_depths-supersample-ba_only"),
				null),
		SPINNING_DIST_5X_UNCENTERED_ALT_LENGTH("Spinning distribution (5x), uncentered, alt grid lengths", "SpinningDistUncenteredAltLengths",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning_along-alt_wc_lengths-supersample-ba_only"),
				null),
		SPINNING_DIST_5X_UNCENTERED_NO_SS("Spinning distribution (5x), uncentered, no supersample", "SpinningDistUncenteredNoSS",
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-five_pt_spinning_along-ba_only"),
				new File(INV_DIR, "2025_10_07-nshm23_pt_src_tests-zoom-five_pt_spinning_along-ba_only"));
		
		public final String name;
		public final String texName;
		public final File mapDir;
		public final File zoomMapDir;

		private Models(String name, String texName, File mapDir, File zoomMapDir) {
			this.name = name;
			this.texName = texName;
			this.mapDir = mapDir;
			this.zoomMapDir = zoomMapDir;
		}

		@Override
		public String getName() {
			return name;
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

}
