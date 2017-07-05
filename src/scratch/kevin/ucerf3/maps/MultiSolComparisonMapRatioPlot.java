package scratch.kevin.ucerf3.maps;

import java.io.File;
import java.util.List;

import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.peter.curves.ProbOfExceed;
import scratch.peter.nshmp.CurveContainer;
import scratch.peter.nshmp.NSHMP_DataUtils;
import scratch.peter.ucerf3.calc.UC33_MapMaker;

public class MultiSolComparisonMapRatioPlot {

	public static void main(String[] args) {
		TestGrid grid = TestGrid.CA_RELM;
		ProbOfExceed pe = ProbOfExceed.PE2IN50;
		Period p = Period.GM0P00;
		double spacing = 0.1;
		
		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/biasi_downsample_tests/2016_05_18-biasi-downsample-pga");
		
		String refName = "ref_branch_orig";
		List<String> compareNames = Lists.newArrayList("ref_branch_downsampled_inverted_both",
				"ref_branch_downsampled_inverted_ends", "ref_branch_downsampled_inverted_starts");
		
		GeoDataSet ref = loadMap(new File(mainDir, refName).getAbsolutePath(), pe, grid, p, spacing);
		
		boolean log = false;
		
		Region disaggReg = new Region(new Location(34.0657901, -117.8049615), 50);
		
		for (String compareName : compareNames) {
			GeoDataSet comp = loadMap(new File(mainDir, compareName).getAbsolutePath(), pe, grid, p, spacing);
			GeoDataSet xyz = GeoDataSetMath.divide(comp, ref);
			File mapDir = new File(mainDir, compareName+"_ratio");
			Preconditions.checkState(mapDir.exists() || mapDir.mkdir());
			UC33_MapMaker.makeRatioPlot(xyz, spacing, grid.bounds(), mapDir.getAbsolutePath(), compareName,
					log, 0.3, !log, false);
			if (disaggReg != null) {
				double maxVal = 0d;
				Location maxLoc = null;
				for (Location loc : xyz.getLocationList()) {
					double val = xyz.get(loc);
					if (val > maxVal && disaggReg.contains(loc)) {
						maxVal = val;
						maxLoc = loc;
					}
				}
				System.out.println("Max val in disagg region of "+maxVal+" at "+maxLoc);
			}
		}
	}
	
	private static final String S = File.separator;
	
	private static GeoDataSet loadMap(String dir, ProbOfExceed pe,
			TestGrid grid, Period p, double spacing) {
		File curves = new File(dir + S + grid + S + p + S + "curves.csv");
		CurveContainer cc = CurveContainer.create(curves, grid, spacing);
		System.out.println("Loaded "+cc.size()+" curves");
		GriddedRegion gr = grid.grid(spacing);
		return NSHMP_DataUtils.extractPE(cc, gr, pe);
	}

}
