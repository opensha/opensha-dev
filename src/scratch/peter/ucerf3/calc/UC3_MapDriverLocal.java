package scratch.peter.ucerf3.calc;

import static org.opensha.sha.earthquake.param.IncludeBackgroundOption.EXCLUDE;
import static org.opensha.sha.earthquake.param.IncludeBackgroundOption.INCLUDE;

import java.io.File;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.geo.LocationList;
import org.opensha.nshmp2.calc.HazardResultWriter;
import org.opensha.nshmp2.calc.HazardResultWriterLocal;
import org.opensha.nshmp2.calc.ThreadedHazardCalc;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import com.google.common.base.Stopwatch;
import com.google.common.io.Files;

/**
 * Add comments here
 * 
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class UC3_MapDriverLocal {

	private static String S = File.separator;

	UC3_MapDriverLocal(String solPath, String outDir, SourceIMR imr, TestGrid grid,
		double spacing, Period p, IncludeBackgroundOption bg) {

		try {
			outDir += grid + S + p + S;
			File out = new File(outDir, "curves");
			Files.createParentDirs(out);
			HazardResultWriter writer = new HazardResultWriterLocal(out, p);
			LocationList locs = grid.grid(spacing).getNodeList();
			ThreadedHazardCalc calc = new ThreadedHazardCalc(solPath, imr, locs, p,
			false, bg, writer, false);
			calc.calculate(null);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
//		String solPathFM31 = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/tree/2013_01_14-UC32-MEAN_BRANCH_AVG_SOL_FM31.zip";
//		String solPathFM32 = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/tree/2013_01_14-UC32-MEAN_BRANCH_AVG_SOL_FM32.zip";

		String sol = "tmp/UC33/src/bravg/FM-DM/UC33brAvg_FM31_ABM.zip";
		
		String outDir = "tmp/tmp";
		
		TestGrid tg = TestGrid.CA_RELM;
		Period p = Period.GM0P00;
		IncludeBackgroundOption bg = INCLUDE;
		
		SourceIMR imr = SourceIMR.WUS_FAULT_14;
		
		outDir += (bg == INCLUDE) ? "all/" : (bg == EXCLUDE) ? "flt/" : "bg/";
		
		Stopwatch sw = Stopwatch.createStarted();
		new UC3_MapDriverLocal(sol, outDir, imr, tg, 0.1, p, bg);
		sw.stop();

		System.out.println(sw.elapsed(TimeUnit.MINUTES));
	}

}
