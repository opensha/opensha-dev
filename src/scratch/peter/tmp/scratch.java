package scratch.peter.tmp;

import static org.opensha.commons.geo.LocationUtils.azimuth;
import static org.opensha.commons.geo.LocationUtils.distanceToLineFast;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.sha.faultSurface.FaultTrace;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class scratch {

	
	public static void main(String[] args) throws Exception {

		GriddedRegion gr = new GriddedRegion(TestGrid.CA_NSHMP.grid(0.1), 0.05, GriddedRegion.ANCHOR_0_0);
		
		System.out.println(gr.getNodeCount());
		
		System.out.println(TestGrid.CA_NSHMP.grid(0.05).getNodeCount());
		System.out.println(TestGrid.CA_NSHMP_E.grid(0.1).getNodeCount());
		System.out.println(TestGrid.CA_NSHMP_N.grid(0.1).getNodeCount());
		System.out.println(TestGrid.CA_NSHMP_NE.grid(0.1).getNodeCount());
				
		int tot = TestGrid.CA_NSHMP.grid(0.1).getNodeCount() + 
				TestGrid.CA_NSHMP_E.grid(0.1).getNodeCount() +
				TestGrid.CA_NSHMP_N.grid(0.1).getNodeCount() +
				TestGrid.CA_NSHMP_NE.grid(0.1).getNodeCount();
		System.out.println(tot);
		
		// test
//		Range<Double> r = Range.open(1.0,23.0);
//		System.out.println(r);
		
//		Period p = Period.GM1P00;
//		DiscretizedFunc f = p.getFunction();
//		f = Utils.getExceedProbabilities(f, -0.23573904887559857, 0.6874658901925597, false, 0.0);
//		f.scale(0.01);
//		System.out.println(f);
//		System.out.println(ProbOfExceed.get(f, ProbOfExceed.PE2IN50));
		
//		double[] seq = DataUtils.buildSequence(0, 10, 1.00000000000001, true);
//		System.out.println(Arrays.toString(seq));
//		double[] seq2 = DataUtils.buildSequence(0, 10, 0.99999999999998, true);
//		System.out.println(Arrays.toString(seq2));
		
//		double[] p1 = {1,2,3,4,5};
//		int idx = firstZeroValue(p1);
//		double[] p2 = Arrays.copyOf(p1,idx);
//		System.out.println(Arrays.toString(p2));
//		GriddedRegion gr = TestGrid.LITTLE_SALMON.grid(0.02);
//		System.out.println(gr.getNodeCount());
		

//		tmp();
		
//		double pp = 6. * (4. / 3.);
//		System.out.println(pp);
//		String dir = "/Users/pmpowers/projects/OpenSHA/tmp/UC3maps/mapsUC32b/UC32littleSalmon/multi-UC32/PE1IN100-sol9";
//		consolidateMaps(dir);
		
//		for (NEHRP_TestCity city : NEHRP_TestCity.values()) {
//			System.out.println(city);
//		}
		
//		UCERF3_FaultSysSol_ERF erf = UC3_CalcUtils.getUC3_ERF(
//			"/Users/pmpowers/projects/OpenSHA/tmp/invSols/conv/2013_02_01-ucerf3p2-convergence_bundle5_sol14.zip",
//			IncludeBackgroundOption.INCLUDE, false, true, 1.0);
//		erf.updateForecast();
		
		//		try {
//			String srcDir = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/conv/";
//			String solPath = srcDir + "2013_02_01-ucerf3p2-convergence_bundle5_COMPOUND_SOL.zip";
//			CompoundFaultSystemSolution cfss = UC3_CalcUtils.getCompoundSolution(solPath);
//			for (LogicTreeBranch br : cfss.getBranches()) {
//				System.out.println(br.buildFileName());
//				
//			}
//			
//			int idx = 0;
//			for (FaultSystemSolution fss : cfss) {
//				System.out.println(fss.getClass());
//				String solName = "2013_02_01-ucerf3p2-convergence_bundle5_sol" + (idx++) + ".zip";
//				SimpleFaultSystemSolution sfss = new SimpleFaultSystemSolution(fss);
//				File solZip = new File(srcDir + solName);
//				sfss.toZipFile(solZip);
//			}
//			
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
	}
	
	private static double tmpRx(FaultTrace trace, Location loc) {
		// avg strike 
		double rX_avgStrk = distanceToLineFast(trace.first(), trace.last(), loc);
		double avgStrkAzStart = azimuth(trace.first(), trace.last());
		double avgStrkAzEnd = azimuth(trace.last(), trace.first());
		
//		Location preStart = location(trace.first(), avgStrkAzStart, 10.0);
//		Location postEnd = location(trace.last(), avgStrkAzEnd, 10.0);
//		
		// could just start with 
		// closest point
		int closeIdx = trace.closestPoint(loc);
		Location closest = trace.get(closeIdx);
		// set segment azimuths away from closets point and using
		// avg strike off the ends of the fault
		double seg1az = (closeIdx == 0) ?
			avgStrkAzStart : azimuth(closest, trace.get(closeIdx - 1));
		double seg2az = (closeIdx == trace.size() - 1) ? 
			avgStrkAzEnd : azimuth(closest, trace.get(closeIdx + 1));
		
		// 
		
//		double SEG1_AZ = LocationUtils.
		return 1;
	}

	private static int firstZeroValue(double[] data) {
		int idx = 0;
		for (double d : data) {
			if (d > 0.0) {
				idx++;
			} else {
				break;
			}
		}
		Preconditions.checkArgument(
			idx > 1, "Curve must have more than two non-zero y-values: " +
				Arrays.toString(data));
		return idx;
	}

	private static void consolidateMaps(String path) throws IOException {
		File srcDir = new File(path);
		File destDir = new File(path, "composite");
		destDir.mkdir();
		File[] mapDirs = srcDir.listFiles();
		for (File mapDir : mapDirs) {
			if (!mapDir.isDirectory() || mapDir.equals(destDir)) continue;
			String mapName = mapDir.getName();
			File srcMap = new File(mapDir, "map.pdf");
			File destMap = new File(destDir, mapName + ".pdf");
			Files.copy(srcMap, destMap);
		}
	}
	
	private static void tmp() throws IOException, DocumentException {
		String srcSol = "tmp/UC33/src/tree/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip";
		String outDir = "tmp/UC33/src/bg";
		String U2_sol = "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU2";
		String U3_sol = "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3";
		
		CompoundFaultSystemSolution cfss = UC3_CalcUtils.getCompoundSolution(srcSol);
		LogicTreeBranch U2_branch = LogicTreeBranch.fromFileName(U2_sol);
		LogicTreeBranch U3_branch = LogicTreeBranch.fromFileName(U3_sol);
		InversionFaultSystemSolution U2_fss = cfss.getSolution(U2_branch);
		InversionFaultSystemSolution U3_fss = cfss.getSolution(U3_branch);
		
		File U2_file = new File(outDir, U2_sol + ".zip");
		File U3_file = new File(outDir, U3_sol + ".zip");
		
		FaultSystemIO.writeSol(U2_fss, U2_file);
		FaultSystemIO.writeSol(U3_fss, U3_file);
	}
}
