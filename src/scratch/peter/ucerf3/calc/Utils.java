package scratch.peter.ucerf3.calc;

import static org.apache.commons.lang3.StringUtils.substringAfterLast;
import static org.apache.commons.lang3.StringUtils.substringBeforeLast;

import java.io.File;
import java.io.IOException;

import org.apache.commons.lang3.ArrayUtils;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.Interpolate;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;

import scratch.peter.curves.ProbOfExceed;

import com.google.common.io.Files;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class Utils {

	private static final String CONSOL_DIR = "consolidated";
	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		String dirName = "/Users/pmpowers/Documents/OpenSHA/RTGM/deaggTmp/TMP/UC3FM3P1/EUREKA";
//		String fileName = "DisaggregationPlot.png";
//		consolidateFiles(fileName, dirName);
		
//		interpolateTest();
//		findMissing();
		combineCurves();
//		GriddedRegion gr = TestGrid.CA_NSHMP.grid(0.05);
//		System.out.println(gr.getNodeCount());
	}
	
	// finds missing locations in big NSHMP13B runs
	private static void findMissing() {
		String path = "tmp/UC33/maps/src/NSHMP14/FLT_E/mean_ucerf3_sol/CA_NSHMP_E/tmp/GM0P20";
		GriddedRegion grid = TestGrid.CA_NSHMP_E.grid(0.1);
		
		String format = "%.3f";
		int count = 0;
		for (Location loc : grid) {
			String fName = String.format(format, loc.getLatitude()) + "_" +
					String.format(format, loc.getLongitude()) + "_prob.txt";
			File f = new File(path, fName);
			if (!f.exists()) {
				count++;
				System.out.println(loc);
			}
		}
		System.out.println(count);
	}
	
	private static void combineCurves() {
		String path = "tmp/UC33/maps/src/NSHMP14/FLT_E/mean_ucerf3_sol/CA_NSHMP_E/tmp/GM0P20";
		UC3_CalcMPJ_Map.aggregateResults(new File(path), Period.GM0P20, false);
	}
	
	
	
	public static void consolidateFiles(
			String fileName,
			String dirName) {
		
		try {
			File dirIn = new File(dirName);
			File dirOut = new File(dirIn, "consolidated");
			dirOut.mkdir();
			File[] dirList = dirIn.listFiles();
			for (File d : dirList) {
				if (!d.isDirectory()) continue;
				if (d.getName().equals(CONSOL_DIR)) continue;
				String idx = substringAfterLast(d.getName(), "-");
				File plotIn = new File(d, fileName);
				String plotOutName = substringBeforeLast(fileName, ".") + "-" +
					idx + "." + substringAfterLast(fileName, ".");
				File plotOut = new File(dirOut, plotOutName);
				Files.copy(plotIn, plotOut);
			}			
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	private static void interpolateTest() {
		double X[] = new double[] {0.0050,0.0070,0.0098,0.0137,0.0192,0.0269,0.0376,0.0527,0.0738,0.103,0.145,0.203,0.284,0.397,0.556,0.778,1.09,1.52,2.13};
		double Y[] = new double[] {0.6075407,0.53996414,0.45535624,0.3607489,0.26581264,0.18253244,0.11834708,0.073424034,0.04487007,0.027603487,0.016950918,0.010716081,0.007054321,0.0048437594,0.0032797337,0.001998431,9.994777E-4,3.8728004E-4,1.0545944E-4};

		ArbitrarilyDiscretizedFunc f = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<X.length; i++) {
			f.set(X[i], Y[i]);
		}
		
		ProbOfExceed pe = ProbOfExceed.PE1IN100;
//		double fVal = f.getFirstInterpolatedX_inLogXLogYDomain(pe.annualRate());
		
		ArrayUtils.reverse(X);
		ArrayUtils.reverse(Y);
		double iVal = Interpolate.findLogLogY(Y, X, pe.annualRate());
		
//		System.out.println(fVal);
		System.out.println(iVal);
		
	}

}
