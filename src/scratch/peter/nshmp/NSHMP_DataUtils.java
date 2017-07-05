/**
 * 
 */
package scratch.peter.nshmp;

import static org.opensha.nshmp2.util.Period.*;
import static org.opensha.sra.rtgm.RTGM.Frequency.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.EnumSet;
import java.util.Formatter;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.NSHMP_Utils;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;

import com.google.common.collect.Maps;

import scratch.peter.curves.ProbOfExceed;

/**
 * Class of static methods for managing, cleaning, and creating USGS NSHMP
 * datasets.
 *
 * @author Peter Powers
 * @version $Id:$
 */
public class NSHMP_DataUtils {

	private static final String NSHMP_SRC_PATH = "/Volumes/Scratch/nshmp-sources/";
	private static final String SHA_SRC_PATH = "/Volumes/Scratch/nshmp-opensha-";
	private static final String SEP = File.separator;
	private static final String RTGM_OUT = "RTGM.dat";
	private static final String RTGM_HEADER = "# lat lon 1hz 5hz\n";
	private static final String RTGM_FORMAT = "%s %s %6.2f %6.2f\n";

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		// national scale data
//		deriveDesignMapRTGM();
//		generateRTGM_NSHMP();
		
		// regional grids
		generateRTGM_SHA();
		
	}
	
	
	
	
	
	
	/*
	 ************************************************************************
	 * Creates a files of lat, lon, prob-rtgm-1hz, prob-rtgm-5hz from the
	 * NSHMP hazard curve datasets:
	 * 
	 * http://earthquake.usgs.gov/hazards/designmaps/datasets/#!nehrp-2009
	 * 
	 * Files are large hence simultaneous streaming through bufferedReaders.
	 ************************************************************************
	 */
	
	private static String dataName = "curves.dat";
	private static String csvName = "curves.csv";
	private static String dataR3 = "DataR3";
	private static String hTool = "HazardTool";
	private static String fortran1 = "FortranLatest";
	private static String fortran2 = "FortranUpdate";
	private static String shaDev = "dev";
	private static String shaTrunk = "trunk";
	private static Set<TestGrid> grids = TestGrid.getLocals();
	private static final double DEF_GM_1HZ = 0.0033;
	private static final double DEF_GM_5HZ = 0.0055;
	
	
	public static void generateRTGM_SHA() {
//		TestGrid tg = TestGrid.LOS_ANGELES;
//		processRTGM_SHA(SHA_SRC_PATH + shaDev + SEP + tg.name() + SEP, tg);
		for (TestGrid grid : grids) {
			processRTGM_SHA(SHA_SRC_PATH + shaDev + SEP + grid.name() + SEP, grid);
			processRTGM_SHA(SHA_SRC_PATH + shaTrunk + SEP + grid.name() + SEP, grid);
		}
	}
	
	private static void processRTGM_SHA(String dataPath, TestGrid tg) {
		File f1hz = new File(dataPath, GM1P00 + SEP + csvName);
		File f5hz = new File(dataPath, GM0P20 + SEP + csvName);
		File out = new File(dataPath, RTGM_OUT);

		GriddedRegion gr = tg.grid(0.1);

		System.out.println("Processing 1hz...");
		CurveContainer cc_1hz = CurveContainer.create(f1hz, tg, 0.1);
		Map<Integer, Double> rtgmMap1hz = calcRTGM(gr, cc_1hz, SA_1P00, DEF_GM_1HZ);

		System.out.println("Processing 5hz...");
		CurveContainer cc_5hz = CurveContainer.create(f5hz, tg, 0.1);
		Map<Integer, Double> rtgmMap5hz = calcRTGM(gr, cc_5hz, SA_0P20, DEF_GM_5HZ);
		
		System.out.println("Writing output...");
		try {
			Formatter formatter = new Formatter(out);
			formatter.format("%s", RTGM_HEADER);
		
			for (int i=0; i<gr.getNodeCount(); i++) {
				Location loc = gr.locationForIndex(i);
				formatter.format(RTGM_FORMAT,
					String.format("%.2f", loc.getLatitude()),
					String.format("%.2f", loc.getLongitude()),
					rtgmMap1hz.get(i), rtgmMap5hz.get(i));
			}
			formatter.flush();
			formatter.close();
		} catch (FileNotFoundException fnfe) {
			fnfe.printStackTrace();
		}
	}

	/**
	 * Utility method to generate RTGM values from hazard curve data sets and
	 * output them to file.
	 */
	public static void generateRTGM_NSHMP() {
		processRTGM_NSHMP(NSHMP_SRC_PATH + dataR3 + SEP);
		processRTGM_NSHMP(NSHMP_SRC_PATH + hTool + SEP);
		processRTGM_NSHMP(NSHMP_SRC_PATH + fortran1 + SEP);
		processRTGM_NSHMP(NSHMP_SRC_PATH + fortran2 + SEP);
	}
		
	private static void processRTGM_NSHMP(String dataPath) {
		File f1hz = new File(dataPath, GM1P00 + SEP + dataName);
		File f5hz = new File(dataPath, GM0P20 + SEP + dataName);
		File out = new File(dataPath, RTGM_OUT);

		GriddedRegion gr = NSHMP_UtilsDev.getNSHMP_Region(0.1);

		System.out.println("Processing 1hz...");
		CurveContainer cc_1hz = CurveContainer.create(f1hz);
		Map<Integer, Double> rtgmMap1hz = calcRTGM(gr, cc_1hz, SA_1P00, DEF_GM_1HZ);

		System.out.println("Processing 5hz...");
		CurveContainer cc_5hz = CurveContainer.create(f5hz);
		Map<Integer, Double> rtgmMap5hz = calcRTGM(gr, cc_5hz, SA_0P20, DEF_GM_5HZ);
		
		System.out.println("Writing output...");
		try {
			Formatter formatter = new Formatter(out);
			formatter.format("%s", RTGM_HEADER);
		
			for (int i=0; i<gr.getNodeCount(); i++) {
				Location loc = gr.locationForIndex(i);
				formatter.format(RTGM_FORMAT,
					String.format("%.2f", loc.getLatitude()),
					String.format("%.2f", loc.getLongitude()),
					rtgmMap1hz.get(i), rtgmMap5hz.get(i));
			}
			formatter.flush();
			formatter.close();
		} catch (FileNotFoundException fnfe) {
			fnfe.printStackTrace();
		}
	}
	
	/*
	 * Returns an index map of rtgm values indexed according to the indices of 
	 * the supplied gridded region. 'defaultGM' is a floor value to be used 
	 */
	private static Map<Integer, Double> calcRTGM(GriddedRegion region,
			CurveContainer cc, Frequency f, double defaultGM) {
		
		// init thread mgr
		int numProc = Runtime.getRuntime().availableProcessors();
		ExecutorService ex = Executors.newFixedThreadPool(numProc);
		CompletionService<RTGM> ecs = new ExecutorCompletionService<RTGM>(ex);

		// process two files linearly
		
		// store results in maps based on location index
		int numLocs = region.getNodeCount();
		
		Map<Integer, Double> rtgm_vals = Maps.newHashMapWithExpectedSize(numLocs);
		
		for (int i=0; i<numLocs; i++) {
			Location loc = region.locationForIndex(i);
			DiscretizedFunc hc = cc.getCurve(loc);
			RTGM rtgm = RTGM.createIndexed(hc, f, 0.8, i);
			ecs.submit(rtgm);
		}
		System.out.println("   Jobs submitted: " + numLocs);
		ex.shutdown();
		
		int bad = 0;
		int good = 0;
		for (int i = 0; i < numLocs; i++) {
			if (i % 10000 == 0) System.out.println("   Jobs completed: " + i);
			
			RTGM rtgm = null;
			try {
				rtgm = ecs.take().get();
			} catch (ExecutionException ee) {
				bad++;
				continue;
				// ignore execution exceptions as they indicate
				// a bad hazard curve
			} catch (InterruptedException ie) {
				ie.printStackTrace();
				ex.shutdownNow();
				return null;
			}
			rtgm_vals.put(rtgm.index(), rtgm.get());
			good++;
		}

		// block while working
		try {
			ex.awaitTermination(48, TimeUnit.HOURS);
		} catch (InterruptedException ie) {
			ex.shutdownNow();
			return null;
		}
		
		// Any missing indices in result map (skipped because hazard curve was
		// bad) should dnow be filled in with 1.0. Good values should be
		// converted to %g
		for (int i = 0; i < numLocs; i++) {
			Double val = rtgm_vals.get(i);
			if (val == null) val = defaultGM;
			val *= 100;
			rtgm_vals.put(i,val);
		}
		
		System.out.println("   Finished good="+good+" bad="+bad+" total="+(good+bad));
		return rtgm_vals;
	}
	
	
	
	
		
	/*
	 ************************************************************************
	 * Creates a file of lat, lon, prob-rtgm-1hz, prob-rtgm-5hz from the
	 * 2009 NEHRP raw data files available here:
	 * 
	 * http://earthquake.usgs.gov/hazards/designmaps/datasets/#!nehrp-2009
	 * 
	 * Files are large hence simultaneous streaming through bufferedReaders.
	 ************************************************************************
	 */
	
	private static String in1hz_rc = "2009_NEHRP-1p0s_Risk_Coefficients.txt";
	private static String in5hz_rc = "2009_NEHRP-0p2s_Risk_Coefficients.txt";
	private static String in1hz2p50 = "2009_NEHRP-1p0s_Uniform_Hazard_Ground_Motion.txt";
	private static String in5hz2p50 = "2009_NEHRP-0p2s_Uniform_Hazard_Ground_Motion.txt";
	
	/**
	 * Utility method that mines 4 huge files (1hz-2%50, 5hz-2%50, 1hz-rc,
	 * 5hz-rc) to compute and output the probabilistic component of RTGM
	 * for the conterminous US.
	 */
	public static void deriveDesignMapRTGM() {
		String dataPath = NSHMP_SRC_PATH + "DesignMap/";
		
		BufferedReader br_p1hz = null;
		BufferedReader br_p5hz = null;
		BufferedReader br_rc1hz = null;
		BufferedReader br_rc5hz = null;
		
		String p1hz_line = null;
		String p5hz_line = null;
		String rc1hz_line = null;
		String rc5hz_line = null;
		try {
			br_p1hz = new BufferedReader(new FileReader(new File(dataPath, in1hz2p50)));
			br_p5hz = new BufferedReader(new FileReader(new File(dataPath, in5hz2p50)));
			br_rc1hz = new BufferedReader(new FileReader(new File(dataPath, in1hz_rc)));
			br_rc5hz = new BufferedReader(new FileReader(new File(dataPath, in5hz_rc)));
			
			Formatter formatter = new Formatter(new File (dataPath, RTGM_OUT));
			formatter.format("%s", RTGM_HEADER);
			
			int count = 0;
			
			while ((p1hz_line = br_p1hz.readLine()) != null) {
				p5hz_line = br_p5hz.readLine();
				rc1hz_line = br_rc1hz.readLine();
				rc5hz_line = br_rc5hz.readLine();
				processLine(formatter, p1hz_line, p5hz_line, rc1hz_line, rc5hz_line);
				count++;
				if (count % 1000000 == 0) System.out.println(count);
			}
			formatter.flush();
			formatter.close();

		} catch (FileNotFoundException fnfe) {
			fnfe.printStackTrace();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		} catch (NullPointerException npe) {
			System.out.println("p1hz_line: " + p1hz_line);
			System.out.println("p5hz_line: " + p5hz_line);
			System.out.println("rc1hz_line: " + rc1hz_line);
			System.out.println("rc5hz_line: " + rc5hz_line);
			npe.printStackTrace();
		}
		IOUtils.closeQuietly(br_p1hz);
		IOUtils.closeQuietly(br_p5hz);
		IOUtils.closeQuietly(br_rc1hz);
		IOUtils.closeQuietly(br_rc5hz);
	}
	
	private static void processLine(Formatter f, String p1hz, String p5hz, String rc1hz,
			String rc5hz) {
		if (p1hz.startsWith("%")) return;
		String[] p1hzVals = StringUtils.split(p1hz);
		String lat = p1hzVals[0];
		String lon = p1hzVals[1];
		if (!validLoc(lat, lon)) return;
		double p1hzVal = NSHMP_Utils.readDouble(p1hz, 2);
		double p5hzVal = NSHMP_Utils.readDouble(p5hz, 2);
		double rc1hzVal = NSHMP_Utils.readDouble(rc1hz, 2);
		double rc5hzVal = NSHMP_Utils.readDouble(rc5hz, 2);
		double rtgm1hz = p1hzVal * rc1hzVal;
		double rtgm5hz = p5hzVal * rc5hzVal;
		f.format(RTGM_FORMAT, lat, lon, rtgm1hz, rtgm5hz);
	}
	
	private static boolean validLoc(String lat, String lon) {
		int latVal = (int) Math.rint((Double.parseDouble(lat) * 100));
		int lonVal = (int) Math.rint((Double.parseDouble(lon) * 100));
		return (latVal % 10 == 0) && (lonVal % 10 == 0);
	}

	
	
	
	/*
	 ************************************************************************
	 * Utility methods
	 ************************************************************************
	 */
	
	
	/**
	 * Utility method to extract probabilities of exceedance from hazard curves
	 * and return them in a GeoDataSet. PE's are calculated fairly quicklly and
	 * can be computed on the fly for the purposes of making maps.
	 * @param cc hazard curve container
	 * @param region gridded region for which curves/PE's should be extracted
	 * @param pe the target PE
	 * @return the PE results in a GeoDataSet
	 */
	public static GeoDataSet extractPE(CurveContainer cc,
			GriddedRegion region, ProbOfExceed pe) {
		GriddedGeoDataSet gDat = new GriddedGeoDataSet(region, true);
		double targetRate = pe.annualRate();
		for (Location loc : region) {
			
			// two possible problems can arise here: (1) the region to extract
			// has points outside the data source region and (2) the data
			// source curve may have uniform very low values and fails
			// interpolation. In either case set ground motion value to 0
			
			double gm = 0;
			DiscretizedFunc f = null;
			try {
				f = cc.getCurve(loc);
			} catch (NullPointerException npe) {
				System.out.println("Missing location? " + loc);
				npe.printStackTrace();
			}
			
			try {
				gm = f.getFirstInterpolatedX_inLogXLogYDomain(targetRate);
			} catch (Exception e) {
				System.out.println("Problem location: " + loc);
				// do nothing; let gm be 0
			}
			gDat.set(loc, gm);
		}
		return gDat;
	}
	
	/**
	 * Flavor of method above that just takes a probability; assumes supplied
	 * curves are probabilistic.
	 * @param cc
	 * @param region
	 * @param prob
	 * @return
	 */
	public static GeoDataSet extractPE(CurveContainer cc,
			GriddedRegion region, double prob) {
		GriddedGeoDataSet gDat = new GriddedGeoDataSet(region, true);
		for (Location loc : region) {
			
			// two possible problems can arise here: (1) the region to extract
			// has points outside the data source region and (2) the data
			// source curve may have uniform very low values and fails
			// interpolation. In either case set ground motion value to 0
			
			double gm = 0;
			DiscretizedFunc f = cc.getCurve(loc);
			try {
				gm = f.getFirstInterpolatedX_inLogXLogYDomain(prob);
			} catch (Exception e) {
				System.out.println("Problem location: " + loc);
				// do nothing; let gm be 0
			}
			gDat.set(loc, gm);
		}
		return gDat;
	}


	/**
	 * Utility method to extract RTGM values from a RTGM data set
	 * and return them in a GeoDataSet.
	 * @param rc RTGM data container
	 * @param region gridded region for which RTGM values should be extracted
	 * @param f the frequency of interest
	 * @return the PE results in a GeoDataSet
	 */
	public static GeoDataSet extractRTGM(RTGM_Container rc,
			GriddedRegion region, Frequency f) {
		GriddedGeoDataSet gDat = new GriddedGeoDataSet(region, true);
		for (Location loc : region) {
			double val = 0;
			try {
				val = rc.getValue(loc, f);
			} catch (Exception e) {
				System.out.println(loc);
				// do nothing; let gm be 0
			}
			gDat.set(loc, val);
		}
		return gDat;
	}
	
}
