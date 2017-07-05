package scratch.peter.nshmp;

import static scratch.peter.nshmp.NSHMP_UtilsDev.toNum;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.sra.rtgm.RTGM.Frequency;

import scratch.peter.nshmp.CurveContainer.CurveFileProcessor_SHA;

import com.google.common.base.CharMatcher;
import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;

/**
 * Class for storint ghrids of RTGM values.
 *
 * @author Peter Powers
 * @version $Id:$
 */
public class RTGM_Container implements Iterable<Location> {

	private GriddedRegion region;
	private Map<Integer, Double> rtgm1hz;
	private Map<Integer, Double> rtgm5hz;
	
	/**
	 * Returns the hazard curve for the supplied location.
	 * @param loc of interest
	 * @param f frequency of interest (1hz or 5hz)
	 * @return the associated hazard curve
	 * @throws IllegalArgumentException if loc is out of range of curve data
	 */
	public Double getValue(Location loc, Frequency f) {
		int idx = region.indexForLocation(loc);
		Preconditions.checkArgument(
			idx != -1, "Location is out of range: " + loc);
		switch (f) {
			case SA_0P20:
				return rtgm5hz.get(idx);
			case SA_1P00:
				return rtgm1hz.get(idx);
			default:
				return null;
		}
	}
	
	/**
	 * Returns the number of curves stored in this container.
	 * @return the container size
	 */
	public int size() {
		return region.getNodeCount();
	}
	
	@Override
	public Iterator<Location> iterator() {
		return region.iterator();
	}

	/**
	 * Creates a curve container for NSHMP national scale datafrom the supplied
	 * data file and region. The supplied file is assumed to be in the standard
	 * format output by NSHMP fortran 'combine' codes.
	 * 
	 * @param f file
	 * @return a new curve container object
	 */
	public static RTGM_Container create(File f) {
		RTGM_FileProcessor rfp = new RTGM_FileProcessor(
			NSHMP_UtilsDev.getNSHMP_Region());
		RTGM_Container curves = null;
		try {
			curves = Files.readLines(f, Charsets.US_ASCII, rfp);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		return curves;
	}
	
	
	/**
	 * Creates a curve container for a localized area from supplied data file
	 * and region. The data locations should match the nodes in the gridded
	 * region. Results are unspecified if the two do not agree. The suplied
	 * file is assumed to be in an OpenSHA csv format.
	 * 
	 * @param f file
	 * @param tg grid for file
	 * @return a new curve container object
	 */
	public static CurveContainer create(File f, TestGrid tg) {
		CurveFileProcessor_SHA cfp = new CurveFileProcessor_SHA(tg.grid(0.1));
		CurveContainer curves = null;
		try {
			curves = Files.readLines(f, Charsets.US_ASCII, cfp);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		return curves;
	}


	static class RTGM_FileProcessor implements LineProcessor<RTGM_Container> {
			
			private Splitter split;
			private RTGM_Container rc;
			
			RTGM_FileProcessor(GriddedRegion region) {
				split = Splitter.on(CharMatcher.WHITESPACE).omitEmptyStrings();
				rc = new RTGM_Container();
				rc.region = region;
				rc.rtgm1hz = Maps.newHashMapWithExpectedSize(region.getNodeCount());
				rc.rtgm5hz = Maps.newHashMapWithExpectedSize(region.getNodeCount());
			}

			@Override
			public RTGM_Container getResult() {
				return rc;
			}

			@Override
			public boolean processLine(String line) throws IOException {
				if (line.startsWith("#")) return true;
				addValues(line);
				return true;
			}
			
			private void addValues(String line) {
				try {
				Iterator<String> it = split.split(line).iterator();
				// read location
				Location loc = new Location(toNum(it.next()), toNum(it.next()));
				int idx = rc.region.indexForLocation(loc);
				rc.rtgm1hz.put(idx, toNum(it.next()));
				rc.rtgm5hz.put(idx, toNum(it.next()));
				} catch (Exception e) {
					System.out.println(line);
					e.printStackTrace();
				}
			}
		}
	

	/**
	 * Test
	 * @param args
	 */
	public static void main(String[] args) {
		File f = new File("/Volumes/Scratch/nshmp-sources/DesignMap/RTGM.dat");
		Stopwatch sw = Stopwatch.createStarted();
		RTGM_Container rc = RTGM_Container.create(f);
		sw.stop();
		System.out.println("time: " + sw.elapsed(TimeUnit.MILLISECONDS));
		System.out.println("size: " + rc.size());
		System.out.println(rc.getValue(NEHRP_TestCity.MEMPHIS.shiftedLocation(), Frequency.SA_0P20));
		System.out.println(rc.getValue(NEHRP_TestCity.MEMPHIS.shiftedLocation(), Frequency.SA_1P00));
		
//		sw.start();
//		CurveContainer cc = RTGM_Container.create(f, TestGrid.MEMPHIS);
//		sw.stop();
//		System.out.println("time: " + sw.elapsedMillis());
//		System.out.println("size: " + cc.size());
//		System.out.println(cc.getCurve(NEHRP_TestCity.MEMPHIS.location()));
	}

}
