package scratch.peter.nshmp;

import static org.opensha.nshmp2.util.SourceRegion.WUS;
import static org.opensha.nshmp2.util.SourceType.FAULT;
import static org.opensha.nshmp2.util.SourceType.*;

import java.awt.geom.Point2D;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.Interpolate;
import org.opensha.commons.util.XMLUtils;
import org.opensha.nshmp2.erf.NSHMP2008;
import org.opensha.nshmp2.erf.NSHMP_ListERF;
import org.opensha.nshmp2.erf.source.FaultERF;
import org.opensha.nshmp2.erf.source.FaultSource;
import org.opensha.nshmp2.erf.source.GridERF;
import org.opensha.nshmp2.erf.source.NSHMP_ERF;
import org.opensha.nshmp2.erf.source.SubductionERF;
import org.opensha.nshmp2.util.NSHMP_Utils;
import org.opensha.nshmp2.util.SourceType;
import org.opensha.sha.magdist.GaussianMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

/**
 * Add comments here
 * 
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class CSEP_Forecast {

	private static final String TEMPLATE_PATH = "tmp/CSEP/csep-forecast-template-m5.xml";
	private static final String DATE_FMT = "yyyy-MM-dd'T'HH:mm:ss'Z'";
	private static final String issueDate, startDate, endDate;
	private static final String VERSION = "3.0";
	private static final String AUTHOR = "USGS";
	private static final double NODE_SIZE = 0.1;
	private static final double MFD_MIN = 5.0;
	private static final double MFD_MAX = 9.0;
	private static final int MFD_NUM = 41;

	private NSHMP_ListERF erf;
	private Document forecastDoc;
	private Map<Integer, Location> locIndexMap;
	private Map<Integer, Element> cellElementMap;
	private Map<Integer, IncrementalMagFreqDist> cellMfdMap;

	static {
		TimeZone tz = TimeZone.getTimeZone("UTC");
		SimpleDateFormat sdf = new SimpleDateFormat(DATE_FMT);
		sdf.setTimeZone(tz);
		GregorianCalendar cal = new GregorianCalendar(tz);
		cal.set(2010, 0, 1, 0, 0, 0);
		issueDate = sdf.format(cal.getTime());
		cal.set(2008, 0, 1, 0, 0, 0);
		startDate = sdf.format(cal.getTime());
		cal.set(2013, 11, 31, 23, 59, 59);
		endDate = sdf.format(cal.getTime());
	}

	/**
	 * Builds a CSEP compatible forecast (XML format) and wirtes it to the
	 * supplied file.
	 * 
	 * @param name
	 * @param erf
	 * @param out
	 * @throws Exception if problem encountered
	 */
	public static void createForecast(String name, NSHMP_ListERF erf, File out)
			throws Exception {
		new CSEP_Forecast(name, erf, out);
	}

	private CSEP_Forecast(String name, NSHMP_ListERF erf, File out) throws Exception {
		System.out.println("Initializing...");
		this.erf = erf;
		Files.createParentDirs(out);
		forecastDoc = XMLUtils.loadDocument(TEMPLATE_PATH);
		Element cellElem = initDocument(name);
		// use index map instead of Locations as keys in latter two maps
		locIndexMap = Maps.newHashMap();
		cellElementMap = Maps.newHashMap();
		cellMfdMap = Maps.newHashMap();
		initCellMaps(cellElem);
		System.out.println("Processing...");
		processForecast();
		checkEmpties();
		System.out.println("Writing MFDs...");
		writeMFDs();
		System.out.println("Writing XML...");
		XMLUtils.writeDocumentToFile(out, forecastDoc);
		System.out.println("Done.");
	}

	/*
	 * Sets header elements of forecast.
	 */
	private Element initDocument(String name) {
		Element fd = forecastDoc.getRootElement().element("forecastData");
		fd.element("modelName").setText(name);
		fd.element("version").setText(VERSION);
		fd.element("author").setText(AUTHOR);
		fd.element("issueDate").setText(issueDate);
		fd.element("forecastStartDate").setText(startDate);
		fd.element("forecastEndDate").setText(endDate);
		fd.element("lastMagBinOpen").setText("1");
		return fd.element("depthLayer");
	}

	/*
	 * Builds maps of Locations indices; the indices ar eused to reference maps
	 * of ELements and MFDs.
	 */
	private void initCellMaps(Element cellWrapper) {
		List cellElements = cellWrapper.elements("cell");
		for (int i = 0; i < cellElements.size(); i++) {
			Element cell = (Element) cellElements.get(i);
			Location loc = new Location(Double.valueOf(cell
				.attributeValue("lat")), Double.valueOf(cell
				.attributeValue("lon")));
			locIndexMap.put(i, loc);
			cellElementMap.put(i, cell);
			cellMfdMap.put(i, null);
		}
	}

	private void processForecast() throws Exception {
		// init thread mgr
		int numProc = Runtime.getRuntime().availableProcessors();
		ExecutorService ex = Executors.newFixedThreadPool(numProc);
		CompletionService<ProcessorResult> ecs = new ExecutorCompletionService<ProcessorResult>(
			ex);

		Set<Integer> indices = locIndexMap.keySet();
		for (Integer idx : indices) {
			Location loc = locIndexMap.get(idx);
			LocationProcessor locProc = new LocationProcessor(erf, loc, idx);
			ecs.submit(locProc);
		}
		ex.shutdown();

		for (int i = 0; i < indices.size(); i++) {
			ProcessorResult pr = ecs.take().get();
			cellMfdMap.put(pr.idx, pr.mfd);
		}
	}

	private static class LocationProcessor implements Callable<ProcessorResult> {

		private NSHMP_ListERF erfList;
		private Location loc;
		private int idx;

		LocationProcessor(NSHMP_ListERF erfList, Location loc, int idx) {
			this.erfList = erfList;
			this.loc = loc;
			this.idx = idx;
		}

		@Override
		public ProcessorResult call() {
			SummedMagFreqDist sum = newSummedForecastMFD();
			for (NSHMP_ERF erf : erfList) {

				SourceType type = erf.getSourceType();

				if (type.equals(GRIDDED)) {
					GridERF gerf = (GridERF) erf;
					IncrementalMagFreqDist gmfd = getOffsetGridMFD(gerf, loc,
						NODE_SIZE);
					// getOffsetGrid may return null mfds
					if (gmfd == null) continue;
					sum.addIncrementalMagFreqDist(gmfd);

				} else if (type.equals(FAULT)) {
					FaultERF ferf = (FaultERF) erf;
					IncrementalMagFreqDist fmfd = getFaultMFD(ferf, loc);
					// getFaultMFD will always return a (possibly empty) sum
					sum.addIncrementalMagFreqDist(fmfd);

				} else if (type.equals(SUBDUCTION)) {
					SubductionERF serf = (SubductionERF) erf;
					IncrementalMagFreqDist fmfd = getSubductionMFD(serf, loc);
					// getSubductionMFD will always return a (possibly empty)
					// sum
					sum.addIncrementalMagFreqDist(fmfd);

				} else {
					throw new UnsupportedOperationException(
						"Invalid Source Type");
				}

			}
			return createResult(idx, sum);
		}

	}

	private static ProcessorResult createResult(int idx,
			IncrementalMagFreqDist mfd) {
		ProcessorResult pr = new ProcessorResult();
		pr.idx = idx;
		pr.mfd = mfd;
		return pr;
	}

	private static class ProcessorResult {
		int idx;
		IncrementalMagFreqDist mfd;
	}

	private void checkEmpties() {
		for (int idx : cellMfdMap.keySet()) {
			double rate = cellMfdMap.get(idx).getCumRate(0);
			if (rate == 0.0) {
				System.out.println("Empty MFD: " + locIndexMap.get(idx));
			}
		}
	}

	/*
	 * Transfers raw MFDs to value fields in forecast XML elements.
	 */
	private void writeMFDs() {
		for (int idx : cellMfdMap.keySet()) {
			IncrementalMagFreqDist mfd = cellMfdMap.get(idx);
			Element e = cellElementMap.get(idx);
			writeMFD(mfd, e);
		}
	}

	private static void writeMFD(IncrementalMagFreqDist mfd, Element e) {
		List bins = e.elements("bin");
		for (int i = 0; i < bins.size(); i++) {
			Element bin = (Element) bins.get(i);
			double mag = Double.valueOf(bin.attributeValue("m"));
			double rate = mfd.getY(mag);
			bin.setText(Double.toString(rate));
		}
	}

	public static IncrementalMagFreqDist getFaultMFD(FaultERF erf, Location loc) {
		SummedMagFreqDist sum = newSummedForecastMFD();
		Site nodeSite = new Site(loc);
		Region nodeRegion = getNodeRegion(loc, NODE_SIZE, NODE_SIZE);

		for (FaultSource fSrc : erf.getSources()) {

			// skip if likely outside node
			double quickDist = fSrc.getMinDistance(nodeSite);
			if (quickDist > 10.0) continue;

			// compute fraction of source inside node
			double nodeWt = fSrc.getSourceSurface()
				.getFractionOfSurfaceInRegion(nodeRegion);
			if (nodeWt == 0.0) continue;

			// mfds are weighted already
			for (IncrementalMagFreqDist mfd : fSrc.getMFDs()) {
				boolean moBalance = mfd.getClass().equals(
					GaussianMagFreqDist.class);
				// no need to clone original as resample does not alter src
				IncrementalMagFreqDist mfdResamp = resample(mfd,
					newIncrForecastMFD(), moBalance, true);
				mfdResamp.scale(nodeWt); // scale by fraction in node
				sum.addIncrementalMagFreqDist(mfdResamp);
			}
		}
		sum.scale(erf.getSourceWeight()); // scale by the erf wt
		return sum;
	}

	public static IncrementalMagFreqDist getSubductionMFD(SubductionERF erf,
			Location loc) {
		SummedMagFreqDist sum = newSummedForecastMFD();
		Site nodeSite = new Site(loc);
		Region nodeRegion = getNodeRegion(loc, NODE_SIZE, NODE_SIZE);

		for (FaultSource fSrc : erf.getSources()) {

			// skip if likely outside node
			double quickDist = fSrc.getMinDistance(nodeSite);
			if (quickDist > 10.0) continue;

			// compute fraction of source inside node
			double nodeWt = fSrc.getSourceSurface()
				.getFractionOfSurfaceInRegion(nodeRegion);
			if (nodeWt == 0.0) continue;

			// mfds are weighted already
			for (IncrementalMagFreqDist mfd : fSrc.getMFDs()) {
				// no need to clone original as resample does not alter src
				IncrementalMagFreqDist mfdResamp = resample(mfd,
					newIncrForecastMFD(), false, true);
				mfdResamp.scale(nodeWt); // scale by fraction in node
				sum.addIncrementalMagFreqDist(mfdResamp);
			}
		}
		sum.scale(erf.getSourceWeight()); // scale by the erf wt
		return sum;
	}

	/**
	 * Returns an MFD for a {@code Location} that is assumed to be at the center
	 * of 4 grid nodes in the supplied {@code erf}. Method returns {@code null}
	 * if any of the 4 nodes are {@code null}.
	 * 
	 * @param erf to process
	 * @param loc to process
	 * @param spacing of erf grid
	 * @return an MFD
	 */
	public static IncrementalMagFreqDist getOffsetGridMFD(GridERF erf,
			Location loc, double spacing) {
		// skip locations outside area of grid with non-zero a-values
//		if (!erf.getBorder().contains(loc)) return null;
		// Update 4/29/2013: the above line causes some edge cases to be skipped
		// and has been commented out in favor of just grabbing 0-rate mfds by
		// processing every gridERF at every location of interest.
		double w = spacing / 2;
		double lat = loc.getLatitude();
		double lon = loc.getLongitude();
		List<Location> locs = Lists.newArrayList(
			new Location(lat - w, lon - w), new Location(lat + w, lon - w),
			new Location(lat + w, lon + w), new Location(lat - w, lon + w));
		SummedMagFreqDist sum = newSummedForecastMFD();
		for (Location cLoc : locs) {
			IncrementalMagFreqDist mfd = erf.getMFD(cLoc);
			// bail if all four surrounding nodes are not available
			// Update 4/29/2013: this is creating some cells with mfds that
			// are truncated at 6.5 because cells in the larger grids (e.g.
			// EXTmap) are masked to accomodate 6.5+ fixed-strike grid sources
			// but the fixed-strike sources are not getting fully apportioned
			// at the edges. Now we use fractional mfds if available.
			if (mfd == null) continue; //return null;
			// no need to clone original as resample does not alter src
			IncrementalMagFreqDist mfdResamp = resample(mfd,
				newIncrForecastMFD(), false, true);
			sum.addIncrementalMagFreqDist(mfdResamp);
//			if (LocationUtils.areSimilar(loc, tmpLoc)) {
//				System.out.println(erf.getName() + " " + tmpLoc);
//				System.out.println(mfdResamp);
//			}
		}
		sum.scale(0.25); // 4 mfds were used so scale to 1/4
		sum.scale(erf.getSourceWeight()); // scale again by the erf wt
		return sum;
	}
	
	private static final Location tmpLoc = new Location(39.85,-120.05);

	/**
	 * Returns a CSEP forecast MFD with all rates initialized to 0.
	 * 
	 * @return a summed MFD that can be used to compile MFDs from different
	 *         sources
	 */
	public static SummedMagFreqDist newSummedForecastMFD() {
		return new SummedMagFreqDist(MFD_MIN, MFD_MAX, MFD_NUM);
	}

	/**
	 * Returns a CSEP forecast MFD with all rates initialized to 0.
	 * 
	 * @return a summed MFD that can be used to compile MFDs from different
	 *         sources
	 */
	public static IncrementalMagFreqDist newIncrForecastMFD() {
		return new IncrementalMagFreqDist(MFD_MIN, MFD_MAX, MFD_NUM);
	}

	/**
	 * MFD resampler. Algorithm determines the {@code dest} magnitudes bins that
	 * span the entire range of {@code src} magnitudes. It determines event
	 * rates via linear or logY interpolation of {@code src} rates, fills out
	 * {@code dest} with rates resampled from {@code src}, and scales
	 * {@code dest} to have the same cumulative rate or Mo rate as original.
	 * 
	 * This could be publicised; compare with resampling in MFD api (I believe
	 * MFDs force rates into nearby bins resulting in odd resampling; NOTE there
	 * is an NSHMP customization to force events into the M=5 bin; M=5.05 rounds
	 * up to M=5.1
	 * 
	 * @param src MFD
	 * @param dest MFD, rate values will be overwritten
	 * @param momentBalance {@code true} for Mo preservation, {@code false} for
	 *        cumulative rate
	 * @param logInterp {@code true} for logY interpolation, {@code false} for
	 *        linear
	 * @return a reference to the supplied {@code dest} MFD
	 */
	public static IncrementalMagFreqDist resample(IncrementalMagFreqDist src,
			IncrementalMagFreqDist dest, boolean momentBalance,
			boolean logInterp) {

		// Handle case where src only has one mag - put rate in closest bin
		// Single mag mfds that are high (e.g. M=9.2 cascadia) will go in
		// highest dest bin (M=9.0)
		if (src.size() == 1) {
			int destIdx = dest.getClosestXIndex(src.getX(0));
			dest.set(destIdx, src.getY(0));
			return dest;
		}

		// src mfd points used as basis for interpolation
		double[] srcMags = new double[src.size()];
		double[] srcRates = new double[src.size()];
		int idx = 0;
		for (Point2D p : src) {
			srcMags[idx] = p.getX();
			srcRates[idx++] = p.getY();
		}

		// iterate target dest indices
		// the slight offset on the low side below guarantess that the CSEP
		// M=5 bin is filled, otherwise the NSHMP M=5.05 events all go into
		// the CSEP M=5.1 bin. This effectively increases the NSHMP mfd by
		// one bin, but the cumulative rate stays the same.

		int minDestIdx = dest
			.getClosestXIndex(src.getMinMagWithNonZeroRate() - 0.000001);
		int maxDestIdx = dest.getClosestXIndex(src.getMaxMagWithNonZeroRate());
		for (int destIdx = minDestIdx; destIdx <= maxDestIdx; destIdx++) {
			// min and max indices are already clamped to dest min max
			double destMag = dest.getX(destIdx);
			double destRate = (logInterp) ? Interpolate.findLogY(srcMags,
				srcRates, destMag) : Interpolate.findY(srcMags, srcRates,
				destMag);
			dest.set(destIdx, destRate);
		}
		if (momentBalance) {
			dest.scaleToTotalMomentRate(src.getTotalMomentRate());
		} else {
			dest.scaleToCumRate(0, src.getCumRate(0));
		}
		return dest;
	}

	/**
	 * Returns a {@code Region} marks the border of grid node identified by the
	 * supplied location (center of the node).
	 * 
	 * @param p node center
	 * @param w node width (in decimal degrees)
	 * @param h node height (in decimal degrees)
	 * @return the node {@code Region}
	 */
	public static Region getNodeRegion(Location p, double w, double h) {
		double halfW = w / 2;
		double halfH = h / 2;
		double nodeLat = p.getLatitude();
		double nodeLon = p.getLongitude();
		LocationList locs = new LocationList();
		locs.add(new Location(nodeLat + halfH, nodeLon + halfW)); // top right
		locs.add(new Location(nodeLat - halfH, nodeLon + halfW)); // bot right
		locs.add(new Location(nodeLat - halfH, nodeLon - halfW)); // bot left
		locs.add(new Location(nodeLat + halfH, nodeLon - halfW)); // top left
		return new Region(locs, BorderType.MERCATOR_LINEAR);
	}

	public static void main(String[] args) {
		File out;
		NSHMP2008 erf;
		String name;

		try {
			out = new File("tmp/CSEP/USGS_NSHMP_CA_2008r3.xml");
			erf = NSHMP2008.createCaliforniaCSEP();
			erf.updateForecast();
			System.out.println(erf);
			name = "USGS NSHMP CA 2008r3";
			createForecast(name, erf, out);

			out = new File("tmp/CSEP/MEAN_UCERF2_2008.xml");
			erf = NSHMP2008.createCalifornia();
			erf.updateForecast();
			System.out.println(erf);
			name = "MEAN UCERF2 2008";
			createForecast(name, erf, out);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

// <cell mask='0'....> if no rates for bin

// <?xml version='1.0' encoding='UTF-8'?>
// <CSEPForecast xmlns='http://www.scec.org/xml-ns/csep/forecast/0.1'>
// <forecastData publicID='smi:org.scec/csep/forecast/1'>
// <modelName>unknown</modelName>
// <version>1.0</version>
// <author>CSEP</author>
// <issueDate>2005-06-18T10:30:00Z</issueDate>
// <forecastStartDate>2006-01-01T00:00:00Z</forecastStartDate>
// <forecastEndDate>2011-01-01T00:00:00Z</forecastEndDate>
// <defaultCellDimension latRange='0.1' lonRange='0.1'/>
// <defaultMagBinDimension>0.1</defaultMagBinDimension>
// <lastMagBinOpen>1</lastMagBinOpen>
// <depthLayer max='30.0' min='0.0'>
// <cell lat='40.15' lon='-125.35'>
// <bin m='5.0'>0.0</bin>
// <bin m='5.1'>0.0</bin>
// <bin m='5.2'>0.0</bin>
// <bin m='5.3'>0.0</bin>
// <bin m='5.4'>0.0</bin>
// <bin m='5.5'>0.0</bin>
// <bin m='5.6'>0.0</bin>
// <bin m='5.7'>0.0</bin>
// <bin m='5.8'>0.0</bin>
// <bin m='5.9'>0.0</bin>
// <bin m='6.0'>0.0</bin>
// <bin m='6.1'>0.0</bin>
// <bin m='6.2'>0.0</bin>
// <bin m='6.3'>0.0</bin>
// <bin m='6.4'>0.0</bin>
// <bin m='6.5'>0.0</bin>
// <bin m='6.6'>0.0</bin>
// <bin m='6.7'>0.0</bin>
// <bin m='6.8'>0.0</bin>
// <bin m='6.9'>0.0</bin>
// <bin m='7.0'>0.0</bin>
// <bin m='7.1'>0.0</bin>
// <bin m='7.2'>0.0</bin>
// <bin m='7.3'>0.0</bin>
// <bin m='7.4'>0.0</bin>
// <bin m='7.5'>0.0</bin>
// <bin m='7.6'>0.0</bin>
// <bin m='7.7'>0.0</bin>
// <bin m='7.8'>0.0</bin>
// <bin m='7.9'>0.0</bin>
// <bin m='8.0'>0.0</bin>
// <bin m='8.1'>0.0</bin>
// <bin m='8.2'>0.0</bin>
// <bin m='8.3'>0.0</bin>
// <bin m='8.4'>0.0</bin>
// <bin m='8.5'>0.0</bin>
// <bin m='8.6'>0.0</bin>
// <bin m='8.7'>0.0</bin>
// <bin m='8.8'>0.0</bin>
// <bin m='8.9'>0.0</bin>
// <bin m='9.0'>0.0</bin>
// </cell>

