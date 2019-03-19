package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;
import org.opensha.sha.simulators.utils.SimulatorUtils;
import org.opensha.sha.util.NEHRP_TestCity;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RSQSimBBP_Config {
	
	public static BBP_PlanarSurface estimateBBP_PlanarSurface(SimulatorEvent event) {
		double length = SimulatorUtils.estimateRuptureLength(event);
		double zTOR = Double.POSITIVE_INFINITY;
		double zBOR = Double.NEGATIVE_INFINITY;
		ArrayList<SimulatorElement> elems = event.getAllElements();
		for (SimulatorElement e : elems) {
			for (Location loc : e.getVertices()) {
				double z = loc.getDepth();
				zTOR = Math.min(zTOR, z);
				zBOR = Math.max(zBOR, z);
			}
		}
		if (zTOR < 0)
			zTOR = 0;
//		System.out.println("zTOR = "+zTOR);
//		System.out.println("zBOR = "+zBOR);
		List<Double> strikes = new ArrayList<>();
		List<Double> rakes = new ArrayList<>();
		List<Double> dips = new ArrayList<>();
		for (SimulatorElement e : elems) {
			FocalMechanism mech = e.getFocalMechanism();
			double strike = mech.getStrike();
			double rake = mech.getRake();
			double dip = mech.getDip();
			if (dip == 90 && (rake == 180 || rake == -180)) {
				
			}
			strikes.add(strike);
			rakes.add(rake);
			dips.add(dip);
		}
		double aveRake = FaultUtils.getInRakeRange(FaultUtils.getAngleAverage(rakes));
		double aveStrike = FaultUtils.getAngleAverage(strikes);
		double aveDip = 0d;
		for (double dip : dips)
			aveDip += dip;
		aveDip /= dips.size();
		FocalMechanism mech = new FocalMechanism(aveStrike, aveDip, aveRake);
		
		double width = (zBOR - zTOR)/Math.sin(Math.toRadians(aveDip));
		
		double aveLat = 0d;
		double aveLon = 0d;
		for (SimulatorElement elem : elems) {
			Location loc = elem.getCenterLocation();
			aveLat += loc.getLatitude();
			aveLon += loc.getLongitude();
		}
		aveLat /= elems.size();
		aveLon /= elems.size();
		
		Location topCenter = new Location(aveLat, aveLon, zTOR);
		
		return new BBP_PlanarSurface(topCenter, length, width, mech);
	}
	
	public static BBP_PlanarSurface planarEquivalentU3Surface(RSQSimCatalog catalog, RSQSimEvent event) {
		return planarEquivalentU3Surface(catalog, event, MIN_SUB_SECT_FRACT, ADJ_WIDTH_MATCH_AREA);
	}
	
	public static BBP_PlanarSurface planarEquivalentU3Surface(RSQSimCatalog catalog, RSQSimEvent event,
			double minFractForInclusion, boolean adjWidthMatchArea) {
		EqkRupture rup = catalog.getMappedSubSectRupture(event, minFractForInclusion);
		RuptureSurface u3Surf = rup.getRuptureSurface();
		
		Location firstLoc = u3Surf.getFirstLocOnUpperEdge();
		Location lastLoc = u3Surf.getLastLocOnUpperEdge();
		
		LocationVector vector = LocationUtils.vector(firstLoc, lastLoc);
		double strike = vector.getAzimuth();
		double rake = rup.getAveRake();
		double dip = u3Surf.getAveDip();
		
		FocalMechanism mech = new FocalMechanism(strike, dip, rake);
		
		double length = vector.getHorzDistance();
		
		LocationVector halfVector = new LocationVector(vector.getAzimuth(), 0.5*length, 0d);
		Location topCenter = LocationUtils.location(firstLoc, halfVector);
		topCenter = new Location(topCenter.getLatitude(), topCenter.getLongitude(), u3Surf.getAveRupTopDepth());
		
		double width = u3Surf.getAveWidth();
		
		if (adjWidthMatchArea) {
			double eventArea = event.getArea()*1e-6; // to km^2
			width = Math.min(eventArea/length, MAX_ADJ_WIDTH);
		}
		
		return new BBP_PlanarSurface(topCenter, length, width, mech);
	}
	
	private static double absAngleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		angleDiff = Math.abs(angleDiff);
		return angleDiff;
	}
	
	public static BBP_SourceFile buildBBP_Source(SimulatorEvent event, BBP_PlanarSurface surface, int seed) {
		ArrayList<SimulatorElement> elems = event.getAllElements();
		Preconditions.checkState(elems.size() > 0, "No elements for event %s", event.getID());
		
		Location hypo = null;
		double earliestTime = Double.POSITIVE_INFINITY;
		for (EventRecord rec : event) {
			List<SimulatorElement> patches = rec.getElements();
			double[] patchTimes = rec.getElementTimeFirstSlips();
			for (int i=0; i<patches.size(); i++) {
				if (patchTimes[i] < earliestTime) {
					earliestTime = patchTimes[i];
					hypo = patches.get(i).getCenterLocation();
				}
			}
		}
		Preconditions.checkNotNull(hypo, "Couldn't detect hypocenter for event %s.",
				event.getID());
		LocationVector hypoVector = LocationUtils.vector(surface.getTopCenter(), hypo);
		double hypoAlongStrike = hypoVector.getHorzDistance();
		if (hypoAlongStrike > surface.getLength())
			hypoAlongStrike = surface.getLength();
		double angleDiff = absAngleDiff(hypoVector.getAzimuth(), surface.getFocalMechanism().getStrike());
		if (angleDiff > 90)
			// it's reverse
			hypoAlongStrike = -hypoAlongStrike;
		double hypoDownDip = Math.abs(hypoVector.getVertDistance());
		if (hypoDownDip > surface.getWidth())
			hypoDownDip = surface.getWidth();
		
		double meanArea = 0d;
		for (SimulatorElement e : elems)
			meanArea += e.getArea()/1000000d; // to km^2
		meanArea /= elems.size();
		double sqSideLen = Math.sqrt(meanArea);
		
//		System.out.println("Mean area: "+meanArea+", dx: "+sqSideLen);
		
		return new BBP_SourceFile(surface, event.getMagnitude(), hypoAlongStrike, hypoDownDip, sqSideLen, sqSideLen, 0.15, seed);
//		return new BBP_SourceFile(event.getMagnitude(), topCenter, length, width, mech, hypoAlongStrike, hypoDownDip,
//				sqSideLen, sqSideLen, 0.15, seed);
	}
	
	static final String EVENT_SRF_DIR_NAME = "event_srfs";
	
	static File getEventSrcFile(RSQSimCatalog catalog, int eventID) {
		File srfDir = new File(catalog.getCatalogDir(), EVENT_SRF_DIR_NAME);
		Preconditions.checkState(srfDir.exists() || srfDir.mkdir());
		return new File(srfDir, "event_"+eventID+".src");
	}
	
	static File getEventSRFFile(RSQSimCatalog catalog, int eventID, SRFInterpolationMode interp, double dt) {
		return getEventSRFFile(catalog, eventID, interp, dt, 1d, false);
	}
	
	static File getEventSRFFile(RSQSimCatalog catalog, int eventID, SRFInterpolationMode interp, double dt,
			double timeScale, boolean scaleVelocities) {
		File srfDir = new File(catalog.getCatalogDir(), EVENT_SRF_DIR_NAME);
		Preconditions.checkState(srfDir.exists() || srfDir.mkdir());
		String prefix = "event_"+eventID+"_"+(float)dt+"s_"+interp.name();
		if (timeScale != 1d) {
			prefix += "_timeScale"+(float)timeScale;
			if (scaleVelocities)
				prefix += "_velScale";
		}
		return new File(srfDir, prefix+".srf");
	}
	
	static File getEventBBPDir(RSQSimCatalog catalog, int eventID, SRFInterpolationMode interp, double dt) {
		return getEventBBPDir(catalog, eventID, interp, dt, 1d, false);
	}
	
	static File getEventBBPDir(RSQSimCatalog catalog, int eventID, SRFInterpolationMode interp, double dt,
			double timeScale, boolean scaleVelocities) {
		File srfDir = new File(catalog.getCatalogDir(), EVENT_SRF_DIR_NAME);
		Preconditions.checkState(srfDir.exists() || srfDir.mkdir());
		String prefix = "event_"+eventID+"_"+(float)dt+"s_"+interp.name();
		if (timeScale != 1d) {
			prefix += "_timeScale"+(float)timeScale;
			if (scaleVelocities)
				prefix += "_velScale";
		}
		return new File(srfDir, prefix+"_bbp");
	}
	
	static final boolean U3_SURFACES = true;
	static VelocityModel VM = VelocityModel.LA_BASIN_500;
	static final Method METHOD = Method.GP;
	static final SRFInterpolationMode SRF_INTERP_MODE = SRFInterpolationMode.ADJ_VEL;
	static final double SRF_DT = 0.05;
	static final double SRF_VERSION = 1.0;
	static final int DEFAULT_SEED = 12345;
	public static final double MIN_SUB_SECT_FRACT = 0.2;
	static final boolean ADJ_WIDTH_MATCH_AREA = true;
	static final double MAX_ADJ_WIDTH = 30;
	static final boolean DO_HF = false;
	
	static final double MAX_DIST = 200d;
	
	static final double SITE_LO_PASS_FREQ = 0.15;
	static final double SITE_HI_PASS_FREQ = 100;
	
	private static final List<BBP_Site> allSites;
	private static final List<BBP_Site> csLAInitialSites;
	static {
		List<BBP_Site> sites = new ArrayList<>();
		
		sites.add(new BBP_Site("USC", new Location(34.0192, -118.286), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("SBSM", new Location(34.064986, -117.29201), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		
		for (NEHRP_TestCity site : NEHRP_TestCity.getShortListCA()) {
			String name = site.name();
			if (name.length() > 10)
				name = name.substring(0, 10);
			sites.add(new BBP_Site(name, site.location(), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		}
		
		allSites = Collections.unmodifiableList(sites);
		
		sites = new ArrayList<>();
		
		sites.add(new BBP_Site("USC", new Location(34.0192, -118.286), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("PAS", new Location(34.148426, -118.17119), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("SBSM", new Location(34.064986, -117.29201), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("WNGC", new Location(34.041824, -118.0653), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("STNI", new Location(33.93088, -118.17881), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("LAPD", new Location(34.557, -118.125), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
//		sites.add(new BBP_Site("s119", new Location(34.55314, -118.72826), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
//		sites.add(new BBP_Site("s279", new Location(34.37809, -118.34757), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
//		sites.add(new BBP_Site("s480", new Location(34.15755, -117.87389), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		sites.add(new BBP_Site("SMCA", new Location(34.00909, -118.48939), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		
		csLAInitialSites = Collections.unmodifiableList(sites);
	}
	
	public static List<BBP_Site> getStandardSites() {
		return allSites;
	}
	
	public static List<BBP_Site> getCyberShakeInitialLASites() {
		return csLAInitialSites;
	}
	
	public static List<BBP_Site> getCAGriddedSites(double spacing) {
		return getGriddedSites(new CaliforniaRegions.RELM_TESTING(), spacing);
	}
	
	public static List<BBP_Site> getSoCalGriddedSites(double spacing) {
		return getGriddedSites(new CaliforniaRegions.RELM_SOCAL(), spacing);
	}
	
	public static List<BBP_Site> getGriddedSites(Region region, double spacing) {
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		System.out.println("Creating "+gridReg.getNodeCount()+" gridded BBP sites");
		
		List<BBP_Site> sites = new ArrayList<>();
		
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			String name = "grid"+i;
			Preconditions.checkState(name.length() <= 10);
			sites.add(new BBP_Site(name, gridReg.getLocation(i), VM.getVs30(), SITE_LO_PASS_FREQ, SITE_HI_PASS_FREQ));
		}
		
		return sites;
	}
	
	public static String siteCleanName(BBP_Site site) {
		for (NEHRP_TestCity city : NEHRP_TestCity.values()) {
			if (city.name().startsWith(site.getName()))
				return city.toString();
		}
		return site.getName();
	}
	
	public static List<BBP_Site> getStandardSites(BBP_SourceFile source) {
		List<BBP_Site> sites = new ArrayList<>();
		
		QuadSurface surf = source.getSurface().getQuadSurface();
		for (BBP_Site site : allSites)
			if (surf.getDistanceJB(site.getLoc()) <= MAX_DIST)
				sites.add(site);
		
		return sites;
	}
	
	public static List<BBP_Site> getStandardSites(File bbpRunDir) {
		List<BBP_Site> sites = new ArrayList<>();
		
		for (BBP_Site site : allSites) {
			for (File file : bbpRunDir.listFiles()) {
				if (file.getName().endsWith(".rd50") && file.getName().contains(site.getName())) {
					sites.add(site);
					break;
				}
			}
		}
		
		return sites;
	}
	
	public static BBP_SourceFile generateBBP_Inputs(RSQSimCatalog catalog, RSQSimEvent event, boolean overwrite)
			throws IOException {
		return generateBBP_Inputs(catalog, event, overwrite, 1d, false);
	}
	
	public static BBP_SourceFile generateBBP_Inputs(RSQSimCatalog catalog, RSQSimEvent event, boolean overwrite,
			double timeScale, boolean scaleVelocities) throws IOException {
		File srcFile = getEventSrcFile(catalog, event.getID());
		
		BBP_SourceFile bbpSource;
		if (!srcFile.exists() || overwrite) {
			BBP_PlanarSurface surface;
			if (U3_SURFACES)
				surface = planarEquivalentU3Surface(catalog, event, MIN_SUB_SECT_FRACT, ADJ_WIDTH_MATCH_AREA);
			else
				surface = estimateBBP_PlanarSurface(event);
			bbpSource = buildBBP_Source(event, surface, DEFAULT_SEED);
			bbpSource.writeToFile(srcFile);
		} else {
			bbpSource = BBP_SourceFile.readFile(srcFile);
		}
		
		File srfFile = getEventSRFFile(catalog, event.getID(), SRF_INTERP_MODE, SRF_DT, timeScale, scaleVelocities);
		if (!srfFile.exists() || overwrite) {
			RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
			if (timeScale != 1d)
				func = func.getTimeScaledFunc(timeScale, scaleVelocities);
			System.out.println("Generating SRF for dt="+(float)SRF_DT+", "+SRF_INTERP_MODE);
			List<SRF_PointData> srf = RSQSimSRFGenerator.buildSRF(func, event.getAllElements(), SRF_DT, SRF_INTERP_MODE);
			SRF_PointData.writeSRF(srfFile, srf, SRF_VERSION);
		}
		
		return bbpSource;
	}
	
	public static void runBBP(RSQSimCatalog catalog, RSQSimEvent event) throws IOException {
		runBBP(catalog, event, null);
	}
	
	public static void runBBP(RSQSimCatalog catalog, RSQSimEvent event, List<BBP_Site> sites) throws IOException {
		runBBP(catalog, event, sites, 1d, false);
	}
	
	public static void runBBP(RSQSimCatalog catalog, RSQSimEvent event, List<BBP_Site> sites,
			double timeScale, boolean scaleVelocities) throws IOException {
		BBP_SourceFile source = generateBBP_Inputs(catalog, event, false, timeScale, scaleVelocities);
		if (sites == null)
			sites = getStandardSites(source);
		File bbpOutputDir = getEventBBPDir(catalog, event.getID(), SRF_INTERP_MODE, SRF_DT, timeScale, scaleVelocities);
		Preconditions.checkState(bbpOutputDir.exists() || bbpOutputDir.mkdir());
		
		File sitesFile = new File(bbpOutputDir, "sites.stl");
		BBP_Site.writeToFile(sitesFile, sites);
		
		BBP_Wrapper bbpWrap = new BBP_Wrapper(VM, METHOD, getEventSrcFile(catalog, event.getID()), null,
				getEventSRFFile(catalog, event.getID(), SRF_INTERP_MODE, SRF_DT, timeScale, scaleVelocities), sitesFile, bbpOutputDir);
		bbpWrap.setDoHF(DO_HF);
		bbpWrap.run();
	}
	
	static VelocityModel detectVM(File bbpDir) throws IOException {
		for (VelocityModel vm : VelocityModel.values())
			if (bbpDir.getName().contains(vm.name()))
				return vm;
		// now try sites
		System.out.println("Fallback to VM detetion from sties...");
		BBP_Site site = BBP_Site.readFile(bbpDir).get(0);
		for (VelocityModel vm : VelocityModel.values())
			if ((float)vm.getVs30() == (float)site.getVs30())
				return vm;
		throw new IllegalStateException("Couldn't detect VM for dir "+bbpDir.getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		int[] ids = { 1670183 };
//		int[] ids = { targetID-2, targetID-1, targetID, targetID+1, targetID+2 };
		double timeScale = 1d;
		boolean scaleVelocities = false;
		List<RSQSimEvent> events = catalog.loader().byIDs(ids);
		RSQSimStateTransitionFileReader trans = catalog.getTransitions();
		for (RSQSimEvent event : events) {
			System.out.println("Running for "+event.getID());
			RSQSimStateTransitionFileReader.printTransitions(event, trans.getTransitions(event));
			List<BBP_Site> sites = BBP_Site.readFile(new File("/data/kevin/bbp/parallel/"
					+ "2018_04_13-rundir2585_1myrs-all-m6.5-skipYears5000-noHF-csLASites/sites.stl"));
			runBBP(catalog, event, sites, timeScale, scaleVelocities);
		}
		System.exit(0);
		events = catalog.loader().minMag(6.5).load();
		int firstZero = -1;
		double firstTime = 0d;
		int lastZero = -1;
		double lastTime = 0d;
		int numZero = 0;
		for (RSQSimEvent event : events) {
			Map<Integer, List<RSQSimStateTime>> eventTrans = trans.getTransitions(event);
			int patches = eventTrans.size();
			int stateTimes = 0;
			for (List<RSQSimStateTime> list : eventTrans.values())
				stateTimes += list.size();
			if (stateTimes == 0) {
				if (firstZero < 0) {
					firstZero = event.getID();
					firstTime = event.getTime();
				}
				lastZero = event.getID();
				lastTime = event.getTime();
				numZero++;
				System.out.println("Event "+event.getID()+" (M="+event.getMagnitude()+") at t="+event.getTime()
					+" has "+stateTimes+" transitions on "+patches+" patches");
			}
		}
		System.out.println("First zero: event "+firstZero+", t="+firstTime);
		System.out.println("Last zero: event "+lastZero+", t="+lastTime);
		System.out.println("Num zero: "+numZero);
	}

}
