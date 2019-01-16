package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTime;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RotatedRupVariabilityConfig {

	private List<Location> siteLocs;
	private List<RSQSimEvent> ruptures;
	private List<SimulatorElement> elements;
	private Map<Integer, RSQSimEvent> idToOrigMap;
	private List<RotationSpec> rotations;
	
	// caches
	private Map<RotationSpec, RSQSimEvent> rotationCache;
	private Map<RSQSimEvent, Location> centroidCache;
	private Table<RSQSimEvent, Double, RSQSimEvent> centroidRotationCache;
	private Table<Location, Double, Map<RSQSimEvent, RSQSimEvent>> siteDistTranslationCache;
	
	public RotatedRupVariabilityConfig(List<Location> siteLocs, List<RSQSimEvent> ruptures, List<SimulatorElement> elements,
			double[] distances, int numSourceToSiteAz, int numSiteToSourceAz) {
		init(siteLocs, ruptures, elements, buildRotations(siteLocs, ruptures, distances, numSourceToSiteAz, numSiteToSourceAz));
	}
	
	private void init(List<Location> siteLocs, List<RSQSimEvent> ruptures, List<SimulatorElement> elements, List<RotationSpec> rotations) {
		this.siteLocs = siteLocs;
		this.ruptures = ruptures;
		// copy to new list so that we can add to it
		this.elements = new ArrayList<>(elements);
		if (ruptures != null) {
			idToOrigMap = new HashMap<>();
			for (RSQSimEvent rupture : ruptures)
				idToOrigMap.put(rupture.getID(), rupture);
		}
		this.rotations = rotations;
	}
	
	public class RotationSpec {
		public final int index;
		public final Location siteLoc;
		public final int eventID;
		public final double distance;
		public final double sourceToSiteAz;
		public final double siteToSourceAz;
		
		public RotationSpec(int index, Location siteLoc, int eventID, double distance, double sourceToSiteAz,
				double siteToSourceAz) {
			this.index = index;
			this.siteLoc = siteLoc;
			this.eventID = eventID;
			this.distance = distance;
			this.sourceToSiteAz = sourceToSiteAz;
			this.siteToSourceAz = siteToSourceAz;
		}
	}
	
	private List<RotationSpec> buildRotations(List<Location> siteLocs, List<RSQSimEvent> ruptures, double[] distances,
			int numSourceToSiteAz, int numSiteToSourceAz) {
		List<RotationSpec> rotations = new ArrayList<>();
		int index = 0;
		double sourceToSiteDeltaAz = 360/(double)numSourceToSiteAz;
		double siteToSourceDeltaAz = 360/(double)numSiteToSourceAz;
		for (Location siteLoc : siteLocs) {
			for (double distance : distances) {
				for (RSQSimEvent rupture : ruptures) {
					int eventID = rupture.getID();
					for (int nSrc=0; nSrc<numSourceToSiteAz; nSrc++) {
						double sourceToSiteAz = nSrc*sourceToSiteDeltaAz;
						for (int nSite=0; nSite<numSiteToSourceAz; nSite++) {
							double siteToSourceAz = nSite*siteToSourceDeltaAz;
							rotations.add(new RotationSpec(index++, siteLoc, eventID, distance, sourceToSiteAz, siteToSourceAz));;
						}
					}
				}
			}
		}
		return rotations;
	}
	
	public List<RotationSpec> getRotations() {
		return rotations;
	}
	
	public synchronized RSQSimEvent getRotatedRupture(RotationSpec rotation) {
		if (rotationCache == null)
			rotationCache = new HashMap<>();
		if (rotationCache.containsKey(rotation))
			return rotationCache.get(rotation);
		Preconditions.checkNotNull(ruptures);
		if (centroidCache == null)
			centroidCache = new HashMap<>();
		if (centroidRotationCache == null)
			centroidRotationCache = HashBasedTable.create();
		if (siteDistTranslationCache == null)
			siteDistTranslationCache = HashBasedTable.create();
		Map<RSQSimEvent, RSQSimEvent> translationCache = siteDistTranslationCache.get(rotation.siteLoc, rotation.distance);
		if (translationCache == null) {
			translationCache = new HashMap<>();
			siteDistTranslationCache.put(rotation.siteLoc, rotation.distance, translationCache);
		}
		
		RSQSimEvent rupture = idToOrigMap.get(rotation.eventID);
		Preconditions.checkNotNull(rupture);
		if (rotation.sourceToSiteAz != 0d) {
			// first rotate the rupture around it's centroid
			RSQSimEvent rotated = centroidRotationCache.get(rupture, rotation.sourceToSiteAz);
			if (rotated == null) {
				// not yet cached
				Location centroid = centroidCache.get(rupture);
				if (centroid == null) {
					centroid = RuptureRotationUtils.calcRuptureCentroid(rupture);
					centroidCache.put(rupture, centroid);
				}
				rotated = RuptureRotationUtils.getRotated(rupture, centroid, rotation.sourceToSiteAz, elements);
				centroidRotationCache.put(rupture, rotation.sourceToSiteAz, rotated);
			}
			rupture = rotated;
		}
		Preconditions.checkNotNull(rupture);
		
		// now translate it to the supplied distance
		RSQSimEvent translated = translationCache.get(rupture);
		if (translated == null) {
			// not yet cached, have to do it
			translated = rupture;
			
			double origDist = -1;
			
			int numTrans = 0;
			double minDist = Double.NaN, pDiff = Double.NaN, absDiff = Double.NaN;
			while (numTrans < 5) { // max 5 translations
				Location closest = null;
				minDist = Double.POSITIVE_INFINITY;
				for (SimulatorElement elem : translated.getAllElements()) {
					for (Vertex v : elem.getVertices()) {
						double elemDist = LocationUtils.horzDistanceFast(rotation.siteLoc, v);
						if (elemDist < minDist) {
							minDist = elemDist;
							closest = v;
						}
					}
				}
				if (origDist < 0)
					origDist = minDist;
				
				pDiff = DataUtils.getPercentDiff(minDist, rotation.distance);
				absDiff = Math.abs(minDist - rotation.distance);
				if (pDiff < 0.5 || absDiff < 0.5)
					break;
				
				LocationVector rupToOrigin = LocationUtils.vector(closest, rotation.siteLoc);
				LocationVector transVector = new LocationVector(rupToOrigin.getAzimuth(),
						rupToOrigin.getHorzDistance()-rotation.distance, 0d);
				translated = RuptureRotationUtils.getTranslated(translated, transVector, elements);
				numTrans++;
			}
			
			Preconditions.checkState(pDiff < 0.5 || absDiff < 0.5,
					"Translation didn't work after %s rounds! target: %s, actual: %s, orig: %s",
					numTrans, rotation.distance, minDist, origDist);
			translationCache.put(rupture, translated);
			rupture = translated;
		}
		
		if (rotation.siteToSourceAz != 0d)
			// rotate it around the site
			rupture = RuptureRotationUtils.getRotated(rupture, rotation.siteLoc, rotation.siteToSourceAz, elements);
		rotationCache.put(rotation, rupture);
		
		return rupture;
	}
	
	public List<SimulatorElement> getElements() {
		return elements;
	}
	
	public synchronized RSQSimEventSlipTimeFunc getRotatedSlipTimeFunc(RSQSimCatalog catalog, RSQSimEvent event) throws IOException {
		RSQSimEvent origEvent = idToOrigMap.get(event.getID());
		Map<Integer, List<RSQSimStateTime>> origTransMap = catalog.getTransitions().getTransitions(event);
		Map<Integer, Double> origSlipVels = catalog.getSlipVelocities();
		
		Map<Integer, List<RSQSimStateTime>> patchTransitionsMap = new HashMap<>();
		Map<Integer, Double> slipVels = new HashMap<>();
		boolean variableSlipSpeed = catalog.isVariableSlipSpeed();
		
		int[] origPatchIDs = origEvent.getAllElementIDs();
		int[] rotPatchIDs = origEvent.getAllElementIDs();
		Preconditions.checkState(origPatchIDs.length == rotPatchIDs.length, "Length mistmatch! %s != %s",
				origPatchIDs.length, rotPatchIDs.length);
		Preconditions.checkState(origPatchIDs.length > 0, "No patch IDs?");
		for (int i=0; i<origPatchIDs.length; i++) {
			int origID = origPatchIDs[i];
			int rotID = rotPatchIDs[i];
			List<RSQSimStateTime> origTrans = origTransMap.get(origID);
			List<RSQSimStateTime> rotTrans = new ArrayList<>();
			for (RSQSimStateTime trans : origTrans)
				rotTrans.add(new RSQSimStateTime(rotID, trans.getStartTime(), trans.getEndTime(), trans.getState(), trans.getVelocity()));
			patchTransitionsMap.put(rotID, rotTrans);
			slipVels.put(rotID, origSlipVels.get(origID));
		}
		
		return new RSQSimEventSlipTimeFunc(patchTransitionsMap, slipVels, variableSlipSpeed);
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		int skipYears = 2000;
		
		List<BBP_Site> bbpSites = RSQSimBBP_Config.getCyberShakeInitialLASites();
		
		double[] distances = BBP_PartBValidationConfig.DISTANCES;
		int numSourceToSiteAz = 10;
		int numSiteToSourceAz = 5;
		
		List<Location> siteLocs = new ArrayList<>();
		for (BBP_Site bbpSite : bbpSites)
			siteLocs.add(bbpSite.getLoc());
		
		System.out.println("Loading ruptures for scenario");
		List<RSQSimEvent> ruptures = BBP_PartBValidationConfig.Scenario.M6p6_VERT_SS_SURFACE.getMatches(catalog, skipYears);
		System.out.println("Loaded "+ruptures.size()+" ruptures");
		
		RotatedRupVariabilityConfig config = new RotatedRupVariabilityConfig(
				siteLocs, ruptures, catalog.getElements(), distances, numSourceToSiteAz, numSiteToSourceAz);
		
		System.out.println("Generated "+config.getRotations().size()+" rotations");
		
		System.out.println("Rotating ruptures. Starting element count: "+catalog.getElements().size());
		for (RotationSpec rotation : config.getRotations()) {
			config.getRotatedRupture(rotation);
		}
		System.out.println("Done rotating ruptures. Final element count: "+config.elements.size());
	}

}
