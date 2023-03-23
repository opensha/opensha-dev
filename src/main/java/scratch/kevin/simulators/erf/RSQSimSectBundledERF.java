package scratch.kevin.simulators.erf;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.InvalidRangeException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.EnumParameter;
import org.opensha.commons.param.impl.FileParameter;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.XMLUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.SimulatorUtils;
import org.opensha.sha.simulators.utils.SimulatorUtils.SimulatorElementIDComparator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Range;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class RSQSimSectBundledERF extends AbstractERF {
	
	private static final int GEOM_LONG_ZONE = 11;
	private static final char GEOM_LAT_ZONE = 'S';
	
	// inputs
	private List<? extends FaultSection> subSects;
	private double catDurationYears;
	private List<SimulatorElement> elements;
	
	// either generated on load or loaded from mappings file
	private List<RSQSimSectBundledSource> sourceList;
	private Map<Integer, RSQSimProbEqkRup> eventIDtoRupMap;
	
	// params
	private static final String MAPPING_FILE_PARAM_NAME = "Mapping File";
	private FileParameter mappingFileParam;
	
	private static final String GEOM_FILE_PARAM_NAME = "Geometry File";
	private FileParameter geomFileParam;
	
	private static final String FAULT_MODEL_PARAM_NAME = "Fault Model";
	private static final FaultModels FAULT_MODEL_DEFAULT = FaultModels.FM3_1;
	private EnumParameter<FaultModels> faultModelParam;
	
	private static final String DEF_MODEL_PARAM_NAME = "Deformation Model";
	private static final DeformationModels DEF_MODEL_DEFAULT = DeformationModels.GEOLOGIC;
	private EnumParameter<DeformationModels> defModelParam;
	
	static final String RUP_SURF_RESOLUTION_PARAM_NAME = "Rupture Surface Resolution";
	private DoubleParameter rupSurfResParam;
	
	static final String SRF_DT_PARAM_NAME = "SRF Time Step";
	private DoubleParameter srfDT;
	
	static final String SPIN_UP_PARAM_NAME = "Spin-Up Time Skipped";
	private DoubleParameter spinUpYears;
	
	// misc
	private Map<Integer, List<FaultSection>> parentSectMappings;
	
	// caches
	private Map<IDPairing, Double> subSectDistsCache = new HashMap<>();
	private double parentRegBuffer = Double.NaN;
	private Map<Integer, Region> parentBufferedRegionCache = new HashMap<>();
	private Location prevElemDistLoc = null;
	private Map<SimulatorElement, Double> elemSiteDistances = null;

	public RSQSimSectBundledERF(List<SimulatorElement> elements, List<RSQSimEvent> events, RupSetFaultModel fm, RupSetDeformationModel dm,
			List<? extends FaultSection> subSects, double minMag, double minFractForInclusion, double sourceBuffer, double dt,
			double skipYears) {
		double surfRes = 0d;
		for (SimulatorElement elem : elements)
			surfRes += elem.getArea();
		surfRes /= elements.size();
		surfRes = Math.sqrt(surfRes); // from area to spacing
		surfRes *= 1e-3; // m -> km
		init(null, null, fm, dm, subSects, elements, dt, skipYears, surfRes);
		adjustableParams.removeParameter(mappingFileParam);
		adjustableParams.removeParameter(geomFileParam);
		
		if (minMag > 0) {
			events = new ArrayList<>(events);
			events.removeIf(e -> e.getMagnitude() < minMag);
		}
		this.catDurationYears = SimulatorUtils.getSimulationDurationYears(events);
		
		List<RSQSimSubSectEqkRupture> allRuptures = buildRuptures(elements, events, minFractForInclusion);
		
		buildSources(allRuptures, sourceBuffer);
	}

	public RSQSimSectBundledERF(File mappingFile, File geometryFile, RupSetFaultModel fm, RupSetDeformationModel dm,
			List<? extends FaultSection> subSects, List<SimulatorElement> elements) throws IOException {
		Preconditions.checkNotNull(mappingFile);
		Preconditions.checkState(geometryFile != null || elements != null);
		init(mappingFile, geometryFile, fm, dm, subSects, elements, Double.NaN, Double.NaN, Double.NaN);
		adjustableParams.removeParameter(mappingFileParam);
		adjustableParams.removeParameter(geomFileParam);
	}

	private RSQSimSectBundledERF() {
		init(null, null, null, null, null, null, Double.NaN, Double.NaN, Double.NaN);
	}
	
	private void init(File mappingFile, File geometryFile, RupSetFaultModel fm, RupSetDeformationModel dm,
			List<? extends FaultSection> subSects, List<SimulatorElement> elements, double dt,
			double skipYears, double surfRes) {
		mappingFileParam = new FileParameter(MAPPING_FILE_PARAM_NAME, mappingFile);
		mappingFileParam.addParameterChangeListener(this);
		adjustableParams.addParameter(mappingFileParam);
		
		geomFileParam = new FileParameter(GEOM_FILE_PARAM_NAME, geometryFile);
		geomFileParam.addParameterChangeListener(this);
		adjustableParams.addParameter(geomFileParam);
		
		if (fm == null)
			fm = FAULT_MODEL_DEFAULT;
		if (fm instanceof FaultModels) {
			faultModelParam = new EnumParameter<FaultModels>(FAULT_MODEL_PARAM_NAME,
					EnumSet.allOf(FaultModels.class), (FaultModels)fm, null);
			faultModelParam.addParameterChangeListener(this);
			adjustableParams.addParameter(faultModelParam);
		} else {
			System.err.println("FM not yet supported: "+fm);
		}
		
		if (dm == null)
			dm = DEF_MODEL_DEFAULT;
		if (dm instanceof DeformationModels) {
			defModelParam = new EnumParameter<DeformationModels>(DEF_MODEL_PARAM_NAME,
					EnumSet.allOf(DeformationModels.class), (DeformationModels)dm, null);
			defModelParam.addParameterChangeListener(this);
			adjustableParams.addParameter(defModelParam);
		} else {
			System.err.println("DM not yet supported: "+dm);
		}
		
		rupSurfResParam = new DoubleParameter(RUP_SURF_RESOLUTION_PARAM_NAME, 0d, 100d);
		rupSurfResParam.getConstraint().setNullAllowed(true);
		adjustableParams.addParameter(rupSurfResParam);
		if (Double.isFinite(surfRes))
			rupSurfResParam.setValue(surfRes);
		else
			rupSurfResParam.setValue(null);
		
		srfDT = new DoubleParameter(SRF_DT_PARAM_NAME, 0d, 1d);
		srfDT.getConstraint().setNullAllowed(true);
		srfDT.setUnits("Seconds");
		if (Double.isFinite(dt))
			srfDT.setValue(dt);
		else
			srfDT.setValue(null);
		adjustableParams.addParameter(srfDT);
		
		spinUpYears = new DoubleParameter(SPIN_UP_PARAM_NAME, 0d, 1000000d);
		spinUpYears.getConstraint().setNullAllowed(true);
		spinUpYears.setUnits("Years");
		if (Double.isFinite(skipYears))
			spinUpYears.setValue(skipYears);
		else
			spinUpYears.setValue(null);
		adjustableParams.addParameter(spinUpYears);
		
		if (subSects == null)
			loadSubSects();
		else
			setSubSects(subSects);
		
		// can be null
		this.elements = elements;
		if (elements != null) {
			double aveArea = 0d;
			for (SimulatorElement e : elements)
				aveArea += e.getArea()*1e-6;
			aveArea /= elements.size();
			rupSurfResParam.setValue(aveArea);
		}
		
		this.timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		this.timeSpan.setDuration(1d);
	}
	
	private void loadSubSects() {
		FaultModels fm = faultModelParam.getValue();
		DeformationModels dm = defModelParam.getValue();
		System.out.println("Loading subsections for "+fm+" "+dm);
		setSubSects(RSQSimUtils.getUCERF3SubSectsForComparison(fm, dm));
	}
	
	private void loadElements() {
		File geomFile = geomFileParam.getValue();
		Preconditions.checkNotNull(geomFile, "No geometry file given!");
		Preconditions.checkState(geomFile.exists(), "Geometry file doesn't exist: %s", geomFile.getAbsolutePath());
		try {
			elements = RSQSimFileReader.readGeometryFile(geomFile, GEOM_LONG_ZONE, GEOM_LAT_ZONE);
			
			double aveArea = 0d;
			for (SimulatorElement e : elements)
				aveArea += e.getArea()*1e-6;
			aveArea /= elements.size();
			rupSurfResParam.setValue(aveArea);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	private void setSubSects(List<? extends FaultSection> subSects) {
		subSectDistsCache.clear();
		parentRegBuffer = Double.NaN;
		parentBufferedRegionCache.clear();
		
		parentSectMappings = new HashMap<>();
		for (FaultSection sect : subSects) {
			List<FaultSection> sectsForParents = parentSectMappings.get(sect.getParentSectionId());
			if (sectsForParents == null) {
				sectsForParents = new ArrayList<>();
				parentSectMappings.put(sect.getParentSectionId(), sectsForParents);
			}
			sectsForParents.add(sect);
		}
		this.subSects = subSects;
	}
	
	@Override
	public void parameterChange(ParameterChangeEvent event) {
		if (event.getParameter() == mappingFileParam) {
			sourceList = null;
			eventIDtoRupMap = null;
		} else if (event.getParameter() == faultModelParam || event.getParameter() == defModelParam) {
			sourceList = null;
			eventIDtoRupMap = null;
			subSects = null;
		} else if (event.getParameter() == geomFileParam) {
			sourceList = null;
			eventIDtoRupMap = null;
			elements = null;
		}
	}

	private List<RSQSimSubSectEqkRupture> buildRuptures(List<SimulatorElement> elements,
			List<RSQSimEvent> events, double minFractForInclusion) {
		List<RSQSimSubSectEqkRupture> ruptures = new ArrayList<>();
		
		Map<Integer, Double> subSectAreas = null;
		if (minFractForInclusion > 0)
			subSectAreas = RSQSimUtils.calcSubSectAreas(elements, subSects);
		
		System.out.print("Building subsection-based ruptures...");
		RSQSimSubSectionMapper mapper = new RSQSimSubSectionMapper(subSects, elements,
				minFractForInclusion, subSectAreas, subSectDistsCache);
		for (RSQSimEvent event : events) {
			RSQSimSubSectEqkRupture rupture = RSQSimUtils.buildSubSectBasedRupture(mapper, event);
			ruptures.add(rupture);
		}
		System.out.println("DONE");
		
		return ruptures;
	}
	
	private double calcEventProbability() {
		double eventRate = 1d/catDurationYears;
		
		double timeSpanDuration = this.getTimeSpan().getDuration();
		return 1d - Math.exp(-eventRate * timeSpanDuration);
	}
	
	private void buildSources(List<RSQSimSubSectEqkRupture> allRuptures, double sourceBuffer) {
		Map<RSQSimSubSectEqkRupture, List<Integer>> eventParentSectsMap = new HashMap<>();
		
		for (RSQSimSubSectEqkRupture rup : allRuptures) {
			HashSet<Integer> eventParentSects = new HashSet<>();
			
			for (FaultSection sect : rup.getSubSections())
				eventParentSects.add(sect.getParentSectionId());
			
			eventParentSectsMap.put(rup, sortedList(eventParentSects));
		}
		
		// bundle by source
		System.out.println("Bundling "+allRuptures.size()+" ruptures into sources...");
		Map<List<Integer>, List<RSQSimSubSectEqkRupture>> sourceBundles = doBundleRuptures(eventParentSectsMap, allRuptures);
		System.out.println("DONE. "+sourceBundles.size()+" sources");
		
		MinMaxAveTracker sourceBundleTrack = new MinMaxAveTracker();
		// bundle by rupture
		System.out.println("Building sources");
		RuptureMagComparator rupMagComparator = new RuptureMagComparator();
		
		double rupProb = calcEventProbability();
		sourceList = new ArrayList<>();
		eventIDtoRupMap = new HashMap<>();
		
		int numWithCulled = 0;
		int totalNumCulled = 0;
		
		for (List<Integer> sourceKey : sourceBundles.keySet()) {
			List<RSQSimSubSectEqkRupture> allSourceRuptures = sourceBundles.get(sourceKey);
			// sort by mag, increasing
			Collections.sort(allSourceRuptures, rupMagComparator);
			
			List<RSQSimProbEqkRup> probRups = new ArrayList<>();
			HashSet<FaultSection> sourceSects = new HashSet<>();
			for (RSQSimSubSectEqkRupture rup : allSourceRuptures) {
				List<SimulatorElement> rupElems = rup.getEvent().getAllElements();
				if (Double.isFinite(sourceBuffer) && sourceBuffer > 0) {
					int prevNum = rupElems.size();
					rupElems = filterElementsByBuffer(rupElems, sourceKey, sourceBuffer);
					int numOutside =  prevNum - rupElems.size();
					if (numOutside > 0) {
						System.out.println("Culled "+numOutside+"/"+rupElems.size()+" elements outside of buffered region");
						numWithCulled++;
						totalNumCulled += numOutside;
					}
				}
				RSQSimProbEqkRup probRup = new RSQSimProbEqkRup(rup, rupProb, rupElems);
				eventIDtoRupMap.put(probRup.getEventID(), probRup);
				probRups.add(probRup);
				for (FaultSection sect : rup.getSubSections())
					sourceSects.add(sect);
			}
			
			sourceBundleTrack.addValue(probRups.size());
			
			List<FaultSection> sortedSourceSects = getSortedRuptureSections(sourceSects);
			
			sourceList.add(new RSQSimSectBundledSource(sortedSourceSects, probRups));
		}
		System.out.println("Source bundle stats: "+sourceBundleTrack);
		if (numWithCulled > 0) {
			double each = (double)totalNumCulled/(double)numWithCulled;
			System.out.println("Culled from "+numWithCulled+"/"+allRuptures.size()+" events ("
					+totalNumCulled+" elements, "+(float)each+" per culled event)");
		}
		
		// sort sources
		Collections.sort(sourceList, new SourceNameComparator());
	}
	
	private static String buildSourceName(Set<String> parentSectNames) {
		if (parentSectNames.size() == 1)
			return parentSectNames.iterator().next();
		
		// map of prefixes to full names that share the prefix
		Map<String, List<String>> prefixes = new HashMap<>();
		// map of prefixes to suffixes for that prefix
		Map<String, List<String>> suffixes = new HashMap<>();
		for (String name : parentSectNames) {
			String prefix, suffix;
			if (name.contains(",")) {
				prefix = name.substring(0, name.indexOf(",")).trim();
				suffix = name.substring(name.indexOf(",")+1, name.length()).trim();
			} else if (name.contains("(") && name.contains(")")) {
				prefix = name.substring(0, name.indexOf("(")).trim();
				suffix = name.substring(name.indexOf("(")+1, name.length()).trim().replace(")", "");
			} else {
				continue;
			}
			if (!prefixes.containsKey(prefix)) {
				prefixes.put(prefix, new ArrayList<>());
				suffixes.put(prefix, new ArrayList<>());
			}
			prefixes.get(prefix).add(name);
			suffixes.get(prefix).add(suffix);
		}
		
		List<String> bundledNames = new ArrayList<>();
		for (String prefix : prefixes.keySet()) {
			List<String> namesForPrefix = prefixes.get(prefix);
			if (namesForPrefix.size() < 2)
				continue;
			for (String name : namesForPrefix)
				parentSectNames.remove(name);
			
			// lets bundle!
			List<String> mySuffixes = suffixes.get(prefix);
			Collections.sort(mySuffixes);
			bundledNames.add(prefix+" ("+commaSpaceJoin.join(mySuffixes)+")");
		}
		// now add any individual names that are left
		bundledNames.addAll(parentSectNames);
		
		// sort
		Collections.sort(bundledNames);
		
		return semicolonSpaceJoin.join(bundledNames);
	}
	
	private List<FaultSection> getSortedRuptureSections(Set<FaultSection> subSects) {
		Map<Integer, List<FaultSection>> subSectsBundledMap = new HashMap<>();
		
		for (FaultSection sect : subSects) {
			Integer parentID = sect.getParentSectionId();
			List<FaultSection> parentSects = subSectsBundledMap.get(parentID);
			if (parentSects == null) {
				parentSects = new ArrayList<>();
				subSectsBundledMap.put(parentID, parentSects);
			}
			parentSects.add(sect);
		}
		
		List<List<FaultSection>> subSectsBundledList = new ArrayList<>();
		SectIDComparator sectComp = new SectIDComparator();
		for (List<FaultSection> parentSects : subSectsBundledMap.values()) {
			Collections.sort(parentSects, sectComp);
			subSectsBundledList.add(parentSects);
		}
		
		if (subSectsBundledList.size() > 1)
			// sort it
			subSectsBundledList = SimulatorFaultSystemSolution.sortRupture(
				this.subSects, subSectsBundledList, subSectDistsCache);
		
		List<FaultSection> rupSects = new ArrayList<>();
		for (List<FaultSection> sects : subSectsBundledList)
			rupSects.addAll(sects);
		Preconditions.checkState(!rupSects.isEmpty());
		
		return rupSects;
	}
	
	private static RuptureSurface buildSubSectSurface(List<FaultSection> sortedRupSects) {
		if (sortedRupSects == null)
			return null;
		double gridSpacing = 1d;

		List<RuptureSurface> rupSurfs = new ArrayList<>();
		for (FaultSection sect : sortedRupSects)
			rupSurfs.add(sect.getFaultSurface(gridSpacing, false, false));

		RuptureSurface surf;
		if (rupSurfs.size() == 1)
			surf = rupSurfs.get(0);
		else
			surf = new CompoundSurface(rupSurfs);
		return surf;
	}
	
	private static final Joiner commaSpaceJoin = Joiner.on(", ");
	private static final Joiner semicolonSpaceJoin = Joiner.on("; ");
	
	private static List<Integer> sortedList(Collection<Integer> ids) {
		List<Integer> list = new ArrayList<>(ids);
		Collections.sort(list);
		return list;
	}
	
	private Map<List<Integer>, List<RSQSimSubSectEqkRupture>> doBundleRuptures(
			Map<RSQSimSubSectEqkRupture, List<Integer>> eventSectsMap,
			List<RSQSimSubSectEqkRupture> rupturesToBundle) {
		Map<List<Integer>, List<RSQSimSubSectEqkRupture>> bundles = new HashMap<>();
		
		for (RSQSimSubSectEqkRupture rup : rupturesToBundle) {
			List<Integer> rupSects = eventSectsMap.get(rup);
			
			List<RSQSimSubSectEqkRupture> bundle = bundles.get(rupSects);
			if (bundle == null) {
				bundle = new ArrayList<>();
				bundles.put(rupSects, bundle);
			}
			bundle.add(rup);
		}
		
		return bundles;
	}
	
	private class RuptureMagComparator implements Comparator<EqkRupture> {

		@Override
		public int compare(EqkRupture o1, EqkRupture o2) {
			return Double.compare(o1.getMag(), o2.getMag());
		}
		
	}
	
	private class SourceNameComparator implements Comparator<ProbEqkSource> {

		@Override
		public int compare(ProbEqkSource o1, ProbEqkSource o2) {
			return o1.getName().compareTo(o2.getName());
		}
		
	}
	
	private class SectIDComparator implements Comparator<FaultSection> {

		@Override
		public int compare(FaultSection o1, FaultSection o2) {
			return Integer.compare(o1.getSectionId(), o2.getSectionId());
		}
		
	}
	
	public static class RSQSimProbEqkRup extends ProbEqkRupture {

		private int eventID;
		private double timeSeconds;
		private List<? extends FaultSection> subSects;
		
		// this can be filtered to remove elements far from the mapped surface
		private List<SimulatorElement> rupElems;
		
		private Range<Double> elemLatRange;
		private Range<Double> elemLonRange;
		private Range<Double> elemDepthRange;

		public RSQSimProbEqkRup(RSQSimSubSectEqkRupture rup, double probability, List<SimulatorElement> rupElems) {
			super(rup.getMag(), rup.getAveRake(), probability, rup.getRuptureSurface(), rup.getHypocenterLocation());
			
			init(rup.getEventID(), rup.getTime(), rup.getSubSections(), rupElems);
		}

		public RSQSimProbEqkRup(double mag, double rake, double probability, Location hypo, int eventID, double timeSeconds,
				List<FaultSection> subSects, List<SimulatorElement> rupElems) {
			super(mag, rake, probability, buildSubSectSurface(subSects), hypo);
			init(eventID, timeSeconds, subSects, rupElems);
		}
		
		private void init(int eventID, double timeSeconds, List<? extends FaultSection> subSects, List<SimulatorElement> rupElems) {
			this.eventID = eventID;
			Preconditions.checkState(eventID >= 0);
			this.timeSeconds = timeSeconds;
			this.subSects = subSects;
			if (rupElems != null) {
				MinMaxAveTracker latTrack = new MinMaxAveTracker();
				MinMaxAveTracker lonTrack = new MinMaxAveTracker();
				MinMaxAveTracker depthTrack = new MinMaxAveTracker();
				for (SimulatorElement elem : rupElems) {
					for (Vertex v : elem.getVertices()) {
						latTrack.addValue(v.getLatitude());
						lonTrack.addValue(v.getLongitude());
						depthTrack.addValue(v.getDepth());
					}
				}
				Range<Double> elemLatRange = Range.closed(latTrack.getMin(), latTrack.getMax());
				Range<Double> elemLonRange = Range.closed(lonTrack.getMin(), lonTrack.getMax());
				Range<Double> elemDepthRange = Range.closed(depthTrack.getMin(), depthTrack.getMax());
				this.rupElems = rupElems;
				this.elemLatRange = elemLatRange;
				this.elemLonRange = elemLonRange;
				this.elemDepthRange = elemDepthRange;
			}
		}

		public int getEventID() {
			return eventID;
		}
		
		public double getTime() {
			return timeSeconds;
		}
		
		public double getTimeYears() {
			return timeSeconds / SimulatorUtils.SECONDS_PER_YEAR;
		}
		
		public List<? extends FaultSection> getSortedSubSects() {
			return subSects;
		}
		
		public List<SimulatorElement> getElements() {
			return rupElems;
		}

		public Range<Double> getElemLatRange() {
			return elemLatRange;
		}

		public Range<Double> getElemLonRange() {
			return elemLonRange;
		}

		public Range<Double> getElemDepthRange() {
			return elemDepthRange;
		}
		
	}
	
	public synchronized double getElementDistance(Location loc, SimulatorElement elem) {
		if (loc != prevElemDistLoc) {
			prevElemDistLoc = loc;
			elemSiteDistances = new HashMap<>();
		}
		Double dist = elemSiteDistances.get(elem);
		if (dist == null) {
			dist = LocationUtils.horzDistanceFast(elem.getCenterLocation(), loc);
			elemSiteDistances.put(elem, dist);
		}
		return dist;
	}
	
	public class RSQSimSectBundledSource extends ProbEqkSource {
		
		private RuptureSurface sourceSurf;
		private List<FaultSection> sortedSourceSects;
		private List<RSQSimProbEqkRup> ruptures;
		private Set<Integer> parentIDs;
		private HashSet<SimulatorElement> sourceElementsSet;

		public RSQSimSectBundledSource(List<FaultSection> sortedSourceSects, List<RSQSimProbEqkRup> ruptures) {
			HashSet<String> parentNames = new HashSet<>();
			parentIDs = new HashSet<>();
			for (FaultSection sect : sortedSourceSects) {
				parentNames.add(sect.getParentSectionName());
				parentIDs.add(sect.getParentSectionId());
			}
			this.name = buildSourceName(parentNames);
			this.sortedSourceSects = sortedSourceSects;
			this.ruptures = ruptures;
			sourceElementsSet = new HashSet<>();
			for (RSQSimProbEqkRup rup : ruptures)
				sourceElementsSet.addAll(rup.getElements());
		}

		@Override
		public LocationList getAllSourceLocs() {
			return getSourceSurface().getEvenlyDiscritizedListOfLocsOnSurface();
		}

		@Override
		public synchronized RuptureSurface getSourceSurface() {
			if (sourceSurf == null)
				sourceSurf = buildSubSectSurface(sortedSourceSects);
			return sourceSurf;
		}

		@Override
		public double getMinDistance(Site site) {
			double minDist = Double.POSITIVE_INFINITY;
			Location loc = site.getLocation();
			for (SimulatorElement elem : sourceElementsSet)
				minDist = Math.min(minDist, getElementDistance(loc, elem));
			return minDist;
		}

		@Override
		public int getNumRuptures() {
			return ruptures.size();
		}

		@Override
		public RSQSimProbEqkRup getRupture(int nRupture) {
			return ruptures.get(nRupture);
		}

		public List<FaultSection> getSortedSourceSects() {
			return sortedSourceSects;
		}
		
	}
	
	private static Map<Integer, RSQSimEvent> eventIDmap(List<RSQSimEvent> events) {
		Map<Integer, RSQSimEvent> eventIDmap = new HashMap<>();
		for (RSQSimEvent e : events)
			eventIDmap.put(e.getID(), e);
		return eventIDmap;
	}
	
	/**
	 * RSQSim events can have stray elements that are far from the main rupture surface. We define the main rupture surface
	 * as a subsection-based surface with a minimum section area requirement. This method filters out stray elements which are
	 * at least surfaceBuffer away from any of the parent fault sections involved in given source.
	 * @param elements
	 * @param source
	 * @param surfaceBuffer
	 * @return
	 */
	private List<SimulatorElement> filterElementsByBuffer(List<SimulatorElement> elements,
			Collection<Integer> parentIDs, double surfaceBuffer) {
		List<SimulatorElement> filteredElems = new ArrayList<>();
		
		for (SimulatorElement elem : elements)
			if (isWithinBuffer(elem.getCenterLocation(), parentIDs, surfaceBuffer))
				filteredElems.add(elem);
		
		Preconditions.checkState(!filteredElems.isEmpty());
		return filteredElems;
	}
	
	public synchronized boolean isWithinBuffer(Location loc, Collection<Integer> parentIDs, double buffer) {
		if (buffer < 0 || buffer > 1000 || !Double.isFinite(buffer))
			return true;
		if (buffer != parentRegBuffer) {
			parentRegBuffer = buffer;
			parentBufferedRegionCache.clear();
		}
		for (Integer parentID : parentIDs) {
			Region parentReg = parentBufferedRegionCache.get(parentID);
			if (parentReg == null) {
				// need to build it
				double spacing = 1d;
				
				List<RuptureSurface> surfs = new ArrayList<>();
				for(FaultSection fltData : parentSectMappings.get(parentID))
					surfs.add(fltData.getFaultSurface(spacing, false, true));
				RuptureSurface compound;
				if (surfs.size() == 1)
					compound = surfs.get(0);
				else
					compound = new CompoundSurface(surfs);
				LocationList trace = compound.getEvenlyDiscritizedUpperEdge();
				// this can have duplicates in it, remove those
				for (int p=trace.size(); --p>0;) {
					Location p1 = trace.get(p);
					Location p2 = trace.get(p-1);
					double dist = LocationUtils.horzDistanceFast(p1, p2);
					if (dist < 0.5*spacing)
						trace.remove(p);
				}
				
				parentReg = new Region(trace, buffer);
				
				parentBufferedRegionCache.put(parentID, parentReg);
			}
			if (parentReg.contains(loc))
				return true;
		}
		return false;
	}
	
	public void writeRupturePointFiles(File mainDir) throws IOException {
		writeRupturePointFiles(this, mainDir);
	}
	
	static void writeRupturePointFiles(AbstractERF erf, File mainDir) throws IOException {
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			writeRupturePointFiles(erf, mainDir, sourceID, false);
		}
	}

	static void writeRupturePointFiles(AbstractERF erf, File mainDir, int sourceID, boolean overwrite) throws IOException {
		File sourceDir = new File(mainDir, sourceID+"");
		Preconditions.checkState(sourceDir.exists() || sourceDir.mkdir(),
				"Couldn't build dir for source %s: %s", sourceID, sourceDir.getAbsolutePath());
		ProbEqkSource source = erf.getSource(sourceID);
		for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
			File rupDir = new File(sourceDir, rupID+"");
			Preconditions.checkState(rupDir.exists() || rupDir.mkdir(),
					"Couldn't build dir for source %s, rup %s: %s", sourceID, rupID, rupDir.getAbsolutePath());
			File pointFile = new File(rupDir, sourceID+"_"+rupID+".txt");
			if (!overwrite && pointFile.exists() && pointFile.length() > 0) {
				System.out.println("Already wrote point file for "+sourceID+", "+rupID);
				continue;
			}
			RSQSimProbEqkRup rup = (RSQSimProbEqkRup)source.getRupture(rupID);
			List<SimulatorElement> eventElems = rup.rupElems;
			
			FileWriter fw = new FileWriter(pointFile);
			fw.write("Probability = "+(float)rup.getProbability()+"\n");
			fw.write("Magnitude = "+(float)rup.getMag()+"\n");
			
			double aveArea = 0d;
			for (SimulatorElement e : eventElems)
				aveArea += e.getArea()*1e-6;
			aveArea /= eventElems.size();
			
			fw.write("AveArea = "+(float)aveArea+"\n");
			fw.write("NumPoints = "+eventElems.size()+"\n");
			
			fw.write("#   Lat         Lon         Depth      Rake    Dip     Strike"+"\n");
			for (SimulatorElement e : eventElems) {
				Location loc = e.getCenterLocation();
				FocalMechanism mech = e.getFocalMechanism();
				fw.write((float)loc.getLatitude()+"    "+(float)loc.getLongitude()+"    "
						+(float)loc.getDepth()+"    "+(float)mech.getRake()+"    "+(float)mech.getDip()
						+"    "+(float)mech.getStrike()+"\n");
			}
			
			fw.close();
		}
	}
	
	public void writeRuptureSRFs(File mainDir, RSQSimCatalog catalog, List<RSQSimEvent> events,
			double dt, SRFInterpolationMode interpMode, double momentPDiffThreshold) throws IOException {
		// most efficient to do it in event ID order (more cache hits on transitions file)
		
		Map<Integer, IDPairing> eventToSourceRupMap = new HashMap<>();
		for (int sourceID=0; sourceID<sourceList.size(); sourceID++) {
			RSQSimSectBundledSource source = sourceList.get(sourceID);
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				RSQSimProbEqkRup rup = source.getRupture(rupID);
				eventToSourceRupMap.put(rup.getEventID(), new IDPairing(sourceID, rupID));
			}
		}
		
		Map<Integer, Double> patchAreas = new HashMap<>();
		for (SimulatorElement elem : catalog.getElements())
			patchAreas.put(elem.getID(), elem.getArea());
		
		SummaryStatistics momentStats = new SummaryStatistics();
		for (RSQSimEvent event : events) {
			IDPairing ids = eventToSourceRupMap.get(event.getID());
			Preconditions.checkNotNull(ids, "No rupture found for event %s with M=%s. %s ruptures, %s mappings",
					event.getID(), event.getMagnitude(), events.size(), eventToSourceRupMap.size());
			
			int sourceID = ids.getID1();
			int rupID = ids.getID2();
			
			writeRuptureSRF(mainDir, catalog, sourceID, rupID, event, dt, interpMode,
					momentPDiffThreshold, patchAreas, momentStats, false);
		}
		
		System.out.println("Transition file event moment % differences:");
		System.out.println("\tMean: "+(float)momentStats.getMean());
		System.out.println("\tMin: "+(float)momentStats.getMin());
		System.out.println("\tMax: "+(float)momentStats.getMax());
	}
	
	static void writeRuptureSRF(File mainDir, RSQSimCatalog catalog, int sourceID, int rupID, RSQSimEvent event,
			double dt, SRFInterpolationMode interpMode, double momentPDiffThreshold, Map<Integer, Double> patchAreas,
			SummaryStatistics momentStats, boolean overwrite) throws IOException {
		File sourceDir = new File(mainDir, sourceID+"");
		Preconditions.checkState(sourceDir.exists() || sourceDir.mkdir());
		File rupDir = new File(sourceDir, rupID+"");
		Preconditions.checkState(rupDir.exists() || rupDir.mkdir());
		File srfFile = new File(rupDir, sourceID+"_"+rupID+"_event"+event.getID()+".srf");
		if (!overwrite && srfFile.exists() && srfFile.length() > 0) {
			System.out.println("Already wrote SRF for "+sourceID+", "+rupID);
			return;
		}
		
		RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
		
		validateTransitions(event, func, patchAreas, momentPDiffThreshold, momentStats);
		
		List<SimulatorElement> eventElems = event.getAllElements();
		List<SRF_PointData> srf = RSQSimSRFGenerator.buildSRF(func, eventElems, dt, interpMode);
		
		SRF_PointData.writeSRF(srfFile, srf, 1d);
	}
	
	private static void validateTransitions(RSQSimEvent event, RSQSimEventSlipTimeFunc slipTimeFunc,
			Map<Integer, Double> patchAreas, double pDiffThreshold, SummaryStatistics momentStats) {
		double listMoment = 0d;
		for (EventRecord rec : event)
			listMoment += rec.getMoment();
		
		double transMoment = 0d;
		
		int[] elems = event.getAllElementIDs();
		double[] slips = event.getAllElementSlips();
		
		double time = slipTimeFunc.getEndTime();
		for (int i=0; i<elems.length; i++) {
			double transSlip = slipTimeFunc.getCumulativeEventSlip(elems[i], time);
			Preconditions.checkState(Double.isFinite(transSlip), "Bad slip=%s for patch %s in event %s. Expected %s.",
					transSlip, elems[i], event.getID(), slips[i]);
			transMoment += FaultMomentCalc.getMoment(patchAreas.get(elems[i]), transSlip);
		}
		double pDiff = DataUtils.getPercentDiff(transMoment, listMoment);
		double absDiff = Math.abs(transMoment - listMoment);
		if (pDiff > pDiffThreshold) {
			System.out.println("WARNING: moment mismatch for event "+event.getID()+" (M="+(float)event.getMagnitude()+") at t="
					+event.getTime()+" ("+(float)event.getTimeInYears()+" yrs)");
			double listMag = MagUtils.momentToMag(listMoment);
			double transMag = MagUtils.momentToMag(transMoment);
			System.out.println("\tList File Moment:  "+(float)listMoment+", M="+(float)listMag);
			System.out.println("\tTrans File Moment:  "+(float)transMoment+", M="+(float)transMag);
			System.out.println("\tMag diff: "+(float)(transMag - listMag));
			System.out.println("\t% Diff: "+(float)pDiff);
			System.out.println("\tAbs Diff: "+(float)absDiff+", M="+(float)MagUtils.momentToMag(absDiff));
		}
		if (momentStats != null)
			momentStats.addValue(pDiff);
	}
	
	private static int sequential_flag = -999;
	
	public void writeMappingBinaryFile(File mappingFile) throws IOException {
		DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(mappingFile)));
		
		out.writeInt(getNumSources()); // write number of sources
		out.writeDouble(catDurationYears); // write catalog duration in years
		
		SimulatorElementIDComparator elemComp = new SimulatorElementIDComparator();
		
		long sequentialIDsSkipped = 0;
		long numElemIDs = 0;
		
		for (int sourceID=0; sourceID<getNumSources(); sourceID++) {
			RSQSimSectBundledSource source = getSource(sourceID);
			List<FaultSection> sectsList = source.getSortedSourceSects();
			
			out.writeInt(sourceID); // write source ID
			out.writeInt(sectsList.size()); // write number of sections for this source
			
			for (FaultSection sect : sectsList)
				out.writeInt(sect.getSectionId()); // section ID
			
			out.writeInt(source.getNumRuptures()); // num ruptures
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				RSQSimProbEqkRup rup = source.getRupture(rupID);
				out.writeInt(rupID); // rup ID
				out.writeInt(rup.getEventID()); // event ID
				out.writeDouble(rup.getTime()); // event time (s)
				out.writeDouble(rup.getMag()); // mag
				out.writeDouble(rup.getAveRake()); // rake
				Location hypo = rup.getHypocenterLocation();
				out.writeDouble(hypo.getLatitude()); // hypo lat
				out.writeDouble(hypo.getLongitude()); // hypo lon
				out.writeDouble(hypo.getDepth()); // hypo depth
				// elements
				List<SimulatorElement> rupElems = new ArrayList<>(rup.getElements());
				rupElems.sort(elemComp);
				int numElems = rupElems.size();
				out.writeInt(numElems); // number of elements
				int curIndex = 0;
				while (curIndex < numElems) {
					int curID = rupElems.get(curIndex).getID();
					out.writeInt(curID);
					
					// now see if there is group of sequential IDs which we can compact
					int numSequential = 0; // number of sequential IDs following this one
					for (int i=curIndex+1; i<numElems; i++) {
						int testID = rupElems.get(i).getID();
						if ((testID - curID) == (i - curIndex))
							numSequential++;
						else
							break;
					}
					if (numSequential > 2) {
						// worth it to compact. if 2 or fewer, then writing [id0] [id1] [id2] is no more efficient than
						// writing [id0] [sequential_flag] [id2]
						out.writeInt(sequential_flag);
						sequentialIDsSkipped += numSequential - 2;
						curIndex += numSequential;
						curID = rupElems.get(curIndex).getID();
						out.writeInt(curID);
					}
					curIndex++;
				}
				numElemIDs += rupElems.size();
				List<? extends FaultSection> rupSects = rup.getSortedSubSects();
				out.writeInt(rupSects.size()); // number of sections
				for (FaultSection sect : rupSects)
					out.writeInt(sect.getSectionId()); // section ID
			}
		}
		
		double percent = 100d*sequentialIDsSkipped/numElemIDs;
		System.out.println("Saved "+sequentialIDsSkipped+"/"+numElemIDs+" elem ID writes ("+(float)percent+" %)");
		
		out.close();
	}
	
	private SimulatorElement getElement(int id, int elementOffset) {
		SimulatorElement elem = elements.get(id - elementOffset);
		Preconditions.checkState(elem.getID() == id);
		return elem;
	}
	
	private void buildSourcesFromMappingFile(File mappingFile) throws IOException {
		DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(mappingFile)));
		
		int numSources = in.readInt();
		System.out.println("Loading "+numSources+" sources");
		sourceList = new ArrayList<>();
		eventIDtoRupMap = new HashMap<>();
		
		catDurationYears = in.readDouble();
		System.out.println("Have "+catDurationYears+" years");
		double rupProb = calcEventProbability();
		
		int elementOffset = elements.get(0).getID();
		int lastID = elements.get(elements.size()-1).getID();
		Preconditions.checkState(lastID == elementOffset + elements.size() - 1, "Element list not sequential");
		
		for (int sourceID=0; sourceID<numSources; sourceID++) {
//			System.out.println("Source "+sourceID);
			int testSourceID = in.readInt();
			Preconditions.checkState(sourceID == testSourceID, "Expected sourceID=%s, have %s", sourceID, testSourceID);
			
			int numSourceSects = in.readInt();
//			System.out.println("Source has "+numSourceSects+" sections");
			List<FaultSection> sourceSects = new ArrayList<>();
			for (int s=0; s<numSourceSects; s++)
				sourceSects.add(subSects.get(in.readInt()));
			
			int numRups = in.readInt();
//			System.out.println("Source has "+numRups+" ruptures");
			List<RSQSimProbEqkRup> sourceRups = new ArrayList<>(numRups);
			for (int rupID=0; rupID<numRups; rupID++) {
				int testRupID = in.readInt();
				Preconditions.checkState(rupID == testRupID, "Expected rupID=%s, have %s (sourceID=%s)",
						rupID, testRupID, sourceID);
				int eventID = in.readInt();
				double time = in.readDouble();
				double mag = in.readDouble();
				double rake = in.readDouble();
				double lat = in.readDouble();
				double lon = in.readDouble();
				double depth = in.readDouble();
				Location hypo = new Location(lat, lon, depth);
				
				int numElems = in.readInt();
				List<SimulatorElement> rupElems = new ArrayList<>(numElems);
				
				while (rupElems.size() < numElems) {
					int id = in.readInt();
					if (id == sequential_flag) {
						int prevID = rupElems.get(rupElems.size()-1).getID();
						int endID = in.readInt();
						Preconditions.checkState(endID > prevID && endID <= lastID);
						for (int i=prevID+1; i<=endID; i++)
							rupElems.add(getElement(i, elementOffset));
					} else {
						Preconditions.checkState(id >= elementOffset && id <= lastID);
						rupElems.add(getElement(id, elementOffset));
					}
				}
				
				int numRupSects = in.readInt();
				List<FaultSection> rupSects = new ArrayList<>();
				for (int s=0; s<numRupSects; s++)
					rupSects.add(subSects.get(in.readInt()));
				sourceRups.add(new RSQSimProbEqkRup(mag, rake, rupProb, hypo, eventID, time, rupSects, rupElems));
			}
			
			RSQSimSectBundledSource source = new RSQSimSectBundledSource(sourceSects, sourceRups);
			sourceList.add(source);
		}
		
		in.close();
	}
	
	@Override
	public int getNumSources() {
		return sourceList.size();
	}

	@Override
	public RSQSimSectBundledSource getSource(int idx) {
		return sourceList.get(idx);
	}

	@Override
	public RSQSimProbEqkRup getRupture(int iSource, int nRupture) {
		return getSource(iSource).getRupture(nRupture);
	}

	@Override
	public void updateForecast() {
		if (subSects == null)
			loadSubSects();
		if (elements == null)
			loadElements();
		if (sourceList == null) {
			File mappingFile = mappingFileParam.getValue();
			Preconditions.checkNotNull(mappingFile, "Sources not loaded and mapping file null");
			try {
				buildSourcesFromMappingFile(mappingFile);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		double rupProb = calcEventProbability();
		for (RSQSimSectBundledSource source : sourceList)
			for (ProbEqkRupture rup : source)
				rup.setProbability(rupProb);
	}

	@Override
	public String getName() {
		return "RSQSim Sect Bundled ERF";
	}

	private static void validateMappings(RSQSimSectBundledERF erf1, RSQSimSectBundledERF erf2) {
		System.out.println("Validating mappings");
		Preconditions.checkState(erf1.getNumSources() == erf2.getNumSources());
		
		SimulatorElementIDComparator elemComp = new SimulatorElementIDComparator();
		for (int sourceID=0; sourceID<erf1.getNumSources(); sourceID++) {
			RSQSimSectBundledSource source1 = erf1.getSource(sourceID);
			RSQSimSectBundledSource source2 = erf2.getSource(sourceID);
			
			Preconditions.checkState(source1.getNumRuptures() == source2.getNumRuptures());
			for (int rupID=0; rupID<source1.getNumRuptures(); rupID++) {
				RSQSimProbEqkRup rup1 = source1.getRupture(rupID);
				RSQSimProbEqkRup rup2 = source2.getRupture(rupID);
				Preconditions.checkState(rup1.getEventID() == rup2.getEventID());
				Preconditions.checkState(rup1.getTime() == rup2.getTime());
				Preconditions.checkState(rup1.getMag() == rup2.getMag());
				Preconditions.checkState(rup1.getAveRake() == rup2.getAveRake());
				Preconditions.checkState(rup1.getProbability() == rup2.getProbability());
				Preconditions.checkState(rup1.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface().size()
						== rup2.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface().size());
				Preconditions.checkState(rup1.elemLatRange.equals(rup2.elemLatRange));
				Preconditions.checkState(rup1.elemLonRange.equals(rup2.elemLonRange));
				Preconditions.checkState(rup1.elemDepthRange.equals(rup2.elemDepthRange));
				
				List<SimulatorElement> elems1 = rup1.getElements();
				elems1.sort(elemComp);
				List<SimulatorElement> elems2 = rup1.getElements();
				elems2.sort(elemComp);
				Preconditions.checkState(elems1.size() == elems2.size());
				for (int i=0; i<elems1.size(); i++)
					Preconditions.checkState(elems1.get(i).getID() == elems2.get(i).getID());
			}
		}
		System.out.println("Done validating");
	}

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog;
		boolean writePoints, writeSRFs, writeMappings, testReadOnly;
		double maxDuration = 0;
		
		double dt = 0.05;
//		double dt = 0.1;
		double skipYears = 5000;
		
		File compMappings = null;
		if (args.length >= 1 && args.length < 8) {
			File catalogDir = new File(args[0]);
			catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
					null, null, null, FaultModels.FM3_1, DeformationModels.GEOLOGIC); // TODO
			
			if (args.length >= 2)
				writePoints = Boolean.parseBoolean(args[1]);
			else
				writePoints = true;
			if (args.length >= 3)
				writeSRFs = Boolean.parseBoolean(args[2]);
			else
				writeSRFs = true;
			if (args.length >= 4)
				writeMappings = Boolean.parseBoolean(args[3]);
			else
				writeMappings = true;
			if (args.length >= 5)
				testReadOnly = Boolean.parseBoolean(args[4]);
			else
				testReadOnly = false;
			if (args.length >= 6)
				dt = Double.parseDouble(args[5]);
			if (args.length >= 7)
				skipYears = Double.parseDouble(args[6]);
		} else {
			System.out.println("Assuming hardcoded. Otherwise usage is:");
			System.out.println("<catalogDir> [writePoints writeSRFs writeMappings testReadOnly]");
			File baseDir = new File("/data/kevin/simulators/catalogs");
			if (!baseDir.exists()) {
				System.err.println("hardcoded dir doesn't exist: "+baseDir.getAbsolutePath());
				System.exit(2);
			}
			catalog = Catalogs.BRUCE_4860_10X.instance(baseDir);
			writePoints = false;
			writeSRFs = false;
			writeMappings = true;
			testReadOnly = false;
//			compMappings = new File(catalog.getCatalogDir(), "erf_mappings.bin.bak");
//			maxDuration = 10000;
			maxDuration = 0;
		}
		
		RupSetFaultModel fm = catalog.getFaultModel();
		RupSetDeformationModel dm = catalog.getDeformationModel();
		List<? extends FaultSection> subSects = catalog.getSubSects();
		List<SimulatorElement> elements = catalog.getElements();
		
		if (testReadOnly) {
			File mappingFile = new File(catalog.getCatalogDir(), "erf_mappings.bin");
			File geomFile = compMappings == null ? null : catalog.getGeomFile();
			RSQSimSectBundledERF erf = new RSQSimSectBundledERF(mappingFile, geomFile, fm, dm, subSects, elements);
			erf.updateForecast();
			System.out.println("Source 0 has "+erf.getNumRuptures(0)+" ruptures");
			if (compMappings != null) {
				RSQSimSectBundledERF erf2 = new RSQSimSectBundledERF(
						compMappings, geomFile, fm, dm, subSects, elements);
				erf2.updateForecast();
				validateMappings(erf, erf2);
			}
		} else {
			double minMag = 6.5;
			double minFractForInclusion = 0.2;
			
			System.out.println("Getting ready to load "+catalog.getName());
			System.out.println("\twritePoints="+writePoints);
			System.out.println("\twriteSRFs="+writeSRFs);
			System.out.println("\twriteMappings="+writeMappings);
			System.out.println("\tdt="+dt+" s");
			System.out.println("\tskipYears="+skipYears+" yrs");
			
			File catalogDir = catalog.getCatalogDir();
			File csDataDir = new File(catalogDir, "cybershake_inputs_"+(float)dt+"s_skip"+(int)skipYears+"yrs");
			Preconditions.checkState(csDataDir.exists() || csDataDir.mkdir());
			
			System.out.println("\toutputDir="+csDataDir.getAbsolutePath());
			Loader loader = catalog.loader().skipYears(skipYears).minMag(minMag);
			try {
				catalog.getTransitions();
				loader.hasTransitions();
			} catch (IOException e) {
				e.printStackTrace();
				System.out.println("WARNING: not forcing events to have transitions");
			}
			if (maxDuration > 0)
				loader.maxDuration(maxDuration);
			List<RSQSimEvent> events = loader.load();
			
			SRFInterpolationMode interpMode = SRFInterpolationMode.ADJ_VEL;
			double srfPointCullDist = 100;
			double momentPDiffThreshold = 2;
			
			RSQSimSectBundledERF erf = new RSQSimSectBundledERF(catalog.getElements(), events,
					fm, dm , subSects, minMag, minFractForInclusion, srfPointCullDist, dt, skipYears);
			erf.updateForecast();
			
			System.out.println("Source 0 has "+erf.getNumRuptures(0)+" ruptures");
			
			if (writePoints) {
				System.out.println("Writing rupture points");
				erf.writeRupturePointFiles(csDataDir);
			}
			if (writeSRFs) {
				System.out.println("Writing SRFs");
				erf.writeRuptureSRFs(csDataDir, catalog, events, dt, interpMode, momentPDiffThreshold);
			}
			if (writeMappings) {
				File mappingFile = new File(catalogDir, "erf_mappings.bin");
				System.out.println("Writing mappings file");
				erf.writeMappingBinaryFile(mappingFile);
				
				File xmlFile = new File(catalogDir, "erf_params.xml");
				File geomFile = catalog.getGeomFile();
				RSQSimSectBundledERF metadataERF = new RSQSimSectBundledERF();
				metadataERF.setParameter(MAPPING_FILE_PARAM_NAME, mappingFile);
				metadataERF.setParameter(GEOM_FILE_PARAM_NAME, geomFile);
				metadataERF.setParameter(FAULT_MODEL_PARAM_NAME, fm);
				metadataERF.setParameter(DEF_MODEL_PARAM_NAME, dm);
				metadataERF.setParameter(RUP_SURF_RESOLUTION_PARAM_NAME, Math.sqrt(catalog.getAveArea()));
				metadataERF.setParameter(SRF_DT_PARAM_NAME, dt);
				metadataERF.setParameter(SPIN_UP_PARAM_NAME, skipYears);
				Document doc = XMLUtils.createDocumentWithRoot();
				Element root = doc.getRootElement();
				metadataERF.toXMLMetadata(root);
				XMLUtils.writeDocumentToFile(xmlFile, doc);
				
				// now test reading it in
				System.out.println("Building from mappings");
				RSQSimSectBundledERF erf2 = new RSQSimSectBundledERF(mappingFile, null, fm, dm , subSects, elements);
				erf2.updateForecast();
				validateMappings(erf, erf2);
			}
		}
	}

}
