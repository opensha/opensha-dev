package scratch.kevin.simulators.ruptures.azimuthal;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.rotation.RuptureRotationUtils;

public class RSQSimAzimuthalSiteConfig extends AzimuthalSiteConfig<RSQSimEvent> {
	
	private RSQSimCatalog catalog;
	private List<RSQSimEvent> events;
	private double buffer;
	private double spacing;
	
	private Map<RSQSimEvent, RSQSimEvent> orientedEvents;
	private Map<RSQSimEvent, QuadSurface> orientedQuads;

	public RSQSimAzimuthalSiteConfig(RSQSimCatalog catalog, Scenario scenario, List<RSQSimEvent> events,
			double buffer, double spacing, boolean positiveX) {
		super(scenario, events);
		this.catalog = catalog;
		this.events = events;
		this.buffer = buffer;
		this.spacing = spacing;
		
		orientedEvents = new HashMap<>();
		orientedQuads = new HashMap<>();
		
		MinMaxAveTracker lenTrack = new MinMaxAveTracker();
		
		for (RSQSimEvent event : events) {
			Location centroid = RuptureRotationUtils.calcRuptureCentroid(event);
			RSQSimEvent oriented = RuptureRotationUtils.getInitiallyOriented(catalog, event, centroid);
			orientedEvents.put(event, oriented);
			QuadSurface quad = RuptureRotationUtils.getIdealizedQuadSurfaceRepresentation(oriented, centroid);
			orientedQuads.put(event, quad);
			lenTrack.addValue(quad.getUpperEdge().getTraceLength());
		}
		
		System.out.println("Lenghts: "+lenTrack);
		double maxLength = Math.min(lenTrack.getMax(), lenTrack.getAverage()*1.3);
		maxLength = spacing*Math.ceil(maxLength/spacing);
		System.out.println("Using "+(float)maxLength+" as max for site gridding");
		
		init(spacing, buffer, maxLength, positiveX);
	}
	
	public RSQSimEvent getOrientedEvent(RSQSimEvent rupture) {
		return orientedEvents.get(rupture);
	}

	@Override
	public List<Location> getRuptureSiteLocs(RSQSimEvent rupture) {
		List<Location> locs = new ArrayList<>();
		QuadSurface quad = orientedQuads.get(rupture);
		Location p1 = quad.getFirstLocOnUpperEdge();
		Location p2 = quad.getLastLocOnUpperEdge();
		EvenlyDiscrXYZ_DataSet gc2xyz = getGC2XYZ();
		for (int i=0; i<gc2xyz.size(); i++) {
			Point2D pt = gc2xyz.getPoint(i);
			
			double rx = pt.getX();
			double ry = pt.getY();
			
			Location loc = SimpleGC2SiteLocationCalc.gc2ToLoc(p1, p2, rx, ry);
			locs.add(loc);
		}
		return locs;
	}

	@Override
	public double getHypocenterDAS(RSQSimEvent rupture) {
		Location hypo = RSQSimUtils.getHypocenter(orientedEvents.get(rupture));
		QuadSurface quad = orientedQuads.get(rupture);
		double[] gc2 = SimpleGC2SiteLocationCalc.locToGC2(
				quad.getFirstLocOnUpperEdge(), quad.getLastLocOnUpperEdge(), hypo);
		return gc2[1]; // ry
	}

	@Override
	public double getHypocenterDDW(RSQSimEvent rupture) {
		Location hypo = RSQSimUtils.getHypocenter(orientedEvents.get(rupture));
		QuadSurface quad = orientedQuads.get(rupture);
		double horzDist = LocationUtils.distanceToLineFast(quad.getFirstLocOnUpperEdge(), quad.getLastLocOnUpperEdge(), hypo);
		double vertDist = Math.max(0d, hypo.getDepth()-quad.getAveRupTopDepth());
		return Math.min(getWidth(rupture), Math.sqrt(horzDist*horzDist + vertDist*vertDist));
	}

	@Override
	public double getLength(RSQSimEvent rupture) {
		QuadSurface quad = orientedQuads.get(rupture);
		return quad.getUpperEdge().getTraceLength();
	}

	@Override
	public double getWidth(RSQSimEvent rupture) {
		QuadSurface quad = orientedQuads.get(rupture);
		return quad.getAveWidth();
	}

	@Override
	public double getDip(RSQSimEvent rupture) {
		QuadSurface quad = orientedQuads.get(rupture);
		return quad.getAveDip();
	}

	@Override
	public double getMag(RSQSimEvent rupture) {
		return rupture.getMagnitude();
	}

	@Override
	public int getID(RSQSimEvent rupture) {
		return rupture.getID();
	}

}
