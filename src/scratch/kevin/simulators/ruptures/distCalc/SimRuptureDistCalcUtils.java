package scratch.kevin.simulators.ruptures.distCalc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.faultSurface.cache.SurfaceDistances;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.Vertex;

import com.google.common.base.Preconditions;

public class SimRuptureDistCalcUtils {
	
	public static LocationElementDistanceCache buildSiteLocDistCache(Location siteLoc) {
		return new LocationElementDistanceCache(siteLoc);
	}
	
	public static class LocationElementDistanceCache {
		
		private Location siteLoc;
		private Map<SimulatorElement, Location[]> elemVertexMap;
		private Map<Location, SurfaceDistances> vertexDistsMap;
		private Map<SimulatorElement, SurfaceDistances> elemDistsMap;

		private LocationElementDistanceCache(Location siteLoc) {
			this.siteLoc = siteLoc;
			
			elemVertexMap = new HashMap<>();
			vertexDistsMap = new HashMap<>();
			elemDistsMap = new HashMap<>();
		}
		
		public SurfaceDistances getDistances(SimulatorElement elem) {
			SurfaceDistances dists = elemDistsMap.get(elem);
			if (dists == null) {
				// calculate it
				Location[] locs = elemVertexMap.get(elem);
				if (locs == null) {
					// need to convert to vanilla locations, with equals and hashCode based on location only and not any
					// vertex metadata. also do some rounding, such that almost identical points are mapped to the same location
					Vertex[] verts = elem.getVertices();
					locs = new Location[verts.length];
					for (int i=0; i<locs.length; i++)
						locs[i] = new Location((float)verts[i].getLatitude(),
								(float)verts[i].getLongitude(), (float)verts[i].getDepth());
					elemVertexMap.putIfAbsent(elem, locs);
				}
				double rJB = Double.POSITIVE_INFINITY;
				double rRup = Double.POSITIVE_INFINITY;
				for (Location loc : locs) {
					SurfaceDistances vertDists = vertexDistsMap.get(loc);
					if (vertDists == null) {
						double vertJB = LocationUtils.horzDistanceFast(siteLoc, loc);
						double vertRup = Math.sqrt(vertJB*vertJB + loc.getDepth()*loc.getDepth());
						vertDists = new SurfaceDistances(vertRup, vertJB, Double.NaN);
						vertexDistsMap.putIfAbsent(loc, vertDists);
					}
					rJB = Double.min(rJB, vertDists.getDistanceJB());
					rRup = Double.min(rRup, vertDists.getDistanceRup());
				}
				dists = new SurfaceDistances(rRup, rJB, Double.NaN);
				elemDistsMap.putIfAbsent(elem, dists);
			}
			return dists;
		}
	}
	
	public static enum DistanceType {
		R_JB("Rjb", "R<sub>JB</sub>") {
			@Override
			protected double calc(SurfaceDistances dists) {
				return dists.getDistanceJB();
			}
		},
		R_RUP("Rrup", "R<sub>Rup</sub>") {
			@Override
			protected double calc(SurfaceDistances dists) {
				return dists.getDistanceRup();
			}
		};
		
		public final String displayName;
		public final String htmlName;

		private DistanceType(String name, String htmlName) {
			this.displayName = name;
			this.htmlName = htmlName;
		}
		
		protected abstract double calc(SurfaceDistances dists);
	}
	
	public static enum Scalar {
		MOMENT("Scalar Moment", "Scalar Moment", "Nm", "Nm") {
			@Override
			protected double calc(SimulatorElement element, double slip) {
				return FaultMomentCalc.getMoment(element.getArea(), slip);
			}
		},
		AREA("Area", "Area", "km^2", "km<sup>2</sup>") {
			@Override
			protected double calc(SimulatorElement element, double slip) {
				return element.getArea()*1e-6;
			}
		};
		
		public final String displayName;
		public final String htmlName;
		public final String units;
		public final String htmlUnits;

		private Scalar(String name, String htmlName, String units, String htmlUnits) {
			this.displayName = name;
			this.htmlName = htmlName;
			this.units = units;
			this.htmlUnits = htmlUnits;
		}
		
		protected abstract double calc(SimulatorElement element, double slip);
	}
	
	public static DiscretizedFunc calcDistScalarFunc(SimulatorEvent event, Location siteLoc,
			LocationElementDistanceCache siteLocDistCache, DistanceType distType, Scalar scalar) {
		List<SimulatorElement> elems = event.getAllElements();
		double[] slips = event.getAllElementSlips();
		
		DiscretizedFunc incrementalFunc = new ArbitrarilyDiscretizedFunc();
		
		for (int i=0; i<slips.length; i++) {
			SimulatorElement elem = elems.get(i);
			double scalarVal = scalar.calc(elem, slips[i]);
			double dist = distType.calc(siteLocDistCache.getDistances(elem));
			int xInd = incrementalFunc.getXIndex(dist);
			if (xInd >= 0)
				incrementalFunc.set(xInd, scalarVal + incrementalFunc.getY(xInd));
			else
				incrementalFunc.set(dist, scalarVal);
		}
		
		double cumulativeVal = 0d;
		DiscretizedFunc cumulativeFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<incrementalFunc.size(); i++) {
			cumulativeVal += incrementalFunc.getY(i);
			cumulativeFunc.set(incrementalFunc.getX(i), cumulativeVal);
		}
		
		return cumulativeFunc;
	}
	
//	public static double calcAveStrikeForDistX(SimulatorEvent event) {
//		ArrayList<SimulatorElement> elems = event.getAllElements();
//		double[] slips = event.getAllElementSlips();
//		
//		for (int i=0; i<elems.size(); i++) {
//			SimulatorElement elem = elems.get(i);
//			double moment = FaultMomentCalc.getMoment(elem.getArea(), slips[i]);
//			
//			elem.getFocalMechanism();
//		}
//	}

}
