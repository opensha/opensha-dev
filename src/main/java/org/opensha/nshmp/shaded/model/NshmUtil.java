package org.opensha.nshmp.shaded.model;

import org.opensha.commons.data.Site;
import org.opensha.sha.util.TectonicRegionType;

import org.opensha.nshmp.shaded.geo.NshmpLocation;
import org.opensha.nshmp.shaded.geo.NshmpLocationList;
import org.opensha.nshmp.shaded.geo.NshmpLocations;

class NshmUtil {

	static NshmpLocation siteToNshmLocation(Site site) {
		return NshmpLocation.create(
				site.getLocation().getLatitude(),
				site.getLocation().getLatitude());
	}

	static org.opensha.commons.geo.Location toOpenShaLocation(NshmpLocation loc) {
		return new org.opensha.commons.geo.Location(loc.latitude, loc.longitude, loc.depth);
	}

	static org.opensha.commons.geo.LocationList toOpenShaLocationList(NshmpLocationList locs) {
		org.opensha.commons.geo.LocationList out = new org.opensha.commons.geo.LocationList();
		locs.forEach(loc -> out.add(toOpenShaLocation(loc)));
		return out;
	}

	static NshmpLocation fromOpenShaLocation(org.opensha.commons.geo.Location loc) {
		return NshmpLocation.create(loc.lon, loc.lat, loc.depth);
	}

	static double distance(Site site, NshmpRuptureSet ruptureSet) {
		NshmpLocation p1 = fromOpenShaLocation(site.getLocation());
		NshmpLocation p2 = ruptureSet.location(p1);
		return NshmpLocations.horzDistanceFast(p1, p2);
	}

	static double distance(Site site, NshmpGridSource source) {
		NshmpLocation p1 = fromOpenShaLocation(site.getLocation());
		NshmpLocation p2 = source.location(p1);
		return NshmpLocations.horzDistanceFast(p1, p2);
	}

	static TectonicRegionType tectonicSettingToType(
			NshmpTectonicSetting setting, NshmpSourceType type) {
		switch (setting) {
		case ACTIVE_CRUST:
			return TectonicRegionType.ACTIVE_SHALLOW;
		case STABLE_CRUST:
			return TectonicRegionType.STABLE_SHALLOW;
		case SUBDUCTION:
			return (type == NshmpSourceType.SLAB || type == NshmpSourceType.INTRASLAB_GRID)
					? TectonicRegionType.SUBDUCTION_SLAB
							: TectonicRegionType.SUBDUCTION_INTERFACE;
		case VOLCANIC:
			return TectonicRegionType.VOLCANIC;
		default:
			throw new UnsupportedOperationException();
		}
	}
}
