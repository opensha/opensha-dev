package gov.usgs.earthquake.nshmp.model;

import org.opensha.commons.data.Site;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.geo.Location;
import gov.usgs.earthquake.nshmp.geo.LocationList;
import gov.usgs.earthquake.nshmp.geo.Locations;

class NshmUtil {

  static Location siteToNshmLocation(Site site) {
    return Location.create(
        site.getLocation().getLatitude(),
        site.getLocation().getLatitude());
  }

  static org.opensha.commons.geo.Location toOpenShaLocation(Location loc) {
    return new org.opensha.commons.geo.Location(loc.latitude, loc.longitude, loc.depth);
  }

  static org.opensha.commons.geo.LocationList toOpenShaLocationList(LocationList locs) {
    org.opensha.commons.geo.LocationList out = new org.opensha.commons.geo.LocationList();
    locs.forEach(loc -> out.add(toOpenShaLocation(loc)));
    return out;
  }

  static Location fromOpenShaLocation(org.opensha.commons.geo.Location loc) {
    return Location.create(loc.lon, loc.lat, loc.depth);
  }

  static double distance(Site site, RuptureSet ruptureSet) {
    Location p1 = fromOpenShaLocation(site.getLocation());
    Location p2 = ruptureSet.location(p1);
    return Locations.horzDistanceFast(p1, p2);
  }

  static double distance(Site site, GridSource source) {
    Location p1 = fromOpenShaLocation(site.getLocation());
    Location p2 = source.location(p1);
    return Locations.horzDistanceFast(p1, p2);
  }

  static TectonicRegionType tectonicSettingToType(
      TectonicSetting setting, SourceType type) {
    switch (setting) {
      case ACTIVE_CRUST:
        return TectonicRegionType.ACTIVE_SHALLOW;
      case STABLE_CRUST:
        return TectonicRegionType.STABLE_SHALLOW;
      case SUBDUCTION:
        return (type == SourceType.SLAB || type == SourceType.INTRASLAB_GRID)
            ? TectonicRegionType.SUBDUCTION_SLAB
            : TectonicRegionType.SUBDUCTION_INTERFACE;
      case VOLCANIC:
        return TectonicRegionType.VOLCANIC;
      default:
        throw new UnsupportedOperationException();
    }
  }
}
