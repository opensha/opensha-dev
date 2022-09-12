package gov.usgs.earthquake.nshmp.model;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.Optional;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import gov.usgs.earthquake.nshmp.model.SystemRuptureSet.SystemSource;

class NshmSource extends ProbEqkSource {

  private final Source delegate;
  private final List<NshmRupture> ruptures;
  private final Optional<List<NshmSurface>> systemSurfaces;

  NshmSource(
      Source delegate,
      double weight, double duration) {

    this.delegate = delegate;
    this.ruptures = delegate.stream()
        .map(rupture -> new NshmRupture(
            rupture.mag(),
            rupture.rake(),
            rupture.rate(),
            weight,
            duration,
            new NshmSurface(rupture.surface())))
        .collect(toList());
    this.systemSurfaces = Optional.empty();
  }

  NshmSource(
      SystemSource delegate, List<NshmSurface> surfaces,
      double weight, double duration) {

    this.delegate = delegate;
    this.systemSurfaces = Optional.of(surfaces);
    NshmRupture rupture = new NshmRupture(
        delegate.magnitude(),
        delegate.rake(),
        delegate.rate(),
        weight,
        duration,
        new CompoundSurface(surfaces));
    this.ruptures = List.of(rupture);
  }

  List<NshmRupture> ruptures() {
    return ruptures;
  }

  @Override
  public double getMinDistance(Site site) {
    if (delegate instanceof SystemSource) {
      return systemSurfaces.orElseThrow().stream()
          .map(NshmSurface::centroid)
          .mapToDouble(p -> LocationUtils.horzDistanceFast(site.getLocation(), p))
          .min()
          .orElseThrow();

    }
    return NshmUtil.distance(site, delegate);
  }

  @Override
  public int getNumRuptures() {
    return ruptures.size();
  }

  @Override
  public ProbEqkRupture getRupture(int index) {
    return ruptures.get(index);
  }

  @Override
  public LocationList getAllSourceLocs() {
    throw new UnsupportedOperationException();
  }

  @Override
  public RuptureSurface getSourceSurface() {
    throw new UnsupportedOperationException();
  }

}
