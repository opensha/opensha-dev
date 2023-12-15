package gov.usgs.earthquake.nshmp.model;

import static java.util.stream.Collectors.toList;

import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import gov.usgs.earthquake.nshmp.Maths;
import gov.usgs.earthquake.nshmp.model.SystemRuptureSet.SystemSource;

public abstract class NshmSource extends ProbEqkSource {

  final Source delegate;

  NshmSource(Source delegate) {
    this.delegate = delegate;
    this.setName(delegate.name());
  }

  public abstract void setDuration(double duration);
  
  public void setName(String name) {
	  this.name = name;
  }
  
  public int getNSHM_ID() {
	  return delegate.id();
  }

  public Source delegate() {
    return delegate;
  }

  @Override
  public double getMinDistance(Site site) {
    return NshmUtil.distance(site, delegate);
  }

  @Override
  public int getNumRuptures() {
    return delegate.size();
  }

  @Override
  public LocationList getAllSourceLocs() {
    throw new UnsupportedOperationException();
  }

  @Override
  public RuptureSurface getSourceSurface() {
    throw new UnsupportedOperationException();
  }
  
  public Source getDelegate() {
	  return delegate;
  }

  static class Fault extends NshmSource {

    private final List<NshmRupture> ruptures;

    Fault(
        Source delegate,
        double weight,
        double duration) {

      super(delegate);
      this.ruptures = delegate.stream()
          .map(rupture -> new NshmRupture(
              rupture.mag(),
              rupture.rake(),
              rupture.rate(),
              weight,
              duration,
              new NshmSurface(rupture.surface())))
          .collect(toList());
    }

    @Override
    public ProbEqkRupture getRupture(int index) {
      return ruptures.get(index);
    }

    @Override
    public void setDuration(double duration) {
      ruptures.forEach(rup -> rup.setProbability(
          Maths.rateToProbability(
              rup.rate * rup.weight,
              duration)));
    }
  }

  static class Point extends NshmSource {

    final double weight;
    double duration;

    Point(
        PointSource delegate,
        double weight,
        double duration) {

      super(delegate);
      this.weight = weight;
      this.duration = duration;
    }

    @Override
    public ProbEqkRupture getRupture(int index) {
      Rupture rupture = ((PointSource) delegate).get(index);
      return new NshmRupture(
          rupture.mag(),
          rupture.rake(),
          rupture.rate(),
          weight,
          duration,
          new NshmSurface(rupture.surface()));
    }

    @Override
    public void setDuration(double duration) {
      this.duration = duration;
    }
  }

  static class System extends NshmSource {

    final List<NshmSurface> surfaces;
    final NshmRupture rupture;

    System(
        SystemSource delegate,
        double weight,
        double duration,
        List<NshmSurface> surfaces) {

      super(delegate);
      this.surfaces = surfaces;
      this.rupture = new NshmRupture(
          delegate.magnitude(),
          delegate.rake(),
          delegate.rate(),
          weight,
          duration,
          new CompoundSurface(surfaces));
    }

    @Override
    public double getMinDistance(Site site) {
      return surfaces.stream()
          .map(NshmSurface::centroid)
          .mapToDouble(p -> LocationUtils.horzDistanceFast(site.getLocation(), p))
          .min()
          .orElseThrow();
    }

    @Override
    public ProbEqkRupture getRupture(int index) {
      return rupture;
    }

    @Override
    public void setDuration(double duration) {
      double p = Maths.rateToProbability(
          rupture.rate * rupture.weight,
          duration);
      rupture.setProbability(p);
    }
  }

}
