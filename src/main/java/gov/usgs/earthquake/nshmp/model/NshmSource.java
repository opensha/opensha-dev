package gov.usgs.earthquake.nshmp.model;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import gov.usgs.earthquake.nshmp.Maths;
import gov.usgs.earthquake.nshmp.model.SystemRuptureSet.SystemRupture;

public abstract class NshmSource extends ProbEqkSource {

  final Object delegate;
  private final int id;

  NshmSource(Object delegate, String name, int id) {
    this.delegate = delegate;
    this.id = id;
    this.setName(name);
  }

  public abstract void setDuration(double duration);

  public void setName(String name) {
    this.name = name;
  }

  public int getNSHM_ID() {
    return id;
  }

  public Object delegate() {
    return delegate;
  }

  @Override
  public LocationList getAllSourceLocs() {
    throw new UnsupportedOperationException();
  }

  @Override
  public RuptureSurface getSourceSurface() {
    throw new UnsupportedOperationException();
  }

  public Object getDelegate() {
    return delegate;
  }

  static class Fault extends NshmSource {

    private final List<NshmRupture> ruptures;

    Fault(
        IterableRuptureSet delegate,
        double weight,
        double duration) {

      super(delegate, delegate.name(), delegate.id());
      this.ruptures = new ArrayList<>();
      for (Rupture rupture : delegate) {
        ruptures.add(new NshmRupture(
            rupture.magnitude(),
            rupture.rake(),
            rupture.rate(),
            weight,
            duration,
            new NshmSurface(rupture.surface())));
      }
    }

    @Override
    public double getMinDistance(Site site) {
      return NshmUtil.distance(site, (IterableRuptureSet) delegate);
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
        GridSource delegate,
        double weight,
        double duration) {

      super(delegate, delegate.name(), delegate.id());
      this.weight = weight;
      this.duration = duration;
    }

    @Override
    public double getMinDistance(Site site) {
      return NshmUtil.distance(site, (GridSource) delegate);
    }

    @Override
    public int getNumRuptures() {
      return ((GridSource) delegate).size();
    }

    @Override
    public ProbEqkRupture getRupture(int index) {
      Rupture rupture = ((GridSource) delegate).get(index);
      return new NshmRupture(
          rupture.magnitude(),
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
        SystemRuptureSet ruptureSet,
        SystemRupture delegate,
        double weight,
        double duration,
        List<NshmSurface> surfaces) {

      super(delegate, ruptureSet.name(), ruptureSet.id());
      this.surfaces = surfaces;
      this.rupture = new NshmRupture(
          delegate.magnitude(),
          delegate.rake(),
          delegate.rate(),
          weight,
          duration,
          CompoundSurface.get(surfaces));
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
    public int getNumRuptures() {
      return 1;
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
