package gov.usgs.earthquake.nshmp.model;

import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;

import gov.usgs.earthquake.nshmp.Maths;

class NshmRupture extends ProbEqkRupture {

  final double rate;
  final double weight;

  NshmRupture(
      double mag,
      double rake,
      double rate,
      double weight,
      double duration,
      RuptureSurface surface) {

    super(
        mag, rake,
        Maths.rateToProbability(rate * weight, duration),
        surface, null);

    this.rate = rate;
    this.weight = weight;
  }
}
