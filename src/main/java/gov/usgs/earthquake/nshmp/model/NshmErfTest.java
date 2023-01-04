package gov.usgs.earthquake.nshmp.model;

import static org.opensha.sha.util.TectonicRegionType.ACTIVE_SHALLOW;

import java.awt.geom.Point2D;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers.ASK_2014_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.TectonicRegionType;

public class NshmErfTest {

  // private static final Path MODEL =
  // Path.of("../nshmp-lib/src/test/resources/model/test-model");
  // private static final Path MODEL = Path.of("../nshm-conus-2018-5.x-maint");
  private static final Path MODEL = Path.of("../nshm-conus");

  static gov.usgs.earthquake.nshmp.geo.Location testLoc =
      gov.usgs.earthquake.nshmp.geo.Location.create(-122, 39.0);

  // static gov.usgs.earthquake.nshmp.geo.Location testLoc =
  // gov.usgs.earthquake.nshmp.geo.Location.create(-110, 37.5);

  public static void main(String[] args) {

    Set<TectonicRegionType> trts = EnumSet.of(ACTIVE_SHALLOW);
    NshmErf erf = new NshmErf(MODEL, trts, true);
    System.out.println("NSHM ERF size: " + erf.getNumSources());
    erf.getTimeSpan().setDuration(50.0);
    erf.updateForecast();

    for (ProbEqkSource src : erf) {
      // Source nshmSrc = ((NshmSource) src).delegate;
      // if (nshmSrc instanceof PointSourceFinite) {
      // PointSourceFinite ptSrc = (PointSourceFinite) nshmSrc;
      // if (ptSrc.loc.equals(testLoc)) {
      // System.out.println(testLoc);
      // System.out.println(ptSrc.mfd);
      // for (Rupture rup : ptSrc) {
      // PointSourceFinite.FiniteSurface surf =
      // (PointSourceFinite.FiniteSurface) rup.surface();
      // System.out.println(surf.mag + " " + surf.zTor + " " + surf.dip());
      // }
      //
      // for (ProbEqkRupture rup : src) {
      // System.out.println(rup.getRuptureSurface().getAveRupTopDepth());
      // }
      // }
      // }

      for (ProbEqkRupture rup : src) {
        LocationList locs = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
        if (locs.size() < 1) {
          System.out.println("Problem rupture: " + src.getName() + " " + locs.size());
          break;
        }
      }
    }

    // calcHazard(erf);
  }

  private static void calcHazard(NshmErf erf) {
    ScalarIMR gmpe = new ASK_2014_Wrapper();
    gmpe.setParamDefaults();
    gmpe.setIntensityMeasure(PGA_Param.NAME);

    Site site = new Site(new Location(34, -118)); // Los Angeles
    // Site site = new Site(new Location(40.75, -111.90)); // Salt lake City

    for (Parameter<?> param : gmpe.getSiteParams()) {
      site.addParameter((Parameter<?>) param.clone());
    }
    site.getParameter(Double.class, Vs30_Param.NAME).setValue(760.0);

    DiscretizedFunc linearXVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
    DiscretizedFunc hazardCurve = new ArbitrarilyDiscretizedFunc();
    for (Point2D pt : linearXVals) {
      hazardCurve.set(Math.log(pt.getX()), 0d);
    }

    HazardCurveCalculator curveCalc = new HazardCurveCalculator();
    curveCalc.getHazardCurve(hazardCurve, site, gmpe, erf);

    System.out.println("DONE");
    System.out.println(hazardCurve);

  }

}
