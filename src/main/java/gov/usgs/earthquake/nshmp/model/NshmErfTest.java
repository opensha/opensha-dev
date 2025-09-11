package gov.usgs.earthquake.nshmp.model;

import java.awt.geom.Point2D;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.Optional;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.RuptureSurface;
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
  private static final Path MODEL = Path.of("../nshm-prvi");

   static gov.usgs.earthquake.nshmp.geo.Location testLoc =
       gov.usgs.earthquake.nshmp.geo.Location.create(-66.117, 18.465);
//  static gov.usgs.earthquake.nshmp.geo.Location testLoc =
//      gov.usgs.earthquake.nshmp.geo.Location.create(-80, 33.2);

  // static gov.usgs.earthquake.nshmp.geo.Location testLoc =
  // gov.usgs.earthquake.nshmp.geo.Location.create(-110, 37.5);

  public static void main(String[] args) {

//    Set<TectonicRegionType> trts =
//        EnumSet.of(TectonicRegionType.STABLE_SHALLOW);
     Set<TectonicRegionType> trts = EnumSet.noneOf(TectonicRegionType.class);

    HazardModel model = HazardModel.load(MODEL);
    NshmErf erf = new NshmErf(model, trts, IncludeBackgroundOption.INCLUDE);
    System.out.println("NSHM ERF size: " + erf.getNumSources());
    erf.getTimeSpan().setDuration(1.0);
    erf.updateForecast();

    System.out.println(Models.mfd(
        model,
        TectonicSetting.SUBDUCTION,
        Optional.of(SourceType.INTERFACE)));

    // for (ProbEqkSource src : erf) {
    //
    // Source nshmSrc = ((NshmSource) src).delegate();
    // System.out.println(nshmSrc.id() + " " + nshmSrc.name());
    //
    // if (nshmSrc instanceof PointSourceFixedStrike) {
    //
    // PointSourceFixedStrike ptSrc = (PointSourceFixedStrike) nshmSrc;
    // if (ptSrc.loc.equals(testLoc)) {
    // System.out.println(testLoc);
    // System.out.println(ptSrc.mfd);
    // for (Rupture rup : ptSrc) {
    // PointSourceFinite.FiniteSurface surf =
    // (PointSourceFinite.FiniteSurface) rup.surface();
    // System.out.println(surf.mag + " " + surf.zTor + " " + surf.dip());
    // System.out.println(surf);
    // }
    // }
    // }
    // }

    for (ProbEqkSource src : erf) {

      // Source nshmSrc = ((NshmSource) src).delegate();
      NshmSource nshmSrc = (NshmSource) src;
      System.out.println(nshmSrc.delegate().id() + " " + nshmSrc.delegate().name());

      System.out.println(nshmSrc.getClass());
      if (nshmSrc instanceof NshmSource.Point) {
        System.out.println(true);

        NshmSource.Point ptSrc = (NshmSource.Point) nshmSrc;
        PointSource nhPtSrc = (PointSource) nshmSrc.delegate();

        // PointSourceFixedStrike ptSrc = (PointSourceFixedStrike) nshmSrc;

        if (nhPtSrc.loc.equals(testLoc)) {

          System.out.println(testLoc);
          System.out.println(nhPtSrc.mfd);

          for (ProbEqkRupture rup : ptSrc) {
            RuptureSurface surface = rup.getRuptureSurface();

            // PointSourceFinite.FiniteSurface surf =
            // (PointSourceFinite.FiniteSurface) rup.surface();
            System.out.println(
                rup.getMag() + " " +
                    surface.getAveRupTopDepth() + " " +
                    surface.getAveDip());
            System.out.println(surface);
          }
          //
          // for (ProbEqkRupture rup : src) {
          // System.out.println(rup.getRuptureSurface().getAveRupTopDepth());
          // }
          // }
          // }

          // for (ProbEqkRupture rup : src) {
          // LocationList locs =
          // rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
          // if (locs.size() < 1) {
          // System.out.println("Problem rupture: " + src.getName() + " " +
          // locs.size());
          // break;
        }
      }
    }

     calcHazard(erf);
  }

  private static void calcHazard(NshmErf erf) {
    ScalarIMR gmpe = new ASK_2014_Wrapper();
    gmpe.setParamDefaults();
    gmpe.setIntensityMeasure(PGA_Param.NAME);

    Site site = new Site(new Location(18.465, -66.117)); // San Juan
//    Site site = new Site(new Location(34, -118)); // Los Angeles
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
