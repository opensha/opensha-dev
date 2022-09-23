package gov.usgs.earthquake.nshmp.model;

import java.awt.geom.Point2D;
import java.nio.file.Path;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers.ASK_2014_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

public class NshmErfTest {

  // private static final Path MODEL =
  // Path.of("../nshmp-lib/src/test/resources/model/test-model");
  private static final Path MODEL = Path.of("../nshm-conus-2018-5.x-maint");

  public static void main(String[] args) {

    NshmErf erf = new NshmErf(MODEL, false, false);
    System.out.println("NSHM ERF size: " + erf.getNumSources());
    erf.getTimeSpan().setDuration(50.0);
    erf.updateForecast();

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
