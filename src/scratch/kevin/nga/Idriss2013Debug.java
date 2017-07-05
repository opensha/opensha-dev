package scratch.kevin.nga;

import java.util.Map;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.rupForecastImpl.FaultRuptureSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.griddedSeis.Point2Vert_FaultPoisSource;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FourPointEvenlyGriddedSurface;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers.Idriss_2014_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.collect.Maps;

public class Idriss2013Debug {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Idriss_2014_Wrapper imr = new Idriss_2014_Wrapper();
		imr.setParamDefaults();
		Parameter<Double> imt = imr.getParameter(PGA_Param.NAME);
		
		Site site = new Site(new Location(34, -120));
		Vs30_Param vsParam = new Vs30_Param();
		vsParam.setValueAsDefault();
		site.addParameter(vsParam);
		Vs30_TypeParam vsTypeParam = new Vs30_TypeParam();
		vsTypeParam.setValueAsDefault();
		site.addParameter(vsTypeParam);
		site.addParameter(new DepthTo1pt0kmPerSecParam());
		site.addParameter(new DepthTo2pt5kmPerSecParam());
		
		Location faultLoc = new Location(34.1, -120.1);
		EvenlyGriddedSurface faultSurf = new FourPointEvenlyGriddedSurface(faultLoc,
				new Location(faultLoc.getLatitude(), faultLoc.getLongitude(), 10),
				new Location(faultLoc.getLatitude()+0.1, faultLoc.getLongitude(), 10),
				new Location(faultLoc.getLatitude()+0.1, faultLoc.getLongitude(), 0));
		
		// create fault rupture
		
		double prob = 0.1;
		double duration = 30d;
		double mag = 7d;
		double rake = 180;

		FaultRuptureSource fltSrc = new FaultRuptureSource(mag, faultSurf, rake, prob, true);
		
		// point source
		BackgroundRupType bgRupType = BackgroundRupType.CROSSHAIR;
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(mag, 1, 0.1);
		double fracStrikeSlip = 1d;
		double fracNormal = 0d;
		double fracReverse = 0d;
		WC1994_MagLengthRelationship magLenRel = new WC1994_MagLengthRelationship();
		double ptSrcCutoff = 6.0;
		double[] DEPTHS = new double[] {5.0, 1.0};
		
		ProbEqkSource pointSource;

		switch (bgRupType) {
		case CROSSHAIR:
			pointSource = new Point2Vert_FaultPoisSource(faultLoc, mfd, magLenRel, duration,
					ptSrcCutoff, fracStrikeSlip, fracNormal,
					fracReverse, true);
			break;
		case FINITE:
			pointSource = new Point2Vert_FaultPoisSource(faultLoc, mfd, magLenRel, duration,
					ptSrcCutoff, fracStrikeSlip, fracNormal,
					fracReverse, false);
			break;
		case POINT:
			Map<FocalMech, Double> mechMap = Maps.newHashMap();
			mechMap.put(FocalMech.STRIKE_SLIP, fracStrikeSlip);
			mechMap.put(FocalMech.REVERSE, fracReverse);
			mechMap.put(FocalMech.NORMAL, fracNormal);
			pointSource = new PointSource13b(faultLoc, mfd, duration, DEPTHS, mechMap);
			break;
		default:
			throw new IllegalStateException("Unknown Background Rup Type: "+bgRupType);
		}
		System.out.println("Point has "+pointSource.getNumRuptures()+" rups");
		
		imr.setAll(fltSrc.getRupture(0), site, imt);
		System.out.println("Fault Based");
		System.out.println("IMR mean: "+Math.exp(imr.getMean())+"\tstd dev: "+(float)imr.getStdDev());
		imr.setAll(pointSource.getRupture(0), site, imt);
		System.out.println("Gridded");
		System.out.println("IMR mean: "+Math.exp(imr.getMean())+"\tstd dev: "+(float)imr.getStdDev());
	}

}
