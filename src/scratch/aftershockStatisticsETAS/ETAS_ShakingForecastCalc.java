package scratch.aftershockStatisticsETAS;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;

public class ETAS_ShakingForecastCalc {
	
	private static double magDelta = 0.1;
	private static double[] depths = { 7, 2 }; // depth of <6.5 and >=6.5, respectively
	
	public static DiscretizedFunc[] calcForecast(GeoDataSet rateModel, double refMag, double maxMag, double b, ScalarIMR gmpe,
			Map<FocalMech, Double> mechWts, double durationYears, SiteData<Double> vs30Provider) {
		// uniform focal mechanism distribution
		List<Map<FocalMech, Double>> mechWtsList = new ArrayList<>();
		mechWtsList.add(mechWts);
		return calcForecast(rateModel, refMag, maxMag, b, gmpe, mechWtsList, durationYears, vs30Provider);
	}
	
	public static DiscretizedFunc[] calcForecast(GeoDataSet rateModel, double refMag, double maxMag, double b, ScalarIMR gmpe,
			List<Map<FocalMech, Double>> mechWts, double durationYears, SiteData<Double> vs30Provider) {
		Preconditions.checkArgument(mechWts.size() == 1 || mechWts.size() == rateModel.size(),
				"Must either have a single forcal mechanism weight set, or one for each node. Expected %s (or 1), got %s",
				rateModel.size(), mechWts.size());
		
		GriddedForecast erf = new GriddedForecast(rateModel, refMag, maxMag, b, mechWts, depths, durationYears);
		erf.updateForecast();
		
		DiscretizedFunc[] curves = new DiscretizedFunc[rateModel.size()];
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1);
		
		SiteDataValueList<Double> vs30Vals = null;
		SiteTranslator siteTrans = null;
		if (vs30Provider != null) {
			siteTrans = new SiteTranslator();
			try {
				vs30Vals = vs30Provider.getAnnotatedValues(rateModel.getLocationList());
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
		for (int i=0; i<rateModel.size(); i++) {
			if (i % 100 == 0)
				System.out.println("Calculating curve for site "+i+"/"+rateModel.size());
			
			Location loc = rateModel.getLocation(i);
			
			Site site = new Site(loc);
			if (vs30Vals != null) {
				// default site params from GMPE
				for (Parameter<?> param : gmpe.getSiteParams())
					site.addParameter((Parameter<?>)param.clone());
				// set Vs30
				SiteDataValue<Double> vs30 = vs30Vals.getValue(i);
				if (vs30Provider.isValueValid(vs30.getValue())) {
					for (Parameter<?> param : site)
						siteTrans.setParameterValue(param, vs30);
				}
			} else {
				// default site params from GMPE
				site.addParameterList(gmpe.getSiteParams());
			}
			
			DiscretizedFunc logCurve = logXVals.deepClone();
			calc.getHazardCurve(logCurve, site, gmpe, erf);
			
			// convert back to linear
			curves[i] = new ArbitrarilyDiscretizedFunc();
			for (int j=0; j<xVals.size(); j++)
				curves[i].set(xVals.getX(j), logCurve.getY(j));
		}
		
		return curves;
	}
	
	public static class GriddedForecast extends AbstractERF {
		
		private GeoDataSet rateModel;
		private double refMag;
		private double maxMag;
		private double b;
		private List<Map<FocalMech, Double>> mechWts;
		private double[] depths;
		private double durationYears;
		
		private int mfdNum;
		
		private List<ProbEqkSource> sources;

		public GriddedForecast(GeoDataSet rateModel, double refMag, double maxMag, double b,
				List<Map<FocalMech, Double>> mechWts, double[] depths, double durationYears) {
			this.rateModel = rateModel;
			this.refMag = refMag;
			this.maxMag = maxMag;
			this.b = b;
			this.mechWts = mechWts;
			this.depths = depths;
			this.durationYears = durationYears;
			
			mfdNum = (int)((maxMag - refMag)/magDelta + 0.5) + 1;
			Preconditions.checkState(mfdNum >= 1, "Ref mag < max mag?");
			
			timeSpan = new TimeSpan(TimeSpan.YEARS, TimeSpan.YEARS);
			
			getTimeSpan().setDuration(durationYears);
		}

		@Override
		public int getNumSources() {
			return sources.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return sources.get(idx);
		}

		@Override
		public void updateForecast() {
			sources = new ArrayList<>();
			
			for (int i=0; i<rateModel.size(); i++) {
				Location loc = rateModel.getLocation(i);
				double rate = rateModel.get(i);
				if (rate == 0)
					continue;
				
				// rate should be total cumulative rate
				GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(b, rate, refMag, maxMag, mfdNum);
				
				Map<FocalMech, Double> mechWtMap;
				if (mechWts.size() == 1)
					mechWtMap = mechWts.get(0);
				else
					mechWtMap = mechWts.get(i);
				
				PointSource13b source = new PointSource13b(loc, mfd, durationYears, depths, mechWtMap);
				sources.add(source);
			}
		}

		@Override
		public String getName() {
			return "Gridded Rate Model Forecast";
		}
		
	}
	
	public static GeoDataSet extractMap(List<Location> locations, DiscretizedFunc[] curves, boolean isProbAt_IML, double level) {
		GeoDataSet geo = new ArbDiscrGeoDataSet(false);
		
		for (Location loc : locations)
			geo.set(loc, 0d);
		populateMap(geo, curves, isProbAt_IML, level);
		
		return geo;
	}
	
	public static GriddedGeoDataSet extractMap(GriddedRegion gridReg, DiscretizedFunc[] curves, boolean isProbAt_IML, double level) {
		GriddedGeoDataSet geo = new GriddedGeoDataSet(gridReg, false);
		
		populateMap(geo, curves, isProbAt_IML, level);
		
		return geo;
	}
	
	private static void populateMap(GeoDataSet geo, DiscretizedFunc[] curves, boolean isProbAt_IML, double level) {
		Preconditions.checkArgument(geo.size() == curves.length);
		
		for (int i=0; i<geo.size(); i++)
			geo.set(i, HazardDataSetLoader.getCurveVal(curves[i], isProbAt_IML, level));
	}
	
	public static void main(String[] args) throws IOException {
		GeoDataSet rateModel = ArbDiscrGeoDataSet.loadXYZFile("/tmp/rateMap.txt", true);
		double refMag = 5d;
		double maxMag = 8.5d;
		double b = 1;
		
		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(PGA_Param.NAME);
		System.out.println(gmpe.getIntensityMeasure().getName());
		
		WaldAllenGlobalVs30 vs30Provider = new WaldAllenGlobalVs30();
		vs30Provider.setActiveCoefficients();
		
		double durationYears = 30d/365d;
		
		Map<FocalMech, Double> mechWts = new HashMap<>();
		mechWts.put(FocalMech.STRIKE_SLIP, 0.5);
		mechWts.put(FocalMech.NORMAL, 0.25);
		mechWts.put(FocalMech.REVERSE, 0.25);
		
		DiscretizedFunc[] curves = calcForecast(rateModel, refMag, maxMag, b, gmpe, mechWts, durationYears, vs30Provider);
		
		GeoDataSet map = extractMap(rateModel.getLocationList(), curves, false, 0.5); // IML with 50% prob
		
		XYZGraphPanel xyzGP = new XYZGraphPanel();
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		
		XYZPlotSpec spec = new XYZPlotSpec(map, cpt, "Spatial Forecast", "Longitude", "Latitude", gmpe.getIntensityMeasure().getName());
		
		xyzGP.drawPlot(spec, false, false, null, null);
		
		JFrame frame = new JFrame("");
		frame.setContentPane(xyzGP);
		frame.setVisible(true);
	}

}
