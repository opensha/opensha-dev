package scratch.aftershockStatisticsETAS;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.aftershockStatisticsETAS.griddedInterpGMPE.DistanceInterpolator;
import scratch.aftershockStatisticsETAS.griddedInterpGMPE.DoubleParameterInterpolator;
import scratch.aftershockStatisticsETAS.griddedInterpGMPE.GriddedInterpGMPE_Calc;

public class ETAS_ShakingForecastCalc {
	
	private static boolean D = false;
	private static double magDelta = 0.1;
	private static double[] depths = { 7, 2 }; // depth of <6.5 and >=6.5, respectively
	
	/**
	 * 
	 * @param calcRegion gridded region for which the the hazard curves will be calculated (one at every grid node)
	 * @param rateModel rate of M>=refMag for the input region
	 * @param refMag reference magnitude for the rate model
	 * @param maxMag Mmax at each grid node
	 * @param b G-R b-balue at each grid node
	 * @param gmpe GMPE, should be initialized by calling setParamDefaults() and have site data set to desired default values
	 * @param mechWts focal mechanism distribution at each grid node
	 * @param maxSourceDist maximum source distance to consider in km. must be finite, as interpolater uses this. typical value is 200 km
	 * @param vs30Provider source of Vs30 data, or null for constant Vs30 (in which case it will use the Vs30 value set in the GMPE)
	 * @return array of hazard curves where the ith element of the array is the curve for the ith node in calcRegion
	 * @throws IOException
	 */
	public static DiscretizedFunc[] calcForecast(GriddedRegion calcRegion, GeoDataSet rateModel, double refMag, double maxMag, double b, ScalarIMR gmpe,
			Map<FocalMech, Double> mechWts, double maxSourceDist, SiteData<Double> vs30Provider) throws IOException {
		return calcForecast( calcRegion,  rateModel,  refMag,  maxMag,  b,  gmpe,
				 mechWts,  maxSourceDist, vs30Provider,  false); 
	}
	
	public static DiscretizedFunc[] calcForecast(GriddedRegion calcRegion, GeoDataSet rateModel, double refMag, double maxMag, double b, ScalarIMR gmpe,
			Map<FocalMech, Double> mechWts, double maxSourceDist, SiteData<Double> vs30Provider, boolean prompt) throws IOException {
		
		double durationYears = 1d; // this must be set to one, because the PSHA codes assume rateModel is annual, but we must give it the total number expected.
		
		GriddedForecast erf = new GriddedForecast(rateModel, refMag, maxMag, b, mechWts, depths, durationYears);
		erf.updateForecast();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		
		
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1);
		
		int numMag = (int)((maxMag - refMag)/magDelta + 0.5) + 1;
		
		// I think this is giving too high of values by assuming everything in-grid is at zero distance. better to assume it is at min distance
//		DistanceInterpolator distInterp = new DistanceInterpolator(true, calcRegion.getSpacing()/2, maxSourceDist, 100);
		DistanceInterpolator distInterp = new DistanceInterpolator(false, calcRegion.getSpacing()/2, maxSourceDist, 100);
		
		DoubleParameterInterpolator vs30Interp = null;
		if (vs30Provider != null)
			vs30Interp = new DoubleParameterInterpolator(
				Vs30_Param.NAME, 180, 760, 20); // matches Wald Allen range
		
		GriddedInterpGMPE_Calc calc = new GriddedInterpGMPE_Calc(gmpe, xVals, b, refMag, maxMag, numMag, distInterp, vs30Interp);
		calc.setPromptForLongCalc(prompt);
		
		// this precalculates to set up the interpolators
		if(D) System.out.println("Setting up interpolators...");
		calc.precalc(durationYears, depths, mechWts);
		
		
		// build sites list
		List<Site> sites = new ArrayList<>();
		List<Double> vs30s = null;
		if (vs30Provider != null)
			vs30s = vs30Provider.getValues(calcRegion.getNodeList());
		for (int i=0; i<calcRegion.getNodeCount(); i++) {
			Site site = new Site(calcRegion.locationForIndex(i));
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>) param.clone());
			if (vs30s != null)
				site.getParameter(Double.class, Vs30_Param.NAME).setValue(vs30s.get(i));
			sites.add(site);
		}
		
		// do interpolated calculation
		if(D) System.out.println("Interpolating...");
		DiscretizedFunc[] curves = calc.calc(rateModel, sites);
		
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
		
		private static <E> List<E> wrapInList(E val) {
			List<E> list = new ArrayList<>();
			list.add(val);
			return list;
		}

		public GriddedForecast(GeoDataSet rateModel, double refMag, double maxMag, double b,
				Map<FocalMech, Double> mechWts, double[] depths, double durationYears) {
			this(rateModel, refMag, maxMag, b, wrapInList(mechWts), depths, durationYears);
		}

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
	
//	public static void main(String[] args) throws IOException {
//		GeoDataSet rateModel = GriddedGeoDataSet.loadXYZFile(new File("/tmp/rateMap.txt"), true);
//		double refMag = 5d;
//		double maxMag = 8.5d;
//		double b = 1;
//		
//		ScalarIMR gmpe = AttenRelRef.BSSA_2014.instance(null);
//		gmpe.setParamDefaults();
//		// for PGA
//		gmpe.setIntensityMeasure(PGA_Param.NAME);
//		// for SA
////		gmpe.setIntensityMeasure(SA_Param.NAME);
////		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 1d);
//		// for PGV
////		gmpe.setIntensityMeasure(PGV_Param.NAME);
//		System.out.println(gmpe.getIntensityMeasure().getName());
//		
//		// Vs30 provider, or null for no Vs30
////		WaldAllenGlobalVs30 vs30Provider = null;
//		WaldAllenGlobalVs30 vs30Provider = new WaldAllenGlobalVs30();
//		// active tectonic coefficients
//		vs30Provider.setActiveCoefficients();
//		// stable coefficients
////		vs30Provider.setStableCoefficients();
//		
////		double durationYears = 30d/365d;
//		
//		Map<FocalMech, Double> mechWts = new HashMap<>();
//		mechWts.put(FocalMech.STRIKE_SLIP, 0.5);
//		mechWts.put(FocalMech.NORMAL, 0.25);
//		mechWts.put(FocalMech.REVERSE, 0.25);
//		
//		// use the rate map region/spacing
//		double calcSpacing = 0.05;
//		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(rateModel.getMaxLat(), rateModel.getMaxLon()),
//				new Location(rateModel.getMinLat(), rateModel.getMinLon())), calcSpacing, null);
//		
//		double maxSourceDist = 200d;
//		
//		DiscretizedFunc[] curves = calcForecast(calcRegion, rateModel, refMag, maxMag, b, gmpe, mechWts,
//				 maxSourceDist, vs30Provider);
//		
//		GeoDataSet map = extractMap(calcRegion.getNodeList(), curves, false, 0.1); // IML with 10% prob. change false--> true and 0.1-->desired IML level for prob exceed IML
//		
//		XYZGraphPanel xyzGP = new XYZGraphPanel();
//		
//		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, map.getMaxZ());
//		
//		XYZPlotSpec spec = new XYZPlotSpec(map, cpt, "Spatial Forecast", "Longitude", "Latitude", gmpe.getIntensityMeasure().getName());
//		
//		Range xRange = new Range(calcRegion.getMinGridLon()-0.5*calcRegion.getLonSpacing(),
//				calcRegion.getMaxGridLon()+0.5*calcRegion.getLonSpacing());
//		Range yRange = new Range(calcRegion.getMinGridLat()-0.5*calcRegion.getLatSpacing(),
//				calcRegion.getMaxGridLat()+0.5*calcRegion.getLatSpacing());
//		
//		xyzGP.drawPlot(spec, false, false, xRange, yRange);
//		
//		JFrame frame = new JFrame("");
//		frame.setContentPane(xyzGP);
//		frame.pack();
//		frame.setVisible(true);
//		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//	}

}
