package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ReturnPeriodUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;

public class HazardCurveMagDependence {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "hazard_curve_mags";
		
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(Presets.BOTH_FM_BRANCH_AVG);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		Site site = new Site(new Location(34.0192, -118.286)); // USC
		
		double[] periods = { 0d, 1d, 2d, 3d };
		
//		EvenlyDiscretizedFunc magBins = new EvenlyDiscretizedFunc(5.1d, 21, 0.2);
		EvenlyDiscretizedFunc magBins = new EvenlyDiscretizedFunc(5.25d, 7, 0.5);
		
		OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
		SiteTranslator trans = new SiteTranslator();
		ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());
		
		for (Parameter<?> param : gmpe.getSiteParams()) {
			trans.setParameterValue(param, datas);
			site.addParameter(param);
		}
		
		double minX = 1e-3;
		double maxX = 1e1;
		
		DiscretizedFunc logXVals = new EvenlyDiscretizedFunc(Math.log(minX), Math.log(maxX), 100);
		DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logXVals)
			xVals.set(Math.exp(pt.getX()), 1d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		cpt = cpt.rescale(0d, magBins.size()-1d);
		
		for (double period : periods) {
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			String xAxisLAbel;
			String myPrefix;
			if (period > 0d) {
				gmpe.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
				xAxisLAbel = (float)period+"s SA (g)";
				myPrefix = prefix+"_"+(float)period+"s";
			} else {
				gmpe.setIntensityMeasure(PGA_Param.NAME);
				xAxisLAbel = "PGA (g)";
				myPrefix = prefix+"_pga";
			}
			
			System.out.println("Doing "+xAxisLAbel);
			
			System.out.println("Calculating full curve");
			calc.getHazardCurve(logXVals, site, gmpe, erf);
			
			DiscretizedFunc curve = xVals.deepClone();
			for (int j=0; j<curve.size(); j++)
				curve.set(j, logXVals.getY(j));
			
			curve.setName("Full Model");
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			for (int m=0; m<magBins.size(); m++) {
				double center = magBins.getX(m);
				double minMag = center - 0.5*magBins.getDelta();
				double maxMag = center + 0.5*magBins.getDelta();
				
				MagBinERF magERF = new MagBinERF(erf, minMag, maxMag);
				System.out.println("Calculating "+magERF.getName()+" curve");
				
				calc.getHazardCurve(logXVals, site, gmpe, magERF);
				
				DiscretizedFunc magCurve = xVals.deepClone();
				for (int j=0; j<magCurve.size(); j++)
					magCurve.set(j, logXVals.getY(j));
				
				magCurve.setName(magERF.getName());
				funcs.add(magCurve);
				Color color = cpt.getColor((float)m).darker();
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Hazard Curve Mag Dependence", xAxisLAbel,
					"Annual Probability of Exceedance");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(24);
			gp.setLegendFontSize(20);
			gp.setBackgroundColor(Color.WHITE);
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			
			Range xRange = new Range(minX, maxX);
			Range yRange = new Range(1e-8, 1e0);
			
			gp.drawGraphPanel(spec, true, true, xRange, yRange);
			
			File file = new File(outputDir, myPrefix);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		}
	}
	
	private static class MagBinERF extends AbstractERF {
		
		private AbstractERF erf;
		private double minMag;
		private double maxMag;

		public MagBinERF(AbstractERF erf, double minMag, double maxMag) {
			this.erf = erf;
			this.minMag = minMag;
			this.maxMag = maxMag;
			
		}

		@Override
		public int getNumSources() {
			return erf.getNumSources();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return new SubsetSource(erf.getSource(idx), minMag, maxMag);
		}

		@Override
		public void updateForecast() {
			erf.updateForecast();
		}

		@Override
		public String getName() {
			return "M∈["+(float)minMag+","+(float)maxMag+")";
		}
		
	}
	
	private static class SubsetSource extends ProbEqkSource {
		
		private ProbEqkSource source;
		private List<Integer> includedIndexes = new ArrayList<>();

		public SubsetSource(ProbEqkSource source, double minMag, double maxMag) {
			this.source = source;
			for (int i=0; i<source.getNumRuptures(); i++) {
				double mag = source.getRupture(i).getMag();
				if (mag >= minMag && mag < maxMag)
					includedIndexes.add(i);
			}
		}

		@Override
		public LocationList getAllSourceLocs() {
			return source.getAllSourceLocs();
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return source.getSourceSurface();
		}

		@Override
		public double getMinDistance(Site site) {
			return source.getMinDistance(site);
		}

		@Override
		public int getNumRuptures() {
			return includedIndexes.size();
		}

		@Override
		public ProbEqkRupture getRupture(int nRupture) {
			return source.getRupture(includedIndexes.get(nRupture));
		}
	}

}
