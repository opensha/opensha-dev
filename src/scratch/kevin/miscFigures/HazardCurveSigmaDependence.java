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
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.mod.ModAttenRelRef;
import org.opensha.sha.imr.mod.ModAttenuationRelationship;
import org.opensha.sha.imr.mod.impl.FixedStdDevMod;
import org.opensha.sha.imr.mod.impl.SimpleScaleMod;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;

public class HazardCurveSigmaDependence {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "sigma_comparison";
		
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(Presets.BOTH_FM_BRANCH_AVG);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();

//		boolean fixed = true;
//		double[] sigmas = {0.7, 0.5, 0.3};
		boolean fixed = false;
		double[] sigmas = { 1d, 0.66, 0.33};
		
		ModAttenuationRelationship modAttenRel = new ModAttenuationRelationship(AttenRelRef.ASK_2014,
				fixed ? ModAttenRelRef.FIXED_STD_DEV : ModAttenRelRef.SIMPLE_SCALE);
		modAttenRel.setParamDefaults();
		Site site = new Site(new Location(34.0192, -118.286)); // USC
		double period = 3d;
		
		modAttenRel.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(modAttenRel.getIntensityMeasure(), period);
		
		OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
		SiteTranslator trans = new SiteTranslator();
		ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());
		
		for (Parameter<?> param : modAttenRel.getSiteParams()) {
			trans.setParameterValue(param, datas);
			site.addParameter(param);
		}
		
		Color[] colors = {Color.BLACK, Color.DARK_GRAY, Color.GRAY};
		
		double minX = 1e-3;
		double maxX = 1e1;
		
		DiscretizedFunc logXVals = new EvenlyDiscretizedFunc(Math.log(minX), Math.log(maxX), 100);
		DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logXVals)
			xVals.set(Math.exp(pt.getX()), 1d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<sigmas.length; i++) {
			if (fixed)
				((FixedStdDevMod)modAttenRel.getCurrentMod()).setStdDev(sigmas[i]);
			else
				((SimpleScaleMod)modAttenRel.getCurrentMod()).setStdDevScaleFactor(sigmas[i]);
			
			calc.getHazardCurve(logXVals, site, modAttenRel, erf);
			
			DiscretizedFunc curve = xVals.deepClone();
			for (int j=0; j<curve.size(); j++)
				curve.set(j, logXVals.getY(j));
			
			if (fixed)
				curve.setName("σ="+(float)sigmas[i]+" ");
			else {
				if (sigmas[i] == 1d)
					curve.setName("Total σ ");
				else
					curve.setName((float)sigmas[i]+"*σ ");
			}
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, colors[i]));
		}
		
//		String title = "Sigma Comparison";
		String title = " ";
		PlotSpec spec = new PlotSpec(funcs, chars, title, (float)period+"s SA (g)", "Annual Probability of Exceedance");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		Range xRange = new Range(minX, maxX);
		Range yRange = new Range(1e-8, 1e0);
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		funcs.get(0).setName("Hazard Curve");
		spec.setPlotElems(funcs.subList(0, 1));
		spec.setChars(chars.subList(0, 1));
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		file = new File(outputDir, prefix+"_main");
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}

}
