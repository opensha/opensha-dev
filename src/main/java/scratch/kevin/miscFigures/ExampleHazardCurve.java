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
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.kevin.util.ReturnPeriodUtils;

public class ExampleHazardCurve {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "hazard_curve";
		
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(Presets.BOTH_FM_BRANCH_AVG);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		Site site = new Site(new Location(34.0192, -118.286)); // USC
		double period = 5d;
		
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
		
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
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		calc.getHazardCurve(logXVals, site, gmpe, erf);
		
		DiscretizedFunc curve = xVals.deepClone();
		for (int j=0; j<curve.size(); j++)
			curve.set(j, logXVals.getY(j));
		
		double twoIn50_poe = ReturnPeriodUtils.calcExceedanceProb(0.02, 50d, 1d);
		System.out.println("2% in 50 means "+twoIn50_poe+" annual");
		System.out.println("2% in 50: "+curve.getFirstInterpolatedX_inLogXLogYDomain(twoIn50_poe));
		double annualProb0p1 = curve.getInterpolatedY_inLogXLogYDomain(0.1);
		System.out.println("0.1g: "+annualProb0p1+" annual prob");
		double prob0p1_50yr = ReturnPeriodUtils.calcExceedanceProb(annualProb0p1, 1d, 50d);
		System.out.println("0.1g: "+prob0p1_50yr+" 50yr prob");
		double annualProb0p2 = curve.getInterpolatedY_inLogXLogYDomain(0.2);
		System.out.println("0.2g: "+annualProb0p2+" annual prob");
		double prob0p2_50yr = ReturnPeriodUtils.calcExceedanceProb(annualProb0p2, 1d, 50d);
		System.out.println("0.2g: "+prob0p2_50yr+" 50yr prob");
		
		funcs.add(curve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Hazard Curve", (float)period+"s SA (g)", "Annual Probability of Exceedance");
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
	}

}
