package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;
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
import org.opensha.commons.util.ReturnPeriodUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.imr.attenRelImpl.ngaw2.IMT;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;

public class GMM_SigmaDists {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "sigma_dists";
		
		double[] mags = { 6.5d, 7.5d, 8.5d };
		double[] dists =  { 25d, 50d, 100d };
		double[] probs = { 0.1, 0.01, 0.001 };
		double[] sigmaScales = { 0.666, 0.333 };
		
		Color[] colors = { Color.BLUE.darker(), Color.RED.darker(), Color.green.darker() };
		
		ASK_2014 gmpe = new ASK_2014();
		
		gmpe.set_dip(90d);
		gmpe.set_fault(FaultStyle.STRIKE_SLIP);
		gmpe.set_IMT(IMT.PGA);
		gmpe.set_vs30(760d);
		gmpe.set_vsInf(true);
		gmpe.set_width(12d);
		gmpe.set_z1p0(Double.NaN);
		gmpe.set_z2p5(Double.NaN);
		gmpe.set_zHyp(6d);
		gmpe.set_zTop(0d);
		
		Range xRange = new Range(1e-3, 1e1);
		double logLower = Math.log10(xRange.getLowerBound());
		double logUpper = Math.log10(xRange.getUpperBound());
		double logDelta = 0.01d;
		int numLog = (int)((logUpper - logLower)/logDelta)+1;
		
		EvenlyDiscretizedFunc hazardCurve = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
		EvenlyDiscretizedFunc[] magCurves = new EvenlyDiscretizedFunc[mags.length];
		for (int i=0; i<mags.length; i++)
			magCurves[i] = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
		EvenlyDiscretizedFunc[] distCurves = new EvenlyDiscretizedFunc[dists.length];
		for (int i=0; i<dists.length; i++)
			distCurves[i] = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
		EvenlyDiscretizedFunc[] sigmaCurves = new EvenlyDiscretizedFunc[sigmaScales.length];
		for (int i=0; i<sigmaScales.length; i++)
			sigmaCurves[i] = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
		// init to 1, non-exceedance
		for (int i=0; i<numLog; i++) {
			hazardCurve.set(i, 1d);
			for (EvenlyDiscretizedFunc curve : magCurves)
				curve.set(i, 1d);
			for (EvenlyDiscretizedFunc curve : distCurves)
				curve.set(i, 1d);
			for (EvenlyDiscretizedFunc curve : sigmaCurves)
				curve.set(i, 1d);
		}
		
		DecimalFormat sigDF = new DecimalFormat("0.00");
		
		for (int m=0; m<mags.length; m++) {
			double mag = mags[m];
			gmpe.set_Mw(mag);
			for (int d=0; d<dists.length; d++) {
				double dist = dists[d];
				gmpe.set_rJB(dist);
				gmpe.set_rRup(dist);
				gmpe.set_rX(dist);
				
				ScalarGroundMotion gm = gmpe.calc();
				
				NormalDistribution normDist = new NormalDistribution(gm.mean(), gm.stdDev());
				
				EvenlyDiscretizedFunc pdf = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
				EvenlyDiscretizedFunc exceed = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
				for (int i=0; i<numLog; i++) {
					double x = Math.log(Math.pow(10, pdf.getX(i)));
					pdf.set(i, normDist.density(x));
					exceed.set(i, 1d-normDist.cumulativeProbability(x));
					hazardCurve.set(i, hazardCurve.getY(i)*Math.pow(1-probs[m], exceed.getY(i)));
					magCurves[m].set(i, magCurves[m].getY(i)*Math.pow(1-probs[m], exceed.getY(i)));
					distCurves[d].set(i, distCurves[d].getY(i)*Math.pow(1-probs[m], exceed.getY(i)));
				}
				
				EvenlyDiscretizedFunc[] scaledPDFs = new EvenlyDiscretizedFunc[sigmaScales.length];
				EvenlyDiscretizedFunc[] scaledExceeds = new EvenlyDiscretizedFunc[sigmaScales.length];
				
				// now scaled sigma curves
				for (int s=0; s<sigmaScales.length; s++) {
					NormalDistribution scaleDist = new NormalDistribution(gm.mean(), sigmaScales[s]*gm.stdDev());
					
					scaledPDFs[s] = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
					scaledExceeds[s] = new EvenlyDiscretizedFunc(logLower, numLog, logDelta);
					
					for (int i=0; i<numLog; i++) {
						double x = Math.log(Math.pow(10, pdf.getX(i)));
						scaledPDFs[s].set(i, scaleDist.density(x));
						scaledExceeds[s].set(i, 1d-scaleDist.cumulativeProbability(x));
						sigmaCurves[s].set(i, sigmaCurves[s].getY(i)*Math.pow(1-probs[m], scaledExceeds[s].getY(i)));
					}
				}
				
//				System.out.println("M"+(float)mag+", R="+(int)dist);
//				System.out.println("PDF");
//				System.out.println(pdf);
//				System.out.println("\nExceed");
//				System.out.println(exceed);
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				funcs.add(toLinear(pdf, false));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				PlotSpec spec = new PlotSpec(funcs, chars, " ", "PGA (g)", "Probability Density");
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(24);
				gp.setLegendFontSize(28);
				gp.setBackgroundColor(Color.WHITE);
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				Range yRange = null;
				
				gp.drawGraphPanel(spec, true, false, xRange, yRange);
				gp.getYAxis().setTickLabelsVisible(false);
				gp.getYAxis().setTickMarksVisible(false);
				gp.getYAxis().setMinorTickMarksVisible(false);
				gp.getPlot().setRangeGridlinesVisible(false);
				
				String magDistPrefix = prefix+"_m"+(float)mag+"_d"+(int)dist;
				
				File file = new File(outputDir, magDistPrefix+"_pdf");
				gp.getChartPanel().setSize(800, 500);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				if (sigmaScales.length > 0 && m == mags.length-1 && d == dists.length-1) {
					// do sigma scaled curves as well
					funcs.get(0).setName("σ="+sigDF.format(gm.stdDev()));
					for (int s=0; s<sigmaScales.length; s++) {
						scaledPDFs[s].setName("σ="+sigDF.format(sigmaScales[s]*gm.stdDev()));
						funcs.add(toLinear(scaledPDFs[s], false));
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[s]));
					}
					spec.setLegendVisible(true);
					
					gp.drawGraphPanel(spec, true, false, xRange, yRange);
					gp.getYAxis().setTickLabelsVisible(false);
					gp.getYAxis().setTickMarksVisible(false);
					gp.getYAxis().setMinorTickMarksVisible(false);
					gp.getPlot().setRangeGridlinesVisible(false);
					
					file = new File(outputDir, magDistPrefix+"_pdf_sigma");
					gp.getChartPanel().setSize(800, 500);
					gp.saveAsPNG(file.getAbsolutePath()+".png");
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				}
				
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				funcs.add(toLinear(exceed, false));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				spec = new PlotSpec(funcs, chars, " ", "PGA (g)", "Exceedance Probability");
				
				gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(24);
				gp.setLegendFontSize(28);
				gp.setBackgroundColor(Color.WHITE);
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				yRange = new Range(0d, 1d);
				
				gp.drawGraphPanel(spec, true, false, xRange, yRange);
				
				file = new File(outputDir, magDistPrefix+"_exceed");
				gp.getChartPanel().setSize(800, 500);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				if (sigmaScales.length > 0 && m == mags.length-1 && d == dists.length-1) {
					// do sigma scaled curves as well
					funcs.get(0).setName("σ="+sigDF.format(gm.stdDev()));
					for (int s=0; s<sigmaScales.length; s++) {
						scaledExceeds[s].setName(scaledPDFs[s].getName());
						funcs.add(toLinear(scaledExceeds[s], false));
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[s]));
					}
					spec.setLegendVisible(true);
					
					gp.drawGraphPanel(spec, true, false, xRange, yRange);
					
					file = new File(outputDir, magDistPrefix+"_exceed_sigma");
					gp.getChartPanel().setSize(800, 500);
					gp.saveAsPNG(file.getAbsolutePath()+".png");
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				}
			}
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		hazardCurve.setName("Total Hazard");
		
		funcs.add(toLinear(hazardCurve, true));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "PGA (g)", "Annual Probability of Exceedance");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		Range yRange = new Range(1e-8, 1e0);
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		File file = new File(outputDir, prefix+"_curve");
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		// now mag curves
		for (int i=0; i<magCurves.length; i++) {
			magCurves[i].setName("M"+(float)mags[i]);
			funcs.add(toLinear(magCurves[i], true));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[i]));
		}
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		file = new File(outputDir, prefix+"_curve_mags");
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		while (funcs.size() > 1) {
			funcs.remove(funcs.size()-1);
			chars.remove(chars.size()-1);
		}
		
		// now dist curves
		for (int i=0; i<distCurves.length; i++) {
			distCurves[i].setName((int)dists[i]+" km");
			funcs.add(toLinear(distCurves[i], true));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[i]));
		}
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		file = new File(outputDir, prefix+"_curve_dists");
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		while (funcs.size() > 1) {
			funcs.remove(funcs.size()-1);
			chars.remove(chars.size()-1);
		}
		
		// now sigma curves
		funcs.get(0).setName("Full σ");
		for (int i=0; i<sigmaCurves.length; i++) {
			sigmaCurves[i].setName(sigDF.format(sigmaScales[i])+"*σ");
			funcs.add(toLinear(sigmaCurves[i], true));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[i]));
		}
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		file = new File(outputDir, prefix+"_curve_sigmas");
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private static DiscretizedFunc toLinear(DiscretizedFunc log10Func, boolean nonExceed) {
		DiscretizedFunc linear = new ArbitrarilyDiscretizedFunc(log10Func.getName());
		
		for (Point2D pt : log10Func)
			if (nonExceed)
				linear.set(Math.pow(10, pt.getX()), 1d-pt.getY());
			else
				linear.set(Math.pow(10, pt.getX()), pt.getY());
		
		return linear;
	}

}
