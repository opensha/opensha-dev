package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

class SupraBValFigure {

	public static void main(String[] args) throws IOException {
		EvenlyDiscretizedFunc refMFD = new EvenlyDiscretizedFunc(6.55, 10, 0.1);
		RupSetScalingRelationship scale = NSHM23_ScalingRelationships.LOGA_C4p2;
		boolean cmlOffset = false;
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		String prefix = "supra_mfd_sweep";
		
		double areaMinRup = Math.pow(10, refMFD.getMinX() - 4.2);
		System.out.println("Area for M"+(float)refMFD.getMinX()+": "+areaMinRup);
		double areaMaxRup = Math.pow(10, refMFD.getMaxX() - 4.2);
		System.out.println("Area for M"+(float)refMFD.getMaxX()+": "+areaMaxRup);
		
		double ddw = Math.sqrt(areaMinRup);
		double len = areaMaxRup/ddw;
		
//		double ddw = 15d;
//		double len = 120d;
//		double areaMinRup = ddw*ddw;
//		double areaMaxRup = len*ddw;
//		System.out.println("Actual mag for smallest rup: "+scale.getMag(areaMinRup*1e6, (len*areaMinRup/areaMaxRup)*1e3, ddw*1e3, ddw*1e3, Double.NaN));
//		System.out.println("Actual mag for largest rup: "+scale.getMag(areaMaxRup*1e6, len*1e3, ddw*1e3, ddw*1e3, Double.NaN));
		
		double area = len*ddw;
		System.out.println("Fault dimensions: "+(float)len+" x "+(float)ddw+" = "+(float)area);
		
		double slipRate = 10 * 1e-3; // mm/yr -> m/yr
		
		double lenMinRup = len*areaMinRup/area;
		double lenMaxRup = len;
		
		double slipInMinRup = scale.getAveSlip(areaMinRup*1e6, lenMinRup*1e3, ddw*1e3, ddw*1e3, Double.NaN);
		double slipInMaxRup = scale.getAveSlip(areaMaxRup*1e6, lenMaxRup*1e3, ddw*1e3, ddw*1e3, Double.NaN);
		
		System.out.println("Slip for M"+(float)refMFD.getMinX()+": "+slipInMinRup);
		System.out.println("Slip for M"+(float)refMFD.getMaxX()+": "+slipInMaxRup);
		
		double moRate = FaultMomentCalc.getMoment(area*1e6, slipRate);
		
		double moMinRup = FaultMomentCalc.getMoment(areaMinRup*1e6, slipInMinRup);
		double moMaxRup = FaultMomentCalc.getMoment(areaMaxRup*1e6, slipInMaxRup);
		
		double rateMinRup = moRate / moMinRup;
		System.out.println("Maximum rate model: "+(float)rateMinRup+" ("+(float)(1d/rateMinRup)+" yrs)");
		double rateMaxRup = moRate / moMaxRup;
		System.out.println("Minimum rate model: "+(float)rateMaxRup+" ("+(float)(1d/rateMaxRup)+" yrs)");
		
		double[] bVals = { 1d, 0.75d, 0.5d, 0.25d, 0d };
		
		List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
		List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> cmlChars = new ArrayList<>();
		
		// max rate model
		IncrementalMagFreqDist maxRateMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		maxRateMFD.set(0, rateMinRup);
		maxRateMFD.setName("Maximum Rate Model");
		incrFuncs.add(maxRateMFD);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		cmlFuncs.add(maxRateMFD);
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		
		// max rate model
		IncrementalMagFreqDist minRateMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		minRateMFD.setName("Minimum Rate Model");
		minRateMFD.set(minRateMFD.size()-1, rateMaxRup);
		incrFuncs.add(minRateMFD);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		cmlFuncs.add(minRateMFD);
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		ArbitrarilyDiscretizedFunc minRateDashed = new ArbitrarilyDiscretizedFunc();
		minRateDashed.set(refMFD.getMinX()-0.45*refMFD.getDelta(), rateMaxRup);
		minRateDashed.set(refMFD.getMaxX()+0.45*refMFD.getDelta(), rateMaxRup);
		cmlFuncs.add(minRateDashed);
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, new Color(0, 0, 255, 127)));
		
		// GR+char
		double fractGR1 = 0.5;
		IncrementalMagFreqDist grPlusChar1 = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		grPlusChar1.setName("1/2 G-R, 1/2 Char");
		double fractGR2 = 2d/3d;
		IncrementalMagFreqDist grPlusChar2 = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		grPlusChar2.setName("2/3 G-R, 1/3 Char");
		GutenbergRichterMagFreqDist grb1 = new GutenbergRichterMagFreqDist(
				refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), moRate, 1d);
		for (int i=0; i<grPlusChar1.size(); i++) {
			grPlusChar1.set(i, fractGR1*grb1.getY(i) + (1d-fractGR1)*minRateMFD.getY(i));
			grPlusChar2.set(i, fractGR2*grb1.getY(i) + (1d-fractGR2)*minRateMFD.getY(i));
		}
		incrFuncs.add(grPlusChar1);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.GREEN.darker()));
		cmlFuncs.add(cml(grPlusChar1, cmlOffset));
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.GREEN.darker()));
		incrFuncs.add(grPlusChar2);
		incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.GREEN.brighter()));
		cmlFuncs.add(cml(grPlusChar2, cmlOffset));
		cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.GREEN.brighter()));
		
		for (double bVal : bVals) {
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
					refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), moRate, bVal);
			double totRate = gr.calcSumOfY_Vals();
			System.out.println("G-R b="+(float)bVal+" model: "+(float)totRate+" ("+(float)(1d/totRate)+" yrs)");
			
			gr.setName("G-R b="+(float)bVal);
			
			incrFuncs.add(gr);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
			cmlFuncs.add(cml(gr, cmlOffset));
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		}
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setTickLabelFontSize(22);
		gp.setAxisLabelFontSize(26);
		
		Range xRange = new Range(refMFD.getMinX() - 0.5*refMFD.getDelta(), refMFD.getMaxX() + 0.5*refMFD.getDelta());
		double maxY = rateMinRup;
		double minY = rateMaxRup;
		for (IncrementalMagFreqDist func : incrFuncs)
			for (Point2D pt : func)
				if (pt.getY() > 0)
					minY = Math.min(minY, pt.getY());
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minY))), Math.pow(10, Math.ceil(Math.log10(maxY))));
		
		PlotSpec spec = new PlotSpec(incrFuncs, incrChars, " ", "Magnitude", "Incremental Rate (1/yr)");
		spec.setLegendInset(true);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_incr", gp, 800, 750, true, true, false);
		
		spec = new PlotSpec(cmlFuncs, cmlChars, " ", "Magnitude", "Cumulative Rate (1/yr)");
		spec.setLegendInset(true);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, true, false);
		
		// now plot tot rate as a function of bval
		xRange = new Range(-1d, 2d);
		yRange = new Range(1e-3, 1e-1);
		EvenlyDiscretizedFunc bRateFunc = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), 1000);
		for (int i=0; i<bRateFunc.size(); i++) {
			double bVal = bRateFunc.getX(i);
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
					refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), moRate, bVal);
			double totRate = gr.calcSumOfY_Vals();
			bRateFunc.set(i, totRate);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet maxFunc = new DefaultXY_DataSet();
		maxFunc.set(xRange.getLowerBound(), rateMaxRup);
		maxFunc.set(xRange.getUpperBound(), rateMaxRup);
		funcs.add(maxFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		DefaultXY_DataSet minFunc = new DefaultXY_DataSet();
		minFunc.set(xRange.getLowerBound(), rateMinRup);
		minFunc.set(xRange.getUpperBound(), rateMinRup);
		funcs.add(minFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		
		funcs.add(bRateFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		List<XYTextAnnotation> anns = new ArrayList<>();
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 20);
		
		for (double bVal : bVals) {
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
					refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), moRate, bVal);
			double totRate = gr.calcSumOfY_Vals();
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			xy.set(xRange.getLowerBound(), totRate);
			xy.set(bVal, totRate);
			xy.set(bVal, yRange.getLowerBound());
			
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			xy = new DefaultXY_DataSet();
			xy.set(bVal, totRate);
			
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.GRAY));
			
			String text = (int)(1d/totRate+0.5)+" yrs  ";
			XYTextAnnotation ann = new XYTextAnnotation(text, bVal, logShift(totRate, 0.01));
			ann.setTextAnchor(TextAnchor.BASELINE_RIGHT);
			ann.setFont(font);
			anns.add(ann);
			
			text = "  b="+(float)bVal;
			ann = new XYTextAnnotation(text, bVal+0.015, totRate);
			ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
			ann.setRotationAnchor(TextAnchor.BASELINE_LEFT);
			ann.setRotationAngle(Math.PI/2d);
			ann.setFont(font);
			anns.add(ann);
		}
		
		XYTextAnnotation ann = new XYTextAnnotation("  Max Rate: "+(int)(1d/rateMinRup+0.5)+" yrs",
				xRange.getLowerBound(), logShift(rateMinRup, 0.01));
		ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
		ann.setFont(font);
		anns.add(ann);
		ann = new XYTextAnnotation("  Min Rate: "+(int)(1d/rateMaxRup+0.5)+" yrs",
				xRange.getLowerBound(), rateMaxRup);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(font);
		anns.add(ann);
		
		spec = new PlotSpec(funcs, chars, " ", "G-R b-value", "Cumulative Rate (1/yr)");
		spec.setPlotAnnotations(anns);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		PlotUtils.setXTick(gp, 0.5);
		
		PlotUtils.writePlots(outputDir, prefix+"_cml_bVal", gp, 800, 750, true, true, false);
		
//		scale.getma
		
//		scale.get
//		double minRateSlip = scale.getAveSlip(slipRate, slipRate, slipRate, slipRate, slipRate)
	}
	
	private static double logShift(double x, double delta) {
		double logX = Math.log10(x);
		logX += delta;
		return Math.pow(10, logX);
	}
	
	private static EvenlyDiscretizedFunc cml(IncrementalMagFreqDist mfd, boolean offset) {
		if (offset)
			return mfd.getCumRateDistWithOffset();
		else
			return mfd.getCumRateDist();
	}

}
