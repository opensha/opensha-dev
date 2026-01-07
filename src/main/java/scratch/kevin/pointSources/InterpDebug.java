package scratch.kevin.pointSources;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.interp.DistanceInterpolator;
import org.opensha.sha.calc.PointSourceOptimizedExceedProbCalc;
import org.opensha.sha.calc.RuptureExceedProbCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import net.mahdilamb.colormap.Colors;

public class InterpDebug {

	public static void main(String[] args) throws IOException {
		RuptureExceedProbCalculator basic = RuptureExceedProbCalculator.BASIC_IMPLEMENTATION;
		RuptureExceedProbCalculator optimized = new PointSourceOptimizedExceedProbCalc();

		ScalarIMR imr = AttenRelRef.WRAPPED_ASK_2014.get();

		// Simple point rupture
		Location sourceLoc = new Location(0d, 0d, 5d);
		PointSurface surf = new PointSurface(sourceLoc);
		EqkRupture rup = new EqkRupture(6.5, 0d, surf, null);

		DistanceInterpolator distInterp = DistanceInterpolator.get();
		double distance = 12.34;
		int binBefore = distInterp.getIndexAtOrBefore(distance);
		Site site = new Site(LocationUtils.location(sourceLoc, 0d, distance));
		site.addParameterList(imr.getSiteParams());
		imr.setSite(site);

		// Default IMT for these tests
		imr.setIntensityMeasure(PGA_Param.NAME);

		// ln(IML) X values
		EvenlyDiscretizedFunc xVals50 = new EvenlyDiscretizedFunc(-3d, 0d, 50);
		
		EvenlyDiscretizedFunc basicCurve = xVals50.deepClone();
		EvenlyDiscretizedFunc optCurve = xVals50.deepClone();

		basic.getExceedProbabilities(imr, rup, basicCurve);
		optimized.getExceedProbabilities(imr, rup, optCurve);
		
		// now get values before and after
		site.setLocation(LocationUtils.location(sourceLoc, 0d, distInterp.getDistance(binBefore)));
		imr.setSite(site);
		EvenlyDiscretizedFunc curveBefore = xVals50.deepClone();
		basic.getExceedProbabilities(imr, rup, curveBefore);
		
		site.setLocation(LocationUtils.location(sourceLoc, 0d, distInterp.getDistance(binBefore+1)));
		imr.setSite(site);
		EvenlyDiscretizedFunc curveAfter = xVals50.deepClone();
		basic.getExceedProbabilities(imr, rup, curveAfter);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		curveBefore.setName("Before ("+(float)distInterp.getDistance(binBefore)+" km)");
		funcs.add(curveBefore);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_blue));
		
		curveAfter.setName("After ("+(float)distInterp.getDistance(binBefore+1)+" km)");
		funcs.add(curveAfter);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_green));
		
		basicCurve.setName("Basic Calc ("+(float)distance+" km)");
		funcs.add(basicCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		optCurve.setName("Interpolated Calc ("+(float)distance+" km)");
		funcs.add(optCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, Colors.tab_orange));
		
		for (int i=0; i<xVals50.size(); i++) {
			double basicY = basicCurve.getY(i);
			double optY = optCurve.getY(i);
			double diff = optY - basicY;
			double pDiff = diff*100d/basicY;
			System.out.println(i+". "+(float)xVals50.getX(i)+"\tbasic="+(float)basicY+"\topt="+(float)optY
					+"\tdiff="+(float)diff+"\t("+pDiff+"%);\t\t["+(float)curveBefore.getY(i)+", "+(float)curveAfter.getY(i)+"]");
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Ln(IML)", "POE");
		plot.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot);
		
		PlotUtils.writePlots(new File("/tmp"), "interp_debug", gp, 5000, 800, true, false, false);
	}

}
