package scratch.kevin;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveWriter;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;
import org.opensha.sha.imr.param.IntensityMeasureParams.MMI_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class BinaryCurvesToMMI_Calc {
	
	public static void writeMMI(File pgaFile, File pgvFile, File mmiFile) throws Exception {
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(MMI_Param.NAME);
		
		writeMMI(pgaFile, pgvFile, mmiFile, xVals);
	}
	
	public static void writeMMI(File pgaFile, File pgvFile, File mmiFile, DiscretizedFunc xVals) throws Exception {
		BinaryHazardCurveReader pgaReader = new BinaryHazardCurveReader(pgaFile.getAbsolutePath());
		BinaryHazardCurveReader pgvReader = new BinaryHazardCurveReader(pgvFile.getAbsolutePath());
		
		List<Location> locs = Lists.newArrayList();
		List<DiscretizedFunc> mmiCurves = Lists.newArrayList();
		
		while (true) {
			int index = locs.size();
			boolean print = index < 10000 && index % 1000 == 0 || index % 10000 == 0;
			if (print) {
				String validateStr = "";
				if (validate_tol_percent > 0)
					validateStr = " (max PDiff: "+(float)maxPDiff+" %)";
				System.out.println("Processing site "+locs.size()+validateStr);
			}
			ArbitrarilyDiscretizedFunc pgaCurve = pgaReader.nextCurve();
			ArbitrarilyDiscretizedFunc pgvCurve = pgvReader.nextCurve();
			
			if (pgaCurve == null) {
				Preconditions.checkState(pgvCurve == null);
				break;
			}
			
			Location loc = pgaReader.currentLocation();
			Preconditions.checkState(loc.equals(pgvReader.currentLocation()));
			
			DiscretizedFunc mmiCurve = getMMI(pgaCurve, pgvCurve, xVals);
			locs.add(loc);
			mmiCurves.add(mmiCurve);
			
			if (debug_plot_interval > 0 && index % debug_plot_interval == 0)
				debugPlot(pgaCurve, pgvCurve, mmiCurve);
		}
		
		BinaryHazardCurveWriter writer = new BinaryHazardCurveWriter(mmiFile);
		writer.writeCurves(mmiCurves, locs);
	}
	
	private static final int num_interp = 1000;
	private static final boolean prob_interp_log = true;
	private static final boolean im_interp_log = false;
	private static final boolean interp_rate_space = false;
	
	private static DiscretizedFunc superSampleEnds(DiscretizedFunc func, int n, int factor, boolean first, boolean last) {
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		
		Preconditions.checkState(n >= 1);
		Preconditions.checkState(factor > 1);
		Preconditions.checkState(first || last);
		
		for (Point2D pt : func)
			ret.set(pt.getX(), 0d);
		
		EvenlyDiscretizedFunc superSample;
		if (first) {
			superSample = new EvenlyDiscretizedFunc(func.getX(0), func.getX(n), n*factor);
			for (Point2D pt : superSample)
				ret.set(pt.getX(), 0d);
		}
		if (last) {
			superSample = new EvenlyDiscretizedFunc(func.getX(func.size()-1-n), func.getX(func.size()-1), n*factor);
			for (Point2D pt : superSample)
				ret.set(pt.getX(), 0d);
		}
		
		return ret;
	}
	
	private static DiscretizedFunc getMMI(DiscretizedFunc pgaCurve, DiscretizedFunc pgvCurve, DiscretizedFunc xVals) {
		DiscretizedFunc calcPGA_curve, calcPGV_curve;
		if (im_interp_log) {
			calcPGA_curve = getLogFunc(pgaCurve);
			calcPGV_curve = getLogFunc(pgvCurve);
		} else {
			calcPGA_curve = pgaCurve;
			calcPGV_curve = pgvCurve;
		}
		
		if (interp_rate_space) {
			calcPGA_curve = probToRate(calcPGA_curve);
			calcPGV_curve = probToRate(calcPGV_curve);
		}
		
		// initialize to zero
		DiscretizedFunc mmiCurve = xVals.deepClone();
		for (int i=0; i<mmiCurve.size(); i++)
			mmiCurve.set(i, 0d);
		
		// first compute really high resolution MMI curve by fixing Y values and grabbing X from each curve
		double upperProb = Math.min(calcPGA_curve.getMaxY(), calcPGV_curve.getMaxY());
		if (upperProb == 0d)
			// zero prob curve
			return mmiCurve;
		double lowerProb = Math.max(minAboveZero(calcPGA_curve), minAboveZero(calcPGV_curve));
		Preconditions.checkState(!Double.isInfinite(lowerProb));
		// ln(prob) here is X, MMI is Y
		DiscretizedFunc highResFunc;
		if (prob_interp_log)
			highResFunc = new EvenlyDiscretizedFunc(Math.log(lowerProb), Math.log(upperProb), num_interp);
		else
			highResFunc = new EvenlyDiscretizedFunc(lowerProb, upperProb, num_interp);
		// now super sample the really high probability reagion
		highResFunc = superSampleEnds(highResFunc, num_interp/100, 200, true, true);
//		System.out.println("Now interpolating on "+highResFunc.size()+" points (orig "+num_interp+")");
		for (int i=0; i<highResFunc.size(); i++) {
			double prob = highResFunc.getX(i);
			if (prob_interp_log)
				prob = Math.exp(prob);
			// now get IMs at this POE from each curve
			double pga = interpX(calcPGA_curve, prob, false);
			double pgv = interpX(calcPGV_curve, prob, false);
			if (im_interp_log) {
				pga = Math.exp(pga);
				pgv = Math.exp(pgv);
			}
			double mmi = Wald_MMI_Calc.getMMI(pga, pgv);
			highResFunc.set(i, mmi);
		}
//		System.out.println("Upper Prob: "+upperProb);
//		System.out.println("Lower Prob: "+lowerProb);
//		System.out.println(highResFunc);
//		new GraphWindow(highResFunc, "High Res Func");
		
		double minMMI = highResFunc.getMinY();
		double maxMMI = highResFunc.getMaxY();
		// now interpolate onto expected X values
		// go backwards so that you can fill in with previous vals
		for (int i=mmiCurve.size(); --i>=0;) {
			double mmi = mmiCurve.getX(i);
//			System.out.println("");
			if (mmi > maxMMI)
				// to the right of the actual values, leave at zero
				continue;
			if (mmi < minMMI) {
				// preserve value to the right
				mmiCurve.set(i, mmiCurve.getY(i+1));
				continue;
			}
			// inerpolate
			double prob = interpX(highResFunc, mmi, false);
			if (prob_interp_log)
				prob = Math.exp(prob);
			if (interp_rate_space)
				prob = rateToProb(prob);
//			System.out.println("MMI: "+mmi+", ln(P)="+Math.log(prob)+", P="+prob);
			mmiCurve.set(i, prob);
		}
		
		if (validate_tol_percent > 0)
			validateMMI(pgaCurve, pgvCurve, mmiCurve, validate_tol_percent);
		
		return mmiCurve;
	}
	
	private static DiscretizedFunc getLogFunc(DiscretizedFunc input) {
		ArbitrarilyDiscretizedFunc lnFunc = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : input)
			lnFunc.set(Math.log(pt.getX()), pt.getY());
		return lnFunc;
	}
	
	private static DiscretizedFunc probToRate(DiscretizedFunc curve) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : curve)
			ret.set(pt.getX(), probToRate(pt.getY()));
		return ret;
	}
	
	private static double probToRate(double prob) {
		return -Math.log(1d - prob);
	}
	
	private static double rateToProb(double rate) {
		return 1 - Math.exp(-rate);
	}
	
	private static double maxPDiff = 0d;
	
	private static void validateMMI(DiscretizedFunc pgaCurve, DiscretizedFunc pgvCurve, DiscretizedFunc mmiCurve, double tolPercent) {
		for (int i=0; i<mmiCurve.size(); i++) {
			double mmi = mmiCurve.getX(i);
			double prob = mmiCurve.getY(i);
//			System.out.println(prob);
			double pga = interpX(pgaCurve, prob, false);
			double pgv = interpX(pgvCurve, prob, false);
			Preconditions.checkState(Doubles.isFinite(pga));
			Preconditions.checkState(Doubles.isFinite(pgv));
			double calcMMI = Wald_MMI_Calc.getMMI(pga, pgv);
			double pDiff = DataUtils.getPercentDiff(mmi, calcMMI);
			maxPDiff = Math.max(maxPDiff, pDiff);
			Preconditions.checkState(pDiff <= tolPercent, "Bad MMI at i=%s. P=%s, PGA=%s, PGV=%s, MMI=%s, calcMMI=%s, pDiff=%s %",
					i, prob, pga, pgv, mmi, calcMMI, pDiff);
		}
	}
	
	private static double interpX(DiscretizedFunc curve, double y, boolean logDomain) {
		double firstY = curve.getY(0);
		if ((float)y == (float)firstY)
			return curve.getX(0);
		double lastY = curve.getY(curve.size()-1);
		if ((float)y == (float)lastY)
			return curve.getX(curve.size()-1);
		if (logDomain)
			return curve.getFirstInterpolatedX_inLogXLogYDomain(y);
		return curve.getFirstInterpolatedX(y);
	}
	
	private static double minAboveZero(DiscretizedFunc curve) {
		double min = Double.POSITIVE_INFINITY;
		for (Point2D pt : curve)
			if (pt.getY() > 0)
				min = Math.min(min, pt.getY());
		return min;
	}
	
	private static int debug_plot_interval = 0;
	private static boolean debug_wait = false;
	private static double validate_tol_percent = 0;
	
	private static void debugPlot(DiscretizedFunc pgaCurve, DiscretizedFunc pgvCurve, DiscretizedFunc mmiCurve) {
		List<PlotSpec> specs = Lists.newArrayList();
		specs.add(buildSpec(pgaCurve, "PGA", Color.BLACK));
		specs.add(buildSpec(pgvCurve, "PGV", Color.BLUE));
		specs.add(buildSpec(mmiCurve, "MMI", Color.GREEN.darker()));
		
		GraphWindow gw = null;
		for (PlotSpec spec : specs) {
			if (gw == null)
				gw = new GraphWindow(spec);
			else
				gw.addTab(spec);
			gw.setXLog(true);
			gw.setYLog(true);
		}
		while (debug_wait && gw.isVisible()) {
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}
	
	private static PlotSpec buildSpec(DiscretizedFunc curve, String imt, Color color) {
		List<DiscretizedFunc> curves = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		curves.add(curve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
		
		return new PlotSpec(curves, chars, imt+" Curve", imt, "Probability of Exceedance");
	}

	public static void main(String[] args) throws Exception {
		debug_plot_interval = 50000;
		debug_wait = false;
//		validate_tol_percent = 0.1;
		validate_tol_percent = 1;
		String[] suffixes = { "_combined", "_gridded", "fault" };
		File baseDir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard/"
				+ "2017_03_01-mojave_m7_fulltd_descendents-NGA2-0.01-site-effects-with-basin");
		
		for (String suffix : suffixes) {
			File pgaFile = new File(baseDir, "results_pga"+suffix+".bin");
			File pgvFile = new File(baseDir, "results_pgv"+suffix+".bin");
			File mmiFile = new File(baseDir, "results_mmi"+suffix+".bin");
			
			System.out.println("Processing "+mmiFile.getName());
			
			writeMMI(pgaFile, pgvFile, mmiFile);
		}
	}

}
