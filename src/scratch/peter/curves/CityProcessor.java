package scratch.peter.curves;

import static com.google.common.base.Preconditions.checkArgument;
import static org.opensha.nshmp2.util.Period.*;
import static scratch.peter.curves.ProbOfExceed.*;
import static org.opensha.sra.rtgm.RTGM.Frequency.*;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.util.NSHMP_Utils;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

class CityProcessor implements Runnable {

	private ScalarIMR imr;
	private ERF erf;
	private Iterable<NEHRP_TestCity> locs;
	private Period per;
	private boolean doRTGM = false;
	
	private String outDir;
	private static final String S = File.separator;
	private static final double TIME = 1;
	
	private List<List<String>> curveData;

	CityProcessor(ScalarIMR imr, ERF erf, Iterable<NEHRP_TestCity> locs,
		Period per, String outDir) {
		this.outDir = outDir;
		this.imr = imr;
		this.erf = erf;
		this.locs = locs;
		this.per = per;
		doRTGM = per == GM1P00 || per == GM0P20;
	}

	@Override
	public void run() {
		System.out.println("Starting: " + toString());
		init();
		HazardCurveCalculator calc = new HazardCurveCalculator();
		calc.setPtSrcDistCorrType(PtSrcDistCorr.Type.NSHMP08);
		for (NEHRP_TestCity loc : locs) {
//			DiscretizedFunc f = per.getLogFunction();
			DiscretizedFunc f = per.getFunction();
			Site site = loc.getSite();
			try {
				f = calc.getHazardCurve(f, site, imr, erf);
//				f = deLog(f);
				
				for (Point2D p : f) {
					f.set(p.getX(), NSHMP_Utils.probToRate(p.getY(), 1));
				}

//				f = calc.getAnnualizedRates(f, TIME);
//				System.out.println(f);
				addResults(loc, f);

			} catch (Exception e) {
				e.printStackTrace();
			}
			System.out.println(loc.toString().substring(0, 6) + " " + per + 
				" " + imr.getShortName());
		}
		writeFiles();
		System.out.println("Finished: " + toString());
	}
	
	private static DiscretizedFunc deLog(DiscretizedFunc f) {
		DiscretizedFunc fOut = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<f.size(); i++) {
			fOut.set(Math.exp(f.getX(i)), f.getY(i));
		}
		return fOut;
	}
	
	private void addResults(NEHRP_TestCity loc, DiscretizedFunc f) {
		double pe2in50 = ProbOfExceed.get(f, PE2IN50);
		double pe10in50 = ProbOfExceed.get(f, PE10IN50);
		double rtgm = (doRTGM) ? getRTGM(f) : 0;
		
		// curve data
		List<String> curveDat = Lists.newArrayList();
		curveDat.add(loc.name());
		curveDat.add(Double.toString(pe2in50));
		curveDat.add(Double.toString(pe10in50));
		curveDat.add(Double.toString(rtgm));
		for (Point2D p : f) {
			curveDat.add(Double.toString(p.getY()));
		}
		curveData.add(curveDat);
	}
	
	
	private double getRTGM(DiscretizedFunc f) {
		Frequency freq = per.equals(GM0P20) ? SA_0P20 : SA_1P00;
		RTGM rtgm = RTGM.create(f, freq, 0.8).call();
		return rtgm.get();
	}

	private void writeFiles() {
		String outDirName = outDir + S + per + S;
		File outDir = new File(outDirName);
		outDir.mkdirs();
		String curveFile = outDirName +  imr.getShortName() + "_curves.csv";
		toCSV(new File(curveFile), curveData);
	}
	
	private static void toCSV(File file, List<List<String>> content) {
		if (file.exists()) file.delete();
		Joiner joiner = Joiner.on(',').useForNull(" ");
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(file, true));
			for (List<String> lineDat : content) {
				String line = joiner.join(lineDat);
				pw.println(line);
			}
			pw.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	@Override
	public String toString() {
		return "  " + imr.getShortName() + " " + per;
	}
	
	private void init() {
		curveData = Lists.newArrayList();
		List<String> curveHeader = Lists.newArrayList();
		curveHeader.add("city");
		curveHeader.add("2in50");
		curveHeader.add("10in50");
		curveHeader.add("rtgm");
		
		for (Double d : per.getIMLs()) {
			curveHeader.add(d.toString());
		}
		curveData.add(curveHeader);
	}

}
