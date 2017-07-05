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
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp2.calc.HazardCalc;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.erf.NSHMP2008;
import org.opensha.nshmp2.imr.NSHMP08_WUS;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

class CityProcessorNSHMP implements Runnable {

	private NSHMP2008 erf;
	private Map<String, Location> locMap;
	private Period per;
	private boolean epi;
	private boolean doRTGM = false;
	
	private String outDir;
	private static final String S = File.separator;
	private static final double TIME = 1;
	
	private List<List<String>> curveData;

	CityProcessorNSHMP(NSHMP2008 erf, Map<String, Location> locMap,
		Period per, boolean epi,  String outDir) {
		this.outDir = outDir;
		this.erf = erf;
		this.locMap = locMap;
		this.per = per;
		this.epi = epi;
		doRTGM = per == GM1P00 || per == GM0P20;
	}

	@Override
	public void run() {
		System.out.println("Starting: " + toString());
		init();
		for (String locName : locMap.keySet()) {
			Site site = new Site(locMap.get(locName));
			HazardCalc calc = HazardCalc.create(erf, site, per, epi);
			HazardResult result = calc.call();
			addResults(locName, result.curve());
			System.out.println(locName + " " + per);
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
	
	private void addResults(String name, DiscretizedFunc f) {
		double pe2in50 = getPE(f, PE2IN50);
		double pe10in50 = getPE(f, PE10IN50);
		double rtgm = getRTGM(f);
		
		// curve data
		List<String> curveDat = Lists.newArrayList();
		curveDat.add(name);
		curveDat.add(Double.toString(pe2in50));
		curveDat.add(Double.toString(pe10in50));
		if (doRTGM) curveDat.add(Double.toString(rtgm));
		for (Point2D p : f) {
			curveDat.add(Double.toString(p.getY()));
		}
		curveData.add(curveDat);
	}
	
	private double getPE(DiscretizedFunc f, ProbOfExceed pe) {
		return f.getFirstInterpolatedX_inLogXLogYDomain(pe.annualRate());
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
		String curveFile = outDirName +  "NSHMP08" + "_curves.csv";
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
		return "  " + "NSHMP08" + " " + per;
	}
	
	private void init() {
		curveData = Lists.newArrayList();
		List<String> curveHeader = Lists.newArrayList();
		curveHeader.add("city");
		curveHeader.add("2in50");
		curveHeader.add("10in50");
		if (doRTGM) curveHeader.add("rtgm");
		
		for (Double d : per.getIMLs()) {
			curveHeader.add(d.toString());
		}
		curveData.add(curveHeader);
	}

}
