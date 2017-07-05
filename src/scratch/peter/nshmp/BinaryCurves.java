package scratch.peter.nshmp;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.Interpolate;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.io.LittleEndianDataInputStream;
import com.google.common.io.LittleEndianDataOutputStream;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Doubles;

/**
 * Add comments here
 *
 * @author Peter Powers
 */
public class BinaryCurves {
	
	// needed becasue SH changed the NSHMP discretization
	static double[] xPGA = new double[] {0.005, 0.007, 0.0098, 0.0137, 0.0192, 0.0269,
	                              0.0376, 0.0527, 0.0738, 0.103, 0.145, 0.203,
	                              0.284, 0.397, 0.556, 0.778, 1.09, 1.52, 2.2,
	                              3.3};
	
	// our 0P30 calcs include an extra value (0.0025), should be identical to 0P20 (5Hz)
	static double[] x0P30 = Doubles.toArray(Period.GM0P20.getIMLs());

	public static void main(String[] args) throws IOException {
//		String file = "/Users/pmpowers/projects/svn/NSHMP/tmp/out/CAdeep_2014.5hz";
//		String file = "tmp/NSHMP-CA-binaries/CA-5Hz.curves";
//		read(new File(file).toURI().toURL());
		
		double[] xx = {1.52, 2.13};
		double[] yy = {1e-13, 1e-16};
		System.out.println(Interpolate.findLogLogY(xx, yy, 2.2));
		
//		Metadata meta = new Metadata();
//		meta.description = "UCERF3.3 Fault Sources";
//		meta.timestamp = (new Timestamp(System.currentTimeMillis())).toString();
//		meta.period = Period.GM0P00.getValue();
//		meta.Xs = xPGA; //Doubles.toArray(Period.GM0P00.getIMLs());
//		
//		String path = "tmp/binaryTest/test.curves";
//		File testFile = new File(path);
//		Files.createParentDirs(testFile);
//		
//		write(null, testFile, meta);
//		
//		read(testFile.toURI().toURL());
//		
	}
	

	public static void read(URL url) throws IOException {
		LittleEndianDataInputStream in = 
				new LittleEndianDataInputStream(url.openStream());

		CurveContainer cc = new CurveContainer();

		// read names 6 * char(128)
		int n = 128;
		for (int i=0; i<6; i++) {
			byte[] nameDat = new byte[n];
			in.read(nameDat, 0, n);
			System.out.println(new String(nameDat));
		}
		float period = in.readFloat();
		int nX = in.readInt();
		Period p = Period.valueForPeriod(period);
		System.out.println("period: " + period);
		System.out.println("nX: " + nX);
		
		// read x-vals real*4 * 20
		for (int i=0; i<20; i++) {
			double val = Precision.round((double) in.readFloat(), 3);
			System.out.println(val);
			// need to read 20 values to advance caret, but only save ones used
//			if (i<nX) cc.xs.add(val);
		}
//		System.out.println("xVals: " + cc.xs);

		// read extras real*4 * 10
		List<Double> extras = Lists.newArrayList();
		for (int i=0; i<10; i++) {
			double val = Precision.round((double) in.readFloat(), 2);
			extras.add(val);
 		}
		System.out.println("extras: " + extras);
		
		for (int i=0; i<1000; i++) {
			System.out.println(in.readFloat());
		}
		
//		System.out.println("C1:");
//		for (int i=0; i<19; i++) {
//			System.out.println(in.readFloat());
//		}
//		System.out.println("C2:");
//		for (int i=0; i<19; i++) {
//			System.out.println(in.readFloat());
//		}
		in.close();
	}
	
	public static void writeUC3(CurveContainer uc3, Period p, double spacing, String desc, File out) throws IOException {
		
		// create WUS cc for supplied curves
		CurveContainer cc = createNSHMP(uc3, spacing, p);
		
		Metadata meta = new Metadata();
		meta.description = desc;
		meta.timestamp = (new Timestamp(System.currentTimeMillis())).toString();
		meta.period = p;
		
		meta.Xs = (p == Period.GM0P00) ? xPGA : (p == Period.GM0P30) ? x0P30 : Doubles.toArray(p.getIMLs());
		
		Files.createParentDirs(out);
		write(cc, spacing, out, meta);
	}

	
	private static void write(CurveContainer cc, double outSpacing, File file, Metadata meta) throws IOException {
		LittleEndianDataOutputStream out = 
				new LittleEndianDataOutputStream(new FileOutputStream(file));

		// write info lines 6 * char(128)
		int n = 128;
		Charset charset = Charsets.US_ASCII;
		byte[] desc = Strings.padEnd(meta.description, n, ' ').getBytes(charset);
		byte[] ts = Strings.padEnd(meta.timestamp, n, ' ').getBytes(charset);
		byte[] dummy = Strings.padEnd("", n, ' ').getBytes(charset);
		out.write(desc);
		out.write(ts);
		for (int i=0; i<4; i++) {
			out.write(dummy);
		}
		
		out.writeFloat((float) meta.period.getValue());
		out.writeInt(meta.Xs.length);
		
		for (int i=0; i<meta.Xs.length; i++) {
			out.writeFloat((float) meta.Xs[i]);
		}
		// pad end of curve with 0s so that 20 values are present, even though
		// fewer may be used
		int pad = 20 - meta.Xs.length;
		for (int i=0; i<pad; i++) {
			out.writeFloat(0.0f);
		}
		
		// grid info
		float empty = -1.0f;
		float lonMin = (float) LON_MIN; // 2
		float lonMax = (float) LON_MAX; // 3
		float spacing = (float) outSpacing; // 4,7
		float latMin = (float) LAT_MIN; // 5
		float latMax = (float) LAT_MAX; // 6
		float nodeCt = 64005f; // 3
		float vs = 760.0f; // 9
		float sedDepth = 2.0f; // 10
		
		out.writeFloat(empty);
		out.writeFloat(lonMin);
		out.writeFloat(lonMax);
		out.writeFloat(spacing);
		out.writeFloat(latMin);
		out.writeFloat(latMax);
		out.writeFloat(spacing);
		out.writeFloat(nodeCt);
		out.writeFloat(vs);
		out.writeFloat(sedDepth);
		
		// may just be testing header
		if (cc == null) {
			out.close();
			return;
		}
		
		// write curves
		int nRows = (int) Math.rint((LAT_MAX - LAT_MIN) / outSpacing) + 1;
		int nCols = (int) Math.rint((LON_MAX - LON_MIN) / outSpacing) + 1;
		for (int i=0; i<nRows; i++) {
			double lat = LAT_MAX - outSpacing * i;
			for (int j=0; j<nCols; j++) {
				double lon = LON_MIN + outSpacing * j;
				Location loc = new Location(lat, lon);
				List<Double> vals = cc.getValues(loc);
				// we compute one too many values for 0.3s; strip first value to
				// bring array in line with 5Hz
				if (meta.period == Period.GM0P30) {
					vals = vals.subList(1, vals.size());
				}
				for (double val : vals) {
					out.writeFloat((float) val);
				}
			}
		}
		out.close();		
	}
	
	public static class Metadata {
		String description;
		String timestamp;
		Period period;
		double[] Xs;
		
		
	}
	
	
	private static final double LON_MIN = -125.0; // 2
	private static final double LON_MAX = -100.0; // 3
//	private static final double SPACING = 0.05; // 4,7
	private static final double LAT_MIN = 24.6; // 5
	private static final double LAT_MAX = 50.0; // 6

	// creates a WUS NSHMP curve container and populates it with the supplied
	// curves; really just handles a fix for large PGA values being slightly
	// different; SH extended the PGA range from
	// ... 1.09 1.52 2.13 to
	// ... 1.09 1.52 2.2 3.3
	public static CurveContainer createNSHMP(CurveContainer shaCC, double spacing, Period p) {
		GriddedRegion nshmpRegion = new GriddedRegion(
			new Location(LAT_MIN, LON_MIN),
			new Location(LAT_MAX, LON_MAX),
			spacing,
			GriddedRegion.ANCHOR_0_0);
		
		if (p != Period.GM0P00) {
			CurveContainer nshmpCC = CurveContainer.create(nshmpRegion, p.getIMLs());
			nshmpCC.union(shaCC);
			return nshmpCC;
		}
		
		// for PGA
		//		- create receiver container with 1 more x-value
		//		- replace 2.13 with extrapolated 2.2 (idx=18)
		// 		- set extrapolated 3.3 (idx=19)
		
		CurveContainer pgaCC = CurveContainer.create(
			TestGrid.CA_NSHMP.grid(spacing), Doubles.asList(xPGA));
		
		List<Double> xSrc = Period.GM0P00.getIMLs();
		double[] xs = {xSrc.get(17), xSrc.get(18)};
		
		for (Location loc : shaCC) {
			
			List<Double> ySrc = shaCC.getValues(loc);
			List<Double> yDest = pgaCC.getValues(loc);
			// copy values at indices 0-17
			for (int i=0; i<18; i++) {
				yDest.set(i, ySrc.get(i));
			}
			if (ySrc.get(17) <= 0.0) continue;
			double[] ys = {ySrc.get(17), ySrc.get(18)};
			double interp2p2 = Interpolate.findLogLogY(xs, ys, 2.2);
			double interp3p3 = Interpolate.findLogLogY(xs, ys, 3.3);
			yDest.set(18, interp2p2);
			yDest.set(19, interp3p3);
		}

		CurveContainer nshmpCC = CurveContainer.create(nshmpRegion,
			Doubles.asList(xPGA));
		nshmpCC.union(pgaCC);
		
		return nshmpCC;
	}
	
}
