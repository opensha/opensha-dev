package org.opensha.nshmp2.calc;

import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.nshmp2.util.Period;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.io.Closeables;
import com.google.common.io.Files;
import com.google.common.io.Flushables;

/**
 * Writer of hazard result
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class HazardResultWriterLocal implements HazardResultWriter {

	private static final Joiner J = Joiner.on(',').useForNull(" ");
	private BufferedWriter writer;

	/**
	 * Creates anew local writer instance.
	 * @param outFile output location
	 * @param period
	 * @throws IOException
	 */
	public HazardResultWriterLocal(File outFile, Period period)
			throws IOException {
		Files.createParentDirs(outFile);
		writer = Files.newWriter(outFile, Charsets.US_ASCII);
		writeCurveHeader(writer, period);
	}

	@Override
	public void write(HazardResult result) throws IOException {
		String resultStr = formatResult(result);
		writer.write(resultStr);
		writer.newLine();
	}

	@Override
	public void close() throws IOException {
		Flushables.flushQuietly(writer);
		Closeables.close(writer, true);
	}

	private static String formatResult(HazardResult result) {
		List<String> dat = Lists.newArrayList();
		Location loc = result.location();
		dat.add(Double.toString(loc.getLatitude()));
		dat.add(Double.toString(loc.getLongitude()));
		for (Point2D p : result.curve()) {
			dat.add(Double.toString(p.getY()));
		}
		return J.join(dat);
	}

	/**
	 * Utility method to write header for a set of curves into a {@code Writer}.
	 * @param writer to add header to
	 * @param period for reference x (ground motion) values
	 * @throws IOException
	 */
	public static void writeCurveHeader(BufferedWriter writer, Period period)
			throws IOException {
		List<String> headerVals = Lists.newArrayList();
		headerVals.add("lat");
		headerVals.add("lon");
		for (Double d : period.getIMLs()) {
			headerVals.add(d.toString());
		}
		writer.write(J.join(headerVals));
		writer.newLine();
	}

}
