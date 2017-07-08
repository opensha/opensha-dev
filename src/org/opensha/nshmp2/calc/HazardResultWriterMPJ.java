package org.opensha.nshmp2.calc;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp2.util.Period;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.io.Closeables;
import com.google.common.io.Files;
import com.google.common.io.Flushables;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class HazardResultWriterMPJ implements HazardResultWriter {

	private static String format = "%.3f";
	private static final Joiner J = Joiner.on(',').useForNull(" ");

	private File outDir;

	/**
	 * Creates anew local writer instance.
	 * @param outDir output location
	 * @throws IOException
	 */
	public HazardResultWriterMPJ(File outDir) throws IOException {
		this.outDir = outDir;
		outDir.mkdirs();
	}

	@Override
	public void write(HazardResult result) throws IOException {
		String fName = createFileName(result.location());
		String resultStr = formatResult(result.curve());
		File file = new File(outDir, fName);
		Files.write(resultStr, file, Charsets.US_ASCII);
	}

	@Override
	public void close() throws IOException {
		// do nothing
	}
	
	private String createFileName(Location loc) {
		return String.format(format, loc.getLatitude()) + "_" +
			String.format(format, loc.getLongitude()) + ".txt";
	}
	
	private static String formatResult(DiscretizedFunc curve) {
		List<String> dat = Lists.newArrayList();
		for (Point2D p : curve) {
			dat.add(Double.toString(p.getY()));
		}
		return J.join(dat);
	}

}
