package scratch.peter.nshmp;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.calc.HazardResultWriter;
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
public class HazardResultWriterMPJ_NSHMP implements HazardResultWriter {

	private static String format = "%.3f";
	private static final Joiner Jc = Joiner.on(',').useForNull(" ");
	private static final Joiner Jt = Joiner.on('\t').useForNull(" ");

	private File outDir;

	/**
	 * Creates anew local writer instance.
	 * @param outDir output location
	 * @throws IOException
	 */
	public HazardResultWriterMPJ_NSHMP(File outDir) throws IOException {
		this.outDir = outDir;
		outDir.mkdirs();
	}

	@Override
	public void write(HazardResult result) throws IOException {
		String probName = createFileName(result.location()) + "_prob.txt";
		String probStr = formatResult(result.curve());
		File probFile = new File(outDir, probName);
		Files.write(probStr, probFile, Charsets.US_ASCII);

		String detName = createFileName(result.location()) +  "_det.txt";
		String detStr = formatResult(result.detData());
		File detFile = new File(outDir, detName);
		Files.write(detStr, detFile, Charsets.US_ASCII);
	}

	@Override
	public void close() throws IOException {
		// do nothing
	}
	
	private String createFileName(Location loc) {
		return String.format(format, loc.getLatitude()) + "_" +
			String.format(format, loc.getLongitude());
	}
	
	private static String formatResult(DiscretizedFunc curve) {
		List<String> dat = Lists.newArrayList();
		for (Point2D p : curve) {
			dat.add(Double.toString(p.getY()));
		}
		return Jc.join(dat);
	}

	private static String formatResult(DeterministicResult detData) {
		return Jt.join(detData.median, detData.mag, detData.rRup, detData.name);
	}
}
