package scratch.peter.nshmp;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.calc.HazardResultWriter;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.base.StandardSystemProperty;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

/**
 * Writier for deterministic hazard results for CA. Recently changed to aggregate 
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class HazardResultWriterMPJ_NSHMP_Det implements HazardResultWriter {

	private static String format = "%.2f";
	private static final Joiner Jt = Joiner.on('\t').useForNull(" ");
	private static final String LF = StandardSystemProperty.LINE_SEPARATOR.value();
	
	private List<Location> locs;
	private List<DeterministicResult> detData;

	/**
	 * Creates anew local writer instance.
	 * @throws IOException
	 */
	public HazardResultWriterMPJ_NSHMP_Det() throws IOException {
		locs = Lists.newArrayList();
		detData = Lists.newArrayList();
	}

	@Override
	public void write(HazardResult result) throws IOException {
		locs.add(result.location());
		detData.add(result.detData());
	}

	/**
	 * Write collected resultes to file identified by the node 'rank'. 
	 * @param outDir
	 * @param id
	 * @throws IOException
	 */
	public void toFile(File outDir, int id) throws IOException {
		String detName = id + "_det.txt";
		File detFile = new File(outDir, detName);
		Files.write("", detFile, Charsets.US_ASCII);
		
		for (int i=0; i<locs.size(); i++) {
			String result = formatResult(locs.get(i), detData.get(i)) + LF;
			Files.append(result, detFile, Charsets.US_ASCII);
		}
	}
	
	@Override
	public void close() throws IOException {
		// do nothing
	}
	
	private static String formatResult(Location loc, DeterministicResult detData) {
		return Jt.join(
			String.format(format, loc.getLatitude()),
			String.format(format, loc.getLongitude()),
			detData.median,
			detData.mag,
			detData.rRup,
			detData.name);
	}
}
