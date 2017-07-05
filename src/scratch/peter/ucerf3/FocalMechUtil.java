package scratch.peter.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.IOUtils;

import com.google.common.base.CharMatcher;
import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

/**
 * Used to clean up some of D. Jacksons focal mech files. Compacts rotation
 * uncertainty file to 3 columns.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class FocalMechUtil {

	static String path = "/Users/pmpowers/Documents/UCERF3/seismicity/focalMechs/data";
	static String fileIn = "HAZTBL_eth_seism.txt";
	static String fileOut = "FocalMechs_Uncert_DJ_2-14-2012.txt";
	// static String fileOut = "FocalMechs_SDR_DJ_2-14-2012";

	static Joiner join = Joiner.on(" ");
	static Splitter split = Splitter.on(CharMatcher.WHITESPACE).omitEmptyStrings();

	public static void main(String[] args) {
		List<String> lines = null;
		try {
			File f = new File(path, fileIn);
			lines = Files.readLines(f, Charsets.US_ASCII);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		List<String> newLines = Lists.newArrayList();
		for (String line : lines) {
			if (line.startsWith("#")) continue;
			if (line.isEmpty()) continue;
			Iterable<String> parts = split.split(line);
//			System.out.println(Iterables.toString(parts));
			String lon = Iterables.get(parts, 0);
			String lat = Iterables.get(parts, 1);
			String unc = Iterables.get(parts, 10);
			String newLine = join.join(lat, lon, unc);
			newLines.add(newLine);
		}

		try {
			File f = new File(path, fileOut);
			Files.write(Joiner.on(IOUtils.LINE_SEPARATOR).join(newLines), f,
				Charsets.US_ASCII);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
}
