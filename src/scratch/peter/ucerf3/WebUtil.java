package scratch.peter.ucerf3;

import static com.google.common.base.Charsets.US_ASCII;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.text.WordUtils;
import org.opensha.commons.geo.Location;

import com.google.common.base.Splitter;
import com.google.common.base.Strings;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

import scratch.peter.ucerf3.calc.UC3_CalcUtils;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class WebUtil {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO do nothing
		addHeaders();

	}
	
	
	/**
	 * Adds HEADER.html to each UCERF3 hazard analysis site directory
	 */
	public static void addHeaders() throws IOException {
		String plotDir = "/Users/pmpowers/projects/UCERF3.3/HazardCurves/Sites/sites/";
		String sitesPath = "tmp/UC33/curvejobs/sites/all.txt";
		File readMe = new File(plotDir + "README.html");
		File header = new File(plotDir + "HEADER.html");
		String headerStr = Files.readFirstLine(header, US_ASCII);

		File f = new File(sitesPath);
		List<String> lines = Files.readLines(f, US_ASCII);
		boolean pbr = false;
		for (String line : lines) {
			if (line.startsWith("# PBR")) pbr = true;
			if (line.startsWith("#")) continue;
			Iterator<String> it = Splitter.on(',').split(line).iterator();
			String name = it.next();
			double lat = Double.parseDouble(it.next());
			double lon = Double.parseDouble(it.next());
			String displayName = pbr ? name : displayName(name);
			String title = String.format(displayName + " (%.2f, %.2f)", lat, lon);
			
			String newHeaderStr = headerStr.replace("####",title);
			
			File siteReadMe = new File(plotDir + name + "/README.html");
			File siteHeader = new File(plotDir + name + "/HEADER.html");
			Files.copy(readMe, siteReadMe);
			Files.write(newHeaderStr, siteHeader, US_ASCII);
		}
	}
	
	private static String displayName(String name) {
		return WordUtils.capitalizeFully(name.replace('_', ' ').toLowerCase());
	}

}
