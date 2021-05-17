package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.text.WordUtils;
import org.opensha.commons.util.FileNameComparator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class PowellMeetingWebOrganizer {

	public static void main(String[] args) throws IOException {
		File etasDir = new File("/var/www/html/ftp/kmilner/ucerf3/etas_results");
		File combLossDir = new File(etasDir, "losses_combined");
		File outputMainDir = new File(etasDir, "2017_04-powell-center");
		File griddedNuclDir = new File(new File(etasDir, "science_paper"), "figure_6");
		String mapDirPrefix = "2017_03_23";
		String title = "UCERF3-ETAS Example Products for April 2017 Powell Center Meeting";
		
		List<String> scenarios = Lists.newArrayList();
		List<String> combinedMainDirNames = Lists.newArrayList();
		List<String> combinedLossDirNames = Lists.newArrayList();
		
		scenarios.add("Haywired M7.1, Fault Based ETAS");
		combinedMainDirNames.add("2016_06_15-haywired_m7-10yr-full_td-no_ert-combined");
		combinedLossDirNames.add("haywired_100k_both_models");
		
		scenarios.add("Haywired M7.1, Pure ETAS (no faults)");
		combinedMainDirNames.add("2017_01_02-haywired_m7-10yr-gridded-only-200kcombined");
		combinedLossDirNames.add(null);
		
		scenarios.add("Northridge M6.7, Fault Based ETAS");
		combinedMainDirNames.add("2017_02_01-northridge-m6.7-10yr-full_td-no_ert-combined");
		combinedLossDirNames.add(null);
		
		scenarios.add("Northridge M6.7, Pure ETAS (no faults)");
		combinedMainDirNames.add("2017_02_01-northridge-m6.7-10yr-gridded-only-combined200k");
		combinedLossDirNames.add(null);
		
		scenarios.add("Mojave M7, Fault Based ETAS");
		combinedMainDirNames.add("2016_02_22-mojave_m7-10yr-full_td-no_ert-combined");
		combinedLossDirNames.add("mojave_100k_both_models");
		
		scenarios.add("Mojave M7, Pure ETAS (no faults)");
		combinedMainDirNames.add("2016_12_03-mojave_m7-10yr-gridded-only");
		combinedLossDirNames.add(null);
		
		scenarios.add("2016 Bombay Beach Swarm, Fault Based ETAS");
		combinedMainDirNames.add("2016_10_27-2016_bombay_swarm-10yr-full_td-no_ert-combined");
		combinedLossDirNames.add("bombay_beach_2016_swarm_100k_both_models");
		
		List<String> anchorNames = Lists.newArrayList();
		for (String name : scenarios)
			anchorNames.add(name.replaceAll(" ", "_").replaceAll("\\(", "").replaceAll("\\)", "")
					.replaceAll(",", "").replaceAll("\\.", ""));
		
		FileWriter fw = new FileWriter(new File(outputMainDir, "index.html"));
		
		fw.write("<!DOCTYPE html>\n");
		fw.write("<html>\n");
		fw.write("\n");
		fw.write("<head>\n");
		fw.write("  <title>"+title+"</title>\n");
		fw.write("</head>\n");
		fw.write("\n");
		fw.write("<body>\n");
		fw.write("<style>\n");
		fw.write("table {\n");
		fw.write("    border-collapse: collapse;\n");
		fw.write("    width: 100%;\n");
		fw.write("}\n");
		fw.write("\n");
		fw.write("td, th {\n");
		fw.write("    border: 1px solid #ddd;\n");
		fw.write("    padding: 8px;\n");
		fw.write("}\n");
		fw.write("\n");
		fw.write("tr:nth-child(even){background-color: #f2f2f2;}\n");
		fw.write("\n");
		fw.write("tr:hover {background-color: #ddd;}\n");
		fw.write("\n");
		fw.write("th {\n");
		fw.write("    padding-top: 12px;\n");
		fw.write("    padding-bottom: 12px;\n");
		fw.write("    text-align: left;\n");
		fw.write("    background-color: #4CAF50;\n");
		fw.write("    color: white;\n");
		fw.write("}\n");
		fw.write("</style>\n");
		fw.write("\n");
		fw.write("<h1>"+title+"</h1>\n");
		fw.write("\n");
		
		// scenario list
		fw.write("<h2>Available Scenarios:</h2>\n");
		fw.write("<ul>\n");
		for (int i=0; i<scenarios.size(); i++)
			fw.write("<li><b><a href=\"#"+anchorNames.get(i)+"\">"+scenarios.get(i)+"</a></b></li>\n");
		fw.write("</ul>\n");
		
		for (int i=0; i<scenarios.size(); i++) {
			String scenario = scenarios.get(i);
			fw.write("<a name=\""+anchorNames.get(i)+"\"/><h2>"+scenario+" <i>(<a href=\"#top\">top</a>)</i></h2>\n");
			File mainDir = new File(etasDir, combinedMainDirNames.get(i));
			File plotsDir = new File(mainDir, "plots");
			File lossDir = null;
			if (combinedLossDirNames.get(i) != null)
				lossDir = new File(combLossDir, combinedLossDirNames.get(i));
			
			fw.write("<table>\n");
			
			writeTableLine(fw, "Product", "Links", true);
			
			String left = "Magnitude Probability Distributions";
			String fullMFDs = getLinks(plotsDir, "consolidated_aftershocks_mag_num_cumulative");
			String weekMFDs = getLinks(plotsDir, "one_week_consolidated_aftershocks_mag_num_cumulative");
			String right = null;
			if (!weekMFDs.isEmpty())
				right = "1 Week: "+weekMFDs;
			if (!fullMFDs.isEmpty()) {
				if (right == null)
					right = "";
				else
					right += "<br>\n";
				right += "10 years: "+fullMFDs;
			}
			writeTableLine(fw, left, right, false);
			
			left = "Magnitude Probability Maps";
			right = null;
			File png = null;
			File pdf = null;
			String scenarioPrefix = scenario.substring(0, scenario.indexOf(" ")).toLowerCase();
			for (File file : griddedNuclDir.listFiles()) {
				String name = file.getName();
				if (scenario.contains("Fault Based")) {
					if (name.toLowerCase().contains("gridded-only"))
						continue;
				} else {
					if (!name.toLowerCase().contains("gridded-only"))
						continue;
				}
				if (name.startsWith(scenarioPrefix) && name.endsWith("1wk_m2.5.png"))
					png = file;
				if (name.startsWith(scenarioPrefix) && name.endsWith("1wk_m2.5.pdf"))
					pdf = file;
			}
			if (png != null)
				right = "1 Week, M>=2.5: "+getLink(png, "[PNG]")+" "+getLink(pdf, "[PDF]");
			writeTableLine(fw, left, right, false);
			
			File curveDir = new File(mainDir, "hazard_maps/curves");
			String[] hazardDurations = { "day", "days_0_3", "week", "month", "year", "full" };
			String[] hazardDurationLabels = { "1 Day", "3 Day", "1 Week", "1 Month", "1 Year", "10 Years" };
			left = "Hazard Curves";
			right = "";
			String[] imtLabels = { "PGA", "PGV", "MMI" };
			if (curveDir.exists()) {
				List<String> curveLines = Lists.newArrayList();
				File[] subDirs = curveDir.listFiles();
				Arrays.sort(subDirs, new FileNameComparator());
				for (File dir : subDirs) {
					
					if (!dir.isDirectory() || dir.getName().startsWith("."))
						continue;
					for (String imt : imtLabels) {
						String name = dir.getName();
						int index = name.indexOf(imt.toLowerCase());
						if (index < 0)
							continue;
						name = name.replaceAll("_", " ").substring(0, index).trim();
						name = WordUtils.capitalize(name);
						String line = name+", "+imt+":";
						for (int d=0; d<hazardDurations.length; d++) {
							String label = hazardDurationLabels[d];
							File map = getPlot(dir, dir.getName(), hazardDurations[d]+".pdf");
							line += " "+getLink(map, label);
							if (hazardDurations[d].equals("days_0_3")) {
								File anim = getPlot(dir.getParentFile(), dir.getName(), ".gif");
								line += " "+getLink(anim, "3 Day Animation");
							}
						}
						curveLines.add(line);
					}
				}
				curveLines.add("<font size=\"2\">"+getLink(curveDir, "Browse all available curves")+"</font>");
				right = Joiner.on("<br>\n").join(curveLines);
			}
			writeTableLine(fw, left, right, false);
			
			left = "Hazard Maps";
			right = "";
			String[] poePrefixes = { "map_mmi_poe_6", "map_mmi_poe_8" };
			String[] poePrefixLabels = { "POE MMI 6", "POE MMI 8" };
			String[] mmiPrefixes = { "map_mmi_p0.5", "map_mmi_p0.25", "map_mmi_p0.1", "map_mmi_p0.01", "map_mmi_p0.001", };
			String[] mmiPrefixLabels = { "MMI @ 50% POE", "MMI @ 25% POE", "MMI @ 10% POE", "MMI @ 1% POE", "MMI @ 0.1% POE", };
			String[] gainPrefixes = { "gain_map_mmi_poe_6", "gain_map_mmi_poe_8" };
			String[] gainPrefixLabels = { "POE MMI 6 Gain", "POE MMI 8 Gain" };
			File mapDir = findMapDir(mainDir, mapDirPrefix);
			if (mapDir != null) {
				List<String> mapLines = Lists.newArrayList();
				mapLines.add("<b>Probability of Exceeding MMI</b>");
				boolean hasGain = false;
				for (int p=0; p<poePrefixes.length; p++) {
					String prefix = poePrefixes[p];
					String line = poePrefixLabels[p]+":";
					for (int d=0; d<hazardDurations.length; d++) {
						File durationDir = new File(mapDir, hazardDurations[d]);
						hasGain = hasGain || hasMatch(durationDir, gainPrefixes[0]);
						String label = hazardDurationLabels[d];
						File map = getPlot(durationDir, prefix, ".pdf");
						line += " "+getLink(map, label);
						if (hazardDurations[d].equals("days_0_3")) {
							File anim = getPlot(mapDir, prefix, ".gif");
							line += " "+getLink(anim, "3 Day Animation");
						}
					}
					mapLines.add(line);
				}
				if (hasGain) {
					mapLines.add("<b>MMI POE GAIN to Long Term UCERF3 Models</b>");
					for (boolean td : new boolean[] { true, false }) {
						for (int p=0; p<poePrefixes.length; p++) {
							String prefix = gainPrefixes[p];
							String line = gainPrefixLabels[p];
							String extra;
							if (td) {
								line += " UCERF3-TD:";
								extra = "u3-td";
							} else {
								line += " UCERF3-TI:";
								extra = "u3-ti";
							}
							for (int d=0; d<hazardDurations.length; d++) {
								File durationDir = new File(mapDir, hazardDurations[d]);
								String label = hazardDurationLabels[d];
								File map = getPlot(durationDir, prefix, ".pdf", extra);
								line += " "+getLink(map, label);
								if (hazardDurations[d].equals("days_0_3")) {
									File anim = getPlot(mapDir, prefix, ".gif", extra);
									line += " "+getLink(anim, "3 Day Animation");
								}
							}
							mapLines.add(line);
						}
					}
				}
				mapLines.add("<b>MMI with fixed Probability of Exceedance</b>");
				for (int p=0; p<mmiPrefixes.length; p++) {
					String prefix = mmiPrefixes[p];
					String line = mmiPrefixLabels[p]+":";
					for (int d=0; d<hazardDurations.length; d++) {
						File durationDir = new File(mapDir, hazardDurations[d]);
						String label = hazardDurationLabels[d];
						File map = getPlot(durationDir, prefix, ".pdf");
						line += " "+getLink(map, label);
						if (hazardDurations[d].equals("days_0_3")) {
							File anim = getPlot(mapDir, prefix, ".gif");
							line += " "+getLink(anim, "3 Day Animation");
						}
					}
					mapLines.add(line);
				}
				mapLines.add("<font size=\"2\">"+getLink(mapDir, "Browse all available maps")+"</font>");
				right = Joiner.on("<br>\n").join(mapLines);
			}
			writeTableLine(fw, left, right, false);
			
			left = "Loss Exceedance Curves";
			right = "";
			if (lossDir != null) {
				String[] lossDurations = { "1day_logY.pdf", "1wk_logY.pdf", "1mo_logY.pdf", "1yr_logY.pdf",
						"10yr_logY.pdf", "exceed_logY.pdf" };
				String[] lossDurationLabels = { "1 Day", "1 Week", "1 Month", "1 Year", "10 Years", "All" };
				for (int d=0; d<lossDurations.length; d++) {
					File plot = getPlot(lossDir, scenarioPrefix, lossDurations[d]);
					if (plot == null || !plot.exists())
						continue;
					if (!right.isEmpty())
						right += " ";
					right += getLink(plot, lossDurationLabels[d]);
				}
			}
			writeTableLine(fw, left, right, false);
			
			left = "Fault Participation Probabilites<br><font size=\"2\">NOTE: Annualized rates, across 10 years.<br>"
					+ "Multiply by 10 for total occurrence rate in simulation duration.</font>";
			String supraPlots = getLinks(plotsDir, "full_children_sect_partic");
			String m6p7Plots = getLinks(plotsDir, "full_children_sect_partic_m6.7");
			String m7p8Plots = getLinks(plotsDir, "full_children_sect_partic_m7.8");
			right = null;
			if (!supraPlots.isEmpty())
				right = "Supra Seis: "+supraPlots;
			if (!m6p7Plots.isEmpty()) {
				if (right == null)
					right = "";
				else
					right += "<br>\n";
				right += "M>=6.7: "+m6p7Plots;
			}
			if (!m7p8Plots.isEmpty()) {
				if (right == null)
					right = "";
				else
					right += "<br>\n";
				right += "M>=7.8: "+m7p8Plots;
			}
			writeTableLine(fw, left, right, false);
			
			fw.write("</table>\n");
		}
		fw.write("</body>\n");
		fw.write("</html>\n");
		fw.close();
	}
	
	private static boolean hasMatch(File dir, String contains) {
		for (File file : dir.listFiles())
			if (file.isFile() && file.getName().contains(contains))
				return true;
		return false;
	}
	
	private static File getPlot(File dir, String prefix, String ext, String... extras) {
		if (!dir.exists())
			return null;
		fileLoop:
		for (File file : dir.listFiles()) {
			if (!file.isFile())
				continue;
			String name = file.getName();
			for (String extra : extras)
				if (!name.contains(extra))
					continue fileLoop;
			if (name.startsWith(prefix) && name.endsWith(ext))
				return file;
		}
		return null;
	}
	
	private static File findMapDir(File baseDir, String mapPrefix) {
		File mapDir = new File(baseDir, "hazard_maps");
		if (!mapDir.exists())
			return null;
		for (File subDir : mapDir.listFiles()) {
			if (subDir.isDirectory() && subDir.getName().startsWith(mapPrefix)) {
				File combDir = new File(subDir, "maps_combined");
				if (combDir.exists())
					return combDir;
				File gridDir = new File(subDir, "maps_gridded");
				if (gridDir.exists())
					return gridDir;
			}
		}
		return null;
	}
	
	private static void writeTableLine(FileWriter fw, String left, String right, boolean header) throws IOException {
		String tag;
		if (header)
			tag = "th";
		else
			tag = "td";
		if (right == null)
			right = "";
		fw.write("  <tr>\n");
		fw.write("    <"+tag+">"+left+"</"+tag+">\n");
		fw.write("    <"+tag+">"+right+"</"+tag+">\n");
		fw.write("  </tr>\n");
	}
	
	private static String getLinks(File dir, String prefix) {
		File png = new File(dir, prefix+".png");
		File pdf = new File(dir, prefix+".pdf");
		if (!png.exists())
			return "";
		Preconditions.checkState(pdf.exists());
		return getLink(png, "[PNG]")+" "+getLink(pdf, "[PDF]");
	}
	
	private static String getLink(File file, String text) {
		if (file == null || !file.exists())
			return "";
		String url = file.getAbsolutePath().replace("/var/www/html", "http://opensha.usc.edu");
		return "<a href=\""+url+"\">"+text+"</a>";
	}
	
	private static String getJSKeyScript() {
		String script = "<script type=\"text/javascript\">";
		script += "\n";
		script += "\ndocument.onkeydown = checkKey;";
		script += "\nfunction checkKey(e) {";
		script += "\n    e = e || window.event;";
		script += "\n    if (e.altKey) {";
		script += "\n      return;";
		script += "\n    }";
		script += "\n    if (e.keyCode == '37') {";
		script += "\n        // left arrow";
		script += "\n       goToPrevAnchor()";
		script += "\n    }";
		script += "\n    if (e.keyCode == '39') {";
		script += "\n        // right arrow";
		script += "\n       goToNextAnchor()";
		script += "\n    }";
		script += "\n}";
		script += "\n";
		// next anchor
		script += "\nfunction goToNextAnchor() {";
		script += "\nvar anchors = document.anchors;";
		script += "\nvar loc = window.location.href.replace(/#.*/,'');";
		script += "\nvar nextAnchorName;";
		script += "\n";
		script += "\n// Get name of the current anchor from the hash";
		script += "\n// if there is one";
		script += "\nvar anchorName = window.location.hash.replace(/#/,'');";
		script += "\n";
		script += "\n// If there is an anchor name...";
		script += "\nif (anchorName) {";
		script += "\n";
		script += "\n  // Find current element in anchor list, then";
		script += "\n  // get next anchor name, or if at last anchor, set to first";
		script += "\n  for (var i=0, iLen=anchors.length; i<iLen; i++) {";
		script += "\n    if (anchors[i].name == anchorName) {";
		script += "\n      nextAnchorName = anchors[++i % iLen].name;";
		script += "\n      break;";
		script += "\n    }";
		script += "\n  }";
		script += "\n}";
		script += "\n";
		script += "\n// If there was no anchorName or no match,";
		script += "\n// set nextAnchorName to first anchor name";
		script += "\nif (!nextAnchorName) {";
		script += "\n  nextAnchorName = anchors[0].name;";
		script += "\n}";
		script += "\n";
		script += "\n// Go to new URL";
		script += "\nwindow.location.href = loc + '#' + nextAnchorName;";
		script += "\n}";
		// prev anchor
		script += "\nfunction goToPrevAnchor() {";
		script += "\nvar anchors = document.anchors;";
		script += "\nvar loc = window.location.href.replace(/#.*/,'');";
		script += "\nvar nextAnchorName;";
		script += "\n";
		script += "\n// Get name of the current anchor from the hash";
		script += "\n// if there is one";
		script += "\nvar anchorName = window.location.hash.replace(/#/,'');";
		script += "\n";
		script += "\n// If there is an anchor name...";
		script += "\nif (anchorName) {";
		script += "\n";
		script += "\n  // Find current element in anchor list, then";
		script += "\n  // get next anchor name, or if at last anchor, set to first";
		script += "\n  for (var i=0, iLen=anchors.length; i<iLen; i++) {";
		script += "\n    if (anchors[i].name == anchorName) {";
		script += "\n      nextAnchorName = anchors[--i % iLen].name;";
		script += "\n      break;";
		script += "\n    }";
		script += "\n  }";
		script += "\n}";
		script += "\n";
		script += "\n// If there was no anchorName or no match,";
		script += "\n// set nextAnchorName to first anchor name";
		script += "\nif (!nextAnchorName) {";
		script += "\n  nextAnchorName = anchors[anchors.length - 1].name;";
		script += "\n}";
		script += "\n";
		script += "\n// Go to new URL";
		script += "\nwindow.location.href = loc + '#' + nextAnchorName;";
		script += "\n}";
		script += "\n</script>";
		return script;
	}

}
