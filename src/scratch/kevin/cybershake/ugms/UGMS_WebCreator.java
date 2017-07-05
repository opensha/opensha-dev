package scratch.kevin.cybershake.ugms;

import java.io.File;
import java.io.FileWriter;
import java.net.URL;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.text.WordUtils;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM_Vs30;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.siteData.impl.CVM_Vs30.CVM;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.SiteInfo2DB;

import com.google.common.base.CaseFormat;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class UGMS_WebCreator {
	
	private static final String GMAPS_API_KEY = "AIzaSyCOfe8NIHLR0Z6l4KzajcDAwxOjlhLlEb4";

	public static void main(String[] args) {
		String prefix = "study_15_4";
		File mainDir = new File("/home/kevin/CyberShake/MCER/mcer_data_products");
		
		boolean reDownloadMaps = false;
		
		boolean[] finalMCER = { true, false };
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		try {
			SiteInfo2DB sites2db = new SiteInfo2DB(db);
			
			List<SiteData<Double>> siteDataProvs = Lists.newArrayList();
			List<String> siteDataProvNames = Lists.newArrayList();
			siteDataProvs.add(new WillsMap2006());
			siteDataProvNames.add("Wills 2006 Vs30 (m/s)");
			siteDataProvs.add(new CVM_Vs30(CVM.CVMS4i26));
			siteDataProvNames.add("CVMS-4.26 Vs30 (m/s)");
			siteDataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
			siteDataProvNames.add("CVMS-4.26 Z1.0 (km)");
			siteDataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			siteDataProvNames.add("CVMS-4.26 Z2.5 (km)");
			
			for (boolean myFinalMCER : finalMCER) {
				for (File setDir : mainDir.listFiles()) {
					String name = setDir.getName();
					if (!setDir.isDirectory() || !name.startsWith(prefix))
						continue;
					
					String formattedName = name;
					// add periods to study names
					for (int i=1; i<formattedName.length()-1; i++) {
						if (formattedName.charAt(i) != '_')
							continue;
						if (Character.isDigit(formattedName.charAt(i-1)) && Character.isDigit(formattedName.charAt(i+1)))
							formattedName = formattedName.substring(0,i)+'.'+formattedName.substring(i+1);
					}
					formattedName = formattedName.replaceAll("_", " ");
					formattedName = WordUtils.capitalize(formattedName);
					
					System.out.println("Handling "+formattedName+" ("+name+")");
					
					List<String> lines = Lists.newArrayList();
					
					lines.add("<h1>"+formattedName+"</h1>");
					lines.add("<i>Note: quickly navigate sites with your keyboard left/right arrows</i>");
					
					File[] subDirs = setDir.listFiles();
					Arrays.sort(subDirs, new FileNameComparator());
					
					List<File> runDirs = Lists.newArrayList();
					List<String> siteNames = Lists.newArrayList();
					
					for (File runDir : subDirs) {
						String runDirName = runDir.getName();
						if (!runDir.isDirectory() || !runDirName.contains("_run"))
							continue;
						System.out.println("Found run: "+runDirName);
						
						runDirs.add(runDir);
						siteNames.add(runDirName.substring(0, runDirName.indexOf("_run")));
					}
					
					lines.add("<h2>Sites</h2>");
					lines.add("<ul>");
					for (int i=0; i<runDirs.size(); i++) {
						String dirName = runDirs.get(i).getName();
						String siteName = siteNames.get(i);
						String dataLink = "(<a href=\""+dirName+"\">all plots</a>)";
						lines.add("<li><a href=\"#"+dirName+"\">"+siteName+"</a> "+dataLink+"</li>");
					}
					lines.add("</ul>");
					
					for (int i=0; i<runDirs.size(); i++) {
						String siteName = siteNames.get(i);
						File runDir = runDirs.get(i);
						File mcerFile = new File(runDir, runDir.getName()+"_RotD100_MCER.png");
						File mcerFinalFile = new File(runDir, runDir.getName()+"_RotD100_final_MCER.png");
						Preconditions.checkState(mcerFile.exists());
						
						String mcerImage;
						if (myFinalMCER) {
							mcerImage = "<img src=\""+runDir.getName()+"/"+mcerFinalFile.getName()+"\" width=\"50%\" height=\"50%\"/>"
									+ "<img src=\""+runDir.getName()+"/"+mcerFile.getName()+"\" width=\"50%\" height=\"50%\"/>";
						} else {
							mcerImage = "<img src=\""+runDir.getName()+"/"+mcerFile.getName()+"\"/>";
						}
						Location loc = sites2db.getSiteFromDB(siteName).createLocation();
						File mapFile = new File(runDir, "site_location_map.png");
						if (!mapFile.exists() || reDownloadMaps)
							FileUtils.downloadURL(new URL(getMiniMap(loc)), mapFile);
						String mapImage = "<img src=\""+runDir.getName()+"/"+mapFile.getName()+"\"/>";
						
						String metadata = "<h3><ul>";
						metadata += "<li><a href=\""+runDir.getName()+"\">All Plots</a></li>";
						metadata += "<li><a href=\""+runDir.getName()+"/disaggregations\">Disaggregations</a></li>";
						metadata += "<li>Location: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"</li>";
						for (int s=0; s<siteDataProvs.size(); s++) {
							SiteData<Double> prov = siteDataProvs.get(s);
							double val = prov.getValue(loc);
							String valStr;
							if (prov.getDataType().equals(SiteData.TYPE_VS30))
								valStr = (int)(val+0.5d)+"";
							else
								valStr = (float)val+"";
							metadata += "<li>"+siteDataProvNames.get(s)+": "+valStr+"</li>";
						}
						metadata += "</ul></h3>";
						
						lines.add("<a name=\""+runDir.getName()+"\"/><h2>"+siteName+"</h2> "
								+"(<a href=\"#top\">top</a>) (<a href=\"../\">parent</a>)");
						lines.add("<table>");
						lines.add("  <tr>");
						lines.add("    <td rowspan=\"2\">"+mcerImage+"</td>");
						lines.add("    <td>"+mapImage+"</td>");
						lines.add("  </tr>");
						lines.add("  <tr>");
						lines.add("    <td>"+metadata+"</td>");
						lines.add("  </tr>");
						lines.add("</table>");
					}
					
					File htmlFile;
					if (myFinalMCER)
						htmlFile = new File(setDir, "final.html");
					else
						htmlFile = new File(setDir, "index.html");
					FileWriter fw = new FileWriter(htmlFile);
					
					for (String line : lines)
						fw.write(line+"\n");
					
					fw.write(getJSKeyScript());
					
					fw.close();
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			db.destroy();
		}
	}
	
	private static String getMiniMap(Location loc) {
		String locStr = loc.getLatitude()+","+loc.getLongitude();
		return "https://maps.googleapis.com/maps/api/staticmap?center="+locStr+"&zoom=9"
				+ "&size=400x300&maptype=roadmap&markers="+locStr+"&key="+GMAPS_API_KEY;
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
