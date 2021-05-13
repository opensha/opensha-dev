package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.DiscretizedFunc;

import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

public class BBP_PartBSimZipLoader extends BBP_SimZipLoader {
	
	private Map<String, DiscretizedFunc[]> rdMap;

	public BBP_PartBSimZipLoader(File file, int numSites) throws ZipException, IOException {
		this(new ZipFile(file), numSites);
	}
	
	private static List<BBP_Site> buildSites(int numSites) {
		List<BBP_Site> sites = new ArrayList<>();
		
		for (int i=0; i<numSites; i++)
			// only name matters here
			sites.add(new BBP_Site("s"+i, null, Double.NaN, Double.NaN, Double.NaN));
		
		return sites;
	}

	public BBP_PartBSimZipLoader(ZipFile zip, int numSites) {
		super(zip, buildSites(numSites));
		rdMap = new HashMap<>();
	}
	
	public Collection<Integer> getEventsIDs(Scenario scenario) {
		HashSet<Integer> ids = new HashSet<>();

		String prefix = scenario.getPrefix();
		Table<BBP_Site, String, Map<FileType, ZipEntry>> table = getEntriesTable();
		for (String dirName : table.columnKeySet()) {
			if (dirName.contains(prefix) && dirName.contains("event_")) {
				String eventStr = dirName.substring(dirName.indexOf("event_"));
				eventStr = eventStr.substring(eventStr.indexOf("_")+1);
				eventStr = eventStr.substring(0, eventStr.indexOf("_"));
				Integer id = Integer.parseInt(eventStr);
				ids.add(id);
			}
		}
		
		return ids;
	}
	
	public boolean hasScenario(Scenario scenario) {
		String prefix = scenario.getPrefix();
		Table<BBP_Site, String, Map<FileType, ZipEntry>> table = getEntriesTable();
		for (String dirName : table.columnKeySet())
			if (dirName.contains(prefix))
				return true;
		return false;
	}
	
	public DiscretizedFunc[] getRotD50(Integer eventID, Scenario scenario, double distance) throws IOException {
		String dirName = MPJ_BBP_PartBSim.getDirName(scenario, eventID, distance);
		if (rdMap.containsKey(dirName))
			return rdMap.get(dirName);
		DiscretizedFunc[] ret = new DiscretizedFunc[sites.size()];
		for (int i=0; i<sites.size(); i++) {
			ret[i] = readRotD50(sites.get(i), dirName);
		}
		rdMap.putIfAbsent(dirName, ret);
		return ret;
	}

}
