package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.DiscretizedFunc;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.ruptures.BBP_CatalogPartBValidationConfig.Scenario;

public class BBP_PartBSimZipLoader extends BBP_SimZipLoader {
	
	private List<BBP_Site> sites;
	
	private Map<String, DiscretizedFunc[]> rdMap;

	public BBP_PartBSimZipLoader(File file, Scenario[] scenarios, int numSites) throws ZipException, IOException {
		this(new ZipFile(file), scenarios, numSites);
	}
	
	private static List<BBP_Site> buildSites(int numSites) {
		List<BBP_Site> sites = new ArrayList<>();
		
		for (int i=0; i<numSites; i++)
			// only name matters here
			sites.add(new BBP_Site("s"+i, null, Double.NaN, Double.NaN, Double.NaN));
		
		return sites;
	}

	public BBP_PartBSimZipLoader(ZipFile zip, Scenario[] scenarios, int numSites) {
		super(zip, buildSites(numSites));
		sites = buildSites(numSites);
		rdMap = new HashMap<>();
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
