package scratch.kevin.simulators.ruptures.azimuthal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.zip.ZipException;

import org.opensha.commons.data.function.DiscretizedFunc;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;

public class AzimuthalZipLoader extends BBP_SimZipLoader {

	private String scenarioPrefix;
	private int siteBinSize;
	private int numSites;

	public AzimuthalZipLoader(File file, int numSites, String scenarioPrefix)
			throws ZipException, IOException {
		super(file, buildSites(numSites));
		this.numSites = numSites;
		
		// detect site binning
		Collection<String> dirNames = getDirNames(sites.get(0));
		siteBinSize = -1;
		for (String dirName : dirNames) {
			dirName = dirName.substring(dirName.indexOf("_event"));
			if (dirName.contains("_s")) {
				dirName = dirName.substring(dirName.indexOf("_s")+2);
				String[] split = dirName.split("_");
				Preconditions.checkState(split.length == 2);
				int startIndex = Integer.parseInt(split[0]);
				int endIndex = Integer.parseInt(split[1]);
				int span = endIndex - startIndex;
				siteBinSize = Integer.max(siteBinSize, span);
			}
		}
		
		if (siteBinSize > 0)
			System.out.println("Sites binned in batches of "+siteBinSize);
		
		this.scenarioPrefix = scenarioPrefix;
	}
	
	public DiscretizedFunc readRotD50(int siteIndex, int eventID) throws IOException {
		String dirName = scenarioPrefix+"_event_"+eventID;
		if (siteBinSize > 0) {
			int binNum = siteIndex / siteBinSize;
			int binStart = binNum * siteBinSize;
			int binEnd = Integer.min(numSites, binStart + siteBinSize);
			dirName += "_s"+binStart+"_"+binEnd;
		}
		return readRotD50(sites.get(siteIndex), dirName);
	}
	
	private static List<BBP_Site> buildSites(int numRuptures) {
		List<BBP_Site> sites = new ArrayList<>();
		for (int i=0; i<numRuptures; i++)
			sites.add(new BBP_Site("s"+i, null, Double.NaN, Double.NaN, Double.NaN));
		return sites;
	}
}
