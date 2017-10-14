package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.DiscretizedFunc;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;

public class BBP_CatalogSimZipLoader extends BBP_SimZipLoader {
	
	public BBP_CatalogSimZipLoader(File file, List<BBP_Site> sites) throws ZipException, IOException {
		super(file, sites);
	}
	
	public BBP_CatalogSimZipLoader(ZipFile zip, List<BBP_Site> sites) throws ZipException, IOException {
		super(zip, sites);
	}
	
	private static String getDirName(int eventID) {
		return "event_"+eventID;
	}
	
	private static int getEventID(String dirName) {
		Preconditions.checkState(dirName.startsWith("event_"));
		return Integer.parseInt(dirName.substring("event_".length()));
	}
	
	public DiscretizedFunc readRotD50(BBP_Site site, int eventID) throws IOException {
		return readRotD50(site, getDirName(eventID));
	}
	
	public DiscretizedFunc readFAS(BBP_Site site, int eventID) throws IOException {
		return readFAS(site, getDirName(eventID));
	}
	
	public DiscretizedFunc[] readAccelSeis(BBP_Site site, int eventID) throws IOException {
		return readAccelSeis(site, getDirName(eventID));
	}
	
	public DiscretizedFunc[] readVelSeis(BBP_Site site, int eventID) throws IOException {
		return readVelSeis(site, getDirName(eventID));
	}
	
	public Collection<Integer> getEventIDs(BBP_Site site) {
		Collection<String> dirNames = super.getDirNames(site);
		List<Integer> ids = new ArrayList<>(dirNames.size());
		for (String dirName : dirNames)
			ids.add(getEventID(dirName));
		return ids;
	}
	
	public static void main(String[] args) throws IOException {
		File file = new File("/home/kevin/bbp/parallel/2017_10_09-rundir2194_long-all-m6.5-skipYears5000-noHF/results.zip");
		List<BBP_Site> sites = BBP_Site.readFile(file.getParentFile());
		
		new BBP_CatalogSimZipLoader(file, sites);
	}

}
