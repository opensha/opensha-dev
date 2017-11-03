package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;

public class BBP_CatalogSimZipLoader extends BBP_SimZipLoader {
	
	private Table<Integer, BBP_Site, DiscretizedFunc> rd50Table;
	private Table<Integer, BBP_Site, DiscretizedFunc[]> rdTable;
	
	public BBP_CatalogSimZipLoader(File file, List<BBP_Site> sites) throws ZipException, IOException {
		this(new ZipFile(file), sites);
	}
	
	public BBP_CatalogSimZipLoader(ZipFile zip, List<BBP_Site> sites) throws ZipException, IOException {
		super(zip, sites);
		rd50Table = HashBasedTable.create();
		if (hasRotD100())
			rdTable = HashBasedTable.create();
	}
	
	private static String getDirName(int eventID) {
		return "event_"+eventID;
	}
	
	private static int getEventID(String dirName) {
		Preconditions.checkState(dirName.startsWith("event_"));
		return Integer.parseInt(dirName.substring("event_".length()));
	}
	
	public boolean contains(BBP_Site site, int eventID) {
		return contains(site, getDirName(eventID));
	}
	
	public DiscretizedFunc readRotD50(BBP_Site site, int eventID) throws IOException {
		if (rdTable != null)
			return readRotD(site, eventID)[0];
		if (rd50Table.contains(eventID, site))
			return rd50Table.get(eventID, site);
		synchronized (rd50Table) {
			// repeat here as could have been populated while waiting for the lock
			if (rd50Table.contains(eventID, site))
				return rd50Table.get(eventID, site);
			DiscretizedFunc spectra = readRotD50(site, getDirName(eventID));
			rd50Table.put(eventID, site, spectra);
			return spectra;
		}
	}
	
	public DiscretizedFunc[] readRotD(BBP_Site site, int eventID) throws IOException {
		if (rdTable.contains(eventID, site))
			return rdTable.get(eventID, site);
		synchronized (rdTable) {
			// repeat here as could have been populated while waiting for the lock
			if (rdTable.contains(eventID, site))
				return rdTable.get(eventID, site);
			DiscretizedFunc[] spectras = readRotD(site, getDirName(eventID));
			rdTable.put(eventID, site, spectras);
			return spectras;
		}
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
			if (!dirName.contains("rup_gen"))
				ids.add(getEventID(dirName));
		return ids;
	}
	
	private Table<Integer, BBP_Site, DiscretizedFunc[]> rd50RupGenTable;
	
	public synchronized DiscretizedFunc[] readRupGenRotD50(BBP_Site site, int eventID) {
		if (rd50RupGenTable == null)
			rd50RupGenTable = HashBasedTable.create();
		DiscretizedFunc[] ret = rd50RupGenTable.get(eventID, site);
		if (ret == null) {
			int cnt = 0;
			String dirName = getDirName(eventID);
			List<DiscretizedFunc> funcs = new ArrayList<>();
			while (true) {
				try {
					DiscretizedFunc rd50 = readRotD50(site, dirName+"/rup_gen_"+(cnt++));
					funcs.add(rd50);
				} catch (Exception e) {
					break;
				}
			}
			Preconditions.checkState(funcs.size() > 0, "No RG files for %s, event %s", site.getName(), eventID);
			ret = funcs.toArray(new DiscretizedFunc[funcs.size()]);
			rd50RupGenTable.put(eventID, site, ret);
		}
		return ret;
	}
	
	public static void main(String[] args) throws IOException {
		File file = new File("/home/kevin/bbp/parallel/2017_10_09-rundir2194_long-all-m6.5-skipYears5000-noHF/results.zip");
		List<BBP_Site> sites = BBP_Site.readFile(file.getParentFile());
		
		new BBP_CatalogSimZipLoader(file, sites);
	}

}
