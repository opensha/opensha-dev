package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.SimulationRotDProvider;

public class BBP_CatalogSimZipLoader extends BBP_SimZipLoader implements SimulationRotDProvider<RSQSimEvent> {
	
	private BiMap<BBP_Site, Site> gmpeSites;
	private Map<Site, BBP_Site> gmpeToBBP;
	private Map<Integer, RSQSimEvent> eventsMap;
	private double durationYears;
	
	private Table<Integer, BBP_Site, DiscretizedFunc> rd50Table;
	private Table<Integer, Site, DiscretizedFunc> rdRatioTable;
	private Table<Integer, BBP_Site, DiscretizedFunc[]> rdTable;
	
	public BBP_CatalogSimZipLoader(File file, List<BBP_Site> sites, BiMap<BBP_Site, Site> gmpeSites,
			Map<Integer, RSQSimEvent> eventsMap) throws ZipException, IOException {
		this(new ZipFile(file), sites, gmpeSites, eventsMap);
	}
	
	public BBP_CatalogSimZipLoader(ZipFile zip, List<BBP_Site> sites, BiMap<BBP_Site, Site> gmpeSites,
			Map<Integer, RSQSimEvent> eventsMap) throws ZipException, IOException {
		super(zip, sites);
		this.gmpeSites = gmpeSites;
		this.eventsMap = eventsMap;
		double minTime = Double.POSITIVE_INFINITY;
		double maxTime = Double.NEGATIVE_INFINITY;
		for (RSQSimEvent event : eventsMap.values()) {
			double t = event.getTimeInYears();
			minTime = Math.min(minTime, t);
			maxTime = Math.max(maxTime, t);
		}
		durationYears = maxTime - minTime;
		gmpeToBBP = gmpeSites.inverse();
		rd50Table = HashBasedTable.create();
		if (hasRotD100()) {
			rdTable = HashBasedTable.create();
			rdRatioTable = HashBasedTable.create();
		}
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
//		File file = new File("/home/kevin/bbp/parallel/2017_10_09-rundir2194_long-all-m6.5-skipYears5000-noHF/results.zip");
//		List<BBP_Site> sites = BBP_Site.readFile(file.getParentFile());
//		
//		new BBP_CatalogSimZipLoader(file, sites);
	}

	@Override
	public DiscretizedFunc getRotD50(Site site, RSQSimEvent rupture, int index) throws IOException {
		Preconditions.checkState(index == 0);
		return readRotD50(gmpeToBBP.get(site), rupture.getID());
	}

	@Override
	public DiscretizedFunc getRotD100(Site site, RSQSimEvent rupture, int index) throws IOException {
		Preconditions.checkState(index == 0);
		return getRotD(site, rupture, index)[1];
	}

	@Override
	public DiscretizedFunc[] getRotD(Site site, RSQSimEvent rupture, int index) throws IOException {
		Preconditions.checkState(index == 0);
		return readRotD(gmpeToBBP.get(site), rupture.getID());
	}

	@Override
	public DiscretizedFunc getRotDRatio(Site site, RSQSimEvent rupture, int index) throws IOException {
		Preconditions.checkState(index == 0);
		if (rdRatioTable.contains(rupture.getID(), site))
			return rdRatioTable.get(rupture.getID(), site);
		synchronized (rdTable) {
			// repeat here as could have been populated while waiting for the lock
			if (rdRatioTable.contains(rupture.getID(), site))
				return rdRatioTable.get(rupture.getID(), site);
			DiscretizedFunc[] spectras = getRotD(site, rupture, index);
			DiscretizedFunc ratio = SimulationRotDProvider.calcRotDRatio(spectras);
			rdRatioTable.put(rupture.getID(), site, ratio);
			return ratio;
		}
	}

	@Override
	public Collection<RSQSimEvent> getRupturesForSite(Site site) {
		List<RSQSimEvent> events = new ArrayList<>();
		for (Integer id : getEventIDs(gmpeToBBP.get(site)))
			events.add(eventsMap.get(id));
		return events;
	}

	@Override
	public double getAnnualRate(RSQSimEvent rupture) {
		return 1d/durationYears;
	}

	@Override
	public String getName() {
		return "RSQSim/BBP";
	}

	@Override
	public double getMagnitude(RSQSimEvent rupture) {
		return rupture.getMagnitude();
	}

	@Override
	public int getNumSimulations(Site site, RSQSimEvent rupture) {
		if (contains(gmpeToBBP.get(site), rupture.getID()))
			return 1;
		return 0;
	}

	@Override
	public double getMinimumCurvePlotRate() {
		return getAnnualRate(null);
	}

	@Override
	public Location getHypocenter(RSQSimEvent rupture, int index) {
		return RSQSimUtils.getHypocenter(rupture);
	}

}
