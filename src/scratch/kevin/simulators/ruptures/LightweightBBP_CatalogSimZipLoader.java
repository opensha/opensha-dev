package scratch.kevin.simulators.ruptures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.SimulationRotDProvider;

public class LightweightBBP_CatalogSimZipLoader extends BBP_SimZipLoader implements SimulationRotDProvider<Integer> {
	
	private Map<Site, BBP_Site> gmpeToBBP;
	
	private double durationYears;
	
	private Map<BBP_Site, List<Integer>> eventIDsMap;
	
	private Table<Integer, BBP_Site, DiscretizedFunc> rd50Table;
	private Table<Integer, Site, DiscretizedFunc> rdRatioTable;
	private Table<Integer, BBP_Site, DiscretizedFunc[]> rdTable;
	
	public LightweightBBP_CatalogSimZipLoader(ZipFile zip, List<BBP_Site> sites, BiMap<BBP_Site, Site> gmpeSites, double durationYears) {
		super(zip, sites);
		
		this.durationYears = durationYears;
		eventIDsMap = new HashMap<>();
		
		for (Cell<BBP_Site, String, Map<FileType, ZipEntry>> cell : super.getEntriesTable().cellSet()) {
			String dirName = cell.getColumnKey();
			if (dirName.startsWith("event_")) {
				int id = Integer.parseInt(dirName.substring(6));
				BBP_Site site = cell.getRowKey();
				List<Integer> eventIDs = eventIDsMap.get(site);
				if (eventIDs == null) {
					eventIDs = new ArrayList<>();
					eventIDsMap.put(site, eventIDs);
				}
				eventIDs.add(id);
			}
		}

		System.out.println("Detected "+eventIDsMap.size()+" RSQSim events/site mappings in zip file");
		
		gmpeToBBP = gmpeSites.inverse();
		rd50Table = HashBasedTable.create();
		if (hasRotD100()) {
			rdTable = HashBasedTable.create();
			rdRatioTable = HashBasedTable.create();
		}
	}

	@Override
	public String getName() {
		return "RSQSim/BBP";
	}

	@Override
	public DiscretizedFunc getRotD50(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		if (rdTable != null)
			return getRotD(site, eventID, index)[0];
		if (rd50Table.contains(eventID, site))
			return rd50Table.get(eventID, site);
		synchronized (rd50Table) {
			// repeat here as could have been populated while waiting for the lock
			if (rd50Table.contains(eventID, site))
				return rd50Table.get(eventID, site);
			BBP_Site bbpSite = gmpeToBBP.get(site);
			DiscretizedFunc spectra = readRotD50(bbpSite, getDirName(eventID));
			rd50Table.put(eventID, bbpSite, spectra);
			return spectra;
		}
	}

	@Override
	public DiscretizedFunc getRotD100(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		return getRotD(site, eventID, index)[1];
	}

	@Override
	public DiscretizedFunc[] getRotD(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		if (rdTable.contains(eventID, site))
			return rdTable.get(eventID, site);
		synchronized (rdTable) {
			// repeat here as could have been populated while waiting for the lock
			if (rdTable.contains(eventID, site))
				return rdTable.get(eventID, site);
			BBP_Site bbpSite = gmpeToBBP.get(site);
			DiscretizedFunc[] spectras = readRotD(bbpSite, getDirName(eventID));
			rdTable.put(eventID, bbpSite, spectras);
			return spectras;
		}
	}

	@Override
	public DiscretizedFunc getRotDRatio(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		if (rdRatioTable.contains(eventID, site))
			return rdRatioTable.get(eventID, site);
		synchronized (rdTable) {
			// repeat here as could have been populated while waiting for the lock
			if (rdRatioTable.contains(eventID, site))
				return rdRatioTable.get(eventID, site);
			DiscretizedFunc[] spectras = getRotD(site, eventID, index);
			DiscretizedFunc ratio = SimulationRotDProvider.calcRotDRatio(spectras);
			rdRatioTable.put(eventID, site, ratio);
			return ratio;
		}
	}
	
	private static String getDirName(int eventID) {
		return "event_"+eventID;
	}

	@Override
	public int getNumSimulations(Site site, Integer eventID) {
		return 1;
	}

	@Override
	public Collection<Integer> getRupturesForSite(Site site) {
		return eventIDsMap.get(gmpeToBBP.get(site));
	}

	@Override
	public double getAnnualRate(Integer rupture) {
		return 1d/durationYears;
	}

	@Override
	public double getMinimumCurvePlotRate(Site site) {
		return getAnnualRate(null);
	}

	@Override
	public double getMagnitude(Integer rupture) {
		return Double.NaN;
	}

	@Override
	public Location getHypocenter(Integer rupture, int index) {
		return null;
	}

}
