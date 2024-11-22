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
import org.opensha.sha.imr.param.IntensityMeasureParams.DurationTimeInterval;

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
	
	private Table<Integer, Site, DiscretizedFunc> rd50Table;
	private Table<Integer, Site, DiscretizedFunc> rdRatioTable;
	private Table<Integer, Site, DiscretizedFunc[]> rdTable;
	private Table<Integer, Site, double[]> pgvTable;
	private Table<Integer, Site, Map<DurationTimeInterval, double[]>> durationTable;
	
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
		
		gmpeToBBP = gmpeSites == null ? null : gmpeSites.inverse();
		rd50Table = HashBasedTable.create();
		if (hasRotD100()) {
			rdTable = HashBasedTable.create();
			rdRatioTable = HashBasedTable.create();
		}
		if (hasPGV())
			pgvTable = HashBasedTable.create();
		if (hasDurations())
			durationTable = HashBasedTable.create();
	}
	
	public void preloadAllRotD(Site site) throws IOException {
		BBP_Site bbpSite = gmpeToBBP.get(site);
		Collection<Integer> rups = getRupturesForSite(site);
		if (rups == null || rups.isEmpty())
			return;
		if (rdTable != null) {
			Map<Integer, DiscretizedFunc[]> spectraMap = new HashMap<>(rups.size());
			for (int eventID : getRupturesForSite(site)) {
				DiscretizedFunc[] spectra = readRotD(bbpSite, getDirName(eventID));
				spectraMap.put(eventID, spectra);
			}
			synchronized (this) {
				rdTable.column(site).putAll(spectraMap);
			}
		} else {
			Map<Integer, DiscretizedFunc> spectraMap = new HashMap<>(rups.size());
			for (int eventID : getRupturesForSite(site)) {
				DiscretizedFunc spectra = readRotD50(bbpSite, getDirName(eventID));
				spectraMap.put(eventID, spectra);
			}
			synchronized (this) {
				rd50Table.column(site).putAll(spectraMap);
			}
		}
	}
	
	@Override
	public String getName() {
		return "RSQSim-BBP";
	}

	@Override
	public DiscretizedFunc getRotD50(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		if (rdTable != null)
			return getRotD(site, eventID, index)[0];
		DiscretizedFunc spectra = rd50Table.get(eventID, site);
		if (spectra != null)
			return spectra;
		BBP_Site bbpSite = gmpeToBBP.get(site);
		spectra = readRotD50(bbpSite, getDirName(eventID));
		synchronized (rd50Table) {
			rd50Table.put(eventID, site, spectra);
		}
		return spectra;
	}

	@Override
	public DiscretizedFunc getRotD100(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		return getRotD(site, eventID, index)[1];
	}

	@Override
	public DiscretizedFunc[] getRotD(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		DiscretizedFunc[] spectras = rdTable.get(eventID, site);
		if (spectras != null)
			return spectras;
		BBP_Site bbpSite = gmpeToBBP.get(site);
		spectras = readRotD(bbpSite, getDirName(eventID));
		synchronized (rdTable) {
			rdTable.put(eventID, site, spectras);
		}
		return spectras;
	}

	@Override
	public DiscretizedFunc getRotDRatio(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		DiscretizedFunc ratio = rdRatioTable.get(eventID, site);
		if (ratio != null)
			return ratio;
		DiscretizedFunc[] spectras = getRotD(site, eventID, index);
		ratio = SimulationRotDProvider.calcRotDRatio(spectras);
		synchronized (rdTable) {
			rdRatioTable.put(eventID, site, ratio);
		}
		return ratio;
	}

	@Override
	public double getPGV(Site site, Integer eventID, int index) throws IOException {
		Preconditions.checkState(index == 0);
		double[] vals = pgvTable.get(eventID, site);
		if (vals != null)
			return vals[0];
		BBP_Site bbpSite = gmpeToBBP.get(site);
		vals = readPGV(bbpSite, getDirName(eventID));
		
		synchronized (pgvTable) {
			pgvTable.put(eventID, site, vals);
		}
		return vals[0];
	}

	@Override
	public double getPGA(Site site, Integer eventID, int index) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public double getDuration(Site site, Integer eventID, DurationTimeInterval interval, int index)
			throws IOException {
		Preconditions.checkState(index == 0);
		Map<DurationTimeInterval, double[]> vals = durationTable.get(eventID, site);
		if (vals != null)
			return vals.get(interval)[3];
		BBP_Site bbpSite = gmpeToBBP.get(site);
		vals = readDurations(bbpSite, getDirName(eventID));
		synchronized (durationTable) {
			durationTable.put(eventID, site, vals);
		}
		return vals.get(interval)[3];
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

	@Override
	public boolean hasPGA() {
		return false;
	}

	@Override
	public double getRake(Integer rupture) {
		return Double.NaN;
	}

}
