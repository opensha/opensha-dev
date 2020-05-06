package scratch.kevin.bbp;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.imr.param.IntensityMeasureParams.DurationTimeInterval;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

public class BBP_SimZipLoader {
	
	protected enum FileType {
		RotD50(".rd50"),
		RotD100(".rd100"),
		RotDPGV(".rdpgv", ".rdvel"),
		ARIAS_DURATION(".ard"),
		VEL_SEIS(".vel.bbp"),
		ACC_SEIS(".acc.bbp"),
		FAS(".fs.col");
		
		private String[] exts;
		private FileType(String... exts) {
			this.exts = exts;
		}
		
		public boolean matches(String name) {
			for (String ext : exts)
				if (name.endsWith(ext))
					return true;
			return false;
		}
	}
	
	public static class BBP_RupGenSimZipLoader extends BBP_SimZipLoader {
		
		private boolean individualSites;

		public BBP_RupGenSimZipLoader(File file, List<BBP_Site> sites) throws ZipException, IOException {
			super(file, sites);
			individualSites = hasIndividualSites(sites);
		}
		
		public BBP_RupGenSimZipLoader(ZipFile zip, List<BBP_Site> sites) throws ZipException, IOException {
			super(zip, sites);
			individualSites = hasIndividualSites(sites);
		}
		
		private boolean hasIndividualSites(List<BBP_Site> sites) {
			String dirName = getEntriesTable().columnKeySet().iterator().next();
			for (BBP_Site site : sites)
				if (dirName.contains(site.getName()))
					return true;
			return false;
		}
		
		public int getNumSims() {
			BBP_Site firstSite = getEntriesTable().rowKeySet().iterator().next();
			return getEntriesTable().row(firstSite).keySet().size();
		}
		
		private String getDirName(int index, BBP_Site site) {
			if (individualSites)
				return "run_"+index+"_"+site.getName();
			return "run_"+index;
		}
		
		public DiscretizedFunc readRotD50(BBP_Site site, int index) throws IOException {
			return readRotD50(site, getDirName(index, site));
		}
		
		public DiscretizedFunc readRotD100(BBP_Site site, int index) throws IOException {
			return readRotD100(site, getDirName(index, site));
		}
		
		public DiscretizedFunc[] readRotD(BBP_Site site, int index) throws IOException {
			return readRotD(site, getDirName(index, site));
		}
		
		public DiscretizedFunc readFAS(BBP_Site site, int index) throws IOException {
			return readFAS(site, getDirName(index, site));
		}
		
		public DiscretizedFunc[] readAccelSeis(BBP_Site site, int index) throws IOException {
			return readAccelSeis(site, getDirName(index, site));
		}
		
		public DiscretizedFunc[] readVelSeis(BBP_Site site, int index) throws IOException {
			return readVelSeis(site, getDirName(index, site));
		}
		
	}
	
	public static class BBP_ShakeMapSimZipLoader extends BBP_SimZipLoader {
		
		private List<BBP_Site> sites;

		public BBP_ShakeMapSimZipLoader(File file, List<BBP_Site> sites) throws ZipException, IOException {
			super(file, sites);
			this.sites = sites;
		}
		
		public BBP_ShakeMapSimZipLoader(ZipFile zip, List<BBP_Site> sites) throws ZipException, IOException {
			super(zip, sites);
			this.sites = sites;
		}
		
		public int getNumSims() {
			return getEntriesTable().rowKeySet().size();
		}
		
		private String getDirName(int index) {
			return "site_"+index;
		}
		
		public DiscretizedFunc readRotD50(int index) throws IOException {
			return readRotD50(sites.get(index), getDirName(index));
		}
		
		public DiscretizedFunc readRotD100(int index) throws IOException {
			return readRotD100(sites.get(index), getDirName(index));
		}
		
		public DiscretizedFunc[] readRotD(int index) throws IOException {
			return readRotD(sites.get(index), getDirName(index));
		}
		
		public double[] readPGV(int index) throws IOException {
			return readPGV(sites.get(index), getDirName(index));
		}
		
		public DiscretizedFunc readFAS(int index) throws IOException {
			return readFAS(sites.get(index), getDirName(index));
		}
		
		public DiscretizedFunc[] readAccelSeis(int index) throws IOException {
			return readAccelSeis(sites.get(index), getDirName(index));
		}
		
		public DiscretizedFunc[] readVelSeis(int index) throws IOException {
			return readVelSeis(sites.get(index), getDirName(index));
		}

		@Override
		protected BBP_Site forEntry(String entryName, List<BBP_Site> sites) {
			if (entryName.startsWith("site_")) {
				if (entryName.contains("/"))
					entryName = entryName.substring(0, entryName.indexOf("/"));
				return sites.get(Integer.parseInt(entryName.substring(entryName.indexOf("_")+1)));
			}
			return null;
		}
		
	}
	
	private Table<BBP_Site, String, Map<FileType, ZipEntry>> entriesTable;
	private HashSet<FileType> allTypes;
	private ZipFile zip;
	protected List<BBP_Site> sites;
	
	public BBP_SimZipLoader(File file, List<BBP_Site> sites) throws ZipException, IOException {
		this(new ZipFile(file), sites);
	}
	
	public BBP_SimZipLoader(ZipFile zip, List<BBP_Site> sites) {
		this.zip = zip;
		this.sites = Collections.unmodifiableList(sites);
		entriesTable = HashBasedTable.create();
		allTypes = new HashSet<>();
		
		System.out.println("Mapping entries...");
		int numFiles = 0;
		int numSims = 0;
		Enumeration<? extends ZipEntry> entries = zip.entries();
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			String name = entry.getName();
			BBP_Site site = forEntry(name, sites);
			if (site != null) {
				for (FileType type : FileType.values()) {
					if (type.matches(name)) {
						Preconditions.checkState(name.contains("/"));
						String dirName = name.substring(0, name.lastIndexOf("/"));
						if (dirName.endsWith("/FAS"))
							dirName = dirName.substring(0, dirName.lastIndexOf("/"));
						Map<FileType, ZipEntry> fileMap = entriesTable.get(site, dirName);
						if (fileMap == null) {
							fileMap = new HashMap<>();
							entriesTable.put(site, dirName, fileMap);
							numSims++;
						}
						fileMap.put(type, entry);
						allTypes.add(type);
						numFiles++;
						break;
					}
				}
			}
			
		}
		System.out.println("Loaded "+numFiles+" files from "+numSims+" sims for "
				+entriesTable.rowKeySet().size()+" sites and "+entriesTable.columnKeySet().size()+" dirs");
	}
	
	protected boolean contains(BBP_Site site, String entryName) {
		return entriesTable.contains(site, entryName);
	}
	
	protected BBP_Site forEntry(String entryName, List<BBP_Site> sites) {
		for (BBP_Site site : sites)
			if (entryName.contains("."+site.getName()+"."))
				return site;
		return null;
	}
	
	protected Table<BBP_Site, String, Map<FileType, ZipEntry>> getEntriesTable() {
		return entriesTable;
	}
	
	private ZipEntry locate(BBP_Site site, String dirName, FileType type) {
		Map<FileType, ZipEntry> fileMap = entriesTable.get(site, dirName);
//		if (fileMap == null) {
//			System.out.println("Debug!");
//			Map<BBP_Site, Map<FileType, ZipEntry>> col = entriesTable.column(dirName);
//			System.out.println("Col null? "+(col == null));
//			System.out.println("Col empty? "+col.isEmpty());
//			if (col != null) {
//				for (BBP_Site key : col.keySet())
//					System.out.println("\t"+key.getName());
//			}
//		}
		Preconditions.checkNotNull(fileMap, "No files for dir %s, site %s", dirName, site.getName());
		ZipEntry entry = fileMap.get(type);
		Preconditions.checkNotNull(entry, "No file of type %s for dir %s, site %s", type.name(), dirName, site.getName());
		return entry;
	}
	
	private List<String> loadFileLines(ZipEntry entry) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
		List<String> lines = new ArrayList<>();
		
		String line;
		while ((line = br.readLine()) != null)
			lines.add(line);
		
		return lines;
	}
	
	public DiscretizedFunc readRotD50(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry;
		try {
			entry = locate(site, dirName, FileType.RotD50);
		} catch (NullPointerException e) {
			// RotD50 also stored in RotD100 files
			entry = locate(site, dirName, FileType.RotD100);
		}
		try {
			return SpectraPlotter.loadRotD50(loadFileLines(entry));
		} catch (RuntimeException e) {
			System.err.println("Exception loading from "+entry.getName());
			throw e;
		}
	}
	
	public DiscretizedFunc readRotD100(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.RotD100);
		try {
			return SpectraPlotter.loadRotD100(loadFileLines(entry));
		} catch (RuntimeException e) {
			System.err.println("Exception loading from "+entry.getName());
			throw e;
		}
	}
	
	public DiscretizedFunc[] readRotD(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.RotD100);
		try {
			return SpectraPlotter.loadRotD(loadFileLines(entry));
		} catch (RuntimeException e) {
			System.err.println("Exception loading from "+entry.getName());
			throw e;
		}
	}
	
	public double[] readPGV(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.RotDPGV);
		try {
			DiscretizedFunc[] spectra = SpectraPlotter.loadRotD(loadFileLines(entry));
			double[] ret = { spectra[0].getY(0), spectra[1].getY(0) };
			return ret;
		} catch (RuntimeException e) {
			System.err.println("Exception loading from "+entry.getName());
			throw e;
		}
	}
	
	/**
	 * 
	 * @param site
	 * @param dirName
	 * @return map from duration interval to array containing [N, E, Z, GEOM] durations
	 * @throws IOException
	 */
	public Map<DurationTimeInterval, double[]> readDurations(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.ARIAS_DURATION);
		try {
			Map<DurationTimeInterval, double[]> ret = new HashMap<>();
			ret.put(DurationTimeInterval.INTERVAL_5_75, new double[4]);
			ret.put(DurationTimeInterval.INTERVAL_5_95, new double[4]);
			ret.put(DurationTimeInterval.INTERVAL_20_80, new double[4]);
			// format:
			// Component  Peak Arias (cm/s)  T5-75 (s)  T5-95 (s)  T20-80 (s)
			// 
			for (String line : loadFileLines(entry)) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#"))
					continue;
				line = line.replaceAll("\t", " ");
				while (line.contains("  "))
					line = line.replaceAll("  ", " ");
				String[] split = line.split(" ");
				Preconditions.checkState(split.length == 5, "Bad arias line: %s", line);
				int arrayIndex;
				if (split[0].equals("N"))
					arrayIndex = 0;
				else if (split[0].equals("E"))
					arrayIndex = 1;
				else if (split[0].equals("Z"))
					arrayIndex = 2;
				else if (split[0].equals("GEOM"))
					arrayIndex = 3;
				else
					throw new IllegalStateException("Unknown compoent for line: "+line);
				double d5_75 = Double.parseDouble(split[2]);
				double d5_95 = Double.parseDouble(split[3]);
				double d20_80 = Double.parseDouble(split[4]);
				ret.get(DurationTimeInterval.INTERVAL_5_75)[arrayIndex] = d5_75;
				ret.get(DurationTimeInterval.INTERVAL_5_95)[arrayIndex] = d5_95;
				ret.get(DurationTimeInterval.INTERVAL_20_80)[arrayIndex] = d20_80;
			}
			return ret;
		} catch (RuntimeException e) {
			System.err.println("Exception loading from "+entry.getName());
			throw e;
		}
	}
	
	public boolean hasRotD50() {
		return allTypes.contains(FileType.RotD50);
	}
	
	public boolean hasRotD100() {
		return allTypes.contains(FileType.RotD100);
	}
	
	public boolean hasPGV() {
		return allTypes.contains(FileType.RotDPGV);
	}
	
	public boolean hasDurations() {
		return allTypes.contains(FileType.ARIAS_DURATION);
	}
	
	public boolean hasFAS() {
		return allTypes.contains(FileType.FAS);
	}
	
	public DiscretizedFunc readFAS(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.FAS);
		return SpectraPlotter.loadRotD50(loadFileLines(entry));
	}
	
	public boolean hasAccelSeis() {
		return allTypes.contains(FileType.ACC_SEIS);
	}
	
	public DiscretizedFunc[] readAccelSeis(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.ACC_SEIS);
		return SeismogramPlotter.loadBBP_Seis(loadFileLines(entry));
	}
	
	public boolean hasVelSeis() {
		return allTypes.contains(FileType.VEL_SEIS);
	}
	
	public DiscretizedFunc[] readVelSeis(BBP_Site site, String dirName) throws IOException {
		ZipEntry entry = locate(site, dirName, FileType.VEL_SEIS);
		return SeismogramPlotter.loadBBP_Seis(loadFileLines(entry));
	}
	
	public Collection<String> getDirNames(BBP_Site site) {
		return entriesTable.row(site).keySet();
	}

}
