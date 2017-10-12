package scratch.kevin.simulators.ruptures;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.DiscretizedFunc;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.SeismogramPlotter;
import scratch.kevin.bbp.SpectraPlotter;

public class BBP_CatalogSimZipLoader {
	
	private enum FileType {
		RotD50(".rd50"),
		VEL_SEIS(".vel.bbp"),
		ACC_SEIS(".acc.bbp"),
		FAS(".fs.col");
		
		private String ext;
		private FileType(String ext) {
			this.ext = ext;
		}
		
		public boolean matches(String name) {
			return name.endsWith(ext);
		}
	}
	
	private Table<BBP_Site, Integer, Map<FileType, ZipEntry>> entriesTable;
	private ZipFile zip;
	
	public BBP_CatalogSimZipLoader(File file, List<BBP_Site> sites) throws ZipException, IOException {
		this(new ZipFile(file), sites);
	}
	
	public BBP_CatalogSimZipLoader(ZipFile zip, List<BBP_Site> sites) {
		this.zip = zip;
		entriesTable = HashBasedTable.create();
		
		System.out.println("Mapping entries...");
		int numFiles = 0;
		int numSims = 0;
		Enumeration<? extends ZipEntry> entries = zip.entries();
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			String name = entry.getName();
			for (BBP_Site site : sites) {
				if (name.contains(site.getName())) {
					for (FileType type : FileType.values()) {
						if (type.matches(name)) {
							Preconditions.checkState(name.contains("/"));
							String eventPortion = name.substring(0, name.indexOf("/"));
							Preconditions.checkState(eventPortion.startsWith("event_"));
							int eventID = Integer.parseInt(eventPortion.substring("event_".length()));
							Map<FileType, ZipEntry> fileMap = entriesTable.get(site, eventID);
							if (fileMap == null) {
								fileMap = new HashMap<>();
								entriesTable.put(site, eventID, fileMap);
								numSims++;
							}
							fileMap.put(type, entry);
							numFiles++;
							break;
						}
					}
					break;
				}
			}
		}
		System.out.println("Loaded "+numFiles+" files from "+numSims+" sims for "
				+entriesTable.rowKeySet().size()+" sites and "+entriesTable.columnKeySet().size()+" events");
	}
	
	private ZipEntry locate(BBP_Site site, int eventID, FileType type) {
		Map<FileType, ZipEntry> fileMap = entriesTable.get(site, eventID);
		Preconditions.checkNotNull(fileMap, "No files for event %s, site %s", eventID, site.getName());
		ZipEntry entry = fileMap.get(type);
		Preconditions.checkNotNull("No file of type %s for event %s, site %s", type.name(), eventID, site.getName());
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
	
	public DiscretizedFunc readRotD50(BBP_Site site, int eventID) throws IOException {
		ZipEntry entry = locate(site, eventID, FileType.RotD50);
		return SpectraPlotter.loadRotD50(loadFileLines(entry));
	}
	
	public DiscretizedFunc readFAS(BBP_Site site, int eventID) throws IOException {
		ZipEntry entry = locate(site, eventID, FileType.FAS);
		return SpectraPlotter.loadRotD50(loadFileLines(entry));
	}
	
	public DiscretizedFunc[] readAccelSeis(BBP_Site site, int eventID) throws IOException {
		ZipEntry entry = locate(site, eventID, FileType.ACC_SEIS);
		return SeismogramPlotter.loadBBP_Seis(loadFileLines(entry));
	}
	
	public DiscretizedFunc[] readVelSeis(BBP_Site site, int eventID) throws IOException {
		ZipEntry entry = locate(site, eventID, FileType.VEL_SEIS);
		return SeismogramPlotter.loadBBP_Seis(loadFileLines(entry));
	}
	
	public Set<Integer> getEventIDs(BBP_Site site) {
		return entriesTable.row(site).keySet();
	}
	
	public static void main(String[] args) throws IOException {
		File file = new File("/home/kevin/bbp/parallel/2017_10_09-rundir2194_long-all-m6.5-skipYears5000-noHF/results.zip");
		List<BBP_Site> sites = BBP_Site.readFile(file.getParentFile());
		
		new BBP_CatalogSimZipLoader(file, sites);
	}

}
