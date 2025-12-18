package scratch.kevin.prvi25;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Stream;

import org.opensha.commons.data.CSVReader;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.util.io.archive.ArchiveInput;
import org.opensha.commons.util.io.archive.CopyAvoidantInMemorySeekableByteChannel;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import scratch.kevin.ZipMD5Compare;

public class GriddedRebuildValidation {

	public static void main(String[] args) throws IOException {
		// validate that nothing changed in gridded seismicity rebuild besides removing duplication

		File inputDir = new File("/project2/scec_608/kmilner/fss_inversions/"
				+ "2025_09_12-prvi25_crustal_branches-dmSample10x");

		File origSLT = new File(inputDir, "results_full_gridded.zip.orig");
		File newSLT = new File(inputDir, "results_full_gridded.zip");

		ArchiveInput origZip = ArchiveInput.getDefaultInput(origSLT);
		ArchiveInput newZip = ArchiveInput.getDefaultInput(newSLT);
		
		origZip = new RemovePointZerosWrapper(origZip);
		newZip = new RemovePointZerosWrapper(newZip);

		// build MD5 sums for everything
		int threads = FaultSysTools.defaultNumThreads();

		// create mappings
		Map<String, List<String>> newToOldMappings = new HashMap<>();
		System.out.println("Mapping old entries to new entries");
		int[] compMappings = {-1, 1, 2, -1, 3, -1, 4, 5, 6, 7, 8, -1};
		List<String> newEntries = new ArrayList<>();
		List<String[]> newEntryComponents = new ArrayList<>();
		for (String entry : newZip.getEntries()) {
			if (entry.endsWith("grid_sources.csv")) {
				newEntries.add(entry);
				newEntryComponents.add(entry.split("/"));
			}
		}
		for (String origEntry : origZip.getEntries()) {
			if (!origEntry.endsWith("grid_sources.csv"))
				continue;
			String[] origPathComponents = origEntry.split("/");
			Preconditions.checkState(origPathComponents.length == compMappings.length,
					"Orig path components length mismatch: %s != %s\n\t%s",
					compMappings.length, origPathComponents.length, origPathComponents);
			boolean mapped = false;
			for (int i=0; i<newEntries.size(); i++) {
				String newEntry = newEntries.get(i);
				String[] newPathComponents = newEntryComponents.get(i);
				boolean match = true;
				for (int j=0; j<compMappings.length; j++) {
					if (compMappings[j] >= 0 && !origPathComponents[j].equals(newPathComponents[compMappings[j]])) {
						match = false;
						break;
					}
				}
				if (match) {
					Preconditions.checkState(!mapped, "Multiply mapped: \n\tOrig:\t%s\n\tThis new:\t%s",
							origPathComponents, newPathComponents);
					List<String> mappings = newToOldMappings.get(newEntry);
					if (mappings == null) {
						mappings = new ArrayList<>();
						newToOldMappings.put(newEntry, mappings);
					}
					mappings.add(origEntry);
					mapped = true;
				}
			}
			Preconditions.checkState(mapped, "No mappings found for:\t%s", origEntry);
		}
		Preconditions.checkState(newToOldMappings.size() == newEntries.size(),
				"Only mapped %s/%s new entries", newToOldMappings.size(), newEntries.size());
		int numMappingsEach = newToOldMappings.values().iterator().next().size();
		for (String newEntry : newToOldMappings.keySet()) {
			int myNum = newToOldMappings.get(newEntry).size();
			Preconditions.checkState(numMappingsEach == myNum,
					"Mapping count mismatch for %s: %s != %s", newEntry, myNum, numMappingsEach);
		}
		System.out.println("Mapped "+numMappingsEach+" original entries for each of "+newEntries.size()+" new entries");

		ExecutorService exec = Executors.newFixedThreadPool(threads);
		System.out.println("Calculating MD5s for new SLT");
		Map<String, String> newMD5s = ZipMD5Compare.loadCalcMD5s(newZip, exec);
		System.out.println("\tFound "+numUniqueGridSoruces(newMD5s)+" unique new grid source MD5s");

		System.out.println("Calculating MD5s for orig SLT");
		Map<String, String> origMD5s = ZipMD5Compare.loadCalcMD5s(origZip, exec);
		System.out.println("\tFound "+numUniqueGridSoruces(origMD5s)+" unique orig grid source MD5s");
		exec.shutdown();

		double absTol = 1e-8;
		double relTol = 1e-5;
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		
		// find the gridded region
		Feature gridRegFeature = Feature.read(new InputStreamReader(new BufferedInputStream(
				newZip.getInputStream("solution_logic_tree/grid_region.geojson"))));
		GriddedRegion gridReg = GriddedRegion.fromFeature(gridRegFeature);

		System.out.println("Comparing MD5s");
		int numPerfectMatches = 0;
		int numMD5Mismatches = 0;
		for (String newEntry : newMD5s.keySet()) {
			if (!newEntry.endsWith("grid_sources.csv"))
				continue;
			String newMD5 = newMD5s.get(newEntry);
			System.out.println("Processing MD5s for:\t"+newMD5+"\t"+newEntry);

			List<String> mappings = newToOldMappings.get(newEntry);
			int numValidated = 0;
			HashSet<String> md5s = new HashSet<>(1);
			GridSourceList newList = null;
			for (String origEntry : mappings) {
				String origMD5 = origMD5s.get(origEntry);
				if (newMD5.equals(origMD5)) {
					numValidated++;
				} else {
					if (!md5s.contains(origMD5)) {
						System.out.println("\tMismatching MD5: "+origMD5+"\t"+origEntry+", will test GridSourceList for statistical equality");
						// first time we've encountered this one
						if (newList == null) {
							// need to load it
							EnumMap<TectonicRegionType, List<List<GriddedRupture>>> newRups =
									GridSourceList.loadGridSourcesCSV(new CSVReader(newZip.getInputStream(newEntry)), gridReg.getNodeList());
							newList = new GridSourceList.Precomputed(gridReg, newRups);
						}
						EnumMap<TectonicRegionType, List<List<GriddedRupture>>> origRups =
								GridSourceList.loadGridSourcesCSV(new CSVReader(origZip.getInputStream(origEntry)), gridReg.getNodeList());
						GridSourceList origList = new GridSourceList.Precomputed(gridReg, origRups);
						
						SolGriddedToleranceCompare.compare(origList, newList, absTol, relTol);
					}
				}
				md5s.add(origMD5);
			}
			System.out.println("\tHave "+numValidated+"/"+mappings.size()+" perfect MD5 matches");
			if (numValidated == mappings.size()) {
				numPerfectMatches++;
			} else {
				numMD5Mismatches++;
			}
//			Preconditions.checkState(numValidated == mappings.size(), "%s of %s didn't match for %s",
//					(mappings.size()-numValidated), mappings.size(), newEntry);
		}
		System.out.println(numPerfectMatches+"/"+newEntries.size()+" perfect matches ("+numMD5Mismatches+" mismatches)");
	}
	
	private static int numUniqueGridSoruces(Map<String, String> md5s) {
		HashSet<String> md5sSet = new HashSet<>();
		for (String entry : md5s.keySet()) {
			if (!entry.endsWith("grid_sources.csv"))
				continue;
			md5sSet.add(md5s.get(entry));
		}
		return md5sSet.size();
	}
	
	/**
	 * Often, the only differences will be rounding on if X.0 is displayed as X.0 or just X. This removes all trailing .0's
	 */
	private static class RemovePointZerosWrapper implements ArchiveInput {
		
		private ArchiveInput wrapped;
		private CopyAvoidantInMemorySeekableByteChannel byteChannel = new CopyAvoidantInMemorySeekableByteChannel();

		public RemovePointZerosWrapper(ArchiveInput wrapped) {
			this.wrapped = wrapped;
		}

		@Override
		public String getName() {
			return "Wrapped "+wrapped.getName();
		}

		@Override
		public void close() throws IOException {
			wrapped.close();
		}

		@Override
		public boolean hasEntry(String name) throws IOException {
			return wrapped.hasEntry(name);
		}

		@Override
		public InputStream getInputStream(String name) throws IOException {
			if (name.endsWith("grid_sources.csv")) {
				BufferedReader bRead = new BufferedReader(new InputStreamReader(wrapped.getInputStream(name)));
				byteChannel.truncate(0l);
				BufferedWriter bWrite = new BufferedWriter(new OutputStreamWriter(byteChannel.getOutputStream()));
				while (true) {
					String line = bRead.readLine();
					if (line == null)
						break;
					line = line.replace(".0,", ",");
					bWrite.write(line);
					bWrite.write('\n');
				}
				bWrite.flush();
				byteChannel.position(0l);
				return byteChannel.getInputStream();
			}
			return wrapped.getInputStream(name);
		}

		@Override
		public Stream<String> entryStream() throws IOException {
			return wrapped.entryStream();
		}
		
	}

}
