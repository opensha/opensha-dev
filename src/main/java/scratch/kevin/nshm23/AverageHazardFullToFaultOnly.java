package scratch.kevin.nshm23;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;

import com.google.common.base.Preconditions;
import com.google.common.io.ByteStreams;

public class AverageHazardFullToFaultOnly {

	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("USAGE: <results_hazard_in.zip> <logic_tree_subset.json> <results_hazard_out.zip>");
			System.exit(1);
		}
		
		String[] mapPrefixes = {
				"map_pga_TWO_IN_50.txt",
				"map_pga_TEN_IN_50.txt",
				"map_1.0s_TWO_IN_50.txt",
				"map_1.0s_TEN_IN_50.txt",
		};
		
		File resultsIn = new File(args[0]);
		Preconditions.checkState(resultsIn.exists());
		
		File treeSubsetFile = new File(args[1]);
		Preconditions.checkState(treeSubsetFile.exists());
		LogicTree<?> subsetTree = LogicTree.read(treeSubsetFile);
		
		File resultsOut = new File(args[2]);
		Preconditions.checkState(!resultsOut.exists(),
				"Output file already exists, delete previous first: %s", resultsOut.getAbsolutePath());
		
		ZipFile	zip = new ZipFile(resultsIn);
		
		ZipEntry fullTreeEntry = zip.getEntry("logic_tree.json");
		BufferedReader treeRead = new BufferedReader(new InputStreamReader(zip.getInputStream(fullTreeEntry)));
		LogicTree<?> fullTree = LogicTree.read(treeRead);
		
		System.out.println("Downsampling from "+fullTree.size()+" to "+subsetTree.size());
		
		ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(resultsOut)));
		
		ZipEntry gridRegEntry = zip.getEntry("gridded_region.geojson");
		Preconditions.checkNotNull(gridRegEntry);
		BufferedInputStream regIS = new BufferedInputStream(zip.getInputStream(gridRegEntry));
		Feature gridRegFeature = Feature.read(new InputStreamReader(regIS));
		GriddedRegion gridReg = GriddedRegion.fromFeature(gridRegFeature);
		
		CompletableFuture<Void> writeFuture = null;
		
		// assume (and verify that) the full tree is in the same order as the subset tree, just deeper
		int fullBranchIndex = 0;
		for (int subsetIndex=0; subsetIndex<subsetTree.size(); subsetIndex++) {
			LogicTreeBranch<?> subsetBranch = subsetTree.getBranch(subsetIndex);
			
			System.out.println("Subset branch "+subsetIndex+"/"+subsetTree.size()+": "+subsetBranch);
			
			List<LogicTreeBranch<?>> fullMatches = new ArrayList<>();
			List<Double> fullWeights = new ArrayList<>();
			double sumWeight = 0d;
			for (int i=fullBranchIndex; i<fullTree.size(); i++) {
				LogicTreeBranch<?> fullBranch = fullTree.getBranch(i);
				boolean match = true;
				for (int l=0; match && l<subsetBranch.size(); l++)
					match = fullBranch.hasValue(subsetBranch.getValue(l));
				if (match) {
					fullBranchIndex = i;
					fullMatches.add(fullBranch);
					double weight = fullTree.getBranchWeight(i);
					fullWeights.add(weight);
					sumWeight += weight;
				} else {
					break;
				}
			}
			System.out.println("\tFound "+fullMatches.size()+" matching full branches");
			Preconditions.checkState(!fullMatches.isEmpty());
			
			List<CompletableFuture<GriddedGeoDataSet>> loadAvgFutures = new ArrayList<>(mapPrefixes.length);
			
			for (String mapPrefix : mapPrefixes)
				loadAvgFutures.add(CompletableFuture.supplyAsync(new LoadAverageSupplier(
						zip, gridReg, mapPrefix, fullMatches, fullWeights, sumWeight)));
			
			GriddedGeoDataSet[] avgMaps = new GriddedGeoDataSet[mapPrefixes.length];
			for (int i=0; i<avgMaps.length; i++)
				avgMaps[i] = loadAvgFutures.get(i).join();
			System.out.println("\tDone loading");
			
			if (writeFuture != null && !writeFuture.isDone()) {
				System.out.println("\tWaiting on previous write");
				writeFuture.join();
			}
			String outDirName = subsetBranch.buildFileName();
			writeFuture = CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					try {
						zout.putNextEntry(new ZipEntry(outDirName+"/"));
						zout.flush();
						zout.closeEntry();
						for (int i=0; i<avgMaps.length; i++) {
							zout.putNextEntry(new ZipEntry(outDirName+"/"+mapPrefixes[i]));
							GriddedGeoDataSet.writeXYZStream(avgMaps[i], zout);
							zout.flush();
							zout.closeEntry();
						}
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
				}
			});
			
			fullBranchIndex++;
		}
		
		if (writeFuture != null && !writeFuture.isDone()) {
			System.out.println("\tWaiting on previous write");
			writeFuture.join();
		}
		
		System.out.println("Copying over extras");
		Enumeration<? extends ZipEntry> entries = zip.entries();
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			
			if (entry.getName().startsWith("level_choice_maps") || entry.getName().startsWith("mean_")) {
				System.out.println("\tCopying: "+entry.getName());
				zout.putNextEntry(new ZipEntry(entry.getName()));
				if (!entry.getName().endsWith("/")) {
					InputStream is = zip.getInputStream(entry);
					ByteStreams.copy(new BufferedInputStream(zip.getInputStream(entry)), zout);
					is.close();
				}
				zout.flush();
				zout.closeEntry();
			}
		}
		
		zip.close();
		
		System.out.println("Writing gridded region");
		zout.putNextEntry(new ZipEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME));
		Feature gridFeature = gridReg.toFeature();
		Feature.write(gridFeature, new OutputStreamWriter(zout));
		zout.closeEntry();
		
		System.out.println("Writing logic tree");
		zout.putNextEntry(new ZipEntry("logic_tree.json"));
		subsetTree.writeToStream(new BufferedOutputStream(zout));
		zout.closeEntry();
		
		zout.close();
	}
	
	private static class LoadAverageSupplier implements Supplier<GriddedGeoDataSet> {
		
		private ZipFile zip;
		private String mapPrefix;
		private List<LogicTreeBranch<?>> fullBranches;
		private List<Double> fullWeights;
		private double sumWeight;
		private GriddedRegion gridReg;

		public LoadAverageSupplier(ZipFile zip, GriddedRegion gridReg, String mapPrefix,
				List<LogicTreeBranch<?>> fullBranches, List<Double> fullWeights, double sumWeight) {
					this.zip = zip;
					this.gridReg = gridReg;
					this.mapPrefix = mapPrefix;
					this.fullBranches = fullBranches;
					this.fullWeights = fullWeights;
					this.sumWeight = sumWeight;
			
		}

		@Override
		public GriddedGeoDataSet get() {
			GriddedGeoDataSet avgMap = new GriddedGeoDataSet(gridReg);
			
			for (int i=0; i<fullBranches.size(); i++) {
				double weight = fullWeights.get(i)/sumWeight;
				String dirName = fullBranches.get(i).buildFileName();
				
				String entryName = dirName+"/"+mapPrefix;
				ZipEntry entry = zip.getEntry(entryName);
				Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
				
				GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
				int index = 0;
				try {
					BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
					String line = bRead.readLine();
					while (line != null) {
						line = line.trim();
						if (!line.startsWith("#")) {
							StringTokenizer tok = new StringTokenizer(line);
							double lon = Double.parseDouble(tok.nextToken());
							double lat = Double.parseDouble(tok.nextToken());
							double val = Double.parseDouble(tok.nextToken());
							Location loc = new Location(lat, lon);
							Preconditions.checkState(LocationUtils.areSimilar(loc, gridReg.getLocation(index)));
							xyz.set(index++, val);
						}
						line = bRead.readLine();
					}
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				Preconditions.checkState(index == gridReg.getNodeCount());
				
				for (int n=0; n<gridReg.getNodeCount(); n++)
					avgMap.set(n, avgMap.get(n) + weight*xyz.get(n));
			}
			return avgMap;
		}
		
	}

}
