package scratch.kevin.sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

public class HazardConvergenceCalcs {

	public static void main(String[] args) throws IOException {
		String treeFileName = "logic_tree_analysis.json";
		String hazardFileName = "results_hazard.zip";
		
		File refMCSDir = new File(PaperPaths.INVS_DIR, "2026_07_17-nshm27-AMSAM-20000samples-mcs");
		LogicTree<?> refMCSTree = LogicTree.read(new File(refMCSDir, treeFileName));
		File refMCSHazardZip = new File(refMCSDir, hazardFileName);
		
		GriddedRegion gridReg = GriddedRegion.fromFeature(Feature.read(new File(refMCSDir, "gridded_region.geojson")));
		
		double period = 0;
		String periodName = "PGA";
		String periodPrefix = "pga";
		
//		double period = 1;
//		String periodName = "1s SA";
//		String periodPrefix = "1s_sa";
		
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		
		System.out.println("Ref has "+refMCSTree.size()+" branches");
		
		Table<SamplingMethod, Integer, List<File>> runDirs = HashBasedTable.create();
		
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 512, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-512samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-512samples-sobol_scrambled-unique_seed")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 1024, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-1024samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-1024samples-sobol_scrambled-unique_seed")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 2048, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-2048samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-2048samples-sobol_scrambled-unique_seed")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 4096, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-4096samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-4096samples-sobol_scrambled-unique_seed")
				));
		
		ModelHazarMaps refMCSMaps = loadMaps(refMCSHazardZip, refMCSTree, gridReg, period, rp);
	}
	
	private static String mapFilePrefix(double period, ReturnPeriods rp) {
		String perStr = period == 0d ? "pga" : (float)period+"s";
		return perStr+"_"+rp.name();
	}
	
	private static GriddedGeoDataSet readMap(GriddedRegion gridReg, InputStream is) throws IOException {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(is));
		String line = bRead.readLine();
		int index = 0;
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
		Preconditions.checkState(index == gridReg.getNodeCount());
		bRead.close();
		return xyz;
	}
	
	private static ModelHazarMaps loadMaps(File hazardZip, LogicTree<?> tree, GriddedRegion gridReg,
			double period, ReturnPeriods rp) throws ZipException, IOException {
		System.out.println("Loading maps from "+hazardZip.getAbsolutePath());
		try (ZipFile zip = new ZipFile(hazardZip)) {
			String suffix = mapFilePrefix(period, rp)+".txt";
			String meanEntryName = "mean_map_"+suffix;
			ZipEntry meanEntry = zip.getEntry(meanEntryName);
			Preconditions.checkNotNull(meanEntry, "Entry doesn't exist in %s: %s", hazardZip.getAbsolutePath(), meanEntryName);
			GriddedGeoDataSet meanMap = readMap(gridReg, zip.getInputStream(meanEntry));
			
			GriddedGeoDataSet[] individual = new GriddedGeoDataSet[tree.size()];
			for (int i=0; i<individual.length; i++) {
				LogicTreeBranch<?> branch = tree.getBranch(i);
				String mapName = branch.buildFileName()+"/map_"+suffix;
				ZipEntry mapEntry = zip.getEntry(mapName);
				Preconditions.checkNotNull(mapEntry, "Entry doesn't exist in %s: %s", hazardZip.getAbsolutePath(), mapName);
				individual[i] = readMap(gridReg, zip.getInputStream(mapEntry));
			}
			
			System.out.println("\tLoaded mean & "+individual.length+" individual");
			
			return new ModelHazarMaps(meanMap, individual);
		}
	}
	
	record ModelHazarMaps(GriddedGeoDataSet mean, GriddedGeoDataSet[] individual) {}

}
