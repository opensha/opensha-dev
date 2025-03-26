package scratch.kevin.pointSources;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.ZipFile;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.util.ComparablePairing;

import com.google.common.base.Preconditions;

public class HazardZipBiggestDifferenceFinder {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions");
		String zipName = "results_hazard_INCLUDE.zip";
		String mapName = "map_pga_TWO_IN_50.txt";
		
		int numLocs = 20;
		
		File mainDir = new File(invsDir, "2025_03_25-nshm23_pt_src_tests-pt_src_corr-analytical_5pt");
		File compDir = new File(invsDir, "2025_03_25-nshm23_pt_src_tests-finite_crosshair");
		
		ZipFile mainZip = new ZipFile(new File(mainDir, zipName));
		ZipFile compZip = new ZipFile(new File(compDir, zipName));
		
		Feature gridRegFeature = Feature.read(new InputStreamReader(mainZip.getInputStream(mainZip.getEntry("gridded_region.geojson"))));
		GriddedRegion gridReg = GriddedRegion.fromFeature(gridRegFeature);
		
		GriddedGeoDataSet mainMap = readMap(mainZip, mapName, gridReg);
		GriddedGeoDataSet compMap = readMap(compZip, mapName, gridReg);
		
		Map<Integer, Double> pDiffs = new HashMap<>(gridReg.getNodeCount());
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			double mainVal = mainMap.get(i);
			double compVal = compMap.get(i);
			if (mainVal == 0d || compVal == 0d)
				continue;
			double pDiff = 100d*(mainVal-compVal)/compVal;
			pDiffs.put(i, pDiff);
		}
		List<Integer> sorted = ComparablePairing.getSortedData(pDiffs);
		
		System.out.println("Sorted % differences:");
		for (int i=0; i<numLocs && i<sorted.size(); i++) {
			Location loc = gridReg.getLocation(sorted.get(i));
			double pDiff = pDiffs.get(sorted.get(i));
			System.out.println((float)loc.lat+", "+(float)loc.lon+":\t"+(float)pDiff+" %");
		}
		System.out.println("...");
		for (int i=sorted.size()-numLocs-1; i<sorted.size(); i++) {
			if (i < 0)
				continue;
			Location loc = gridReg.getLocation(sorted.get(i));
			double pDiff = pDiffs.get(sorted.get(i));
			System.out.println((float)loc.lat+", "+(float)loc.lon+":\t"+(float)pDiff+" %");
		}
	}
	
	private static GriddedGeoDataSet readMap(ZipFile zip, String mapName, GriddedRegion gridReg) throws IOException {
		BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(zip.getEntry(mapName))));
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
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
		return xyz;
	}

}
