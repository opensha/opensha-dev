package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.XYZClosestPointFinder;

public class CSVSimMeshConvert {
	
	private static void handleFile(String dir, String fName) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(new File(dir + File.separator + fName), true);
		
		System.out.println("loaded "+csv.getNumRows()+" from "+fName);
		
		GeoDataSet pgv = new ArbDiscrGeoDataSet(true);
		GeoDataSet pga = new ArbDiscrGeoDataSet(true);
		GeoDataSet sa03 = new ArbDiscrGeoDataSet(true);
		GeoDataSet sa10 = new ArbDiscrGeoDataSet(true);
		
		System.out.println("populating GeoDataSet's");
		for (int i=0; i<csv.getNumRows(); i++) {
			List<String> line = csv.getLine(i);
			double lat = Double.parseDouble(line.get(0));
			double lon = Double.parseDouble(line.get(1));
			double pgvVal = Double.parseDouble(line.get(2));
			double pgaVal = Double.parseDouble(line.get(3));
			double sa03Val = Double.parseDouble(line.get(4));
			double sa10Val = Double.parseDouble(line.get(5));
			Location loc = new Location(lat, lon);
			
			pgv.set(loc, pgvVal);
			pga.set(loc, pgaVal);
			sa03.set(loc, sa03Val);
			sa10.set(loc, sa10Val);
		}
		
		double spacing = 0.016667;
//		double tolerance = spacing * 2;
		double tolerance = 10; // KM
		
		Location topLeft = new Location(pgv.getMaxLat(), pgv.getMinLon());
		Location bottomRight = new Location(pgv.getMinLat(), pgv.getMaxLon());
		System.out.println("Creating gridded region with corners:\n"+topLeft+"\n"+bottomRight);
		GriddedRegion region = new GriddedRegion(topLeft, bottomRight, spacing, null);
		
		XYZClosestPointFinder xyz = new XYZClosestPointFinder(pgv);
		
		GriddedGeoDataSet pgvGridded = new GriddedGeoDataSet(region, true);
		GriddedGeoDataSet pgaGridded = new GriddedGeoDataSet(region, true);
		GriddedGeoDataSet sa03Gridded = new GriddedGeoDataSet(region, true);
		GriddedGeoDataSet sa10Gridded = new GriddedGeoDataSet(region, true);
		
		int tot = region.getNodeCount();
		System.out.println("matching to "+tot+" grid points");
		int nans = 0;
		int cnt = 0;
		for (Location gridLoc : region.getNodeList()) {
			Location dataLoc = xyz.getClosestLoc(gridLoc, tolerance);
			
			if (dataLoc == null) {
				pgvGridded.set(gridLoc, Double.NaN);
				pgaGridded.set(gridLoc, Double.NaN);
				sa03Gridded.set(gridLoc, Double.NaN);
				sa10Gridded.set(gridLoc, Double.NaN);
				nans++;
			} else {
				pgvGridded.set(gridLoc, pgv.get(dataLoc));
				pgaGridded.set(gridLoc, pga.get(dataLoc));
				sa03Gridded.set(gridLoc, sa03.get(dataLoc));
				sa10Gridded.set(gridLoc, sa10.get(dataLoc));
			}
			cnt++;
			if (cnt % 1000 == 0)
				System.out.println("processed "+cnt+"/"+tot+" points ("+nans+" nans, "+
						((float)((double)nans / (double)cnt) * 100f)+" %)");
		}
		System.out.println("processed "+cnt+"/"+tot+" points ("+nans+" nans, "+
				((float)((double)nans / (double)cnt) * 100f)+" %)");
		
		System.out.println("writing results");
		GriddedGeoDataSet.writeXYZFile(pgvGridded, dir+File.separator+"pgv_map_data.txt");
		GriddedGeoDataSet.writeXYZFile(pgaGridded, dir+File.separator+"pga_map_data.txt");
		GriddedGeoDataSet.writeXYZFile(sa03Gridded, dir+File.separator+"sa03_map_data.txt");
		GriddedGeoDataSet.writeXYZFile(sa10Gridded, dir+File.separator+"sa10_map_data.txt");
		
		System.out.println("DONE");
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		String dir = "/home/kevin/OpenSHA/hope_2011_01/final2_3";
//		String dir = "/home/kevin/OpenSHA/hope_2011_01/g6d3";
//		String dir = "/home/kevin/OpenSHA/hope_2011_01/g7d1_final";
//		String dir = "/home/kevin/OpenSHA/hope_2011_01/v1d3_final";
//		String dir = "/home/kevin/OpenSHA/hope_2011_01/final2_4";
		String dir = "/home/kevin/OpenSHA/hope_2011_01/M8_final";
//		String fName = "final2_3.csv";
//		String fName = "g6d3.csv";
//		String fName = "g7d1_final.csv";
//		String fName = "v1d3_final.csv";
//		String fName = "final2_4.csv";
		String fName = "M8_final.csv";
		handleFile(dir, fName);
	}

}
