package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.FaultUtils.AngleAverager;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.SeismicityRegions;

import com.google.common.base.Preconditions;

public class SlabDepthCSVConverter {

	public static void main(String[] args) throws IOException {
		File inputDir = new File("/tmp/slab_depths");
		File outputDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/seismicity/depths");
		
		Map<SeismicityRegions, File> inputDepthCSVs = new HashMap<>();
		Map<SeismicityRegions, File> inputStrikeCSVs = new HashMap<>();
		
		inputDepthCSVs.put(SeismicityRegions.CAR_INTERFACE, new File(inputDir, "apdf_car_interface_gk_adN4.csv"));
		inputDepthCSVs.put(SeismicityRegions.CAR_INTRASLAB, new File(inputDir, "apdf_car_intraslab_gk_adN4.csv"));
		inputDepthCSVs.put(SeismicityRegions.MUE_INTERFACE, new File(inputDir, "apdf_mue_interface_gk_adN3.csv"));
		inputDepthCSVs.put(SeismicityRegions.MUE_INTRASLAB, new File(inputDir, "apdf_mue_intraslab_gk_adN2.csv"));
		
		inputStrikeCSVs.put(SeismicityRegions.CAR_INTERFACE, new File(inputDir, "strike_car_interface_gk_adN4.csv"));
		inputStrikeCSVs.put(SeismicityRegions.MUE_INTERFACE, new File(inputDir, "strike_mue_interface_gk_adN3.csv"));
		
		for (SeismicityRegions seisReg : inputDepthCSVs.keySet()) {
			System.out.println("Doing "+seisReg);
			Region reg = seisReg.load();
			GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
			GriddedGeoDataSet depthData = new GriddedGeoDataSet(gridReg);
			// initialized to nan
			for (int i=0; i<depthData.size(); i++)
				depthData.set(i, Double.NaN);
			
			File depthFile = inputDepthCSVs.get(seisReg);
			System.out.println("Reading "+depthFile.getName());
			CSVFile<String> depthCSV = CSVFile.readFile(depthFile, true);
			for (int row=1; row<depthCSV.getNumRows(); row++) {
				double lat = depthCSV.getDouble(row, 1);
				double lon = depthCSV.getDouble(row, 0);
				double depth = depthCSV.getDouble(row, 2);
				if (Double.isFinite(depth)) {
					Location loc = new Location(lat, lon);
					int index = gridReg.indexForLocation(loc);
					if (index < 0)
						System.err.println("WARNING: skipping location outside of region: "+loc);
					else
						depthData.set(index, depth);
				}
			}
			
			System.out.println("Depth range: "+(float)depthData.getMinZ()+" "+(float)depthData.getMaxZ());
			
			for (int i=0; i<depthData.size(); i++)
				Preconditions.checkState(Double.isFinite(depthData.get(i)),
						"No value found for index=%s, loc=%s", i, depthData.getLocation(i));
			
			File strikeFile = inputStrikeCSVs.get(seisReg);
			GriddedGeoDataSet strikeData = null;
			if (strikeFile != null) {
				AngleAverager strikeAvg = new AngleAverager();
				System.out.println("Reading "+strikeFile.getName());
				CSVFile<String> strikeCSV = CSVFile.readFile(strikeFile, true);
				strikeData = new GriddedGeoDataSet(gridReg);
				// initialized to nan
				for (int i=0; i<strikeData.size(); i++)
					strikeData.set(i, Double.NaN);
				for (int row=1; row<strikeCSV.getNumRows(); row++) {
					double lat = strikeCSV.getDouble(row, 1);
					double lon = strikeCSV.getDouble(row, 0);
					String strikeStr = strikeCSV.get(row, 2).toLowerCase().trim();
					double strike = strikeStr.equals("nan") ? Double.NaN : Double.parseDouble(strikeStr);
					if (Double.isFinite(strike)) {
						Location loc = new Location(lat, lon);
						int index = gridReg.indexForLocation(loc);
						strikeAvg.add(strike, 1d);
						if (index < 0)
							System.err.println("WARNING: skipping location outside of region: "+loc);
						else
							strikeData.set(index, strike);
					}
				}
				System.out.println("Average strike: "+(float)strikeAvg.getAverage());
				
//				for (int i=0; i<strikeData.size(); i++)
//					Preconditions.checkState(Double.isFinite(strikeData.get(i)),
//							"No value found for index=%s, loc=%s", i, strikeData.getLocation(i));
				int numMissing = 0;
				for (int i=0; i<strikeData.size(); i++)
					if (!Double.isFinite(strikeData.get(i)))
						numMissing++;
				if (numMissing > 0) {
					System.err.println("WARNING: filling in "+numMissing+" missing strike values");
					// figure out lat/lon gridding
					int numLat = gridReg.getNumLatNodes();
					int numLon = gridReg.getNumLonNodes();
					int[][] gridIndexes = new int[numLon][numLat];
					for (int i=0; i<numLon; i++)
						for (int j=0; j<numLat; j++)
							gridIndexes[i][j] = -1;
					int[] latIndexes = new int[gridReg.getNodeCount()];
					int[] lonIndexes = new int[gridReg.getNodeCount()];
					for (int i=0; i<gridReg.getNodeCount(); i++) {
						Location loc = gridReg.getLocation(i);
						int latIndex = gridReg.getLatIndex(loc);
						Preconditions.checkState(latIndex < numLat, "bad latIndex=%s for numLat=%s", latIndex, numLat);
						int lonIndex = gridReg.getLonIndex(loc);
						Preconditions.checkState(lonIndex < numLon, "bad lonIndex=%s for numLon=%s", lonIndex, numLon);
						gridIndexes[lonIndex][latIndex] = i;
						latIndexes[i] = latIndex;
						lonIndexes[i] = lonIndex;
					}
					
					// use nearest in the same column (same lon)
					for (int i=0; i<strikeData.size(); i++) {
						if (Double.isNaN(strikeData.get(i))) {
							int latIndex = latIndexes[i];
							int lonIndex = lonIndexes[i];
							int maxNumAway = Integer.MAX_VALUE;
							double closestValue = Double.NaN;
							for (int testLatIndex=0; testLatIndex<numLat; testLatIndex++) {
								if (gridIndexes[lonIndex][testLatIndex] >= 0 && Double.isFinite(strikeData.get(gridIndexes[lonIndex][testLatIndex]))) {
									// we have a point at this location
									int numAway = latIndex - testLatIndex;
									if (numAway < 0)
										numAway = -numAway;
									Preconditions.checkState(numAway > 0);
									if (numAway < maxNumAway) {
										maxNumAway = numAway;
										closestValue = strikeData.get(gridIndexes[lonIndex][testLatIndex]);
									}
								}
							}
							System.out.println("\tClosest in-column value for "+gridReg.locationForIndex(i)+" is "+maxNumAway+" cells away: "+closestValue);
							strikeData.set(i, closestValue);
						}
					}
				}
			}
			CSVFile<String> outCSV = new CSVFile<>(true);
			outCSV.addLine("Grid Index", "Latitude", "Longitude", "Depth (km)", "Strike");
			
			for (int i=0; i<depthData.size(); i++) {
				Location loc = depthData.getLocation(i);
				List<String> line = new ArrayList<>(outCSV.getNumCols());
				line.add(i+"");
				line.add((float)loc.lat+"");
				line.add((float)loc.lon+"");
				line.add((float)depthData.get(i)+"");
				if (strikeData == null)
					line.add("");
				else
					line.add((float)strikeData.get(i)+"");
				outCSV.addLine(line);
			}
			
			outCSV.writeToFile(new File(outputDir, seisReg.name()+".csv"));
		}
	}

}
