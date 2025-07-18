package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class SeisPDFCheckAndInterfaceDepthCutNormalize {

	public static void main(String[] args) throws IOException {
		File pdfsMainDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/seismicity/spatial_seis_pdfs/2025_07_18");
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		for (PRVI25_SeismicityRegions seisReg : PRVI25_SeismicityRegions.values()) {
			File regDir = new File(pdfsMainDir, seisReg.name());
			
			System.out.println("Region: "+seisReg);
			
			Region region = seisReg.load();
			GriddedRegion gridReg = new GriddedRegion(region, 0.1d, GriddedRegion.ANCHOR_0_0);
			
			boolean isInterface = seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE || seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE;
			
			GriddedGeoDataSet depths = null;
			GriddedGeoDataSet strikes = null;
			GriddedGeoDataSet dips = null;
			if (seisReg != PRVI25_SeismicityRegions.CRUSTAL) {
				depths = PRVI25_GridSourceBuilder.loadSubductionDepths(seisReg);
				if (isInterface) {
					strikes = PRVI25_GridSourceBuilder.loadSubductionStrikes(seisReg);
					dips = PRVI25_GridSourceBuilder.loadSubductionDips(seisReg);
				}

//				double minLatValidDepth = Double.POSITIVE_INFINITY;
//				double maxLatValidDepth = Double.NEGATIVE_INFINITY;
//				double minLonValidDepth = Double.POSITIVE_INFINITY;
//				double maxLonValidDepth = Double.NEGATIVE_INFINITY;
//				for (int i=0; i<depths.size(); i++) {
//					if (Double.isFinite(depths.get(i))) {
//						Location loc = depths.getLocation(i);
//						minLatValidDepth = Math.min(minLatValidDepth, loc.lat);
//						maxLatValidDepth = Math.max(maxLatValidDepth, loc.lat);
//						minLonValidDepth = Math.min(minLonValidDepth, loc.lon);
//						maxLonValidDepth = Math.max(maxLonValidDepth, loc.lon);
//					}
//				}
//
//				System.out.println("\tLat range with valid depths: ["+(float)minLatValidDepth+", "+(float)maxLatValidDepth+"]");
//				System.out.println("\tLon range with valid depths: ["+(float)minLonValidDepth+", "+(float)maxLonValidDepth+"]");

				MinMaxAveTracker depthRange = new MinMaxAveTracker();
				MinMaxAveTracker strikeRange = new MinMaxAveTracker();
				MinMaxAveTracker dipRange = new MinMaxAveTracker();
				
				for (int i=0; i<gridReg.getNodeCount(); i++) {
					Location loc = gridReg.getLocation(i);
					int depthIndex = depths.indexOf(loc);
					Preconditions.checkState(depthIndex >= 0, "No depth point for location: %s. %s", i, loc);
					Location depthLoc = depths.getLocation(depthIndex);
					Preconditions.checkState(LocationUtils.areSimilar(loc, depthLoc),
							"Depth point doesn't match: %s. %s != %s. %s", i, loc, depthIndex, depthLoc);
					Preconditions.checkState(depthIndex >= 0, "No depth data for location: %s. %s", i, loc);
					double depth = depths.get(depthIndex);
					Preconditions.checkState(Double.isFinite(depth), "Bad depth for loc %s: %s", loc, depth);
					depthRange.addValue(depth);
					if (isInterface) {
						double strike = strikes.get(depthIndex);
						Preconditions.checkState(Double.isFinite(strike), "Bad strike for loc %s: %s", loc, strike);
						strikeRange.addValue(strike);
						double dip = dips.get(depthIndex);
						Preconditions.checkState(Double.isFinite(dip), "Bad depth for loc %s: %s", loc, dip);
						dipRange.addValue(dip);
					}
				}
				System.out.println("\tMapped depths:\t"+depthRange);
				if (isInterface) {
					System.out.println("\tMapped strikes:\t"+strikeRange);
					System.out.println("\tMapped dips:\t"+dipRange);
				}
			}
			
			double sumInterfaceMask = 0d;
			int numInterfaceMasks = 0;
			
			for (PRVI25_DeclusteringAlgorithms decluster : PRVI25_DeclusteringAlgorithms.values()) {
				if (decluster.getNodeWeight(null) == 0d)
					continue;
				for (PRVI25_SeisSmoothingAlgorithms smooth : PRVI25_SeisSmoothingAlgorithms.values()) {
					if (smooth.getNodeWeight(null) == 0d)
						continue;
					File testOrigFile = new File(regDir, decluster.getFilePrefix()+"_"+smooth.getFilePrefix()+"_original.csv");
					File inputFile;
					if (testOrigFile.exists())
						inputFile = testOrigFile;
					else
						inputFile = new File(regDir, decluster.getFilePrefix()+"_"+smooth.getFilePrefix()+".csv");
					System.out.println("\tParsing "+inputFile.getName());
					
					CSVFile<String> csv = CSVFile.readFile(inputFile, true);
					
					GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
					double sum = 0d;
					int numMapped = 0;
					
					for (int row=0; row<csv.getNumRows(); row++) {
						double lon = csv.getDouble(row, 0);
						double lat = csv.getDouble(row, 1);
						double val = csv.getDouble(row, 2);
						sum += val;
						Location loc = new Location(lat, lon);
						int gridIndex = gridReg.indexForLocation(loc);
						if (gridIndex >= 0) {
							numMapped++;
							Preconditions.checkState(xyz.get(gridIndex) == 0d);
							xyz.set(gridIndex, val);
						}
					}
					double sumMapped = xyz.getSumZ();
					System.out.println("\t\ttotWeight="+(float)sum+";\tmappedWeight="+(float)sumMapped+"; mapping results:");
					System.out.println("\t\t\t\t"+numMapped+"/"+csv.getNumRows()+" ("
							+pDF.format((double)numMapped/(double)csv.getNumRows())+") of locations from input CSV mapped");
					System.out.println("\t\t\t\t"+numMapped+"/"+gridReg.getNodeCount()+" ("
							+pDF.format((double)numMapped/(double)gridReg.getNodeCount())+") of gridded region mapped");
					Preconditions.checkState(Precision.equals(sumMapped, 1d, 0.01),
							"PDF (%s) doesn't sum to 1 when mapped to region: sum=%s, sumMapped=%s", inputFile.getName(), (float)sum, (float)sumMapped);
					
					if (isInterface) {
						float maxDepth = seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE ?
								(float)PRVI25_GridSourceBuilder.INTERFACE_CAR_MAX_DEPTH : (float)PRVI25_GridSourceBuilder.INTERFACE_MUE_MAX_DEPTH;
						double maskedWeight = 0d;
						int maskedLocs = 0;
						for (int i=0; i<xyz.size(); i++) {
							if (xyz.get(i) > 0d) {
								Location loc = xyz.getLocation(i);
								int depthIndex = depths.indexOf(loc);
								Preconditions.checkState(depthIndex >= 0);
								Preconditions.checkState(LocationUtils.areSimilar(loc, depths.getLocation(depthIndex)));
								double depth = depths.get(depthIndex);
								if ((float)depth > maxDepth) {
									maskedLocs++;
									maskedWeight += xyz.get(i);
									xyz.set(i, 0d);
								}
							}
						}
						
						sumInterfaceMask += maskedWeight;
						numInterfaceMasks++;
						System.out.println("\t\tMasking out "+(float)maskedWeight+" across "+maskedLocs);
						
						File origFileCopy = new File(regDir, decluster.getFilePrefix()+"_"+smooth.getFilePrefix()+"_original.csv");
						if (!origFileCopy.exists())
							Files.copy(inputFile, origFileCopy);
						
						xyz.scale(1d/xyz.getSumZ());
						
						CSVFile<String> outCSV = new CSVFile<>(true);
						for (int i=0; i<xyz.size(); i++) {
							Location loc = xyz.getLocation(i);
							outCSV.addLine((float)loc.lon+"", (float)loc.lat+"", xyz.get(i)+"");
						}
						
						outCSV.writeToFile(new File(regDir, decluster.getFilePrefix()+"_"+smooth.getFilePrefix()+".csv"));
					}
				}
				if (isInterface) {
					System.out.println("\tAverage mask: "+(float)(sumInterfaceMask/(double)numInterfaceMasks));
				}
			}
			
			System.out.println();
			System.out.println();
		}
		
	}

}
