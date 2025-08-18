package scratch.kevin.nshm23.wrapper;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

public class RewriteExternalHazard {

	public static void main(String[] args) throws IOException {
		// utility to convert nshmp-haz runs to our hazard results zip file format
		
		File invsDir = new File("/data/kevin/nshm23/batch_inversions");
		
//		File extCalcDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
//				+ "prvi-t2-2003-ERF-2025-vB1-GMMs-2025conf-0p01-vs760-20241213-49f58ecb02d600/vs30-760");
//		File outputDir = new File(invsDir, "2024_12_13-nshmp-haz-external-prvi-2b1-prvi25gmms-vs760");
		
		File extCalcDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
				+ "prvi-t14-2003-ERF-2025-GMMs-v1.7.10-2025conf-0p01ext-vs760-20250812-6aa514571a45ae/vs30-760");
		File outputDir = new File(invsDir, "2025_08_15-nshmp-haz-external-prvi03-t14-prvi25gmms-vs760");
		
		GriddedRegion gridReg = new GriddedRegion(
				PRVI25_RegionLoader.loadPRVI_MapExtents(),
				0.01, GriddedRegion.ANCHOR_0_0);

		double[] periods = {0d, 0.2, 1d, 5d};
		
//		File extCalcDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
//				+ "conus-6b4-nshmp-haz-lib1415-grid_smooth_OFF-0p1-vs760-20240307-087a4937807eb5/vs30-760");
//		File outputDir = new File(invsDir, "2024_03_08-nshmp-haz-external-conus-6b4-ask2014-vs760-grid_smooth_OFF");
//				+ "conus-6b4-nshmp-haz-grid_smooth_optimize_OFF-0p1-vs760-20240213-195bddb0d73730/vs30-760");
//		File outputDir = new File(invsDir, "2024_02_15-nshmp-haz-external-conus-6b4-ask2014-vs760-grid_smooth_optimize_OFF");
//				+ "conus-6b4-nshmp-haz-0p1-vs760-20240208-afaff93cd918f5/vs30-760");
//		File outputDir = new File(invsDir, "2024_02_08-nshmp-haz-external-conus-6b4-ask2014-vs760");
//				+ "conus-6b4-nshmp-haz-0p1-vs760-20240228-ef0f647c24c8d1/vs30-760");
//		File outputDir = new File(invsDir, "2024_02_28-nshmp-haz-external-conus-6b4-ask2014-vs760");
		
//		GriddedRegion gridReg = new GriddedRegion(
//				NSHM23_RegionLoader.loadFullConterminousUS(),
//				0.1, GriddedRegion.ANCHOR_0_0);
		
//		double[] periods = {0d, 1d};
		ReturnPeriods[] rps = {ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50};
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(
				new File(outputDir, "results_hazard.zip"))));
		
		for (double period : periods) {
			String extDirName;
			String mapPrefix;
			if (period == 0d) {
				extDirName = "PGA";
				mapPrefix = "map_pga";
			} else {
				extDirName = "SA"+((float)period+"").replace('.', 'P');
				mapPrefix = "map_"+(float)period+"s";
			}
			File extSubDir = new File(extCalcDir, extDirName);
			Preconditions.checkState(extSubDir.exists(), "Ext period dir doesn't exist: %s", extSubDir.getAbsolutePath());
			
			CSVFile<String> extCSV = CSVFile.readFile(new File(extSubDir, "map.csv"), true);
			
			for (ReturnPeriods rp : rps) {
				double years = rp.returnPeriod;
				double minDiff = Double.POSITIVE_INFINITY;
				int closestCol = -1;
				for (int col=2; col<extCSV.getNumCols(); col++) {
					double csvYears = extCSV.getDouble(0, col);
					double diff = Math.abs(csvYears - years);
					if (diff < minDiff) {
						closestCol = col;
						minDiff = diff;
					}
				}
				System.out.println("Matched "+rp+" with CSV column "+extCSV.get(0, closestCol));
				Preconditions.checkState(minDiff < 5d, "RP mimatch! |ours (%s = %s) - theirs (%s)| = %s",
						rp, rp.returnPeriod, extCSV.get(0, closestCol), minDiff);
				
				int numSkipped = 0;
				int numMatched = 0;
				GriddedGeoDataSet extXYZ = new GriddedGeoDataSet(gridReg);
				for (int i=0; i<extXYZ.size(); i++)
					extXYZ.set(i, Double.NaN);
				for (int row=1; row<extCSV.getNumRows(); row++) {
					double lon = extCSV.getDouble(row, 0);
					double lat = extCSV.getDouble(row, 1);
					double val = extCSV.getDouble(row, closestCol);
					Location loc = new Location(lat, lon);
					int locIndex = gridReg.indexForLocation(loc);
					if (locIndex < 0) {
						numSkipped++;
						continue;
					} else {
						Preconditions.checkState(Double.isNaN(extXYZ.get(locIndex)));
						extXYZ.set(locIndex, val);
						numMatched++;
					}
				}
				System.out.println("Filled in "+numMatched+"/"+extXYZ.size()+" points (skipped "+numSkipped+")");
				
				String destName = mapPrefix+"_"+rp.name()+".txt";
				ZipEntry entry = new ZipEntry(destName);
				zout.putNextEntry(entry);
				GriddedGeoDataSet.writeXYZStream(extXYZ, zout);
				zout.flush();
			}
		}
		// write the region
		ZipEntry entry = new ZipEntry("gridded_region.geojson");
		zout.putNextEntry(entry);
		OutputStreamWriter writer = new OutputStreamWriter(zout);
		Feature.write(gridReg.toFeature(), writer);
		writer.flush();
		zout.close();
	}

}
