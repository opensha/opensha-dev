package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;

public class WrapperCompareToExternal {

	public static void main(String[] args) throws IOException {
		File wrapperCalcDir = new File("/data/kevin/nshm23/batch_inversions/2024_02_08-nshm23-wrapped-conus-hazard-ask2014-0.1deg");
//		File wrapperCalcDir = new File("/data/kevin/nshm23/batch_inversions/2024_02_13-nshm23-wrapped-conus-hazard-wrapedask2014-0.1deg");
		double[] periods = {0d, 1d};
		ReturnPeriods[] rps = {ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50};
		
		Region region = NSHM23_RegionLoader.loadFullConterminousUS();
		GriddedRegion gridReg = new GriddedRegion(region, 0.1, GriddedRegion.ANCHOR_0_0);
		
		File extCalcDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
//				+ "conus-6b4-nshmp-haz-0p1-vs760-20240208-afaff93cd918f5/vs30-760");
//		File outputDir = new File(wrapperCalcDir, "hazard_comparison_nsmp_haz");
				+ "conus-6b4-nshmp-haz-grid_smooth_OFF-0p1-vs760-20240213-bbc02463b2e1ff/vs30-760");
		File outputDir = new File(wrapperCalcDir, "hazard_comparison_nsmp_haz_no_smooth");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		for (ReturnPeriods rp : rps)
			System.out.println(rp+" is "+rp.returnPeriod+" years");
		
		ZipFile wrapperZip = new ZipFile(new File(wrapperCalcDir, "results_hazard.zip"));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setFaultSections(NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.WUS_FM_v3));
		mapMaker.setSectOutlineChar(null);
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));

		Color transparent = new Color(255, 255, 255, 0);
//		CPT pDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance(), 10d, 50d);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(transparent);
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2d, 0.2d);
		diffCPT.setNanColor(transparent);
		
		for (double period : periods) {
			String extDirName;
			String perLabel;
			String mapPrefix;
			if (period == 0d) {
				extDirName = "PGA";
				mapPrefix = "map_pga";
				perLabel = "PGA";
			} else {
				perLabel = (float)period+"s SA";
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
				
				String mapFileName = mapPrefix+"_"+rp.name()+".txt";
				ZipEntry mapEntry = wrapperZip.getEntry(mapFileName);
				InputStream is = wrapperZip.getInputStream(mapEntry);
				Preconditions.checkNotNull(is, "IS is null for %s", wrapperZip);
				
				GriddedGeoDataSet wrapperXYZ = new GriddedGeoDataSet(gridReg, false);
				BufferedReader bRead = new BufferedReader(new InputStreamReader(wrapperZip.getInputStream(mapEntry)));
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
						wrapperXYZ.set(index++, val);
					}
					line = bRead.readLine();
				}
				Preconditions.checkState(index == gridReg.getNodeCount());
				
				GriddedGeoDataSet pDiff = WUS_HazardChangePageGen.mapPDiff(wrapperXYZ, extXYZ);
				GriddedGeoDataSet diff = WUS_HazardChangePageGen.mapDiff(wrapperXYZ, extXYZ);
				
				String prefix = mapPrefix+"_"+rp.name();
				String label = perLabel+", "+rp.label;
				
				mapMaker.plotXYZData(diff, diffCPT, "Wrapper - NSHMP-Haz, "+label);
				mapMaker.plot(outputDir, prefix+"_diff", " ");
				mapMaker.plotXYZData(pDiff, pDiffCPT, "Wrapper vs NSHMP-Haz, % Difference, "+label);
				mapMaker.plot(outputDir, prefix+"_pDiff", " ");
			}
		}
		
		wrapperZip.close();
	}

}
