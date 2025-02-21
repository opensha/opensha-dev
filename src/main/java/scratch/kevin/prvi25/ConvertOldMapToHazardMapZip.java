package scratch.kevin.prvi25;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class ConvertOldMapToHazardMapZip {
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Downloads/prvi_2003_static-total/map");
		
		List<File> inputFiles = new ArrayList<>();
		List<String> outputNames = new ArrayList<>();
		
		inputFiles.add(new File(dir, "map_prvi_2003_static-total_vs760_PGA_00475yrs.gmt"));
		outputNames.add("map_pga_TEN_IN_50.txt");
		
		inputFiles.add(new File(dir, "map_prvi_2003_static-total_vs760_PGA_02475yrs.gmt"));
		outputNames.add("map_pga_TWO_IN_50.txt");
		
		inputFiles.add(new File(dir, "map_prvi_2003_static-total_vs760_SA0P2_00475yrs.gmt"));
		outputNames.add("map_0.2s_TEN_IN_50.txt");
		
		inputFiles.add(new File(dir, "map_prvi_2003_static-total_vs760_SA0P2_02475yrs.gmt"));
		outputNames.add("map_0.2s_TWO_IN_50.txt");
		
		inputFiles.add(new File(dir, "map_prvi_2003_static-total_vs760_SA1P0_00475yrs.gmt"));
		outputNames.add("map_1.0s_TEN_IN_50.txt");
		
		inputFiles.add(new File(dir, "map_prvi_2003_static-total_vs760_SA1P0_02475yrs.gmt"));
		outputNames.add("map_1.0s_TWO_IN_50.txt");
		
		Region reg = new Region(new Location(17.5, -67.5), new Location(19, -64.5));
		GriddedRegion gridReg = new GriddedRegion(reg, 0.01, GriddedRegion.ANCHOR_0_0);
		
		File outputFile = new File(dir, "results_hazard.zip");
		ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		for (int p=0; p<inputFiles.size(); p++) {
			File inputFile = inputFiles.get(p);
			String outputName = outputNames.get(p);
			
			System.out.println(inputFile.getName()+" -> "+outputName);
			
			int numSkipped = 0;
			int numSet = 0;
			
			GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg);
			// set to NaN
			xyz.scale(Double.NaN);
			
			// these CSVs have annoying padding extra empty columns
			for (String line : Files.readLines(inputFile, Charset.defaultCharset())) {
				line = line.trim();
				if (line.startsWith("#") || line.isEmpty())
					continue;
				Preconditions.checkState(line.contains(","));
				String[] split = line.split(",");
				double lon = Double.parseDouble(split[0]);
				double lat = Double.parseDouble(split[1]);
				double val = Double.parseDouble(split[2]);
				
				Location loc = new Location(lat, lon);
				int index = gridReg.indexForLocation(loc);
				if (index < 0) {
					numSkipped++;
					continue;
				}
				Preconditions.checkState(Double.isNaN(xyz.get(index)));
				xyz.set(index, val);
				numSet++;
			}
			System.out.println("\tSkipped "+numSkipped+"/"+xyz.size());
			System.out.println("\tSet "+numSet+"/"+xyz.size());
			
			zout.putNextEntry(new ZipEntry(outputName));
			ArbDiscrGeoDataSet.writeXYZStream(xyz, zout);
			zout.flush();
			zout.closeEntry();
		}
		
		zout.putNextEntry(new ZipEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME));
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(zout));
		Feature.write(gridReg.toFeature(), writer);
		writer.flush();
		zout.closeEntry();
		
		zout.close();
	}

}
