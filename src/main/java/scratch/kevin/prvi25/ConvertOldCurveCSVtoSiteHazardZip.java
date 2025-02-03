package scratch.kevin.prvi25;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.CSVFile;

public class ConvertOldCurveCSVtoSiteHazardZip {
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/tmp");
		
		List<File> inputFiles = new ArrayList<>();
		List<String> outputSuffixes = new ArrayList<>();
		
		inputFiles.add(new File(dir, "curves_pga_2023_total_mean.csv"));
		outputSuffixes.add("pga.csv");
		
		inputFiles.add(new File(dir, "curves_sa_0.2_2023_total_mean.csv"));
		outputSuffixes.add("sa_0.2.csv");
		
		inputFiles.add(new File(dir, "curves_sa_1_2023_total_mean.csv"));
		outputSuffixes.add("sa_1.0.csv");
		
		File outputFile = new File(dir, "prvi03_original_site_hazard.zip");
		ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		for (int p=0; p<inputFiles.size(); p++) {
			File inputFile = inputFiles.get(p);
			String suffix = outputSuffixes.get(p);
			
			// these CSVs have annoying padding extra empty columns
			CSVFile<String> inCSV = CSVFile.readFile(inputFile, false);
			int numRealCols = 0;
			for (String val : inCSV.getLine(0)) {
				if (val.isBlank())
					break;
				numRealCols++;
			}
			
			List<String> outHeader = new ArrayList<>();
			outHeader.add("Site Name");
			outHeader.add("Branch Index");
			outHeader.add("Branch Weight");
			for (int i=3; i<numRealCols; i++)
				outHeader.add(inCSV.get(0, i));
			for (int row=1; row<inCSV.getNumRows(); row++) {
				String siteName = inCSV.get(row, 0);
				CSVFile<String> outCSV = new CSVFile<>(true);
				outCSV.addLine(outHeader);
				List<String> line = new ArrayList<>();
				line.add(siteName);
				line.add("0"); // branch index
				line.add("1.0"); // branch weight
				for (int i=3; i<numRealCols; i++) {
					double rate = inCSV.getDouble(row, i);
					double prob = 1d - Math.exp(-rate);
					line.add((float)prob+"");
				}
				outCSV.addLine(line);
				
				String sitePrefix = siteName.replaceAll("\\W+", "_");
				zout.putNextEntry(new ZipEntry(sitePrefix+"_"+suffix));
				outCSV.writeToStream(zout);
				zout.flush();
				zout.closeEntry();
			}
		}
		
		zout.close();
	}

}
