package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;

public class PortfolioVs30Writer {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/eal/");
		
		File inputFile = new File(dir, "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio.csv");
		int latCol = 3;
		int lonCol = 4;
		int vsCol = 9;
		int delCol = 10;
		
		CSVFile<String> input = CSVFile.readFile(inputFile, true);
		LocationList locs = new LocationList();
		for (int row=1; row<input.getNumRows(); row++)
			locs.add(new Location(input.getDouble(row, latCol), input.getDouble(row, lonCol)));
		
		WillsMap2015 wills = new WillsMap2015();
		WaldAllenGlobalVs30 wald = new WaldAllenGlobalVs30();

		ArrayList<Double> willsVals = wills.getValues(locs);
		ArrayList<Double> waldVals = wald.getValues(locs);
		
		List<String> header = input.getLine(0);
		if (delCol >= 0) {
			header = new ArrayList<>(header);
			header.remove(delCol);
		}
		
		CSVFile<String> willsCSV = new CSVFile<>(true);
		CSVFile<String> waldCSV = new CSVFile<>(true);
		
		willsCSV.addLine(header);
		waldCSV.addLine(header);
		
		for (int i=0; i<locs.size(); i++) {
			List<String> line = input.getLine(i+1);
			
			List<String> willsLine = new ArrayList<>(line);
			List<String> waldLine = new ArrayList<>(line);
			
			double willsVal = willsVals.get(i);
			double waldVal = waldVals.get(i);
			if (!wills.isValueValid(willsVal)) {
				System.err.println("Warning, Wills="+willsVal+", defaulting to Wald="+waldVal
						+" for "+locs.get(i));
				willsVal = waldVal;
			}
			willsLine.set(vsCol, (float)willsVal+"");
			waldLine.set(vsCol, (float)waldVal+"");
			
			if (delCol >= 0) {
				willsLine.remove(delCol);
				waldLine.remove(delCol);
			}
			
			willsCSV.addLine(willsLine);
			waldCSV.addLine(waldLine);
		}
		
		willsCSV.writeToFile(new File(dir, inputFile.getName().replaceAll(".csv", "-wills2015.csv")));
		waldCSV.writeToFile(new File(dir, inputFile.getName().replaceAll(".csv", "-wald2007.csv")));
	}

}
