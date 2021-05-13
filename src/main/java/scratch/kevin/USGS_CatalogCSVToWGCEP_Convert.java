package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.opensha.commons.data.CSVFile;

public class USGS_CatalogCSVToWGCEP_Convert {

	public static void main(String[] args) throws IOException, ParseException {
		File csvFile = new File("/tmp/2.5_week.csv");
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		FileWriter fw = new FileWriter(new File("/tmp/2.5_week_WGCEPformat.txt"));
		
		for (int row=csv.getNumRows(); --row>=1;) {
			Date date = inputDF.parse(csv.get(row, 0));
			double lat = Double.parseDouble(csv.get(row, 1));
			double lon = Double.parseDouble(csv.get(row, 2));
			double depth = Double.parseDouble(csv.get(row, 3));
			double mag = Double.parseDouble(csv.get(row, 4));
			
			String line = outputDF.format(date)+" "+lat+" "+lon+" "+depth+" "+mag+" -1 -1 -1";
			System.out.println(line);
			fw.write(line+"\n");
		}
		fw.close();
	}
	
	private static final SimpleDateFormat inputDF = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'");
	private static final SimpleDateFormat outputDF = new SimpleDateFormat("yyyy MM dd HH mm ss.SSS");

}
