package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.util.FileUtils;

public class SandarshStationListCreate {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		CSVFile<String> file = CSVFile.readFile(new File("/tmp/CS_Gridded_Sites.txt"), true, 3);
//		
//		LocationList locs = new LocationList();
//		
//		WillsMap2006 wills = new WillsMap2006();
//		
//		for (int i=0; i<file.getNumRows(); i++) {
//			List<String> line = file.getLine(i);
//			double lat = Double.parseDouble(line.get(1));
//			double lon = Double.parseDouble(line.get(2));
//			
//			locs.add(new Location(lat, lon));
//		}
//		
//		ArrayList<Double> vals = wills.getValues(locs);
//		
//		CSVFile<String> out = new CSVFile<String>(true);
//		
//		for (int i=0; i<file.getNumRows(); i++) {
//			List<String> line = file.getLine(i);
//			line.add(vals.get(i)+"");
//			out.addLine(line);
//		}
//		
//		out.writeToFile(new File("/tmp/CS_Gridded_Sites_withVs30.txt"));
		
		LocationList locs = new LocationList();
		ArrayList<String[]> lines = new ArrayList<String[]>();
		
		for (String line : FileUtils.loadFile("/tmp/hw_stats.txt")) {
			if (line.length() < 2)
				continue;
			
			StringTokenizer tok = new StringTokenizer(line);
			
			String[] srs = new String[tok.countTokens()];
			for (int i=0; i<srs.length; i++)
				srs[i] = tok.nextToken();
			
			lines.add(srs);
			
			double lon = Double.parseDouble(srs[0]);
			double lat = Double.parseDouble(srs[1]); 
			
			locs.add(new Location(lat, lon));
		}
		
		WillsMap2006 wills = new WillsMap2006();
		ArrayList<Double> vals = wills.getValues(locs);
		
		FileWriter out = new FileWriter(new File("/tmp/hw_stats.new.txt"));
		
		for (int i=0; i<vals.size(); i++) {
			String[] line = lines.get(i);
			for (int j=0; j<3; j++) {
				if (j>0)
					out.write("\t");
				out.write(line[j]);
			}
			double val = vals.get(i);
			if (Double.isNaN(val))
				val = 180d;
			out.write("\t"+val);
			for (int j=3; j<line.length; j++) {
				out.write("\t"+line[j]);
			}
			out.write("\n");
		}
		
		out.close();
	}

}
