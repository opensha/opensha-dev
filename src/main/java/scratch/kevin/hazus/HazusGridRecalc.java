package scratch.kevin.hazus;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.sha.calc.hazus.parallel.HazusJobWriter;

public class HazusGridRecalc {
	
	private static void writeCSV(LocationList locs, File file) throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		for (Location loc : locs) {
			ArrayList<String> line = new ArrayList<String>();
			line.add(loc.getLatitude()+"");
			line.add(loc.getLongitude()+"");
			csv.addLine(line);
		}
		csv.writeToFile(file);
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File origFile = new File("/home/kevin/OpenSHA/hazus/05grid.csv");
		File newFile = new File("/home/kevin/OpenSHA/hazus/05GridV2.csv");
		
		LocationList origLocs = HazusJobWriter.loadCSV(origFile);
		LocationList newLocs = HazusJobWriter.loadCSV(newFile);
		
		LocationList uniqueNew = new LocationList();
		LocationList uniqueOld = new LocationList();
		
		for (Location newLoc : newLocs) {
			if (origLocs.contains(newLoc))
				continue;
			uniqueNew.add(newLoc);
		}
		for (Location origLoc : origLocs) {
			if (newLocs.contains(origLoc))
				continue;
			uniqueOld.add(origLoc);
		}
		
		writeCSV(uniqueNew, new File("/home/kevin/OpenSHA/hazus/grid_fix/new.csv"));
		writeCSV(uniqueOld, new File("/home/kevin/OpenSHA/hazus/grid_fix/old.csv"));
		
		System.out.println("Unique new: "+uniqueNew.size());
		System.out.println("Unique old: "+uniqueOld.size());
	}

}
