package scratch.kevin.ucerf3;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.FileUtils;

public class FaultRuptureMatch {
	
	private static class Quake {
		private int year;
		private int month;
		private int day;
		private double mag;
		private double lat;
		private double lon;
		
		private String name;
		private int id;
		
		public Quake(int year, int month, int day, double mag, double lat,
				double lon) {
			super();
			this.year = year;
			this.month = month;
			this.day = day;
			this.mag = mag;
			this.lat = lat;
			this.lon = lon;
		}
		
		public Quake(int year, int month, int day, double mag, double lat,
				double lon, String name, int id) {
			this(year, month, day, mag, lat, lon);
			this.name = name;
			this.id = id;
		}

		@Override
		public String toString() {
			String ret = "Quake [year=" + year + ", month=" + month + ", day=" + day
					+ ", mag=" + mag + ", lat=" + lat + ", lon=" + lon;
			if (name != null)
				ret += ", name=" + name+ ", id=" + id;
			ret += "]";
			
			return ret;
		}
	}
	
	private static List<Quake> loadTXT(File file) throws FileNotFoundException, IOException {
		ArrayList<Quake> quakes = new ArrayList<Quake>();
		for (String line : FileUtils.loadFile(file.getAbsolutePath())) {
			line = line.trim();
			while (line.contains("  "))
				line = line.replaceAll("  ", " ");
			String[] strs = line.split(" ");
			
//			System.out.println(line);
			
			int year = Integer.parseInt(strs[0]);
			int month = Integer.parseInt(strs[1]);
			int day = Integer.parseInt(strs[2]);
			
			double lat = Double.parseDouble(strs[6]);
			double lon = Double.parseDouble(strs[7]);
			double mag = Double.parseDouble(strs[9]);
			
			quakes.add(new Quake(year, month, day, mag, lat, lon));
		}
		
		return quakes;
	}
	
	private static List<Quake> loadExcel(File file) throws FileNotFoundException, IOException {
		ArrayList<Quake> quakes = new ArrayList<Quake>();
		
		POIFSFileSystem fs = new POIFSFileSystem(new BufferedInputStream(new FileInputStream(file)));
		HSSFWorkbook wb = new HSSFWorkbook(fs);
		HSSFSheet sheet = wb.getSheetAt(0);
		
		for (int i=1; i<=sheet.getLastRowNum(); i++) {
			HSSFRow row = sheet.getRow(i);
			int id = (int)row.getCell(0).getNumericCellValue();
			String name = row.getCell(1).getStringCellValue().trim();
			int year;
			try {
				year = (int)row.getCell(2).getNumericCellValue();
			} catch (Exception e) {
				continue;
			}
			String monDayStr;
			try {
				monDayStr = row.getCell(3).getStringCellValue();
			} catch (Exception e) {
				monDayStr = (int)row.getCell(3).getNumericCellValue() + "";
			}
			if (monDayStr == null || monDayStr.length()<3)
				continue;
			int month;
			int day;
			if (monDayStr.length() == 3) {
				month = Integer.parseInt(monDayStr.substring(0, 1));
				day = Integer.parseInt(monDayStr.substring(1));
			} else {
				month = Integer.parseInt(monDayStr.substring(0, 2));
				day = Integer.parseInt(monDayStr.substring(2));
			}
			
			double mag = row.getCell(5).getNumericCellValue();
			
			double lat;
			double lon;
			try {
				lat = row.getCell(20).getNumericCellValue();
				lon = row.getCell(21).getNumericCellValue();
			} catch (Exception e) {
				continue;
			}
			
			Quake quake = new Quake(year, month, day, mag, lat, lon, name, id);
//			System.out.println(quake);
			quakes.add(quake);
		}
		
		return quakes;
	}
	
	public static void printMatches(List<Quake> excel, List<Quake> txt) {
		ArrayList<Quake> noMatches = new ArrayList<FaultRuptureMatch.Quake>();
		
		int numMatches = 0;
		for (Quake txtQuake : txt) {
			Location txtLoc = new Location(txtQuake.lat, txtQuake.lon);
			ArrayList<Quake> matches = new ArrayList<FaultRuptureMatch.Quake>();
			ArrayList<Double> distances = new ArrayList<Double>();
			for (Quake excelQuake : excel) {
				if (txtQuake.year != excelQuake.year)
					continue;
				if (txtQuake.month != excelQuake.month)
					continue;
				if (txtQuake.day != excelQuake.day)
					continue;
				
				if (Math.abs(txtQuake.mag - excelQuake.mag) > 0.5)
					continue;
				
				Location excelLoc = new Location(excelQuake.lat, excelQuake.lon);
				double dist = LocationUtils.horzDistance(txtLoc, excelLoc);
				if (dist > 50)
					continue;
				
				matches.add(excelQuake);
				distances.add(dist);
			}
//			
			if (matches.isEmpty()) {
//				System.out.println("\t(none)");
				noMatches.add(txtQuake);
			} else {
				System.out.println("\n"+txtQuake);
				numMatches++;
				for (int i=0; i<matches.size(); i++) {
					Quake quake = matches.get(i);
					double dist = distances.get(i);
					System.out.println("\t* ID: "+quake.id+", dist: "+(float)dist+" (KM), name: "+quake.name+", "+quake);
				}
			}
		}
		System.out.println("Num matches: "+numMatches+"/"+txt.size());
		System.out.println();
		System.out.println("Non Matched Quakes:");
		for (Quake quake : noMatches)
			System.out.println(quake);
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		File dir = new File("D:\\Documents\\temp");
		File txtFile = new File(dir, "large_EQs-1.txt");
		File exelFile = new File(dir, "EQ.V8.xls");
		
		System.out.println("----- TXT Quakes -----");
		List<Quake> txtQuakes = loadTXT(txtFile);
		
		System.out.println("\n----- Excel Quakes -----");
		List<Quake> excelQuakes = loadExcel(exelFile);
		
		printMatches(excelQuakes, txtQuakes);
	}

}

