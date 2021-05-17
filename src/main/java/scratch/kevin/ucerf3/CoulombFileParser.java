package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.StringTokenizer;
import java.util.TimeZone;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.SimpleFaultData;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;

import com.google.common.base.Preconditions;

public class CoulombFileParser {
	
	/**
	 * Loads a the finite fault surface from a Coulomb input file.
	 * 
	 * @param file
	 * @param gridSpacing
	 * @return
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static RuptureSurface loadFaultSurface(File file, double gridSpacing) throws FileNotFoundException, IOException {
		ArrayList<String> lines = FileUtils.loadFile(file.getAbsolutePath());
		
		Preconditions.checkArgument(gridSpacing > 0, "gridSpacing must be positive");
		
		double zeroLon = Double.NaN;
		double zeroLat = Double.NaN;
		
		// first go through and find the zero points
		for (String line : lines) {
			if (line.contains("zero lon")) {
				String[] split = line.trim().split(" ");
				zeroLon = Double.parseDouble(split[split.length-1]);
			} else if (line.contains("zero lat")) {
				String[] split = line.trim().split(" ");
				zeroLat = Double.parseDouble(split[split.length-1]);
			}
		}
		
		Preconditions.checkState(!Double.isNaN(zeroLat), "zero lat field not found in file!");
		Preconditions.checkState(!Double.isNaN(zeroLon), "zero lon field not found in file!");
		
		Location origin = new Location(zeroLat, zeroLon);
		
		ArrayList<EvenlyGriddedSurface> surfs = new ArrayList<EvenlyGriddedSurface>();
		
		// now read the fault surfaces
		int cnt = 0;
		boolean reading = false;
		for (String line : lines) {
			if (line.contains("X-start") && line.contains("Y-start")) {
				reading = true;
				continue;
			}
			if (reading) {
				line = line.trim();
				if (line.startsWith("xxx")) {
					// skip
					continue;
				}
				if (line.isEmpty() || line.contains("Grid Parameters")) {
					reading = false;
					continue;
				}
				StringTokenizer tok = new StringTokenizer(line);
				tok.nextToken(); // index
				
				// offsets, all in KM
				double xStart =	Double.parseDouble(tok.nextToken());
				double yStart =	Double.parseDouble(tok.nextToken());
				double xFin =	Double.parseDouble(tok.nextToken());
				double yFin =	Double.parseDouble(tok.nextToken());
				tok.nextToken(); // Kode
				tok.nextToken(); // rt.lat
				tok.nextToken(); // reverse
				double dip =	Double.parseDouble(tok.nextToken());
				double top =	Double.parseDouble(tok.nextToken());
				double bottom =	Double.parseDouble(tok.nextToken());
				
				Location start = getRelativeLocation(origin, xStart, yStart);
				start = new Location(start.getLatitude(), start.getLongitude(), top);
//				System.out.println("Start: "+start);
				Location end = getRelativeLocation(origin, xFin, yFin);
				end = new Location(end.getLatitude(), end.getLongitude(), top);
//				System.out.println("end: "+end);
				
				FaultTrace trace = new FaultTrace("Coulomb Fault "+(cnt++)+" ("+file.getName()+")");
				trace.add(start);
				trace.add(end);
				
				SimpleFaultData sfd = new SimpleFaultData(dip, bottom, top, trace);
				
				surfs.add(new StirlingGriddedSurface(sfd, gridSpacing));
			}
		}
		
		Preconditions.checkState(!surfs.isEmpty(), "no fault surfaces found in file!");
		
		if (surfs.size() == 1)
			return surfs.get(0);
		
		return new CompoundSurface(surfs);
	}
	
	private static Location getRelativeLocation(Location origin, double xOffset, double yOffset) {
		// first move vertically/North (or South if negative) (azimuth = 0)
		Location loc = LocationUtils.location(origin, 0, yOffset);
		// then move horizontally/East (or West if negative) (azimuth = PI/2)
		loc = LocationUtils.location(loc, Math.PI/2d, xOffset);
		
		return loc;
	}
	
	private static ObsEqkRupture loadObsEqkRup(File file, Location hypocenter,
			long date, double mag) throws FileNotFoundException, IOException {
		
		System.out.println("Loading rupture: "+file.getName());
		RuptureSurface surf = loadFaultSurface(file, 1d);
		
		ObsEqkRupture rup = new ObsEqkRupture(file.getName(), date, hypocenter, mag);
		rup.setRuptureSurface(surf);
		
		return rup;
	}
	
	private static long getDate(int year, int month, int day, int hours, int mins, double secs) {
		GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("UTC"));
		cal.set(year, month-1, day, hours, mins, (int)(secs+0.5));
		System.out.println("Parsed date: "+new SimpleDateFormat().format(cal.getTime()));
		return cal.getTimeInMillis();
	}
	
	public static ObsEqkRupList loadObsEqkRups(File dir) throws IOException {
		ObsEqkRupList rups = new ObsEqkRupList();
		
		// hardcoded filenames since the files themselves don't have the needed info
		
		// missing
//		rups.add(loadObsEqkRup(
//				new File(dir, "1976.11.26, Mw=6.7, W of Trinidad, CA.inp"),
//				)
		
		rups.add(loadObsEqkRup(
				new File(dir, "1980.11.08, Mw=7.3, W of Eureka, CA.inp"),
				new Location(41.084, -124.616), getDate(1980, 11, 8, 10, 27, 33.2), 7.3));
		
		// no match (closest is mag 5.5
//		rups.add(loadObsEqkRup(
//				new File(dir, "1983.08.24, Mw=6.3, off Cape Mendocino, CA.inp"),
//				);
		
		rups.add(loadObsEqkRup(
				new File(dir, "1984.09.10, Mw=6.6, Mendocino Fault Zone.inp"),
				new Location(40.504, -125.130), getDate(1984, 9, 10, 3, 14, 28.1), 6.6));
		
		rups.add(loadObsEqkRup(
				new File(dir, "1987.07.31, Mw=6.0, off Cape Mendocino, CA.inp"),
				new Location(40.416, -124.383), getDate(1987, 7, 31, 23, 56, 57.830), 6.0));
		
		// no match
//		rups.add(loadObsEqkRup(
//				new File(dir, "1991.07.13, Mw=6.8, bent model.inp"),
//				);
		
		// no match
//		rups.add(loadObsEqkRup(
//				new File(dir, "1991.08.16, Mw=6.3, off Crescent City, CA.inp"),
//				);
		
		rups.add(loadObsEqkRup(
				new File(dir, "1991.08.17, Mw=6.1, Honeydew, CA.inp"),
				new Location(40.252, -124.286), getDate(1991, 8, 17, 19, 29, 40.00), 6.1));
		
		// no match
//		rups.add(loadObsEqkRup(Mw=7.1, off Trinidad, CA.inp
//				new File(dir, "1991.08.17, Mw=6.1, Honeydew, CA.inp"),
//				);
		
		rups.add(loadObsEqkRup(
				new File(dir, "1992.04.25, Mw=6.9, Cape Mendocino, CA, Oppenheimer et al 1993 displacement model.inp"),
				new Location(40.335, -124.229), getDate(1992, 4, 25, 18, 6, 5.18), 6.9));
		
		rups.add(loadObsEqkRup(
				new File(dir, "1992.04.26, Mw=6.5, off Cape Mendocino, CA.inp"),
				new Location(40.432, -124.566), getDate(1992, 4, 26, 7, 41, 40.09), 6.5));
		
		rups.add(loadObsEqkRup(
				new File(dir, "1992.04.26, Mw=6.6, off Cape Mendocino, CA.inp"),
				new Location(40.383, -124.555), getDate(1992, 4, 26, 11, 18, 25.98), 6.6));
		
		// no match, and duplicate?
//		1994.09.01, Mw=7.0, Mendocino Fracture Zone, bilateral rupture.inp
//		1994.09.01, Mw=7.0, Mendocino Fracture Zone, unilateral eastward rupture.inp
		
		// no match
//		1995.02.19, Mw=6.6, off Cape Mendocino, CA.inp
		
		rups.add(loadObsEqkRup(
				new File(dir, "2000.03.16, Mw=5.9, Mendocino Fault Zone.inp"),
				new Location(40.389, -125.239), getDate(2000, 3, 16, 15, 19, 56.38), 5.9));
		
		// no match
//		2005.06.15, Mw=7.2, off Eureka, CA (Shao and Ji finite fault model).inp
//		2005.06.17, Mw=6.6, off Cape Mendocino, CA.inp
//		2008.11.28, Mw=5.9, Mendocino Fault Zone, CA.inp
		
		rups.add(loadObsEqkRup(
				new File(dir, "2010.01.10, Mw=6.5, off Ferndale, CA (Dreger et al variable slip model with UCB epicenter).inp"),
				new Location(40.652, -124.692), getDate(2010, 1, 10, 0, 27, 39.320), 6.5));
		
		
		
		return rups;
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/rollins_stein_finite");
		
		loadObsEqkRups(dir);
		System.exit(0);
		
		File outputDir = new File(dir, "gridded");
		outputDir.mkdir();
		
		double gridSpacing = 1d;
		
		int cnt = 0;
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				continue;
			if (!file.getName().endsWith(".inp"))
				continue;
			
			System.out.println("Loading: "+file.getName());
			
			RuptureSurface surf = loadFaultSurface(file, gridSpacing);
			
			File outputFile = new File("/tmp/surf"+(cnt++)+".txt");
			System.out.println("Writing: "+outputFile.getAbsolutePath());
			CatalogWriter.writeFiniteSurfaceFile(surf, outputFile, 0.5*gridSpacing);
		}
	}

}
