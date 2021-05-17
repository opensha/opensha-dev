package scratch.ned.GK_Declustering;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.TimeZone;

import org.jfree.data.Range;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.IncrementalMagFreqDist;


/**
 * This reads a USGS NSHMP catalog (at least as they were defined in 2020).
 * 
 * This defines event IDs as the nth event in the catalog (starting at 0)
 * @author field
 *
 */
public class Hauksson2012_CatalogParser {
	
	public static ObsEqkRupList loadCatalog(File file) throws IOException {
		ObsEqkRupList rups = new ObsEqkRupList();
		
		BufferedReader in = new BufferedReader(new FileReader(file));
		
		// used if 10 column input without IDs
		int startID = 1;
		
		String line;
		String[] split;
		int eventID = -1;
		
		while (in.ready()) {
			line = in.readLine();
			line = line.trim();
			while (line.contains("  "))
				line = line.replaceAll("  ", " ");
			split = line.split(" ");

			int year			= Integer.parseInt(split[0]);
			int month			= Integer.parseInt(split[1]);
			int date			= Integer.parseInt(split[2]);
			int hourOfDay		= Integer.parseInt(split[3]);
			int minute			= Integer.parseInt(split[4]);
			double second			= Double.parseDouble(split[5]);
			// sixth element is SCSN cuspid (up to 9 digits)
			double latitude		= Double.parseDouble(split[7]);
			double longitude	= Double.parseDouble(split[8]);
			double depth		= Double.parseDouble(split[9]);
			double mag			= Double.parseDouble(split[10]);
			
			eventID += 1;
						
			GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("GMT-0:00"));
			cal.clear();
			cal.set(year, month-1, date, hourOfDay, minute, (int)second);
			
			Location hypoLoc = new Location(latitude, longitude, depth);
						
			rups.add(new ObsEqkRupture(eventID+"", cal.getTimeInMillis(), hypoLoc, mag));
		}
		
		in.close();
				
		return rups;
	}
	

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File file = new File("/Users/field/MiscDocs/MuellerGK_FortranCode/wmm.c2");
		ObsEqkRupList rupList = loadCatalog(file);
		for (ObsEqkRupture rup : rupList) {
//			System.out.println(rup);
		}
		System.out.println("numEvents = "+rupList.size());
	}
	
	

}
