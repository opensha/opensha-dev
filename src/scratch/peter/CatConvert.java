package scratch.peter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.TimeZone;

import org.apache.commons.io.IOUtils;
import org.opensha.commons.eq.cat.CatTools;
import org.opensha.commons.eq.cat.Catalog;
import org.opensha.commons.eq.cat.MutableCatalog;
import org.opensha.commons.eq.cat.filters.ExtentsFilter;
import org.opensha.commons.eq.cat.io.Reader_SCEDC;

public class CatConvert {
	
	private static Calendar cal = new GregorianCalendar(TimeZone.getTimeZone("PDT"));
	//
	// fjk

	//
	private String test;
	
	// El-Mayor-Cucapah.cat
	//		start: 01-01-2010
	//		  end: 12-09-2010

	//
	public static void main(String[] args) {

		Catalog catOut = null;
		try {
			File in = new File("tmp/cat_gen/El-Mayor-Cucapah.cat");
			File out = new File("tmp/cat_gen/El-Mayor-Cucapah_Mgt3.js");
			MutableCatalog cat = new MutableCatalog(in, new Reader_SCEDC(20000));
			System.out.println(" Source Size: " + cat.size());
			cal.set(2010,0,1,0,0,0);
			long start = cal.getTimeInMillis();
			cal.set(2010,11,31,0,0,0);
			long end = cal.getTimeInMillis();
			ExtentsFilter filter = new ExtentsFilter();
			filter.setDates(start, end);
			filter.setMagnitudes(3, 8);
			int[] idx = filter.process(cat);
			catOut = cat.deriveCatalog(idx, new MutableCatalog());
			System.out.println("Trimmed Size: " + catOut.size());
			BufferedWriter br = new BufferedWriter(new FileWriter(out));
			br.write("var cat='");
			br.write(CatTools.toJSON(catOut));
			br.write("';");
			IOUtils.closeQuietly(br);
			//System.out.println(CatTools.toJSON(catOut));
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	 
	
	
}
