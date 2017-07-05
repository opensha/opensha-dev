package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.sql.Date;
import java.util.List;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;

public class FirstRuptureFinder {

	public static void main(String[] args) {
		File intputFile = new File(args[0]);
		
		long minOrigTime = Long.MAX_VALUE;
		ETAS_EqkRupture firstRup = null;
		
		long maxOrigTime = Long.MIN_VALUE;
		ETAS_EqkRupture lastRup = null;
		
		for (List<ETAS_EqkRupture> catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(intputFile, -1)) {
			ETAS_EqkRupture rup = catalog.get(0);
			long orig = rup.getOriginTime();
			if (orig < minOrigTime) {
				minOrigTime = orig;
				firstRup = rup;
			}
			rup = catalog.get(catalog.size()-1);
			if (orig > maxOrigTime) {
				maxOrigTime = orig;
				lastRup = rup;
			}
		}
		
		System.out.println("First rupture occured at "+minOrigTime);
		System.out.println("Date object: "+new Date(minOrigTime));
		System.out.println("ETAS Catalog:");
		System.out.println("\t"+ETAS_CatalogIO.EVENT_FILE_HEADER);
		System.out.println("\t"+ETAS_CatalogIO.getEventFileLine(firstRup));
		System.out.println();
		System.out.println("Last rupture occured at "+maxOrigTime);
		System.out.println("Date object: "+new Date(maxOrigTime));
		System.out.println("ETAS Catalog:");
		System.out.println("\t"+ETAS_CatalogIO.EVENT_FILE_HEADER);
		System.out.println("\t"+ETAS_CatalogIO.getEventFileLine(lastRup));
	}

}
