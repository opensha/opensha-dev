package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.zip.ZipException;

import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ETAS_RandomizeBinary {

	public static void main(String[] args) throws ZipException, IOException {
		if (args.length != 4) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ETAS_RandomizeBinary.class)
						+" <input> <output> <startYear> <durationYears>");
			System.exit(2);
		}
		
		File inputFile = new File(args[0]);
		File outputFile = new File(args[1]);
		double startYear = Double.parseDouble(args[2]);
		double duration = Double.parseDouble(args[3]);
		
		long ot = Math.round((startYear-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		double durationMillis = duration * ProbabilityModelsCalc.MILLISEC_PER_YEAR;
		
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(inputFile);
		
		CatalogComparator comp = new CatalogComparator();
		
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			for (ETAS_EqkRupture rup : catalog) {
				long time = ot + (long)(Math.random()*durationMillis);
				rup.setOriginTime(time);
			}
			Collections.sort(catalog, comp);
		}
		
		ETAS_CatalogIO.writeCatalogsBinary(outputFile, catalogs);
	}
	
	private static class CatalogComparator implements Comparator<ETAS_EqkRupture> {

		@Override
		public int compare(ETAS_EqkRupture o1, ETAS_EqkRupture o2) {
			return new Long(o1.getOriginTime()).compareTo(o2.getOriginTime());
		}
		
	}

}
