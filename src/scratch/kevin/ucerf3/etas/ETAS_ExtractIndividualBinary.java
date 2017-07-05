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
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.BinarayCatalogsIterable;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ETAS_ExtractIndividualBinary {

	public static void main(String[] args) throws ZipException, IOException {
		if (args.length < 3 || args.length > 4) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ETAS_ExtractIndividualBinary.class)
						+" <input> <output-dir> <num> [<startIndex>]");
			System.exit(2);
		}
		
		File inputFile = new File(args[0]);
		File outputDir = new File(args[1]);
		Preconditions.checkState((outputDir.exists() && outputDir.isDirectory()) || outputDir.mkdir(),
				"Output directory doesn't exist and can't be created, or exists and isn't a directory!");
		int num = Integer.parseInt(args[2]);
		int startIndex = 0;
		if (args.length == 4)
			startIndex = Integer.parseInt(args[3]);
		
		BinarayCatalogsIterable it = ETAS_CatalogIO.getBinaryCatalogsIterable(inputFile, 0d);
		
		int digits = ((it.getNumCatalogs()-1)+"").length();
		
		int index = 0;
		int count = 0;
		for (List<ETAS_EqkRupture> catalog : it) {
			if (index < startIndex) {
				index++;
				continue;
			}
			String countStr = count+"";
			while (countStr.length() < digits)
				countStr = "0"+countStr;
			File outputFile = new File(outputDir, "catalog_"+countStr+".bin");
			System.out.println("Writing "+outputFile.getAbsolutePath());
			ETAS_CatalogIO.writeCatalogBinary(outputFile, catalog);
			index++;
			count++;
			if (count == num)
				break;
		}
	}

}
