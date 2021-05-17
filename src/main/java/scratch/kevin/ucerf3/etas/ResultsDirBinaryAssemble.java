package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;

public class ResultsDirBinaryAssemble {

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.err.println("USAGE: <results-dir> <output-file> [<min-mag>]");
			System.exit(0);
		}
		File resultsDir = new File(args[0]);
		Preconditions.checkState(resultsDir.exists() && resultsDir.isDirectory());
		File outputFile = new File(args[1]);
		double minMag = args.length > 2 ? Double.parseDouble(args[2]) : 0d;
		
		List<ETAS_Catalog> catalogs = new ArrayList<>();
		for (File subDir : resultsDir.listFiles()) {
			File binFile = new File(subDir, "simulatedEvents.bin");
			File asciiFile = new File(subDir, "simulatedEvents.txt");
			if (binFile.exists()) {
				System.out.println("Loading "+binFile.getAbsolutePath());
				ETAS_Catalog catalog = ETAS_CatalogIO.loadCatalogBinary(binFile, minMag);
				System.out.println("\tLoaded "+catalog.size()+" events");
				catalogs.add(catalog);
			} else if (asciiFile.exists()) {
				if (ETAS_Launcher.isAlreadyDone(subDir)) {
					System.out.println("Loading "+asciiFile.getAbsolutePath());
					ETAS_Catalog catalog = ETAS_CatalogIO.loadCatalog(asciiFile, minMag);
					System.out.println("\tLoaded "+catalog.size()+" events");
					catalogs.add(catalog);
				} else {
					System.out.println("Skipping "+subDir.getName()+" (in progress...)");
				}
			}
		}
		System.out.println("Loaded "+catalogs.size()+" catalogs, writing to "+outputFile.getAbsolutePath());
		ETAS_CatalogIO.writeCatalogsBinary(outputFile, catalogs);
		System.out.println("DONE");
		System.exit(0);
	}

}
