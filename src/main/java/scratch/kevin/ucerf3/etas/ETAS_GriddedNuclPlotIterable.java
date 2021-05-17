package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.util.ClassUtils;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.BinarayCatalogsIterable;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;

public class ETAS_GriddedNuclPlotIterable {

	public static void main(String[] args) throws IOException, GMT_MapException {
		if (args.length < 5 || args.length > 6) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ETAS_GriddedNuclPlotIterable.class)
					+" <bin-file> <output-dir> <duration> <name> <prefix> [<maxOT>]");
			System.exit(2);
		}
		File catFile = new File(args[0]);
		File outputDir = new File(args[1]);
		
		BinarayCatalogsIterable catalogs = ETAS_CatalogIO.getBinaryCatalogsIterable(catFile, 0d);
		int numCatalogs = catalogs.getNumCatalogs();
		
		double duration = Double.parseDouble(args[2]);
		String name = args[3];
		String prefix = args[4];
		
		double[] mags = { 2.5 };
		
		long maxOT = Long.MIN_VALUE;
		if (args.length == 6)
			maxOT = Long.parseLong(args[5]);
		
		ETAS_MultiSimAnalysisTools.plotCubeNucleationRates(catalogs, numCatalogs,
				duration, maxOT, outputDir, name, prefix, mags, true);
	}

}
