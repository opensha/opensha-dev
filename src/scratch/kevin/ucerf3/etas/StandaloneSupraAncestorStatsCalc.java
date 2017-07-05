package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.util.List;

import org.opensha.commons.util.ClassUtils;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;

public class StandaloneSupraAncestorStatsCalc {

	public static void main(String[] args) {
		if (args.length < 1 || args.length > 2) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(StandaloneSupraAncestorStatsCalc.class)
				+" <cat-bin-file> [<output-dir>]");
			System.exit(2);
		}
		
		File catalogsFile = new File(args[0]);
		File outputDir = null;
		if (args.length == 2)
			outputDir = new File(args[1]);		
		
		try {
			Iterable<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.getBinaryCatalogsIterable(catalogsFile, 0d);
			
			ETAS_MultiSimAnalysisTools.calcSupraAncestorStats(catalogs, outputDir);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}

}
