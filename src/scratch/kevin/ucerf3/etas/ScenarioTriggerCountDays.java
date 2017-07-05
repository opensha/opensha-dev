package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ScenarioTriggerCountDays {

	public static void main(String[] args) throws IOException {
		int startYear = 2012;
		long ot = Math.round((startYear-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		double minMag = 7;
		int maxDays = 3;
		
		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_02_24-bombay_beach_m4pt8-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont");
		
		System.out.println(dir.getName());
		doIt(new File(dir, "results_m4.bin"), ot, minMag, maxDays);
		
		
		dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_02_24-bombay_beach_m4pt8-10yr-no_ert-subSeisSupraNucl-gridSeisCorr-noSpont");
		
		System.out.println();
		System.out.println(dir.getName());
		doIt(new File(dir, "results_m4.bin"), ot, minMag, maxDays);
	}
	
	private static void doIt(File binFile, long ot, double minMag, int maxDays) throws IOException {
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(binFile);
		int triggerParentID = 0;
		List<List<ETAS_EqkRupture>> childrenCatalogs = Lists.newArrayList();
		for (List<ETAS_EqkRupture> catalog : catalogs)
			childrenCatalogs.add(ETAS_SimAnalysisTools.getChildrenFromCatalog(catalog, triggerParentID));
		
		ETAS_MultiSimAnalysisTools.calcNumWithMagAbove(childrenCatalogs, ot, minMag, triggerParentID, maxDays);
	}

}
