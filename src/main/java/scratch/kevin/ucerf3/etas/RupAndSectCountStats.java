package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.BitSet;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;

public class RupAndSectCountStats {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2026_05_27-FSS_Rupture_201887_M7p8_Start_2026_10_15_1_yr_kCOV_1p5_MaxPtSrcM_6");
		File configFile = new File(dir, "config.json");
		File binFile = new File(dir, "results_m5_preserve_chain.bin");
		
		int targetRupID = 218331;
		
		String[] targetSects = {
				"Raymond",
				"Hollywood"
		};
		
		ETAS_Config config = ETAS_Config.readJSON(configFile);
		ETAS_Launcher launcher = new ETAS_Launcher(config, false);
		FaultSystemSolution sol = launcher.checkOutFSS();
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		int[] parentIDs = new int[targetSects.length];
		BitSet[] targetSectRups = new BitSet[targetSects.length];
		BitSet targetCorups = null;
		for (int i=0; i<targetSects.length; i++) {
			int parentID = FaultSectionUtils.findParentSectionID(sol.getRupSet().getFaultSectionDataList(), targetSects[i]);
			Preconditions.checkState(parentID >= 0);
			parentIDs[i] = parentID;
			targetSectRups[i] = new BitSet(rupSet.getNumRuptures());
			for (int rupIndex : rupSet.getRupturesForParentSection(parentID))
				targetSectRups[i].set(rupIndex);
			if (targetCorups == null)
				targetCorups = (BitSet)targetSectRups[i].clone();
			else
				targetCorups.and(targetSectRups[i]);
		}
		
		int numWithRup = 0;
		int numWithCorupture = 0;
		int numWithAny = 0;
		int[] numWithEach = new int[targetSects.length];
		int numCatalogs = 0;
		
		for (ETAS_Catalog catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(binFile, 0d)) {
			boolean hasRup = false;
			boolean hasCorup = false;
			boolean hasAny = false;
			boolean[] hasIndv = new boolean[targetSects.length];
			for (ETAS_EqkRupture rup : catalog) {
				int fssIndex = rup.getFSSIndex();
				if (fssIndex < 0)
					continue;
				if (fssIndex == targetRupID)
					hasRup = true;
				if (targetCorups.get(fssIndex))
					hasCorup = true;
				for (int i=0; i<targetSects.length; i++) {
					if (targetSectRups[i].get(fssIndex)) {
						hasAny = true;
						hasIndv[i] = true;
					}
				}
			}
			if (hasRup)
				numWithRup++;
			if (hasCorup)
				numWithCorupture++;
			if (hasAny)
				numWithAny++;
			for (int i=0; i<targetSects.length; i++)
				if (hasIndv[i])
					numWithEach[i]++;
			numCatalogs++;
		}
		
		DecimalFormat pDF = new DecimalFormat("0.0000%");
		
		System.out.println("Processed "+numCatalogs+" catalogs");
		System.out.println(numWithRup+"/"+numCatalogs+" ("
				+pDF.format((double)numWithRup/(double)numCatalogs)+") have rupture "+targetRupID);
		System.out.println(numWithAny+"/"+numCatalogs+" ("
				+pDF.format((double)numWithAny/(double)numCatalogs)
				+") have rupture on "+Joiner.on(" or ").join(targetSects));
		System.out.println(numWithCorupture+"/"+numCatalogs+" ("
				+pDF.format((double)numWithCorupture/(double)numCatalogs)
				+") have corupture of "+Joiner.on(" and ").join(targetSects));
		for (int i=0; i<targetSects.length; i++)
			System.out.println(numWithEach[i]+"/"+numCatalogs+" ("
					+pDF.format((double)numWithEach[i]/(double)numCatalogs)+") have a rupture on "+targetSects[i]);
	}

}
