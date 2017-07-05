package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import scratch.UCERF3.inversion.BatchPlotGen;

public class MisfitComparison {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		Map<String, Double> shortMisfits = BatchPlotGen.loadMisfitsFile(
				new File("/tmp/FM2_1_UC2ALL_Shaw09Mod_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU2_VarPaleo1_VarSectNuclMFDWt0.01_sol.zip.misfits"));
		Map<String, Double> longMisfits = BatchPlotGen.loadMisfitsFile(
				new File("/tmp/FM2_1_UC2ALL_Shaw09Mod_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU2_VarPaleo1_run00_sol.zip.misfits"));
		
		for (String key : shortMisfits.keySet()) {
			Double shortMisfit = Math.sqrt(shortMisfits.get(key));
			Double longMisfit = Math.sqrt(longMisfits.get(key));
			double p_imp = (shortMisfit - longMisfit) / shortMisfit;
			p_imp *= 100;
			System.out.println(key+":\t "+shortMisfit+" => "+longMisfit+"\t("+(float)p_imp+" % improvement)");
		}
		
//		double[] runLengths = { 3, 4, 6, 500d/60d };
		double[] runLengths = { 3, 4, 6, 8 };
		double availableHours = 6*24 + 8; // midnight tues night through 8am next tues morn
		int availableNodes = 60;
		
		for (double runLength : runLengths) {
			double padRunLength = runLength + 0.5; // for processing time
			int runs = 0;
			for (double hours=0; hours+padRunLength<availableHours; hours+=padRunLength)
				runs += availableNodes;
			System.out.println((int)runLength+" hrs: "+runs+" runs");
		}
	}

}
