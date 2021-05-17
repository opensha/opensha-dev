package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.utils.FaultSystemIO;

public class RupSmoothnessPlotTest {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File solFile = new File("/tmp/comp_plots/comp_plots_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
		
		CommandLineInversionRunner.writeRupPairingSmoothnessPlot(
				sol, "comp_plots", solFile.getParentFile());
	}

}
