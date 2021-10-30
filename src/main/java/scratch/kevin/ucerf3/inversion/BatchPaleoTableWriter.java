package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoProbabilityModel;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.inversion.UCERF2_ComparisonSolutionFetcher;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoFitPlotter;
import scratch.UCERF3.utils.paleoRateConstraints.U3PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

import com.google.common.base.Preconditions;

public class BatchPaleoTableWriter {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		Preconditions.checkArgument(args.length == 1, "Must specify directory.");
		File dir = new File(args[0]);
		Preconditions.checkArgument(dir.exists() && dir.isDirectory(),
				"Directory doesn't exist or isn't directory.");
		
		FaultSystemSolution ucerf2Sol = UCERF2_ComparisonSolutionFetcher
				.getUCERF2Solution(FaultModels.FM2_1);
		List<U3AveSlipConstraint> ucerf2AveSlipConstraints;
		List<U3PaleoRateConstraint> ucerf2PaleoConstraints;
		try {
			ucerf2AveSlipConstraints = U3AveSlipConstraint.load(
					ucerf2Sol.getRupSet().getFaultSectionDataList());
			ucerf2PaleoConstraints = UCERF3_PaleoRateConstraintFetcher
					.getConstraints(ucerf2Sol.getRupSet().getFaultSectionDataList());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		PaleoProbabilityModel paleoProbModel = null;
		try {
			paleoProbModel = UCERF3_PaleoProbabilityModel.load();
		} catch (IOException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		
		handleDir(dir, ucerf2Sol, ucerf2AveSlipConstraints, ucerf2PaleoConstraints, paleoProbModel);
	}
	
	private static void handleDir(File dir, FaultSystemSolution ucerf2Sol,
			List<U3AveSlipConstraint> ucerf2AveSlipConstraints,
			List<U3PaleoRateConstraint> ucerf2PaleoRateConstraints,
			PaleoProbabilityModel paleoProbModel) throws ZipException, IOException, DocumentException {
		for (File file : dir.listFiles()) {
			if (file.isDirectory()) {
				handleDir(file, ucerf2Sol, ucerf2AveSlipConstraints,
						ucerf2PaleoRateConstraints, paleoProbModel);
				continue;
			}
			
			// make sure it's a solution file
			if (!file.getName().endsWith("_sol.zip"))
				continue;
			
			File paleoDir = new File(dir, "paleo_fault_based");
			if (!paleoDir.exists())
				paleoDir.mkdir();
			
			InversionFaultSystemSolution sol = U3FaultSystemIO.loadInvSol(file);
			
			List<U3AveSlipConstraint> aveSlipConstraints =
					U3AveSlipConstraint.load(sol.getRupSet().getFaultSectionDataList());
			
			List<U3PaleoRateConstraint> paleoRateConstraints =
					CommandLineInversionRunner.getPaleoConstraints(
							sol.getRupSet().getFaultModel(), sol.getRupSet());
			
			PaleoFitPlotter.writeTables(paleoDir, sol, aveSlipConstraints, paleoRateConstraints,
					ucerf2Sol, ucerf2AveSlipConstraints, ucerf2PaleoRateConstraints, paleoProbModel);
		}
	}

}
