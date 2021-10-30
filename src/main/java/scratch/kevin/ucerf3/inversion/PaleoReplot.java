package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.dom4j.DocumentException;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.paleoRateConstraints.U3PaleoRateConstraint;

public class PaleoReplot {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File(args[0]);
		
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				continue;
			if (!file.getName().endsWith("_sol.zip"))
				continue;
			System.out.println("Working on: "+file.getName());
			InversionFaultSystemSolution sol = U3FaultSystemIO.loadInvSol(file);
			String prefix = file.getName().substring(0, file.getName().indexOf("_sol.zip"));
			
			U3LogicTreeBranch branch = U3LogicTreeBranch.fromFileName(file.getName());
			
			ArrayList<U3PaleoRateConstraint> paleoRateConstraints = CommandLineInversionRunner.getPaleoConstraints(branch.getValue(FaultModels.class),
					sol.getRupSet());
			
			CommandLineInversionRunner.writePaleoPlots(paleoRateConstraints, null, sol, dir, prefix);
		}
	}

}
