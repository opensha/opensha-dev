package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;

import com.google.common.collect.Lists;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class VarsPaleoBoundsPlot {
	
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File dir = new File("/tmp/zeng_ref_vars_paleo");
		
		InversionFaultSystemRupSet rupSet = null;
		
		List<InversionFaultSystemSolution> sols = Lists.newArrayList();
		
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				continue;
			if (!file.getName().endsWith("_sol.zip"))
				continue;
			
			if (file.getName().startsWith("FM3_2"))
				continue;
			
			System.out.println("Loading..."+file.getName());
			
			if (rupSet == null)
				rupSet = FaultSystemIO.loadInvRupSet(file);
			InversionFaultSystemSolution sol = FaultSystemIO.loadInvSol(file);
			
			sols.add(sol);
		}
		
		AverageFaultSystemSolution.writePaleoBoundsPlot(dir, "zeng_ref_vars_lowpaleo", sols);
	}

}
