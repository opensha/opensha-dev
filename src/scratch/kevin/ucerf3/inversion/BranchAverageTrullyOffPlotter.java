package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.collect.Lists;

import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class BranchAverageTrullyOffPlotter {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File u3p2File = new File("/tmp/2013_01_14-stampede_3p2_production_runs_fm3p1_dm_scale_subset_MEAN_BRANCH_AVG_SOL.zip");
		File u3p3File = new File("/tmp/2013_05_01-ucerf3p3-proposed-subset-hpcc-salmonfix_COMPOUND_SOL_MEAN_BRANCH_AVG_SOL.zip");
		
		InversionFaultSystemSolution u3p2Sol = FaultSystemIO.loadInvSol(u3p2File);
		InversionFaultSystemSolution u3p3Sol = FaultSystemIO.loadInvSol(u3p3File);
		
		SummedMagFreqDist u3p2TrulyOffMFD = new SummedMagFreqDist(0.05, 100, 0.1);
		SummedMagFreqDist u3p3TrulyOffMFD = new SummedMagFreqDist(0.05, 100, 0.1);
		
		GridSourceProvider u3p2Grid = u3p2Sol.getGridSourceProvider();
		GridSourceProvider u3p3Grid = u3p3Sol.getGridSourceProvider();
		
		for (int i=0; i<u3p2Grid.size(); i++) {
			IncrementalMagFreqDist u3p2MFD = u3p2Grid.getNodeUnassociatedMFD(i);
			IncrementalMagFreqDist u3p3MFD = u3p3Grid.getNodeUnassociatedMFD(i);
			
			if (u3p2MFD != null)
				u3p2TrulyOffMFD.addIncrementalMagFreqDist(u3p2MFD);
				
			if (u3p3MFD != null)
				u3p3TrulyOffMFD.addIncrementalMagFreqDist(u3p3MFD);
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(u3p2TrulyOffMFD);
		funcs.add(u3p3TrulyOffMFD);
		
		new GraphWindow(funcs, "Truly Off Fault MFDs");
	}

}
