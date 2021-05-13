package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.Region;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.RELM_RegionUtils;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF2_MFD_ConstraintFetcher;

public class BulkAvgMFDCalc {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_06_27-ref-char-unconst");
		String prefix = "FM3_1_NEOK_EllB_DsrUni_CharUnconst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3";
		InversionFaultSystemRupSet rupSet = FaultSystemIO.loadInvRupSet(
				new File(dir, "FM3_1_NEOK_EllB_DsrUni_CharUnconst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_run000_sol.zip"));
		AverageFaultSystemSolution avgSol = AverageFaultSystemSolution.fromDirectory(rupSet, dir, prefix);
		
		Region region = RELM_RegionUtils.getGriddedRegionInstance();
		UCERF2_MFD_ConstraintFetcher ucerf2Fetch = new UCERF2_MFD_ConstraintFetcher();
		
		int digits = ((avgSol.getNumSolutions()-1)+"").length();
		
		for (int i=0; i<avgSol.getNumSolutions(); i++) {
			System.out.println("Plotting "+i);
			double[] rates = avgSol.getRates(i);
			InversionFaultSystemSolution sol = new InversionFaultSystemSolution(rupSet, rates);
			String numStr = i+"";
			while (numStr.length() < digits)
				numStr = "0"+numStr;
			String indPrefix = prefix+"_run"+numStr;
			IncrementalMagFreqDist totalMFD = rupSet.getInversionTargetMFDs().getTotalTargetGR();
			IncrementalMagFreqDist targetMFD = rupSet.getInversionTargetMFDs().getOnFaultSupraSeisMFD();
			IncrementalMagFreqDist solutionMFD = sol.calcNucleationMFD_forRegion(null, // null since we want everything
					totalMFD.getMinX(), 9.05, 0.1, true);
			double minMag = rupSet.getMinMag();
			double maxMag = avgSol.getInversionConfiguration().getMFDTransitionMag();
			double e = 0;
			for (int j=0; j<targetMFD.size(); j++) {
				double x = targetMFD.getX(j);
				double targetY = targetMFD.getY(j);
				if (x < minMag || x > maxMag)
					continue;
				double solutionY = solutionMFD.getY(x);
				double e_part = (solutionY - targetY) / targetY * 10d;
				System.out.println(x+":\ttarget="+targetY+"\tsolution="+solutionY+"\te_part="+e_part+"\tenergy="+(e_part*e_part));
				if (targetY == 0)
					continue;
				e += e_part * e_part;
			}
			System.out.println("MISFIT: "+e);
			CommandLineInversionRunner.writeMFDPlot(sol, dir, indPrefix, totalMFD,
					targetMFD, region, ucerf2Fetch);
		}
	}

}
