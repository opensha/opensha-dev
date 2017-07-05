package scratch.kevin.ucerf3.eal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.FaultSystemIO;

public class InternInterestingRuptureFinder {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		// true mean FSS which includes rupture mapping information. this must be the exact file used to calulate EALs
		File trueMeanSolFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		
//		int numToWrite = 1000;

		// dollars from CEA proxy portfolio
//		File dataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_01_15-ucerf3-eal-calc-NGA2s-2013");
		File dataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller");
		String label = "$ (Billions)";
		final double scale = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		File outputFile = new File("/tmp/rup_losses.csv");
		
		// fatalities
//		File dataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-fatality-smaller");
//		String label = "Fatalities";
//		final double scale = 1d;
//		File outputFile = new File("/tmp/rup_fatalities.csv");
		
		// IMR for which EAL data has already been computed
		AttenRelRef attenRelRef = AttenRelRef.BSSA_2014;

		// Fault model of interest
		FaultModels fm = FaultModels.FM3_1;
		
		// Compound fault system solution
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		// Branch averaged FSS
		FaultSystemSolution baSol = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
						+ "COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		double minMag = 6.7;
		double maxMag = 7.1;
		
		UCERF3_BranchAvgLossFetcher dollarFetch = new UCERF3_BranchAvgLossFetcher(trueMeanSolFile, cfss, dataDir);
		
		DiscretizedFunc[] faultLosses = dollarFetch.getFaultLosses(attenRelRef, fm, true);
		
		List<Rupture> rups = Lists.newArrayList();
		
		for (int i = 0; i < faultLosses.length; i++) {
			double mag = baSol.getRupSet().getMagForRup(i);
			if (mag < minMag || mag > maxMag)
				continue;
			double rate = baSol.getRateForRup(i);
			DiscretizedFunc faultLoss = faultLosses[i];
			double loss = calcMeanLoss(faultLoss)*scale;
			
			Rupture rup = new Rupture(i, mag, loss, rate);
			rups.add(rup);
		}
		
		Collections.sort(rups);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Index", "Mag", "Loss: "+label, "Rate", "Loss*Rate", "First Section", "Last Section");
		for (int i = 0; i < rups.size(); i++) {
			Rupture rup = rups.get(i);
			List<FaultSectionPrefData> sects = baSol.getRupSet().getFaultSectionDataForRupture(rup.id);
			csv.addLine(rup.id+"", rup.mag+"", rup.loss+"", rup.rate+"", (rup.loss*rup.rate)+"",
					sects.get(0).getName(), sects.get(sects.size()-1).getName());
		}
		csv.writeToFile(outputFile);
	}
	
	private static final double calcMeanLoss(DiscretizedFunc faultLoss) {
		double loss = 0;
		for (Point2D pt : faultLoss)
			// x=loss, y=weight
			loss += pt.getX()*pt.getY();
		return loss;
	}
	
	static class Rupture implements Comparable<Rupture> {
		int id;
		double mag;
		double loss;
		double rate;
		
		public Rupture(int id, double mag, double loss, double rate) {
			super();
			this.id = id;
			this.mag = mag;
			this.loss = loss;
			this.rate = rate;
		}

		@Override
		public int compareTo(Rupture o) {
			return Double.compare(o.loss, loss);
		}
	}

}
