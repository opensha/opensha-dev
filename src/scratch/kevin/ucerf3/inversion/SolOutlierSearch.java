package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.RegionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Joiner;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.utils.FaultSystemIO;

public class SolOutlierSearch {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws RuntimeException 
	 * @throws GMT_MapException 
	 */
	public static void main(String[] args) throws IOException, DocumentException, GMT_MapException, RuntimeException {
		File dir = new File(args[0]);
		File saveDir = new File(args[1]);
		File referenceSolFile = new File(args[2]);
		FaultSystemSolution referenceSol = FaultSystemIO.loadSol(referenceSolFile);
		Region region = new CaliforniaRegions.RELM_TESTING();
		
		Region bayRegion = new CaliforniaRegions.SF_BOX_GRIDDED();
		Region laRegion = new CaliforniaRegions.LA_BOX_GRIDDED();
		
		String urlBase = "http://opensha.usc.edu/ftp/kmilner/ucerf3/"+saveDir.getName()+"/";
		
		ArrayList<DiffRecord> recs = new ArrayList<SolOutlierSearch.DiffRecord>();
		
		double minMag = 6.7;
		double maxMag = 10;
		
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				continue;
			if (!file.getName().endsWith("_sol.zip"))
				continue;
			System.out.println("Working on: "+file.getName());
			FaultSystemSolution sol = FaultSystemIO.loadSol(file);
			String prefix = file.getName().substring(0, file.getName().indexOf("_sol.zip"));
			
			double[] newVals = sol.calcParticRateForAllSects(minMag, maxMag);
			double[] refVals = referenceSol.calcParticRateForAllSects(minMag, maxMag);
			
			FaultBasedMapGen.plotParticipationRates(sol, region, saveDir, prefix, false, minMag, maxMag);
			double diff = FaultBasedMapGen.plotParticipationDiffs(sol, referenceSol, region, saveDir, prefix,
					false, minMag, maxMag);
			FaultBasedMapGen.plotParticipationRatios(sol, referenceSol, region, saveDir, prefix, false, minMag, maxMag, true);
			
			String rateURL = urlBase+prefix+"_partic_rates_"+(float)+minMag+"+.png";
			String diffURL = urlBase+prefix+"_ref_partic_diff_"+(float)+minMag+"+.png";
			
			// now calc the diffs for regions;
			double[] diffs = new double[3];
			diffs[0] = diff;
			for (int i=0; i<newVals.length; i++) {
				FaultTrace trace = sol.getRupSet().getFaultSectionData(i).getFaultTrace();
				trace = FaultUtils.resampleTrace(trace, 10);
				double myDiff = newVals[i] - refVals[i];
				diffs[1] += myDiff * RegionUtils.getFractionInside(bayRegion, trace);
				diffs[2] += myDiff * RegionUtils.getFractionInside(laRegion, trace);
			}
			
			recs.add(new DiffRecord(diffs, prefix, rateURL, diffURL));
		}
		Collections.sort(recs);
		
		
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Name", "Total Diff", "Bay Area Diff", "LA Diff", "Rate Diff URL");
		for (DiffRecord rec : recs) {
			csv.addLine(rec.prefix, rec.diffs[0]+"", rec.diffs[1]+"", rec.diffs[2]+"", rec.diffURL);
		}
		csv.writeToFile(new File(saveDir, "partic_diffs.csv"));
		FileWriter fw = new FileWriter(new File(saveDir, "partic_diffs.txt"));
		for (int row=0; row<csv.getNumRows(); row++) {
			fw.write(Joiner.on('\t').join(csv.getLine(row))+"\n");
		}
		fw.close();
	}
	
	private static class DiffRecord implements Comparable<DiffRecord> {
		private double[] diffs;
		private String prefix, rateURL, diffURL;
		public DiffRecord(double[] diffs, String prefix, String rateURL, String diffURL) {
			this.diffs = diffs;
			this.prefix = prefix;
			this.rateURL = rateURL;
			this.diffURL = diffURL;
		}

		@Override
		public int compareTo(DiffRecord o) {
			return -Double.compare(diffs[0], o.diffs[0]);
		}
		
	}

}
