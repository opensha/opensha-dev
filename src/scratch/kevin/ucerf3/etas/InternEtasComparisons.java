package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class InternEtasComparisons {

	public static void main(String[] args) throws IOException {
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined/results_descendents_m5_preserve.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined/results_m5_preserve.bin");
		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2017_04_13-parkfield-10yr-full_td-no_ert-combined/results_descendents_m4_preserve.bin");
		double targetMag = 7d;
		
		List<? extends List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(catalogFile, targetMag);
		
		double[] times = { 1d / (365.25 * 24), 1d / 365.25, 7d / 365.25, 30 / 365.25, 1d, 10d };
		
		int[] counts = new int[times.length];
		
		long ot = 1325419200000l;
		
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			// these are already filtered for the magnitude of interest
			if (catalog.isEmpty())
				continue;
			
			long firstMatch = catalog.get(0).getOriginTime();
			long delta = firstMatch - ot;
			double deltaYears = (double)delta / ProbabilityModelsCalc.MILLISEC_PER_YEAR;
			for (int i=0; i<times.length; i++)
				if (deltaYears <= times[i])
					counts[i]++;
		}
		
		System.out.println("Time (yr)\tTime Label\tProb");
		for (int i=0; i<times.length; i++) {
			String label = getTimeLabel(times[i]);
			double prob = (double)counts[i]/catalogs.size();
			System.out.println((float)times[i]+"\t"+label+"\t"+(float)prob);
		}
	}
	
	private static String getTimeLabel(double time) {
		String label;
		if (time < 1d) {
			int days = (int) (time * 365.25 + 0.5);
			if (days == 30)
				label = "1 mo";
			else if (days == 7)
				label = "1 wk";
			else if (time < 1d / 365.25) {
				int hours = (int) (time * 365.25 * 24 + 0.5);
				label = hours + " hr";
			} else
				label = days + " d";
		} else {
			label = (int) time + " yr";
		}
		return label;
	}

}
