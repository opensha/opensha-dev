package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;

class FaultParticipationPlots {

	public static void main(String[] args) throws IOException {
		File solFile = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-6, -1);
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setNanColor(Color.GRAY);
		
		double[] minMags = { 0, 6, 6.7, 7.5, 7.8 };
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/partic_plot");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(sol.getRupSet(), NSHM23_RegionLoader.loadFullConterminousWUS());
		
		DecimalFormat oDF = new DecimalFormat("0.#");
		
		for (double minMag : minMags) {
			String label, prefix;
			
			if (minMag > 0d) {
				label = "M≥"+oDF.format(minMag);
				prefix = "m"+oDF.format(minMag);
			} else {
				label = "Supra-Seismogenic";
				prefix = "supra_seis";
			}
			prefix = "paric_"+prefix;
			label = "Log10 "+label+" Participation Rate (/yr)";
			
			double[] rates = sol.calcParticRateForAllSects(minMag, Double.POSITIVE_INFINITY);
			for (int i=0; i<rates.length; i++) {
				if (rates[i] == 0d)
					rates[i] = Double.NaN;
				else
					rates[i] = Math.log10(rates[i]);
			}
			
			mapMaker.plotSectScalars(rates, cpt, label);
			
			mapMaker.plot(outputDir, prefix, " ");
		}
	}

}
