package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

public class FaultParticipationPlots {
	
	public static void main(String[] args) throws IOException {
		File solFile = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2024_11_19-prvi25_crustal_branches-dmSample5x/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-5, -1);
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setNanColor(Color.GRAY);
		
		double[] minMags = { 0, 6, 6.7, 7, 7.5, 7.8 };
		
//		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/partic_plot");
		File outputDir = new File("/tmp");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(sol.getRupSet(), PRVI25_RegionLoader.loadPRVI_ModelBroad());
		
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
