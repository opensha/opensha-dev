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

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class FaultParticipationPlots {
	
	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(CRUSTAL_SOL_SUPRA_ONLY);
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-5, -1);
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setNanColor(Color.GRAY);
		cpt.setLog10(true);
		
		double[] minMags = { 0, 6, 6.7, 7, 7.5, 7.8 };
		
		File outputDir = new File(FIGURES_DIR, "partic_plot");
//		File outputDir = new File("/tmp");
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
			label = label+" Participation Rate (/yr)";
			
			double[] rates = sol.calcParticRateForAllSects(minMag, Double.POSITIVE_INFINITY);
			for (int i=0; i<rates.length; i++) {
				if (rates[i] == 0d)
					rates[i] = Double.NaN;
//				else
//					rates[i] = Math.log10(rates[i]);
			}
			
			mapMaker.plotSectScalars(rates, cpt, label);
			
			mapMaker.plot(outputDir, prefix, " ");
		}
	}

}
