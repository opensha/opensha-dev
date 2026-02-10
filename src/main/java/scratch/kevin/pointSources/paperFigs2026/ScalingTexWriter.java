package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import scratch.kevin.latex.LaTeXUtils;

public class ScalingTexWriter {

	public static void main(String[] args) throws IOException {
		FileWriter fw = new FileWriter(new File(FIGURES_DIR, "scaling_stats.tex"));
		
		List<Function<Double, Double>> mls = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		mls.add(LEONARD_DIP_SLIP_ML);
		prefixes.add("LeonardDip");
		mls.add(LEONARD_STRIKE_SLIP_ML);
		prefixes.add("LeonardSS");
		mls.add(WC94_RLD_ML);
		prefixes.add("WCRLD");
		mls.add(WC94_SRL_ML);
		prefixes.add("WCSRL");
		
		double[] mags = {4d, 5d, 6d, 7.05};
		String[] magNames = {"Four", "Five", "Six", "Sevenish"};
		
		DecimalFormat df = new DecimalFormat("0");
		
		for (int i=0; i<mls.size(); i++) {
			Function<Double, Double> ml = mls.get(i);
			
			for (int m=0; m<mags.length; m++) {
				String name = prefixes.get(i)+"M"+magNames[m];
				double len = ml.apply(mags[m]); 
				
				fw.write(LaTeXUtils.defineValueCommand(name, df.format(len))+"\n");
			}
		}
		
		fw.close();
	}

}
