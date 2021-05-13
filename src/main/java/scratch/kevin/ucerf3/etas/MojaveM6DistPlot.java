package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.cpt.CPT;

import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class MojaveM6DistPlot {

	public static void main(String[] args) throws IOException {
		String dirPrefix = "U3_ETAS_KevinTestMojave_";
		int[] dists = { 0, 1, 2, 4, 6, 8, 12, 16 };
		String mfdFileName = "Diagnostics_Scenario/Scenario_ExpPrimMFD_Cum.txt";
		String yLabel = "Primary Expected Cumulative Trigger Rate";
		double[] mags = { 6d, 6.5d, 7d };
		
		ArbitrarilyDiscretizedFunc[] funcs = new ArbitrarilyDiscretizedFunc[mags.length];
		for (int i=0; i<mags.length; i++)
			funcs[i] = new ArbitrarilyDiscretizedFunc();
		
		DecimalFormat df = new DecimalFormat("0.0");
		String[] searches = new String[mags.length];
		for (int i=0; i<mags.length; i++)
			searches[i] = df.format(mags[i]);
		
		for (int dist : dists) {
			File file = new File(dirPrefix+dist+"km/"+mfdFileName);
			if (!file.exists()) {
				System.out.println(dist+"km not yet done, skipping");
				continue;
			}
			
			boolean[] dones = new boolean[mags.length];
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				line = line.trim();
				for (int i=0; i<mags.length; i++) {
					if (!dones[i] && line.startsWith(searches[i])) {
						double val = Double.parseDouble(line.split("\t")[1].trim());
						System.out.println(dist+"km, M≥"+(float)mags[i]+": "+val);
						funcs[i].set((double)dist, val);
						dones[i] = true;
						// should check if all done and break but inexpensive and lazy
					}
				}

			}
		}
		
		List<DiscretizedFunc> funcsList = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		CPT cpt = new CPT(mags[0], mags[mags.length-1], Color.BLUE, Color.RED);
		for (int i=0; i<mags.length; i++) {
			funcs[i].setName("M≥"+(float)mags[i]);
			funcsList.add(funcs[i]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, cpt.getColor((float)mags[i])));
		}
		PlotSpec spec = new PlotSpec(funcsList, chars, "Mojave Trigger Dist Decay",
				"Fault Perpendicular Distance (km)", yLabel);
		spec.setLegendVisible(true);
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		gw.saveAsPNG(new File("/tmp/mojave_dist_trigger.png").getAbsolutePath());
		gw.saveAsPDF(new File("/tmp/mojave_dist_trigger.pdf").getAbsolutePath());
	}

}
