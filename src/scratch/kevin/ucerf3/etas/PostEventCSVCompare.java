package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;

public class PostEventCSVCompare {

	public static void main(String[] args) throws IOException {
		File partialCSVfile = new File("/tmp/one_week_mag_num_cumulative.csv");
		File finalCSVfile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2018_04_05-santa_cruz_m5p3-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont/"
				+ "plots/one_week_mag_num_cumulative.csv");
		File outputDir = new File("/tmp");
		String prefix = "santa_cruz_sim_comparison";
		
		CSVFile<String> partialCSV = CSVFile.readFile(partialCSVfile, true);
		CSVFile<String> finalCSV = CSVFile.readFile(finalCSVfile, true);
		
		ArbitrarilyDiscretizedFunc partialMean = new ArbitrarilyDiscretizedFunc();
		partialMean.setName("Partial Mean");
		ArbitrarilyDiscretizedFunc finalMean = new ArbitrarilyDiscretizedFunc();
		finalMean.setName("Final Mean");
		
		ArbitrarilyDiscretizedFunc partialFract = new ArbitrarilyDiscretizedFunc();
		partialFract.setName("Partial Fract With ≥ Mag");
		ArbitrarilyDiscretizedFunc finalFract = new ArbitrarilyDiscretizedFunc();
		finalFract.setName("Final Fract With ≥ Mag");
		
		double maxY = 0d;
		double minY = Double.POSITIVE_INFINITY;
		
		for (int row=1; row<partialCSV.getNumRows(); row++) {
			double mag = Double.parseDouble(partialCSV.get(row, 0));
			double mean1 = Double.parseDouble(partialCSV.get(row, 1));
			double mean2 = Double.parseDouble(finalCSV.get(row, 1));
			double fract1 = Double.parseDouble(partialCSV.get(row, 10));
			double fract2 = Double.parseDouble(finalCSV.get(row, 10));
			
			if (fract1 > 0 && fract1 < minY)
				minY = fract1;
			if (fract2 > 0 && fract2 < minY)
				minY = fract2;
			maxY = Math.max(maxY, Math.max(mean1, mean2));
			
			partialMean.set(mag, mean1);
			finalMean.set(mag, mean2);
			
			partialFract.set(mag, fract1);
			finalFract.set(mag, fract2);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars= new ArrayList<>();
		
		funcs.add(partialMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		
		funcs.add(finalMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		funcs.add(partialFract);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		funcs.add(finalFract);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Partial vs Final Simulation Comparison", "Magnitude", "Num");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setUserBounds(partialMean.getMinX(), partialMean.getMaxX(), minY, maxY);
		gp.setLegendFontSize(20);

		ETAS_MultiSimAnalysisTools.setFontSizes(gp, 10);

		gp.drawGraphPanel(spec, false, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(outputDir, prefix + ".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix + ".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(outputDir, prefix + ".txt").getAbsolutePath());
	}

}
