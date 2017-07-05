package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.collect.Lists;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;

public class NonzeroWaterlevelPlot {
	
	private static ArbitrarilyDiscretizedFunc loadCounts(File dir) throws IOException {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		double[] nonZeros = null;
		int cnt = 0;
		for (File file : dir.listFiles()) {
			if (!file.getName().endsWith("noMinRates.bin"))
				continue;
			double[] rates = MatrixIO.doubleArrayFromFile(file);
			
			if (nonZeros == null)
				nonZeros = new double[rates.length];
			cnt++;
			for (int i=0; i<rates.length; i++)
				if (rates[i] > 0)
					nonZeros[i] = 1d;
			
			func.set((double)cnt, StatUtils.sum(nonZeros));
		}
		func.setInfo(nonZeros.length+" rups");
		return func;
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
//		File binsDir = new File("/home/kevin/OpenSHA/UCERF3/inversions/2013_05_09-ucerf3p3-convergence-test/bins");
//		File startUCERF2Dir = new File("/home/kevin/OpenSHA/UCERF3/nomins_startucerf2");
		
//		ArbitrarilyDiscretizedFunc refBranch = loadCounts(binsDir);
//		refBranch.setName("UCERF3.3 Reference Branch");
//		ArbitrarilyDiscretizedFunc startUCERF2Func = loadCounts(startUCERF2Dir);
//		startUCERF2Func.setName("Starting with UCERF2");
		
//		int numRups = Integer.parseInt(refBranch.getInfo().split(" ")[0]);
		
		File solFile = new File("/tmp/branch_avg_avg/mean_noMins.zip");
		AverageFaultSystemSolution avgSol = FaultSystemIO.loadAvgInvSol(solFile);
		EvenlyDiscretizedFunc refBranch = new EvenlyDiscretizedFunc(1d, avgSol.getNumSolutions(), 1d);
		int numRups = avgSol.getRupSet().getNumRuptures();
		double[] nonZeros = new double[numRups];
		for (int i=0; i<avgSol.getNumSolutions(); i++) {
			double[] rates = avgSol.getRates(i);
			for (int r=0; r<rates.length; r++) {
				if (rates[r] > 0)
					nonZeros[r] = 1;
			}
			refBranch.set(i, StatUtils.sum(nonZeros));
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(refBranch);
//		GraphWindow gw = new GraphWindow(funcs, "Num Non Zeros");
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		gp.setBackgroundColor(Color.WHITE);
		gp.setUserBounds(refBranch.getMinX(), refBranch.getMaxX(), 0, numRups);
		gp.drawGraphPanel("# Runs", "# Ruptures", funcs, chars,
				"Ruptures Above Waterlevel");
		File file = new File("/tmp/compound_rups_above_waterlevel");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");

		file = new File(file.getAbsolutePath()+"_small");
		gp.getChartPanel().setSize(500, 400);
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}

}
