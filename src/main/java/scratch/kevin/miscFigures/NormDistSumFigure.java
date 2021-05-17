package scratch.kevin.miscFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import com.google.common.base.Preconditions;

public class NormDistSumFigure {

	public static void main(String[] args) throws IOException {
		int[] allNumDists =  { 10, 50, 100, 1000 };
		File outputDir = new File("/tmp");
		
		RandomGenerator rng = new Well19937c(12345l);
		NormalDistribution normDist = new NormalDistribution(rng, 0d, 1d);
		NormalDistribution meanDist = new NormalDistribution(rng, 0d, 1d);
//		RealDistribution sqDist = new NormalDistribution(rng, 1d, 0.25d);
//		RealDistribution sqDist = new ExponentialDistribution(1d);
		
		double minX = -3;
		double maxX = 3;
		int numX = 500;
//		double inverseCumAccuracy = (maxX - minX)/(double)numX;
//		inverseCumAccuracy /= 10;
		
		for (int numDists : allNumDists) {
			System.out.println("Num Dists: "+numDists);
			
			NormalDistribution[] randDists = new NormalDistribution[numDists];
			NormalDistribution[] uniformDists = new NormalDistribution[numDists];
			NormalDistribution[] spacedDists = new NormalDistribution[numDists];
			
			double deltaInvP = 1d/(double)numDists;

			double sumRandMean = 0d;
			double sumSpacedMean = 0d;
			double sumSigmaSq = 0d;
			for (int i=0; i<numDists; i++) {
				uniformDists[i] = normDist;
				
				double mean, sigma;
//				if (i < numDists - 1) {
					mean = meanDist.sample();
//					sigma = Double.NaN;
//					while (Double.isNaN(sigma))
//						// force it to be positive
//						sigma = Math.sqrt(sqDist.sample());
//				} else {
//					System.out.println("Sum mean before last: "+sumMean);
//					mean = -sumMean;
//					System.out.println("Sum sigma^2 before last: "+sumSigmaSq);
//					sigma = Math.sqrt(numDists - sumSigmaSq);
//				}
//				System.out.println(i+". mean="+mean+"\tsigma="+sigma);
				sigma = Math.sqrt(1d/(double)numDists);
//				sigma = 1d;
				sumRandMean += mean;
				sumSigmaSq += sigma*sigma;
				randDists[i] = new NormalDistribution(rng, mean, sigma);
				
				// now spaced
				double p = 0.5*deltaInvP + i*deltaInvP;
				double spacedMean = normDist.inverseCumulativeProbability(p);
//				System.out.println("p="+p+"\tmean="+spacedMean);
				spacedDists[i] = new NormalDistribution(rng, spacedMean, sigma);
			}
			System.out.println("Sum rand mean: "+sumRandMean);
			System.out.println("Sum spaced mean: "+sumSpacedMean);
			System.out.println("Sum sigma^2: "+sumSigmaSq);

			System.out.println("Plotting uniform");
			Range yRange = plotStackedHists(minX, maxX, numX, uniformDists, outputDir,
					"stacked_uniform_dists_"+numDists, null);
			System.out.println("Plotting random");
			plotStackedHists(minX, maxX, numX, randDists, outputDir,
					"stacked_rand_dists_"+numDists, yRange);
			System.out.println("Plotting spaced");
			plotStackedHists(minX, maxX, numX, spacedDists, outputDir,
					"stacked_spaced_dists_"+numDists, yRange);
			System.out.println("DONE");
		}
	}
	
	private static final DecimalFormat df = new DecimalFormat("0.00");
	
	private static Range plotStackedHists(double minX, double maxX, int numX,
			NormalDistribution[] dists, File outputDir, String prefix, Range yRange)
					throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		HistogramFunction sumFunc = new HistogramFunction(minX, maxX, numX);
		double delta = sumFunc.getDelta();
		double halfDelta = 0.5*delta;
		
		for (int i=0; i<dists.length; i++) {
			String name = "μ="+df.format(dists[i].getMean())
				+", σ="+df.format(dists[i].getStandardDeviation());
//			System.out.println("Processing dist "+i+": "+name);
			for (int j=0; j<numX; j++) {
				double x = sumFunc.getX(j);
				double x0 = x-halfDelta;
				double x1 = x+halfDelta;
				double prob = dists[i].probability(x0, x1);
				Preconditions.checkState(prob >= 0, "Bad prob=%s for x=[%s %s]", prob, x0, x1);
				double density = prob/delta;
				density /= (double)dists.length;
				sumFunc.add(j, density);
			}
			DiscretizedFunc distFunc = sumFunc.deepClone();
			distFunc.setName(name);
			funcs.add(distFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
		}
		
		sumFunc.setName("Sum");
		funcs.add(sumFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "", "X", "Probability Density");
//		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(minX, maxX);
		if (yRange == null)
			yRange = new Range(0d, sumFunc.getMaxY()*1.25);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getYAxis().setTickLabelsVisible(false);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		return yRange;
	}

}
