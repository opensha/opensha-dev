package scratch.kevin.nshm26;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.TruncatedNormalDistribution;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.apache.commons.statistics.distribution.ContinuousDistribution.Sampler;
import org.apache.commons.statistics.distribution.CorrTruncatedNormalDistribution;
import org.apache.commons.statistics.distribution.NormalDistribution;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

import net.mahdilamb.colormap.Colors;

public class SamplerTest {

	public static void main(String[] args) throws IOException {
//		double mean = 1d;
//		double sd = 0.1;
//		double lower = 0.7;
//		double upper = 1.3;
		
//		double mean = -1d;
//		double sd = 0.1;
//		double lower = -1.3;
//		double upper = -0.7;

//		double mean = 1d;
//		double sd = 0.1;
//		double lower = 0.7;
//		double upper = 1.3;
		
		double mean = 1d;
		double sd = 0.1;
		double lower = -0.1;
		double upper = 1.3;
		
//		double mean = -1d;
//		double sd = 0.1;
//		double lower = -1d;
//		double upper = -0.7;
		
//		double mean = -0.1d;
//		double sd = 0.1;
//		double lower = 0.01;
//		double upper = 1d;
		
		int samples = 10000000;
		
		EvenlyDiscretizedFunc hist1 = HistogramFunction.getEncompassingHistogram(lower-0.1, upper+0.1, 0.01);
		
		ContinuousDistribution dist = TruncatedNormalDistribution.of(mean, sd, lower, upper);
		
		System.out.println("Reported dist params:");
		System.out.println("\tmean: "+dist.getMean());
		System.out.println("\tvar: "+dist.getVariance());
		System.out.println("\tsqrt(var): "+Math.sqrt(dist.getVariance()));
		System.out.println("\tlower: "+dist.getSupportLowerBound());
		System.out.println("\tupper: "+dist.getSupportUpperBound());
		
		Sampler sampler = dist.createSampler(RandomSource.XO_RO_SHI_RO_128_PP.create(123456l));
		for (int i=0; i<samples; i++) {
			double value = sampler.sample();
			hist1.add(hist1.getClosestXIndex(value), 1d);
		}
		hist1.scale(1d/(samples*hist1.getDelta()));
		
		EvenlyDiscretizedFunc hist2 = new EvenlyDiscretizedFunc(hist1.getMinX(), hist1.size(), hist1.getDelta());
		ContinuousDistribution corrDist = CorrTruncatedNormalDistribution.of(mean, sd, lower, upper);
		Sampler corrSampler = corrDist.createSampler(RandomSource.XO_RO_SHI_RO_128_PP.create(123456l));
		for (int i=0; i<samples; i++) {
			double value = corrSampler.sample();
			hist2.add(hist2.getClosestXIndex(value), 1d);
		}
		hist2.scale(1d/(samples*hist2.getDelta()));
		
		ArbitrarilyDiscretizedFunc pdfDensity = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<hist1.size(); i++) {
			double x = hist1.getX(i);
			double density = dist.density(x);
			if (density > 0)
				pdfDensity.set(x, density);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		hist1.setName("Current Rejection Sampler");
		funcs.add(hist1);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		Color corrHistColor = Colors.tab_orange;
		corrHistColor = new Color(corrHistColor.getRed(), corrHistColor.getGreen(), corrHistColor.getBlue(), 127);
		hist2.setName("Corrected Rejection Sampler");
		funcs.add(hist2);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, corrHistColor));
		
		pdfDensity.setName("PDF Density");
		funcs.add(pdfDensity);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Colors.tab_blue));
		
		double maxY = 0d;
		for (DiscretizedFunc func : funcs)
			maxY = Math.max(maxY, func.getMaxY());
		
		String title = "μ="+(float)mean+", σ="+(float)sd+", bounds=["+(float)lower+", "+(float)upper+"]";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, "X", "Density");
		plot.setLegendInset(RectangleAnchor.TOP_LEFT);
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		
		gp.setLegendFontSize(16);
		
		gp.drawGraphPanel(plot, false, false, new Range(hist1.getMinX(), hist1.getMaxX()), new Range(0d, maxY*1.2));
		
		String prefix = "trunc_dist_"+(float)mean+"_"+(float)sd+"_"+(float)lower+"_"+(float)upper;
		
		PlotUtils.writePlots(new File("/tmp"), prefix, gp, 700, 450, true, false, false);
	}

}
