package scratch.kevin.nshm27;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.TruncatedNormalDistribution;
import org.apache.commons.statistics.distribution.UniformContinuousDistribution;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeLevel.ContinuousDistributionBinnedLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.ContinuousDistributionSampledLevel;
import org.opensha.commons.logicTree.LogicTreeNode.SimpleValuedNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

import com.google.common.collect.Range;

import net.mahdilamb.colormap.Colors;

public class DistSampleCountTests {

	public static void main(String[] args) throws IOException {
//		ContinuousDistribution dist = TruncatedNormalDistribution.of(7.6, 0.134, 7.15, 8.05);
		ContinuousDistribution dist = TruncatedNormalDistribution.of(8, 0.2, 7.45, 8.55);
		int precisionScale = 1;
		int bins = 3;
		
////		ContinuousDistribution dist = UniformContinuousDistribution.of(0d, 1d);
////		ContinuousDistribution dist = UniformContinuousDistribution.of(0.5d, 1d);
//		ContinuousDistribution dist = TruncatedNormalDistribution.of(0.75, 0.5, 0d, 1.5);
//		int precisionScale = 2;
		
//		ContinuousDistribution dist = UniformContinuousDistribution.of(1000d, 1500d);
//		int precisionScale = 2;
		
		int samples = 10000;
		
		ContinuousDistributionSampledLevel level = new ContinuousDistributionSampledLevel(
				"Test level", "Level", dist, precisionScale, "Node ", "Node", "Node");
		
		level.build(123456789l, samples);
				
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (SimpleValuedNode<Double> node : level.getNodes())
			track.addValue(node.getValue());
		
		double delta = 1d/Math.pow(10, precisionScale);
		
		System.out.println("Values: "+track);
//		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(track.getMin(), track.getMax(), delta);
		EvenlyDiscretizedFunc hist = new EvenlyDiscretizedFunc(track.getMin(), track.getMax(), 1 + (int)(track.getLength()/delta + 0.5));
		
		for (SimpleValuedNode<Double> node : level.getNodes())
			hist.add(hist.getClosestXIndex(node.getValue()), 1d);
		
		hist.scale(1d/(samples*hist.getDelta()));
		
		EvenlyDiscretizedFunc pdfDensity = new EvenlyDiscretizedFunc(track.getMin()-0.5*delta, track.getMax()+0.5*delta, 1000);
		
		for (int i=0; i<pdfDensity.size(); i++) {
			double x = pdfDensity.getX(i);
			double density = dist.density(x);
			if (density > 0)
				pdfDensity.set(i, density);
		}
		
		EvenlyDiscretizedFunc pdfCoarseDensity = new EvenlyDiscretizedFunc(hist.getMinX(), hist.getMaxX(), hist.size());
		for (int i=0; i<pdfCoarseDensity.size(); i++) {
			double x = pdfCoarseDensity.getX(i);
			double density = dist.density(x);
			if (density > 0)
				pdfCoarseDensity.set(i, density);
		}
		
		double histSum = hist.calcSumOfY_Vals();
		for (int i=0; i<hist.size(); i++) {
			double middle = hist.getX(i);
			double start = middle - 0.5*delta;
			double end = middle + 0.5*delta;
			double prob = dist.cumulativeProbability(end) - dist.cumulativeProbability(start);
			System.out.println("Bin "+i+"; center="+(float)middle+", ["+(float)start+", "+(float)end+"]");
			System.out.println("\tSampled weight:\t"+(float)(hist.getY(i)/histSum));
			System.out.println("\tIntegrated prob:\t"+(float)prob);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_blue));
		
		funcs.add(pdfDensity);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		funcs.add(pdfCoarseDensity);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
		
		for (int i=0; i<=3; i++) {
			double x = dist.inverseCumulativeProbability(i/3d);
			double y = dist.density(x);
			DefaultXY_DataSet third = new DefaultXY_DataSet(x, 0, x, y);
			
			funcs.add(third);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_green));
			
			third = new DefaultXY_DataSet(x, y);
			
			funcs.add(third);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Colors.tab_green));
		}
		
		ContinuousDistributionBinnedLevel binned = bins > 0 ? level.toBinnedLevel(bins) : level.toBinnedLevel();
		List<? extends SimpleValuedNode<Range<Double>>> binNodes = binned.getNodes();
		for (int i=0; i<binNodes.size(); i++) {
			SimpleValuedNode<Range<Double>> node = binNodes.get(i);
			System.out.println("Bin node "+i+": "+node.getName()+" (wt="+(float)node.getNodeWeight()+")");
			if (i == binNodes.size()-1)
				break;
			double x = node.getValue().upperEndpoint();
			double y = dist.density(x);
			DefaultXY_DataSet edge = new DefaultXY_DataSet(x, 0, x, y);
			
			funcs.add(edge);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Colors.tab_orange));
			
			edge = new DefaultXY_DataSet(x, y);
			
			funcs.add(edge);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Colors.tab_orange));
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, samples+" samples", "X", "Density");
		
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		
		gp.drawGraphPanel(plot);
		
		PlotUtils.writePlots(new File("/tmp"), "dist_test", gp, 800, 650, true, false, false);
	}

}
