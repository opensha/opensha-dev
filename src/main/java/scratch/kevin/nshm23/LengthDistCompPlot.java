package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.TablesAndPlotsGen;

public class LengthDistCompPlot {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_09_28-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File outputDir = new File(invDir, "misc_plots");
		String prefix = "wells_2013_length_dist_compare";
		
		File resultsFile = new File(invDir, "results.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
//		tree = tree.sample(10, false);
		
		double[] lengths = slt.forBranch(tree.getBranch(0)).getRupSet().getLengthForAllRups();
		
		Range xRange = new Range(0d, 500d);
		HistogramFunction meanHist = HistogramFunction.getEncompassingHistogram(
				xRange.getLowerBound(), xRange.getUpperBound(), 50);
		
		double sumWeight = 0d;
		HistogramFunction lowerHist = new HistogramFunction(meanHist.getMinX(), meanHist.size(), meanHist.getDelta());
		HistogramFunction upperHist = new HistogramFunction(meanHist.getMinX(), meanHist.size(), meanHist.getDelta());
		for (int index=0; index<tree.size(); index++) {
			LogicTreeBranch<?> branch = tree.getBranch(index);
			System.out.println("Processing branch "+index+"/"+tree.size());
			double weight = tree.getBranchWeight(branch);
			sumWeight += weight;
			
			double[] rates = slt.loadRatesForBranch(branch);
			Preconditions.checkState(rates.length == lengths.length);
			
			HistogramFunction myHist = new HistogramFunction(meanHist.getMinX(), meanHist.size(), meanHist.getDelta());
			for (int r=0; r<rates.length; r++) {
				if (rates[r] > 0) {
					double len = lengths[r]*1e-3; // m -> km
					if (len > myHist.getMaxX()+0.5*myHist.getDelta())
						continue;
					int lenIndex = myHist.getClosestXIndex(len);
					myHist.add(lenIndex, rates[r]);
				}
			}
			myHist.normalizeBySumOfY_Vals();
			for (int i=0; i<myHist.size(); i++) {
				double val = myHist.getY(i);
				if (val > 0d) {
					meanHist.add(i, weight*myHist.getY(i));
					if (val > upperHist.getY(i))
						upperHist.set(i, val);
					if (val < lowerHist.getY(i) || lowerHist.getY(i) == 0d)
						lowerHist.set(i, val);
				}
			}
		}
		meanHist.scale(1d/sumWeight);
		
		HistogramFunction data = TablesAndPlotsGen.loadSurfaceRupData();
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		data.setName("Wells (2013), Observed");
		funcs.add(data);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.LIGHT_GRAY));
		
		meanHist.setName("NSHN23 WUS, Mean");
		funcs.add(meanHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, Color.BLACK));
		
		lowerHist.setName("NSHN23 WUS, Min & Max");
		funcs.add(lowerHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		upperHist.setName(null);
		funcs.add(upperHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Rupture Length Distribution", "Rupture Length (km)",
				"Fraction of Earthquakes");
		spec.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, false, xRange, new Range(0d, 1d));
		
		PlotUtils.writePlots(outputDir, prefix, gp, 1000, 800, true, true, true);
		
		gp.drawGraphPanel(spec, false, true, xRange, new Range(1e-3, 1d));
		
		PlotUtils.writePlots(outputDir, prefix+"_log", gp, 1000, 800, true, true, true);
	}

}
