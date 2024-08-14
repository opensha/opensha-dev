package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class MagExceedProbFigure {

	public static void main(String[] args) throws IOException {
		File sltFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/full_logic_tree.zip");
		File outputDir = new File("/tmp");
		
		String prefix = "mag_exceedance";
		
		SolutionLogicTree slt = SolutionLogicTree.load(sltFile);
		LogicTree<?> tree = slt.getLogicTree();
		
//		tree = tree.sample(10, false);
		
		double[] durations = { 30d, 15d, 5d, 1d };
		Color[] colors = { Color.RED.darker(), Color.GREEN.darker(), Color.BLUE.darker(), Color.CYAN.darker()};
		
		double[] percentiles = { 0d, 0.16d, 0.84d, 1d };
		String percentileStr = "p0,16,84,100";
		
		Range xRange = new Range(5d, 8.25d);
		
		EvenlyDiscretizedFunc refIncrFunc = FaultSysTools.initEmptyMFD(8.5);
		EvenlyDiscretizedFunc refCmlFunc = null;
		
		ArbDiscrEmpiricalDistFunc[][] dists = new ArbDiscrEmpiricalDistFunc[durations.length][refIncrFunc.size()];
		for (int d=0; d<durations.length; d++)
			for (int i=0; i<dists[d].length; i++)
				dists[d][i] = new ArbDiscrEmpiricalDistFunc();
		
		for (LogicTreeBranch<?> branch : tree) {
			double weight = tree.getBranchWeight(branch);
			FaultSystemSolution sol = slt.forBranch(branch);
			IncrementalMagFreqDist fullMFD = sol.calcTotalNucleationMFD(
					refIncrFunc.getMinX(), refIncrFunc.getMaxX(), refIncrFunc.getDelta());
			GridSourceProvider prov = sol.getGridSourceProvider();
			for (int i=0; i<prov.getNumLocations(); i++) {
				IncrementalMagFreqDist mfd = prov.getMFD(i);
				if (mfd == null)
					continue;
				int indexOffset = refIncrFunc.getClosestXIndex(mfd.getMinX());
				Preconditions.checkState(mfd.getMinX() >= refIncrFunc.getMinX());
				for (int j=0; j<mfd.size(); j++) {
					double y = mfd.getY(j);
					if (y > 0d)
						fullMFD.add(j+indexOffset, y);
				}
			}
			
			EvenlyDiscretizedFunc cmlMFD = fullMFD.getCumRateDistWithOffset();
			if (refCmlFunc == null)
				refCmlFunc = cmlMFD;
			
			for (int d=0; d<durations.length; d++)
				for (int i=0; i<refIncrFunc.size(); i++)
					dists[d][i].set(rateToProb(cmlMFD.getY(i), durations[d]), weight);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		List<EvenlyDiscretizedFunc> durMeans = new ArrayList<>();
		List<EvenlyDiscretizedFunc[]> durPercentiles = new ArrayList<>();
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		for (int d=0; d<durations.length; d++) {
			double duration = durations[d];
			Color color = colors[d];
			
			Color transColor = new Color(color.getRed(), color.getGreen(), color.getBlue(), 40);
			
			EvenlyDiscretizedFunc mean = new EvenlyDiscretizedFunc(refCmlFunc.getMinX(), refCmlFunc.size(), refCmlFunc.getDelta());
			EvenlyDiscretizedFunc median = new EvenlyDiscretizedFunc(refCmlFunc.getMinX(), refCmlFunc.size(), refCmlFunc.getDelta());
			EvenlyDiscretizedFunc[] percentileFuncs = new EvenlyDiscretizedFunc[percentiles.length];
			for (int p=0; p<percentileFuncs.length; p++)
				percentileFuncs[p] = new EvenlyDiscretizedFunc(refCmlFunc.getMinX(), refCmlFunc.size(), refCmlFunc.getDelta());
			
			for (int i=0; i<refCmlFunc.size(); i++) {
				ArbDiscrEmpiricalDistFunc dist = dists[d][i];
				
				mean.set(i, dist.getMean());
				median.set(i, dist.getMedian());
				for (int p=0; p<percentiles.length; p++) {
					double val;
					if (percentiles[p] == 0d)
						val = dist.getMinX();
					else if (percentiles[p] == 1d)
						val = dist.getMaxX();
					else
						val = dist.getInterpolatedFractile(percentiles[p]);
					percentileFuncs[p].set(i, val);
				}
			}
			
			durMeans.add(mean);
			durPercentiles.add(percentileFuncs);
			
			mean.setName(oDF.format(duration)+"-year");
			funcs.add(mean);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
			
			Preconditions.checkState(percentiles.length % 2 == 0);
			for (int i=0; i<percentiles.length/2; i++) {
				EvenlyDiscretizedFunc lower = percentileFuncs[i];
				EvenlyDiscretizedFunc upper = percentileFuncs[percentileFuncs.length - (1+i)];
				UncertainArbDiscFunc uncert = new UncertainArbDiscFunc(median, lower, upper);
				funcs.add(uncert);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, transColor));
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Magnitude", "Probability of Exceedance");
		spec.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		gp.getPlotPrefs().scaleFontSizes(1.3d);
		
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, false, false, xRange, new Range(0d, 1d));
		
		PlotUtils.setXTick(gp, 0.5);
		PlotUtils.setYTick(gp, 0.1);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 800, true, true, false);
		
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		header.add("Magnitude");
		for (double duration : durations) {
			header.add(oDF.format(duration)+"-year Mean");
			for (double percentile : percentiles) {
				if (percentile == 0d)
					header.add("Minimum");
				else if (percentile == 1d)
					header.add("Maximum");
				else
					header.add("p"+oDF.format(percentile*100d));
			}
		}
		csv.addLine(header);
		
		for (int i=0; i<refCmlFunc.size(); i++) {
			double mag = refCmlFunc.getX(i);
			if ((float)mag < (float)xRange.getLowerBound())
				continue;
			List<String> line = new ArrayList<>(header.size());
			line.add((float)mag+"");
			for (int d=0; d<durations.length; d++) {
				line.add((float)durMeans.get(d).getY(i)+"");
				for (EvenlyDiscretizedFunc func : durPercentiles.get(d))
					line.add((float)func.getY(i)+"");
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, prefix+".csv"));
	}
	
	private static double rateToProb(double rate, double duration) {
		return 1d - Math.exp(-rate*duration);
	}

}
