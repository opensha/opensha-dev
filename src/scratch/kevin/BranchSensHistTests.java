package scratch.kevin;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotSpec;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.analysis.BranchSensitivityHistogram;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class BranchSensHistTests {

	public static void main(String[] args) throws ZipException, IOException {
		BranchSensitivityHistogram hists = new BranchSensitivityHistogram("Ratio");
		
		// use CFSS for branch
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
						"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		APrioriBranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		List<LogicTreeBranch> branches = Lists.newArrayList(cfss.getBranches());
		Collections.sort(branches);
		
		// first create "data"
		double[] vals = new double[branches.size()];
		for (int i=0; i<vals.length; i++)
//			vals[i] = Math.pow(Math.random(), 2d);
			vals[i] = Math.random();
		double mean = StatUtils.mean(vals);
		
		for (int i=0; i<branches.size(); i++) {
			LogicTreeBranch branch = branches.get(i);
			double weight = weightProv.getWeight(branch);
			hists.addValues(branch, vals[i]/mean, weight);
		}
		
		File outputDir = new File("/tmp/hist_test");
		if (!outputDir.exists())
			outputDir.mkdir();
		
		// write the ratio hists
		Map<String, PlotSpec> histPlots = hists.getStackedHistPlots(true, 0d, 21, 0.1);
		for (String categoryName : histPlots.keySet()) {
			PlotSpec spec = histPlots.get(categoryName);
			
			mean = hists.calcMean(categoryName);
			double stdDev = hists.calcStdDev(categoryName);
			
			System.out.println(categoryName+": mean="+mean+", sigma="+stdDev);
			
			XYTextAnnotation ann = new XYTextAnnotation("StdDev="+new DecimalFormat("0.00").format(stdDev), 0.05, 0.95);
			ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
			ann.setTextAnchor(TextAnchor.TOP_LEFT);
			spec.setPlotAnnotations(Lists.newArrayList(ann));
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			CommandLineInversionRunner.setFontSizes(gp);
			
			gp.drawGraphPanel(spec, false, false, new Range(0d, 2d), new Range(0d, 1d));
			
			File file = new File(outputDir, categoryName+"_hists");
			gp.getChartPanel().setSize(1000, 600);
			gp.saveAsPDF(file.getAbsolutePath() + ".pdf");
			gp.saveAsPNG(file.getAbsolutePath() + ".png");
		}
	}
}
