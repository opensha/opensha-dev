package scratch.kevin.ucerf3.eal;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

import scratch.kevin.ucerf3.eal.branches.U3_EAL_LogicTreeBranch;

public class LECCurvePlotter {

	public static void main(String[] args) throws IOException {
		File covFile = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "2020_09_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC-covModel/"
				+ "all_branch_lec_results.csv");
		File meanLossFile = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "2020_04_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC/"
				+ "all_branch_lec_results.csv");
		File outputDir = new File("/tmp");
		String prefix = "lec_compare";
		
		// billions
//		Range xRange = new Range(0d, 200d);
		Range xRange = new Range(1e-3, 200d);
		Range yRange = new Range(1e-7, 1e1);
//		Range yRange = new Range(1e-7, 1e0);
		boolean xLog = true;
		boolean yLog = true;
		
		Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> covLECs = UCERF3_LEC_TreeTrimmer.loadLECs(covFile);
		Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> meanLECs = UCERF3_LEC_TreeTrimmer.loadLECs(meanLossFile);
		
		U3_EAL_LogicTreeBranch singleBranch = covLECs.rowKeySet().iterator().next();
		System.out.println("Single branch: "+singleBranch.buildFileName());
		DiscretizedFunc meanSingleBranchCurve = meanLECs.row(singleBranch).values().iterator().next();
		DiscretizedFunc covSingleBranchCurve = covLECs.row(singleBranch).values().iterator().next();
		LossCOV_Model model = LossCOV_Model.PORTER_POWER_LAW_2020_09_01;
		System.out.println("Plotting single branch convolution");
		testPlotConvolve(meanSingleBranchCurve, covSingleBranchCurve, model, outputDir,
				prefix+"_single_convolve", xLog, yLog, xRange, yRange);
		
		System.out.println("Building fractile calcs...");
		FractileCurveCalculator covFractileCalc = fractileCalc(covLECs);
		FractileCurveCalculator meanFractileCalc = fractileCalc(meanLECs);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		System.out.println("Calculating fractiles/means...");
		DiscretizedFunc upperFunc = toBillions(covFractileCalc.getFractile(0.975));
		DiscretizedFunc lowerFunc = toBillions(covFractileCalc.getFractile(0.025));
		DiscretizedFunc covMean = toBillions(covFractileCalc.getMeanCurve());
		DiscretizedFunc mean = toBillions(meanFractileCalc.getMeanCurve());
		
		System.out.println("Plotting mean convolution");
		testPlotConvolve(fromBillions(mean), fromBillions(covMean), model, outputDir,
				prefix+"_convolve", xLog, yLog, xRange, yRange);
		
		System.out.println("Plotting...");
		
		UncertainArbDiscDataset covUncertain = new UncertainArbDiscDataset(covMean, lowerFunc, upperFunc);
		
		covUncertain.setName("All Branch 95% Range");
		funcs.add(covUncertain);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 0, 255, 80)));
		
		covMean.setName("COV Model Mean LEC");
		funcs.add(covMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		mean.setName("COV=0 Mean LEC");
		funcs.add(mean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Loss Eceeedance Curves",
				"Loss (Billion $)", "Annual Probability of Exceedance");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(22);
		gp.setBackgroundColor(Color.WHITE);
//		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, xLog, yLog, xRange, yRange);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
		
		System.out.println("DONE");
	}
	
	private static FractileCurveCalculator fractileCalc(Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> table) {
		XY_DataSetList functionList = new XY_DataSetList();
		List<Double> relativeWts = new ArrayList<>();
		for (Cell<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> cell : table.cellSet()) {
			functionList.add(cell.getValue());
			relativeWts.add(cell.getColumnKey());
		}
		return new FractileCurveCalculator(functionList, relativeWts);
	}
	
	private static DiscretizedFunc toBillions(XY_DataSet xy) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xy)
			ret.set(pt.getX()*1e-6, pt.getY());
		return ret;
	}
	
	private static DiscretizedFunc fromBillions(XY_DataSet xy) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xy)
			ret.set(pt.getX()*1e6, pt.getY());
		return ret;
	}
	
	private static void testPlotConvolve(DiscretizedFunc meanCurve, DiscretizedFunc covCurve, LossCOV_Model model,
			File outputDir, String prefix, boolean xLog, boolean yLog, Range xRange, Range yRange) throws IOException {
		DiscretizedFunc convolved = new ArbitrarilyDiscretizedFunc();
		for (int i=1; i<meanCurve.size(); i++)
			convolved.set(meanCurve.getX(i), 0d);
		
		for (int i=1; i<meanCurve.size(); i++) {
			double loss0 = meanCurve.getX(i-1);
			double loss1 = meanCurve.getX(i);
			double prob0 = meanCurve.getY(i-1);
			double prob1 = meanCurve.getY(i);
			double probInBin = prob0-prob1;
			if (probInBin == 0d)
				continue;
			double loss = 0.5*(loss0+loss1);
			DiscretizedFunc condExceed = model.calcLossExceedanceProbs(convolved, loss);
			for (int j=0; j<convolved.size(); j++)
				convolved.set(j, convolved.getY(j)+condExceed.getY(j)*probInBin);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		covCurve = toBillions(covCurve);
		covCurve.setName("COV Model LEC");
		funcs.add(covCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		meanCurve = toBillions(meanCurve);
		meanCurve.setName("COV=0 LEC");
		funcs.add(meanCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));

		convolved = toBillions(convolved);
		convolved.setName("COV=0 LEC, Convolved");
		funcs.add(convolved);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Loss Eceeedance Curves",
				"Loss (Billion $)", "Annual Probability of Exceedance");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(22);
		gp.setBackgroundColor(Color.WHITE);
//		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, xLog, yLog, xRange, yRange);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
	}

}
