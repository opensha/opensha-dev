package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder.ParkfieldSelectionCriteria;

public class ParkfieldRIHists {

	public static void main(String[] args) throws IOException {
		File u3ResultsFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2021_11_30-u3_branches-orig_calcs-5h/results.zip");
		HistogramFunction histU3 = initHist();
		double riU3 = calcBranchParkfieldHist(u3ResultsFile, histU3);
		
		File ingredientsFile = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/results.zip");
		HistogramFunction histIngred = initHist();
		double riIngred = calcBranchParkfieldHist(ingredientsFile, histIngred);
		
		File modelFile = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/results.zip");
		HistogramFunction histModel = initHist();
		double riModel = calcBranchParkfieldHist(modelFile, histModel);
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DecimalFormat riDF = new DecimalFormat("0.0");
		
		histU3.setName("UCERF3, RI="+riDF.format(riU3)+" yrs");
		funcs.add(histU3);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		histIngred.setName("NSHM23 Methodology, RI="+riDF.format(riIngred)+" yrs");
		funcs.add(histIngred);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
		histModel.setName("NSHM23 Model, RI="+riDF.format(riModel)+" yrs");
		funcs.add(histModel);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		
		double maxY = 0.5d;
		for (XY_DataSet func : funcs)
			maxY = Math.max(maxY, Math.ceil(10d*1.2*func.getMaxY())/10d);
		
		Range yRange = new Range(0d, Math.min(1d, maxY));
		Range xRange = new Range(histU3.getMinX(), histU3.getMaxX());
		
		DefaultXY_DataSet dataFunc = new DefaultXY_DataSet();
		double dataRI = 1d/NSHM23_ConstraintBuilder.PARKFIELD_RATE.bestEstimate;
		// same fractional
		double fractRI = dataRI*NSHM23_ConstraintBuilder.PARKFIELD_RATE.getPreferredStdDev()/NSHM23_ConstraintBuilder.PARKFIELD_RATE.bestEstimate;
		dataFunc.set(dataRI, 0d);
		dataFunc.set(dataRI, yRange.getUpperBound());
		dataFunc.setName("Bakun et al. (2005), "+riDF.format(dataRI)+" ± "+riDF.format(fractRI)+" yrs");
		funcs.add(0, dataFunc);
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Parkfield M~6 Recurrence Interval (yrs)", "Weighted-Fraction");
		spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setAxisLabelFontSize(26);
		gp.setTickLabelFontSize(24);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.setXTick(gp, 1d);
		
		PlotUtils.writePlots(outputDir, "parkfield_hist", gp, 800, 650, true, true, false);
	}
	
	private static HistogramFunction initHist() {
		return new HistogramFunction(20.0, 61, 0.25);
	}
	
	private static double calcBranchParkfieldHist(File resultsFile, HistogramFunction hist) throws IOException {
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		Map<String, List<Integer>> parkRupsCache = new HashMap<>();
		
		double rateSum = 0d;
		double weightSum = 0d;
		LogicTree<?> tree = slt.getLogicTree();
//		tree = tree.sample(100, false);
		for (LogicTreeBranch<?> branch : tree) {
			RupSetFaultModel fm = branch.requireValue(RupSetFaultModel.class);
			RupSetDeformationModel dm = branch.requireValue(RupSetDeformationModel.class);
			RupSetScalingRelationship scale = branch.requireValue(RupSetScalingRelationship.class);
			
			String rupsKey = fm.getFilePrefix()+"_"+dm.getFilePrefix()+"_"+scale.getFilePrefix();
			if (!parkRupsCache.containsKey(rupsKey)) {
				ParkfieldSelectionCriteria criteria = NSHM23_InvConfigFactory.getParkfieldSelectionCriteria(fm);
				FaultSystemRupSet rupSet = slt.forBranch(branch).getRupSet();
				List<Integer> rups = NSHM23_ConstraintBuilder.findParkfieldRups(rupSet, criteria);
				parkRupsCache.put(rupsKey, rups);
			}
			
			List<Integer> parkRups = parkRupsCache.get(rupsKey);
			
			double[] rates = slt.loadRatesForBranch(branch);
			
			double parkRate = 0d;
			for (int rupIndex : parkRups)
				parkRate += rates[rupIndex];
			
			double weight = tree.getBranchWeight(branch);
			rateSum += parkRate*weight;
			weightSum += weight;
			
			double parkRI = 1d/parkRate;
			System.out.println(branch+": "+parkRI);
			hist.add(hist.getClosestXIndex(parkRI), weight);
		}
		
		hist.normalizeBySumOfY_Vals();
		
		double avgRate = rateSum/weightSum;
		double avgRI = 1d/avgRate;
		System.out.println(resultsFile.getAbsolutePath()+": aveRI="+avgRI);
		return avgRI;
	}

}
