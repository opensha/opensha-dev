package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.calc.recurInterval.BPT_DistCalc;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

import com.google.common.collect.Lists;

public class PaleoOpenIntervalRIPlot {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution baSol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/"
				+ "scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = baSol.getRupSet();
		
		List<FaultSectionPrefData> subSects = rupSet.getFaultSectionDataList();
		
		ArrayList<PaleoRateConstraint> constraints = UCERF3_PaleoRateConstraintFetcher.getConstraints(
				subSects);
//		PaleoProbabilityModel paleoProbs = UCERF3_PaleoProbabilityModel.load();
		PaleoProbabilityModel paleoProbs = null;
		// filter out constraints without OI data
		for (int i=constraints.size(); --i>=0;) {
			PaleoRateConstraint constr = constraints.get(i);
			if (subSects.get(constr.getSectionIndex()).getDateOfLastEvent() == Long.MIN_VALUE)
				constraints.remove(i);
		}
		System.out.println(constraints.size()+" constraints have open intervals");
		
		double aperiodicity = 0.3;
		double duration = 30d;
		
		DefaultXY_DataSet safXY = new DefaultXY_DataSet();
		DefaultXY_DataSet sjXY = new DefaultXY_DataSet();
		DefaultXY_DataSet otherXY = new DefaultXY_DataSet();
		
		HashSet<Integer> safParents = new HashSet<Integer>();
		HashSet<Integer> sjParents = new HashSet<Integer>();
		Map<String, List<Integer>> namedFaultsMap = FaultModels.FM3_1.getNamedFaultsMapAlt();
		for (String fName : namedFaultsMap.keySet()) {
			if (fName.toLowerCase().contains("andreas"))
				safParents.addAll(namedFaultsMap.get(fName));
			else if (fName.toLowerCase().contains("jacinto"))
				sjParents.addAll(namedFaultsMap.get(fName));
		}
		
//		Date oiRelativeDate = new Date(2014, 0, 1);
//		Date oiRelativeDate = new Date();
		long oiRelativeDate = Math.round((2015.0-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		
//		Range range = new Range(1e1, 1e4);
//		boolean log = true;
		Range range = new Range(0, 400);
		boolean log = false;
		boolean contour = true;
		boolean lighter = true;
		
		Region reg = new CaliforniaRegions.RELM_SOCAL();
		
		List<ConfidenceInterval> intervals = Lists.newArrayList();
		
		for (int i=0; i<constraints.size(); i++) {
			PaleoRateConstraint constr = constraints.get(i);
			if (reg != null && !reg.contains(constr.getPaleoSiteLoction()))
				continue;
			FaultSectionPrefData subSect = subSects.get(constraints.get(i).getSectionIndex());
			int s = subSect.getSectionId();
//			double rate = 0d;
//			for (int r : rupSet.getRupturesForSection(s)) {
//				if (paleoProbs == null)
//					rate += baSol.getRateForRup(r);
//				else
//					rate += paleoProbs.getProbPaleoVisible(rupSet, r, s)*baSol.getRateForRup(r);
//			}
			double rate = constr.getMeanRate();
			
			long oi = oiRelativeDate - subSect.getDateOfLastEvent();
			double oiYears = (double)oi/(double)ProbabilityModelsCalc.MILLISEC_PER_YEAR;
			
			double meanRI = 1d/rate;
			
			if (meanRI > range.getUpperBound() || oiYears > range.getUpperBound())
				continue;
			
//			if (oiYears > meanRI)
//				System.out.println
//			
			if (meanRI < range.getUpperBound() && oiYears < range.getUpperBound())
				System.out.println(constraints.get(i).getPaleoSiteName()+": oi="+oiYears+", meanRI="+meanRI
					+"\t"+subSect.getName());
			
			if (safParents.contains(subSect.getParentSectionId()))
				safXY.set(meanRI, oiYears);
			else if (sjParents.contains(subSect.getParentSectionId()))
				sjXY.set(meanRI, oiYears);
			else
				otherXY.set(meanRI, oiYears);
			
			intervals.add(new ConfidenceInterval(meanRI, oiYears, true,
//					1d/constr.getLower95ConfOfRate(), 1d/constr.getUpper95ConfOfRate()));
					1d/constr.getLower68ConfOfRate(), 1d/constr.getUpper68ConfOfRate()));
			if (subSect.getName().contains("Coachella")) {
				intervals.add(new ConfidenceInterval(meanRI, oiYears, false,
						oiYears-40d, oiYears+40d));
			}
		}
		
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		for (ConfidenceInterval interval : intervals) {
			funcs.addAll(interval.generateFuncs());
			chars.addAll(interval.generateChars());
		}
		
		ArbitrarilyDiscretizedFunc diag = new ArbitrarilyDiscretizedFunc();
		for (double exp = -5; exp<10d; exp++)
			diag.set(Math.pow(10, exp), Math.pow(10, exp));
		funcs.add(diag);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		
		float symbolThick = 15f;
		funcs.add(safXY);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, symbolThick, Color.BLACK));
		funcs.add(sjXY);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, symbolThick, Color.BLACK));
		funcs.add(otherXY);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_DIAMOND, symbolThick*0.7f, Color.BLACK));
		
		String title = "UCERF3 Paleo Sites";
		String xLabel = "Mean Recurrence (yr)";
		String yLabel = "Open Interval (yr)";
		
//		PlotSpec spec = new PlotSpec(funcs, chars, title, xLabel, yLabel);
		
		// calc BPT
		double gridSpacing = 1d;
		double minVal = range.getLowerBound()+0.5*gridSpacing;
		int num = (int)((range.getUpperBound()-range.getLowerBound())/gridSpacing);
		EvenlyDiscrXYZ_DataSet bptXYZ = new EvenlyDiscrXYZ_DataSet(num, num, minVal, minVal, gridSpacing);
		for (int xInd=0; xInd<num; xInd++) {
			for (int yInd=0; yInd<num; yInd++) {
				double mean = bptXYZ.getX(xInd);
				double timeSinceLast = bptXYZ.getY(yInd);
				double bptProb = BPT_DistCalc.getCondProb(mean, aperiodicity, timeSinceLast, duration);
				bptXYZ.set(xInd, yInd, bptProb);
			}
		}
		CPT bptCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		if (contour) {
			CPT newCPT = new CPT();
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0.025, 20, 0.05);
			for (int i=0; i<func.size(); i++) {
				float start = (float)(func.getX(i)-0.025d);
				float end = (float)(func.getX(i)+0.025d);
				Color c = bptCPT.getColor((float)func.getX(i));
				CPTVal val = new CPTVal(start, c, end, c);
				newCPT.add(val);
			}
			bptCPT = newCPT;
		}
		if (lighter) {
			for (CPTVal val : bptCPT) {
				val.minColor = brighter(val.minColor);
				val.maxColor = brighter(val.maxColor);
			}
		}
		String zLabel = (int)duration+"yr BPT Prob";
		
		XYZPlotSpec xyzSpec = new XYZPlotSpec(bptXYZ, bptCPT, title, xLabel, yLabel, zLabel);
		
		xyzSpec.setXYElems(funcs);
		xyzSpec.setXYChars(chars);
		
		String outPrefix = "/tmp/paleo_oi_compare";
		if (paleoProbs != null)
			outPrefix += "_paleo_probs";
		
		PlotPreferences prefs = PlotPreferences.getDefault();
		prefs.setTickLabelFontSize(18);
		prefs.setAxisLabelFontSize(20);
		prefs.setPlotLabelFontSize(21);
		prefs.setBackgroundColor(Color.WHITE);
		
		XYZGraphPanel gp = new XYZGraphPanel(prefs);
		
		gp.drawPlot(xyzSpec, log, log, range, range);
//		gp.setUserBounds(range, range);
//		gp.drawGraphPanel(spec, log, log);
		
		gp.getChartPanel().setSize(600, 600);
		gp.saveAsPDF(outPrefix+".pdf");
		gp.saveAsPNG(outPrefix+".png");
	}
	
	private static Color brighter(Color c) {
		double r = (c.getRed() + 255d)/2d;
		double g = (c.getGreen() + 255d)/2d;
		double b = (c.getBlue() + 255d)/2d;
		return new Color((int)r, (int)g, (int)b);
	}
	
	private static class ConfidenceInterval {
		private double centerX;
		private double centerY;
		private boolean isX;
		private double min;
		private double max;
		
		private static PlotCurveCharacterstics pChar =
				new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY);
		private static List<PlotCurveCharacterstics> chars = Lists.newArrayList(pChar, pChar, pChar);
		
		private static double tickWidth = 20d;
		
		public ConfidenceInterval(double centerX, double centerY, boolean isX, double min, double max) {
			this.centerX = centerX;
			this.centerY = centerY;
			this.isX = isX;
			this.min = min;
			this.max = max;
		}
		
		public List<XY_DataSet> generateFuncs() {
			List<XY_DataSet> funcs = Lists.newArrayList();
			DefaultXY_DataSet interval = new DefaultXY_DataSet();
			DefaultXY_DataSet minTick = new DefaultXY_DataSet();
			DefaultXY_DataSet maxTick = new DefaultXY_DataSet();
			if (isX) {
				interval.set(min, centerY);
				interval.set(max, centerY);
				minTick.set(min, centerY+tickWidth*0.5);
				minTick.set(min, centerY-tickWidth*0.5);
				maxTick.set(max, centerY+tickWidth*0.5);
				maxTick.set(max, centerY-tickWidth*0.5);
			} else {
				interval.set(centerX, min);
				interval.set(centerX, max);
				minTick.set(centerX+tickWidth*0.5, min);
				minTick.set(centerX-tickWidth*0.5, min);
				maxTick.set(centerX+tickWidth*0.5, max);
				maxTick.set(centerX-tickWidth*0.5, max);
			}
			funcs.add(interval);
			funcs.add(minTick);
			funcs.add(maxTick);
			return funcs;
		}
		
		public List<PlotCurveCharacterstics> generateChars() {
			return chars;
		}
	}

}
