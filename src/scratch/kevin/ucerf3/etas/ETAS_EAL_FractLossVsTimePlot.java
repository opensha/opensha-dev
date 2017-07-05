package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.eal.UCERF3_BranchAvgLossFetcher;

public class ETAS_EAL_FractLossVsTimePlot {

	public static void main(String[] args) throws IOException, DocumentException {
		File trueMeanSolFile = new File("dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		// Branch averaged FSS
		FaultSystemSolution baSol = FaultSystemIO.loadSol(
				new File("dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultModels fm = FaultModels.FM3_1;
		
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		File ealMainDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		
		File[] ealDataDirs = { new File(ealMainDir, "2014_05_28-ucerf3-99percent-wills-smaller"),
				new File(ealMainDir, "2016_06_06-ucerf3-90percent-wald")};
		
		File etasCatalogsFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_06_15-haywired_m7-10yr-full_td-no_ert-combined/results_descendents_m5.bin");
		
		File outputDir = new File(etasCatalogsFile.getParentFile(), "loss_results");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		if (ealDataDirs.length == 1)
			outputDir = new File(outputDir, ealDataDirs[0].getName());
		else
			outputDir = new File(outputDir, "combined");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		// IMR for which EAL data has already been computed
		Map<AttenRelRef, Double> imrWeightsMap = Maps.newHashMap();
//		imrWeightsMap.put(AttenRelRef.BSSA_2014, 1d);
		
		imrWeightsMap.put(AttenRelRef.CB_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.CY_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.ASK_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.BSSA_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.IDRISS_2014, 0.12);
		
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(etasCatalogsFile, 5d);
		
		double minTimeYears = 1d/365.25/24d/60d/60d;
		double maxTimeYears = 10d;
		int numTimeBins = 500;
		
		String valueLabel = "$ (Billions)";
		// portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		double thousandsToBillions = 1d/1e6;
		double inflationScalar = 1d/0.9d;
		double valueScale = thousandsToBillions*inflationScalar;
		
		EvenlyDiscretizedFunc logTimeFunc = new EvenlyDiscretizedFunc(Math.log10(minTimeYears),
				Math.log10(maxTimeYears), numTimeBins);
		EvenlyDiscretizedFunc cumLogTimeFunc = new EvenlyDiscretizedFunc(Math.log10(minTimeYears),
				Math.log10(maxTimeYears), numTimeBins);
		
		long minOT = Long.MAX_VALUE;
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			if (!catalog.isEmpty() && catalog.get(0).getOriginTime() < minOT)
				minOT = catalog.get(0).getOriginTime();
		}
		
		double[] meanFaultLosses = null;
		DiscretizedFunc[] meanGriddedLosses = null;
		GridSourceProvider gridProv = baSol.getGridSourceProvider();
		GriddedRegion region = gridProv.getGriddedRegion();
		double ealWeightEach = 1d/(double)ealDataDirs.length;
		for (AttenRelRef ref : imrWeightsMap.keySet()) {
			double gmpeWeight = imrWeightsMap.get(ref);
			for (File ealDataDir : ealDataDirs) {
				double weight = gmpeWeight*ealWeightEach;
				UCERF3_BranchAvgLossFetcher lossFetch =
						new UCERF3_BranchAvgLossFetcher(trueMeanSolFile, cfss, ealDataDir);
				DiscretizedFunc[] myFaultLosses = lossFetch.getFaultLosses(ref, fm, true);
				DiscretizedFunc[] myGriddedLosses = lossFetch.getGriddedMagLossDists(ref, region);
				
				if (meanFaultLosses == null) {
					meanFaultLosses = new double[myFaultLosses.length];
					meanGriddedLosses = new DiscretizedFunc[myGriddedLosses.length];
					for (int i=0; i<meanGriddedLosses.length; i++)
						meanGriddedLosses[i] = new ArbitrarilyDiscretizedFunc();
				}
				
				for (int i=0; i<meanFaultLosses.length; i++) {
					DiscretizedFunc condLossFunc = myFaultLosses[i];
					if (condLossFunc == null || condLossFunc.size() == 0)
						continue;
					double sumY = condLossFunc.calcSumOfY_Vals();
					Preconditions.checkState((float)sumY == 1f, "Cond losses don't sum to 1: %s", sumY);
					double meanLoss = 0d;
					for (Point2D pt : condLossFunc)
						meanLoss += pt.getX()*pt.getY();
					meanFaultLosses[i] += meanLoss*weight;
				}
				
				for (int i=0; i<meanGriddedLosses.length; i++) {
					DiscretizedFunc magLossFunc = myGriddedLosses[i];
					for (Point2D pt : magLossFunc) {
						int magIndex = UCERF3_BranchAvgLossFetcher.getMatchingXIndexFloatPrecision(pt.getX(), meanGriddedLosses[i]);
						if (magIndex >= 0)
							meanGriddedLosses[i].set(magIndex, meanGriddedLosses[i].getY(magIndex) + pt.getY()*weight);
						else
							meanGriddedLosses[i].set(pt.getX(), pt.getY()*weight);
					}
				}
			}
		}
		
//		ETAS_CatalogEALCalculator etasCalc = new ETAS_CatalogEALCalculator(lossFetch, baSol, fm, catalogs);
		
		System.out.println("Calculating fractional losses");
		double totLoss = 0d;
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			for (ETAS_EqkRupture rup : catalog) {
				double rupLoss;
				if (rup.getFSSIndex() >= 0) {
					rupLoss = meanFaultLosses[rup.getFSSIndex()];
				} else {
					int nodeIndex = ETAS_CatalogEALCalculator.calcNodeIndex(rup, region);
					rupLoss = meanGriddedLosses[nodeIndex].getClosestYtoX(rup.getMag());
				}
				rupLoss *= valueScale; // convert to billions
				long deltaOT = rup.getOriginTime() - minOT;
				double deltaYears = (double)deltaOT / ProbabilityModelsCalc.MILLISEC_PER_YEAR;
				int timeIndex;
				if (deltaYears == 0)
					timeIndex = 0;
				else
					timeIndex = logTimeFunc.getClosestXIndex(Math.log10(deltaYears));
				// now incremental
				logTimeFunc.add(timeIndex, rupLoss);
				for (int i=timeIndex; i<logTimeFunc.size(); i++)
					cumLogTimeFunc.add(i, rupLoss);
				totLoss += rupLoss;
			}
		}
		System.out.println("Preparing plot");
		logTimeFunc.scale(1d/catalogs.size());
		ArbitrarilyDiscretizedFunc linearTimeFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc cumLinearTimeFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<logTimeFunc.size(); i++) {
			double x = Math.pow(10, logTimeFunc.getX(i));
			linearTimeFunc.set(x, logTimeFunc.getY(i));
			cumLinearTimeFunc.set(x, cumLogTimeFunc.getY(i));
		}
		
		System.out.println("1 year val: "+linearTimeFunc.getInterpolatedY(1d));
		System.out.println("10 year val: "+linearTimeFunc.getInterpolatedY(10d));
		
		// now normalize to one
		ArbitrarilyDiscretizedFunc normalizedCumTimeFunc = cumLinearTimeFunc.deepClone();
		normalizedCumTimeFunc.scale(1d/normalizedCumTimeFunc.getMaxY());
		
		// and build density for later
		ArbitrarilyDiscretizedFunc densityFunc = linearTimeFunc.deepClone();
		double logDelta = logTimeFunc.getDelta();
		for (int i=0; i<densityFunc.size(); i++) {
			double logX = logTimeFunc.getX(i);
			double minLogX = logX - logDelta;
			double maxLogX = logX + logDelta;
			double width = Math.pow(10, maxLogX) - Math.pow(10, minLogX);
			densityFunc.set(i, densityFunc.getY(i)/width);
		}
		
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(normalizedCumTimeFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		List<XYTextAnnotation> annotations = addTimeLables(0d, 1d, false, funcs, chars);
		
		PlotSpec spec = new PlotSpec(funcs, chars, "", "Time (years)", "Fraction of Losses Before");
		spec.setPlotAnnotations(annotations);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setUserBounds(minTimeYears, maxTimeYears, 0d, 1d);
		
		ETAS_MultiSimAnalysisTools.setFontSizes(gp);
		
		File file = new File(outputDir, "fract_losses_before");
		
		gp.drawGraphPanel(spec, true, false);
		gp.getChartPanel().setSize(1000, 800);
		System.out.println("Writing TXT");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
		System.out.println("Writing PNG");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		System.out.println("Writing PDF");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		funcs.clear();
		chars.clear();;
		
		// now plot density
		funcs.add(densityFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		// add long term line
		XY_DataSet longTerm = new DefaultXY_DataSet();
		longTerm.set(minTimeYears, 4.09975);
		longTerm.set(maxTimeYears, 4.09975);
		funcs.add(longTerm);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		
		double minY = Math.min(1d, densityFunc.getMinY()*0.9);
		double maxY = densityFunc.getMaxY()*1.2;
		
		annotations = addTimeLables(minY, maxY, true, funcs, chars);
		
		spec = new PlotSpec(funcs, chars, "", "Time (years)", "Average Loss Rate Density");
		spec.setPlotAnnotations(annotations);
		
		gp = new HeadlessGraphPanel();
		gp.setUserBounds(minTimeYears, maxTimeYears, minY, maxY);
		
		ETAS_MultiSimAnalysisTools.setFontSizes(gp);
		
		file = new File(outputDir, "loss_rate_density");
		
		gp.drawGraphPanel(spec, true, true);
		gp.getChartPanel().setSize(1000, 800);
		System.out.println("Writing TXT");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
		System.out.println("Writing PNG");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		System.out.println("Writing PDF");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private static List<XYTextAnnotation> addTimeLables(double minY, double maxY, boolean logY,
			List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars) {
List<XYTextAnnotation> annotations = Lists.newArrayList();
		
		List<Double> times = Lists.newArrayList();
		List<String> timeLabels = Lists.newArrayList();
		
		times.add(1d/365.25/24d/60d/60d);
		timeLabels.add("1 s");
		
		times.add(1d/365.25/24d/60d);
		timeLabels.add("1 m");
		
		times.add(1d/365.25/24d);
		timeLabels.add("1 hr");
		
		times.add(1d/365.25);
		timeLabels.add("1 d");
		
		times.add(7d/365.25);
		timeLabels.add("1 wk");
		
		times.add(30/365.25);
		timeLabels.add("1 mo");
		
		times.add(1d);
		timeLabels.add("1 yr");
		
		times.add(10d);
		timeLabels.add("10 yr");
		
		double mainAnnY, lastAnnY;
		if (logY) {
			double deltaY = Math.log10(maxY) - Math.log10(minY);
			mainAnnY = Math.pow(10, Math.log10(minY) + deltaY*0.95);
			lastAnnY = Math.pow(10, Math.log10(minY) + deltaY*0.9);
		} else {
			double deltaY = maxY - minY;
			mainAnnY = minY + deltaY*0.95;
			lastAnnY = minY + deltaY*0.9;
		}
		
		for (int i=0; i<times.size(); i++) {
			double time = times.get(i);
			String label = timeLabels.get(i);
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			xy.set(time, minY);
			xy.set(time, maxY);
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			
			XYTextAnnotation ann = new XYTextAnnotation(label, time, mainAnnY);
			if (i == (times.size()-1)) {
				ann.setTextAnchor(TextAnchor.TOP_RIGHT);
				ann.setY(lastAnnY); // put it below
			} else {
				ann.setTextAnchor(TextAnchor.TOP_LEFT);
			}
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
			annotations.add(ann);
		}
		
		return annotations;
	}

}
