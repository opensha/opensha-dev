package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FSSRupsInRegionCache;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.etas.ConditionalSectTriggerMPDCalc;

public class TomRanchoMirageCalcs {
	
	private static final double minX = 5.05d;
	private static final double maxX = 9.05d;
	private static final double delta = 0.1d;
	private static final int num = (int) ((maxX - minX) / delta + 1);

	public static void main(String[] args) throws IOException, DocumentException {
		// use this one because it has rupture MFDs, so calculations will include mag uncertainty
//		File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/UCERF3_ERF/"
//				+ "cached_FM3_1_dep100.0_depMean_rakeMean.zip");
		File solFileForTD = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/UCERF3_ERF/"
				+ "cached_FM3_1_dep100.0_depMean_rakeMean.zip"); // has rup mag prob dists
		File solFileForETAS = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"); // correct IDs for ETAS sims
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFileForTD);
		
//		File etasFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined/results_descendents_m5_preserve.bin");
//		String etasScenarioLabel = "mojave_m7";
//		String etasScenarioTitle = "Following Mojave M7";
		File etasFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_10_27-2016_bombay_swarm-10yr-full_td-no_ert-combined/results.bin");
		String etasScenarioLabel = "2016_bb_swarm";
		String etasScenarioTitle = "Following 2016 Bombay Beach Swarm";
		double etasDuration = 7d/365.25;
		String etasDurationLabel = "1 week";
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/tom_rancho_mirage");
		
		int coachellaID = 295;
		Region soCalRegion = new CaliforniaRegions.RELM_SOCAL();
		
//		double[] durations = { 1d, 10d, 30d };
		double[] durations = {};
		int startYear = 2017;
		
		FaultSystemSolutionERF tdERF = new FaultSystemSolutionERF(sol);
		tdERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		tdERF.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
		tdERF.setParameter(HistoricOpenIntervalParam.NAME, startYear-1875d);
		tdERF.getTimeSpan().setStartTime(startYear);
		
		FaultSystemSolutionERF tiERF = new FaultSystemSolutionERF(sol);
		tiERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		tiERF.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
		
		FSSRupsInRegionCache rupsCache = new FSSRupsInRegionCache();
		
		for (double duration : durations) {
			System.out.println("Calculating "+duration+" yr");
			tdERF.getTimeSpan().setDuration(duration);
			tiERF.getTimeSpan().setDuration(duration);
			
			tdERF.updateForecast();
			tiERF.updateForecast();
			
			// these are cumulative
			EvenlyDiscretizedFunc tdFaultMPD = calcParentSectMPD(tdERF, coachellaID);
			EvenlyDiscretizedFunc tiFaultMPD = calcParentSectMPD(tiERF, coachellaID);
			
			SummedMagFreqDist tdRegionMFD = ERF_Calculator
					.getParticipationMagFreqDistInRegion(tdERF, soCalRegion, minX,
							num, delta, true, rupsCache);
			// to cumulative probs
			EvenlyDiscretizedFunc tdRegionMPD = FaultSysSolutionERF_Calc.calcProbsFromSummedMFD(
					tdRegionMFD.getCumRateDistWithOffset(), duration);
			SummedMagFreqDist tiRegionMFD = ERF_Calculator
					.getParticipationMagFreqDistInRegion(tiERF, soCalRegion, minX,
							num, delta, true, rupsCache);
			// to cumulative probs
			EvenlyDiscretizedFunc tiRegionMPD = FaultSysSolutionERF_Calc.calcProbsFromSummedMFD(
					tiRegionMFD.getCumRateDistWithOffset(), duration);
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(tiFaultMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			tiFaultMPD.setName("Coachella-TI");
			
			funcs.add(tdFaultMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			tdFaultMPD.setName("Coachella-TD");
			
			funcs.add(tiRegionMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
			tiRegionMPD.setName("SoCal-TI");
			
			funcs.add(tdRegionMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
			tdRegionMPD.setName("SoCal-TD");
			
			PlotSpec spec = new PlotSpec(funcs, chars, (int)duration+" yr Magnitude-Probability Distributions",
					"Magnitude", (int)duration+" yr Probability");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setLegendFontSize(18);
			gp.setBackgroundColor(Color.WHITE);
			
//			gp.setUserBounds(tdFaultMPD.getMinX(), tdFaultMPD.getMaxX(), 0d, 1d);
			gp.setUserBounds(6.5, 8.5, 0d, 1d);
			gp.drawGraphPanel(spec);
			gp.getChartPanel().setSize(800, 600);
			
			String prefix = "mpds_"+(int)duration+"yr";
			
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
		}
		
		if (etasFile != null) {
			System.out.println("Loading ETAS catalogs");
			List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(etasFile);
			long maxOT = ConditionalSectTriggerMPDCalc.calcMaxOTforDuration(catalogs, etasDuration);
			
			System.out.println("Calculating ETAS");
			FaultSystemSolution solForETAS = FaultSystemIO.loadSol(solFileForETAS);
			EvenlyDiscretizedFunc etasMPD = ConditionalSectTriggerMPDCalc.calcCumulativeMPD(
					catalogs, solForETAS, coachellaID, 6.5, 8.5, 0.1, maxOT);
			System.out.println("Max OT: "+maxOT);
			
			System.out.println("Calculating long term comparisons");
			tdERF.getTimeSpan().setDuration(etasDuration);
			tiERF.getTimeSpan().setDuration(etasDuration);
			
			tdERF.updateForecast();
			tiERF.updateForecast();
			
			// these are cumulative
			EvenlyDiscretizedFunc tdMPD = calcParentSectMPD(tdERF, coachellaID);
			EvenlyDiscretizedFunc tiMPD = calcParentSectMPD(tiERF, coachellaID);
			
			System.out.println("TD:\n"+tdMPD);
			System.out.println("ETAS:\n"+etasMPD);
			
			EvenlyDiscretizedFunc etasCombMPD = new EvenlyDiscretizedFunc(etasMPD.getMinX(), etasMPD.getMaxX(), etasMPD.size());
			for (int i=0; i<etasCombMPD.size(); i++) {
				double x = etasMPD.getX(i);
				double etasProb = etasMPD.getY(i);
				int compIndex = tdMPD.getClosestXIndex(x);
				Preconditions.checkState(compIndex >= 0, "No TD index for x="+x);
				double tdProb = tdMPD.getY(compIndex);
				
				double combProb = 1d - (1d - etasProb)*(1d - tdProb);
				etasCombMPD.set(i, combProb);
			}
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(tiMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			tiMPD.setName("Coachella-TI");
			
			funcs.add(tdMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			tdMPD.setName("Coachella-TD");
			
			funcs.add(etasCombMPD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			etasCombMPD.setName("Coachella-ETAS");
			
			PlotSpec spec = new PlotSpec(funcs, chars, etasDurationLabel+" "+etasScenarioTitle,
					"Magnitude", etasDurationLabel+" Probability");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setLegendFontSize(18);
			gp.setBackgroundColor(Color.WHITE);
			
//			gp.setUserBounds(tdFaultMPD.getMinX(), tdFaultMPD.getMaxX(), 0d, 1d);
			gp.setUserBounds(6.5, 8.25, 1e-7, 1e-1);
//			gp.drawGraphPanel(spec);
			gp.drawGraphPanel(spec, false, true);
			gp.getChartPanel().setSize(800, 600);
			
			String prefix = "mpds_etas_"+etasScenarioLabel+"_"+etasDurationLabel.replaceAll(" ", "_");
			
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
			
			System.out.println("ETAS "+etasDurationLabel+" M7 Probabilities "+etasScenarioTitle+":");
			System.out.println("\tTI Prob: "+(float)tiMPD.getY(7d));
			System.out.println("\tTD Prob: "+(float)tdMPD.getY(7d));
			System.out.println("\tETAS Prob: "+(float)etasMPD.getY(7d));
			System.out.println("\tETAS+TD Prob: "+(float)etasCombMPD.getY(7d));
		}
	}
	
	private static EvenlyDiscretizedFunc calcParentSectMPD(FaultSystemSolutionERF erf, int parentID) {
		EvenlyDiscretizedFunc probFunc = new EvenlyDiscretizedFunc(minX-delta*0.5, num, delta);
		
		FaultSystemRupSet rupSet = erf.getSolution().getRupSet();
		
		HashSet<Integer> rupIndexes = new HashSet<Integer>(rupSet.getRupturesForParentSection(parentID));
		
		List<List<Double>> probs = Lists.newArrayList();
		for (int i=0; i<probFunc.size(); i++)
			probs.add(new ArrayList<Double>());
		
		for (int sourceID=0; sourceID<erf.getNumFaultSystemSources(); sourceID++) {
			int fssIndex = erf.getFltSysRupIndexForSource(sourceID);
			if (!rupIndexes.contains(fssIndex))
				continue;
			
			ProbEqkSource source = erf.getSource(sourceID);
			
			for (ProbEqkRupture rup : source) {
				double mag = rup.getMag();
				double prob = rup.getProbability();
				if (prob == 0)
					continue;
				for (int i=0; i<probFunc.size(); i++) {
					if (mag >= probFunc.getX(i)) {
						probs.get(i).add(prob);
					}
				}
			}
		}
		
		for (int i=0; i<probFunc.size(); i++) 
			probFunc.set(i, FaultSysSolutionERF_Calc.calcSummedProbs(probs.get(i)));
		
		return probFunc;
	}

}
