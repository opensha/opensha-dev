package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.TimeZone;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;

public class LongTermDebug {

	public static void main(String[] args) throws IOException {
		File simsDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/");
		
		List<String> names = new ArrayList<>();
		List<File> binFiles = new ArrayList<>();
		List<Long> startTimes = new ArrayList<>();
		List<Double> durations = new ArrayList<>();
		List<Color> colors = new ArrayList<>();
		
		double[] plotDurations = { 50d, 100d, 500d };
		
		double minMag = 5d;
		double[] deltas = { 1d, 5d, 10d };
		File outputDir = new File("/tmp/etas_debug");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		names.add("BSSA Paper (2016)");
		binFiles.add(new File(simsDir, "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m5_preserve.bin"));
		startTimes.add(1325419200000l);
		durations.add(1000d);
		colors.add(Color.BLACK);
		
		GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("UTC"));
		cal.clear();
		cal.set(2012, 0, 1);
		long newStart = cal.getTimeInMillis();
		
		names.add("2019 (RELAXED)");
		binFiles.add(new File(simsDir, "2019_10_11-Start2012_500yr_Spontaneous_HistoricalCatalog/results_m5_preserve_chain.bin"));
		startTimes.add(newStart);
		durations.add(500d);
		colors.add(Color.BLUE);
		
		names.add("2021 (STRICT)");
		binFiles.add(new File(simsDir, "2021_03_02-Start2012_500yr_CompletenessSTRICT_Spontaneous_HistoricalCatalog/results_m5_preserve_chain.bin"));
		startTimes.add(newStart);
		durations.add(500d);
		colors.add(Color.RED);
		
//		names.add("Nov 2018 JAR (RELAXED)");
//		binFiles.add(new File(simsDir, "2021_03_03-Start2012_500yr_CompletenessSTRICT_Spontaneous_HistoricalCatalog_OLDJAR/results_m5_preserve_chain.bin"));
//		startTimes.add(newStart);
//		durations.add(500d);
//		colors.add(Color.GREEN);
		
//		names.add("Aug 31 2018 Jar");
//		binFiles.add(new File(simsDir, "2021_03_03-Start2012_500yr_CompletenessSTRICT_Spontaneous_HistoricalCatalog_OLDJAR2/results_m5_preserve_chain.bin"));
//		startTimes.add(newStart);
//		durations.add(500d);
//		colors.add(Color.ORANGE);
		
		names.add("Aug 15 2018 Jar");
		binFiles.add(new File(simsDir, "2021_03_03-Start2012_500yr_CompletenessSTRICT_Spontaneous_HistoricalCatalog_OLDJAR3/results_m5_preserve_chain.bin"));
		startTimes.add(newStart);
		durations.add(500d);
		colors.add(Color.MAGENTA);
		
		names.add("Old Launcher Current Jar");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14/results_m4_preserve.bin"));
		startTimes.add(newStart);
		durations.add(500d);
		colors.add(Color.CYAN);
		
		names.add("2016 REDO w/ 1 thread");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-REDO-2016-JAR/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(500d);
		colors.add(Color.GREEN);
		
//		names.add("2016 REDO w/ 2 threads");
//		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-REDO-2016-JAR-2thread/results_m5_ip.bin"));
//		startTimes.add(newStart);
//		durations.add(500d);
//		colors.add(Color.YELLOW.darker());
		
		names.add("2016 REDO-Scale w/ 1 thread");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-REDO-2016-SCALE-JAR/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(500d);
		colors.add(Color.GRAY);
		
		names.add("2016-09 SVN");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-2016-09-SCALE-JAR-2thread/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(200d);
		colors.add(Color.ORANGE);
		
		names.add("2017-03 SVN");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-2017-03-SCALE-JAR-2thread/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(200d);
		colors.add(Color.YELLOW.darker());
		
		names.add("2017-09 END SVN");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-2017-09-SCALE-JAR-2thread/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(200d);
		colors.add(Color.GREEN.darker());
		
		names.add("2018-03 GIT");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-2018-03-GIT-JAR-2thread/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(200d);
		colors.add(Color.RED.darker());
		
		names.add("2017-09 GIT");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-2017-09-GIT-JAR-2thread/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(200d);
		colors.add(Color.BLUE.darker());
		
		names.add("2017-12 GIT");
		binFiles.add(new File(simsDir, "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-2017-12-GIT-JAR-2thread/results_m5.bin"));
		startTimes.add(newStart);
		durations.add(200d);
		colors.add(Color.MAGENTA.darker());
		
		List<List<EvenlyDiscretizedFunc>> spontFuncs = new ArrayList<>();
		List<List<EvenlyDiscretizedFunc>> triggerFuncs = new ArrayList<>();
		List<List<EvenlyDiscretizedFunc>> triggerPrimaryFuncs = new ArrayList<>();
		for (int i=0; i<deltas.length; i++) {
			spontFuncs.add(new ArrayList<>());
			triggerFuncs.add(new ArrayList<>());
			triggerPrimaryFuncs.add(new ArrayList<>());
		}
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<binFiles.size(); i++) {
			List<ETAS_Catalog> catalogs = ETAS_CatalogIO.loadCatalogsBinary(binFiles.get(i), minMag);
			double duration = durations.get(i);
			long startTime = startTimes.get(i);
			
			String name = names.get(i);
			
			double rateEach = 1d/catalogs.size();
			
			long numEvents = 0;
			
			for (int j = 0; j < deltas.length; j++) {
				double delta = deltas[j];
				EvenlyDiscretizedFunc spontFunc = HistogramFunction.getEncompassingHistogram(0d, duration, delta);
				spontFunc.setName(names.get(i));
				EvenlyDiscretizedFunc triggerFunc = new EvenlyDiscretizedFunc(spontFunc.getMinX(), spontFunc.getMaxX(), spontFunc.size());
				triggerFunc.setName(names.get(i));
				EvenlyDiscretizedFunc triggerPrimaryFunc = new EvenlyDiscretizedFunc(spontFunc.getMinX(), spontFunc.getMaxX(), spontFunc.size());
				triggerPrimaryFunc.setName(names.get(i));
				
				for (ETAS_Catalog catalog : catalogs) {
					if (j == 0)
						numEvents += catalog.size();
					for (ETAS_EqkRupture rup : catalog) {
						long time = rup.getOriginTime();
						long tDelta = time - startTime;
						Preconditions.checkState(tDelta >= 0l);
						double relTime = (double)tDelta/ProbabilityModelsCalc.MILLISEC_PER_YEAR;
						Preconditions.checkState((float)relTime <= (float)duration+1f, "%s: relTime=%s, duration=%s", name, relTime, duration);
						if (rup.getParentID() >= 0) {
							triggerFunc.add(triggerFunc.getClosestXIndex(relTime), rateEach);
							if (rup.getGeneration() == 1)
								triggerPrimaryFunc.add(triggerPrimaryFunc.getClosestXIndex(relTime), rateEach);
						} else {
							spontFunc.add(spontFunc.getClosestXIndex(relTime), rateEach);
						}
					}
				}
				
				spontFunc.scale(1d/delta);
				triggerFunc.scale(1d/delta);
				triggerPrimaryFunc.scale(1d/delta);
				spontFuncs.get(j).add(spontFunc);
				triggerFuncs.get(j).add(triggerFunc);
				triggerPrimaryFuncs.get(j).add(triggerPrimaryFunc);
			}
			
			double totalRate = (double)numEvents/(duration*catalogs.size());
			System.out.println(name+" total event rate:\t"+totalRate);
			
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(i)));
		}
		
		for (int i=0; i<deltas.length; i++) {
			PlotSpec spontSpec = new PlotSpec(spontFuncs.get(i), chars, "ETAS Debug", "Time Into Simulation (yr)", "Spontaneous Rate");
			spontSpec.setLegendVisible(true);
			PlotSpec triggerSpec = new PlotSpec(triggerFuncs.get(i), chars, "ETAS Debug", "Time Into Simulation (yr)", "Trigger Rate");
			triggerSpec.setLegendVisible(false);
			PlotSpec triggerPrimarySpec = new PlotSpec(triggerPrimaryFuncs.get(i), chars, "ETAS Debug", "Time Into Simulation (yr)", "Triger Primary Rate");
			triggerPrimarySpec.setLegendVisible(false);
			
			List<EvenlyDiscretizedFunc> combFuncs = new ArrayList<>();
			for (int j=0; j<binFiles.size(); j++) {
				EvenlyDiscretizedFunc spontFunc = spontFuncs.get(i).get(j);
				EvenlyDiscretizedFunc triggerFunc = triggerFuncs.get(i).get(j);
				Preconditions.checkState(spontFunc.size() == triggerFunc.size());
				EvenlyDiscretizedFunc combFunc = new EvenlyDiscretizedFunc(spontFunc.getMinX(), spontFunc.getMaxX(), spontFunc.size());
				for (int k=0; k<combFunc.size(); k++)
					combFunc.set(k, spontFunc.getY(k)+triggerFunc.getY(k));
				combFuncs.add(combFunc);
			}
			PlotSpec combSpec = new PlotSpec(combFuncs, chars, "ETAS Debug", "Time Into Simulation (yr)", "Total Event Rate");
			combSpec.setLegendVisible(false);
			
			List<PlotSpec> specs = new ArrayList<>();
			specs.add(spontSpec);
			specs.add(triggerSpec);
			specs.add(triggerPrimarySpec);
			specs.add(combSpec);
			
			List<Range> yRanges = new ArrayList<>();
			for (PlotSpec spec : specs) {
				MinMaxAveTracker yTrack = new MinMaxAveTracker();
				for (DiscretizedFunc func : spec.getPlotFunctionsOnly()) {
					yTrack.addValue(func.getMinY());
					yTrack.addValue(func.getMaxY());
				}
				yRanges.add(new Range(Math.floor(yTrack.getMin()), Math.ceil(yTrack.getMax())));
			}
			
			for (double plotDur : plotDurations) {
				Range xRange = new Range(0d, plotDur);
				List<Range> xRanges = new ArrayList<>();
				xRanges.add(xRange);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				
				String prefix = "compare_"+(int)plotDur+"yr_bin_"+(int)deltas[i]+"yr";
				
				gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
				gp.getChartPanel().setSize(1000, 1400);
				gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
				gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			}
		}
	}

}
