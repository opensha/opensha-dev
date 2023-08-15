package scratch.kevin.simulators.ruptures.multifaultSeparate;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimEqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.BBP_CatalogSimZipLoader;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

public class MultifaultSeparatePageGen {

	public static void main(String[] args) throws IOException {
		File bbpBaseDir = new File("/data/kevin/bbp/parallel");
		File markdownDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs/");
		
		File bbpDir = new File(bbpBaseDir,
				"2023_07_31-rundir5413_multifault_separate-all-maxDist500-noHF-vmLA_BASIN_500-cs500Sites");
		
		File catalogDir = new File("/data/kevin/simulators/catalogs/bruce/rundir5413_multifault_separate");
		File outputDir = new File(markdownDir, "rundir5413/bbp_LA_BASIN_500/multifault_separate");
		RSQSimCatalog catalog = new RSQSimCatalog(catalogDir, "Test Catalog", FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		
		AttenRelRef gmmRef = AttenRelRef.ASK_2014;
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		CSVFile<String> mappingCSV = CSVFile.readFile(new File(catalogDir, "event_mappings.csv"), false);
		
		VelocityModel vm = VelocityModel.LA_BASIN_500;

		File fullZip = new File(bbpDir, "results_rotD.zip");
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		int[] siteCounts = new int[sites.size()];
		int anySiteCount = 0;
		
		HashBiMap<BBP_Site, Site> sitesBBPtoGMPE = HashBiMap.create();
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site site : sites) {
			Site gmpeSite = site.buildGMPE_Site(vm);
			gmpeSite.setName(RSQSimBBP_Config.siteCleanName(site));
			gmpeSites.add(gmpeSite);
			sitesBBPtoGMPE.put(site, gmpeSite);
		}
		
		double[] periods = { 2d, 3d, 5d, 10d };
		double[] maxMoFracts = { 1d, 0.9, 0.75, 0.5 };

		Table<Quantity, Double, DefaultXY_DataSet[]> periodScatters = HashBasedTable.create();
		Table<Quantity, Double, DefaultXY_DataSet[]> gmmPeriodScatters = HashBasedTable.create();
		Map<Double, DefaultXY_DataSet[]> gmmFullRupPeriodScatters = new HashMap<>();
		
		
		Quantity[] quantities = Quantity.values();
		for (Quantity q : quantities) {
			for (double fract : maxMoFracts) {
				DefaultXY_DataSet[] scatters = new DefaultXY_DataSet[periods.length];
				for (int p=0; p<periods.length; p++)
					scatters[p] = new DefaultXY_DataSet();
				periodScatters.put(q, fract, scatters);
				
				scatters = new DefaultXY_DataSet[periods.length];
				for (int p=0; p<periods.length; p++)
					scatters[p] = new DefaultXY_DataSet();
				gmmPeriodScatters.put(q, fract, scatters);
			}
		}
		for (double fract : maxMoFracts) {
			DefaultXY_DataSet[] scatters = new DefaultXY_DataSet[periods.length];
			for (int p=0; p<periods.length; p++)
				scatters[p] = new DefaultXY_DataSet();
			gmmFullRupPeriodScatters.put(fract, scatters);
		}
		
		List<CSVFile<String>> siteCSVs = new ArrayList<>(sites.size());
		CSVFile<String> allSitesCSV = new CSVFile<>(false);
		List<String> siteHeader = new ArrayList<>();
		siteHeader.add("Site Name");
		siteHeader.add("Event ID");
		siteHeader.add("Num Sub-Events");
		for (int i=0; i<2; i++) {
			String prefix = i == 0 ? "Full" : "Sub-Event "+i;
			siteHeader.add(prefix+" Moment");
			siteHeader.add(prefix+" Magnitude");
			siteHeader.add(prefix+" rRup");
			siteHeader.add(prefix+" rJB");
			for (double period : periods)
				siteHeader.add(prefix+" "+oDF.format(period)+"s RotD50");
		}
		allSitesCSV.addLine(siteHeader);
		for (int s=0; s<sites.size(); s++) {
			CSVFile<String> siteCSV = new CSVFile<>(false);
			siteCSV.addLine(siteHeader);
			siteCSVs.add(siteCSV);
		}
		
		List<RSQSimEvent> events = catalog.loader().load();
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		Map<Integer, RSQSimEqkRupture> eventRupsMap = new HashMap<>();
		for (RSQSimEvent event : events) {
			eventsMap.put(event.getID(), event);
			eventRupsMap.put(event.getID(), catalog.getEqkRupture(event));
		}
		
		System.out.println("Calculating with GMM: "+gmmRef.getShortName());
		ExecutorService exec = Executors.newFixedThreadPool(Integer.min(sites.size(),
				Integer.min(30, Runtime.getRuntime().availableProcessors())));
		List<Future<Map<RSQSimEvent, DiscretizedFunc>>> gmmFutures = new ArrayList<>();
		for (Site site : gmpeSites)
			gmmFutures.add(exec.submit(new SiteGMM_Calc(site, gmmRef, events, eventRupsMap, periods)));
		Table<Site, RSQSimEvent, DiscretizedFunc> gmmEventFuncs = HashBasedTable.create();
		for (int s=0; s<sites.size(); s++) {
			Site site = gmpeSites.get(s);
			try {
				Map<RSQSimEvent, DiscretizedFunc> gms = gmmFutures.get(s).get();
				for (RSQSimEvent event : gms.keySet())
					gmmEventFuncs.put(site, event, gms.get(event));
				System.out.println("\tDone with "+site.getName());
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		exec.shutdown();
		
		BBP_CatalogSimZipLoader bbpLoader = new BBP_CatalogSimZipLoader(fullZip, sites, sitesBBPtoGMPE, eventsMap);
		
		for (int mappingRow=1; mappingRow<mappingCSV.getNumRows(); mappingRow++) {
			int col = 0;
			List<String> line = mappingCSV.getLine(mappingRow);
			int origID = Integer.parseInt(line.get(col++));
			double origMag = Double.parseDouble(line.get(col++));
			int fullID = Integer.parseInt(line.get(col++));
			
			int numParts = Integer.parseInt(line.get(col++));
			List<RSQSimEvent> eventParts = new ArrayList<>(numParts);
			
			for (int i=0; i<numParts; i++) {
				int partID = Integer.parseInt(line.get(col++));
				double partMag = Double.parseDouble(line.get(col++));
				Preconditions.checkState(partMag <= origMag, "partMag=%s but origMag=%s", partMag, origMag);
				RSQSimEvent event = eventsMap.get(partID);
				Preconditions.checkState((float)partMag == (float)event.getMagnitude(),
						"Mag mismatch: %s != %s", (float)partMag, (float)event.getMagnitude());
				eventParts.add(event);
			}
			
			RSQSimEvent fullEvent = eventsMap.get(fullID);
			double fullMoment = moment(fullEvent);
			double maxMoFract = 0d;
			for (RSQSimEvent event : eventParts) {
				double moment = moment(event);
				maxMoFract = Math.max(maxMoFract, moment/fullMoment);
			}
			System.out.println("Processing event "+mappingRow+"/"+mappingCSV.getNumRows()+" with origID="+origID
					+", mag="+origMag+", fullID="+fullID+", and "+numParts+" parts with maxMoFract="+(float)maxMoFract);
			Preconditions.checkState((float)maxMoFract <= 1f);
			
			boolean any = false;
			
			for (int s=0; s<sites.size(); s++) {
				Site site = gmpeSites.get(s);
				BBP_Site bbpSite = sites.get(s);
				if (!bbpLoader.contains(bbpSite, fullID))
					continue;
				any = true;
				siteCounts[s]++;
				DiscretizedFunc fullRD50 = bbpLoader.getRotD50(site, fullEvent, 0);
				List<String> siteLine = new ArrayList<>();
				siteLine.add(site.getName());
				siteLine.add(origID+"");
				siteLine.add(numParts+"");
				
				// binned by period
				List<List<Double>> partGMs = new ArrayList<>();
				for (int p=0; p<periods.length; p++)
					partGMs.add(new ArrayList<>(numParts));
				
				for (int i=-1; i<numParts; i++) {
					RSQSimEvent event = i < 0 ? fullEvent : eventParts.get(i);
					RSQSimEqkRupture rup = eventRupsMap.get(event.getID());
					DiscretizedFunc rd50 = bbpLoader.contains(bbpSite, event.getID()) ? bbpLoader.getRotD50(site, event, 0) : null;
					
					siteLine.add((float)moment(event)+"");
					siteLine.add((float)event.getMagnitude()+"");
					siteLine.add((float)rup.getRuptureSurface().getDistanceRup(bbpSite.getLoc())+"");
					siteLine.add((float)rup.getRuptureSurface().getDistanceJB(bbpSite.getLoc())+"");
					if (rd50 == null) {
						for (int p=0; p<periods.length; p++) {
							siteLine.add("");
							if (i >= 0)
								partGMs.get(p).add(0d);
						}
					} else {
						for (int p=0; p<periods.length; p++) {
							double gm = rd50.getInterpolatedY(periods[p]);
							if (i >= 0)
								partGMs.get(p).add(gm);
							siteLine.add(gm+"");
						}
					}
				}
				
				for (int p=0; p<periods.length; p++) {
					double fullGM = fullRD50.getInterpolatedY(periods[p]);
					List<Double> periodPartGMs = partGMs.get(p);
					for (Quantity q : quantities) {
						double val = q.calc(periodPartGMs);
						for (int f=0; f<maxMoFracts.length; f++)
							if ((float)maxMoFract <= (float)maxMoFracts[f])
								periodScatters.get(q, maxMoFracts[f])[p].set(fullGM, val);
					}
					
					// now with the GMM
					List<Double> gmmPeriodPartGMs = new ArrayList<>(eventParts.size());
					for (RSQSimEvent event : eventParts) {
						double gm = gmmEventFuncs.get(site, event).getY(periods[p]);
						gmmPeriodPartGMs.add(gm);
					}
					for (Quantity q : quantities) {
						double val = q.calc(gmmPeriodPartGMs);
						for (int f=0; f<maxMoFracts.length; f++)
							if ((float)maxMoFract <= (float)maxMoFracts[f])
								gmmPeriodScatters.get(q, maxMoFracts[f])[p].set(fullGM, val);
					}
					double gmmFullVal = gmmEventFuncs.get(site, fullEvent).getY(periods[p]);
					for (int f=0; f<maxMoFracts.length; f++)
						if ((float)maxMoFract <= (float)maxMoFracts[f])
							gmmFullRupPeriodScatters.get(maxMoFracts[f])[p].set(fullGM, gmmFullVal);
				}
				siteCSVs.get(s).addLine(siteLine);
				allSitesCSV.addLine(siteLine);
			}
			if (any)
				anySiteCount++;
		}
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Multifault Ground Motions");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Sites and CSV Files");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("Site Name", "Latitude", "Longitude", "Count", "CSV File");
		
		for (int s=0; s<sites.size(); s++) {
			CSVFile<String> csv = siteCSVs.get(s);
			if (csv.getNumRows() > 1) {
				File file = new File(resourcesDir, sites.get(s).getName()+".csv");
				csv.writeToFile(file);
				BBP_Site site = sites.get(s);
				table.addLine(site.getName(), (float)site.getLoc().lat, (float)site.getLoc().lon, siteCounts[s],
						"["+file.getName()+"]("+resourcesDir.getName()+"/"+file.getName()+")");
			}
		}
		File allCSVFile = new File(resourcesDir, "all_sites.csv");
		table.addLine("__All Sites__", "_N/A_", "_N/A_", anySiteCount,
				"["+allCSVFile.getName()+"]("+resourcesDir.getName()+"/"+allCSVFile.getName()+")");
		allSitesCSV.writeToFile(allCSVFile);
		
		lines.addAll(table.build());
		lines.add("");
		
		for (int i=0; i<maxMoFracts.length; i++) {
			String fractStr, fractPrefix;
			if (maxMoFracts[i] == 1d) {
				fractStr = "All Events";
				fractPrefix = "all";
			} else {
				fractStr = "Max Mo Ratio: "+(float)maxMoFracts[i];
				fractPrefix = "fract"+(float)maxMoFracts[i];
			}
			lines.add("## Scatter Plots, "+fractStr);
			lines.add(topLink); lines.add("");
			
			for (int q=0; q<quantities.length; q++) {
				lines.add("### "+quantities[q].name+" Scatter, "+fractStr);
				lines.add(topLink); lines.add("");

				table = MarkdownUtils.tableBuilder();
				table.addLine("Period", "Linear-Linear", "Log-Log");
				
				DefaultXY_DataSet[] scatters = periodScatters.get(quantities[q], maxMoFracts[i]);
				
				for (int p=-1; p<periods.length; p++) {
					String prefix = "scatter_"+fractPrefix+"_"+quantities[q].name();
					
					System.out.println("Plotting: "+prefix);
					
					DefaultXY_DataSet scatter;
					double period;
					String label;
					if (p < 0) {
						label = "All";
						scatter = union(scatters);
						prefix += "_all";
						period = Double.NaN;
					} else {
						label = oDF.format(periods[p])+"s";
						scatter = scatters[p];
						prefix += "_"+oDF.format(periods[p])+"s";
						period = periods[p];
					}
					
					writeScatter(resourcesDir, prefix, period, scatter, quantities[q].axisLabel, Color.BLACK);
					
					table.addLine("__"+label+"__", "![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+".png)",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+"_log.png)");
				}
				
				lines.addAll(table.build());
				lines.add("");
				
				lines.add("### "+gmmRef.getShortName()+", "+quantities[q].name+" Scatter, "+fractStr);
				lines.add(topLink); lines.add("");
				
				scatters = gmmPeriodScatters.get(quantities[q], maxMoFracts[i]);
				
				for (int p=-1; p<periods.length; p++) {
					String prefix = "scatter_"+gmmRef.name()+"_"+fractPrefix+"_"+quantities[q].name();
					String gmmFullPrefix = "scatter_"+gmmRef.name()+"_"+fractPrefix+"_full";
					
					System.out.println("Plotting: "+prefix);
					
					table = MarkdownUtils.tableBuilder();
					table.addLine("", "Linear-Linear", "Log-Log");
					
					DefaultXY_DataSet scatter;
					double period;
					if (p < 0) {
						lines.add("__All Periods__");
						lines.add("");
						scatter = union(scatters);
						prefix += "_all";
						gmmFullPrefix += "_all";
						period = Double.NaN;
					} else {
						lines.add("__"+oDF.format(periods[p])+"s__");
						lines.add("");
						scatter = scatters[p];
						prefix += "_"+oDF.format(periods[p])+"s";
						gmmFullPrefix += "_"+oDF.format(periods[p])+"s";
						period = periods[p];
					}
					
					writeScatter(resourcesDir, prefix, period, scatter,
							gmmRef.getShortName()+" "+quantities[q].axisLabel, Color.RED.darker());
					writeLogDiffScatter(resourcesDir, prefix+"_diff", period, scatter,
							gmmRef.getShortName()+" "+quantities[q].axisLabel, Color.RED.darker());
					
					table.addLine(gmmRef.getShortName()+" "+quantities[q].name,
							"![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+".png)",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+"_log.png)");
					
					
					if (q == 0) {
						// plot the full GMM rupture scatter (will reuse later)
						DefaultXY_DataSet[] fullScatters = gmmFullRupPeriodScatters.get(maxMoFracts[i]);
						DefaultXY_DataSet fullScatter = p < 0 ? union(fullScatters) : fullScatters[p];
						writeScatter(resourcesDir, gmmFullPrefix, period, fullScatter,
								gmmRef.getShortName()+" Full Rupture", Color.BLUE.darker());
						writeLogDiffScatter(resourcesDir, gmmFullPrefix+"_diff", period, fullScatter,
								gmmRef.getShortName()+" Full Rupture", Color.BLUE.darker());
					}
					table.addLine(gmmRef.getShortName()+" Full Rupture",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+gmmFullPrefix+".png)",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+gmmFullPrefix+"_log.png)");
					table.addLine("Residuals",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+"_diff.png)",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+gmmFullPrefix+"_diff.png)");
					
					lines.addAll(table.build());
					lines.add("");
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static double moment(RSQSimEvent event) {
		double ret = 0d;
		for (EventRecord rec : event)
			ret += rec.getMoment();
		return ret;
	}
	
	private enum Quantity {
		MAX("Maximum", "Max(Sub-Event GMs)") {

			@Override
			public double calc(List<Double> subEventGMs) {
				double max = Double.NEGATIVE_INFINITY;
				for (double gm : subEventGMs)
					max = Math.max(max, gm);
				return max;
			}
			
		},
		SUM_SQ("Sum Of Squares", "sqrt(GM_1^2 + ... + GM_N^2)") {

			@Override
			public double calc(List<Double> subEventGMs) {
				double sumSq = 0d;
				for (double gm : subEventGMs)
					sumSq += gm*gm;
				return Math.sqrt(sumSq);
			}
			
		};
//		SUM("Sum", "GM_1 + ... + GM_N") {
//
//			@Override
//			public double calc(List<Double> subEventGMs) {
//				double sum = 0d;
//				for (double gm : subEventGMs)
//					sum += gm;
//				return sum;
//			}
//			
//		};
		
		private String name;
		private String axisLabel;

		private Quantity(String name, String axisLabel) {
			this.name = name;
			this.axisLabel = axisLabel;
		}
		
		public abstract double calc(List<Double> subEventGMs);
	}
	
	private static DefaultXY_DataSet union(DefaultXY_DataSet[] datas) {
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (DefaultXY_DataSet data : datas)
			for (Point2D pt : data)
				ret.set(pt);
		return ret;
	}
	
	private static void writeScatter(File outputDir, String prefix, double period, DefaultXY_DataSet scatter,
			String yAxisLabel, Color color) throws IOException {
		double min = Math.min(scatter.getMinX(), scatter.getMinY());
		double minNonZero = Double.POSITIVE_INFINITY;
		for (Point2D pt : scatter) {
			if (pt.getX() > 0d)
				minNonZero = Math.min(pt.getX(), minNonZero);
			if (pt.getY() > 0d)
				minNonZero = Math.min(pt.getY(), minNonZero);
		}
		double max = Math.max(scatter.getMaxX(), scatter.getMaxY());
		if (!Double.isFinite(minNonZero)) {
			minNonZero = 1e-2;
			max = 1e1;
		}
		
		String title = Double.isFinite(period) ? oDF.format(period)+"s SA" : "All Periods";

		color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 80);
		
		for (boolean log : new boolean[] {false, true}) {
			Range range;
			if (log)
				range = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))),
						Math.pow(10, Math.ceil(Math.log10(max))));
			else
				range = new Range(0d, max*1.05);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
			oneToOne.set(range.getLowerBound(), range.getLowerBound());
			oneToOne.set(range.getUpperBound(), range.getUpperBound());
			
			funcs.add(scatter);
//			chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, color));
			chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, color));
			
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Full Ground Motion (g)", yAxisLabel+" (g)");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, log, log, range, range);
			
			PlotUtils.writePlots(outputDir, prefix+(log ? "_log" : ""), gp, 800, false, true, false, false);
		}
	}
	
	private static void writeLogDiffScatter(File outputDir, String prefix, double period, DefaultXY_DataSet scatter,
			String title, Color color) throws IOException {
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(-2.99d, 2.99d, 0.1);
		// shift by a half bin so that we have one centered on zero
		hist = new HistogramFunction(hist.getMinX()-0.5*hist.getDelta(), hist.size()+1, hist.getDelta());
		
		for (Point2D pt : scatter) {
			double logRef = Math.log(pt.getX());
			double logComp = Math.log(pt.getY());
			
			double diff = logComp - logRef;
			hist.add(hist.getClosestXIndex(diff), 1);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
		
		Range xRange = new Range(hist.getMinX()-0.5*hist.getDelta(), hist.getMaxX() + 0.5*hist.getDelta());
		Range yRange = new Range(0d, hist.getMaxY()*1.05);
		
		DefaultXY_DataSet vertLine = new DefaultXY_DataSet();
		vertLine.set(0d, yRange.getLowerBound());
		vertLine.set(0d, yRange.getUpperBound());
		
		funcs.add(vertLine);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Log Residual vs Full Simulated GM (g)", "Count");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		title += ", " + (Double.isFinite(period) ? oDF.format(period)+"s SA" : "All Periods");
		
		gp.drawGraphPanel(spec, false, false, xRange, null);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 700, true, true, false);
	}
	
	private static class SiteGMM_Calc implements Callable<Map<RSQSimEvent, DiscretizedFunc>> {
		
		private Site site;
		private AttenRelRef gmmRef;
		private List<RSQSimEvent> events;
		private Map<Integer, RSQSimEqkRupture> eventRupsMap;
		private double[] periods;

		public SiteGMM_Calc(Site site, AttenRelRef gmmRef, List<RSQSimEvent> events,
				Map<Integer, RSQSimEqkRupture> eventRupsMap, double[] periods) {
			this.site = site;
			this.gmmRef = gmmRef;
			this.events = events;
			this.eventRupsMap = eventRupsMap;
			this.periods = periods;

		}

		@Override
		public Map<RSQSimEvent, DiscretizedFunc> call() throws Exception {
			ScalarIMR gmm = gmmRef.get();
			
			gmm.setIntensityMeasure(SA_Param.NAME);
			Map<RSQSimEvent, DiscretizedFunc> ret = new HashMap<>(events.size());
			
			gmm.setSite(site);
			
			for (RSQSimEvent event : events) {
				RSQSimEqkRupture rup = eventRupsMap.get(event.getID());
				gmm.setEqkRupture(rup);
				
				double[] yVals = new double[periods.length];
				
				for (int p=0; p<periods.length; p++) {
					SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), periods[p]);
					
					yVals[p] = Math.exp(gmm.getMean());
				}
				
				ret.put(event, new LightFixedXFunc(periods, yVals));
			}
			return ret;
		}
		
	}

}
