package scratch.kevin.simulators.ruptures.subduction;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Supplier;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenRelSupplier;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_AttenRelSupplier;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimEqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import gov.usgs.earthquake.nshmp.gmm.Gmm;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.MultiRupGMPE_ComparePageGen;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_CatalogSimZipLoader;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.RSQSimGeographicMapMaker;

public class CrustalAndSubductionPageGen {
	
	private enum Type {
		FULL("Full"),
		CRUSTAL("Crustal"),
		SUBDUCTION("Subduction");
		
		private String prefix;

		private Type(String prefix) {
			this.prefix = prefix;
		}
	};
	
	private enum Quantity {
		MAX("Maximum", "Max(Subduction, Crustal)") {

			@Override
			public double calc(double subductionGM, double crustalGM) {
				return Double.max(subductionGM, crustalGM);
			}
			
		},
		CRUSTAL("Crustal-Only GM", "Crustal-Only") {

			@Override
			public double calc(double subductionGM, double crustalGM) {
				return crustalGM;
			}
			
		},
		SUBDUCTION("Subduction-Only GM", "Subduction-Only") {

			@Override
			public double calc(double subductionGM, double crustalGM) {
				return subductionGM;
			}
			
		},
		SUM_SQ("Sum Of Squares", "sqrt(Subduction^2 + Crustal^2)") {

			@Override
			public double calc(double subductionGM, double crustalGM) {
				return Math.sqrt(subductionGM*subductionGM + crustalGM*crustalGM);
			}
			
		},
		SUM("Sum", "Subduction + Crustal") {

			@Override
			public double calc(double subductionGM, double crustalGM) {
				return subductionGM + crustalGM;
			}
			
		};
		
		private String name;
		private String axisLabel;

		private Quantity(String name, String axisLabel) {
			this.name = name;
			this.axisLabel = axisLabel;
		}
		
		public abstract double calc(double subductionGM, double crustalGM);
	}
	
	private static final String DEBUG_SITE_NAME = "WLG";
	private static final int DEBUG_EVENT_ID = 152653;

	public static void main(String[] args) throws IOException {
		RSQSimCatalog fullCatalog = Catalogs.BRUCE_5597.instance();
		RSQSimCatalog crustalCatalog = Catalogs.BRUCE_5597_CRUSTAL.instance();
		RSQSimCatalog subductionCatalog = Catalogs.BRUCE_5597_SUB.instance();
		
		File bbpBaseDir = new File("/data/kevin/bbp/parallel");
		
		File fullBBPdir = new File(bbpBaseDir,
//				"2023_03_30-rundir5566-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
				"2023_11_29-rundir5597-all-m6.5-skipYears2000-noHF-vmCENTRAL_JAPAN-standardSitesNZ-griddedSitesNZ");
		File crustalBBPdir = new File(bbpBaseDir,
//				"2023_06_27-rundir5566_crustal-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
//				"2023_06_27-rundir5566_crustal_corupture-all-m5.0-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
				"2023_11_29-rundir5597_crustal_corupture-all-m6.5-skipYears2000-noHF-vmCENTRAL_JAPAN-standardSitesNZ-griddedSitesNZ");
		File subductionBBPdir = new File(bbpBaseDir,
//				"2023_06_27-rundir5566_subduction-all-m6.5-skipYears5000-maxDist500-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
//				"2023_06_27-rundir5566_subduction_corupture-all-m5.0-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
				"2023_11_29-rundir5597_subduction_corupture-all-m6.5-skipYears2000-noHF-vmCENTRAL_JAPAN-standardSitesNZ-griddedSitesNZ");
		
//		AttenRelSupplier subductionGMM = AttenRelRef.ZHAO_2006;
		AttenRelSupplier subductionGMM = new NSHMP_AttenRelSupplier(
				Gmm.AG_20_GLOBAL_INTERFACE, "AG2020 (Global)");
		AttenRelSupplier crustalGMM = AttenRelRef.ASK_2014;
		
		Quantity[] mapQuantities = {Quantity.SUM_SQ};
		double mapPeriod = 10d;
		int maxMaps = 20;
		
		HashSet<String> highlightSiteNames = new HashSet<>();
		
		highlightSiteNames.add("WLG");
		
		File markdownDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs/");
		
		File markdownCatalogDir = new File(markdownDir, fullCatalog.getCatalogDir().getName());
		Preconditions.checkState(markdownCatalogDir.exists() || markdownCatalogDir.mkdir());
		File outputDir = new File(markdownCatalogDir, "crustal_subduction_gms");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		double minMag = 6.5;

		File fullZip = new File(fullBBPdir, "results_rotD.zip");
		File curstalZip = new File(crustalBBPdir, "results_rotD.zip");
		File subductionZip = new File(subductionBBPdir, "results_rotD.zip");
		
		List<BBP_Site> sites = BBP_Site.readFile(crustalBBPdir);
		
		HashBiMap<BBP_Site, Site> sitesBBPtoGMPE = HashBiMap.create();
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site site : sites) {
			Site gmpeSite = site.buildGMPE_Site(vm);
			gmpeSite.setName(RSQSimBBP_Config.siteCleanName(site));
			gmpeSites.add(gmpeSite);
			sitesBBPtoGMPE.put(site, gmpeSite);
		}
		
		List<RSQSimEvent> events = fullCatalog.loader().minMag(minMag).load();
		
		Map<Integer, RSQSimEvent> crustalMap = crustalCatalog.loader().minMag(minMag).loadMap();
		Map<Integer, RSQSimEvent> subductionMap = subductionCatalog.loader().minMag(minMag).loadMap();
		
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		for (RSQSimEvent event : events)
			if (crustalMap.containsKey(event.getID()) && subductionMap.containsKey(event.getID()))
				eventsMap.put(event.getID(), event);

		BBP_CatalogSimZipLoader fullLoader = new BBP_CatalogSimZipLoader(fullZip, sites, sitesBBPtoGMPE, eventsMap);
		BBP_CatalogSimZipLoader crustalLoader = new BBP_CatalogSimZipLoader(curstalZip, sites, sitesBBPtoGMPE, crustalMap);
		BBP_CatalogSimZipLoader subductionLoader = new BBP_CatalogSimZipLoader(subductionZip, sites, sitesBBPtoGMPE, subductionMap);
		
		double[] periods = { 2d, 3d, 5d, 10d };
//		double[] periods = { 2d, 3d, 5d };
		
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

		MultiRupGMPE_ComparePageGen.checkTransSiteParams(crustalGMM, gmpeSites);
		MultiRupGMPE_ComparePageGen.checkTransSiteParams(subductionGMM, gmpeSites);
		System.out.println("Calculating GMMs:");
		System.out.println("Full rup, crustal GMM");
		Table<Integer, Site, DiscretizedFunc> fullRupCrustalGMM = calcGMM(
				fullCatalog, gmpeSites, fullLoader, periods, crustalGMM, exec);
		System.out.println("Full rup, subduction GMM");
		Table<Integer, Site, DiscretizedFunc> fullRupSubductionGMM = calcGMM(
				fullCatalog, gmpeSites, fullLoader, periods, subductionGMM, exec);
		System.out.println("Partial rup, crustal GMM");
		Table<Integer, Site, DiscretizedFunc> partialRupCrustalGMM = calcGMM(
				crustalCatalog, gmpeSites, crustalLoader, periods, crustalGMM, exec);
		System.out.println("Partial rup, subduction GMM");
		Table<Integer, Site, DiscretizedFunc> partialRupSubductionGMM = calcGMM(
				subductionCatalog, gmpeSites, subductionLoader, periods, subductionGMM, exec);
		System.out.println("Done calculating GMMs");
		
		double[] minFracts = { 0d, 0.2, 0.5 };
		
		Quantity[] quantities = Quantity.values();
		
		// now build combined gmms
		List<Table<Integer, Site, DiscretizedFunc>> combinedGMMs = new ArrayList<>();
		List<Table<Integer, Site, DiscretizedFunc>> combinedSimFuncs = new ArrayList<>();
		for (Quantity q : quantities) {
			// do from GMMs now
			Table<Integer, Site, DiscretizedFunc> table = HashBasedTable.create();
			for (int eventID : partialRupCrustalGMM.rowKeySet()) {
				for (Site site : partialRupCrustalGMM.row(eventID).keySet()) {
					DiscretizedFunc crustalGMs = partialRupCrustalGMM.get(eventID, site);
					DiscretizedFunc subductionGMs = partialRupSubductionGMM.get(eventID, site);
					if (crustalGMs != null && subductionGMs != null) {
						double[] yVals = new double[periods.length];
						for (int p=0; p<periods.length; p++)
							yVals[p] = q.calc(subductionGMs.getY(p), crustalGMs.getY(p));
						table.put(eventID, site, new LightFixedXFunc(periods, yVals));
					}
				}
			}
			combinedGMMs.add(table);
			
			// will fill in from simulations later 
			table = HashBasedTable.create();
			combinedSimFuncs.add(table);
		}
		
		DefaultXY_DataSet[][][] scatters = new DefaultXY_DataSet[minFracts.length][quantities.length][periods.length];
		Table<Integer, Site, double[]> eventSitePeriodFracts = HashBasedTable.create();
		
		for (int i=0; i<minFracts.length; i++)
			for (int q=0; q<quantities.length; q++)
				for (int p=0; p<periods.length; p++)
					scatters[i][q][p] = new DefaultXY_DataSet();
		
		Type[] types = Type.values();
		
		CSVFile<String> combCSV = null;
		
		HashSet<Integer> uniqueEvents = new HashSet<>();
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Crustal and Subduction Co-Rupture Ground Motions");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Sites and CSV Files");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("Site Name", "Latitude", "Longitude", "Record Count", "CSV File");
		
		List<BBP_Site> highlightSites = new ArrayList<>();
		
		for (BBP_Site site : sites) {
			CSVFile<String> csv = new CSVFile<>(true);
			
			if (highlightSiteNames.contains(site.getName()))
				highlightSites.add(site);
			
			List<String> header = new ArrayList<>();
			header.add("Site Name");
			header.add("Event ID");
			
			for (Type type : types)
				header.add(type.prefix+" Moment");
			for (Type type : types)
				header.add(type.prefix+" Magnitude");
			for (Type type : types)
				header.add(type.prefix+" rRup");
			for (Type type : types)
				header.add(type.prefix+" rJB");
			for (double period : periods)
				for (Type type : types)
					header.add(type.prefix+" "+oDF.format(period)+"s RotD50");
			
			csv.addLine(header);
			
			if (combCSV == null) {
				combCSV = new CSVFile<>(true);
				combCSV.addLine(header);
			}
			
			int matches = 0;
			
			Site gmmSite = sitesBBPtoGMPE.get(site);
			
			for (RSQSimEvent fullEvent : events) {
				int eventID = fullEvent.getID();
				RSQSimEvent crustalEvent = crustalMap.get(eventID);
				RSQSimEvent subductionEvent = subductionMap.get(eventID);
				
				if (crustalEvent != null && crustalLoader.contains(site, eventID)
						&& subductionEvent != null && subductionLoader.contains(site, eventID)) {
					// exists for both
					matches++;
					uniqueEvents.add(eventID);
					
					Map<Type, RSQSimEvent> typeEventMap = new EnumMap<>(Type.class);
					Map<Type, RSQSimEqkRupture> typeEqkRupMap = new EnumMap<>(Type.class);
					Map<Type, DiscretizedFunc> typeRDMap = new EnumMap<>(Type.class);
					for (Type type : types) {
						RSQSimEvent event;
						RSQSimCatalog catalog;
						BBP_CatalogSimZipLoader loader;
						switch (type) {
						case FULL:
							event = fullEvent;
							catalog = fullCatalog;
							loader = fullLoader;
							break;
						case CRUSTAL:
							event = crustalMap.get(eventID);
							catalog = crustalCatalog;
							loader = crustalLoader;
							break;
						case SUBDUCTION:
							event = subductionMap.get(eventID);
							catalog = subductionCatalog;
							loader = subductionLoader;
							break;

						default:
							throw new IllegalStateException();
						}
						typeEventMap.put(type, event);
						typeEqkRupMap.put(type, catalog.getEqkRupture(event));
						typeRDMap.put(type, loader.readRotD50(site, eventID));
						
					}
					
					List<String> line = new ArrayList<>(header.size());
					line.add(site.getName());
					line.add(eventID+"");
					for (Type type : types)
						line.add((float)moment(typeEventMap.get(type))+"");
					for (Type type : types)
						line.add((float)typeEventMap.get(type).getMagnitude()+"");
					for (Type type : types)
						line.add((float)typeEqkRupMap.get(type).getRuptureSurface().getDistanceRup(site.getLoc())+"");
					for (Type type : types)
						line.add((float)typeEqkRupMap.get(type).getRuptureSurface().getDistanceJB(site.getLoc())+"");
					for (double period : periods)
						for (Type type : types)
							line.add(typeRDMap.get(type).getInterpolatedY(period)+"");
					
					csv.addLine(line);
					combCSV.addLine(line);
					
					double[] siteEventFracts = new double[periods.length];
					
					for (int p=0; p<periods.length; p++) {
						double fullGM = typeRDMap.get(Type.FULL).getInterpolatedY(periods[p]);
						double crustalGM = typeRDMap.get(Type.CRUSTAL).getInterpolatedY(periods[p]);
						double subductionGM = typeRDMap.get(Type.SUBDUCTION).getInterpolatedY(periods[p]);
						Preconditions.checkState(fullGM > 0, "Bad fullGM=%s for event %s, p=%s, site %s",
								fullGM, eventID, periods[p], site.getName());
						Preconditions.checkState(crustalGM > 0, "Bad crustalGM=%s for event %s, p=%s, site %s",
								crustalGM, eventID, periods[p], site.getName());
						Preconditions.checkState(subductionGM > 0, "Bad subductionGM=%s for event %s, p=%s, site %s",
								subductionGM, eventID, periods[p], site.getName());
						
						double fract = crustalGM > subductionGM ? subductionGM/crustalGM : crustalGM/subductionGM;
						siteEventFracts[p] = fract;
						
						for (int q=0; q<quantities.length; q++) {
							double combVal = quantities[q].calc(subductionGM, crustalGM);
							DiscretizedFunc func = combinedSimFuncs.get(q).get(eventID, gmmSite);
							if (func == null) {
								func = new LightFixedXFunc(periods, new double[periods.length]);
								combinedSimFuncs.get(q).put(eventID, gmmSite, func);
							}
							func.set(p, combVal);
							for (int i=0; i<minFracts.length; i++)
								if (fract >= minFracts[i])
									scatters[i][q][p].set(fullGM, combVal);
						}
					}
					
					eventSitePeriodFracts.put(eventID, gmmSite, siteEventFracts);
				}
			}
			System.out.println("Found "+matches+" matches for "+site.getName());
			if (matches > 0) {
				File csvFile = new File(resourcesDir, site.getName()+".csv");
				csv.writeToFile(csvFile);
				table.addLine(site.getName(), (float)site.getLoc().lat, (float)site.getLoc().lon, matches+" events",
						"["+csvFile.getName()+"]("+resourcesDir.getName()+"/"+csvFile.getName()+")");
			}
		}
		
		System.out.println("Found "+(combCSV.getNumRows()-1)+" records for "+uniqueEvents.size()+" events");
		File csvFile = new File(resourcesDir, "all_sites.csv");
		table.addLine("__All Sites__", "_N/A_", "_N/A_", (combCSV.getNumRows()-1)+" records for "+uniqueEvents.size()+" events",
				"["+csvFile.getName()+"]("+resourcesDir.getName()+"/"+csvFile.getName()+")");
		combCSV.writeToFile(csvFile);
		
		lines.addAll(table.build());
		lines.add("");

		String crustalShortName = crustalGMM.getShortName();
		String subductionShortName = subductionGMM.getShortName();
		
		for (int i=0; i<minFracts.length; i++) {
			String fractStr, fractPrefix;
			if (minFracts[i] == 0d) {
				fractStr = "All Events";
				fractPrefix = "all";
			} else {
				fractStr = "Min Ratio: "+(float)minFracts[i];
				fractPrefix = "fract"+(float)minFracts[i];
			}
			lines.add("## Scatter Plots, "+fractStr);
			lines.add(topLink); lines.add("");
			
			for (int q=0; q<quantities.length; q++) {
				lines.add("### "+quantities[q].name+" Scatter, "+fractStr);
				lines.add(topLink); lines.add("");

				table = MarkdownUtils.tableBuilder();
				table.addLine("Period", "Linear-Linear", "Log-Log");

				Table<Integer, Site, DiscretizedFunc> combinedRupComps = combinedGMMs.get(q);
				Table<Integer, Site, DiscretizedFunc> combinedSimRupComps = combinedSimFuncs.get(q);
				
				for (int p=-1; p<periods.length; p++) {
					String prefix = "scatter_"+fractPrefix+"_"+quantities[q].name();
					
					DefaultXY_DataSet scatter;
					double period;
					String label;
					if (p < 0) {
						label = "All";
						scatter = union(scatters[i][q]);
						prefix += "_all";
						period = Double.NaN;
					} else {
						label = oDF.format(periods[p])+"s";
						scatter = scatters[i][q][p];
						prefix += "_"+oDF.format(periods[p])+"s";
						period = periods[p];
					}
					
					writeScatter(resourcesDir, prefix, period, scatter, quantities[q].axisLabel);
					
					table.addLine("__"+label+"__", "![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+".png)",
							"![Scatter Plot]("+resourcesDir.getName()+"/"+prefix+"_log.png)");
				}
				
				lines.addAll(table.build());
				lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				table.addLine("Period", "Combined Simulated Residuals", "Combined Empirical Residuals");
				
				for (int p=-1; p<periods.length; p++) {
					String prefix = "residual_"+fractPrefix+"_"+quantities[q].name();
					
					Double limitPeriod;
					String simulatedTitle = quantities[q].name+" vs Simulated Full";
					String gmmTitle = "GMM "+quantities[q].name+" vs Simulated Full";
					String label;
					if (p < 0) {
						label = "All";
						prefix += "_all";
						limitPeriod = null;
						simulatedTitle += ", All Periods";
						gmmTitle += ", All Periods";
					} else {
						label = oDF.format(periods[p])+"s";
						prefix += "_"+oDF.format(periods[p])+"s";
						limitPeriod = periods[p];
						simulatedTitle += ", "+label+" SA";
						gmmTitle += ", "+label+" SA";
					}

					writeResidualPlot(resourcesDir, prefix+"_gmm", Color.GREEN.darker(), gmmTitle,
							"Full Simulated", "Combined Empirical GMM", periods,
							limitPeriod, minFracts[i], fullLoader, eventSitePeriodFracts, combinedRupComps, sitesBBPtoGMPE.inverse());
					writeResidualPlot(resourcesDir, prefix+"_sim", Color.GRAY, simulatedTitle,
							"Full Simluated", "Combined Simulated", periods,
							limitPeriod, minFracts[i], fullLoader, eventSitePeriodFracts, combinedSimRupComps, sitesBBPtoGMPE.inverse());
					
					table.addLine("__"+label+"__", "![Simulated Hist]("+resourcesDir.getName()+"/"+prefix+"_sim.png)",
							"![Empirical Hist]("+resourcesDir.getName()+"/"+prefix+"_gmm.png)");
				}
				
				lines.addAll(table.build());
				lines.add("");
			}
			
			lines.add("### GMM Residuals, "+fractStr);
			lines.add(topLink); lines.add("");
			
			for (boolean useSubduction : new boolean[] {false,true}) {
				String gmPrefix;
				Table<Integer, Site, DiscretizedFunc> fullRupComps;
				Table<Integer, Site, DiscretizedFunc> partialRupComps;
				BBP_CatalogSimZipLoader partialSimLoader;
				String partialName;
				Color color;
				String gmmName;
				if (useSubduction) {
					lines.add("### Subduction GMM Residuals, "+fractStr);
					lines.add("");
					gmPrefix = fractPrefix+"_subduction_gmm";
					partialName = "Subduction Rupture Portion";
					lines.add("Residuals between simulated ground motions and the subduction GMM, "+subductionShortName);
					fullRupComps = fullRupSubductionGMM;
					partialRupComps = partialRupSubductionGMM;
					color = Color.BLUE.darker();
					gmmName = subductionShortName;
					partialSimLoader = subductionLoader;
				} else {
					lines.add("### Crustal GMM Residuals, "+fractStr);
					lines.add("");
					gmPrefix = fractPrefix+"_crustal_gmm";
					partialName = "Crustal Rupture Portion";
					lines.add("Residuals between simulated ground motions and the crustal GMM, "+crustalShortName);
					fullRupComps = fullRupCrustalGMM;
					partialRupComps = partialRupCrustalGMM;
					color = Color.RED.darker();
					gmmName = crustalShortName;
					partialSimLoader = crustalLoader;
				}
				lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				table.addLine("Period", partialName, "Full Rupture");
				
				for (int p=-1; p<periods.length; p++) {
					String prefix = gmPrefix;
					
					Double myPeriod;
					String fullTitle;
					String partialTitle;
					String label;
					if (p < 0) {
						label = "All";
						prefix += "_all";
						myPeriod = null;
						fullTitle = "Full Rupture, All Periods";
						partialTitle = "Partial Rupture, All Periods";
					} else {
						label = oDF.format(periods[p])+"s";
						prefix += "_"+oDF.format(periods[p])+"s";
						myPeriod = periods[p];
						fullTitle = "Full Rupture, "+label;
						partialTitle = partialName+", "+label;
					}
					
					writeResidualPlot(resourcesDir, prefix+"_partial", color, partialTitle,
							"Partial Simulated", gmmName, periods, myPeriod, minFracts[i],
							partialSimLoader, eventSitePeriodFracts, partialRupComps, sitesBBPtoGMPE.inverse());
					writeResidualPlot(resourcesDir, prefix+"_full", color, fullTitle,
							"Full Simulated", gmmName, periods, myPeriod, minFracts[i],
							fullLoader, eventSitePeriodFracts, fullRupComps, sitesBBPtoGMPE.inverse());
					
					table.addLine("__"+label+"__", "![Partial Hist]("+resourcesDir.getName()+"/"+prefix+"_partial.png)",
							"![Full Plot]("+resourcesDir.getName()+"/"+prefix+"_full.png)");
				}
				
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		if (maxMaps > 0) {
			lines.add("## Event Map GMM Residual");
			lines.add(topLink); lines.add("");
			
			lines.add("Debug residual maps of simulations vs each GMM, all for "+oDF.format(mapPeriod)+"s");
			lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumns(crustalShortName+" Residuals", subductionShortName+" Residuals");
			for (Quantity q : mapQuantities)
				table.addColumn("Combined GMM Residuals, "+q.name);
			table.finalizeLine();
			
			List<RSQSimEvent> mapEvents = new ArrayList<>();
			for (RSQSimEvent event : eventsMap.values()) {
				if (fullRupCrustalGMM.containsRow(event.getID()))
					mapEvents.add(event);
			}
			Collections.sort(mapEvents);
			if (mapEvents.size() > maxMaps) {
				System.out.println("Will plot "+maxMaps+"/"+mapEvents.size()+" events");
				mapEvents = mapEvents.subList(0, maxMaps);
			} else {
				System.out.println("Will plot all "+mapEvents.size()+" events");
			}
			for (RSQSimEvent event : mapEvents) {
				System.out.println("Plotting maps for event "+event.getID());
				table.initNewLine();
				
				String prefix = "event_"+event.getID();
				plotResidualMap(resourcesDir, prefix+"_crustal", event, fullLoader, fullRupCrustalGMM.row(event.getID()),
						 mapPeriod, sitesBBPtoGMPE.inverse(), "Log Residuals: "+crustalShortName+" - Full Simulated");
				table.addColumn("![Crustal Map]("+resourcesDir.getName()+"/"+prefix+"_crustal.png)");
				plotResidualMap(resourcesDir, prefix+"_subduction", event, fullLoader, fullRupSubductionGMM.row(event.getID()),
						 mapPeriod, sitesBBPtoGMPE.inverse(), "Log Residuals: "+subductionShortName+" - Full Simulated");
				table.addColumn("![Subduction Map]("+resourcesDir.getName()+"/"+prefix+"_subduction.png)");
				for (Quantity q : mapQuantities) {
					int qIndex = -1;
					for (int i=0; i<quantities.length; i++)
						if (quantities[i] == q)
							qIndex = i;
					Preconditions.checkState(qIndex >= 0);
					String qPrefix = prefix+"_"+q.name();
					plotResidualMap(resourcesDir, qPrefix, event, fullLoader, combinedGMMs.get(qIndex).row(event.getID()),
							 mapPeriod, sitesBBPtoGMPE.inverse(), "Log Residuals: Combined GMM ("+q.name+") - Full Simulated");
					table.addColumn("!["+q.name+" Map]("+resourcesDir.getName()+"/"+qPrefix+".png)");
				}
				
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
		
		exec.shutdown();
	}
	
	private static DefaultXY_DataSet union(DefaultXY_DataSet[] datas) {
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (DefaultXY_DataSet data : datas)
			for (Point2D pt : data)
				ret.set(pt);
		return ret;
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static double moment(RSQSimEvent event) {
		double ret = 0d;
		for (EventRecord rec : event)
			ret += rec.getMoment();
		return ret;
	}
	
	private static void writeScatter(File outputDir, String prefix, double period, DefaultXY_DataSet scatter,
			String yAxisLabel) throws IOException {
		double min = Math.min(scatter.getMinX(), scatter.getMinY());
		double max = Math.max(scatter.getMaxX(), scatter.getMaxY());
		
		String title = Double.isFinite(period) ? oDF.format(period)+"s SA" : "All Periods";
		
		for (boolean log : new boolean[] {false, true}) {
			Range range;
			if (log)
				range = new Range(Math.pow(10, Math.floor(Math.log10(min))),
						Math.pow(10, Math.ceil(Math.log10(max))));
			else
				range = new Range(0d, max*1.05);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
			oneToOne.set(range.getLowerBound(), range.getLowerBound());
			oneToOne.set(range.getUpperBound(), range.getUpperBound());
			
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
			
			funcs.add(scatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Full Ground Motion (g)", yAxisLabel+" (g)");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, log, log, range, range);
			
			PlotUtils.writePlots(outputDir, prefix+(log ? "_log" : ""), gp, 800, false, true, false, false);
		}
	}
	
	private static Table<Integer, Site, DiscretizedFunc> calcGMM(RSQSimCatalog catalog,
			List<Site> sites, BBP_CatalogSimZipLoader simLoader, double[] periods, Supplier<ScalarIMR> gmmRef,
			ExecutorService exec) {
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		
		Table<Integer, Site, Future<DiscretizedFunc>> futures = HashBasedTable.create();
		
		for (Site site : sites) {
			for (RSQSimEvent event : simLoader.getRupturesForSite(site)) {
				GMM_CalcCallable call = new GMM_CalcCallable(gmmDeque, gmmRef, site,
						event.getID(), catalog.getEqkRupture(event), periods);
				futures.put(event.getID(), site, exec.submit(call));
			}
		}
		
		Table<Integer, Site, DiscretizedFunc> ret = HashBasedTable.create();
		try {
			for (Cell<Integer, Site, Future<DiscretizedFunc>> cell : futures.cellSet()) {
				DiscretizedFunc gms = futures.get(cell.getRowKey(), cell.getColumnKey()).get();
				ret.put(cell.getRowKey(), cell.getColumnKey(), gms);
			}
		} catch (InterruptedException | ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		System.out.println("Calculated GMs for "+ret.rowKeySet().size()+" events, "+ret.size()+" event/site pairs");
		
		return ret;
	}
	
	private static class GMM_CalcCallable implements Callable<DiscretizedFunc> {
		
		private ArrayDeque<ScalarIMR> gmmDeque;
		private Supplier<ScalarIMR> gmmRef;
		private Site site;
		private int eventID;
		private EqkRupture rup;
		private double[] periods;

		public GMM_CalcCallable(ArrayDeque<ScalarIMR> gmmDeque, Supplier<ScalarIMR> gmmRef, Site site, int eventID, EqkRupture rup, double[] periods) {
			this.gmmDeque = gmmDeque;
			this.gmmRef = gmmRef;
			this.site = site;
			this.eventID = eventID;
			this.rup = rup;
			this.periods = periods;
		}

		@Override
		public DiscretizedFunc call() throws Exception {
			ScalarIMR gmm = null;
			synchronized (gmmDeque) {
				if (!gmmDeque.isEmpty())
					gmm = gmmDeque.pop();
			}
			if (gmm == null) {
				gmm = gmmRef.get();
				gmm.setParamDefaults();
				gmm.setIntensityMeasure(SA_Param.NAME);
			}
			
			Parameter<?> im = gmm.getIntensityMeasure();
			gmm.setAll(rup, site, im);
			double[] vals = new double[periods.length];
			boolean D = eventID == DEBUG_EVENT_ID && site.getName().equals(DEBUG_SITE_NAME);
			for (int p=0; p<periods.length; p++) {
				SA_Param.setPeriodInSA_Param(im, periods[p]);
				vals[p] = Math.exp(gmm.getMean());
			}
			if (D) {
				System.out.println("GMM Debug for "+eventID+", "+site.getName()+", "+gmm.getShortName());
				System.out.println("\tMag: "+rup.getMag());
				System.out.println("\tDip: "+rup.getRuptureSurface().getAveDip());
				System.out.println("\tRake: "+rup.getAveRake());
				System.out.println("\tArea: "+rup.getRuptureSurface().getArea());
				System.out.println("\trRup: "+gmm.getParameter(DistanceRupParameter.NAME).getValue());
				System.out.print("\tGMs: ");
				for (double val : vals)
					System.out.print((float)Math.log(val)+", ");
				System.out.println();
				System.out.flush();
			}
			
			synchronized (gmmDeque) {
				gmmDeque.push(gmm);
			}
			return new LightFixedXFunc(periods, vals);
		}
		
	}
	
	private static void writeResidualPlot(File outputDir, String prefix, Color color, String title,
			String simName, String modelName,
			double[] allPeriods, Double limitPeriod, double minFract,
			BBP_CatalogSimZipLoader simLoader, Table<Integer, Site, double[]> eventSitePeriodFracts,
			Table<Integer, Site, DiscretizedFunc> gmmResults, Map<Site, BBP_Site> sitesGMPEtoBBP) throws IOException {
		List<Double> residuals = new ArrayList<>();
		for (Cell<Integer, Site, double[]> cell : eventSitePeriodFracts.cellSet()) {
			int eventID = cell.getRowKey();
			Site site = cell.getColumnKey();
			BBP_Site bppSite = sitesGMPEtoBBP.get(site);
			double[] fracts = cell.getValue();
			Preconditions.checkState(fracts.length == allPeriods.length);
			boolean D = limitPeriod == null && eventID == DEBUG_EVENT_ID && site.getName().equals(DEBUG_SITE_NAME);
			if (D) System.out.println("Residual debug for "+eventID+", "+site.getName()+", "+simName+" vs "+modelName+";\tTitle: "+title);
			for (int p=0; p<fracts.length; p++) {
				if (fracts[p] >= minFract && (limitPeriod == null || limitPeriod.doubleValue() == allPeriods[p])) {
					DiscretizedFunc gmmFunc = gmmResults.get(eventID, site);
					Preconditions.checkNotNull(gmmFunc);
					DiscretizedFunc simFunc = simLoader.readRotD50(bppSite, eventID);
					Preconditions.checkNotNull(simFunc);
					double logSim = Math.log(simFunc.getY(allPeriods[p]));
					double logGMM = Math.log(gmmFunc.getY(allPeriods[p]));
					double residual = logGMM - logSim;
					if (D)
						System.out.println("\tT="+(float)allPeriods[p]+": "
								+(float)logGMM+" - "+(float)logSim+" = "+(float)residual);
					residuals.add(residual);
				}
			}
		}
		
		double[] residualArray = Doubles.toArray(residuals);
		double mean = StatUtils.mean(residualArray);
		double stdDev = Math.sqrt(StatUtils.variance(residualArray));
		
//		Range xRange = new Range(-3, 3);
		Range xRange = new Range(-2, 2);
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(
				xRange.getLowerBound()+0.01, xRange.getUpperBound()-0.01, 0.1);
		for (double val : residuals)
			hist.add(hist.getClosestXIndex(val), 1d);
		
		Range yRange = new Range(0d, hist.getMaxY()*1.1);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
		meanLine.set(mean, yRange.getLowerBound());
		meanLine.set(mean, yRange.getUpperBound());
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
		
		funcs.add(meanLine);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Log Residual: "+modelName+" - "+simName, "Count");
		
		DecimalFormat df = new DecimalFormat("0.00");
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		XYTextAnnotation meanAnn = new XYTextAnnotation("Mean: "+df.format(mean),
				xRange.getUpperBound()-0.1, yRange.getUpperBound()*0.975);
		meanAnn.setFont(font);
		meanAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(meanAnn);
		XYTextAnnotation sdAnn = new XYTextAnnotation("Std. Dev: "+df.format(stdDev),
				xRange.getUpperBound()-0.1, yRange.getUpperBound()*0.925);
		sdAnn.setFont(font);
		sdAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(sdAnn);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 700, true, true, false);
	}
	
	private static void plotResidualMap(File outputDir, String prefix, RSQSimEvent event,
			BBP_CatalogSimZipLoader simLoader, Map<Site, DiscretizedFunc> modelGMs,
			double period, Map<Site, BBP_Site> sitesGMPEtoBBP, String label) throws IOException {
		Region region = new Region(new Location(-33.5, 165), new Location(-47, 179));
		RSQSimGeographicMapMaker mapMaker = new RSQSimGeographicMapMaker(region, PoliticalBoundariesData.loadNZOutlines());
		
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setWritePDFs(false);
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		MinMaxAveTracker depTrack = new MinMaxAveTracker();
		List<Double> depths = new ArrayList<>();
		for (SimulatorElement elem : event.getAllElements()) {
			for (Location loc : elem.getVertices()) {
				latTrack.addValue(loc.getLatitude());
				lonTrack.addValue(loc.getLongitude());
				depTrack.addValue(loc.getDepth());
			}
			depths.add(elem.getAveDepth());
		}
//		System.out.println("Element locations:");
//		System.out.println("\tLatitude:\t"+latTrack);
//		System.out.println("\tLongitude:\t"+lonTrack);
//		System.out.println("\tDepth:\t"+depTrack);
		
//		mapMaker.plotEvent(event, Color.LIGHT_GRAY, null, Float.NaN);
		CPT depthCPT = new CPT(0d, depTrack.getMax(), Color.BLACK, Color.LIGHT_GRAY);
		mapMaker.plotEventFillScalars(event, depths, depthCPT, "Depth (km)");
		mapMaker.plotEventHypocenter(Color.GREEN.darker());
		
		CPT residCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-2d, 2d);
		
		List<Location> siteLocs = new ArrayList<>();
		List<Double> siteResiduals = new ArrayList<>();
		System.out.println("Have "+modelGMs.size()+" GMM sites");
		for (Site site : modelGMs.keySet()) {
			BBP_Site bbpSite = sitesGMPEtoBBP.get(site);
			
			if (simLoader.contains(bbpSite, event.getID())) {
				double simVal = Math.log(simLoader.readRotD50(bbpSite, event.getID()).getY(period));
				double gmVal = Math.log(modelGMs.get(site).getY(period));
				double residual = gmVal - simVal;
				siteLocs.add(site.getLocation());
				siteResiduals.add(residual);
			}
		}
		System.out.println("Plotting "+siteResiduals.size()+" site residuals");
		
		mapMaker.plotScatterScalars(siteLocs, siteResiduals, residCPT, label);
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_TRIANGLE, 8f, PlotSymbol.TRIANGLE, Color.BLACK);
		
//		mapMaker.plotScatters(siteLocs, null);
		
		mapMaker.plot(outputDir, prefix, "Event "+event.getID()+", M"+oDF.format(event.getMagnitude())+", "+oDF.format(period)+"s SA");
	}

}
