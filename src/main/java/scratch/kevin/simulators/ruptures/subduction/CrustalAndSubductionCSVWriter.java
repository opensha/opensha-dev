package scratch.kevin.simulators.ruptures.subduction;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimEqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBiMap;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_CatalogSimZipLoader;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

public class CrustalAndSubductionCSVWriter {
	
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
		SUM("Sum", "Subduction + Crustal)") {

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

	public static void main(String[] args) throws IOException {
		RSQSimCatalog fullCatalog = Catalogs.BRUCE_5566.instance();
		RSQSimCatalog crustalCatalog = Catalogs.BRUCE_5566_CRUSTAL.instance();
		RSQSimCatalog subductionCatalog = Catalogs.BRUCE_5566_SUB.instance();
		
		File bbpBaseDir = new File("/data/kevin/bbp/parallel");
		
		File fullBBPdir = new File(bbpBaseDir,
				"2023_03_30-rundir5566-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
		File crustalBBPdir = new File(bbpBaseDir,
//				"2023_06_27-rundir5566_crustal-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
				"2023_06_27-rundir5566_crustal_corupture-all-m5.0-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
		File subductionBBPdir = new File(bbpBaseDir,
//				"2023_06_27-rundir5566_subduction-all-m6.5-skipYears5000-maxDist500-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
				"2023_06_27-rundir5566_subduction_corupture-all-m5.0-skipYears5000-noHF-vmLA_BASIN_500-standardSitesNZ-griddedSitesNZ");
		
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
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		for (RSQSimEvent event : events)
			eventsMap.put(event.getID(), event);
		
		Map<Integer, RSQSimEvent> crustalMap = crustalCatalog.loader().minMag(minMag).loadMap();
		Map<Integer, RSQSimEvent> subductionMap = subductionCatalog.loader().minMag(minMag).loadMap();

		BBP_CatalogSimZipLoader fullLoader = new BBP_CatalogSimZipLoader(fullZip, sites, sitesBBPtoGMPE, eventsMap);
		BBP_CatalogSimZipLoader crustalLoader = new BBP_CatalogSimZipLoader(curstalZip, sites, sitesBBPtoGMPE, eventsMap);
		BBP_CatalogSimZipLoader subductionLoader = new BBP_CatalogSimZipLoader(subductionZip, sites, sitesBBPtoGMPE, eventsMap);
		
		double[] periods = { 2d, 3d, 5d, 10d };
		
		double[] minFracts = { 0d, 0.2, 0.5 };
		
		Quantity[] quantities = Quantity.values();
		
		DefaultXY_DataSet[][][] scatters = new DefaultXY_DataSet[minFracts.length][quantities.length][periods.length];
		
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
		
		for (BBP_Site site : sites) {
			CSVFile<String> csv = new CSVFile<>(true);
			
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
						
						for (int i=0; i<minFracts.length; i++) {
							if (fract >= minFracts[i]) {
								for (int q=0; q<quantities.length; q++)
									scatters[i][q][p].set(fullGM, quantities[q].calc(subductionGM, crustalGM));
							}
						}
					}
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
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
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

}
