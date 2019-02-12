package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

public class BBP_PartBValidationPageGen {
	
	private RSQSimCatalog catalog;
	private BBP_PartBSimZipLoader loader;
	private int skipYears;
	private int numSites;
	private boolean randomAz;
	private double vs30;

	public BBP_PartBValidationPageGen(RSQSimCatalog catalog, BBP_PartBSimZipLoader loader, int skipYears, int numSites,
			boolean randomAz, double vs30) {
		this.catalog = catalog;
		this.loader = loader;
		this.skipYears = skipYears;
		this.numSites = numSites;
		this.randomAz = randomAz;
		this.vs30 = vs30;
	}
	
	public void generatePage(File outputDir) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		List<String> lines = new ArrayList<>();
		
		lines.add("# "+catalog.getName()+" BBP Part B Validation");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		lines.add("## Background Information");
		lines.add(topLink); lines.add("");
		lines.add("");
		lines.add("This page reproduces the SCEC BroadBand Platform \"Part B\" validation exercise as defined in:");
		lines.add("");
		lines.add("*Goulet, C. A., Abrahamson, N. A., Somerville, P. G., & Wooddell, K. E. (2014). The SCEC broadband platform "
				+ "validation exercise: Methodology for code validation in the context of seismic‐hazard analyses. "
				+ "Seismological Research Letters, 86(1), 17-26.* [(link)](https://pubs.geoscienceworld.org/ssa/srl/article/86/1/17/315438/"
				+ "the-scec-broadband-platform-validation-exercise)");
		lines.add("");
		lines.add("The goal of this exercise was to validate BBP simulation methods (both rupture generation and ground motion simulation) "
				+ "against the NGA-West GMPEs (the original study used NGA-West1, we use NGA-West2) for scenario ruptures where the NGA "
				+ "relations are well constrained:");
		lines.add("");
		lines.add("* M 5.5, 45°-dipping reverse, Ztor = 6 km");
		lines.add("* M 6.2, vertical strike slip, Ztor = 4 km");
		lines.add("* M 6.6, vertical strike slip with a surface rupture");
		lines.add("* M 6.6, 45°-dipping reverse, Ztor = 3 km");
		lines.add("");
		lines.add("50 rupture realizations were generated for each scenario with randomly distributed hypocenters, and 40 sites were distributed "
				+ "at random azimuths on the footwall side of the faults at rupture distances of 20 and 50 km. Resultant ground motions were "
				+ "compared against an evaluation criterion which \"was established so as to be wide enough to limit a pass/fail grade for each "
				+ "scenario considered.\" This criterion only applies at periods up to "
				+ optionalDigitDF.format(BBP_PartBValidationConfig.BBP_MAX_ACCEPTANCE_PERIOD)+"s, becase data above this period \"are fairly "
						+ "sparse and cannot provide a reliable constraint.\"");
		lines.add("");
		lines.add("A method is said to pass the test if the median RotD50 value is within the evaluation criteria at every spectral period. "
				+ "\"Departure from that range is a definite sign that the model is not consistent with our current dataset and is a sign of "
				+ "potential issues with the simulations.\"");
		lines.add("");
		lines.add("## RSQSim BBP Part B Methodology");
		lines.add(topLink); lines.add("");
		lines.add("We reproduce the Part B experiment using RSQSim as the rupture generator, coupled with the Graves & Pitarka ground motion "
				+ "simulation method. While the original BBP Part B validation exercise was prescriptive (the magnitude and fault surface for "
				+ "each scenario was an input to the rupture generators), we can't prescribe RSQSim ruptures. Instead we search catalogs for "
				+ "ruptures which are very similar to the BBP Part B scenarios, and distribute sites around those ruptures. This algorithm "
				+ "can be a little tricky for non-rectangular dipping ruptures. The specific matching criteria for each scenario along with "
				+ "example plots of ruptures and sites can be found under each scenario.");
		lines.add("");
		lines.add("*NOTE: We only currently only consider the larger M6.6 scenarios, and spectral periods "
				+optionalDigitDF.format(BBP_PartBValidationConfig.BBP_MIN_ACCEPTANCE_PERIOD)+"s or larger*");
		
		Color[] distColors = {Color.BLUE.darker(), Color.GREEN.darker()};
		String[] colorNames = {"Blue", "Geen"};
		double[] distances = BBP_PartBValidationConfig.OFFICIAL_DISTANCES;
		Preconditions.checkState(distColors.length == distances.length);
		
		lines.add("");
		lines.add("## Results Summary Table");
		lines.add("");
		int resultsTableIndex = lines.size();
		TableBuilder resultsTable = MarkdownUtils.tableBuilder();
		resultsTable.initNewLine();
		resultsTable.addColumn("Scenario");
		for (double distance : distances)
			resultsTable.addColumn(optionalDigitDF.format(distance)+" km");
		resultsTable.finalizeLine();
		lines.add("");
		
		for (Scenario scenario : Scenario.values()) {
			if (!loader.hasScenario(scenario))
				continue;
			
			System.out.println("Doing scenario: "+scenario.getName());
			
			lines.add("## "+scenario.getName());
			lines.add(topLink); lines.add("");
			lines.add("### "+scenario.getShortName()+" RSQSim Rupture Match Criteria");
			lines.add(topLink); lines.add("");
			List<RSQSimEvent> matches = scenario.getMatches(catalog, skipYears);
			String[] criteria = scenario.getMatchCriteria();
			lines.add(matches.size()+" events in the catalog match the following criteria:");
			lines.add("");
			for (String criterion : criteria)
				lines.add("* "+criterion);
			lines.add("");
			String distsStr = null;
			for (int i=0; i<distColors.length; i++) {
				if (i == 0)
					distsStr = "";
				else
					distsStr += ", ";
				distsStr += optionalDigitDF.format(distances[i])+" km sites in "+colorNames[i];
			}
			lines.add("Example matches ("+distsStr+"):");
			lines.add("");
			TableBuilder builder = MarkdownUtils.tableBuilder();
			builder.initNewLine();
			for (int i=0; i<5 && i<matches.size(); i++) {
				RSQSimEvent event = matches.get(i);
				
				String prefix = scenario.getPrefix()+"_match_"+i+"_event_"+event.getID();
				BBP_PartBValidationConfig.plotEventAndSites(catalog, event, distColors, numSites, randomAz, resourcesDir, prefix);
				
				builder.addColumn("![Event "+event.getID()+"](resources/"+prefix+".png)");
				new File(resourcesDir, prefix+".pdf").delete();
			}
			builder.finalizeLine();
			lines.addAll(builder.build());
			
			resultsTable.initNewLine();
			resultsTable.addColumn(scenario.getShortName());
			
			for (double distance : distances) {
				String labelName = scenario.getShortName()+" "+optionalDigitDF.format(distance)+" km Results";
				String anchorName = MarkdownUtils.getAnchorName(labelName);
				lines.add("### "+labelName);
				lines.add(topLink); lines.add("");
				
				String prefix = scenario.getPrefix()+"_"+optionalDigitDF.format(distance)+"km";
				List<ValidationResult> results = calcPlotScenarioResults(scenario, distance, matches, resourcesDir, prefix);
				boolean pass = true;
				for (ValidationResult result : results)
					pass = pass && result.passes();
				
				if (pass) {
					lines.add("Result: **PASS**");
					resultsTable.addColumn("**[PASS](#"+anchorName+")**");
					System.out.println(scenario+", "+(float)distance+" km: PASS");
				} else {
					lines.add("Result: **FAIL**");
					resultsTable.addColumn("**[FAIL](#"+anchorName+")**");
					System.out.println(scenario+", "+(float)distance+" km: FAIL");
				}
				lines.add("");
				lines.add("![Acceptance Plot](resources/"+prefix+".png)");
				lines.add("");
				builder = MarkdownUtils.tableBuilder();
				builder.initNewLine();
				builder.addColumn("**Period**");
				for (ValidationResult result : results)
					builder.addColumn(optionalDigitDF.format(result.period)+"s");
				builder.finalizeLine();
				builder.initNewLine();
				builder.addColumn("**Lower Bound**");
				for (ValidationResult result : results)
					builder.addColumn(valDF.format(result.lowerBound));
				builder.finalizeLine();
				builder.initNewLine();
				builder.addColumn("**Sim Median**");
				for (ValidationResult result : results)
					if (result.passes())
						builder.addColumn("**"+valDF.format(result.simMedian)+"**");
					else
						builder.addColumn("*"+valDF.format(result.simMedian)+"*");
				builder.finalizeLine();
				builder.initNewLine();
				builder.addColumn("**NGA-W2 Median**");
				for (ValidationResult result : results)
					builder.addColumn(valDF.format(result.dataMedian));
				builder.finalizeLine();
				builder.initNewLine();
				builder.addColumn("**Upper Bound**");
				for (ValidationResult result : results)
					builder.addColumn(valDF.format(result.upperBound));
				builder.finalizeLine();
				lines.addAll(builder.build());
				lines.add("");
				Map<String, String> comparisonFigs = scenario.getPublishedComparisonURLs(distance);
				if (comparisonFigs != null) {
					lines.add("#### "+scenario.getShortName()+" "+optionalDigitDF.format(distance)+" km Comparisons");
					lines.add(topLink); lines.add("");
					builder = MarkdownUtils.tableBuilder();
					List<String> keys = new ArrayList<>(comparisonFigs.keySet());
					Collections.sort(keys);
					builder.initNewLine();
					for (String key : keys)
						builder.addColumn(key);
					builder.finalizeLine();
					builder.initNewLine();
					for (String key : keys)
						builder.addColumn("!["+key+"]("+comparisonFigs.get(key)+")");
					builder.finalizeLine();
					lines.addAll(builder.build());
					lines.add("");
				}
			}
			resultsTable.finalizeLine();
		}
		
		lines.addAll(resultsTableIndex, resultsTable.build());
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private class ValidationResult {
		private double period, lowerBound, upperBound, dataMedian, simMedian, simSD;

		public ValidationResult(double period, double lowerBound, double upperBound, double dataMedian,
				double simMedian, double simSD) {
			super();
			this.period = period;
			this.lowerBound = lowerBound;
			this.upperBound = upperBound;
			this.dataMedian = dataMedian;
			this.simMedian = simMedian;
			this.simSD = simSD;
		}
		
		public String toString() {
			String str = "Period "+optionalDigitDF.format(period)+"s: ";
			if (passes())
				str += "PASS\n";
			else
				str += "FAIL\n";
			str += "\tmedian: "+valDF.format(simMedian)+"\n";
			str += "\tcriterion: ["+valDF.format(lowerBound)+", "+valDF.format(upperBound)+"]";
			return str;
		}
		
		public boolean passes() {
			return simMedian >= lowerBound && simMedian <= upperBound;
		}
	}
	
	private static DecimalFormat valDF = new DecimalFormat("0.0000");
	
	private List<ValidationResult> calcPlotScenarioResults(Scenario scenario, double distance,
			List<RSQSimEvent> events, File resourcesDir, String prefix) throws IOException {
		SummaryStatistics[] lnPeriodStats = null;
		List<List<Double>> lnVals = new ArrayList<>();
		double[] periods = null;
		for (RSQSimEvent event : events) {
			DiscretizedFunc[] rd50s = loader.getRotD50(event.getID(), scenario, distance);
			for (DiscretizedFunc rd50 : rd50s) {
				if (lnPeriodStats == null) {
					lnPeriodStats = new SummaryStatistics[rd50.size()];
					lnVals = new ArrayList<>();
					periods = new double[rd50.size()];
					for (int p=0; p<periods.length; p++) {
						lnPeriodStats[p] = new SummaryStatistics();
						lnVals.add(new ArrayList<>());
						periods[p] = rd50.getX(p);
					}
				} else {
					Preconditions.checkState(lnPeriodStats.length == rd50.size());
				}
				for (int p=0; p<rd50.size(); p++) {
					double lnVal = Math.log(rd50.getY(p));
					lnPeriodStats[p].addValue(lnVal);
					lnVals.get(p).add(lnVal);
				}
			}
		}
		
		DiscretizedFunc gmpeMean = scenario.getMeanPrediction(vs30, distance);
		UncertainArbDiscDataset gmpeBounds = scenario.getAcceptanceCriteria(vs30, distance);
		DiscretizedFunc gmpeUpper = gmpeBounds.getUpper();
		DiscretizedFunc gmpeLower = gmpeBounds.getLower();
		
		double whiskerWidthLog10 = 0.01;
		double boxWidthLog10 = 0.02;
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(gmpeMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		gmpeUpper.setName(gmpeBounds.getName());
		funcs.add(gmpeBounds.getUpper());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		
		gmpeLower.setName(null);
		funcs.add(gmpeLower);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		
		PlotCurveCharacterstics whiskerChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK);
		PlotCurveCharacterstics boxChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE);
		PlotCurveCharacterstics medianChar = new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 4f, Color.RED);
		
		DiscretizedFunc medianFunc = new ArbitrarilyDiscretizedFunc();
		medianFunc.setName("Simulation Median");
		funcs.add(medianFunc);
		chars.add(medianChar);

		boolean firstMedian = true;
		boolean firstWhisker = true;
		
		List<ValidationResult> results = new ArrayList<>();
		
		for (int p=0; p<periods.length; p++) {
			double period = periods[p];
			if (period < BBP_PartBValidationConfig.BBP_MIN_ACCEPTANCE_PERIOD)
				continue;
			
			double logSD = lnPeriodStats[p].getStandardDeviation();
			double logMedian = DataUtils.median(Doubles.toArray(lnVals.get(p)));
			
			double median = Math.exp(logMedian);
			if (period <= BBP_PartBValidationConfig.BBP_MAX_ACCEPTANCE_PERIOD) {
				double acceptMin = gmpeLower.getInterpolatedY_inLogXLogYDomain(period);
				double acceptMax = gmpeUpper.getInterpolatedY_inLogXLogYDomain(period);
				double dataMedian = gmpeMean.getInterpolatedY_inLogXLogYDomain(period);
				ValidationResult result = new ValidationResult(period, acceptMin, acceptMax, dataMedian, median, logSD);
				System.out.println(result);
				results.add(result);
			}
			
			medianFunc.set(period, median);

			double plusSD = Math.exp(logMedian+logSD);
			double minusSD = Math.exp(logMedian-logSD);
			double max = Math.exp(lnPeriodStats[p].getMax());
			double min = Math.exp(lnPeriodStats[p].getMin());
			
			// +/- SD box
			double boxLeft = getXlogDiff(period, -0.5*boxWidthLog10);
			double boxRight = getXlogDiff(period, 0.5*boxWidthLog10);
			// bottom line
			funcs.add(line(boxLeft, minusSD, boxRight, minusSD));
			if (firstMedian) {
				funcs.get(funcs.size()-1).setName("± σ");
				firstMedian = false;
			}
			chars.add(boxChar);
			
			// left line
			funcs.add(line(boxLeft, minusSD, boxLeft, plusSD));
			chars.add(boxChar);
			
			// top line
			funcs.add(line(boxLeft, plusSD, boxRight, plusSD));
			chars.add(boxChar);
			
			// right line
			funcs.add(line(boxRight, minusSD, boxRight, plusSD));
			chars.add(boxChar);
			
			// whisker lines
			funcs.add(line(period, min, period, max));
			if (firstWhisker) {
				funcs.get(funcs.size()-1).setName("Extrema");
				firstWhisker = false;
			}
			chars.add(whiskerChar);
			
			double whiskerLeft = getXlogDiff(period, -0.5*whiskerWidthLog10);
			double whiskerRight = getXlogDiff(period, 0.5*whiskerWidthLog10);
			funcs.add(line(whiskerLeft, min, whiskerRight, min));
			chars.add(whiskerChar);
			funcs.add(line(whiskerLeft, max, whiskerRight, max));
			chars.add(whiskerChar);
		}
		
		// add it again without a name so it's on top
		DiscretizedFunc medianFunc2 = medianFunc.deepClone();;
		medianFunc2.setName(null);
		funcs.add(medianFunc2);
		chars.add(medianChar);
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		for (XY_DataSet func : funcs) {
			minY = Math.min(minY, func.getMinY());
			maxY = Math.max(maxY, func.getMaxY());
		}
		
		Range xRange = new Range(Math.pow(10, Math.log10(BBP_PartBValidationConfig.BBP_MIN_ACCEPTANCE_PERIOD)-0.05),
				Math.pow(10, Math.log10(periods[periods.length-1])+0.05));
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minY))), Math.pow(10, Math.ceil(Math.log10(maxY))));
		
		PlotSpec spec = new PlotSpec(funcs, chars, scenario.getName()+", "+optionalDigitDF.format(distance)+" km", "Period (s)", "RotD50 Sa (g)");
		spec.setLegendVisible(true);
		List<XYTextAnnotation> anns = new ArrayList<>();
		// add counts
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 26);

		double logRangeY = Math.log10(yRange.getUpperBound()) - Math.log10(yRange.getLowerBound());
		XYTextAnnotation eventCountAnn = new XYTextAnnotation(events.size()+" Events  ",
				xRange.getUpperBound(), Math.pow(10, Math.log10(yRange.getLowerBound()) + 0.95*logRangeY));
		eventCountAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		eventCountAnn.setFont(font);
		anns.add(eventCountAnn);

		XYTextAnnotation simCountAnn = new XYTextAnnotation(lnVals.get(0).size()+" Records  ",
				xRange.getUpperBound(), Math.pow(10, Math.log10(yRange.getLowerBound()) + 0.88*logRangeY));
		simCountAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		simCountAnn.setFont(font);
		anns.add(simCountAnn);
		
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setLegendFontSize(18);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
//		gp.getYAxis().setStandardTickUnits(getYTickUnits());
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(resourcesDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(resourcesDir, prefix+".pdf").getAbsolutePath());
		
		return results;
	}
	
	private static double getXlogDiff(double center, double deltaLog10) {
		return Math.pow(10, Math.log10(center)+deltaLog10);
	}
	
	private static XY_DataSet line(double x1, double y1, double x2, double y2) {
		XY_DataSet ret = new DefaultXY_DataSet();
		ret.set(x1, y1);
		ret.set(x2, y2);
		return ret;
	}

	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_2310.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance(baseDir);
		
		System.out.println("Catalog: "+catalog.getName());
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File bbpZipFile = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		int numSites = -1;
		int skipYears = 0;
		boolean randomAz = false;
		VelocityModel vm = null;
		Arrays.sort(allBBPDirs, new FileNameComparator());
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-partB")
					&& name.contains("sites") && name.contains("-vm")) {
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					bbpDir = dir;
					bbpZipFile = zipFile;
					String sitesStr = name.substring(0, name.indexOf("sites"));
					sitesStr = sitesStr.substring(sitesStr.lastIndexOf("-")+1);
					numSites = Integer.parseInt(sitesStr);
					if (name.contains("-skipYears")) {
						String yearsStr = name.substring(name.indexOf("-skipYears")+"-skipYears".length());
						if (yearsStr.contains("-"))
							yearsStr = yearsStr.substring(0, yearsStr.indexOf("-"));
						skipYears = Integer.parseInt(yearsStr);
					}
					randomAz = name.contains("randomAz");
					String vmStr = name.substring(name.indexOf("-vm")+3);
					if (vmStr.contains("-"))
						vmStr = vmStr.substring(0, vmStr.indexOf("-"));
					vm = VelocityModel.valueOf(vmStr);
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("\tInput file: "+bbpZipFile.getName());
		System.out.println("\tNum sites: "+numSites);
		
		BBP_PartBSimZipLoader bbpLoader = new BBP_PartBSimZipLoader(bbpZipFile, numSites);
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File partBDir = new File(catalogOutputDir, "bbp_part_b");
		Preconditions.checkState(partBDir.exists() || partBDir.mkdir());
		
		BBP_PartBValidationPageGen pageGen = new BBP_PartBValidationPageGen(
				catalog, bbpLoader, skipYears, numSites, randomAz, vm.getVs30());
		
		pageGen.generatePage(partBDir);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
