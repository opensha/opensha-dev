package scratch.kevin.simulators.multiFault;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.laughTest.AbstractLaughTest;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.inversion.laughTest.MinSectsPerParentFilter;
import scratch.UCERF3.inversion.laughTest.MinSectsPerParentFilter.CleanupFilter;
import scratch.UCERF3.inversion.laughTest.MinSectsPerParentFilter.ContinualFilter;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;

public class RSQSimU3RupturePageGen {

	public static void main(String[] args) throws IOException, DocumentException, GMT_MapException, RuntimeException {
		File catalogsBaseDir = new File("/data/kevin/simulators/catalogs");
		File mainOutputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");

//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(catalogsBaseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(catalogsBaseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_3271.instance(catalogsBaseDir);
		RSQSimCatalog catalog = Catalogs.JG_tuneBase1m.instance(catalogsBaseDir);
		
		File catalogDir = catalog.getCatalogDir();
		
		File catalogOutputDir = new File(mainOutputDir, catalogDir.getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File multiFaultDir = new File(catalogOutputDir, "multi_fault");
		Preconditions.checkState(multiFaultDir.exists() || multiFaultDir.mkdir());
		
		File resourcesDir = new File(multiFaultDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		File fssDir = new File(catalogDir, "fss");
		Preconditions.checkState(fssDir.exists() || fssDir.mkdir());
		
		String catalogName = catalog.getName();
		String catalogType = "RSQSim";
		String catalogTypeFileName = "rsqsim";
		
		double minMag = 6.5;
		int skipYears = 5000;
		boolean includeNumSects = false;
		
		List<String> lines = new ArrayList<>();
		
		FaultBasedMapGen.MAP_LABEL_SIZE = 24;
		FaultBasedMapGen.MAP_LABEL_TICK_SIZE = 20;
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		FaultBasedMapGen.FAULT_THICKNESS = 4d;
		
		double minFractForInclusion = catalog.getMinSubSectFractForInclusion();
		
		// header
		lines.add("# Multi Fault Rupture Comparisons");
		lines.add("");
		lines.add("*Subsections participate in a rupture if at least "
				+(float)(minFractForInclusion*100d)+" % of its area ruptures*");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		FaultModels fm = catalog.getFaultModel();
		DeformationModels dm = catalog.getDeformationModel();
		
		String catParams = "m"+(float)minMag+"_skip"+skipYears+"_sectArea"+(float)minFractForInclusion;
		File solFile = new File(fssDir, "rsqsim_sol_"+catParams+".zip");
		FaultSystemSolution sol;
		if (solFile.exists()) {
			System.out.println("Loading solution from: "+solFile.getAbsolutePath());
			sol = FaultSystemIO.loadSol(solFile);
		} else {
			System.out.println("Loading events from: "+catalogDir.getAbsolutePath());
			Loader loader = catalog.loader().minMag(minMag).skipYears(skipYears);
			sol = catalog.buildSolution(loader, minMag);
			System.out.println("Writing solution to: "+solFile.getAbsolutePath());
			FaultSystemIO.writeSol(sol, solFile);
		}
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		System.out.println(rupSet.getNumRuptures()+" unique ruptures");
		
		LaughTestFilter filter = LaughTestFilter.getDefault();
		if (!includeNumSects)
			filter.setMinNumSectInRup(0);
		
		File scratchDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/");
		DeformationModelFetcher fetch = new DeformationModelFetcher(fm, dm, scratchDir, 0.1);
		
		List<FaultSectionPrefData> datas = fetch.getSubSectionList();
		
		Map<IDPairing, Double> distances = fetch.getSubSectionDistanceMap(1000d);
		Map<IDPairing, Double> azimuths = fetch.getSubSectionAzimuthMap(distances.keySet());
		Map<Integer, Double> rakesMap = new HashMap<Integer, Double>();
		for (FaultSectionPrefData data : rupSet.getFaultSectionDataList())
			rakesMap.put(data.getSectionId(), data.getAveRake());
		boolean applyGarlockPintoMtnFix = true;
		
		CoulombRates coulombRates = null;
		if (filter.getCoulombFilter() != null) {
			try {
				coulombRates = CoulombRates.loadUCERF3CoulombRates(fm);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
		// we want actual catalog rupture counts before binning into U3 style ruptures
		// find the smallest rate, which will be 1/catLen, then numRups = solRate/minRate
		double minRate = StatUtils.min(sol.getRateForAllRups());
		
		List<List<Integer>> sectionConnectionsListList = SectionClusterList.computeCloseSubSectionsListList(
				datas, distances, filter.getMaxJumpDist(), coulombRates);
		
		List<AbstractLaughTest> tests = filter.buildLaughTests(azimuths, distances, rakesMap, coulombRates, applyGarlockPintoMtnFix,
				sectionConnectionsListList, rupSet.getFaultSectionDataList());
		
		// doesn't come with jump dist filter by default (it is included explicitly in generation for UCERF3)
		tests.add(0, new JumpDistFilter(distances, 5d));
		
		if (includeNumSects) {
			// replace the separate min sects per parent filters with a single one
			tests.add(1, new CombinedMinSectsFilter(removeByClass(tests, ContinualFilter.class),
				removeByClass(tests, CleanupFilter.class)));
		}
		
		Color[] colors = { Color.DARK_GRAY, Color.RED, Color.BLUE, Color.GREEN.darker(), Color.CYAN, Color.ORANGE };
		
		int allPassCount = 0;
		int[] failCounts = new int[tests.size()];
		int[] onlyFailCounts = new int[tests.size()];
		int[] erredCounts = new int[tests.size()];
		
		int tot = 0;
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			int numCatalogOccurances = (int)Math.round(sol.getRateForRup(r)/minRate);
			tot += numCatalogOccurances;
			Preconditions.checkState(numCatalogOccurances >= 1);
			List<FaultSectionPrefData> rupture = rupSet.getFaultSectionDataForRupture(r);
			boolean allPass = true;
			int onlyFailureIndex = -1;
			for (int t=0; t<tests.size(); t++) {
				AbstractLaughTest test = tests.get(t);
				boolean subPass;
				try {
					subPass = test.doesRupturePass(rupture);
				} catch (Exception e) {
					if (erredCounts[t] == 0) {
						System.err.println("First exception for "+test.getName()+":");
						e.printStackTrace();
					}
					erredCounts[t] += numCatalogOccurances;
					subPass = true; // do not fail on error
				}
				if (!subPass && allPass) {
					// this is the first failure
					onlyFailureIndex = t;
				} else if (!subPass) {
					// failed more than 1
					onlyFailureIndex = -1;
				}
				allPass = subPass && allPass;
				if (!subPass)
					failCounts[t] += numCatalogOccurances;
			}
			if (allPass)
				allPassCount += numCatalogOccurances;
			if (onlyFailureIndex >= 0)
				onlyFailCounts[onlyFailureIndex] += numCatalogOccurances;
		}
		System.out.println("Passed all filters: "+countStats(allPassCount, tot));
		for (int t=0; t<tests.size(); t++) {
			System.out.println(tests.get(t).getName());
			System.out.println("\tFailed: "+countStats(failCounts[t], tot));
			System.out.println("\tOnly Failure: "+countStats(onlyFailCounts[t], tot));
			System.out.println("\tErred: "+countStats(erredCounts[t], tot));
		}
		
		// now plot
		double dx = 1d;
		double buffer = 0.2*dx;
		double deltaEachSide = (dx - buffer)/2d;
		float thickness = 80f;
		double maxY = 50;
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		
		List<PlotElement> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(new DefaultXY_DataSet(new double[] {0d, 1d}, new double[] {0d, 0d}));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, Color.WHITE));
		
		List<XYAnnotation> anns = new ArrayList<>();
		
		for (int i=0; i<tests.size(); i++) {
			double x = i*dx + 0.5*dx;
			double percentFailed = 100d*failCounts[i]/tot;
			double percentOnly = 100d*onlyFailCounts[i]/tot;
			double percentErred = 100d*erredCounts[i]/tot;
			
			Color c = colors[i % colors.length];
			
//			funcs.add(vertLine(x, 0, percentFailed));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, c));
			anns.add(filledBox(x-deltaEachSide, 0, x+deltaEachSide, percentFailed, c));
			
			if (percentOnly > 0) {
//				funcs.add(vertLine(x, 0, percentOnly));
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, darker(c)));
				anns.add(filledBox(x-deltaEachSide, 0, x+deltaEachSide, percentOnly, darker(c)));
			}
			
			String title = tests.get(i).getShortName();
			
			if (percentErred > 0) {
//				funcs.add(vertLine(x, percentFailed, percentFailed + percentErred));
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, Color.LIGHT_GRAY));
				anns.add(emptyBox(x-deltaEachSide, percentFailed, x+deltaEachSide, percentFailed + percentErred,
						PlotLineType.DASHED, Color.LIGHT_GRAY, 2f));
				title += "*";
			}
			
			XYTextAnnotation ann = new XYTextAnnotation(title, x, maxY*0.95);
			ann.setTextAnchor(TextAnchor.TOP_CENTER);
			ann.setPaint(c);
			ann.setFont(font);
			
			anns.add(ann);
			
			ann = new XYTextAnnotation(percentDF.format(percentFailed/100d), x, percentFailed+0.6);
			ann.setTextAnchor(TextAnchor.BOTTOM_CENTER);
			ann.setPaint(Color.BLACK);
			ann.setFont(font);
			
			anns.add(ann);
		}
		
		XYTextAnnotation ann = new XYTextAnnotation(
				percentDF.format((double)allPassCount/tot)+" passed all", dx*0.25, maxY*0.88);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setPaint(Color.BLACK);
		ann.setFont(font);
		
		anns.add(ann);
		
		String title = "Rupture Plausibility Filters, M≥"+(float)minMag+", SectArea≥"+(float)minFractForInclusion;
		PlotSpec spec = new PlotSpec(funcs, chars, title, " ", "Percent Failed");
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = new File(resourcesDir, "filters_"+catParams).getAbsolutePath();
		
		gp.drawGraphPanel(spec, false, false, new Range(0, tests.size()*dx), new Range(0, maxY));
		gp.getXAxis().setTickLabelsVisible(false);
//		gp.getXAxis().setvisi
		gp.getChartPanel().setSize(1000, 500);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		
		lines.add("## Plausibility Filter Comparisons");
//		lines.add(topLink); lines.add("");
		lines.add("");
		
		lines.add("### Rupture Failure Percentages");
		lines.add(topLink); lines.add("");
		lines.add("");
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/filters_"+catParams+".png)");
		
		// now jumps
		lines.add("## 1km Jump Count");
		lines.add(topLink); lines.add("");
		lines.add("");
		System.out.println("Plotting num jumps");
		RSQSimRupJumpCompare.plotFixedJumpDist(catalog.getU3CompareSol(), distances, sol, catalogName, minMag, 1d, resourcesDir);
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/jumps_1.0km.png)");
		lines.add("");
		
		// cumulant mag
		System.out.println("Plotting cumulant mag");
		calcCumulantMedianMag(catalog.getU3CompareSol(), sol, catalogType, resourcesDir);
		
		lines.add("## Cumulant Magnitude");
		lines.add(topLink); lines.add("");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3", "Difference");
		File rsPlot = new File(resourcesDir, "mag_cumulant_medians_"+catalogType.toLowerCase()+".png");
		Preconditions.checkState(rsPlot.exists());
		File u3Plot = new File(resourcesDir, "mag_cumulant_medians_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		File diffPlot = new File(resourcesDir, "mag_cumulant_medians_diff.png");
		Preconditions.checkState(diffPlot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")",
				"![Difference]("+resourcesDir.getName()+"/"+diffPlot.getName()+")");
		lines.addAll(table.build());
		lines.add("");
		
		System.out.println("Plotting connectivity clusters");
		plotConnectivity(catalog.getU3CompareSol().getRupSet(), resourcesDir, "connectivity_ucerf3", "UCERF3 Connectivity");
		plotConnectivity(sol.getRupSet(), resourcesDir, "connectivity_"+catalogTypeFileName, catalogName+" Connectivity");
		
		lines.add("## Fault Connectivity");
		lines.add(topLink); lines.add("");
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3");
		rsPlot = new File(resourcesDir, "connectivity_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "connectivity_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, multiFaultDir);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
	}
	
	private static Color darker(Color c) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
//		r += (255-r)/2;
//		g += (255-g)/2;
//		b += (255-b)/2;
		r /= 2;
		g /= 2;
		b /= 2;
		return new Color(r, g, b);
	}
	
	private static DefaultXY_DataSet vertLine(double x, double y0, double y1) {
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.set(x, y0);
		line.set(x, y1);
		return line;
	}
	
	private static XYBoxAnnotation filledBox(double x0, double y0, double x1, double y1, Color c) {
		XYBoxAnnotation ann = new XYBoxAnnotation(x0, y0, x1, y1, null, null, c);
		return ann;
	}
	
	private static XYBoxAnnotation emptyBox(double x0, double y0, double x1, double y1,
			PlotLineType lineType, Color c, float thickness) {
		Stroke stroke = lineType.buildStroke(thickness);
		XYBoxAnnotation ann = new XYBoxAnnotation(x0, y0, x1, y1, stroke, c, null);
		return ann;
	}
	
	private static final DecimalFormat percentDF = new DecimalFormat("0.00%");
	private static String countStats(int count, int tot) {
		return count+"/"+tot+" ("+percentDF.format((double)count/(double)tot)+")";
	}
	
	private static class JumpDistFilter extends AbstractLaughTest {
		
		private Map<IDPairing, Double> distances;
		private double maxJumpDist;
		
		public JumpDistFilter(Map<IDPairing, Double> distances, double maxJumpDist) {
			this.distances = distances;
			this.maxJumpDist = maxJumpDist;
		}

		@Override
		public String getShortName() {
			return "JumpDist";
		}

		@Override
		public String getName() {
			return "Maximum Jump Dist";
		}

		@Override
		public boolean doesLastSectionPass(List<FaultSectionPrefData> rupture, List<IDPairing> pairings,
				List<Integer> junctionIndexes) {
			if (junctionIndexes.isEmpty())
				return true;
			IDPairing pair = pairings.get(junctionIndexes.size()-1);
			return distances.get(pair) <= maxJumpDist;
		}

		@Override
		public boolean isContinueOnFaulure() {
			return false;
		}

		@Override
		public boolean isApplyJunctionsOnly() {
			return true;
		}
	}
	
	private static class CombinedMinSectsFilter extends AbstractLaughTest {
		
		private ContinualFilter continualFilter;
		private CleanupFilter cleanupFilter;

		public CombinedMinSectsFilter(MinSectsPerParentFilter.ContinualFilter continualFilter,
				MinSectsPerParentFilter.CleanupFilter cleanupFilter) {
			this.continualFilter = continualFilter;
			this.cleanupFilter = cleanupFilter;
		}

		@Override
		public String getShortName() {
			return "SectsPerParent";
		}

		@Override
		public String getName() {
			return "Min Sects Per Parent";
		}

		@Override
		public boolean doesLastSectionPass(List<FaultSectionPrefData> rupture, List<IDPairing> pairings,
				List<Integer> junctionIndexes) {
			boolean passContinual = continualFilter.doesLastSectionPass(rupture, pairings, junctionIndexes);
			boolean passCleanup = cleanupFilter.doesLastSectionPass(rupture, pairings, junctionIndexes);
//			if (!junctionIndexes.isEmpty() && junctionIndexes.get(junctionIndexes.size()-1) == pairings.size()-1)
//				return cleanupFilter.doesLastSectionPass(rupture, pairings, junctionIndexes);
			int i = rupture.size()-1;
			boolean junction = i > 0 &&
					rupture.get(i).getParentSectionId() != rupture.get(i-1).getParentSectionId();
			return passContinual && (!junction || passCleanup);
		}

		@Override
		public boolean isContinueOnFaulure() {
			return true;
		}

		@Override
		public boolean isApplyJunctionsOnly() {
			return false;
		}
		
	}
	
	private static <E extends AbstractLaughTest> E removeByClass(List<AbstractLaughTest> tests, Class<E> clazz) {
		for (int i=tests.size(); --i>=0;) {
			AbstractLaughTest test = tests.get(i);
			if (clazz.isInstance(test))
				return (E)tests.remove(i);
		}
		throw new IllegalStateException();
	}
	
	static void calcCumulantMedianMag(FaultSystemSolution u3Sol, FaultSystemSolution rsSol,
			String catalogName, File outputDir) throws IOException, GMT_MapException, RuntimeException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Sect Index", "Sect Name", catalogName, "UCERF3");
		
		FaultSystemRupSet u3RupSet = u3Sol.getRupSet();
		FaultSystemRupSet rsRupSet = rsSol.getRupSet();
		
		Preconditions.checkState(u3RupSet.getNumSections() == rsRupSet.getNumSections());
		
		List<Double> rsVals = new ArrayList<>();
		List<Double> u3Vals = new ArrayList<>();
		
		for (int s=0; s<u3RupSet.getNumSections(); s++) {
			List<String> line = new ArrayList<>();
			line.add(s+"");
			line.add(u3RupSet.getFaultSectionData(s).getName());
			double rsVal = calcCumulantMedianMag(rsSol, s);
			rsVals.add(rsVal);
			line.add(rsVal+"");
			double u3Val = calcCumulantMedianMag(u3Sol, s);
			u3Vals.add(u3Val);
			line.add(u3Val+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, "mag_cumulant_medians.csv"));
		
		List<LocationList> faults = new ArrayList<>();
		for (FaultSectionPrefData fault : u3RupSet.getFaultSectionDataList())
			faults.add(fault.getFaultTrace());
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(6d,  8.5d);
		Region reg = new CaliforniaRegions.RELM_TESTING();
		FaultBasedMapGen.makeFaultPlot(cpt, faults, Doubles.toArray(u3Vals), reg, outputDir,
				"mag_cumulant_medians_ucerf3", false, false, "UCERF3 Mag Cumulant Median");
		FaultBasedMapGen.makeFaultPlot(cpt, faults, Doubles.toArray(rsVals), reg, outputDir,
				"mag_cumulant_medians_"+catalogName.toLowerCase(), false, false, catalogName+" Mag Cumulant Median");
		
//		CPT diffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1d, 1d);
		CPT diffCPT = new CPT(-1, 1d,
				new Color(0, 0, 140), new Color(0, 60, 200 ), new Color(0, 120, 255),
				Color.WHITE,
				new Color(255, 120, 0), new Color(200, 60, 0), new Color(140, 0, 0));
		diffCPT.setBelowMinColor(diffCPT.getMinColor());
		diffCPT.setAboveMaxColor(diffCPT.getMaxColor());
		double[] diffVals = new double[faults.size()];
		
		for (int i=0; i<diffVals.length; i++)
			diffVals[i] = rsVals.get(i) - u3Vals.get(i);
		
		FaultBasedMapGen.makeFaultPlot(diffCPT, faults, diffVals, reg, outputDir,
				"mag_cumulant_medians_diff", false, true, catalogName+"-U3 Mag Cumulant Median");
	}
	
	private static double calcCumulantMedianMag(FaultSystemSolution sol, int s) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(5d, 9d, (int)((9d-5d)/0.01) + 1);
		for (int r : sol.getRupSet().getRupturesForSection(s)) {
			double mag = sol.getRupSet().getMagForRup(r);
			int i = func.getClosestXIndex(mag);
			for (int x=i; x<func.size(); x++)
				func.add(x, MagUtils.magToMoment(mag));
		}
		if (func.calcSumOfY_Vals() == 0)
			return Double.NaN;
		func.scale(1d/func.getMaxY());
		return func.getFirstInterpolatedX(0.5);
	}
	
	static void plotConnectivity(FaultSystemRupSet rupSet, File outputDir, String prefix, String title)
			throws IOException, GMT_MapException, RuntimeException {
		List<HashSet<Integer>> clusters = new ArrayList<>();
		Map<Integer, Integer> sectIndexToClusterIndexMap = new HashMap<>();
		
		for (int s=0; s<rupSet.getNumSections(); s++)
			processClusterRecursive(rupSet, s, clusters.size(), clusters, sectIndexToClusterIndexMap);
		
		System.out.println("Detected "+clusters.size()+" clusters for "+prefix);
		
		CPT refCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		refCPT = refCPT.rescale(0d, 1d);
		// list of values for each discrete color, initially sorted from first color to last
		List<Double> colorValues = Lists.newArrayList();
//		for (CPTVal cptVal : refCPT)
//			colorValues.add((double)cptVal.start);
//		colorValues.add((double)refCPT.get(refCPT.size()-1).end);
		for (double v=0; v<=1d; v+=0.1)
			colorValues.add(v);
		// now sorted from last color to first
		Collections.reverse(colorValues);
		refCPT.setNanColor(Color.GRAY);
		
		// sort from smallest to largest
		List<Integer> sizes = Lists.newArrayList();
		for (HashSet<Integer> cluster : clusters)
			sizes.add(cluster.size());
		clusters = ComparablePairing.getSortedData(sizes, clusters);
		// now reverse, largest to smallest
		Collections.reverse(clusters);
		
		if (clusters.size() > colorValues.size())
			clusters = clusters.subList(0, colorValues.size());
		
		List<LocationList> faults = new ArrayList<>();
		List<Double> values = new ArrayList<>();
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			faults.add(rupSet.getFaultSectionData(s).getFaultTrace());
			double val = Double.NaN;
			for (int i=0; i<clusters.size(); i++) {
				if (clusters.get(i).contains(s))
					val = colorValues.get(i);
			}
			values.add(val);
		}
		
		Region reg = new CaliforniaRegions.RELM_TESTING();
		FaultBasedMapGen.makeFaultPlot(refCPT, faults, Doubles.toArray(values), reg, outputDir,
				prefix, false, false, title+" ("+colorValues.size()+" largest)");
	}
	
	private static void processClusterRecursive(FaultSystemRupSet rupSet, int sect, int clusterIndex, List<HashSet<Integer>> clusters,
			Map<Integer, Integer> sectIndexToClusterIndexMap) {
		if (sectIndexToClusterIndexMap.containsKey(sect))
			// we've already done this one
			return;
		if (clusters.size() == clusterIndex)
			clusters.add(new HashSet<>());
		clusters.get(clusterIndex).add(sect);
		sectIndexToClusterIndexMap.put(sect, clusterIndex);
		for (int r : rupSet.getRupturesForSection(sect)) {
			for (int sect2 : rupSet.getSectionsIndicesForRup(r)) {
				processClusterRecursive(rupSet, sect2, clusterIndex, clusters, sectIndexToClusterIndexMap);
			}
		}
	}

}
