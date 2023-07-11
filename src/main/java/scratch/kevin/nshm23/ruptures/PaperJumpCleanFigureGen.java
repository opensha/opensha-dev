package scratch.kevin.nshm23.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.JumpAzimuthsPlot;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.RupHistogramPlots;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectMaxValuesPlot;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.RupHistogramPlots.HistScalar;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.RupHistogramPlots.HistScalarValues;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.RupHistogramPlots.RakeType;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.Table;
import com.google.common.io.Files;

import scratch.UCERF3.utils.U3FaultSystemIO;

public class PaperJumpCleanFigureGen {
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Documents/papers/2021_UCERF4_Plausibility/figures");
		
		String prefix = "figure_3";
		
		File rupSetFileA = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip");
//		File rupSetFileB = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5.zip");
		File rupSetFileB = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5_plus_coulomb_conns.zip");
		File rupSetFileC = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_plausibleMulti15km_adaptive6km_direct_"
				+ "cmlRake360_jumpP0.001_slipP0.05incrCapDist_cff0.75IntsPos_comb2Paths_cffFavP0.01_cffFavRatioN2P0.5_sectFractGrow0.1.zip");
		
		FaultSystemRupSet[] rupSets = {
				FaultSystemRupSet.load(rupSetFileA),
				FaultSystemRupSet.load(rupSetFileB),
				FaultSystemRupSet.load(rupSetFileC)
		};
		
		SectionDistanceAzimuthCalculator distAzCalc = null;
		List<Table<RakeType, RakeType, List<Double>>> rupSetJumpTables = new ArrayList<>();
		for (FaultSystemRupSet rupSet : rupSets) {
			if (distAzCalc == null)
				distAzCalc = rupSet.getModule(SectionDistanceAzimuthCalculator.class);
			if (distAzCalc == null && rupSet.hasModule(PlausibilityConfiguration.class))
				distAzCalc = rupSet.getModule(PlausibilityConfiguration.class).getDistAzCalc();
			if (distAzCalc == null)
				distAzCalc = new SectionDistanceAzimuthCalculator(rupSet.getFaultSectionDataList());
			if (!rupSet.hasModule(SectionDistanceAzimuthCalculator.class))
				rupSet.addModule(distAzCalc);
			// add connection search
			if (!rupSet.hasModule(RuptureConnectionSearch.class))
				rupSet.addModule(new RuptureConnectionSearch(rupSet, distAzCalc, 100d, false));
			
			if (!rupSet.hasAvailableModule(ClusterRuptures.class)) {
				if (rupSet.getNumRuptures() == 253706)
					rupSet.addModule(ClusterRuptures.singleStranged(rupSet));
				else
					rupSet.addAvailableModule(new Callable<ClusterRuptures>() {
						
						@Override
						public ClusterRuptures call() throws Exception {
							return ClusterRuptures.instance(rupSet, rupSet.requireModule(RuptureConnectionSearch.class));
						}
						
					}, ClusterRuptures.class);
			}
			rupSetJumpTables.add(JumpAzimuthsPlot.calcJumpAzimuths(rupSet));
		}
		
		plotAzimuths(dir, prefix+"a", "UCERF3 Jump Azimuths", false, rupSetJumpTables.get(0));
		plotAzimuths(dir, prefix+"b", "RSQSim Jump Azimuths", false, rupSetJumpTables.get(1));
		plotAzimuths(dir, prefix+"c", "Proposed Model Jump Azimuths", true, rupSetJumpTables.get(2));
		
		// now plot all with labels, e.g., for use int talks
		plotAzimuths(new File("/tmp"), "jump_az_u3", "UCERF3 Jump Azimuths", true, rupSetJumpTables.get(0));
		plotAzimuths(new File("/tmp"), "jump_az_rsqsim", "RSQSim Jump Azimuths", true, rupSetJumpTables.get(1));
		plotAzimuths(new File("/tmp"), "jump_az_coulomb", "Proposed Model Jump Azimuths", true, rupSetJumpTables.get(2));
		
		PlotPreferences prefs = PlotUtils.getDefaultFigurePrefs();
		prefs.setTickLabelFontSize(24);
		prefs.setAxisLabelFontSize(26);
		prefs.setLegendFontSize(24);
		HistScalarValues valsA = new HistScalarValues(HistScalar.LENGTH, rupSets[0], null,
				rupSets[0].requireModule(ClusterRuptures.class).getAll(), distAzCalc);
		HistScalarValues valsB = new HistScalarValues(HistScalar.LENGTH, rupSets[2], null,
				rupSets[2].requireModule(ClusterRuptures.class).getAll(), distAzCalc);
		
		prefix = "figure_14";
		HashSet<UniqueRupture> uniquesA = new HashSet<>();
		for (ClusterRupture rup : rupSets[0].requireModule(ClusterRuptures.class))
			uniquesA.add(rup.unique);
		HashSet<UniqueRupture> uniquesB = new HashSet<>();
		for (ClusterRupture rup : rupSets[2].requireModule(ClusterRuptures.class))
			uniquesB.add(rup.unique);
		RupHistogramPlots.PLOT_PREFS_DEFAULT = prefs;
		RupHistogramPlots.TITLES = false;
		RupHistogramPlots.plotRuptureHistogram(dir, prefix+"a", valsA, valsB, uniquesB, Color.BLUE, false, false);
		RupHistogramPlots.plotRuptureHistogram(dir, prefix+"b", valsB, valsA, uniquesA, Color.RED, false, false);
		
		prefix = "figure_15";
		GeographicMapMaker.PLOT_PREFS_DEFAULT = prefs;
		SectMaxValuesPlot.plotScalarMaxMapView(rupSets[0], dir, prefix+"a", " ", valsA, valsB,
				new CaliforniaRegions.RELM_TESTING(), Color.BLUE, false, false);
		SectMaxValuesPlot.plotScalarMaxMapView(rupSets[2], dir, prefix+"b", " ", valsB, valsA,
				new CaliforniaRegions.RELM_TESTING(), Color.RED, false, false);
		new File(dir, prefix+"a_hist.png").delete();
		new File(dir, prefix+"b_hist.png").delete();
//		plotScalarMaxMapView(rupSets[0], dir, prefix+"a", " ", valsA, new CaliforniaRegions.RELM_TESTING());
//		plotScalarMaxMapView(rupSets[2], dir, prefix+"b", " ", valsB, new CaliforniaRegions.RELM_TESTING());
		
		prefix = "figure_17";
		FaultSystemSolution u3Sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_03_24-u3_branches-FM3_1-2000ip/results_FM3_1_branch_averaged.zip"));
		FaultSystemSolution newSol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_03_25-u3_branches-coulomb-FM3_1-2000ip/results_FM3_1_branch_averaged.zip"));
		
		plotJumpsOverDistance(dir, prefix, u3Sol, newSol);
	}
	
	private static final RakeType sourceRake = null;
	
	private static List<RakeType> destRakes = new ArrayList<>();
	static {
		destRakes.add(null);
		for (RakeType type : RakeType.values())
			destRakes.add(type);
	}
	
	private static void plotAzimuths(File outputDir, String prefix, String title, boolean labels,
			Table<RakeType, RakeType, List<Double>> azTable) throws IOException {
		Map<RakeType, List<Double>> azMap = JumpAzimuthsPlot.getAzimuthsFrom(sourceRake, azTable);
		
		Range xRange = new Range(-180d, 180d);
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(xRange);
		
		List<Range> yRanges = new ArrayList<>();
		List<PlotSpec> specs = new ArrayList<>();
		
		for (int i=0; i<destRakes.size(); i++) {
			RakeType destRake = destRakes.get(i);
			
			HistogramFunction hist = HistogramFunction.getEncompassingHistogram(-179d, 179d, 15d);
			for (RakeType oRake : azMap.keySet()) {
				if (destRake != null && destRake != oRake)
					continue;
				for (double azDiff : azMap.get(oRake)) {
					hist.add(hist.getClosestXIndex(azDiff), 1d);
				}
			}

			Color color;
			String label;
			if (destRake == null) {
				color = Color.DARK_GRAY;
				label = "Any";
			} else {
				color = destRake.color;
				label = destRake.name;
				label = label.replace("Right", "R").replace("Left", "L");
				label = label.replace("Lateral", "L").replace("SS", "S-S");
			}
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
			
			double maxY = Math.max(1.1*hist.getMaxY(), 1d);
			Range yRange = new Range(0d, maxY);
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Azimuthal Difference", "Count");
			
			if (labels) {
				XYTextAnnotation ann = new XYTextAnnotation("To "+label, 175, maxY*0.975);
				ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 26));
				ann.setTextAnchor(TextAnchor.TOP_RIGHT);
				spec.addPlotAnnotation(ann);
			}
			
			specs.add(spec);
			yRanges.add(yRange);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(24);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(26);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		CombinedDomainXYPlot plot = (CombinedDomainXYPlot)gp.getPlot();
		for (XYPlot subPlot : (List<XYPlot>)plot.getSubplots())
			((NumberAxis)subPlot.getRangeAxis()).setNumberFormatOverride(new DecimalFormat("0E0"));
		
		PlotUtils.setTick(plot.getDomainAxis(), 90d);
		
		File pngFile = new File(outputDir, prefix+".png");
		File pdfFile = new File(outputDir, prefix+".pdf");
		gp.getChartPanel().setSize(650, 1000);
		gp.saveAsPNG(pngFile.getAbsolutePath());
		gp.saveAsPDF(pdfFile.getAbsolutePath());
	}
	
	private static void plotJumpsOverDistance(File outputDir, String prefix,
			FaultSystemSolution u3Sol, FaultSystemSolution newSol) throws IOException {
		// Table 3 from Biasi & Wesnousky (2016)
		DiscretizedFunc bw2016 = new ArbitrarilyDiscretizedFunc();
		
		double totNumBW = 76;
		bw2016.set(0d, 26d/totNumBW);
		bw2016.set(1d, 29d/totNumBW);
		bw2016.set(2d, 10d/totNumBW);
		bw2016.set(3d, 6d/totNumBW);
		bw2016.set(4d, 2d/totNumBW);
		bw2016.set(5d, 2d/totNumBW);
		bw2016.set(6d, 0d);
		bw2016.set(7d, 0d);
		bw2016.set(8d, 0d);
		bw2016.set(9d, 1d/totNumBW);
		
		DiscretizedFunc u3Func = calcJumpDistFunc(u3Sol, 0d, 1f);
		u3Func.setName("UCERF3");
		DiscretizedFunc newFunc = calcJumpDistFunc(newSol, 0d, 1f);
		newFunc.setName("Proposed Model w/ UCERF3 Inversion");
		bw2016.setName("Biasi & Wesnousky (2016), Table 3");
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(bw2016);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 8f, Color.GREEN.darker()));
		
		funcs.add(u3Func);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_TRIANGLE, 8f, Color.BLUE.darker()));
		
		funcs.add(newFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_SQUARE, 8f, Color.RED.darker()));
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setLegendFontSize(20);
		gp.setAxisLabelFontSize(26);
		gp.setTickLabelFontSize(24);
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Number of Jumps ≥1 km", "Fraction of Supra-Seis Ruptures");
		spec.setLegendInset(true);
		
		gp.drawGraphPanel(spec, false, false, new Range(0d, 5d), new Range(0d, 1d));
		
		PlotUtils.setTick(gp.getPlot().getDomainAxis(), 1d);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 650, true, true, false);
	}
	
	private static DiscretizedFunc calcJumpDistFunc(FaultSystemSolution sol, double minMag, float jumpDist) {
		DiscretizedFunc solFunc = new ArbitrarilyDiscretizedFunc();
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		ClusterRuptures clusterRups = rupSet.requireModule(ClusterRuptures.class);

		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			double mag = rupSet.getMagForRup(r);

			if (mag < minMag)
				continue;
			
			ClusterRupture rup = clusterRups.get(r);
			int jumpsOverDist = 0;
			for (Jump jump : rup.getJumpsIterable()) {
				if ((float)jump.distance >= jumpDist)
					jumpsOverDist++;
			}

			double rate = sol.getRateForRup(r);
			
			while (jumpsOverDist >= solFunc.size())
				solFunc.set((double)solFunc.size(), 0d);
			
			solFunc.set(jumpsOverDist, solFunc.getY(jumpsOverDist) + rate);
		}
		
		solFunc.scale(1d/sol.getTotalRateForAllFaultSystemRups());
		
		return solFunc;
	}

}
