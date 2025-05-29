package scratch.kevin.prvi25;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.LineIntegralResult;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.VectorComponent;
import scratch.kevin.prvi25.figures.SlipRateFigures;

public class LogicTreeLineIntegralCalc {
	
	private static double minLat = SlipRateFigures.CRUSTAL_FAULT_MAP_REG.getMinLat();
	private static double maxLat = SlipRateFigures.CRUSTAL_FAULT_MAP_REG.getMaxLat();
	private static double minLon = SlipRateFigures.CRUSTAL_FAULT_MAP_REG.getMinLon();
	private static double maxLon = SlipRateFigures.CRUSTAL_FAULT_MAP_REG.getMaxLon();
	private static double delta = 0.1d;

	public static void main(String[] args) throws IOException {
		System.setProperty("java.awt.headless", "true");
		Preconditions.checkState(args.length == 2 || args.length == 3,
				"USAGE: <solution logic tree file> <output dir>");
		File sltFile = new File(args[0]);
		SolutionLogicTree slt = SolutionLogicTree.load(sltFile);
		File outputDir = new File(args[1]);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdirs(),
				"Couldn't create output dir: "+outputDir.getAbsolutePath());
		
		List<Location> startLocs = new ArrayList<>();
		List<Location> endLocs = new ArrayList<>();
		for (double lon=minLon; (float)lon<=(float)maxLon; lon += delta) {
			startLocs.add(new Location(minLat, lon));
			endLocs.add(new Location(maxLat, lon));
		}
		
		VectorComponent[] components = VectorComponent.values();
		
		FaultSystemLineIntegralCalculator calc = null;
		FaultSystemRupSet prevRupSet = null;
		ExecutorService exec = Executors.newFixedThreadPool(Integer.max(5, FaultSysTools.defaultNumThreads()));
		List<Future<DiscretizedFunc[]>> futures = new ArrayList<>();
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		
		LogicTree<?> tree = slt.getLogicTree();
		for (int i=0; i<tree.size(); i++) {
			System.out.println("Processing branch "+i+"/"+tree.size());
			LogicTreeBranch<?> branch = tree.getBranch(i);
			
			FaultSystemSolution sol = slt.forBranch(branch, true);
			if (calc == null || prevRupSet != sol.getRupSet()) {
				prevRupSet = sol.getRupSet();
				calc = new FaultSystemLineIntegralCalculator(prevRupSet);
			}
			
			FaultSystemLineIntegralCalculator myCalc = calc;
			
			futures.add(exec.submit(new Callable<DiscretizedFunc[]>() {
				
				@Override
				public DiscretizedFunc[] call() {
					List<LineIntegralResult> lis = new ArrayList<>(startLocs.size());
					for (int l=0; l<startLocs.size(); l++)
						lis.add(myCalc.calcSolutionLineIntegral(sol, startLocs.get(l), endLocs.get(l), false));
					
					DiscretizedFunc[] ret = new DiscretizedFunc[components.length];
					for (int c=0; c<components.length; c++)
						ret[c] = myCalc.buildIntegralFunction(false, lis, components[c]);
					return ret;
				}
			}));
			
		}
		
		List<List<DiscretizedFunc>> componentBranchResults = new ArrayList<>();
		for (int c=0; c<components.length; c++)
			componentBranchResults.add(new ArrayList<>());
		List<Double> branchWeights = new ArrayList<>();
		
		System.out.println("Waiting on any remaining calc futures");
		for (int i=0; i<tree.size(); i++) {
			double weight = tree.getBranchWeight(i);
			Future<DiscretizedFunc[]> future = futures.get(i);
			try {
				DiscretizedFunc[] funcs = future.get();
				for (int c=0; c<components.length; c++)
					componentBranchResults.get(c).add(funcs[c]);
				branchWeights.add(weight);
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		for (int c=0; c<components.length; c++) {
			System.out.println("Writing results for "+components[c]);
			List<DiscretizedFunc> branchResults = componentBranchResults.get(c);
			
			System.out.println("Writing CSV");
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>();
			header.add("Branch Index");
			header.add("Branch Weight");
			for (Location loc : startLocs)
				header.add((float)loc.lon+"");
			csv.addLine(header);
			for (int b=0; b<branchResults.size(); b++) {
				List<String> line = new ArrayList<>(header.size());
				line.add(b+"");
				line.add(branchWeights.get(b)+"");
				for (Point2D pt : branchResults.get(b))
					line.add(pt.getY()+"");
				csv.addLine(line);
			}
			
			String prefix = "branch_line_integrals_"+components[c].name();
			csv.writeToFile(new File(outputDir, prefix+".csv"));
			
			System.out.println("Building plot"); 
			writePlot(outputDir, prefix, branchResults, branchWeights, components[c]);
		}
		
		exec.shutdown();
	}
	
	public static PlotSpec writePlot(File outputDir, String prefix, List<DiscretizedFunc> branchResults,
			List<Double> branchWeights, VectorComponent component) throws IOException {
		return writePlot(outputDir, prefix, branchResults, branchWeights, component, "Average", null, null, null);
	}
	
	public static PlotSpec writePlot(File outputDir, String prefix, List<DiscretizedFunc> branchResults,
			List<Double> branchWeights, VectorComponent component, String meanName, List<DiscretizedFunc> extraFuncs,
			List<PlotCurveCharacterstics> extraChars, Range yRange) throws IOException {
		int numLocs = branchResults.get(0).size();
		List<LightFixedXFunc> normCDFs = new ArrayList<>(numLocs);
		DiscretizedFunc mean = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<numLocs; i++) {
			double weightedSum = 0d;
			double weightSum = 0d;
			List<Double> vals = new ArrayList<>(branchResults.size());
			for (int b=0; b<branchResults.size(); b++) {
				double val = branchResults.get(b).getY(i);
				double weight = branchWeights.get(b);
				weightedSum += weight*val;
				weightSum += weight;
				vals.add(val);
			}
			
			normCDFs.add(ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(vals, branchWeights));
			mean.set(branchResults.get(0).getX(i), weightedSum/weightSum);
		}
		
		double[] fractiles = {0d, 0.025, 0.16, 0.84, 0.975d, 1d};
		String fractileLabel = "p[0,2.5,16,84,97.5,100]";
		
		ArbitrarilyDiscretizedFunc[] fractileFuncs = new ArbitrarilyDiscretizedFunc[fractiles.length];
		for (int f=0; f<fractileFuncs.length; f++) {
			fractileFuncs[f] = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<normCDFs.size(); i++) {
				LightFixedXFunc normCDF = normCDFs.get(i);
				double x = mean.getX(i);
				double y;
				if (normCDF.size() == 1)
					y = 0;
				else
					y = ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(normCDF, fractiles[f]);
				fractileFuncs[f].set(x, y);
			}
		}
		
		int cnt = 0;
		ArbitrarilyDiscretizedFunc incrMin = fractileFuncs[cnt++];
		ArbitrarilyDiscretizedFunc incrP025 = fractileFuncs[cnt++];
		ArbitrarilyDiscretizedFunc incrP16 = fractileFuncs[cnt++];
		ArbitrarilyDiscretizedFunc incrP84 = fractileFuncs[cnt++];
		ArbitrarilyDiscretizedFunc incrP975 = fractileFuncs[cnt++];
		ArbitrarilyDiscretizedFunc incrMax = fractileFuncs[cnt++];
		UncertainArbDiscFunc bounds = new UncertainArbDiscFunc(
				getAvg(incrMin, incrMax), incrMin, incrMax, null);
		UncertainArbDiscFunc bounds95 = new UncertainArbDiscFunc(
				getAvg(incrP025, incrP975), incrP025, incrP975, null);
		UncertainArbDiscFunc bounds68 = new UncertainArbDiscFunc(
				getAvg(incrP16, incrP84), incrP16, incrP84, null);
		bounds.setName(fractileLabel);
		bounds95.setName(null);
		bounds68.setName(null);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Color distColor = Colors.tab_grey;
		Color meanColor = Color.BLACK;
		Color transColor = new Color(distColor.getRed(), distColor.getGreen(), distColor.getBlue(), 60);
		PlotCurveCharacterstics minMaxChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, transColor);
		
		funcs.add(bounds);
		chars.add(minMaxChar);
		
		funcs.add(bounds95);
		chars.add(minMaxChar);
		
		funcs.add(bounds68);
		chars.add(minMaxChar);
		
		if (extraFuncs != null) {
			funcs.addAll(extraFuncs);
			chars.addAll(extraChars);
		}
		
		mean.setName(meanName);
		funcs.add(mean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, meanColor));
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Longitude",
				component.label+" Summed Rate (mm/yr)");
//		plot.setLegendInset(RectangleAnchor.TOP_RIGHT, 0.975, 0.975, 0.9, false);
//		plot.setLegendInset(RectangleAnchor.TOP_LEFT);
		plot.setLegendVisible(true);
		
		if (yRange == null)
			yRange = FaultSystemLineIntegralCalculator.getPlotYRange(plot);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		int width = 900;
		int height = 650;
		
		Range xRange = new Range(minLon, maxLon);
		
		gp.drawGraphPanel(plot, false, false, xRange, yRange);

		PlotUtils.writePlots(outputDir, prefix, gp, width, height, true, true, false);
		
		return plot;
	}
	
	private static DiscretizedFunc getAvg(DiscretizedFunc min, DiscretizedFunc max) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<min.size(); i++)
			ret.set(min.getX(i), 0.5*(min.getY(i) + max.getY(i)));
		return ret;
	}
}
