package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.VectorComponent;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator;
import scratch.kevin.prvi25.LogicTreeLineIntegralCalc;

public class CombinedLogicTreeLineIntegralsPlot {

	public static void main(String[] args) throws IOException {
//		File crustalLIDir = new File(CRUSTAL_DIR, "line_integrals");
//		File subLIDir = new File(SUBDUCTION_DIR, "line_integrals");
		File crustalLIDir = new File(new File(INV_DIR, "2025_05_09-prvi25_crustal_branches-dmSample10x"), "line_integrals");
		File subLIDir = new File(new File(INV_DIR, "2025_05_09-prvi25_subduction_branches"), "line_integrals");
		File outputDir = new File(FIGURES_DIR, "logic_tree_line_integrals");
//		File crustalLIDir = new File(new File(INV_DIR, "2025_01_17-prvi25_crustal_branches-dmSample10x"), "line_integrals");
//		File subLIDir = new File(new File(INV_DIR, "2025_01_17-prvi25_subduction_branches"), "line_integrals");
//		File outputDir = new File(FIGURES_DIR, "logic_tree_line_integrals_old");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<DiscretizedFunc> crustalPerpFuncs = new ArrayList<>();
		List<DiscretizedFunc> crustalParallelFuncs = new ArrayList<>();
		List<DiscretizedFunc> crustalFullFuncs = new ArrayList<>();
		List<Double> crustalWeights = new ArrayList<>();
		
		System.out.println("Loading crustal");
		loadFuncs(crustalLIDir, crustalPerpFuncs, crustalParallelFuncs, crustalFullFuncs, crustalWeights);
		
		List<DiscretizedFunc> subPerpFuncs = new ArrayList<>();
		List<DiscretizedFunc> subParallelFuncs = new ArrayList<>();
		List<DiscretizedFunc> subFullFuncs = new ArrayList<>();
		List<Double> subWeights = new ArrayList<>();
		
		System.out.println("Loading subduction");
		loadFuncs(subLIDir, subPerpFuncs, subParallelFuncs, subFullFuncs, subWeights);
		
		List<DiscretizedFunc> combPerpFuncs = new ArrayList<>();
		List<DiscretizedFunc> combParallelFuncs = new ArrayList<>();
		List<DiscretizedFunc> combFullFuncs = new ArrayList<>();
		List<Double> combWeights = new ArrayList<>();
		System.out.println("Building combined functions");
		
		double[] xVals = new double[crustalParallelFuncs.get(0).size()];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = crustalParallelFuncs.get(0).getX(i);
		
		for (int i=0; i<crustalWeights.size(); i++) {
			DiscretizedFunc crustalPerpFunc = crustalPerpFuncs.get(i);
			DiscretizedFunc crustalParallelFunc = crustalParallelFuncs.get(i);
			double crustalWeight = crustalWeights.get(i);
			Preconditions.checkState(crustalPerpFunc.size() == xVals.length);
			
			for (int j=0; j<subWeights.size(); j++) {
				DiscretizedFunc subPerpFunc = subPerpFuncs.get(j);
				DiscretizedFunc subParallelFunc = subParallelFuncs.get(j);
				double subWeight = subWeights.get(j);
				Preconditions.checkState(subPerpFunc.size() == xVals.length);

				double weight = crustalWeight * subWeight;
				
				double[] perpYs = new double[xVals.length];
				double[] parallelYs = new double[xVals.length];
				double[] fullYs = new double[xVals.length];
				
				for (int k=0; k<xVals.length; k++) {
					perpYs[k] = crustalPerpFunc.getY(k) + subPerpFunc.getY(k);
					parallelYs[k] = crustalParallelFunc.getY(k) + subParallelFunc.getY(k);
					fullYs[k] = Math.sqrt(perpYs[k]*perpYs[k] + parallelYs[k]*parallelYs[k]);
				}

				combPerpFuncs.add(new LightFixedXFunc(xVals, perpYs));
				combParallelFuncs.add(new LightFixedXFunc(xVals, parallelYs));
				combFullFuncs.add(new LightFixedXFunc(xVals, fullYs));
				combWeights.add(weight);
			}
		}
		
		System.out.println("Combined model has "+crustalWeights.size()+" x "+subWeights.size()+" = "+combWeights.size()+" branches");
		
		DiscretizedFunc crustalMeanPerp = null;
		DiscretizedFunc crustalMeanParallel = null;
		DiscretizedFunc crustalMeanFull = null;
		DiscretizedFunc subMeanPerp = null;
		DiscretizedFunc subMeanParallel = null;
		DiscretizedFunc subMeanFull = null;
		
		for (int i=0; i<3; i++) {
			String prefix;
			List<DiscretizedFunc> perpFuncs;
			List<DiscretizedFunc> parallelFuncs;
			List<DiscretizedFunc> fullFuncs;
			List<Double> weights;
			Range fixedYRange;
			switch (i) {
			case 0:
				prefix = "crustal";
				perpFuncs = crustalPerpFuncs;
				parallelFuncs = crustalParallelFuncs;
				fullFuncs = crustalFullFuncs;
				weights = crustalWeights;
				fixedYRange = new Range(0d, 20d);
				break;
			case 1:
				prefix = "subduction";
				perpFuncs = subPerpFuncs;
				parallelFuncs = subParallelFuncs;
				fullFuncs = subFullFuncs;
				weights = subWeights;
				fixedYRange = new Range(0d, 10d);
				break;
			case 2:
				prefix = "combined";
				perpFuncs = combPerpFuncs;
				parallelFuncs = combParallelFuncs;
				fullFuncs = combFullFuncs;
				weights = combWeights;
				fixedYRange = new Range(0d, 20d);
				break;

			default:
				throw new IllegalStateException();
			}
			List<PlotSpec> plots = new ArrayList<>();
			List<Range> xRanges = List.of(new Range(xVals[0], xVals[xVals.length-1]));
			List<Range> yRanges = new ArrayList<>();
			for (VectorComponent comp : VectorComponent.values()) {
				String myPrefix = prefix+"_"+comp.name();
				System.out.println("Plotting "+myPrefix);
				List<DiscretizedFunc> funcs;
				List<DiscretizedFunc> extraFuncs = null;
				String meanName = i < 2 ? "Average" : "Combined Average";
				if (comp == VectorComponent.FULL_HORIZONTAL) {
					funcs = fullFuncs;
					if (i == 2)
						extraFuncs = List.of(crustalMeanFull, subMeanFull);
				} else if (comp == VectorComponent.PARALLEL) {
					funcs = parallelFuncs;
					if (i == 2)
						extraFuncs = List.of(crustalMeanParallel, subMeanParallel);
				} else if (comp == VectorComponent.PERPENDICULAR) {
					funcs = perpFuncs;
					if (i == 2)
						extraFuncs = List.of(crustalMeanPerp, subMeanPerp);
				} else {
					throw new IllegalStateException();
				}
				List<PlotCurveCharacterstics> extraChars = null;
				if (extraFuncs != null) {
					 extraChars = new ArrayList<>();
					extraFuncs.get(0).setName("Crustal");
					extraChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_blue));
					extraFuncs.get(1).setName("Subduction interface");
					extraChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_orange));
				}
				PlotSpec plot = LogicTreeLineIntegralCalc.writePlot(outputDir, myPrefix, funcs, weights, comp,
						meanName, extraFuncs, extraChars, comp == VectorComponent.FULL_HORIZONTAL ? fixedYRange : null);
				plot.setLegendVisible(plots.isEmpty());
				if (comp == VectorComponent.FULL_HORIZONTAL)
					plot.setYAxisLabel("|Full Horizontal| (mm/yr)");
				else if (comp == VectorComponent.PARALLEL)
					plot.setYAxisLabel("North-South (mm/yr)");
				else if (comp == VectorComponent.PERPENDICULAR)
					plot.setYAxisLabel("East-West (mm/yr)");
				plot.setYAxisLabel(plot.getYAxisLabel().replace("Summed Rate ", ""));
				plots.add(plot);
				yRanges.add(FaultSystemLineIntegralCalculator.getPlotYRange(plot));
				List<? extends PlotElement> origFuncs = plot.getPlotElems();
				DiscretizedFunc meanFunc = (DiscretizedFunc)origFuncs.get(origFuncs.size()-1);
				if (i == 0) {
					if (comp == VectorComponent.FULL_HORIZONTAL)
						crustalMeanFull = meanFunc;
					else if (comp == VectorComponent.PARALLEL)
						crustalMeanParallel = meanFunc;
					else if (comp == VectorComponent.PERPENDICULAR)
						crustalMeanPerp = meanFunc;
				} else if (i == 1) {
					if (comp == VectorComponent.FULL_HORIZONTAL)
						subMeanFull = meanFunc;
					else if (comp == VectorComponent.PARALLEL)
						subMeanParallel = meanFunc;
					else if (comp == VectorComponent.PERPENDICULAR)
						subMeanPerp = meanFunc;
				}
			}
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plots, false, false, xRanges, yRanges);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 1200, true, true, false);
		}
	}
	
	private static void loadFuncs(File dir, List<DiscretizedFunc> perpFuncs, List<DiscretizedFunc> parallelFuncs,
			List<DiscretizedFunc> fullFuncs, List<Double> weights)
			throws IOException {
		for (int c=0; c<3; c++) {
			VectorComponent component;
			List<DiscretizedFunc> funcs;
			switch (c) {
			case 0:
				component = VectorComponent.FULL_HORIZONTAL;
				funcs = fullFuncs;
				break;
			case 1:
				component = VectorComponent.PARALLEL;
				funcs = parallelFuncs;
				break;
			case 2:
				component = VectorComponent.PERPENDICULAR;
				funcs = perpFuncs;
				break;

			default:
				throw new IllegalStateException();
			}
			
			File file = new File(dir, "branch_line_integrals_"+component.name()+".csv");
			System.out.println("Reading "+file.getName());
			CSVFile<String> csv = CSVFile.readFile(file, true);
			
			double[] xVals = new double[csv.getNumCols()-2];
			
			for (int i=0; i<xVals.length; i++)
				xVals[i] = csv.getDouble(0, i+2);
			
			for (int row=1; row<csv.getNumRows(); row++) {
				Preconditions.checkState(csv.getInt(row, 0) == row-1);
				if (c == 0)
					weights.add(csv.getDouble(row, 1));
				double[] yVals = new double[xVals.length];
				for (int i=0; i<xVals.length; i++)
					yVals[i] = csv.getDouble(row, i+2);
				
				funcs.add(new LightFixedXFunc(xVals, yVals));
			}
		}
	}

}
