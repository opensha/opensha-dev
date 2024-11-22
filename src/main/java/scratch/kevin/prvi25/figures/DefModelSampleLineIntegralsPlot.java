package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.LineIntegralResult;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.VectorComponent;

public class DefModelSampleLineIntegralsPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "prvi_dm_sample_line_integrals";
		
		LogicTree<?> randTree = LogicTree.read(new File(
				"/home/kevin/OpenSHA/nshm23/batch_inversions/2024_10_24-prvi25_crustal_branches-dmSample5x/logic_tree.json"));
		
		PRVI25_CrustalFaultModels fm = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1;
		VectorComponent[] components = VectorComponent.values();
		
		double[] fractiles = {0d, 0.025, 0.16, 0.84, 0.975d, 1d};
		String fractileLabel = "Sample p[0,2.5,16,84,97.5,100]";
		
		boolean plotOrig = true;
		
		CPT randColorCPT = GMT_CPT_Files.CATEGORICAL_TAB10.instance();
		int numIndv = Integer.min(randColorCPT.size(), randTree.size());

		List<? extends FaultSection> avgRates = PRVI25_CrustalDeformationModels.GEOLOGIC_DIST_AVG.build(fm);
		List<? extends FaultSection> origRates = PRVI25_CrustalDeformationModels.GEOLOGIC.build(fm);
		List<? extends FaultSection> upperRates = PRVI25_CrustalDeformationModels.GEOLOGIC_HIGH.build(fm);
		List<? extends FaultSection> lowerRates = PRVI25_CrustalDeformationModels.GEOLOGIC_LOW.build(fm);
		
		DiscretizedFunc[] avgIntegrals = calcSummedLineIntegrals(avgRates, components);
		for (int i=1; i<avgIntegrals.length; i++)
			Preconditions.checkState(avgIntegrals[0].size() == avgIntegrals[i].size());
		DiscretizedFunc[] origIntegrals = plotOrig ? calcSummedLineIntegrals(origRates, components) : null;
		DiscretizedFunc[] upperIntegrals = calcSummedLineIntegrals(upperRates, components);
		DiscretizedFunc[] lowerIntegrals = calcSummedLineIntegrals(lowerRates, components);
		
		List<DiscretizedFunc[]> randIntegralsList = new ArrayList<>(randTree.size());
		for (LogicTreeBranch<?> branch : randTree) {
			RupSetDeformationModel dm = branch.requireValue(RupSetDeformationModel.class);
			System.out.println("Processing "+dm);
			List<? extends FaultSection> branchSects = dm.build(fm);
			DiscretizedFunc[] randIntegrals = calcSummedLineIntegrals(branchSects, components);
			randIntegralsList.add(randIntegrals);
			Preconditions.checkState(randIntegrals[0].size() == avgIntegrals[0].size());
		}
		
//		Random rand = new Random(numIndv * randIntegralsList.size());
		Random rand = new Random(123456l);
		List<Integer> randIndexes = new ArrayList<>();
		while (randIndexes.size() < numIndv) {
			int index = rand.nextInt(randIntegralsList.size());
			if (!randIndexes.contains(index))
				randIndexes.add(index);
		}
		
		for (int c=0; c<components.length; c++) {
			System.out.println("Building distributions for "+components[c]);
			
			DiscretizedFunc avgIntegral = avgIntegrals[c];
			
			// build cdfs
			LightFixedXFunc[] cdfs = new LightFixedXFunc[avgIntegrals[0].size()];
			for (int x=0; x<cdfs.length; x++) {
				double xVal = avgIntegral.getX(x);
				double[] values = new double[randIntegralsList.size()];
				for (int i=0; i<values.length; i++) {
					DiscretizedFunc randIntegral = randIntegralsList.get(i)[c];
					Preconditions.checkState((float)xVal == (float)randIntegral.getX(x));
					values[i] = randIntegral.getY(x);
				}
				cdfs[x] = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(values, null);
			}
			
			ArbitrarilyDiscretizedFunc[] fractileFuncs = new ArbitrarilyDiscretizedFunc[fractiles.length];
			for (int f=0; f<fractileFuncs.length; f++) {
				fractileFuncs[f] = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<cdfs.length; i++) {
					double x = avgIntegral.getX(i);
					double y;
					if ((float)fractiles[f] <= (float)cdfs[i].getMinY())
						y = cdfs[i].getMinX();
					else if (fractiles[f] == 1)
						y = cdfs[i].getMaxX();
					else
						y = cdfs[i].getFirstInterpolatedX(fractiles[f]);
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
			
			for (boolean individual : new boolean[] {false,true}) {
				if (individual && numIndv <= 0)
					continue;
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				Color lowerColor = Colors.tab_blue;
				Color upperColor = Colors.tab_orange;
//				Color distColor = Colors.tab_red;
				Color distColor = Colors.tab_grey;
				Color origColor = Colors.tab_green;
				Color transColor = new Color(distColor.getRed(), distColor.getGreen(), distColor.getBlue(), 30);
				PlotCurveCharacterstics minMaxChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, transColor);
				
				funcs.add(bounds);
				chars.add(minMaxChar);
				
				funcs.add(bounds95);
				chars.add(minMaxChar);
				
				funcs.add(bounds68);
				chars.add(minMaxChar);
				
				String plotPrefix = prefix+"_"+components[c].name();
				if (individual) {
					lowerIntegrals[c].setName("Original Extrema");
					funcs.add(lowerIntegrals[c]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					
					upperIntegrals[c].setName(null);
					funcs.add(upperIntegrals[c]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					
					plotPrefix += "_indv";
					for (int i=0; i<randIndexes.size(); i++) {
						DiscretizedFunc func = randIntegralsList.get(randIndexes.get(i))[c];
						func.setName(i == 0 ? randIndexes.size()+" Random Examples" : null);
						funcs.add(func);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f,
								randColorCPT == null ? Color.BLACK : randColorCPT.get(i % randColorCPT.size()).minColor));
					}
				} else {
					lowerIntegrals[c].setName("All Low");
					funcs.add(lowerIntegrals[c]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, lowerColor));
					
					if (plotOrig) {
						origIntegrals[c].setName("Preferred");
						funcs.add(origIntegrals[c]);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, origColor));
					}
					
					avgIntegral.setName("Average");
					funcs.add(avgIntegral);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					
					upperIntegrals[c].setName("All High");
					funcs.add(upperIntegrals[c]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, upperColor));
				}
				
				PlotSpec plot = new PlotSpec(funcs, chars, "Deformation Model Line Integrals", "Longitude",
						components[c].label+" Summed Rate (mm/yr)");
				plot.setLegendInset(RectangleAnchor.TOP_RIGHT, 0.975, 0.975, 0.9, false);
				
				Range yRange = FaultSystemLineIntegralCalculator.getPlotYRange(plot);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				int width = 900;
				int height = 650;
				
				Range xRange = new Range(minLon, maxLon);
				
				gp.drawGraphPanel(plot, false, false, xRange, yRange);

				PlotUtils.writePlots(outputDir, plotPrefix, gp, width, height, true, true, false);
				
				FaultSystemLineIntegralCalculator calc = new FaultSystemLineIntegralCalculator(avgRates, true);
				
				List<LineIntegralResult> mapIntegrals = calcLineIntegrals(calc, 1);
				
				GeographicMapMaker mapMaker = calc.buildMapPlot(false, false, 10d, mapIntegrals);
				mapMaker.setCPTLocation(RectangleEdge.TOP);
				PlotSpec map = mapMaker.buildPlot(plot.getTitle());
				
				gp.drawGraphPanel(List.of(map, plot), false, false, List.of(xRange), List.of(mapMaker.getYRange(), yRange));

				PlotUtils.writePlots(outputDir, plotPrefix, gp, width, true, true, true, false);
				
				double bigScalar = 1.5;
				gp.getPlotPrefs().scaleFontSizes(bigScalar);
				width = (int)(width*bigScalar);
				height = (int)(height*bigScalar);
				
				for (PlotCurveCharacterstics pChar : chars)
					pChar.setLineWidth((float)(pChar.getLineWidth()*bigScalar));
				
				gp.drawGraphPanel(plot, false, false, xRange, yRange);
				
				PlotUtils.writePlots(outputDir, plotPrefix+"_big", gp, width, height, true, true, false);
			}
		}
	}
	
	private static DiscretizedFunc getAvg(DiscretizedFunc min, DiscretizedFunc max) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<min.size(); i++)
			ret.set(min.getX(i), 0.5*(min.getY(i) + max.getY(i)));
		return ret;
	}
	
	private static double minLat = 16;
	private static double maxLat = 21;
	private static double minLon = -71d;
	private static double maxLon = -60.5d;
	
	private static List<LineIntegralResult> calcLineIntegrals(FaultSystemLineIntegralCalculator calc, double delta) {
		List<LineIntegralResult> integrals = new ArrayList<>();
		
		double minLon = DefModelSampleLineIntegralsPlot.minLon;
		if (delta == 1d)
			minLon += 1d;
		
		for (double lon=minLon; (float)lon<=(float)maxLon; lon += delta)
			integrals.add(calc.calcLineIntegral(new Location(minLat, lon), new Location(maxLat, lon)));
		
		return integrals;
	}
	
	private static DiscretizedFunc[] calcSummedLineIntegrals(List<? extends FaultSection> sects, VectorComponent[] components) {
		
		FaultSystemLineIntegralCalculator calc = new FaultSystemLineIntegralCalculator(sects, true);
		
		List<LineIntegralResult> integrals = calcLineIntegrals(calc, 0.1);
		
		DiscretizedFunc[] ret = new DiscretizedFunc[components.length];
		for (int c=0; c<ret.length; c++)
			ret[c] = calc.buildIntegralFunction(false, integrals, components[c]);
		return ret;
	}

}
