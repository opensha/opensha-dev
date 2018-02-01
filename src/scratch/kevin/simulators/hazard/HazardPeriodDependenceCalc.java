package scratch.kevin.simulators.hazard;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.util.MarkdownUtils;

public class HazardPeriodDependenceCalc {

	public static void main(String[] args) throws Exception {
		File hazardJobDir = new File("/home/kevin/Simulators/hazard/");
		File catalogsBaseDir = new File("/data/kevin/simulators/catalogs");
		File mainOutputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2388.instance(catalogsBaseDir);
		File catAllIMTsDir = new File(hazardJobDir, "2018_01_05-bruce2388-m6.5-sectArea0.2-skip5000yr-all_imts-8xPoints/curves_all_imts");
		
		String catalogName = "RSQSim";
		
		File u3AllIMTsDir = new File(hazardJobDir, "ucerf-comparisons/ucerf3-m6.5-8xPoints/curves_all_imts");
		
		double[] rps = { 2500 };
		
		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param saParam = (SA_Param) gmpe.getIntensityMeasure();
		List<Double> allPeriods = saParam.getPeriodParam().getAllowedDoubles();
		Collections.sort(allPeriods);
		allPeriods.add(0, 0d); // PGA is first
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		HashSet<Integer> nehrpGridIndexes = new HashSet<>();
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
			Location loc = city.getSite().getLocation();
			nehrpGridIndexes.add(gridReg.indexForLocation(loc));
		}
		
		ArbitrarilyDiscretizedFunc[] meanFuncs = new ArbitrarilyDiscretizedFunc[rps.length];
		ArbitrarilyDiscretizedFunc[] stdDevFuncs = new ArbitrarilyDiscretizedFunc[rps.length];
		ArbitrarilyDiscretizedFunc[] nehrpMeanFuncs = new ArbitrarilyDiscretizedFunc[rps.length];
		ArbitrarilyDiscretizedFunc[] medianFuncs = new ArbitrarilyDiscretizedFunc[rps.length];
		ArbitrarilyDiscretizedFunc[] madFuncs = new ArbitrarilyDiscretizedFunc[rps.length];
		for (int i=0; i<rps.length; i++) {
			meanFuncs[i] = new ArbitrarilyDiscretizedFunc();
			stdDevFuncs[i] = new ArbitrarilyDiscretizedFunc();
			nehrpMeanFuncs[i] = new ArbitrarilyDiscretizedFunc();
			medianFuncs[i] = new ArbitrarilyDiscretizedFunc();
			madFuncs[i] = new ArbitrarilyDiscretizedFunc();
		}
		
		for (int i=0; i<allPeriods.size(); i++) {
			String binFileName = "imrs"+(i+1)+".bin";
			double period = allPeriods.get(i);
			
			File simBinFile = new File(catAllIMTsDir, binFileName);
			Preconditions.checkState(catAllIMTsDir.exists());
			
			System.out.print("Loading simulation bin file for t="+(float)period+"...");
			BinaryHazardCurveReader reader = new BinaryHazardCurveReader(simBinFile.getAbsolutePath());
			Map<Location, ArbitrarilyDiscretizedFunc> simCurves = reader.getCurveMap();
			System.out.println("DONE");
			
			File u3BinFile = new File(u3AllIMTsDir, binFileName);
			Preconditions.checkState(u3AllIMTsDir.exists());
			
			System.out.print("Loading U3 bin file for t="+(float)period+"...");
			reader = new BinaryHazardCurveReader(u3BinFile.getAbsolutePath());
			Map<Location, ArbitrarilyDiscretizedFunc> u3Curves = reader.getCurveMap();
			System.out.println("DONE");
			
			List<List<Double>> ratioVals = new ArrayList<>();
			List<List<Double>> nehrpRatios = new ArrayList<>();
			for (int r=0; r<rps.length; r++) {
				ratioVals.add(new ArrayList<>());
				nehrpRatios.add(new ArrayList<>());
			}
			
			System.out.println("Calculating for each site");
			for (int s=0; s<gridReg.getNodeCount(); s++) {
				Location loc = gridReg.getLocation(s);
				ArbitrarilyDiscretizedFunc simCurve = simCurves.get(loc);
				ArbitrarilyDiscretizedFunc u3Curve = u3Curves.get(loc);
				
				boolean nehrp = nehrpGridIndexes.contains(s);
				
				for (int r=0; r<rps.length; r++) {
					double prob = 1d/rps[r];
					double simVal = HazardDataSetLoader.getCurveVal(simCurve, false, prob);
					double u3Val = HazardDataSetLoader.getCurveVal(u3Curve, false, prob);
					
					if (simVal == 0 || u3Val == 0 || Double.isNaN(simVal) || Double.isNaN(u3Val))
						continue;
					
					double ratio = Math.log(simVal/u3Val);
					Preconditions.checkState(Double.isFinite(ratio), "Bad ratio: Ln(%s/%s) = %s", simVal, u3Val, ratio);
					
					ratioVals.get(r).add(ratio);
					if (nehrp)
						nehrpRatios.get(r).add(ratio);
				}
			}
			
			for (int r=0; r<rps.length; r++) {
				double[] ratioArray = Doubles.toArray(ratioVals.get(r));
				double[] nehrpRatioArray = Doubles.toArray(nehrpRatios.get(r));
				
				// calc MAD
				double mad, median;
				if (ratioVals.isEmpty()) {
					mad = Double.NaN;
					median = Double.NaN;
				} else {
					median = DataUtils.median(ratioArray);
					double[] deviations = new double[ratioArray.length];
					for (int j=0; j<ratioArray.length; j++)
						deviations[j] = Math.abs(ratioArray[j] - median);
					mad = DataUtils.median(deviations);
				}
				
				double mean = StatUtils.mean(ratioArray);
				double stdDev = Math.sqrt(StatUtils.variance(ratioArray));
				
				meanFuncs[r].set(period, mean);
				stdDevFuncs[r].set(period, stdDev);
				medianFuncs[r].set(period, median);
				madFuncs[r].set(period, mad);
				
				double nehrpMean = StatUtils.mean(nehrpRatioArray);
				nehrpMeanFuncs[r].set(period, nehrpMean);
			}
			System.out.println("Done with t="+(float)period);
		}
		
		System.out.println("Building plots");
		
		File catOutDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catOutDir.exists() || catOutDir.mkdir());
		File catHazardOutDir = new File(catOutDir, "hazard_t_dependence");
		Preconditions.checkState(catHazardOutDir.exists() || catHazardOutDir.mkdir());
		File resourcesDir = new File(catHazardOutDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		// header
		lines.add("# Hazard T-Dependence");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		
		for (int r=0; r<rps.length; r++) {
			lines.add("## "+rpDF.format(rps[r])+"yr");
			lines.add("");
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			nehrpMeanFuncs[r].setName("NEHRP City Mean");
			funcs.add(nehrpMeanFuncs[r]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
			
			funcs.add(meanFuncs[r]);
			meanFuncs[r].setName("Mean");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
			
			funcs.add(stdDevFuncs[r]);
			stdDevFuncs[r].setName("StdDev");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
			
			funcs.add(medianFuncs[r]);
			medianFuncs[r].setName("Median");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.MAGENTA));
			
			funcs.add(madFuncs[r]);
			madFuncs[r].setName("MAD");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.ORANGE.darker()));
			
			PlotSpec spec = new PlotSpec(funcs, chars, rpDF.format(rps[r])+" T Dependence", "T [sec]", "Ln("+catalogName+"/UCERF3)");
			spec.setLegendVisible(true);
			
			PlotPreferences plotPrefs = PlotPreferences.getDefault();
			plotPrefs.setTickLabelFontSize(18);
			plotPrefs.setAxisLabelFontSize(20);
			plotPrefs.setPlotLabelFontSize(21);
			plotPrefs.setBackgroundColor(Color.WHITE);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			gp.drawGraphPanel(spec, true, false, new Range(0.01, 10d), new Range(-0.4, 0.4));
			gp.getChartPanel().setSize(800, 600);
			String prefix = rpDF.format(rps[r])+"yr";
			gp.saveAsTXT(new File(resourcesDir, prefix+".txt").getAbsolutePath());
			gp.saveAsPNG(new File(resourcesDir, prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(resourcesDir, prefix+".pdf").getAbsolutePath());
			
			lines.add("![prefix]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
		}
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, catHazardOutDir);
		
		catalog.writeMarkdownSummary(catOutDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
		
		System.out.println("DONE");
	}
	
	private static final DecimalFormat rpDF = new DecimalFormat("0.#");

}
