package scratch.kevin.nshm23.timJonEpistemic;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

public class CurveDebugPlot {
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/data/kevin/nshm23/batch_inversions/2023_08_23-tim_jon_nshm23_v3_epistemic_calcs_downsampled");
		File resultsDir = new File(dir, "results");
		String prefix = "curves_Berkeley_sa_0.2s";
		File outputDir = new File(dir, "debug");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		AttenRelRef[] gmpeRefs = {
				AttenRelRef.ASK_2014,
				AttenRelRef.BSSA_2014,
				AttenRelRef.CB_2014,
				AttenRelRef.CY_2014
		};
		Map<AttenRelRef, PlotLineType> gmmLineTypes = new HashMap<>();
		gmmLineTypes.put(AttenRelRef.ASK_2014, PlotLineType.SOLID);
		gmmLineTypes.put(AttenRelRef.BSSA_2014, PlotLineType.DASHED);
		gmmLineTypes.put(AttenRelRef.CB_2014, PlotLineType.DOTTED);
		gmmLineTypes.put(AttenRelRef.CY_2014, PlotLineType.DOTTED_AND_DASHED);
		
//		String search = "p50";
		String search = "p97.5";
//		String search = "MEAN";
		
		File[] resultsFiles = resultsDir.listFiles();
		
		for (boolean measured : new boolean[] {false,true}) {
			System.out.println("Processing for measured="+measured);
			
			Map<AttenRelRef, DiscretizedFunc> gmmAvgCurves = new HashMap<>();
			Map<AttenRelRef, Integer> gmmAvgCounts = new HashMap<>();
			Map<Double, DiscretizedFunc> vs30AvgCurves = new HashMap<>();
			Map<Double, Integer> vs30AvgCounts = new HashMap<>();
			
			CPT vs30CPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
			if (measured)
				vs30CPT = vs30CPT.rescale(500d, 800d);
			else
				vs30CPT = vs30CPT.rescale(280d, 1400d);
			
			List<DiscretizedFunc> allCurves = new ArrayList<>();
			List<PlotCurveCharacterstics> allChars = new ArrayList<>();
			
			for (File file : resultsFiles) {
				String name = file.getName();
				if (!name.startsWith(prefix) || !name.endsWith(".csv"))
					continue;
				if (measured && !name.contains("measured"))
					continue;
				if (!measured && !name.contains("inferred"))
					continue;
				AttenRelRef gmm = null;
				for (AttenRelRef ref : gmpeRefs) {
					if (name.contains(ref.name())) {
						gmm = ref;
						break;
					}
				}
				Preconditions.checkNotNull(gmm, "GMM not found for %s", name);
				int vsIndex = name.indexOf("_vs30_");
				Preconditions.checkState(vsIndex >= 0);
				String vsStr = name.substring(vsIndex+6);
				vsStr = vsStr.substring(0, vsStr.indexOf("_"));
				double vs30 = Double.valueOf(vsStr);
				System.out.println("Reading "+name+", GMM: "+gmm.name()+", vs30="+vs30);
				
				CSVFile<String> csv = CSVFile.readFile(file, true);
				int matchRow = -1;
				for (int row=csv.getNumRows(); --row>=0;) {
					if (csv.get(row, 0).equals(search)) {
						matchRow = row;
						break;
					}
				}
				Preconditions.checkState(matchRow > 0);
				DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
				for (int col=2; col<csv.getNumCols(); col++)
					curve.set(csv.getDouble(0, col), csv.getDouble(matchRow, col));
				
				allCurves.add(curve);
				allChars.add(new PlotCurveCharacterstics(gmmLineTypes.get(gmm), 1f, vs30CPT.getColor((float)vs30)));
				
				averageInCurve(gmm, curve, gmmAvgCurves, gmmAvgCounts);
				averageInCurve(vs30, curve, vs30AvgCurves, vs30AvgCounts);
			}
			
			if (allCurves.isEmpty())
				continue;
			
			normalizeInCurves(gmmAvgCurves, gmmAvgCounts);
			normalizeInCurves(vs30AvgCurves, vs30AvgCounts);
			
			String plotPrefix = prefix;
			if (measured)
				plotPrefix += "_measured";
			else
				plotPrefix += "_inferred";
			plotPrefix += "_"+search;
			
			curvePlot(allCurves, allChars, false, outputDir, plotPrefix+"_all");
			
			List<DiscretizedFunc> gmmCurves = new ArrayList<>();
			List<PlotCurveCharacterstics> gmmChars = new ArrayList<>();
			for (AttenRelRef ref : gmpeRefs) {
				DiscretizedFunc curve = gmmAvgCurves.get(ref);
				Preconditions.checkNotNull(curve);
				curve.setName(ref.getShortName());
				
				gmmCurves.add(curve);
				gmmChars.add(new PlotCurveCharacterstics(gmmLineTypes.get(ref), 3f, Color.BLACK));
			}
			
			curvePlot(gmmCurves, gmmChars, true, outputDir, plotPrefix+"_gmm_avg");

			List<DiscretizedFunc> vs30Curves = new ArrayList<>();
			List<PlotCurveCharacterstics> vs30Chars = new ArrayList<>();
			List<Double> vs30s = new ArrayList<>(vs30AvgCounts.keySet());
			Collections.sort(vs30s);
			for (double vs30 : vs30s) {
				DiscretizedFunc curve = vs30AvgCurves.get(vs30);
				curve.setName((int)vs30+" (m/s)");
				
				vs30Curves.add(curve);
				vs30Chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, vs30CPT.getColor((float)vs30)));
			}
			
			curvePlot(vs30Curves, vs30Chars, true, outputDir, plotPrefix+"_vs30_avg");
		}
	}
	
	private static <E> void normalizeInCurves(Map<E, DiscretizedFunc> curvesMap, Map<E, Integer> countsMap) {
		for (E key : curvesMap.keySet()) {
			DiscretizedFunc curve = curvesMap.get(key);
			int count = countsMap.get(key);
			curve.scale(1d/(double)count);
		}
	}
	
	private static <E> void averageInCurve(E key, DiscretizedFunc curve,
			Map<E, DiscretizedFunc> curvesMap, Map<E, Integer> countsMap) {
		if (curvesMap.containsKey(key)) {
			DiscretizedFunc avgCurve = curvesMap.get(key);
			Preconditions.checkState(avgCurve.size() == curve.size());
			for (int i=0; i<curve.size(); i++) {
				Preconditions.checkState((float)curve.getX(i) == (float)avgCurve.getX(i));
				avgCurve.set(i, avgCurve.getY(i)+curve.getY(i));
			}
			countsMap.put(key, countsMap.get(key)+1);
		} else {
			curvesMap.put(key, curve.deepClone());
			countsMap.put(key, 1);
		}
	}
	
	private static void curvePlot(List<DiscretizedFunc> curves, List<PlotCurveCharacterstics> chars, boolean legend,
			File outputDir, String prefix) throws IOException {
		PlotSpec spec = new PlotSpec(curves, chars, " ", "X", "Y");
		if (legend)
			spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, true, true, new Range(1e-2, 1e1), new Range(1e-6, 1e0));
		
		PlotUtils.writePlots(outputDir, prefix, gp, 1000, 850, true, false, false);
	}

}
