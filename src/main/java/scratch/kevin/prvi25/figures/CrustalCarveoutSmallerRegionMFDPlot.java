package scratch.kevin.prvi25.figures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_MaxMagOffFault;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import net.mahdilamb.colormap.Colors;

public class CrustalCarveoutSmallerRegionMFDPlot {

	public static void main(String[] args) throws IOException {
		List<Region> regions = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		regions.add(PRVI25_SeismicityRegions.CRUSTAL.load());
		names.add("Full Model Region");
		prefixes.add("full_region");
		
		regions.add(PRVI25_RegionLoader.loadPRVI_MapExtents());
		names.add("Smaller Map Region");
		prefixes.add("smaller_map_region");
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_10_24-prvi25_crustal_branches-dmSample5x/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(2.55, sol.getRupSet().getMaxMag());
		
		NSHM23_MaxMagOffFault mMax = NSHM23_MaxMagOffFault.MAG_7p6;
		
		DecimalFormat fractDF = new DecimalFormat("0.###");
		
		for (int r=0; r<regions.size(); r++) {
			Region region = regions.get(r);
			IncrementalMagFreqDist solMFD = sol.calcNucleationMFD_forRegion(
					region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size(), false);
			IncrementalMagFreqDist leftoverMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			IncrementalMagFreqDist approxAvgLeftoverMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			
			GriddedGeoDataSet pdf = PRVI25_SeisSmoothingAlgorithms.AVERAGE.loadXYZ(
					PRVI25_SeismicityRegions.CRUSTAL, PRVI25_DeclusteringAlgorithms.AVERAGE);
			
			double fractN = 0d;
			for (int i=0; i<pdf.size(); i++)
				if (region.contains(pdf.getLocation(i)))
					fractN += pdf.get(i);
			
			boolean remapped = !region.equals(PRVI25_SeismicityRegions.CRUSTAL.load());
			
			System.out.println("fractN="+(float)fractN+" for "+names.get(r));
			
			IncrementalMagFreqDist lower = PRVI25_RegionalSeismicity.LOW.build(PRVI25_SeismicityRegions.CRUSTAL, refMFD, mMax.getMaxMagOffFault());
			IncrementalMagFreqDist upper = PRVI25_RegionalSeismicity.HIGH.build(PRVI25_SeismicityRegions.CRUSTAL, refMFD, mMax.getMaxMagOffFault());
			IncrementalMagFreqDist pref = PRVI25_RegionalSeismicity.PREFFERRED.build(PRVI25_SeismicityRegions.CRUSTAL, refMFD, mMax.getMaxMagOffFault());
			IncrementalMagFreqDist origPref = pref.deepClone();
			if (remapped) {
				lower.scale(fractN);
				upper.scale(fractN);
				pref.scale(fractN);
			}
			double lowerWeight = PRVI25_RegionalSeismicity.LOW.getNodeWeight(null);
			double upperWeight = PRVI25_RegionalSeismicity.HIGH.getNodeWeight(null);
			double prefWeight = PRVI25_RegionalSeismicity.PREFFERRED.getNodeWeight(null);
			UncertainBoundedIncrMagFreqDist dist = new UncertainBoundedIncrMagFreqDist(pref, lower, upper, UncertaintyBoundType.TWO_SIGMA);
			
			for (int i=0; i<refMFD.size(); i++) {
				leftoverMFD.set(i, Math.max(0, pref.getY(i) - solMFD.getY(i)));
				double approxAvg = lowerWeight*Math.max(0d, lower.getY(i) - solMFD.getY(i));
				approxAvg += upperWeight*Math.max(0d, upper.getY(i) - solMFD.getY(i));
				approxAvg += prefWeight*Math.max(0d, pref.getY(i) - solMFD.getY(i));
				approxAvgLeftoverMFD.set(i, approxAvg);
			}
			
			List<IncrementalMagFreqDist> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (remapped)
				dist.setName("Remapped (fractN="+fractDF.format(fractN)+") "+dist.getDefaultBoundName());
			else
				dist.setName("Full "+dist.getDefaultBoundName());
			funcs.add(dist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Colors.tab_lightpurple));
			
			if (remapped)
				pref.setName("Remapped (fractN="+fractDF.format(fractN)+") Preferred Rate");
			else
				pref.setName("Full Preferred Rate");
			funcs.add(pref);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_purple));
			
			if (remapped) {
				origPref.setName("Original Full Preferred Rate");
				funcs.add(origPref);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Colors.tab_purple));
			}
			
			solMFD.setName("Inversion Supra-Seis");
			funcs.add(solMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
			
			leftoverMFD.setName("Leftover");
			funcs.add(leftoverMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
			
			approxAvgLeftoverMFD.setName("Approx. Branch Avg. Leftover");
			funcs.add(approxAvgLeftoverMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_olive));
			
			PlotSpec spec = new PlotSpec(funcs, chars, names.get(r), "Magnitude", "Incremental Rate (1/yr)");
			spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, true, new Range(5d, 8d), new Range(1e-6, 1e1));
			
			PlotUtils.writePlots(new File("/tmp"), "mfds_"+prefixes.get(r), gp, 800, 800, true, true, false);
		}
	}

}
