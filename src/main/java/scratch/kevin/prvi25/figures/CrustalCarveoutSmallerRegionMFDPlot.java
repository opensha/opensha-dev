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
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import net.mahdilamb.colormap.Colors;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class CrustalCarveoutSmallerRegionMFDPlot {

	public static void main(String[] args) throws IOException {
		List<Region> regions = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		String fullName = "Full Model Region";
		regions.add(PRVI25_SeismicityRegions.CRUSTAL.load());
		names.add(fullName);
		prefixes.add("full_region");
		
		regions.add(PRVI25_RegionLoader.loadPRVI_MapExtents());
		names.add("Hazard Assessment Region");
		prefixes.add("smaller_map_region");
		
		FaultSystemSolution sol = FaultSystemSolution.load(CRUSTAL_SOL_SUPRA_ONLY);
		
		boolean includeCarveout = false;
		
		File outputDir = new File("/tmp");
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(2.55, sol.getRupSet().getMaxMag());
		
		NSHM23_MaxMagOffFault mMax = NSHM23_MaxMagOffFault.MAG_7p9;
		PRVI25_SeismicityRateEpoch epoch = PRVI25_SeismicityRateEpoch.RECENT_SCALED;
		
		DecimalFormat fractDF = new DecimalFormat("0.###");
		
		for (int r=0; r<regions.size(); r++) {
			Region region = regions.get(r);
			IncrementalMagFreqDist solMFD = sol.calcNucleationMFD_forRegion(
					region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size(), false);
			
			GriddedGeoDataSet pdf = PRVI25_SeisSmoothingAlgorithms.AVERAGE.loadXYZ(
					PRVI25_SeismicityRegions.CRUSTAL, PRVI25_DeclusteringAlgorithms.AVERAGE);
			
			double fractN = 0d;
			for (int i=0; i<pdf.size(); i++)
				if (region.contains(pdf.getLocation(i)))
					fractN += pdf.get(i);
			
			boolean remapped = !region.equals(PRVI25_SeismicityRegions.CRUSTAL.load());
			
			System.out.println("fractN="+(float)fractN+" for "+names.get(r));
			
			IncrementalMagFreqDist lower = PRVI25_CrustalSeismicityRate.LOW.build(epoch, refMFD, mMax.getMaxMagOffFault());
			IncrementalMagFreqDist upper = PRVI25_CrustalSeismicityRate.HIGH.build(epoch, refMFD, mMax.getMaxMagOffFault());
			IncrementalMagFreqDist pref = PRVI25_CrustalSeismicityRate.PREFFERRED.build(epoch, refMFD, mMax.getMaxMagOffFault());
			IncrementalMagFreqDist origPref = pref.deepClone();
			if (remapped) {
				lower.scale(fractN);
				upper.scale(fractN);
				pref.scale(fractN);
			}
			double lowerWeight = PRVI25_CrustalSeismicityRate.LOW.getNodeWeight(null);
			double upperWeight = PRVI25_CrustalSeismicityRate.HIGH.getNodeWeight(null);
			double prefWeight = PRVI25_CrustalSeismicityRate.PREFFERRED.getNodeWeight(null);
			UncertainBoundedIncrMagFreqDist dist = new UncertainBoundedIncrMagFreqDist(pref, lower, upper, UncertaintyBoundType.TWO_SIGMA);
			
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
				origPref.setName(fullName+" Preferred Rate");
				funcs.add(origPref);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Colors.tab_purple));
			}
			
			solMFD.setName("Inversion Supra-Seis");
			funcs.add(solMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
			
			if (includeCarveout) {
				IncrementalMagFreqDist leftoverMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				IncrementalMagFreqDist approxAvgLeftoverMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				
				for (int i=0; i<refMFD.size(); i++) {
					leftoverMFD.set(i, Math.max(0, pref.getY(i) - solMFD.getY(i)));
					double approxAvg = lowerWeight*Math.max(0d, lower.getY(i) - solMFD.getY(i));
					approxAvg += upperWeight*Math.max(0d, upper.getY(i) - solMFD.getY(i));
					approxAvg += prefWeight*Math.max(0d, pref.getY(i) - solMFD.getY(i));
					approxAvgLeftoverMFD.set(i, approxAvg);
				}
				
				leftoverMFD.setName("Leftover");
				funcs.add(leftoverMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
				
				approxAvgLeftoverMFD.setName("Approx. Branch Avg. Leftover");
				funcs.add(approxAvgLeftoverMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_olive));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, names.get(r), "Magnitude", "Incremental Rate (1/yr)");
			spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, true, new Range(5d, 8d), new Range(1e-6, 1e1));
			
			PlotUtils.writePlots(outputDir, "mfds_"+prefixes.get(r), gp, 800, 800, true, true, false);
		}
	}

}
