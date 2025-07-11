package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class BulgeContributionMFDPlot {

	public static void main(String[] args) throws IOException {
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2024_12_12-prvi25_crustal_branches-dmSample5x/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		FaultSystemSolution sol = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED);
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/markdown/inversions/"
//				+ "2025_01_16-prvi25-GEOLOGIC_DIST_AVG_NSHM23_Avg_SupraB0.5_AvgSeg/solution.zip"));
		File outputDir = new File("/tmp");
		
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(2.55, 8.05);
		
		IncrementalMagFreqDist fullMFD = sol.calcNucleationMFD_forRegion(null, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
		
		PRVI25_SeismicityRateEpoch epoch = PRVI25_SeismicityRateEpoch.DEFAULT;
		UncertainBoundedIncrMagFreqDist obsBounds = PRVI25_CrustalSeismicityRate.loadRateModel(epoch).getBounded(refMFD, 8.05);
		// find the largest bulge;
		double maxRatio = 0d;
		double sortMag = Double.NaN;
		int sortMagIndex = -1;
		for (int i=0; i<fullMFD.size(); i++) {
			double mag = fullMFD.getX(i);
			Preconditions.checkState((float)mag == (float)obsBounds.getX(i));
			double obsRate = obsBounds.getY(i);
			double invRate = fullMFD.getY(i);
			if (invRate > obsRate) {
				double ratio = invRate / obsRate;
				if (ratio > maxRatio) {
					maxRatio = ratio;
					sortMag = mag;
					sortMagIndex = i;
				}
			}
		}
		Preconditions.checkState(sortMagIndex >= 0, "No bulge");
		System.out.println("Largest bulge is for M"+(float)sortMag+"; ratio="+(float)maxRatio);
		
		CPT colors = GMT_CPT_Files.CATEGORICAL_TAB10_NOGRAY.instance();
		
		for (boolean remapped : new boolean[] {false,true}) {
			Map<Integer, IncrementalMagFreqDist> parentParticMFDs = new HashMap<>();
			Map<Integer, Double> parentParticSortRates = new HashMap<>();
			
			if (remapped) {
				HashSet<Integer> processedParents = new HashSet<>();
				Map<String, String> remappings = FaultTableWriter.getNameRemappings();
				Map<String, HashSet<Integer>> remappedRupIDs = new HashMap<>();
				for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
					int parentID = sect.getParentSectionId();
					if (!processedParents.contains(parentID)) {
						processedParents.add(parentID);
						String name = sect.getParentSectionName();
						for (String prefix : remappings.keySet()) {
							if (name.startsWith(prefix)) {
								name = remappings.get(prefix);
								break;
							}
						}
						List<Integer> rups = sol.getRupSet().getRupturesForParentSection(parentID);
						if (!remappedRupIDs.containsKey(name))
							remappedRupIDs.put(name, new HashSet<>());
						remappedRupIDs.get(name).addAll(rups);
					}
				}
				int fakeID = 0;
				RupMFDsModule rupMFDs = sol.getModule(RupMFDsModule.class);
				for (String name : remappedRupIDs.keySet()) {
					HashSet<Integer> rups = remappedRupIDs.get(name);
					
					IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
					for (int rupIndex : rups) {
						DiscretizedFunc rupMFD = rupMFDs != null ? rupMFDs.getRuptureMFD(rupIndex) : null;
						if (rupMFD == null) {
							mfd.add(mfd.getClosestXIndex(sol.getRupSet().getMagForRup(rupIndex)), sol.getRateForRup(rupIndex));
						} else {
							for (Point2D pt : rupMFD)
								mfd.add(mfd.getClosestXIndex(pt.getX()), pt.getY());
						}
					}
					name = name.replace("\\", "");
					mfd.setName(name);
					parentParticMFDs.put(fakeID, mfd);
					parentParticSortRates.put(fakeID, mfd.getY(sortMagIndex));
					fakeID++;
				}
			} else {
				for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
					if (!parentParticMFDs.containsKey(sect.getParentSectionId())) {
						IncrementalMagFreqDist parentMFD = sol.calcParticipationMFD_forParentSect(
								sect.getParentSectionId(), refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
						parentMFD.setName(sect.getParentSectionName());
						parentParticMFDs.put(sect.getParentSectionId(), parentMFD);
						parentParticSortRates.put(sect.getParentSectionId(), parentMFD.getY(sortMagIndex));
					}
				}
			}
			
			List<Integer> sortedParents = ComparablePairing.getSortedData(parentParticSortRates);
			Collections.reverse(sortedParents);
			
			sortedParents = sortedParents.subList(0, colors.size());
			
			List<IncrementalMagFreqDist> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			obsBounds.setName("Observed (& 95% bounds)");
			funcs.add(obsBounds);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
			
			fullMFD.setName("Total On-Fault");
			funcs.add(fullMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			IncrementalMagFreqDist obsLower = obsBounds.getLower();
			obsLower.setName(null);
			funcs.add(obsLower);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
			
			IncrementalMagFreqDist obsUpper = obsBounds.getUpper();
			obsUpper.setName(null);
			funcs.add(obsUpper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
			
			System.out.println("Sorted faults (remapped="+remapped+")");
			for (int i=0; i<sortedParents.size(); i++) {
				IncrementalMagFreqDist mfd = parentParticMFDs.get(sortedParents.get(i));
				System.out.println("\t"+mfd.getName());
				funcs.add(mfd);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(i).minColor));
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
//			plot.setLegendVisible(true);
//			plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			plot.setLegendInset(RectangleAnchor.TOP, 0.5, 0.98, 0.9, false);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.getPlotPrefs().setLegendFontSize(18);
			
			gp.drawGraphPanel(plot, false, true, new Range(5d, 8d), new Range(1e-6, 1e0));
			
			String prefix = "mfd_bulge_contributions";
			if (remapped)
				prefix += "_remapped";
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
		}
	}

}
