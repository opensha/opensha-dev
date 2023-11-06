package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateSegmentationConstraint.RateCombiner;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationCalculator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class SegResultPlots {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/seg_results");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File modelDir = new File(invsDir, "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		FaultSystemSolution fullSol = FaultSystemSolution.load(
				new File(modelDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		FaultSystemSolution noneSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_None.zip"));
		
		FaultSystemSolution lowSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_LowSeg.zip"));
		
		FaultSystemSolution middleSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_MidSeg.zip"));
		
		FaultSystemSolution highSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_HighSeg.zip"));
		
		FaultSystemSolution classicSol = FaultSystemSolution.load(
				new File(modelDir, "node_branch_averaged/SegModel_Classic.zip"));
		
		ClusterRuptures cRups = fullSol.getRupSet().requireModule(ClusterRuptures.class);
		PlausibilityConfiguration config = fullSol.getRupSet().requireModule(PlausibilityConfiguration.class);
		ClusterConnectionStrategy connStrat = config.getConnectionStrategy();
		
		FaultSystemSolution[] sols = {
				fullSol,
				noneSol,
				lowSol,
				middleSol,
				highSol,
				classicSol
		};
		String[] prefixes = {
				"full_ba",
				"none",
				"low",
				"middle",
				"high",
				"classic"
		};
		
		SegmentationCalculator.WRITE_PDFS = true;
		
		for (int i=0; i<sols.length; i++) {
			FaultSystemSolution sol = sols[i];
			FaultSystemRupSet rupSet = sol.getRupSet();
			// overwrite the segmentation choice so that the colors used in the dist-depend plot are uniform
			sol.requireModule(LogicTreeBranch.class).setValue(NSHM23_SegmentationModels.NONE);
			SegmentationCalculator segCalc = new SegmentationCalculator(sol, cRups.getAll(),
					connStrat, config.getDistAzCalc(), new double[] { 0d });
			String prefix = prefixes[i];
			
			segCalc.setLegendFontSize(26);
			
			segCalc.plotDistDependComparison(outputDir, "dist_dependence_"+prefix, true, null, " ");
			
			if (i == 0) {
				// write an empty plot for Ned
				SegmentationCalculator emptySegCalc = new SegmentationCalculator(sol, cRups.getAll(),
						connStrat, config.getDistAzCalc(), new double[] { 10d });
				emptySegCalc.plotDistDependComparison(outputDir, "dist_dependence_empty", true, null, " ");
			}
			
			// now annotated map
			GeographicMapMaker mapMaker = new RupSetMapMaker(rupSet,
					RupSetMapMaker.buildBufferedRegion(rupSet.getFaultSectionDataList()));
			mapMaker.setWriteGeoJSON(false);
			mapMaker.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(210, 210, 210)));
//			plotter.setJumpLineThickness(4f);
			
			CPT cpt = SegmentationCalculator.getConnectionFractCPT();
			
			mapMaker.clearJumpScalars();
			
			String label = "Passthrough Rate";
			mapMaker.plotJumpScalars(segCalc.calcJumpPassthroughs(0, RateCombiner.MIN), cpt, label);
			
			Location creepingLoc = calcCreepingLoc(rupSet);
			Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
			Stroke baseStroke = PlotLineType.SOLID.buildStroke(2f);
//			Color lineColor = new Color(0, 0, 0, 127);
			Color lineColor = Color.BLACK;
			
			double annX = -121.3;
			double annY = 35.6;
			double lineX = creepingLoc.getLongitude();
			double lineY = creepingLoc.getLatitude();
			XYTextAnnotation ann = new XYTextAnnotation("Creeping section", annX, annY);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			ann.setFont(annFont);
			mapMaker.addAnnotation(ann);
			mapMaker.addAnnotation(new XYLineAnnotation(lineX, lineY, annX, annY-0.04, baseStroke, lineColor));
			
			annX = -110.88;
			annY = 40.3;
			lineX = -111.78;
			lineY = 40.61;
			ann = new XYTextAnnotation("Wasatch", annX, annY);
			ann.setTextAnchor(TextAnchor.TOP_LEFT);
			ann.setFont(annFont);
			mapMaker.addAnnotation(ann);
			mapMaker.addAnnotation(new XYLineAnnotation(lineX, lineY, annX-0.06, annY-0.01, baseStroke, lineColor));
			
			mapMaker.plot(outputDir, "passthrough_map_"+prefix, " ", 800);
		}
	}
	
	private static Location calcCreepingLoc(FaultSystemRupSet rupSet) {
		int parentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Andreas", "Creeping");
		Preconditions.checkState(parentID >= 0);
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (sect.getParentSectionId() == parentID) {
				for (Location loc : sect.getFaultTrace()) {
					latTrack.addValue(loc.lat);
					lonTrack.addValue(loc.lon);
				}
			}
		}
		
		return new Location(0.5*(latTrack.getMin() + latTrack.getMax()), 0.5*(lonTrack.getMin() + lonTrack.getMax()));
	}

}
