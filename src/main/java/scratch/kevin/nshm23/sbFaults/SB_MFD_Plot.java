package scratch.kevin.nshm23.sbFaults;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

public class SB_MFD_Plot {
	
	public static void main(String[] args) throws IOException {
		File sbDir = new File("/home/kevin/OpenSHA/2025_sb_faults");
		Region region = Region.fromFeature(Feature.read(new File(sbDir, "region.geojson")));
		
		File mfdOutputDir = new File(sbDir, "regional_mfds");
		Preconditions.checkState(mfdOutputDir.exists() || mfdOutputDir.mkdir());
		
		File origDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3");
		File nodeBADir = new File(origDir, "node_branch_averaged");
		
		FaultSystemSolution baSol = FaultSystemSolution.load(new File(origDir, "results_WUS_FM_v3_branch_averaged_mod_gridded.zip"));
		PlotCurveCharacterstics baChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK);
		
		List<List<File>> nodeSolFiles = new ArrayList<>();
		List<List<String>> nodeSolLabels = new ArrayList<>();
		List<String> nodePrefixes = new ArrayList<>();
		List<String> nodeTitles = new ArrayList<>();
		
		nodeSolFiles.add(null);
		nodeSolLabels.add(null);
		nodePrefixes.add("mfds");
		nodeTitles.add(" ");
		
		List<File> segSolFiles = new ArrayList<>();
		List<String> segLabels = new ArrayList<>();
		nodeSolFiles.add(segSolFiles);
		nodeSolLabels.add(segLabels);
		nodePrefixes.add("mfds_seg");
		nodeTitles.add("Segmentation model");
		
		segSolFiles.add(new File(nodeBADir, "SegModel_Classic.zip"));
		segLabels.add("Classic");
		
		segSolFiles.add(new File(nodeBADir, "SegModel_High.zip"));
		segLabels.add("High");
		
		segSolFiles.add(new File(nodeBADir, "SegModel_Middle.zip"));
		segLabels.add("Middle");
		
		segSolFiles.add(new File(nodeBADir, "SegModel_Low.zip"));
		segLabels.add("Low");
		
		segSolFiles.add(new File(nodeBADir, "SegModel_None.zip"));
		segLabels.add("None");
		
		List<File> dmSolFiles = new ArrayList<>();
		List<String> dmLabels = new ArrayList<>();
		nodeSolFiles.add(dmSolFiles);
		nodeSolLabels.add(dmLabels);
		nodePrefixes.add("mfds_dm");
		nodeTitles.add("Deformation model");
		
		dmSolFiles.add(new File(nodeBADir, "DM_EVANS.zip"));
		dmLabels.add("Evans");
		
		dmSolFiles.add(new File(nodeBADir, "DM_GEOLOGIC.zip"));
		dmLabels.add("Geologic");
		
		dmSolFiles.add(new File(nodeBADir, "DM_POLLITZ.zip"));
		dmLabels.add("Pollitz");
		
		dmSolFiles.add(new File(nodeBADir, "DM_SHEN_BIRD.zip"));
		dmLabels.add("Shen-Bird");
		
		dmSolFiles.add(new File(nodeBADir, "DM_ZENG.zip"));
		dmLabels.add("Zeng");
		
		List<File> bSolFiles = new ArrayList<>();
		List<String> bLabels = new ArrayList<>();
		nodeSolFiles.add(bSolFiles);
		nodeSolLabels.add(bLabels);
		nodePrefixes.add("mfds_b_value");
		nodeTitles.add("On-fault b-value");
		
		bSolFiles.add(new File(nodeBADir, "SupraB_SupraB1.0.zip"));
		bLabels.add("b=1");
		
		bSolFiles.add(new File(nodeBADir, "SupraB_SupraB0.75.zip"));
		bLabels.add("b=0.75");
		
		bSolFiles.add(new File(nodeBADir, "SupraB_SupraB0.5.zip"));
		bLabels.add("b=0.5");
		
		bSolFiles.add(new File(nodeBADir, "SupraB_SupraB0.25.zip"));
		bLabels.add("b=0.25");
		
		bSolFiles.add(new File(nodeBADir, "SupraB_SupraB0.0.zip"));
		bLabels.add("b=0");
		
		List<File> scaleFiles = new ArrayList<>();
		List<String> scaleLabels = new ArrayList<>();
		nodeSolFiles.add(scaleFiles);
		nodeSolLabels.add(scaleLabels);
		nodePrefixes.add("mfds_scale");
		nodeTitles.add("Scaling relationship");
		
		for (NSHM23_ScalingRelationships scale : NSHM23_ScalingRelationships.values()) {
			if (scale.getNodeWeight(null) > 0) {
				scaleFiles.add(new File(nodeBADir, "Scale_"+scale.getFilePrefix()+".zip"));
				scaleLabels.add(scale.getShortName());
			}
		}
		
		BranchSectNuclMFDs baNuclMFDs = baSol.requireModule(BranchSectNuclMFDs.class);
		
		double[] regionalNuclFracts = baSol.getRupSet().getFractSectsInsideRegion(region, false);
		
		EvenlyDiscretizedFunc[] cmlMFDFractiles = baNuclMFDs.calcCumulativeFractiles(regionalNuclFracts, SolMFDPlot.standardFractiles);
		
		Range magRange = new Range(6d, 8.25d);
		Range rateRange = new Range(1e-5, 1e0);
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(magRange.getLowerBound()+0.01, magRange.getUpperBound()-0.01);
		
//		Color refColor = Colors.tab_red;
//		Color transColor = new Color(refColor.getRed(), refColor.getGreen(), refColor.getBlue(), 40);
		Color transColor = new Color(0, 0, 0, 40);
		PlotCurveCharacterstics minMaxChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, transColor);
		
		IncrementalMagFreqDist baMFD = baSol.calcNucleationMFD_forRegion(region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size(), false);
		baMFD.setName("Branch-averaged");
		EvenlyDiscretizedFunc baCmlMFD = baMFD.getCumRateDistWithOffset();

		boolean doGridded = false;
		EvenlyDiscretizedFunc baCmlWithGridded = null;
		if (doGridded) {
			SummedMagFreqDist baWithGridded = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			baWithGridded.addIncrementalMagFreqDist(baMFD);
			GridSourceProvider gridProv = baSol.getGridSourceProvider();
			for (int l=0; l<gridProv.getNumLocations(); l++)
				if (region.contains(gridProv.getLocation(l)))
					baWithGridded.addIncrementalMagFreqDist(gridProv.getMFD(l));
			baWithGridded.setName("Branch-averaged (fault + gridded)");
			baCmlWithGridded = baWithGridded.getCumRateDistWithOffset();
		}
		
		for (int n=0; n<nodeSolFiles.size(); n++) {
			List<File> baFiles = nodeSolFiles.get(n);
			List<String> baLabels = nodeSolLabels.get(n);
			String prefix = nodePrefixes.get(n);
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();

			funcs.add(baCmlMFD);
			chars.add(baChar);
			
			if (doGridded) {
				funcs.add(baCmlWithGridded);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, baChar.getLineWidth(), baChar.getColor()));
			}
			
			for (UncertainArbDiscFunc func : SolMFDPlot.processCmlFractiles(cmlMFDFractiles, magRange.getLowerBound())) {
				if (funcs.size() == 1)
					func.setName(SolMFDPlot.fractileLabel);
				else
					func.setName(null);
				funcs.add(func);
				chars.add(minMaxChar);
			}
			
			if (baFiles != null) {
				CPT rainbow = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, baFiles.size()-1).reverse();
				
				for (int s=0; s<baFiles.size(); s++) {
					FaultSystemSolution sol = FaultSystemSolution.load(baFiles.get(s));
					IncrementalMagFreqDist mfd = sol.calcNucleationMFD_forRegion(region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size(), false);
					mfd.setName(baLabels.get(s));
					EvenlyDiscretizedFunc cmlMFD = mfd.getCumRateDistWithOffset();
					
					funcs.add(cmlMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, rainbow.getColor(s)));
				}
				
				// again on top
				baCmlMFD = baCmlMFD.deepClone();
				baCmlMFD.setName(null);
				funcs.add(baCmlMFD);
				chars.add(baChar);
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, nodeTitles.get(n), "Magnitude", "Cumulative Rate (1/yr)");
			plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plot, false, true, magRange, rateRange);
			
			PlotUtils.writePlots(mfdOutputDir, prefix, gp, 1000, 850, true, true, false);
		}
	}

}
