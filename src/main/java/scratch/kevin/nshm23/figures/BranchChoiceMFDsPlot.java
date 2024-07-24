package scratch.kevin.nshm23.figures;

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
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

class BranchChoiceMFDsPlot {

	public static void main(String[] args) throws IOException {
		File solDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/");
		File baFile = new File(solDir, "results_WUS_FM_v3_branch_averaged_gridded.zip");
		
		File nodesDir = new File(solDir, "node_branch_averaged");
		File plotsDir = new File(solDir, "misc_plots");
		File outputDir = new File(plotsDir, "branch_choice_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution baSol = FaultSystemSolution.load(baFile);
		
		Range xRange = new Range (6d, 8.5d);
		Range incrYRange = new Range(1e-5, 2e0);
		Range cmlYRange = new Range(1e-5, 2e0);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);
		
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		List<LogicTreeNode[]> nodesList = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		nodesList.add(new NSHM23_SegmentationModels[] {
				NSHM23_SegmentationModels.CLASSIC,
				NSHM23_SegmentationModels.HIGH,
				NSHM23_SegmentationModels.MID,
				NSHM23_SegmentationModels.LOW,
				NSHM23_SegmentationModels.NONE
		});
		prefixes.add("SegModel");
		
		nodesList.add(new SupraSeisBValues[] {
				SupraSeisBValues.B_1p0,
				SupraSeisBValues.B_0p75,
				SupraSeisBValues.B_0p5,
				SupraSeisBValues.B_0p25,
				SupraSeisBValues.B_0p0
		});
		prefixes.add("SupraB");
		
		IncrementalMagFreqDist baMFD = baSol.calcNucleationMFD_forRegion(
				reg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
		
		UncertainBoundedIncrMagFreqDist dataBounds = NSHM23_RegionalSeismicity.getRemapped(reg,
				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
//		UncertainBoundedIncrMagFreqDist dataBounds = NSHM23_RegionalSeismicity.getBounded(seisReg, refMFD, refMFD.getMaxX());
		dataBounds.setName("Observed");
		dataBounds.setBoundName("95% Bounds");
		UncertainBoundedIncrMagFreqDist dataForBounds = dataBounds.deepClone();
		dataForBounds.setName(dataBounds.getBoundName());
		
		// add observed bounds
		Color obsColor = new Color(125, 80, 145); // "indigo"
		EvenlyDiscretizedFunc dataCumulative = Regional_MFD_Plots.getCmlAsFakeIncr(dataBounds);
		dataCumulative = dataCumulative.deepClone();
		dataCumulative.setName(dataBounds.getName());
		
		EvenlyDiscretizedFunc upperCumulative = Regional_MFD_Plots.getCmlAsFakeIncr(dataBounds.getUpper());
		EvenlyDiscretizedFunc lowerCumulative = Regional_MFD_Plots.getCmlAsFakeIncr(dataBounds.getLower());
		Preconditions.checkState(dataCumulative.size() == upperCumulative.size());
		for (int i=0; i<dataCumulative.size(); i++) {
			upperCumulative.set(i, Math.max(dataCumulative.getY(i), upperCumulative.getY(i)));
			lowerCumulative.set(i, Math.max(0, Math.min(dataCumulative.getY(i), lowerCumulative.getY(i))));
		}
		
		UncertainArbDiscFunc cmlBounded = new UncertainArbDiscFunc(dataCumulative, lowerCumulative, upperCumulative);
		cmlBounded.setName(dataBounds.getBoundName());
		
		for (int i=0; i<nodesList.size(); i++) {
			List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
			
			List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> cmlChars = new ArrayList<>();
			
			incrFuncs.add(dataBounds);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, obsColor));
			cmlFuncs.add(dataCumulative);
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, obsColor));
			
			incrFuncs.add(dataForBounds);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60)));
			cmlFuncs.add(cmlBounded);
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60)));
			
			LogicTreeNode[] nodes = nodesList.get(i);
			String prefix = prefixes.get(i);
			
			CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().reverse().rescale(0d, nodes.length-1);
			
			int startCopyIndex = incrFuncs.size();
			for (int j=0; j<nodes.length; j++) {
				Color color = cpt.getColor((float)j);
				
				LogicTreeNode node = nodes[j];
				
				File solFile = new File(nodesDir, prefix+"_"+node.getFilePrefix()+".zip");
				
				FaultSystemSolution sol = FaultSystemSolution.load(solFile);
				
				IncrementalMagFreqDist mfd = sol.calcNucleationMFD_forRegion(
						reg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
				
				mfd.setName(node.getShortName());
				incrFuncs.add(mfd);
				incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
				cmlFuncs.add(mfd.getCumRateDistWithOffset());
				cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
			}
			int endCopyIndex = incrFuncs.size();
			
			baMFD.setName("Average");
			incrFuncs.add(baMFD);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			cmlFuncs.add(baMFD.getCumRateDistWithOffset());
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			
			for (int j=startCopyIndex; j<endCopyIndex; j++) {
				Color color = cpt.getColor((float)(j-startCopyIndex));
				Color transColor = new Color(color.getRed(), color.getGreen(), color.getBlue(), 127);
				IncrementalMagFreqDist mfd = incrFuncs.get(j).deepClone();
				mfd.setName(null);
				incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, transColor));
				cmlFuncs.add(mfd.getCumRateDistWithOffset());
				cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, transColor));
			}
			
			PlotSpec incrSpec = new PlotSpec(incrFuncs, incrChars, " ", "Magnitude", "Incremental Rate (/yr)");
			incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			PlotSpec cmlSpec = new PlotSpec(cmlFuncs, cmlChars, " ", "Magnitude", "Cumulative Rate (/yr)");
			cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.setPlotLabelFontSize(26);
			gp.setAxisLabelFontSize(28);
			gp.setTickLabelFontSize(24);
			
			gp.drawGraphPanel(incrSpec, false, true, xRange, incrYRange);
			
			prefix += "_mfds";
			
			PlotUtils.writePlots(outputDir, prefix+"_incr", gp, 800, 750, true, true, false);
			
			gp.drawGraphPanel(cmlSpec, false, true, xRange, cmlYRange);
			
			PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, true, false);
		}
	}

}
