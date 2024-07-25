package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeBranchAverageBuilder;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class ObservedMFDComparisons {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File nodeDir = new File(invDir, "node_branch_averaged");
		File meanSolFile = new File(invDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File ltFile = new File(invDir, "logic_tree.json");
		
		File outputDir = new File(invDir, "observed_mfd_comparison");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		AnalysisRegions[] aRegs = {
				AnalysisRegions.CONUS_U3_RELM,
				AnalysisRegions.CONUS_IMW,
				AnalysisRegions.CONUS_PNW
		};
		
		String pathPrefix = "/data/erf/nshm23/seismicity/observed";
		String[] csvNames = {
				"ucerf_observed.csv",
				"imw_observed.csv",
				"pnw_observed.csv"
		};
		
		Region wusReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		IncrementalMagFreqDist[] incrMFDs = new IncrementalMagFreqDist[aRegs.length];
		EvenlyDiscretizedFunc[] cmlMFDs = new EvenlyDiscretizedFunc[incrMFDs.length];
		
		Range magRage = new Range(5d, 8.5d);
		Range yRange = new Range(1e-6, 1e1);
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.45);
		
		for (int i=0; i<incrMFDs.length; i++) {
			System.out.println("Loading "+csvNames[i]);
			CSVFile<String> csv = CSVFile.readStream(
					NSHM23_RegionalSeismicity.class.getResourceAsStream(pathPrefix+"/"+csvNames[i]), true);
			incrMFDs[i] = loadIncrMFD(csv, refMFD);
			cmlMFDs[i] = loadCmlMFD(csv, incrMFDs[i].getCumRateDistWithOffset());
			System.out.println(cmlMFDs[i]);
		}
		
		IncrementalMagFreqDist totIncrMFD = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (IncrementalMagFreqDist mfd : incrMFDs)
			for (int i=0; i<mfd.size(); i++)
				totIncrMFD.set(i, totIncrMFD.getY(i)+mfd.getY(i));
		EvenlyDiscretizedFunc totCmlMFD = totIncrMFD.getCumRateDistWithOffset();
		System.out.println("Total:\n"+totCmlMFD);
		
		UncertainBoundedIncrMagFreqDist[] estIncrMFDs = new UncertainBoundedIncrMagFreqDist[aRegs.length];
		UncertainArbDiscFunc[] estCmlMFDs = new UncertainArbDiscFunc[aRegs.length];
		
		for (int i=0; i<incrMFDs.length; i++) {
			estIncrMFDs[i] = NSHM23_RegionalSeismicity.getRemapped(aRegs[i].load(),
					NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
			EvenlyDiscretizedFunc cmlMean = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFDs[i]);
			EvenlyDiscretizedFunc cmlUpper = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFDs[i].getUpper());
			EvenlyDiscretizedFunc cmlLower = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFDs[i].getLower());
			estCmlMFDs[i] = new UncertainArbDiscFunc(cmlMean, cmlLower, cmlUpper, estIncrMFDs[i].getBoundType());
		}
		
		UncertainBoundedIncrMagFreqDist estTotIncrMFD = NSHM23_RegionalSeismicity.getRemapped(wusReg,
				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
		EvenlyDiscretizedFunc estCmlMean = Regional_MFD_Plots.getCmlAsFakeIncr(estTotIncrMFD);
		EvenlyDiscretizedFunc estCmlUpper = Regional_MFD_Plots.getCmlAsFakeIncr(estTotIncrMFD.getUpper());
		EvenlyDiscretizedFunc estCmlLower = Regional_MFD_Plots.getCmlAsFakeIncr(estTotIncrMFD.getLower());
		UncertainArbDiscFunc estTotCmlMFD = new UncertainArbDiscFunc(estCmlMean, estCmlLower, estCmlUpper, estTotIncrMFD.getBoundType());
		
		FaultSystemSolution baSol = FaultSystemSolution.load(meanSolFile);
		LogicTree<?> tree = LogicTree.read(ltFile);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Observed MFD Comparisons");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		Color estColor = new Color(125, 80, 145); // "indigo"
		Color estTransColor = new Color(estColor.getRed(), estColor.getGreen(), estColor.getBlue(), 60);
		
		List<LogicTreeLevel<?>> levels = new ArrayList<>();
		levels.add(null);
		for (LogicTreeLevel<?> level : tree.getLevels())
			levels.add(level);
		
		for (int l=0; l<levels.size(); l++) {
			LogicTreeLevel<?> level = levels.get(l);
			
			List<FaultSystemSolution> nodeSols = new ArrayList<>();
			List<LogicTreeNode> nodes = new ArrayList<>();
			String prefix;
			String label;
			if (level != null) {
				prefix = MPJ_LogicTreeBranchAverageBuilder.levelPrefix(level);
				for (LogicTreeNode node : level.getNodes()) {
					File solFile = new File(nodeDir, prefix+"_"+node.getFilePrefix()+".zip");
					if (solFile.exists()) {
						nodeSols.add(FaultSystemSolution.load(solFile));
						nodes.add(node);
					}
				}
				if (nodes.isEmpty())
					continue;
				
				System.out.println("Processing for "+level.getName()+" with "+nodes.size()+" node solutions");
				label = level.getName();
			} else {
				System.out.println("Processing for mean solution");
				label = "Branch-Averaged Solution";
				prefix = "avg_sol";
			}
			
			lines.add("## "+label);
			lines.add(topLink); lines.add("");
			

			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			table.addLine("Incremental", "Cumulative");
			
			for (int r=-1; r<aRegs.length; r++) {
				Region reg;
				IncrementalMagFreqDist obsIncr;
				EvenlyDiscretizedFunc obsCml;
				UncertainBoundedIncrMagFreqDist estIncr;
				UncertainArbDiscFunc estCml;
				String regPrefix;
				
				if (r < 0) {
					reg = wusReg;
					regPrefix = prefix+"_WUS";
					obsIncr = totIncrMFD;
					obsCml = totCmlMFD;
					estIncr = estTotIncrMFD;
					estCml = estTotCmlMFD;
				} else {
					reg = aRegs[r].load();
					regPrefix = prefix+"_"+aRegs[r].name();
					obsIncr = incrMFDs[r];
					obsCml = cmlMFDs[r];
					estIncr = estIncrMFDs[r];
					estCml = estCmlMFDs[r];
				}
				
				String regName = reg.getName();
				
				List<DiscretizedFunc> incrFuncs = new ArrayList<>();
				List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				estIncr.setName("Estimated");
				incrFuncs.add(estIncr);
				estCml.setName(estIncr.getName());
				cmlFuncs.add(estCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, estColor));
				
				UncertainBoundedIncrMagFreqDist estForBounds = estIncr.deepClone();
				estForBounds.setName(estForBounds.getBoundName());
				incrFuncs.add(estForBounds);
				UncertainArbDiscFunc estCmlForBounds = estCml.deepClone();
				estCmlForBounds.setName(estCmlForBounds.getBoundName());
				cmlFuncs.add(estCmlForBounds);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 3f, estTransColor));
				
				obsIncr.setName("Observed");
				incrFuncs.add(obsIncr);
				obsCml.setName(obsIncr.getName());
				cmlFuncs.add(obsCml);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.BLACK));
				
				// ba curve
				IncrementalMagFreqDist baIncr = baSol.calcNucleationMFD_forRegion(
						reg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
				baIncr.setName("Branch-Averaged Supra-Seis");
				incrFuncs.add(baIncr);
				EvenlyDiscretizedFunc baCml = baIncr.getCumRateDistWithOffset();
				baCml.setName(baIncr.getName());
				cmlFuncs.add(baCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.DARK_GRAY));
				
				GridSourceProvider gridProv = baSol.getGridSourceProvider();
				if (gridProv != null) {
					IncrementalMagFreqDist baTotIncr = baIncr.deepClone();
					for (int i=0; i<gridProv.getNumLocations(); i++) {
						Location loc = gridProv.getGriddedRegion().getLocation(i);
						if (reg.contains(loc)) {
							IncrementalMagFreqDist gridMFD = gridProv.getMFD(i);
							if (gridMFD != null) {
								for (Point2D pt : gridMFD)
									baTotIncr.add(baTotIncr.getClosestXIndex(pt.getX()), pt.getY());
							}
						}
					}
					
					baTotIncr.setName("Branch-Averaged Total");
					incrFuncs.add(baTotIncr);
					EvenlyDiscretizedFunc baTotCml = baTotIncr.getCumRateDistWithOffset();
					baTotCml.setName(baTotIncr.getName());
					cmlFuncs.add(baTotCml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.DARK_GRAY));
				}
				
				if (nodes.size() > 1) {
					CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, nodes.size()-1d);
					
					for (int n=0; n<nodes.size(); n++) {
						FaultSystemSolution nodeSol = nodeSols.get(n);
						IncrementalMagFreqDist nodeIncr = nodeSol.calcNucleationMFD_forRegion(
								reg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
						LogicTreeNode node = nodes.get(n);
						nodeIncr.setName(node.getShortName());
						incrFuncs.add(nodeIncr);
						EvenlyDiscretizedFunc nodeCml = nodeIncr.getCumRateDistWithOffset();
						nodeCml.setName(nodeIncr.getName());
						cmlFuncs.add(nodeCml);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt.getColor((float)n)));
					}
				}
				
				PlotSpec incrSpec = new PlotSpec(incrFuncs, chars,
						reg.getName(), "Magnitude", "Incremental Rate (/yr)");
				incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				PlotSpec cmlSpec = new PlotSpec(cmlFuncs, chars,
						reg.getName(), "Magnitude", "Cumulative Rate (/yr)");
				cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.getPlotPrefs().scaleFontSizes(1.2d);

				gp.drawGraphPanel(incrSpec, false, true, magRage, yRange);
				PlotUtils.writePlots(resourcesDir, regPrefix+"_incr", gp, 950, 800, true, true, true);
				
				gp.drawGraphPanel(cmlSpec, false, true, magRage, yRange);
				PlotUtils.writePlots(resourcesDir, regPrefix+"_cml", gp, 950, 800, true, true, true);
				table.addLine("![MFD]("+resourcesDir.getName()+"/"+regPrefix+"_incr.png)",
						"![MFD]("+resourcesDir.getName()+"/"+regPrefix+"_cml.png)");
			}
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final double OBS_YEARS = 2023d-1930d;
	private static final double OBS_SCALAR = 1d/OBS_YEARS;
	
	private static IncrementalMagFreqDist loadIncrMFD(CSVFile<String> csv, EvenlyDiscretizedFunc refMFD) {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		for (int row=1; row<csv.getNumRows(); row++) {
			double mag = csv.getDouble(row, 0);
			double rate = csv.getDouble(row, 3);
			mfd.set(mfd.getClosestXIndex(mag+0.01), rate*OBS_SCALAR);
		}
		
		return mfd;
	}
	
	private static IncrementalMagFreqDist loadCmlMFD(CSVFile<String> csv, EvenlyDiscretizedFunc refMFD) {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		for (int row=1; row<csv.getNumRows(); row++) {
			double mag = csv.getDouble(row, 0);
			double rate = csv.getDouble(row, 4);
			mfd.set(mfd.getClosestXIndex(mag+0.01), rate*OBS_SCALAR);
		}
		
		int minXNonZero = -1;
		for (int i=mfd.size(); --i>=0;) {
			if (mfd.getY(i) > 0) {
				minXNonZero = i;
			} else {
				if (minXNonZero >= 0)
					// we've encountered nonzero, and now back to zero
					break;
			}
		}
		for (int i=0; i<minXNonZero; i++)
			mfd.set(i, mfd.getY(minXNonZero));
		
		return mfd;
	}

}
