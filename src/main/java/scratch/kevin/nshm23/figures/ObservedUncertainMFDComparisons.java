package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
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
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class ObservedUncertainMFDComparisons {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/data/kevin/nshm23/batch_inversions/"
//				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File nodeDir = new File(invDir, "node_branch_averaged");
		File meanSolFile = new File(invDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File ltFile = new File(invDir, "logic_tree.json");
		
		File outputDir = new File(invDir, "observed_mfd_uncertainty_comparison");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		String pathPrefix = "/data/erf/nshm23/seismicity/observed";
		
		Region wusReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		Range magRage = new Range(6d, 8.5d);
		Range yRange = new Range(1e-6, 2e0);
		EvenlyDiscretizedFunc refIncrMFD = FaultSysTools.initEmptyMFD(8.45);
		
		int endYear = 2023;
		int[] startYears = { 1850, 1900, 1930 };
		
		int maxZeroBins = 2;
		
		ArbitrarilyDiscretizedFunc[] obsMFDs = new ArbitrarilyDiscretizedFunc[startYears.length];
		UncertainArbDiscFunc[] bounds95 = new UncertainArbDiscFunc[startYears.length];
		UncertainArbDiscFunc[] bounds68 = new UncertainArbDiscFunc[startYears.length];
		
		for (int i=0; i<startYears.length; i++) {
			String path = pathPrefix+"/WUS-Other-"+startYears[i]+"-2023-direct-rates-v2.txt";
			
			InputStream stream = NSHM23_RegionLoader.class.getResourceAsStream(path);
			BufferedReader read = new BufferedReader(new InputStreamReader(stream));
			
			ArbitrarilyDiscretizedFunc p2p5 = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc p16 = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc obs = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc p84 = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc p97p5 = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc middle = new ArbitrarilyDiscretizedFunc();
			
			int numZeroBins = 0;
			String line = read.readLine();
			while (line != null) {
				line = line.trim();
				String[] split = line.split(" ");
				if (split.length == 8) {
					// might be a data line
					try {
						int index = 0;
						double mag = Double.parseDouble(split[index++]);
						int count = Integer.parseInt(split[index++]);
						index++; // duration
						double rate = Double.parseDouble(split[index++]);
						double lower2p5 = Double.parseDouble(split[index++]);
						double lower16 = Double.parseDouble(split[index++]);
						double upper84 = Double.parseDouble(split[index++]);
						double upper97p5 = Double.parseDouble(split[index++]);
						
						if (count == 0 && maxZeroBins == 0)
							break;
						
						obs.set(mag, count>0 ? rate : 0d);
						p2p5.set(mag, lower2p5);
						p16.set(mag, lower16);
						p84.set(mag, upper84);
						p97p5.set(mag, upper97p5);
						middle.set(mag, 0.5*(lower16+upper84));
						
						if (count == 0) {
							numZeroBins++;
							if (numZeroBins == maxZeroBins)
								break;
						}
					} catch (NumberFormatException e) {}
				} 
				
				line = read.readLine();
			}
			
			read.close();
			
			obsMFDs[i] = obs;
			bounds95[i] = new UncertainArbDiscFunc(middle, p2p5, p97p5, UncertaintyBoundType.CONF_95);
			bounds68[i] = new UncertainArbDiscFunc(middle, p16, p84, UncertaintyBoundType.CONF_68);
		}
		
		UncertainBoundedIncrMagFreqDist estIncrMFD = NSHM23_RegionalSeismicity.getRemapped(wusReg,
				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refIncrMFD, refIncrMFD.getMaxX());
		EvenlyDiscretizedFunc cmlMean = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFD);
		EvenlyDiscretizedFunc cmlUpper = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFD.getUpper());
		EvenlyDiscretizedFunc cmlLower = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFD.getLower());
		UncertainArbDiscFunc estCmlMFD = new UncertainArbDiscFunc(cmlMean, cmlLower, cmlUpper, estIncrMFD.getBoundType());
		
		FaultSystemSolution baSol = FaultSystemSolution.load(meanSolFile);
		LogicTree<?> tree = LogicTree.read(ltFile);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Observed MFD Comparisons");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		PlotCurveCharacterstics estPrefChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY);
		PlotCurveCharacterstics estBoundsChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY);
		
		Color obsColor = new Color(125, 80, 145); // "indigo"
		Color obsTransColor = new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 60);
		
		List<LogicTreeLevel<?>> levels = new ArrayList<>();
		levels.add(null);
		for (LogicTreeLevel<?> level : tree.getLevels())
			levels.add(level);
		
		// ba curve
		IncrementalMagFreqDist baIncr = baSol.calcNucleationMFD_forRegion(
				wusReg, refIncrMFD.getMinX(), refIncrMFD.getMaxX(), refIncrMFD.getDelta(), false);
		EvenlyDiscretizedFunc baCml = baIncr.getCumRateDistWithOffset();
		
		GridSourceProvider gridProv = baSol.getGridSourceProvider();
		EvenlyDiscretizedFunc baTotCml = null;
		if (gridProv != null) {
			IncrementalMagFreqDist baTotIncr = baIncr.deepClone();
			for (int i=0; i<gridProv.size(); i++) {
				Location loc = gridProv.getGriddedRegion().getLocation(i);
				if (wusReg.contains(loc)) {
					IncrementalMagFreqDist gridMFD = gridProv.getMFD(i);
					if (gridMFD != null) {
						for (Point2D pt : gridMFD)
							baTotIncr.add(baTotIncr.getClosestXIndex(pt.getX()), pt.getY());
					}
				}
			}
			baTotCml = baTotIncr.getCumRateDistWithOffset();
		}
		
		
		for (int l=0; l<levels.size(); l++) {
			LogicTreeLevel<?> level = levels.get(l);
			
			List<FaultSystemSolution> nodeSols = new ArrayList<>();
			List<LogicTreeNode> nodes = new ArrayList<>();
			List<EvenlyDiscretizedFunc> nodeFuncs = new ArrayList<>();
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
				if (nodes.size() < 2)
					continue;
				
				System.out.println("Processing for "+level.getName()+" with "+nodes.size()+" node solutions");
				label = level.getName();
				
				for (int n=0; n<nodes.size(); n++) {
					FaultSystemSolution nodeSol = nodeSols.get(n);
					IncrementalMagFreqDist nodeIncr = nodeSol.calcNucleationMFD_forRegion(
							wusReg, refIncrMFD.getMinX(), refIncrMFD.getMaxX(), refIncrMFD.getDelta(), false);
					LogicTreeNode node = nodes.get(n);
					EvenlyDiscretizedFunc nodeCml = nodeIncr.getCumRateDistWithOffset();
					nodeCml.setName(node.getShortName());
					nodeFuncs.add(nodeCml);
				}
			} else {
				System.out.println("Processing for mean solution");
				label = "Branch-Averaged Solution";
				prefix = "avg_sol";
			}
			
			lines.add("## "+label);
			lines.add(topLink); lines.add("");
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			table.initNewLine();
			for (int startYear : startYears)
				table.addColumn(MarkdownUtils.boldCentered("Years: "+startYear+"-"+endYear));
			table.finalizeLine();
			
			table.initNewLine();
			for (int y=0; y<startYears.length; y++) {
				String yearPrefix = prefix+"_"+startYears[y]+"_"+endYear;
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				estCmlMFD.setName("Estimated");
				funcs.add(estCmlMFD);
				chars.add(estPrefChar);
				DiscretizedFunc estUpper = estCmlMFD.getUpper();
				estUpper = estUpper.deepClone();
				estUpper.setName("Estimated "+estCmlMFD.getBoundName());
				funcs.add(estUpper);
				chars.add(estBoundsChar);
				DiscretizedFunc estLower = estCmlMFD.getLower();
				estLower = estLower.deepClone();
				estLower.setName(null);
				funcs.add(estLower);
				chars.add(estBoundsChar);
				
				obsMFDs[y].setName("Observed");
				funcs.add(obsMFDs[y]);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
				bounds95[y].setName("Observed 95% and 68% Bounds");
				funcs.add(bounds95[y]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, obsTransColor));
				bounds68[y].setName(null);
				funcs.add(bounds68[y]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, obsTransColor));
				
				
				baCml.setName("Branch-Averaged Supra-Seis");
				funcs.add(baCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.DARK_GRAY));
				
				if (baTotCml != null) {
					baTotCml.setName("Branch-Averaged Total");
					funcs.add(baTotCml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.DARK_GRAY));
				}
				
				if (nodes.size() > 1) {
					CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, nodes.size()-1d);
					
					for (int n=0; n<nodeFuncs.size(); n++) {
						EvenlyDiscretizedFunc nodeCml = nodeFuncs.get(n);
						funcs.add(nodeCml);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt.getColor((float)n)));
					}
				}
				
				ArbitrarilyDiscretizedFunc obsClone = obsMFDs[y].deepClone();
				obsClone.setName(null);
				funcs.add(obsClone);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
				
				PlotSpec cmlSpec = new PlotSpec(funcs, chars,
						wusReg.getName(), "Magnitude", "Cumulative Rate (/yr)");
				cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.getPlotPrefs().scaleFontSizes(1.2d);
				
				gp.drawGraphPanel(cmlSpec, false, true, magRage, yRange);
				PlotUtils.writePlots(resourcesDir, yearPrefix+"_cml", gp, 950, 800, true, true, true);
				table.addColumn("![MFD]("+resourcesDir.getName()+"/"+yearPrefix+"_cml.png)");
			}
			table.finalizeLine();
			
			table = table.wrap(2, 0);
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}

}
