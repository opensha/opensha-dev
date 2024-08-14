package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedDiscretizedFunc;
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
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.NSHM23_BaseRegion;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Table;

import gov.usgs.earthquake.nshmp.model.HazardModel;

public class ObservedUncertainMFDComparisons {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/data/kevin/nshm23/batch_inversions/"
//				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//		new NSHM23_InvConfigFactory.NSHM23_V2(); // set seis values to V2
		File nodeDir = new File(invDir, "node_branch_averaged");
		File meanSolFile = new File(invDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File modelDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.b.3");
		File ltFile = new File(invDir, "logic_tree.json");

//		boolean wus = true;
//		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
//		File outputDir = new File(invDir, "observed_mfd_uncertainty_comparison");
//		int[] startYears = { 1850, 1900, 1930 };
//		String dataPrefix = "WUS-Other-";
		
		boolean wus = false;
		Region region = AnalysisRegions.CONUS_EAST.load();
		File outputDir = new File(invDir, "observed_mfd_uncertainty_comparison_ceus");
		int[] startYears = { 1800, 1850 };
		String dataPrefix = "CEUS-Other-";
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		String pathPrefix = "/data/erf/nshm23/seismicity/observed";
		
		
		Range magRage = new Range(6d, 8.5d);
//		Range yRange = wus ? new Range(1e-6, 2e0) : new Range(1e-6, 2e-1);
		Range yRange = new Range(1e-6, 2e0);
		EvenlyDiscretizedFunc refIncrMFD = FaultSysTools.initEmptyMFD(8.45);
		
		int endYear = 2023;
		
		int maxZeroBins = 2;
		
		ArbitrarilyDiscretizedFunc[] obsMFDs = new ArbitrarilyDiscretizedFunc[startYears.length];
		UncertainArbDiscFunc[] bounds95 = new UncertainArbDiscFunc[startYears.length];
		UncertainArbDiscFunc[] bounds68 = new UncertainArbDiscFunc[startYears.length];
		
		for (int i=0; i<startYears.length; i++) {
			String path = pathPrefix+"/"+dataPrefix+startYears[i]+"-"+endYear+"-direct-rates-v2.txt";
			System.out.println("Loading "+path);
			
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
						double duration = Double.parseDouble(split[index++]);
						double rate = Double.parseDouble(split[index++]);
						double lower2p5 = Double.parseDouble(split[index++]);
						double lower16 = Double.parseDouble(split[index++]);
						double upper84 = Double.parseDouble(split[index++]);
						double upper97p5 = Double.parseDouble(split[index++]);
						
						if (count == 0 && maxZeroBins == 0)
							break;
						
						// originally using Andy's rate which is (1+obs)/duration
//						obs.set(mag, count>0 ? rate : 0d);
						// now using direct obs/duration (maximum liklihood)
						obs.set(mag, (double)count/duration);
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
		
		UncertainBoundedIncrMagFreqDist estIncrMFD = NSHM23_RegionalSeismicity.getRemapped(region,
				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refIncrMFD, refIncrMFD.getMaxX());
		EvenlyDiscretizedFunc cmlMean = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFD);
		EvenlyDiscretizedFunc cmlUpper = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFD.getUpper());
		EvenlyDiscretizedFunc cmlLower = Regional_MFD_Plots.getCmlAsFakeIncr(estIncrMFD.getLower());
		UncertainArbDiscFunc estCmlMFD = new UncertainArbDiscFunc(cmlMean, cmlLower, cmlUpper, estIncrMFD.getBoundType());
		
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
		
		DiscretizedFunc modelSupraCml;
		DiscretizedFunc modelTotCml = null;
		UncertainBoundedDiscretizedFunc modelTotBounds = null;
		if (wus) {
			FaultSystemSolution baSol = FaultSystemSolution.load(meanSolFile);	
			
			LogicTree<?> tree = LogicTree.read(ltFile);
			for (LogicTreeLevel<?> level : tree.getLevels())
				levels.add(level);
			
			IncrementalMagFreqDist modelSupraIncr = baSol.calcNucleationMFD_forRegion(
					region, refIncrMFD.getMinX(), refIncrMFD.getMaxX(), refIncrMFD.getDelta(), false);
			modelSupraCml = modelSupraIncr.getCumRateDistWithOffset();
			
			GridSourceProvider gridProv = baSol.getGridSourceProvider();
			if (gridProv != null) {
				IncrementalMagFreqDist baTotIncr = modelSupraIncr.deepClone();
				for (int i=0; i<gridProv.getNumLocations(); i++) {
					Location loc = gridProv.getGriddedRegion().getLocation(i);
					if (region.contains(loc)) {
						IncrementalMagFreqDist gridMFD = gridProv.getMFD(i);
						if (gridMFD != null) {
							for (Point2D pt : gridMFD)
								baTotIncr.add(baTotIncr.getClosestXIndex(pt.getX()), pt.getY());
						}
					}
				}
				modelTotCml = baTotIncr.getCumRateDistWithOffset();
			}
			
			BranchRegionalMFDs branchMFDs = baSol.getModule(BranchRegionalMFDs.class);
			RegionsOfInterest rois = baSol.getRupSet().getModule(RegionsOfInterest.class);
			
			if (rois != null && branchMFDs != null && branchMFDs.hasRegionalMFDs()
					&& modelTotCml != null && branchMFDs.hasMFDs(MFDType.SUM)) {
				// see if we can add uncertainties
				int exactMatchRegIndex = -1;
				int imwIndex = -1;
				int pnwIndex = -1;
				int u3Index = -1;
				
				List<Region> solRegions = rois.getRegions();
				System.out.println("Looking for compatable region for "+region.getName());
				for (int i=0; i<solRegions.size(); i++) {
					Region solRegion = solRegions.get(i);
					boolean match = solRegion.equalsRegion(region);
					System.out.println("\t"+i+". "+solRegion.getName()+" ? "+match);
					if (match)
						exactMatchRegIndex = i;
					if (solRegion.equals(AnalysisRegions.CONUS_IMW.load()))
						imwIndex = i;
					if (solRegion.equals(AnalysisRegions.CONUS_PNW.load()))
						pnwIndex = i;
					if (solRegion.equals(AnalysisRegions.CONUS_U3_RELM.load()))
						u3Index = i;
				}
				double[] fracts = {0.025, 0.5, 0.975};
				
				if (exactMatchRegIndex >= 0) {
					System.out.println("Matched region with index "+exactMatchRegIndex);
					
					EvenlyDiscretizedFunc[] bounds = branchMFDs.calcRegionalCumulativeFractiles(
							MFDType.SUM, exactMatchRegIndex, fracts);
					modelTotBounds = new UncertainArbDiscFunc(bounds[1], bounds[0], bounds[2], UncertaintyBoundType.CONF_95);
				} else if (imwIndex >= 0 && pnwIndex >= 0 && u3Index >= 0) {
					System.out.println("Summing WUS sub-regions");
					
					EvenlyDiscretizedFunc[] bounds1 = branchMFDs.calcRegionalCumulativeFractiles(
							MFDType.SUM, imwIndex, fracts);
					EvenlyDiscretizedFunc[] bounds2 = branchMFDs.calcRegionalCumulativeFractiles(
							MFDType.SUM, pnwIndex, fracts);
					EvenlyDiscretizedFunc[] bounds3 = branchMFDs.calcRegionalCumulativeFractiles(
							MFDType.SUM, u3Index, fracts);
					
					EvenlyDiscretizedFunc[] bounds = new EvenlyDiscretizedFunc[fracts.length];
					for (int n=0; n<bounds.length; n++) {
						bounds[n] = new EvenlyDiscretizedFunc(
								bounds1[n].getMinX(), bounds1[n].size(), bounds1[n].getDelta());
						for (int i=0; i<bounds[n].size(); i++)
							bounds[n].set(i, bounds1[n].getY(i) + bounds2[n].getY(i) + bounds3[n].getY(i));
					}
					
					modelTotBounds = new UncertainArbDiscFunc(bounds[1], bounds[0], bounds[2], UncertaintyBoundType.CONF_95);
				} else {
					System.out.println("No matching region found!");
				}
			}
		} else {
			// conus east
			NSHM23_BaseRegion baseReg = AnalysisRegions.CONUS_EAST;
			Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW, TectonicRegionType.STABLE_SHALLOW);
			HazardModel model = HazardModel.load(modelDir.toPath());
			Table<NSHM23_BaseRegion, MFDType, IncrementalMagFreqDist> mfds = Regional_MFD_Plots.calcModelMFDs(
					model, trts, new NSHM23_BaseRegion[] {baseReg}, refIncrMFD);
			modelSupraCml = mfds.get(baseReg, MFDType.SUPRA_ONLY).getCumRateDistWithOffset();
			modelTotCml = mfds.get(baseReg, MFDType.SUM).getCumRateDistWithOffset();
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
							region, refIncrMFD.getMinX(), refIncrMFD.getMaxX(), refIncrMFD.getDelta(), false);
					
					GridSourceProvider gridProv = nodeSol.getGridSourceProvider();
					if (gridProv != null) {
						// add gridded seismicity
						for (int i=0; i<gridProv.getNumLocations(); i++) {
							Location loc = gridProv.getGriddedRegion().getLocation(i);
							if (region.contains(loc)) {
								IncrementalMagFreqDist gridMFD = gridProv.getMFD(i);
								if (gridMFD != null) {
									for (Point2D pt : gridMFD)
										nodeIncr.add(nodeIncr.getClosestXIndex(pt.getX()), pt.getY());
								}
							}
						}
					}
					
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
				
				estCmlMFD.setName("Extrapolated");
				funcs.add(estCmlMFD);
				chars.add(estPrefChar);
				DiscretizedFunc estUpper = estCmlMFD.getUpper();
				estUpper = estUpper.deepClone();
				estUpper.setName("Extrapolated "+estCmlMFD.getBoundName());
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
//				if (modelTotBounds == null) {
//					// use both 95 and 86
//					bounds95[y].setName("Observed 95% and 68% bounds");
//					funcs.add(bounds95[y]);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, obsTransColor));
//					bounds68[y].setName(null);
//					funcs.add(bounds68[y]);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, obsTransColor));
//				} else {
					// only use 95
					bounds95[y].setName("Observed 95% bounds");
					funcs.add(bounds95[y]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, obsTransColor));
//				}
				
				modelSupraCml.setName("Model Supra-Seis");
				funcs.add(modelSupraCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.DARK_GRAY));
				
				if (modelTotCml != null) {
					Color color = Color.DARK_GRAY;
					modelTotCml.setName("Model Total");
					funcs.add(modelTotCml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
					
					if (modelTotBounds != null) {
						Color transColor = new Color(color.getRed(), color.getGreen(), color.getBlue(), 80);
						modelTotBounds.setName("Model Total "+((UncertainBoundedDiscretizedFunc)modelTotBounds).getBoundName());
						funcs.add(modelTotBounds);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 3f, transColor));
					}
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
						region.getName(), "Magnitude", "Cumulative Rate (/yr)");
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
