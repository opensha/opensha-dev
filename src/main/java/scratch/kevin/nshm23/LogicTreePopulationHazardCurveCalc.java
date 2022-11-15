package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedDiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.binFile.GeolocatedRectangularBinaryMesh2DCalculator;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.hazard.LogicTreeHazardCompare;
import org.opensha.sha.earthquake.faultSysSolution.hazard.SiteLogicTreeHazardPageGen;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;
import com.google.common.io.LittleEndianDataInputStream;

public class LogicTreePopulationHazardCurveCalc {
	
	static DecimalFormat popDF = new DecimalFormat("0");
	static {
		popDF.setGroupingUsed(true);
		popDF.setGroupingSize(3);
	}

	public static void main(String[] args) throws IOException {
		Preconditions.checkArgument(args.length == 5 || args.length == 8,
				"USAGE: <primary-results-zip> <primary-hazard-zip> <primary-name> [<comparison-results-zip> "
				+ "<comparison-hazard-zip> <comparison-name>] <output-dir> <population-data-file>");
		
		int cnt = 0;
		File resultsFile = new File(args[cnt++]);
		File hazardFile = new File(args[cnt++]);
		String mainName = args[cnt++];
		String compName;
		File compResultsFile, compHazardFile;
		if (args.length > 5) {
			if (args[cnt].toLowerCase().trim().equals("null")) {
				compResultsFile = null;
				cnt++;
			} else {
				compResultsFile = new File(args[cnt++]);
			}
			compHazardFile = new File(args[cnt++]);
			compName = args[cnt++];
		} else {
			compResultsFile = null;
			compHazardFile = null;
			compName = null;
		}
		File outputDir = new File(args[cnt++]);
		File popDataFile = new File(args[cnt++]);
		
		SolutionLogicTree solTree = SolutionLogicTree.load(resultsFile);
		
		ReturnPeriods[] rps = ReturnPeriods.values();
		double[] periods = { 0d, 1d };
		double spacing = -1; // detect
//		double spacing = 0.1;
//		double spacing = 0.25;
		
		LogicTree<?> tree = solTree.getLogicTree();
//		if (currentWeights)
//			tree.setWeightProvider(new BranchWeightProvider.CurrentWeights());
		
		LogicTreeHazardCompare mapper = null;
		LogicTreeHazardCompare comp = null;
		int exit = 0;
		try {
			mapper = new LogicTreeHazardCompare(solTree, tree,
					hazardFile, rps, periods, spacing);
			
			LogicTree<?> compTree = null;
			if (compHazardFile != null) {
				SolutionLogicTree compSolTree;
				if (compResultsFile == null) {
					// just have an average hazard result, possibly an external ERF
					compSolTree = null;
					compTree = null;
				} else {
					compSolTree = SolutionLogicTree.load(compResultsFile);
					compTree = compSolTree.getLogicTree();
//					if (compCurrentWeights)
//						compTree.setWeightProvider(new BranchWeightProvider.CurrentWeights());
				}
				comp = new LogicTreeHazardCompare(compSolTree, compTree,
						compHazardFile, rps, periods, spacing);
			}
			
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
			
			buildPopulationCurveReport(mapper, mainName, tree, comp, compName, compTree, outputDir, popDataFile, rps, periods);
		} catch (Exception e) {
			e.printStackTrace();
			exit = 1;
		} finally {
			if (mapper != null)
				mapper.close();
			if (comp != null)
				comp.close();
		}
		System.exit(exit);
	}
	
	private static void buildPopulationCurveReport(LogicTreeHazardCompare mapper, String name, LogicTree<?> tree,
			LogicTreeHazardCompare comp, String compName, LogicTree<?> compTree, File outputDir, File popDataFile,
			ReturnPeriods[] rps, double[] periods) throws IOException {
		GriddedGeoDataSet populationXYZ = null;
		
		List<String> lines = new ArrayList<>();
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		lines.add("# "+name+" Population Hazard");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "_[(top)](#table-of-contents)_";
		
		Color primaryColor = Color.RED;
		Color compColor = Color.BLUE;
		
		double popScalar = 1e-6;
		String popUnits = "(Millions)";
		
		for (double period : periods) {
			String perLabel, perPrefix;
			if (period == 0d) {
				perLabel = "PGA (g)";
				perPrefix = "pga";
			} else {
				perLabel = (float)period+"s SA";
				perPrefix = (float)period+"s";
			}
			
			DiscretizedFunc xVals = null;
			
			for (ReturnPeriods rp : rps) {
				GriddedGeoDataSet[] maps = mapper.loadMaps(rp, period);
				Preconditions.checkNotNull(maps);
				for (int i=0; i<maps.length; i++)
					Preconditions.checkNotNull(maps[i], "map %s is null", i);
				GriddedRegion gridReg = maps[0].getRegion();
				
				if (xVals == null) {
					DiscretizedFunc rawXVals = mapper.getMapper().getXVals(period);
					xVals = new ArbitrarilyDiscretizedFunc();
					for (int i=1; i<rawXVals.size(); i++) {
						double x1 = rawXVals.getX(i-1);
						double x2 = rawXVals.getX(i);
						if (x1 == 0d)
							continue;
						EvenlyDiscretizedFunc logVals = new EvenlyDiscretizedFunc(Math.log10(x1), Math.log10(x2), 50);
						for (int j=0; j<logVals.size(); j++)
							xVals.set(Math.pow(10, logVals.getX(j)), 0d);
					}
				}
				
				if (populationXYZ == null) {
					// load it
					populationXYZ = loadMapPopulationData(popDataFile, gridReg);
					double maxZ = populationXYZ.getMaxZ();
					CPT popCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(1d, maxZ);
					popCPT.add(0, new CPTVal(0f, Color.WHITE, 1f, Color.WHITE));
					popCPT.setBelowMinColor(Color.WHITE);
					SolHazardMapCalc hazMapper = mapper.getMapper();
					hazMapper.plotMap(resourcesDir, "population_data", populationXYZ, popCPT, "Population Data",
							"Population Per "+(float)gridReg.getSpacing()+" Degree Grid Cell");
					
					lines.add("## Population Data");
					lines.add(topLink); lines.add("");
					
					// TODO don't hardcode?
					lines.add("2021 CONUS population data (nighttime) from [LandScan](https://landscan.ornl.gov), "
							+ "doi: [10.48690/1527701](https://doi.org/10.48690/1527701)");
					lines.add("");
					lines.add("Population data mapped to gridded region with "+(float)gridReg.getSpacing()+" degree "
							+ "grid cells. Total population in region: "+popDF.format(populationXYZ.getSumZ()));
					lines.add("");
					lines.add("![Population data map]("+resourcesDir.getName()+"/population_data.png)");
					lines.add("");
					
					// scale
					populationXYZ.scale(popScalar);
				}
				
				String label = perLabel+", "+rp.label;
				String prefix = perPrefix+"_"+rp.name();
				
				System.out.println(label);
				
				lines.add("## "+label);
				lines.add(topLink);
				lines.add("");
				
				GriddedGeoDataSet meanMap = mapper.buildMean(maps);
				
				DiscretizedFunc meanCurve = calcPopExceedanceCurve(meanMap, populationXYZ, xVals);
				meanCurve.setName(name);
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				funcs.add(meanCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, primaryColor));
				
				if (comp != null) {
					GriddedGeoDataSet[] compMaps = comp.loadMaps(rp, period);
					
					GriddedGeoDataSet compMeanMap = comp.buildMean(compMaps);
					
					DiscretizedFunc compMeanCurve = calcPopExceedanceCurve(compMeanMap, populationXYZ, xVals);
					compMeanCurve.setName(compName);
					
					funcs.add(compMeanCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, compColor));
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, "Population Hazard Curve", label, "Population Exceeding "+popUnits);
				spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				Range linearYRange, logYRange;
				
				double sumZ = populationXYZ.getSumZ();
				
				linearYRange = new Range(0d, sumZ*1.05);
				double maxY = Math.pow(10, Math.ceil(Math.log10(populationXYZ.getSumZ())));
				double minY;
				if (maxY > 1e-1)
					minY = 1e-3;
				else
					minY = Math.pow(10, Math.log10(maxY - 2));
				logYRange = new Range(minY, maxY);
				
				double minX = Double.POSITIVE_INFINITY;
				double maxX = 0d;
				for (DiscretizedFunc func : funcs) {
					for (Point2D pt : func) {
						if (pt.getX() > 1e-2 && (float)pt.getY() < (float)logYRange.getUpperBound()
								&& (float)pt.getY() > (float)logYRange.getLowerBound()) {
							minX = Math.min(minX, pt.getX());
							maxX = Math.max(maxX, pt.getX());
						}
					}
				}
				minX = Math.pow(10, Math.floor(Math.log10(minX)));
				maxX = Math.pow(10, Math.ceil(Math.log10(maxX)));
				Range xRange = new Range(minX, maxX);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				// write linear curve
				gp.drawGraphPanel(spec, true, false, xRange, linearYRange);
				PlotUtils.writePlots(resourcesDir, prefix, gp, 1000, 900, true, true, false);
				gp.drawGraphPanel(spec, true, true, xRange, logYRange);
				PlotUtils.writePlots(resourcesDir, prefix+"_log", gp, 1000, 900, true, true, false);
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine("Linear Curves", "Log10 Curves");
				table.initNewLine();
				table.addColumn("![Population Hazard Curve]("+resourcesDir.getName()+"/"+prefix+".png)");
				table.addColumn("![Population Hazard Curve]("+resourcesDir.getName()+"/"+prefix+"_log.png)");
				table.finalizeLine();
				lines.addAll(table.build());
				lines.add("");
				
				// now logic tree
				List<LogicTreeLevel<?>> levels = new ArrayList<>();
				List<List<LogicTreeNode>> levelNodes = new ArrayList<>();
				List<List<DiscretizedFunc>> levelNodeCurves = new ArrayList<>();
				for (int l=0; l<tree.getLevels().size(); l++) {
					List<LogicTreeNode> nodes = new ArrayList<>();
					List<DiscretizedFunc> nodeMeanCurves = new ArrayList<>();
					
					LogicTreeLevel<?> level = tree.getLevels().get(l);
					for (LogicTreeNode node : level.getNodes()) {
						List<GriddedGeoDataSet> nodeMaps = new ArrayList<>();
						List<Double> nodeWeights = new ArrayList<>();
						for (int i=0; i<tree.size(); i++) {
							LogicTreeBranch<?> branch = tree.getBranch(i);
							if (branch.hasValue(node)) {
								nodeMaps.add(maps[i]);
								nodeWeights.add(tree.getBranchWeight(i));
							}
						}
						
						if (!nodeMaps.isEmpty()) {
							GriddedGeoDataSet nodeMeanMap = mapper.buildMean(nodeMaps, nodeWeights);
							
							DiscretizedFunc nodeMeanCurve = calcPopExceedanceCurve(nodeMeanMap, populationXYZ, xVals);
							
							nodes.add(node);
							nodeMeanCurves.add(nodeMeanCurve);
						}
					}
					
					if (nodes.size() > 1) {
						levels.add(level);
						levelNodes.add(nodes);
						levelNodeCurves.add(nodeMeanCurves);
					}
				}
				if (!levels.isEmpty()) {
					lines.add("### "+label+", Logic Tree Comparisons");
					lines.add(topLink); lines.add("");
					
					for (int l=0; l<levels.size(); l++) {
						LogicTreeLevel<?> level = levels.get(l);
						List<LogicTreeNode> nodes = levelNodes.get(l);
						List<DiscretizedFunc> nodeCurves = levelNodeCurves.get(l);
						
						lines.add("#### "+level.getName()+", "+label);
						lines.add(topLink); lines.add("");
						
						table = MarkdownUtils.tableBuilder();
						table.addLine("Linear Curves", "Log10 Curves");
						
						funcs = new ArrayList<>();
						chars = new ArrayList<>();
						
						funcs.add(meanCurve);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
						
						List<Color> colors = SiteLogicTreeHazardPageGen.getNodeColors(nodes.size());
						
						for (int i=0; i<nodes.size(); i++) {
							DiscretizedFunc nodeCurve = nodeCurves.get(i);
							nodeCurve.setName(nodes.get(i).getShortName());
							
							funcs.add(nodeCurve);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors.get(i)));
						}
						
						spec = new PlotSpec(funcs, chars, level.getName(), label, "Population Exceeding "+popUnits);
						spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
						
						String ltPrefix = prefix+"_"+level.getShortName().replaceAll("\\W+", "_");
						
						table.initNewLine();
						
						gp.drawGraphPanel(spec, true, false, xRange, linearYRange);
						PlotUtils.writePlots(resourcesDir, ltPrefix, gp, 1000, 900, true, true, false);
						table.addColumn("!["+level.getName()+"]("+resourcesDir.getName()+"/"+ltPrefix+".png)");
						
						gp.drawGraphPanel(spec, true, true, xRange, logYRange);
						PlotUtils.writePlots(resourcesDir, ltPrefix+"_log", gp, 1000, 900, true, true, false);
						table.addColumn("!["+level.getName()+"]("+resourcesDir.getName()+"/"+ltPrefix+"_log.png)");
						
						table.finalizeLine();
						
						// now difference and ratio
						
						List<DiscretizedFunc> diffFuncs = new ArrayList<>();
						List<DiscretizedFunc> logRatioFuncs = new ArrayList<>();
						
						diffFuncs.add(flatCurve(xRange, 0d, name));
						logRatioFuncs.add(flatCurve(xRange, 0d, name));
						
						double maxDiff = 1e-3;
						
						for (int n=0; n<nodes.size(); n++) {
							DiscretizedFunc diffFunc = new ArbitrarilyDiscretizedFunc();
							DiscretizedFunc logRatioFunc = new ArbitrarilyDiscretizedFunc();
							
							DiscretizedFunc nodeFunc = nodeCurves.get(n);
							for (int i=0; i<xVals.size(); i++) {
								double x = xVals.getX(i);
								if (i >= nodeFunc.size() && i >= meanCurve.size()) {
									diffFunc.set(x, 0d);
									logRatioFunc.set(x, 0d);
									break;
								}
								double nodeVal, meanVal;
								if (i >= nodeFunc.size()) {
									// zero for this node, nonzero for mean
									nodeVal = 0d;
									meanVal = meanCurve.getY(i);
								} else if (i >= meanCurve.size()) {
									// nonzero for this node, nzero for mean
									nodeVal = nodeFunc.getY(i);
									meanVal = 0d;
								} else {
									nodeVal = nodeFunc.getY(i);
									meanVal = meanCurve.getY(i);
								}
								
								maxDiff = Math.max(maxDiff, Math.abs(nodeVal - meanVal));
								
								diffFunc.set(x, nodeVal - meanVal);
								if (nodeVal == meanVal)
									logRatioFunc.set(x, 0d);
								else
									logRatioFunc.set(x, Math.log10(nodeVal/meanVal));
							}
							
							String nodeName = nodes.get(n).getShortName();
							
							diffFunc.setName(nodeName);
							diffFuncs.add(diffFunc);
							
							logRatioFunc.setName(nodeName);
							logRatioFuncs.add(logRatioFunc);
						}
						
						Range diffRange = new Range(maxDiff*-1.05, maxDiff*1.05);
						
						table.addLine(MarkdownUtils.boldCentered("Difference"), MarkdownUtils.boldCentered("Ratio"));
						
						spec = new PlotSpec(diffFuncs, chars, level.getName(), label, "Population Exceeding, Choice - Mean "+popUnits);
						spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
						
						table.initNewLine();
						
						gp.drawGraphPanel(spec, true, false, xRange, diffRange);
						PlotUtils.writePlots(resourcesDir, ltPrefix+"_diff", gp, 1000, 800, true, true, false);
						table.addColumn("!["+level.getName()+"]("+resourcesDir.getName()+"/"+ltPrefix+"_diff.png)");
						
						spec = new PlotSpec(logRatioFuncs, chars, level.getName(), label, "Population Exceeding, Log10(Choice / Mean)");
						spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
						
						gp.drawGraphPanel(spec, true, false, xRange, new Range(-1d, 1d));
						PlotUtils.writePlots(resourcesDir, ltPrefix+"_ratio", gp, 1000, 800, true, true, false);
						table.addColumn("!["+level.getName()+"]("+resourcesDir.getName()+"/"+ltPrefix+"_ratio.png)");
						
						table.finalizeLine();
						
						lines.addAll(table.build());
						lines.add("");
					}
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static GriddedGeoDataSet loadMapPopulationData(File popDataFile, GriddedRegion gridReg) throws IOException {
		Preconditions.checkState(popDataFile.exists());
		File hdrFile, fltFile;
		String fName = popDataFile.getName();
		Preconditions.checkState(fName.contains("."), "Population file has no extension: %s", popDataFile);
		String prefix = fName.substring(0, fName.lastIndexOf('.'));
		if (fName.toLowerCase().endsWith(".hdr")) {
			// header file supplied
			hdrFile = popDataFile;
			fltFile = new File(popDataFile.getParentFile(), prefix+".flt");
			if (!fltFile.exists())
				// try upper case
				fltFile = new File(popDataFile.getParentFile(), prefix+".FLT");
			Preconditions.checkState(fltFile.exists(),
					"Corresponding float file not found (tried upper and lower case extensions): %s", fltFile);
		} else {
			// assume it's a float file 
			fltFile = popDataFile;
			hdrFile = new File(popDataFile.getParentFile(), prefix+".hdr");
			if (!hdrFile.exists())
				// try upper case
				hdrFile = new File(popDataFile.getParentFile(), prefix+".HDR");
			Preconditions.checkState(hdrFile.exists(),
					"Corresponding header file not found (tried upper and lower case extensions): %s", hdrFile);
		}
		
		System.out.println("Reading population data header from "+hdrFile);
		GeolocatedRectangularBinaryMesh2DCalculator meshCalc = GeolocatedRectangularBinaryMesh2DCalculator.readHDR(hdrFile);
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		
		long ny = meshCalc.getNY();
		long nx = meshCalc.getNX();
		long nxy = nx*ny;
		System.out.println("Population mesh size: "+nx+" x "+ny+" = "+nxy);
		
		System.out.println("Reading population data from "+fltFile);
		BufferedInputStream bin = new BufferedInputStream(new FileInputStream(fltFile));
		LittleEndianDataInputStream din = new LittleEndianDataInputStream(bin);
		
		long fileLocsUsed = 0l;
		BitSet gridLocsUsed = new BitSet(gridReg.getNodeCount());
		
		double sumMappedPop = 0d;
		double sumTotalPop = 0d;
		
		for (long index=0; index<nxy; index++) {
			float val = din.readFloat();
			
			if (val >= 0f) {
				sumTotalPop += val;
				
				Location loc = meshCalc.getLocationForMeshIndex(index);
				int locIndex = gridReg.indexForLocation(loc);
				if (locIndex >= 0) {
					sumMappedPop += val;
					
					fileLocsUsed++;
					gridLocsUsed.set(locIndex);
					xyz.set(locIndex, xyz.get(locIndex)+val);
				}
			}
		}
		
		din.close();
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		
		System.out.println("Used "+fileLocsUsed+"/"+nxy+" ("+pDF.format((double)fileLocsUsed/(double)nxy)+") population data indexes");
		int gridNodesSet = gridLocsUsed.cardinality();
		System.out.println("Set "+gridNodesSet+"/"+xyz.size()+" ("+pDF.format((double)gridNodesSet/(double)xyz.size())+") grid locations");
		System.out.println("Total population in dataset: "+popDF.format(sumTotalPop));
		System.out.println("Population in grid region: "+popDF.format(sumMappedPop)+" ("+pDF.format(sumMappedPop/sumTotalPop)+")");
		
		return xyz;
	}
	
	private static DiscretizedFunc calcPopExceedanceCurve(GriddedGeoDataSet map, GriddedGeoDataSet populationXYZ,
			DiscretizedFunc xVals) {
		Preconditions.checkState(populationXYZ.size() == map.size());
		
		double[] xValArray = new double[xVals.size()];
		for (int i=0; i<xValArray.length; i++)
			xValArray[i] = xVals.getX(i);
		double[] yVals = new double[xVals.size()];
		
		for (int i=0; i<map.size(); i++) {
			double val = map.get(i);
			double pop = populationXYZ.get(i);
			if (val > 0d && pop > 0d) {
				for (int j=0; j<xValArray.length; j++) {
					if (val >= xValArray[j])
						yVals[j] += pop;
					else
						break;
				}
			}
		}
		
//		return new LightFixedXFunc(xValArray, yVals);
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		
		for (int i=0; i<yVals.length; i++) {
			ret.set(xValArray[i], yVals[i]);
			if (yVals[i] == 0d)
				break;
		}
		
		return ret;
	}
	
	private static DiscretizedFunc flatCurve(Range xRange, double y, String name) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		ret.set(xRange.getLowerBound(), y);
		ret.set(xRange.getUpperBound(), y);
		return ret;
	}

}
