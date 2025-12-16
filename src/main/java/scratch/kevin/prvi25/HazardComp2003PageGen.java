package scratch.kevin.prvi25;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.MarkdownUtils.TableTextAlignment;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

import scratch.kevin.latex.LaTeXUtils;
import scratch.kevin.nshm23.figures.MethodsAndIngredientsHazChangeFigures;
import scratch.kevin.nshm23.figures.WUS_HazardChangePageGen;
import scratch.kevin.prvi25.figures.MapSourceTypeDisagg;
import scratch.kevin.prvi25.figures.PRVI_Paths;

public class HazardComp2003PageGen {
	
	private static enum CompType {
			TOTAL("Full Model", "HazardChangeTotal"),
			CRUSTAL("Crustal-Only (Fault & Gridded)", "HazardChangeCrustal"),
			CRUSTAL_FAULTS("Crustal Fault-Only", "HazardChangeCrustalFault"),
			CRUSTAL_GRIDDED("Crustal Gridded-Only", "HazardChangeCrustalGridded"),
			INTERFACE("Interface-Only", "HazardChangeInterface"),
			SLAB("Intraslab-Only", "HazardChangeSlab");
		
		private String label;
		private String texPrefix;

		private CompType(String label, String texPrefix) {
			this.label = label;
			this.texPrefix = texPrefix;
		}
	};
	
	public static void main(String[] args) throws IOException {
//		String imtName = "PGA";
//		String imtDir = "PGA";
//		double period = 0d;
		
//		String imtName = "0.2s SA";
//		String imtDir = "SA0P2";
//		double period = 0.2d;
		
//		String imtName = "1s SA";
//		String imtDir = "SA1P0";
//		double period = 1d;
		
		String imtName = "5s SA";
		String imtDir = "SA5P0";
		double period = 5d;
		
//		ReturnPeriods[] rps = ReturnPeriods.values();
		ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50, ReturnPeriods.TEN_IN_50 };
		
		String name25 = "NSHM25";
		String name03 = "NSHM03";
		
		File combinedDir = PRVI_Paths.COMBINED_DIR;
		File crustalDir = PRVI_Paths.CRUSTAL_DIR;
		File subDir = PRVI_Paths.SUBDUCTION_DIR;
//		File combinedDir = new File(PRVI_Paths.INV_DIR, "2024_12_12-prvi25_crustal_subduction_combined_branches");
//		File crustalDir = new File(PRVI_Paths.INV_DIR, "2024_12_12-prvi25_crustal_branches-dmSample5x");
//		File subDir = new File(PRVI_Paths.INV_DIR, "2024_12_12-prvi25_subduction_branches");
		
		File inputDirFull25 = new File(PRVI_Paths.INV_DIR, combinedDir.getName()+"-ba_only-vs760");
		File inputDirCrustal25 = new File(PRVI_Paths.INV_DIR, crustalDir.getName()+"-ba_only-vs760");
		File inputDirInterface25 = new File(PRVI_Paths.INV_DIR, subDir.getName()+"-ba_only-INTERFACE_only-vs760");
		File inputDirSlab25 = new File(PRVI_Paths.INV_DIR, subDir.getName()+"-ba_only-SLAB_only-vs760");
		File outputDir = new File(combinedDir, "nshm03_erf_comparisons_"+imtDir);
		double sourceSpacing = 0.01;
		File sourcesDir03 = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
				+ "prvi-t14-2003-ERF-2025-GMMs-v1.7.15-2025conf-0p01ext-vs760-20250822-c4785428e56e5c/vs30-760/"+imtDir+"/source");

		System.out.println("Output dir: "+outputDir.getAbsolutePath());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		GriddedRegion mapReg = new GriddedRegion(PRVI25_RegionLoader.loadPRVI_MapExtents(),
				sourceSpacing, GriddedRegion.ANCHOR_0_0);
		GriddedGeoDataSet landMask = MapSourceTypeDisagg.buildLandMask(mapReg);
		
		File texFile = new File(resourcesDir, "stats.tex");
		FileWriter texFW = new FileWriter(texFile);
		
		DiscretizedFunc[] crustalFault25 = loadRegularCurves(inputDirCrustal25, IncludeBackgroundOption.EXCLUDE, mapReg, period);
		DiscretizedFunc[] crustalGrid25 = loadRegularCurves(inputDirCrustal25, IncludeBackgroundOption.ONLY, mapReg, period);
		DiscretizedFunc[] interface25 = loadRegularCurves(inputDirInterface25, IncludeBackgroundOption.INCLUDE, mapReg, period);
		DiscretizedFunc[] slab25 = loadRegularCurves(inputDirSlab25, IncludeBackgroundOption.ONLY, mapReg, period);
		DiscretizedFunc[] full25 = loadRegularCurves(inputDirFull25, IncludeBackgroundOption.INCLUDE, mapReg, period);
		
		DiscretizedFunc[] crustalFault03 = add(
				loadExtCurves(new File(sourcesDir03, "FAULT/curves.csv"), mapReg),
				loadExtCurves(new File(sourcesDir03, "ZONE/curves.csv"), mapReg)
				);
		DiscretizedFunc[] crustalGrid03 = loadExtCurves(new File(sourcesDir03, "GRID/curves.csv"), mapReg);
		DiscretizedFunc[] interface03 = loadExtCurves(new File(sourcesDir03, "INTERFACE/curves.csv"), mapReg);
		DiscretizedFunc[] slab03 = loadExtCurves(new File(sourcesDir03, "SLAB/curves.csv"), mapReg);
		DiscretizedFunc[] full03 = add(crustalFault03, crustalGrid03, interface03, slab03);
		
		// convert to probabilities
		crustalFault03 = ratesToProbs(crustalFault03);
		crustalGrid03 = ratesToProbs(crustalGrid03);
		interface03 = ratesToProbs(interface03);
		slab03 = ratesToProbs(slab03);
		full03 = ratesToProbs(full03);
		
		// interpolate 03 onto our x-value spacing
		crustalFault03 = interpolateToMatch(crustalFault03, full25[0]);
		crustalGrid03 = interpolateToMatch(crustalGrid03, full25[0]);
		interface03 = interpolateToMatch(interface03, full25[0]);
		slab03 = interpolateToMatch(slab03, full25[0]);
		full03 = interpolateToMatch(full03, full25[0]);
		
		// these maps substitute a single nshm03 ingredient at a time, keeping everything else at 25
		DiscretizedFunc[] withCrustalFaults03 = add(crustalFault03, crustalGrid25, interface25, slab25);
		DiscretizedFunc[] withCrustalGrid03 = add(crustalFault25, crustalGrid03, interface25, slab25);
		DiscretizedFunc[] withCrustal03 = add(crustalFault03, crustalGrid03, interface25, slab25);
		DiscretizedFunc[] withInterface03 = add(crustalFault25, crustalGrid25, interface03, slab25);
		DiscretizedFunc[] withSlab03 = add(crustalFault25, crustalGrid25, interface25, slab03);
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(List.of(), mapReg);
		mapMaker.setPoliticalBoundaryChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
//		mapMaker.setDefaultPlotWidth(1000);
		
		Color transparent = new Color(255, 255, 255, 0);
		
//		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-1.5, 0.5);
		hazCPT.setLog10(true);
		hazCPT.setNanColor(transparent);
		
//		CPT pDiffCPT = MethodsAndIngredientsHazChangeFigures.getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance(), 10d, 50d);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(transparent);
		
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2d, 0.2d);
		diffCPT.setNanColor(transparent);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Hazard Comparisons, "+name25+" vs "+name03);
		lines.add("");
		
		lines.add("This page compares "+name25+" and "+name03+" hazard maps due to the ERF. "
				+ "Each individual ERF contributor to hazard changes are plotted separately, but GMMs are held constant.");
		lines.add("");
		
		String csvStr = "Download map CSVs:";
		for (ReturnPeriods rp : rps) {
			System.out.println("Building CSV for "+rp);
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Location Index", "Latitude", "Longitude", name25, name03,
					name03+" Crustal Faults", name03+" Crustal Gridded", name03+" Crustal",
					name03+" Interface", name03+" Slab");
			
			GriddedGeoDataSet[] maps = {
					curvestoMap(full25, mapReg, rp),
					curvestoMap(full03, mapReg, rp),
					curvestoMap(withCrustalFaults03, mapReg, rp),
					curvestoMap(withCrustalGrid03, mapReg, rp),
					curvestoMap(withCrustal03, mapReg, rp),
					curvestoMap(withInterface03, mapReg, rp),
					curvestoMap(withSlab03, mapReg, rp)
			};
			
			for (int i=0; i<maps[0].size(); i++) {
				List<String> line = new ArrayList<>(csv.getNumCols());
				
				line.add(i+"");
				Location loc = maps[0].getLocation(i);
				line.add((float)loc.getLatitude()+"");
				line.add((float)loc.getLongitude()+"");
				for (GriddedGeoDataSet map : maps)
					line.add((float)map.get(i)+"");
				
				csv.addLine(line);
			}
			
			File csvFile = new File(resourcesDir, "maps_"+rp.name()+".csv");
			csv.writeToFile(csvFile);
			
			csvStr += " [__"+csvFile.getName()+"__]("+resourcesDir.getName()+"/"+csvFile.getName()+")";
		}
		lines.add(csvStr);
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Change Summaries");
		lines.add(topLink); lines.add("");
		int summaryTableIndex = lines.size();
		lines.add("");
		TableBuilder[] rpTables = new TableBuilder[rps.length];
		for (int r=0; r<rpTables.length; r++) {
			rpTables[r] = MarkdownUtils.tableBuilder().textAlign(TableTextAlignment.CENTER).initNewLine();
			for (CompType type : CompType.values())
				rpTables[r].addColumn("__"+type.label+"__");
			rpTables[r].finalizeLine().initNewLine();
		}
		
		for (CompType type : CompType.values()) {
			DiscretizedFunc[] curves25, curves03;
			String label = type.label+" Comparison";
			String texPrefix = type.texPrefix;
			String description;
			String mapLabelAdd;
			
			switch (type) {
			case TOTAL:
				curves25 = full25;
				curves03 = full03;
				description = "Full model comparison, including fault, gridded, and subduction model changes.";
				mapLabelAdd = "";
				break;
			case CRUSTAL:
				curves25 = full25;
				curves03 = withCrustal03;
				description = "Crustal-only comparison, holding subduction sources constant (using those from "+name25+").";
				mapLabelAdd = "Crustal, ";
				break;
			case CRUSTAL_FAULTS:
				curves25 = full25;
				curves03 = withCrustalFaults03;
				description = "Crustal fault-only comparison, holding gridded and subduction sources constant (using those from "+name25+").";
				mapLabelAdd = "Crustal Faults, ";
				break;
			case CRUSTAL_GRIDDED:
				curves25 = full25;
				curves03 = withCrustalGrid03;
				description = "Crustal Gridded-only comparison, holding crustal fault and subduction sources constant (using those from "+name25+").";
				mapLabelAdd = "Crustal Gridded, ";
				break;
			case INTERFACE:
				curves25 = full25;
				curves03 = withInterface03;
				description = "Subduction interface-only comparison, crustal and slab sources constant (using those from "+name25+").";
				mapLabelAdd = "Interface, ";
				break;
			case SLAB:
				curves25 = full25;
				curves03 = withSlab03;
				description = "Subduction intraslab-only comparison, crustal and interface sources constant (using those from "+name25+").";
				mapLabelAdd = "Intraslab, ";
				break;

			default:
				throw new IllegalStateException();
			}
			
			lines.add("## "+label);
			lines.add(topLink); lines.add("");
			
			lines.add(description);
			lines.add("");
			
			for (int r=0; r<rps.length; r++) {
				ReturnPeriods rp = rps[r];
				System.out.println("Plotting "+label+", "+rp.label);
				
				lines.add("### "+label+", "+rp.label);
				lines.add(topLink); lines.add("");
				
				String hazLabel = imtName+", "+rp.label;
				String prefix = type.name()+"_"+rp.name();
				
				GriddedGeoDataSet map25 = curvestoMap(curves25, mapReg, rp);
				GriddedGeoDataSet map03 = curvestoMap(curves03, mapReg, rp);
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				
				table.addLine(name25, name03);
				
				table.initNewLine();
				
				mapMaker.plotXYZData(map25, hazCPT, mapLabelAdd+name25+", "+hazLabel+" (g)");
				mapMaker.plot(resourcesDir, prefix+"_nshm25", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_nshm25.png)");
				mapMaker.plotXYZData(map03, hazCPT, mapLabelAdd+name03+", "+hazLabel+" (g)");
				mapMaker.plot(resourcesDir, prefix+"_nshm03", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_nshm03.png)");
				
				table.finalizeLine();
				
				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				GriddedGeoDataSet pDiff = WUS_HazardChangePageGen.mapPDiff(map25, map03);
				GriddedGeoDataSet diff = WUS_HazardChangePageGen.mapDiff(map25, map03);
				
				table.initNewLine();
				
				mapMaker.plotXYZData(pDiff, pDiffCPT, mapLabelAdd+name25+" vs "+name03+", % Change, "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_pDiff.png)");
				mapMaker.plotXYZData(diff, diffCPT, mapLabelAdd+name25+" - "+name03+", "+hazLabel+" (g)");
				mapMaker.plot(resourcesDir, prefix+"_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_diff.png)");
				
				table.finalizeLine();
				
				rpTables[r].addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_pDiff.png)");
				
				lines.addAll(table.build());
				lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				table.addLine("", "Average % Change", "% Change Range", "Average Difference", "Difference Range");
				
				for (boolean land : new boolean[] {false,true}) {
					for (boolean abs : new boolean[] {false,true}) {
						MinMaxAveTracker diffStats = mapStats(diff, abs, land ? landMask : null);
						MinMaxAveTracker pDiffStats = mapStats(pDiff, abs, land ? landMask : null);
						String tableLabel = land ? "On land" : "Full map";
						if (abs)
							tableLabel += ", absolute value";
						table.addLine("__"+tableLabel+"__",
								withPercentAndPlusIfPositive(pDiffStats.getAverage()),
								"["+withPercentAndPlusIfPositive(pDiffStats.getMin())+", "+withPercentAndPlusIfPositive(pDiffStats.getMax())+"]",
								twoDigits.format(diffStats.getAverage()),
								"["+twoDigits.format(diffStats.getMin())+", "+twoDigits.format(diffStats.getMax())+"]");
						
						if (land && rp == ReturnPeriods.TWO_IN_50) {
							String myPrefix = texPrefix;
							if (abs)
								myPrefix += "Abs";
							texFW.write(LaTeXUtils.defineValueCommand(myPrefix+"Average", withPercent(pDiffStats.getAverage()))+"\n");
							texFW.write(LaTeXUtils.defineValueCommand(myPrefix+"Min", withPercent(pDiffStats.getMin()))+"\n");
							texFW.write(LaTeXUtils.defineValueCommand(myPrefix+"Max", withPercent(pDiffStats.getMax()))+"\n");
						}
					}
				}
				
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		List<String> rpSummaryLines = new ArrayList<>();
		for (int r=0; r<rps.length; r++) {
			rpSummaryLines.add("### "+rps[r].label+" Summary");
			rpSummaryLines.add(topLink); rpSummaryLines.add("");
			rpTables[r].finalizeLine();
			rpSummaryLines.addAll(rpTables[r].wrap(2, 0).build());
			rpSummaryLines.add("");
		}
		lines.addAll(summaryTableIndex, rpSummaryLines);
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
		
		texFW.close();
	}
	
	private static final DecimalFormat oneDigit = new DecimalFormat("0.0");
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	
	private static String withPercent(double percent) {
		return oneDigit.format(percent)+"%";
	}
	
	private static String withPercentAndPlusIfPositive(double percent) {
		String ret = withPercent(percent);
		if (percent > 0)
			ret = "+"+ret;
		return ret;
	}
	
	private static MinMaxAveTracker mapStats(GriddedGeoDataSet map, boolean abs, GriddedGeoDataSet mask) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (int i=0; i<map.size(); i++) {
			if (mask != null && mask.get(i) != 1d)
				continue;
			double val = map.get(i);
			if (abs)
				val = Math.abs(val);
			if (Double.isFinite(val))
				track.addValue(val);
		}
		return track;
	}
	
	private static DiscretizedFunc[] loadRegularCurves(File runDir, IncludeBackgroundOption bgOp,
			GriddedRegion gridReg, double period) throws IOException {
		String hazardPrefix = "hazard_"+(float)gridReg.getSpacing()+"deg";
		hazardPrefix += "_grid_seis_";
		hazardPrefix += bgOp.name();
		File resultsDir = new File(runDir, "results");
		File subDir = new File(resultsDir, hazardPrefix);
		Preconditions.checkState(subDir.exists(), "Doesn't exist: %s", subDir.getAbsolutePath());
		
		File curvesFile = new File(subDir, SolHazardMapCalc.getCSV_FileName("curves", period));
		if (!curvesFile.exists())
			curvesFile = new File(curvesFile.getAbsolutePath()+".gz");
		Preconditions.checkState(curvesFile.exists(), "Doesn't exist: %s", curvesFile.getAbsolutePath());
		
		return SolHazardMapCalc.loadCurvesCSV(CSVFile.readFile(curvesFile, true), gridReg);
	}
	
	private static DiscretizedFunc[] loadExtCurves(File csvFile, GriddedRegion gridReg) throws IOException {
		System.out.println("Loading external curves from "+csvFile.getAbsolutePath());
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		double[] xVals = new double[csv.getNumCols()-2];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = csv.getDouble(0, i+2);
		
		DiscretizedFunc[] curves = new DiscretizedFunc[gridReg.getNodeCount()];
		
		int numSkipped = 0;
		int numFilled = 0;
		for (int row=1; row<csv.getNumRows(); row++) {
			double lon = csv.getDouble(row, 0);
			double lat = csv.getDouble(row, 1);
			Location loc = new Location(lat, lon);
			int index = gridReg.indexForLocation(loc);
			if (index < 0) {
				numSkipped++;
			} else {
				numFilled++;
				double[] yVals = new double[xVals.length];
				for (int i=0; i<yVals.length; i++)
					yVals[i] = csv.getDouble(row, i+2);
				curves[index] = new LightFixedXFunc(xVals, yVals);
			}
		}
		System.out.println("\tFilled in "+numFilled+"/"+gridReg.getNodeCount()+" curves (skipped "+numSkipped+")");
		
		return curves;
	}
	
	private static double[] getXVals(DiscretizedFunc[] curves) {
		for (DiscretizedFunc curve : curves) {
			if (curve != null) {
				double[] xVals = new double[curve.size()];
				for (int i=0; i<xVals.length; i++)
					xVals[i] = curve.getX(i);
				return xVals;
			}
		}
		throw new IllegalStateException("All curves are null");
	}
	
	private static DiscretizedFunc[] ratesToProbs(DiscretizedFunc[] curves) {
		DiscretizedFunc[] ret = new DiscretizedFunc[curves.length];
		double[] xVals = getXVals(curves);
		for (int n=0; n<curves.length; n++) {
			DiscretizedFunc curve = curves[n];
			if (curve == null)
				continue;
			double[] probs = new double[xVals.length];
			for (int i=0; i<probs.length; i++)
				probs[i] = 1d - Math.exp(-curve.getY(i));
			ret[n] = new LightFixedXFunc(xVals, probs);
		}
		return ret;
	}
	
	private static DiscretizedFunc[] interpolateToMatch(DiscretizedFunc[] curves, DiscretizedFunc targetSpacing) {
		DiscretizedFunc[] ret = new DiscretizedFunc[curves.length];
		double[] xVals = new double[targetSpacing.size()];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = targetSpacing.getX(i);
		for (int n=0; n<curves.length; n++) {
			DiscretizedFunc curve = curves[n];
			if (curve == null)
				continue;
			double[] interpY = new double[xVals.length];
			for (int i=0; i<interpY.length; i++) {
				double x = xVals[i];
				if (x < curve.getMinX())
					interpY[i] = curve.getY(0);
				else if (x > curve.getMaxX())
					interpY[i] = 0d;
				else
					interpY[i] = curve.getInterpolatedY_inLogXDomain(x);
			}
			ret[n] = new LightFixedXFunc(xVals, interpY);
		}
		return ret;
	}
	
	private static DiscretizedFunc[] add(DiscretizedFunc[]... allCurves) {
		DiscretizedFunc[] ret = new DiscretizedFunc[allCurves[0].length];
		
		double[] xVals = getXVals(allCurves[0]);
		
		for (int i=0; i<ret.length; i++) {
			double[] yVals = new double[xVals.length];
			for (DiscretizedFunc[] curves : allCurves) {
				DiscretizedFunc curve = curves[i];
				if (curve == null || yVals == null) {
					yVals = null;
					continue;
				}
				Preconditions.checkState(curve.size() == xVals.length,
						"X-value count mismatch: %s vs %s; minX: %s vs %s; maxX: %s vs %s",
						curve.size(), xVals.length, (float)curve.getMinX(), (float)xVals[0],
						(float)curve.getMaxX(), (float)xVals[xVals.length-1]);
				for (int j=0; j<xVals.length; j++) {
					Preconditions.checkState((float)xVals[j] == (float)curve.getX(j));
					yVals[j] += curve.getY(j);
				}
			}
			if (yVals != null)
				ret[i] = new LightFixedXFunc(xVals, yVals);
		}
		
		return ret;
	}
	
	private static GriddedGeoDataSet curvestoMap(DiscretizedFunc[] curves, GriddedRegion gridRegion, ReturnPeriods rp) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(gridRegion);
		for (int i=0; i<ret.size(); i++)
			ret.set(i, Double.NaN);
		
		for (int i=0; i<curves.length; i++) {
			DiscretizedFunc curve = curves[i];
			double val;
			if (curve == null)
				val = Double.NaN;
			else if (rp.oneYearProb > curve.getMaxY())
				val = 0d;
			else if (rp.oneYearProb < curve.getMinY())
				// saturated
				val = curve.getMaxX();
			else
				val = curve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
			ret.set(i, val);
		}
		
		return ret;
	}

}
