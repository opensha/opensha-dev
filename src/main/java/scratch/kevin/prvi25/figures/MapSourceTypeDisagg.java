package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.Interpolate;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.latex.LaTeXUtils;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class MapSourceTypeDisagg {
	
	enum MapType {
		COMBINED("Complete model", Colors.tab_brown),
		CRUSTAL("Crustal sources", Colors.tab_blue),
		SUBDUCTION("All subduction sources", Colors.tab_red),
		SUBDUCTION_INTERFACE("Subduction interface sources", Colors.tab_pink),
		SUBDUCTION_SLAB("Subduction intraslab sources", Colors.tab_green, IncludeBackgroundOption.ONLY);
		
		private String name;
		private Color color;
		private IncludeBackgroundOption[] bgOps;
		
		private MapType(String name, Color color) {
			this(name, color, IncludeBackgroundOption.INCLUDE, IncludeBackgroundOption.EXCLUDE, IncludeBackgroundOption.ONLY);
		}

		private MapType(String name, Color color, IncludeBackgroundOption... bgOps) {
			this.name = name;
			this.color = color;
			this.bgOps = bgOps;
		}
	}

	public static void main(String[] args) throws IOException {
		EnumMap<MapType, File> dirs = new EnumMap<>(MapType.class);
		
		double[] periods = {0d, 0.2d, 1d, 5d};
		ReturnPeriods[] rps = SolHazardMapCalc.MAP_RPS;
		
		String suffix = "-vs760";
//		String suffix = "-vs260";
		
		dirs.put(MapType.COMBINED, new File(INV_DIR, COMBINED_DIR.getName()+"-ba_only"+suffix));
		dirs.put(MapType.CRUSTAL, new File(INV_DIR, CRUSTAL_DIR.getName()+"-ba_only"+suffix));
		dirs.put(MapType.SUBDUCTION, new File(INV_DIR, SUBDUCTION_DIR.getName()+"-ba_only-both_fms"+suffix));
		dirs.put(MapType.SUBDUCTION_INTERFACE, new File(INV_DIR, SUBDUCTION_DIR.getName()+"-ba_only-INTERFACE_only"+suffix));
		dirs.put(MapType.SUBDUCTION_SLAB, new File(INV_DIR, SUBDUCTION_DIR.getName()+"-ba_only-SLAB_only"+suffix));
		
		File texFile = new File(new File(FIGURES_DIR, "hazard_map_disagg"), "disagg_stats.tex");
		FileWriter texFW = null;
		if (suffix.contains("760"))
			texFW = new FileWriter(texFile);
		double texPeriod = 1d;
		ReturnPeriods texRP = ReturnPeriods.TWO_IN_50;
		
		double debugPeriod = 5d;
		ReturnPeriods debugRP = ReturnPeriods.TWO_IN_50;
//		Location debugLoc = new Location(17.75, -64.7);
		Location debugLoc = null;
		
		double cptMax = 100; // percent
		double discreteDelta = 5d;
		
		MapType refType = MapType.COMBINED;
		IncludeBackgroundOption refBGType = IncludeBackgroundOption.INCLUDE;
		File refDir = dirs.get(refType);
		
		GriddedRegion gridReg = GriddedRegion.fromFeature(Feature.read(new File(refDir, "gridded_region.json")));
		int debugIndex = debugLoc == null ? -1 : gridReg.indexForLocation(debugLoc);
		
		double annX = gridReg.getMinLon() + 0.95*(gridReg.getMaxLon() - gridReg.getMinLon());
		double annY1 = gridReg.getMinLat() + 0.95*(gridReg.getMaxLat() - gridReg.getMinLat());
		double annY2 = gridReg.getMinLat() + 0.85*(gridReg.getMaxLat() - gridReg.getMinLat());
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 28);
		
		File reportDir = new File(refDir, "hazard_map_disagg");
		Preconditions.checkState(reportDir.exists() || reportDir.mkdir());
		File resourcesDir = new File(reportDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		XY_DataSet[] polBounds = PoliticalBoundariesData.loadDefaultOutlines(gridReg);
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg, polBounds);
		mapMaker.setWritePDFs(true);
		mapMaker.setWriteGeoJSON(false);
		
		GriddedGeoDataSet mask = null;
		mask = new GriddedGeoDataSet(gridReg);
		for (XY_DataSet polBound : polBounds) {
			LocationList list = new LocationList();
			for (Point2D pt : polBound)
				list.add(new Location(pt.getY(), pt.getX()));
			Region reg = new Region(list, BorderType.MERCATOR_LINEAR);
			double halfLatSpacing = gridReg.getLatSpacing()*0.5;
			double halfLonSpacing = gridReg.getLonSpacing()*0.5;
			for (int i=0; i<gridReg.getNodeCount(); i++) {
				if (mask.get(i) == 1d)
					continue;
				Location center = gridReg.getLocation(i);
				Location[] testLocs = {
						center,
						new Location(center.lat+halfLatSpacing, center.lon+halfLonSpacing),
						new Location(center.lat+halfLatSpacing, center.lon-halfLonSpacing),
						new Location(center.lat-halfLatSpacing, center.lon+halfLonSpacing),
						new Location(center.lat-halfLatSpacing, center.lon-halfLonSpacing),
				};
				for (Location loc : testLocs) {
					if (reg.contains(loc)) {
						mask.set(i, 1d);
						break;
					}
				}
			}
		}
		System.out.println("Mask kept "+(float)(100d*mask.getSumZ()/gridReg.getNodeCount())+" % of locations");
		
		lines.add("# PRVI Map Disaggregations");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		GeoDataSet[][] refMaps = new GeoDataSet[periods.length][rps.length];
		ZipFile refZip = new ZipFile(new File(refDir, "results_hazard_INCLUDE.zip"));
		System.out.println("Loading reference hazard maps");
		for (int p=0; p<periods.length; p++) {
			for (int r=0; r<rps.length; r++) {
				String mapFileName = MPJ_LogicTreeHazardCalc.mapPrefix(periods[p], rps[r])+".txt";
				ZipEntry entry = refZip.getEntry(mapFileName);
				Preconditions.checkNotNull(entry, "%s doesn't exist in %s", mapFileName, refZip.getName());
				GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
				BufferedReader bRead = new BufferedReader(new InputStreamReader(refZip.getInputStream(entry)));
				String line = bRead.readLine();
				int index = 0;
				while (line != null) {
					line = line.trim();
					if (!line.startsWith("#")) {
						StringTokenizer tok = new StringTokenizer(line);
						double lon = Double.parseDouble(tok.nextToken());
						double lat = Double.parseDouble(tok.nextToken());
						double val = Double.parseDouble(tok.nextToken());
						Location loc = new Location(lat, lon);
						Preconditions.checkState(LocationUtils.areSimilar(loc, gridReg.getLocation(index)));
						xyz.set(index++, val);
					}
					line = bRead.readLine();
				}
				Preconditions.checkState(index == gridReg.getNodeCount());
				refMaps[p][r] = xyz;
			}
		}
		refZip.close();
		
		System.out.println("Loading exceedance rate maps");
		Table<MapType, IncludeBackgroundOption, GeoDataSet[][]> exceedRatesTable = HashBasedTable.create();
		GeoDataSet[][] refExceedRates = null;
		double maxFract = 0d;
		for (MapType type : MapType.values()) {
			File dir = dirs.get(type);
			File resultsDir = new File(dir, "results");
			Preconditions.checkState(resultsDir.exists());
			for (IncludeBackgroundOption bgOp : type.bgOps) {
				File hazardDir = new File(resultsDir, "hazard_"+(float)gridReg.getSpacing()+"deg_grid_seis_"+bgOp.name());
				Preconditions.checkState(hazardDir.exists());
				System.out.println("Loading "+type+", "+bgOp);
				GeoDataSet[][] exceedRates = new GeoDataSet[periods.length][rps.length];
				exceedRatesTable.put(type, bgOp, exceedRates);
				double myMaxFract = 0d;
				for (int p=0; p<periods.length; p++) {
					for (int r=0; r<rps.length; r++) {
						String curvesFileName = SolHazardMapCalc.getCSV_FileName("curves", periods[p]);
						File curvesFile = new File(hazardDir, curvesFileName);
						if (!curvesFile.exists())
							curvesFile = new File(hazardDir, curvesFileName+".gz");
						CSVFile<String> csv = CSVFile.readFile(curvesFile, true);
						GriddedGeoDataSet map = new GriddedGeoDataSet(gridReg);
						DiscretizedFunc[] curves = SolHazardMapCalc.loadCurvesCSV(csv, gridReg, true);
						for (int i=0; i<curves.length; i++) {
							boolean debug = debugIndex == i && (float)debugPeriod == (float)periods[p] && rps[r] == debugRP;
							if (curves[i] == null) {
								map.set(i, 0d);
							} else {
								double iml = refMaps[p][r].get(i);
								double prob = curves[i].getInterpolatedY_inLogXLogYDomain(iml);
								double rate = -Math.log(1 - prob);
								map.set(i, rate);
								if (refExceedRates != null) {
									double fract = rate/refExceedRates[p][r].get(i);
									if (debug) {
										System.out.println("DEBUG p="+(float)periods[p]+", "+rps[r]);
										System.out.println("\t"+(float)rate+"/"+(float)refExceedRates[p][r].get(i)+" = "+(float)fract);
									}
									myMaxFract = Math.max(myMaxFract, fract);
								}
								if (debug) {
									System.out.println("CURVE DATA");
									System.out.println(curves[i]);
//									for (int j=0; j<curves[i].size(); j++)
//										System.out.println((float)curves[i].getX(j)+", "+(float)curves[i].getY(j));
									System.out.println();
								}
							}
						}
						exceedRates[p][r] = map;
						
					}
				}
				if (bgOp == refBGType && type == refType) {
					refExceedRates = exceedRates;
				} else {
					System.out.println("Max Fract: "+(float)myMaxFract);
					maxFract = Math.max(myMaxFract, maxFract);
				}
			}
		}
		
		System.out.println("Done loading rates");
		
		if (debugIndex >= 0) {
			for (int p=0; p<periods.length; p++) {
				if (periods[p] == debugPeriod) {
					for (int r=0; r<rps.length; r++) {
						if (rps[r] == debugRP) {
							System.out.println("Debug contrib sums");
							double refExceed = refExceedRates[p][r].get(debugIndex);
							double faultExeed = exceedRatesTable.get(MapType.COMBINED, IncludeBackgroundOption.EXCLUDE)[p][r].get(debugIndex)/refExceed;
							double gridExeed = exceedRatesTable.get(MapType.COMBINED, IncludeBackgroundOption.ONLY)[p][r].get(debugIndex)/refExceed;
							System.out.println("Combined fault + gridded = "+(float)faultExeed+" + "+(float)gridExeed
									+" = "+(float)(faultExeed + gridExeed));
							// calculate fault by itself
							double crustalFaultExceed = exceedRatesTable.get(MapType.CRUSTAL, IncludeBackgroundOption.EXCLUDE)[p][r].get(debugIndex)/refExceed;
							double interfaceFaultExceed = exceedRatesTable.get(MapType.SUBDUCTION_INTERFACE, IncludeBackgroundOption.EXCLUDE)[p][r].get(debugIndex)/refExceed;
							System.out.println("crustal fault + interface fault = "+(float)crustalFaultExceed+" + "+(float)interfaceFaultExceed
									+" = "+(float)(crustalFaultExceed + interfaceFaultExceed));
							// calculate gridded by itself
							double crustalGridExceed = exceedRatesTable.get(MapType.CRUSTAL, IncludeBackgroundOption.ONLY)[p][r].get(debugIndex)/refExceed;
							double interfaceGridExceed = exceedRatesTable.get(MapType.SUBDUCTION_INTERFACE, IncludeBackgroundOption.ONLY)[p][r].get(debugIndex)/refExceed;
							double slabGridExceed = exceedRatesTable.get(MapType.SUBDUCTION_SLAB, IncludeBackgroundOption.ONLY)[p][r].get(debugIndex)/refExceed;
							System.out.println("crustal grid + interface grid + slab grid = "+(float)crustalGridExceed+" + "+(float)interfaceGridExceed+" + "+(float)slabGridExceed
									+" = "+(float)(crustalGridExceed + interfaceGridExceed + slabGridExceed));
						}
					}
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		for (int p=0; p<periods.length; p++) {
			String perLabel, perPrefix;
			if (periods[p] == 0d) {
				perLabel = "PGA";
				perPrefix = "pga";
			} else {
				perLabel = oDF.format(periods[p])+"s SA";
				perPrefix = oDF.format(periods[p])+"s";
			}
			for (int r=0; r<rps.length; r++) {
				String rpLabel = rps[r].label;
				String rpPrefix = rps[r].name();
				
				String perRPLabel = perLabel+", "+rpLabel;
				String perRPPrefix = perPrefix+"_"+rpPrefix;
				
				lines.add("## "+perRPLabel);
				lines.add(topLink); lines.add("");
				
				for (MapType compType : MapType.values()) {
					IncludeBackgroundOption[] compBGTypes;
					if (compType == refType) {
						compBGTypes = new IncludeBackgroundOption[] { IncludeBackgroundOption.EXCLUDE, IncludeBackgroundOption.ONLY };
					} else {
						compBGTypes = compType.bgOps;
					}
					
					System.out.println("Plotting "+compType+", "+perRPLabel);
					
					
					System.out.println("--- COLORS ----");
					Color maxColor = compType.color;
					System.out.println("Original color: "+maxColor);
					float[] whiteLab = Colors.toLab(Color.WHITE);
					System.out.println("WHITE LAB: "+whiteLab[0]+" "+whiteLab[1]+" "+whiteLab[2]);
					float[] lab = Colors.toLab(maxColor);
					// reset to the same perceptual luminence
					System.out.println("ORIG LAB: "+lab[0]+" "+lab[1]+" "+lab[2]);
					lab[0] = 20;
					System.out.println("MOD LAB: "+lab[0]+" "+lab[1]+" "+lab[2]);
					// try to convert to 30 luminence
					float[] rgb = Colors.LabToRGB(lab[0], lab[1], lab[2]);
					System.out.println("RGB: "+rgb[0]+" "+rgb[1]+" "+rgb[2]);
					maxColor = new Color(
							Float.max(0f, Float.min(1f, rgb[0])),
							Float.max(0f, Float.min(1f, rgb[1])),
							Float.max(0f, Float.min(1f, rgb[2])));
					float[] backToLab = Colors.toLab(maxColor);
					System.out.println("Back to LAB: "+backToLab[0]+" "+backToLab[1]+" "+backToLab[2]);
					CPT cpt = new CPT(0, cptMax, Color.WHITE, maxColor);
					cpt.setNanColor(Color.WHITE);
					cpt.setAboveMaxColor(cpt.getMaxColor());
					if (discreteDelta > 0d)
						cpt = cpt.asDiscrete(discreteDelta, false);
					System.out.println("------------------");
					
					String includeLabel, excludeLabel, onlyLabel, texPrefix;
					if (compType == MapType.COMBINED) {
						includeLabel = null;
						excludeLabel = "All Fault Sources";
						onlyLabel = "All Gridded Sources";
						texPrefix = "DisaggCombined";
					} else if (compType == MapType.CRUSTAL) {
						includeLabel = "All Crustal Sources";
						excludeLabel = "Crustal Fault Sources";
						onlyLabel = "Crustal Gridded Sources";
						texPrefix = "DisaggCrustal";
					} else if (compType == MapType.SUBDUCTION) {
						includeLabel = "All Subduction Sources";
						excludeLabel = "Subduction Interface Fault Sources";
						onlyLabel = "Subduction Gridded Sources (Interface & Slab)";
						texPrefix = "DisaggSubduction";
					} else if (compType == MapType.SUBDUCTION_INTERFACE) {
						includeLabel = "All Subduction Interface Sources";
						excludeLabel = "Subduction Interface Fault Sources";
						onlyLabel = "Subduction Interface Gridded Sources";
						texPrefix = "DisaggInterface";
					} else if (compType == MapType.SUBDUCTION_SLAB) {
						includeLabel = null;
						excludeLabel = null;
						onlyLabel = "Subduction Intraslab Sources";
						texPrefix = "DisaggSlab";
					} else {
						throw new IllegalStateException();
					}
					List<String> labels = new ArrayList<>();
					
					for (IncludeBackgroundOption bgOp : compBGTypes) {
						String label;
						switch (bgOp) {
						case INCLUDE:
							label = includeLabel;
							break;
						case EXCLUDE:
							label = excludeLabel;
							break;
						case ONLY:
							label = onlyLabel;
							break;

						default:
							throw new IllegalStateException();
						}
						labels.add(label);
					}
					
					lines.add("### "+compType.name+", "+perRPLabel);
					lines.add(topLink); lines.add("");
					
					TableBuilder table = MarkdownUtils.tableBuilder();
					table.addLine(labels);
					
					table.initNewLine();
					for (int i=0; i<compBGTypes.length; i++) {
						String label = labels.get(i);
						String prefix = compType.name()+"_"+compBGTypes[i].name()+"_"+perRPPrefix;
						
						GeoDataSet rateMap = exceedRatesTable.get(compType, compBGTypes[i])[p][r];
						GeoDataSet refMap = refExceedRates[p][r];
						GeoDataSet contribMap = new GriddedGeoDataSet(gridReg);
						int numWith = 0;
						double min = Double.POSITIVE_INFINITY;
						double max = 0d;
						double sum = 0d;
						for (int n=0; n<contribMap.size(); n++) {
							if (mask != null && mask.get(n) != 1d) {
								contribMap.set(n, Double.NaN);
							} else {
								double percent = 100d*rateMap.get(n)/refMap.get(n);
								min = Math.min(min, percent);
								max = Math.max(max, percent);
								sum += percent;
								contribMap.set(n, percent);
								numWith++;
							}
						}
						
						mapMaker.clearAnnotations();
						XYTextAnnotation rangeAnn = new XYTextAnnotation("Range: ["+(int)(min+0.5)+"%, "+(int)(max+0.5)+"%]", annX, annY1);
						rangeAnn.setFont(annFont);
						rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
						mapMaker.addAnnotation(rangeAnn);
						XYTextAnnotation avgAnn = new XYTextAnnotation("Average: "+(int)((sum/numWith)+0.5)+"%", annX, annY2);
						avgAnn.setFont(annFont);
						avgAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
						mapMaker.addAnnotation(avgAnn);
						
						mapMaker.plotXYZData(contribMap, cpt, "% Contribution, "+perRPLabel);
						mapMaker.plot(resourcesDir, prefix, label);
						
						table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+".png)");
						
						if (periods[p] == texPeriod && rps[r] == texRP && texFW != null) {
							String myTexPrefix = texPrefix;
							if (compBGTypes[i] == IncludeBackgroundOption.EXCLUDE)
								myTexPrefix += "Fault";
							else if (compBGTypes[i] == IncludeBackgroundOption.ONLY)
								myTexPrefix += "Gridded";
							System.out.println("Writing tex for T="+(float)periods[p]+"s, "+rps[r]+" with prefix: "+myTexPrefix);
							texFW.write(LaTeXUtils.defineValueCommand(myTexPrefix+"Min",
									LaTeXUtils.numberAsPercent(min, 0))+"\n");
							texFW.write(LaTeXUtils.defineValueCommand(myTexPrefix+"Max",
									LaTeXUtils.numberAsPercent(max, 0))+"\n");
							texFW.write(LaTeXUtils.defineValueCommand(myTexPrefix+"Avg",
									LaTeXUtils.numberAsPercent(sum/numWith, 0))+"\n");
						}
					}
					table.finalizeLine();
					
					lines.addAll(table.build());
				}
			}
		}
		
		System.out.println("DONE");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, reportDir);
		
		if (texFW != null)
			texFW.close();
	}

}
