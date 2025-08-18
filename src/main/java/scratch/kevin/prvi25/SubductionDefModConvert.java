package scratch.kevin.prvi25;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.GeoJSON_Type;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.faultSysSolution.util.minisections.MinisectionMappings;
import org.opensha.sha.earthquake.faultSysSolution.util.minisections.MinisectionSlipRecord;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCouplingModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;
import com.google.gson.Gson;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.figures.PRVI_Paths;
import scratch.kevin.prvi25.figures.PRVI_SubductionSubSectPlots;

public class SubductionDefModConvert {
	
	private static boolean WRITE_SECT_INDEXES = false;

	public static void main(String[] args) throws IOException {
		String version = "v5";
		File fmDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/fault_models/subduction/"+version);
		File dmDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/def_models/subduction/"+version);
		File inputsDir = new File(fmDir, "inputs");
		
//		// disable writing outputs (only plots)
//		fmDir = null;
//		dmDir = null;
		
//		File inputFile = new File(fmDir, "PRVI_sub_v2_13May2024.geojson");
//		File fmOutputFile = new File(fmDir, "PRVI_sub_v2_fault_model.geojson");
//		File fullOutputFile = new File(dmDir, "PRVI_sub_v2_full_minisections.txt");
//		File partOutputFile = new File(dmDir, "PRVI_sub_v2_partial_minisections.txt");
		
		PRVI25_SubductionFaultModels[] fms = {
				PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE, PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL};
		PRVI25_SubductionDeformationModels[] dms = {
				PRVI25_SubductionDeformationModels.FULL, PRVI25_SubductionDeformationModels.PARTIAL};
		PRVI25_SubductionCouplingModels[] couplings = {
				PRVI25_SubductionCouplingModels.LOW, PRVI25_SubductionCouplingModels.PREFERRED,
				PRVI25_SubductionCouplingModels.HIGH};
		
		boolean interpolate = true;
		boolean interpSymmetry = true;
		boolean interpAcrossNH_PR = true;
		double fractUncert = 0.1;
		
//		File plotOutputDir = new File("/home/kevin/Documents/papers/2024_PRVI_Subduction/figures/def_model");
		File plotOutputDir = new File(PRVI_Paths.FIGURES_DIR, "sub_dm");
//		File plotOutputDir = new File("/tmp/dm_plots");
		Preconditions.checkState(plotOutputDir.exists() || plotOutputDir.mkdir());
		
		List<PRVI25_SubductionFaultModels> debugPlotFMs = new ArrayList<>();
		List<PRVI25_SubductionDeformationModels> debugPlotDMs = new ArrayList<>();
		List<PRVI25_SubductionCouplingModels> debugPlotCouplings = new ArrayList<>();
		List<String> debugPlotPrefixes = new ArrayList<>();
		
		for (PRVI25_SubductionFaultModels fm : fms) {
			String inGeomPrefix = fm == PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE ? "LargePolys" : "SmallPolys";
			String outGeomPrefix = fm.getFilePrefix();
			String geomTitle = fm.getShortName()+" Fault Model";
			for (PRVI25_SubductionDeformationModels dm : dms) {
				String inDMPrefix = dm == PRVI25_SubductionDeformationModels.FULL ? "FullRate" : "PartRate";
				String outDMPrefix = dm.getFilePrefix();
				
				File fmOutputFile;
				if (dm == PRVI25_SubductionDeformationModels.FULL && fmDir != null)
					fmOutputFile = new File(fmDir, "PRVI_sub_"+version+"_fault_model_"+outGeomPrefix+".geojson");
				else
//					fmOutputFile = new File(fmDir, "PRVI_sub_"+version+"_fault_model_"+outGeomPrefix+"_dup.geojson");
					fmOutputFile = null;
				
				File inputFile = new File(inputsDir, "PRVI_subduction_"+inDMPrefix+"_"+inGeomPrefix+"_LowPrefHigh.geojson");
				System.out.println("Processing "+inputFile.getName());
				
				String ratePropPrefix;
				String rakePropName;
				if (dm == PRVI25_SubductionDeformationModels.FULL) {
					ratePropPrefix = "RateFull";
					rakePropName = "RakeRateF";
				} else {
					ratePropPrefix = "RatePart";
					rakePropName = "RakeRateP";
				}
				
				for (PRVI25_SubductionCouplingModels coupling : couplings) {
					String outCouplingPrefix = coupling.getFilePrefix();
					String ratePropName;
					switch (coupling) {
					case HIGH:
						ratePropName = ratePropPrefix+"High";
						break;
					case PREFERRED:
						ratePropName = ratePropPrefix+"Pref";
						break;
					case LOW:
						ratePropName = ratePropPrefix+"Low";
						break;

					default:
						throw new IllegalStateException("Unrecognized coupling: "+coupling);
					}
					
					String rateTitle = dm.getName()+", "+coupling.getShortName()+" Coupling";
					
					File dmOutputFile;
					if (dmDir == null)
						// plot only
						dmOutputFile = null;
					else
						dmOutputFile = new File(dmDir, "PRVI_sub_"+version+"_"
								+outGeomPrefix+"_"+outDMPrefix+"_"+outCouplingPrefix+"_rate_minisections.txt");
					
					try {
						Gson inGSON = FeatureCollection.buildGson(DepthSerializationType.DEPTH_M);
						BufferedReader read = new BufferedReader(new FileReader(inputFile));
						List<Feature> inFeatures = inGSON.fromJson(read, FeatureCollection.class).features;
						read.close();
						
						String prefix = "geom_"+inGeomPrefix+"_"+inDMPrefix+"_"+coupling.getShortName();
						String fmTitle = geomTitle;
						String fmDmTitle = geomTitle+", "+rateTitle;
						debugPlotFMs.add(fm);
						debugPlotDMs.add(dm);
						debugPlotCouplings.add(coupling);
						debugPlotPrefixes.add(prefix);
						debugPlotSects(inFeatures, plotOutputDir, prefix, fmTitle, fmDmTitle,
								ratePropName, rakePropName);
						
						Map<String, int[]> stitchInIDs = new HashMap<>();
						Map<String, Integer> stitchOutIDs = new HashMap<>();
						
						String name0 = "Northern Hispaniola Subduction";
						int[] inIDs0 = {7501};
						int outID0 = 7501;
						stitchInIDs.put(name0, inIDs0);
						stitchOutIDs.put(name0, outID0);
						
						String name1 = "PRVI - Lesser Antilles Subduction";
						int[] inIDs1 = {7502, 7503, 7504, 7505, 7506, 7507, 7508};
						int outID1 = 7500;
						stitchInIDs.put(name1, inIDs1);
						stitchOutIDs.put(name1, outID1);
						
						String name2 = "Muertos Thrust";
						int[] inIDs2 = {7553, 7552, 7551};
						int outID2 = 7550;
						stitchInIDs.put(name2, inIDs2);
						stitchOutIDs.put(name2, outID2);
						
						Map<String, int[]> interpInIDs = null;
						
						if (interpolate) {
							interpInIDs = new HashMap<>();
							interpInIDs.put("Hispaniola - PRVI - Lesser Antilles", Ints.concat(inIDs0, inIDs1));
							interpInIDs.put("Muertos", inIDs2);
						}
						
						Map<String, Feature[]> stitchFeatureSets = new HashMap<>();
						for (String name : stitchInIDs.keySet())
							stitchFeatureSets.put(name, new Feature[stitchInIDs.get(name).length]);
						
						for (Feature feature : inFeatures) {
							int id = feature.properties.getInt("id", -1);
							Preconditions.checkState(id >= 0);
							boolean found = false;
							for (String name : stitchInIDs.keySet()) {
								int[] stitchIDs = stitchInIDs.get(name);
								for (int i=0; i<stitchIDs.length; i++) {
									if (id == stitchIDs[i]) {
										Preconditions.checkState(!found, "%s is in multiple sets! 2nd encountered is %s", id, name);
										stitchFeatureSets.get(name)[i] = feature;
										found = true;
									}
								}
							}
							Preconditions.checkState(found);
						}
						
						List<Feature> outFeatures = new ArrayList<>();
						
						Map<Integer, List<MinisectionSlipRecord>> minisectionRecsMap = new HashMap<>();
						Map<Integer, List<MinisectionSlipRecord>> minisectionRecsMapToOrigIDs = new HashMap<>();
						
						Map<Integer, String> shorNamesMap = new HashMap<>();
						
						for (String name : stitchInIDs.keySet()) {
							int outID = stitchOutIDs.get(name);
							Feature[] features = stitchFeatureSets.get(name);
							
							FaultTrace stitchedUpperTrace = new FaultTrace(name);
							FaultTrace stitchedLowerTrace = new FaultTrace(name);
							
							System.out.println("Building "+name);
							
							List<MinisectionSlipRecord> minisectionRecs = new ArrayList<>();
							
							int prevID = -1;
							
//							List<Double> minisectionDASs = new ArrayList<>();
//							ArbitrarilyDiscretizedFunc dasSlipFunc = new ArbitrarilyDiscretizedFunc();
//							ArbitrarilyDiscretizedFunc dasSlipUncertFunc = new ArbitrarilyDiscretizedFunc();
//							ArbitrarilyDiscretizedFunc dasRakeFunc = new ArbitrarilyDiscretizedFunc();
//							double runningDAS = 0d;
							
							for (Feature feature : features) {
								Preconditions.checkState(feature.geometry.type == GeoJSON_Type.MultiLineString);
								MultiLineString geometry = (MultiLineString)feature.geometry;
								Preconditions.checkState(geometry.lines.size() == 2);
								FaultTrace upper = new FaultTrace(null);
								upper.addAll(geometry.lines.get(0));
								FaultTrace lower = new FaultTrace(null);
								lower.addAll(geometry.lines.get(1));
								
								int sectID = feature.properties.getInt("id", -1);
								String shortName = feature.properties.getString("FaultName");
								String sectName = shortName+" ("+feature.properties.getString("FaultDesc")+")";
								System.out.println("Processing "+sectID+". "+sectName);
								
								shorNamesMap.put(sectID, shortName);
								
								// check that lower is in the same direction
								double upperStrike = upper.getAveStrike();
								double lowerSrike = lower.getAveStrike();
								System.out.println("\tUpper trace: depth="+depthStr(upper)+"; strike="+oneDigit.format(upperStrike));
								System.out.println("\tLower trace: depth="+depthStr(lower)+"; strike="+oneDigit.format(lowerSrike));
								double diff = FaultUtils.getAbsAngleDiff(upperStrike, lowerSrike);
								if (diff > 90d) {
									System.out.println("\tREVERSING LOWER (diff="+oneDigit.format(diff)+")");
									lower.reverse();
									geometry = new MultiLineString(List.of(upper, lower));
								}
								
								// now check Aki & Richards
								Location avgUpper = avgTraceLoc(upper);
								Location avgLower = avgTraceLoc(lower);
								double azUpperToLower = LocationUtils.azimuth(avgUpper, avgLower);
								double perfectRight = upperStrike + 90d;
								while (perfectRight > 360d)
									perfectRight -= 360d;
								double arDiff = FaultUtils.getAbsAngleDiff(azUpperToLower, perfectRight);
//							System.out.println("\tARdiff="+oneDigit.format(arDiff));
								if (arDiff > 90d) {
									System.out.println("REVERSING BOTH due to Aki & Richards violation; upper + 90 is "
											+oneDigit.format(perfectRight)+", upper to lower is "
											+oneDigit.format(azUpperToLower)+", diff is "+oneDigit.format(arDiff));
									upper.reverse();
									lower.reverse();
									geometry = new MultiLineString(List.of(upper, lower));
								}
								
								if (!stitchedUpperTrace.isEmpty()) {
									// make sure we're congiguous
									double upperDist = LocationUtils.linearDistance(stitchedUpperTrace.last(), upper.first());
									Preconditions.checkState(upperDist < 0.01d, "Upper trace not contiguous: %s, %s, dist=%s, %s. %s; prev=%s",
											stitchedUpperTrace.last(), upper.first(), upperDist, sectID, sectName, prevID);
									double lowerDist = LocationUtils.linearDistance(stitchedLowerTrace.last(), lower.first());
									Preconditions.checkState(lowerDist < 0.01d, "Lower trace not contiguous: %s, %s, dist=%s, %s. %s; prev=%s",
											stitchedLowerTrace.last(), lower.first(), lowerDist, sectID, sectName, prevID);
									
									// add all but the duplicated first point
									for (int i=1; i<upper.size(); i++)
										stitchedUpperTrace.add(upper.get(i));
									for (int i=1; i<lower.size(); i++)
										stitchedLowerTrace.add(lower.get(i));
								} else {
									// add all
									stitchedUpperTrace.addAll(upper);
									stitchedLowerTrace.addAll(lower);
								}
								
								FeatureProperties props = feature.properties;
								double slip = props.getDouble(ratePropName, Double.NaN);
								double slipUncert = slip*fractUncert;
								double rake = props.getDouble(rakePropName, Double.NaN);
								
								List<MinisectionSlipRecord> myMinis = new ArrayList<>(upper.size()-1);
								
								for (int i=1; i<upper.size(); i++) {
									Location startLoc = upper.get(i-1);
									Location endLoc = upper.get(i);
									int minisectionID = minisectionRecs.size();
									MinisectionSlipRecord mini = new MinisectionSlipRecord(outID, minisectionID, startLoc, endLoc, rake, slip, slipUncert);
									minisectionRecs.add(mini);
									myMinis.add(mini);
								}
								
								minisectionRecsMapToOrigIDs.put(sectID, myMinis);
								
								prevID = sectID;
							}
							
							// build stitched feature
							MultiLineString geometry = new MultiLineString(List.of(stitchedUpperTrace, stitchedLowerTrace));
							FeatureProperties props = new FeatureProperties();
							props.set(GeoJSONFaultSection.FAULT_ID, outID);
							props.set(GeoJSONFaultSection.FAULT_NAME, name);
							props.set("PrimState", "PR");
							outFeatures.add(new Feature(outID, geometry, props));
							
							minisectionRecsMap.put(outID, minisectionRecs);
						}
						
						if (interpolate) {
							System.out.println("Interpolating");
							
							// build original subsections and calculate original moment rates
							List<FaultSection> fullSects = new ArrayList<>();
							for (Feature feature : outFeatures)
								fullSects.add(GeoJSONFaultSection.fromFeature(feature));
							List<FaultSection> subSects = SubSectionBuilder.buildSubSects(fullSects, 2, 0.5, 30);
							MinisectionMappings mappings = new MinisectionMappings(fullSects, subSects);
							mappings.mapDefModelMinisToSubSects(minisectionRecsMap);
							Map<String, Double> origMoRates = new HashMap<>();
							for (FaultSection sect : subSects) {
								String name = sect.getParentSectionName();
								Double moRate = origMoRates.get(name);
								if (moRate == null)
									moRate = 0d;
								origMoRates.put(name, moRate+sect.calcMomentRate(false));
							}
							
							for (String name : interpInIDs.keySet()) {
								int[] ids = interpInIDs.get(name);
								System.out.println("Interpolating "+name);
								
								// figure out section bounds in terms of DAS
								double runningDAS = 0d;
								List<MinisectionSlipRecord> interpMinis = new ArrayList<>();
								List<Double> miniMidDASs = new ArrayList<>();
								
								List<Double> sectStartDASs = new ArrayList<>();
								List<Double> sectMidDASs = new ArrayList<>();
								List<String> sectShortNames = new ArrayList<>();
								List<Double> sectEndDASs = new ArrayList<>();
								List<MinisectionSlipRecord> sectRefMinis = new ArrayList<>();
								
								for (int id : ids) {
									List<MinisectionSlipRecord> minis = minisectionRecsMapToOrigIDs.get(id);
									Preconditions.checkNotNull(minis);
									
									double start = runningDAS;
									
									for (MinisectionSlipRecord mini : minis) {
										double len = LocationUtils.linearDistanceFast(mini.startLoc, mini.endLoc);
										runningDAS += len;
										
										interpMinis.add(mini);
										miniMidDASs.add(runningDAS - 0.5*len);
									}
									double end = runningDAS;
									double middle = 0.5*(start + end);
									runningDAS = end;
									
									sectRefMinis.add(minis.get(0));
									sectStartDASs.add(start);
									sectMidDASs.add(middle);
									sectEndDASs.add(end);
									sectShortNames.add(shorNamesMap.get(id));
								}
								
								System.out.println("Total DAS for "+name+": "+runningDAS);
								
								// build slip/rake as a function of DAS
								ArbitrarilyDiscretizedFunc dasSlipFunc = new ArbitrarilyDiscretizedFunc();
								ArbitrarilyDiscretizedFunc dasSlipUncertFunc = new ArbitrarilyDiscretizedFunc();
								ArbitrarilyDiscretizedFunc dasRakeFunc = new ArbitrarilyDiscretizedFunc();
								
								List<double[]> interpRanges = new ArrayList<>();
								
								for (int i=0; i<sectRefMinis.size(); i++) {
									MinisectionSlipRecord mini = sectRefMinis.get(i);
									double start = sectStartDASs.get(i);
									double middle = sectMidDASs.get(i);
									double end = sectEndDASs.get(i);
									
//									System.out.println(i+". "+mini.parentID);
									
									if (i == 0) {
										// first, start at the beginning not the middle
										Preconditions.checkState(start == 0d);
										dasSlipFunc.set(0d, mini.slipRate);
										dasSlipUncertFunc.set(0d, mini.slipRateStdDev);
										dasRakeFunc.set(0d, mini.rake);
									} else if (!interpAcrossNH_PR && i > 0 && mini.parentID == 7500 && sectRefMinis.get(i-1).parentID == 7501) {
										// skip interpolation between NH and PR
										interpRanges.add(null);
										
										MinisectionSlipRecord prev = sectRefMinis.get(i-1);
										
										double prevEnd = start - 1e-6;
										dasSlipFunc.set(prevEnd, prev.slipRate);
										dasSlipUncertFunc.set(prevEnd, prev.slipRateStdDev);
										dasRakeFunc.set(prevEnd, prev.rake);
										dasSlipFunc.set(start, mini.slipRate);
										dasSlipUncertFunc.set(start, mini.slipRateStdDev);
										dasRakeFunc.set(start, mini.rake);
									} else if (interpSymmetry) {
										// enforce symmetry
										double prevHalfLen = sectEndDASs.get(i-1) - sectMidDASs.get(i-1);
										double myHalfLen = middle - start;
										double interpLen = Math.min(prevHalfLen, myHalfLen);
										interpRanges.add(new double[] {start - interpLen, start + interpLen});
										if (myHalfLen > prevHalfLen) {
											// my portion is longer, truncate
											double interpX = start + prevHalfLen;
											dasSlipFunc.set(interpX, mini.slipRate);
											dasSlipUncertFunc.set(interpX, mini.slipRateStdDev);
											dasRakeFunc.set(interpX, mini.rake);
										} else {
											// previous portion was longer, truncate
											MinisectionSlipRecord prev = sectRefMinis.get(i-1);
											double interpX = start - myHalfLen;
											dasSlipFunc.set(interpX, prev.slipRate);
											dasSlipUncertFunc.set(interpX, prev.slipRateStdDev);
											dasRakeFunc.set(interpX, prev.rake);
										}
									} else {
										interpRanges.add(new double[] {sectMidDASs.get(i-1), middle});
									}
									
									// add the middle
									dasSlipFunc.set(middle, mini.slipRate);
									dasSlipUncertFunc.set(middle, mini.slipRateStdDev);
									dasRakeFunc.set(middle, mini.rake);
									
									if (i == sectRefMinis.size()-1) {
										// last, extend to the end
										dasSlipFunc.set(end, mini.slipRate);
										dasSlipUncertFunc.set(end, mini.slipRateStdDev);
										dasRakeFunc.set(end, mini.rake);
									}
								}
								
								// apply interpolation
								for (int i=0; i<interpMinis.size(); i++) {
									MinisectionSlipRecord mini = interpMinis.get(i);
									double middle = miniMidDASs.get(i);
									
									double slip = dasSlipFunc.getInterpolatedY(middle);
									double slipUncert = dasSlipUncertFunc.getInterpolatedY(middle);
									double rake = dasRakeFunc.getInterpolatedY(middle);
									
									MinisectionSlipRecord modMini = new MinisectionSlipRecord(
											mini.parentID, mini.minisectionID, mini.startLoc, mini.endLoc,
											rake, slip, slipUncert);
									
									// find the match
									List<MinisectionSlipRecord> origRecs = minisectionRecsMap.get(mini.parentID);
									Preconditions.checkState(mini.minisectionID < origRecs.size());
									MinisectionSlipRecord orig = origRecs.get(mini.minisectionID);
									Preconditions.checkState(orig.equals(mini));
									origRecs.set(mini.minisectionID, modMini);
								}
								
								// plot them
								String interpPrefix = prefix+"_interp_"+name.replaceAll("\\W+", "_");
								debugWriteInterpFuncsOrders(plotOutputDir, interpPrefix, name,
										dasSlipFunc, dasSlipUncertFunc, dasRakeFunc, sectStartDASs,
										sectMidDASs, interpRanges, sectShortNames);
							}
							
							// calculate modified moment rates
							mappings.mapDefModelMinisToSubSects(minisectionRecsMap);
							Map<String, Double> modMoRates = new LinkedHashMap<>();
							for (FaultSection sect : subSects) {
								String name = sect.getParentSectionName();
								Double moRate = modMoRates.get(name);
								if (moRate == null)
									moRate = 0d;
								modMoRates.put(name, moRate+sect.calcMomentRate(false));
							}
							
							System.out.println("Interpolation Mo-Rate changes for "+fmDmTitle);
							for (String name : modMoRates.keySet()) {
								double orig = origMoRates.get(name);
								double mod = modMoRates.get(name);
								System.out.println("\t"+name+":\t"+moDF.format(orig)+" -> "+moDF.format(mod)+" ("+pDF.format((mod-orig)/orig)+")");
							}
						}

						if (fmOutputFile != null)
							FeatureCollection.write(new FeatureCollection(outFeatures), fmOutputFile);
						
						if (dmOutputFile != null)
							MinisectionSlipRecord.writeMinisectionsFile(dmOutputFile, minisectionRecsMap);
					} catch (Exception e) {
						System.err.println("WARNING: Exception with "+inputFile.getName()+", skipping");
						e.printStackTrace();
					}
				}
			}
		}
		
		List<String> debugPlotTypes = new ArrayList<>();
		List<String> debugPlotSuffixes = new ArrayList<>();
		
		debugPlotTypes.add("Trace Orders");
		debugPlotSuffixes.add("_order");
		
		debugPlotTypes.add("Depths");
		debugPlotSuffixes.add("_depth");
		
		debugPlotTypes.add("Slip Rates");
		debugPlotSuffixes.add("_slips");
		
		debugPlotTypes.add("Rake Angles");
		debugPlotSuffixes.add("_rakes");
		
		List<String> lines = new ArrayList<>();
		lines.add("# PRVI Subduction FM/DM Debug");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		for (int i=0; i<debugPlotTypes.size(); i++) {
			lines.add("## "+debugPlotTypes.get(i));
			lines.add(topLink); lines.add("");
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			String suffix = debugPlotSuffixes.get(i);
			for (int d=0; d<debugPlotPrefixes.size(); d++) {
				String prefix = debugPlotPrefixes.get(d);
				
				table.addLine("![Plot]("+prefix+suffix+".png)");
			}
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, plotOutputDir);
	}
	
	private static final DecimalFormat moDF = new DecimalFormat("0.00E0");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static String depthStr(FaultTrace trace) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (Location loc : trace)
			track.addValue(loc.depth);
		return "avg="+oneDigit.format(track.getAverage())+" ["+oneDigit.format(track.getMin())+", "+oneDigit.format(track.getMax())+"]";
	}
	private static DecimalFormat oneDigit = new DecimalFormat("0.0");
	
	private static Location avgTraceLoc(FaultTrace trace) {
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		MinMaxAveTracker depTrack = new MinMaxAveTracker();
		for (Location loc : trace) {
			latTrack.addValue(loc.lat);
			lonTrack.addValue(loc.lon);
			depTrack.addValue(loc.depth);
		}
		return new Location(latTrack.getAverage(), lonTrack.getAverage(), depTrack.getAverage());
	}
	
	private static void debugPlotSects(List<Feature> features, File outputDir, String prefix, String fmTitle,
			String fmDmTitle, String slipProp, String rakeProp) throws IOException {
		List<FaultSection> sects = new ArrayList<>();
		List<Location> directionLocs = new ArrayList<>();
		List<Double> directionFracts = new ArrayList<>();
		List<Double> depths = new ArrayList<>();
		List<Double> slips = new ArrayList<>();
//		List<Double> slipUncerts = new ArrayList<>();
		List<Double> rakes = new ArrayList<>();
		CPT orderCPT = new CPT(0d, 1d, Color.RED, Color.GREEN);
		
		for (Feature feature : features) {
			Feature copy = new Feature(feature.properties.getNumber("id"), feature.geometry, feature.properties);
			GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(copy);
			sects.add(sect);
			Preconditions.checkState(feature.geometry.type == GeoJSON_Type.MultiLineString);
			MultiLineString multi = (MultiLineString)feature.geometry;
			for (LocationList line : multi.lines) {
				FaultTrace trace = new FaultTrace(null);
				trace.addAll(line);
				trace = FaultUtils.resampleTrace(trace, (int)(trace.getTraceLength()+0.5));
				for (int l=0; l<trace.size(); l++) {
					directionLocs.add(trace.get(l));
					directionFracts.add((double)l/(double)(trace.size()-1));
					depths.add(trace.get(l).depth);
				}
			}
			double slip = feature.properties.getDouble(slipProp, Double.NaN);
			slips.add(slip);
			sect.setAveSlipRate(slip);
//			slipUncerts.add(feature.properties.getDouble(slipUncertProp, Double.NaN));
			double rake = feature.properties.getDouble(rakeProp, Double.NaN);
			rakes.add(rake);
			sect.setAveRake(rake);
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(PRVI_SubductionSubSectPlots.plotReg, sects);
		
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 0.3f);
		mapMaker.plotScatterScalars(directionLocs, directionFracts, orderCPT, "First ----------> Last");

		List<Location> centerLocs = new ArrayList<>();
		for (FaultSection sect : sects) {
			RuptureSurface surf = sect.getFaultSurface(1d);
			LocationList locs = surf.getEvenlyDiscritizedListOfLocsOnSurface();
			double avgLat = locs.stream().mapToDouble(L->L.lat).average().getAsDouble();
			double avgLon = locs.stream().mapToDouble(L->L.lon).average().getAsDouble();
			centerLocs.add(new Location(avgLat, avgLon));
		}
		
		
		if (WRITE_SECT_INDEXES) {
			Font font = new Font(Font.SANS_SERIF, Font.BOLD, 16);
			
			for (int i = 0; i < sects.size(); i++) {
				FaultSection sect = sects.get(i);
				Location centerLoc = centerLocs.get(i);
				XYTextAnnotation idAnn = new XYTextAnnotation(sect.getSectionId()+"", centerLoc.lon, centerLoc.lat);
				idAnn.setFont(font);
				idAnn.setTextAnchor(TextAnchor.CENTER);
				mapMaker.addAnnotation(idAnn);
			}
		}
		
		mapMaker.plot(outputDir, prefix+"_order", fmTitle);
		
		// add crustal fault lines
		List<LocationList> crustalTraces = new ArrayList<>();
		for (FaultSection sect : PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1.getFaultSections())
			crustalTraces.add(sect.getFaultTrace());
		mapMaker.plotLines(crustalTraces, Color.BLACK, 1f);
		
		CPT depthCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance();
		depthCPT = depthCPT.rescale(0d, 50d);
		mapMaker.plotScatterScalars(directionLocs, depths, depthCPT, "Depth (km)");
		mapMaker.plot(outputDir, prefix+"_depth", fmTitle);
		mapMaker.clearLines();
		
		mapMaker.clearScatters();
		
		mapMaker.setFillSurfaces(true);
		
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 8d);
//		CPT slipUncertCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 2d);
		CPT rakeCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, 90d);
		
//		mapMaker.plotSectScalars(slipUncerts, slipUncertCPT, "Slip Rate Uncertainty (mm/yr)");
//		mapMaker.plot(outputDir, prefix+"_slip_uncerts", fmDmTitle);
		
		if (!WRITE_SECT_INDEXES) {
			// add vectors for rake and slip
			List<LocationList> arrows = new ArrayList<>();
			for (int s=0; s<sects.size(); s++) {
				arrows.addAll(PRVI_SubductionSubSectPlots.buildRakeArrows(sects.get(s), slipCPT.getMaxValue()));
			}
			System.out.println("Plotting "+arrows.size()+" arrow lines for "+prefix);
//			mapMaker.plotLines(arrows, Color.BLACK, 2f);
			mapMaker.plotArrows(arrows, 20d, Colors.tab_red.darker(), 2f);
			mapMaker.setFillArrowheads(true);
			
			mapMaker.plotSectScalars(slips, slipCPT, null);
			mapMaker.plot(outputDir, prefix+"_slips_no_cpt", fmDmTitle);
		}
		
		mapMaker.plotSectScalars(slips, slipCPT, "Slip Rate (mm/yr)");
		mapMaker.plot(outputDir, prefix+"_slips", fmDmTitle);
		mapMaker.plotSectScalars(rakes, rakeCPT, "Rake Angle");
		mapMaker.plot(outputDir, prefix+"_rakes", fmDmTitle);
		// now add annotations
		mapMaker.addAnnotations(PRVI_SubductionSubSectPlots.getLabelAnns(mapMaker));
		mapMaker.plotSectScalars(slips, slipCPT, "Slip Rate (mm/yr)");
		mapMaker.plot(outputDir, prefix+"_slips_ann", fmDmTitle);
	}
	
	private static void debugWriteInterpFuncsOrders(File outputDir, String prefix, String name,
			ArbitrarilyDiscretizedFunc dasSlipFunc, ArbitrarilyDiscretizedFunc dasSlipUncertFunc,
			ArbitrarilyDiscretizedFunc dasRakeFunc, List<Double> sectStartDASs,
			List<Double> sectMidDASs, List<double[]> interpRanges, List<String> sectShortNames) throws IOException {
		
		List<PlotSpec> specs = new ArrayList<>(3);
		List<Range> yRanges = new ArrayList<>(3);

		PlotCurveCharacterstics startChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK);
		PlotCurveCharacterstics middleChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY);
		
		for (int i=0; i<3; i++) {
			ArbitrarilyDiscretizedFunc inFunc;
			String yAxisLabel;
			switch (i) {
			case 0:
				inFunc = dasSlipFunc;
				yAxisLabel = "Slip Rate (mm/yr)";
				break;
			case 1:
				inFunc = dasSlipUncertFunc;
				yAxisLabel = "Slip Rate Uncertainty (mm/yr)";
				break;
			case 2:
				inFunc = dasRakeFunc;
				yAxisLabel = "Rake Angle";
				break;

			default:
				throw new IllegalStateException();
			}
			
			double maxY = inFunc.getMaxY()*1.05;
			yRanges.add(new Range(0d, maxY));
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			for (int s=0; s<sectMidDASs.size(); s++) {
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				double middle = sectMidDASs.get(s);
				xy.set(middle, 0d);
				xy.set(middle, maxY);
				
				funcs.add(xy);
				chars.add(middleChar);
				
				if (s == 0)
					xy.setName("Section Middles");
			}
			
			for (int s=1; s<sectStartDASs.size(); s++) {
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				double start = sectStartDASs.get(s);
				xy.set(start, 0d);
				xy.set(start, maxY);
				
				funcs.add(xy);
				chars.add(startChar);
				
				if (s == 1)
					xy.setName("Section Ends");
			}
			
			inFunc.setName("Interpolated");
			funcs.add(inFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
			
			DefaultXY_DataSet midFunc = new DefaultXY_DataSet();
			for (double das : sectMidDASs)
				midFunc.set(das, inFunc.getInterpolatedY(das));
			
			midFunc.setName("Raw Values");
			funcs.add(midFunc);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 9f, Color.GRAY));

			// add interpolation annotations
			int interpAlpha = 40;
			Color interpColor1 = Colors.tab_orange;
			Color interpColor2 = Colors.tab_blue;
			interpColor1 = new Color(interpColor1.getRed(), interpColor1.getGreen(), interpColor1.getBlue(), interpAlpha);
			interpColor2 = new Color(interpColor2.getRed(), interpColor2.getGreen(), interpColor2.getBlue(), interpAlpha);

			double interpMinY = 0d;
			double interpMaxY = maxY;

			List<XYAnnotation> anns = new ArrayList<>();

			for (double[] interpRange : interpRanges) {
				if (interpRange != null) {
					XYBoxAnnotation boxAnn = new XYBoxAnnotation(interpRange[0], interpMinY, interpRange[1], interpMaxY,
							null, null, anns.size() % 2 == 0 ? interpColor1 : interpColor2);
					anns.add(boxAnn);
				}
			}
			
			// annotate interpolation ranges
			// use darker color
			interpColor1 = new Color(interpColor1.getRed(), interpColor1.getGreen(), interpColor1.getBlue(), Integer.min(150, interpAlpha*4));
			interpColor2 = new Color(interpColor2.getRed(), interpColor2.getGreen(), interpColor2.getBlue(), Integer.min(150, interpAlpha*4));
			EvenlyDiscretizedFunc fakeHist = new EvenlyDiscretizedFunc(0d, 1d, 10);
			fakeHist.setName(" ");
			funcs.add(fakeHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, interpColor1));
			fakeHist = fakeHist.deepClone();
			fakeHist.setName("Interpolation Regions");
			funcs.add(fakeHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, interpColor2));
			
			PlotSpec spec = new PlotSpec(funcs, chars, name, "Distance Along Strike (km)", yAxisLabel);
			if (i == 0)
//				spec.setLegendInset(RectangleAnchor.TOP_LEFT);
				spec.setLegendVisible(true);

			spec.setPlotAnnotations(anns);
			
			specs.add(spec);
		}
		
		int width = 1000;
		int height = width;
		
		// on second though, we don't care about the uncertainty one
		specs.remove(1);
		yRanges.remove(1);
		height = 800;
		
		// add names
		List<Integer> specWeights = new ArrayList<>();
		for (int i=0; i<specs.size(); i++)
			specWeights.add(100);
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		List<XYTextAnnotation> anns = new ArrayList<>();
		double maxY = 10;
		for (int s=0; s<sectStartDASs.size(); s++) {
			String shortName = sectShortNames.get(s);
			
			double mid = sectMidDASs.get(s);
			
			XYTextAnnotation ann = new XYTextAnnotation(shortName, mid, 0.5*maxY);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			ann.setTextAnchor(TextAnchor.CENTER);
			anns.add(ann);
			
			if (s > 0) {
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				double start = sectStartDASs.get(s);
				xy.set(start, 0d);
				xy.set(start, maxY);
				
				funcs.add(xy);
				chars.add(startChar);
			}
		}
		PlotSpec nameSpec = new PlotSpec(funcs, chars, " ", specs.get(0).getXAxisLabel(), "");
		nameSpec.setPlotAnnotations(anns);
		specs.add(1, nameSpec);
		yRanges.add(1, new Range(0d, maxY));
		specWeights.add(1, 10);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		gp.drawGraphPanel(specs, false, false, List.of(new Range(0d, dasSlipFunc.getMaxX())), yRanges);
		
		PlotUtils.setSubPlotWeights(gp, Ints.toArray(specWeights));
		
		XYPlot subPlot = PlotUtils.getSubPlots(gp).get(1);
		subPlot.setDomainGridlinesVisible(false);
		subPlot.setRangeGridlinesVisible(false);
		subPlot.getRangeAxis().setTickLabelsVisible(false);
		
		PlotUtils.writePlots(outputDir, prefix, gp, width, height, true, true, false);
	}
}
