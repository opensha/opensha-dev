package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.kevin.latex.LaTeXUtils;
import scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.Models;

public class HazardMapFigures {

	public static void main(String[] args) throws IOException {
		double[] periods = { 0d, 1d };
		ReturnPeriods[] rps = { ReturnPeriods.TWO_IN_50 };
		
		double[] texPeriods = { 0d, 1d };
		ReturnPeriods texRP = ReturnPeriods.TWO_IN_50;
		
		boolean replot = false;
		
		String mapZipName = "results_hazard_INCLUDE.zip";
		
		File hazardDir = new File(FIGURES_DIR, "hazard_calcs");
		Preconditions.checkState(hazardDir.exists() || hazardDir.mkdir());
		
		Map<Models, List<Models>> comparisons = new HashMap<>();
		
		Table<Models, Models, String> compLabels = HashBasedTable.create();
		
		Models[] models = Models.values();
		Models prevModel = null;
		for (Models model : models) {
			if (model != REF_FINITE_MODEL && model.ordinal() <= PROPOSED_DIST_CORR_MODEL.ordinal())
				// this is a distance correction/point surface representation test
				// include a comparison with our reference model
				compAdd(comparisons, model, REF_FINITE_MODEL);
			if (model.ordinal() >= PROPOSED_DIST_CORR_MODEL.ordinal()) {
				// this is a further grid property test
				// include a reference with original grid props, but proposed dist corr
				compAdd(comparisons, model, PROPOSED_FULL_MODEL);
				if (prevModel != PROPOSED_FULL_MODEL)
					// also add the incremental change from just this compared to the previous
					compAdd(comparisons, model, prevModel);
			}
			prevModel = model;
		}
		
		// published vs full proposed
		compAdd(comparisons, Models.AS_PUBLISHED, PROPOSED_FULL_MODEL);
		
		// alt random and realization count tests
		compAdd(comparisons, Models.FINITE_1X_UNCENTERED, Models.FINITE_1X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_2X_UNCENTERED, Models.FINITE_2X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_5X_UNCENTERED, Models.FINITE_5X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_20X_UNCENTERED, Models.FINITE_20X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_50X_UNCENTERED, Models.FINITE_50X_UNCENTERED_ALT_RAND);
		compAdd(comparisons, Models.FINITE_100X_UNCENTERED, Models.FINITE_100X_UNCENTERED_ALT_RAND);
		
		compAdd(comparisons, PROPOSED_DIST_CORR_MODEL, Models.AS_PUBLISHED);
		compLabels.put(PROPOSED_DIST_CORR_MODEL, Models.AS_PUBLISHED, "Proposed corrections");
		compAdd(comparisons, PROPOSED_FULL_MODEL, Models.AS_PUBLISHED);
		compLabels.put(PROPOSED_FULL_MODEL, Models.AS_PUBLISHED, "Proposed properties & corrections");
		
//		// + improved Rrup and Rx
//		compAdd(comparisons, Models.SPINNING_AVG_M6, Models.AS_PUBLISHED);
//		// + lower magnitude-threshold
//		compAdd(comparisons, Models.SPINNING_AVG, Models.SPINNING_AVG_M6);
//		// + uncentered
//		compAdd(comparisons, Models.SPINNING_AVG_UNCENTERED, Models.SPINNING_AVG);
//
//		// effect of multiple realizations
//		compAdd(comparisons, Models.FINITE_20X_CENTERED, Models.SPINNING_AVG);
//		// effect of uncentering
//		compAdd(comparisons, Models.FINITE_20X_UNCENTERED, Models.FINITE_20X_CENTERED);
//		
//		// new model vs higher resolution
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED, Models.FINITE_50X_UNCENTERED);
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED, Models.FINITE_100X_UNCENTERED);
//		
//		// modeling choice comparisons
//		// M3 vs M4
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_M3, Models.SPINNING_DIST_5X_UNCENTERED_M4);
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_M3, Models.SPINNING_DIST_5X_UNCENTERED);
//		// M4 vs M5
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_M4, Models.SPINNING_DIST_5X_UNCENTERED);
//		// alt depth profile (interpolated, not shelved)
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_ALT_DEPTH, Models.SPINNING_DIST_5X_UNCENTERED);
//		// alt WC94 lengths
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_ALT_LENGTH, Models.SPINNING_DIST_5X_UNCENTERED);
//		// no SS
//		compAdd(comparisons, Models.SPINNING_DIST_5X_UNCENTERED_NO_SS, Models.SPINNING_DIST_5X_UNCENTERED);
		
		ZipFile[] modelZips = new ZipFile[models.length];
		ZipFile[] modelZoomZips = new ZipFile[models.length];
		for (int i=0; i<models.length; i++) {
			Models model = models[i];
			File mapFile = new File(model.mapDir, mapZipName);
			if (mapFile.exists())
				modelZips[i] = new ZipFile(mapFile);
			File zoomFile = new File(model.zoomMapDir, mapZipName);
			if (zoomFile.exists())
				modelZoomZips[i] = new ZipFile(zoomFile);
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(FULL_GRID_REG);
		GeographicMapMaker zoomMapMaker = new GeographicMapMaker(ZOOM_GRID_REG);
		
		List<? extends FaultSection> allSects = NSHM23_FaultModels.WUS_FM_v3.getFaultSections();
		
		for (GeographicMapMaker map : List.of(mapMaker, zoomMapMaker)) {
			map.setFaultSections(allSects);
			map.setSectOutlineChar(null);
			map.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
//			map.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(127, 127, 127, 100)));
			map.setWriteGeoJSON(false);
		}
//		mapMaker.setSectOutlineChar(null);
//		zoomMapMaker.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(127, 127, 127, 100)));
		
		Color trans = new Color(255, 255, 255, 0);
		
		CPT hazardCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		hazardCPT.setLog10(true);
		hazardCPT.setNanColor(trans);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-20, 20);
		pDiffCPT.setNanColor(trans);
		CPT pDiffTightCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-5, 5);
		pDiffTightCPT.setNanColor(trans);
		double maxForTight = 6d;
		
		GriddedGeoDataSet landMask = buildLandMask(FULL_GRID_REG);
		
		for (int p=0; p<periods.length; p++) {
			String perPrefix, perLabel, perUnits;
			String perTexPrefix;
			if (periods[p] == 0d) {
				perPrefix = "pga";
				perLabel = "PGA";
				perUnits = "(g)";
				perTexPrefix = "PGA";
			} else {
				perPrefix = (float)periods[p]+"s";
				perLabel = (float)periods[p]+"s SA";
				perUnits = "(g)";
				Preconditions.checkState(periods[p] == 1d);
				perTexPrefix = "SAOne";
			}
			
			for (ReturnPeriods rp : rps) {
				GriddedGeoDataSet[] maps = new GriddedGeoDataSet[models.length];
				GriddedGeoDataSet[] zoomMaps = new GriddedGeoDataSet[models.length];
				
				String entryName = "map_"+perPrefix+"_"+rp.name()+".txt";
				
				System.out.println("Loading maps for "+perLabel+", "+rp.label);
				
				for (int i=0; i<models.length; i++) {
					if (modelZips[i] != null)
						maps[i] = loadXYZ(modelZips[i], FULL_GRID_REG, entryName, landMask);
					if (modelZoomZips[i] != null)
						zoomMaps[i] = loadXYZ(modelZoomZips[i], ZOOM_GRID_REG, entryName, null);
				}
				
				String mapLabel = perLabel+" "+perUnits+", "+rp.label;
				String mapDiffLabel = "% Change, "+perLabel+", "+rp.label;
				
				File subDir = new File(hazardDir, perPrefix+"_"+rp.name());
				Preconditions.checkState(subDir.exists() || subDir.mkdir());
				
				File mapsDir = new File(subDir, "maps");
				Preconditions.checkState(mapsDir.exists() || mapsDir.mkdir());
				File mapCompDir = new File(subDir, "comparisons");
				Preconditions.checkState(mapCompDir.exists() || mapCompDir.mkdir());
				
				File zoomMapsDir = new File(subDir, "maps_zoom");
				Preconditions.checkState(zoomMapsDir.exists() || zoomMapsDir.mkdir());
				File zoomMapCompDir = new File(subDir, "comparisons_zoom");
				Preconditions.checkState(zoomMapCompDir.exists() || zoomMapCompDir.mkdir());
				
				FileWriter texFW = null;
				boolean plotMaps = false;
				if (Doubles.contains(texPeriods, periods[p]) && rp == texRP) {
					plotMaps = true;
					texFW = new FileWriter(new File(subDir, "hazard_change_stats.tex"));
				}
				
				CSVFile<String> compCSV = new CSVFile<>(true);
				compCSV.addLine("Model", "Comparison Model", "Mean Change", "Mean Absolute Change", "Median Absolute Change", "Maximum Change");
				
				for (int i=0; i<models.length; i++) {
					Models model = models[i];
					if (maps[i] == null) {
						System.err.println("Skipping "+models[i].getName()+", calculations not found");
						continue;
					}
					System.out.println("Plotting maps for "+models[i].getName());
					if (plotMaps) {
						if (replot || !new File(mapsDir, model.name()+".pdf").exists()) {
							mapMaker.plotXYZData(maps[i], hazardCPT, model.getName()+", "+mapLabel);
							mapMaker.plot(mapsDir, model.name(), " ");
						}
						if (zoomMaps[i] != null && (replot || !new File(zoomMapsDir, model.name()+".pdf").exists())) {
							zoomMapMaker.plotXYZData(zoomMaps[i], hazardCPT, model.getName()+", "+mapLabel);
							zoomMapMaker.plot(zoomMapsDir, model.name(), " ");
						}
					}
					
					List<Models> comps = comparisons.get(model);
					if (comps != null) {
						for (Models comp : comps) {
							int compIndex = comp.ordinal();
							System.out.println("\tVs "+comp.getName());
							if (maps[compIndex] == null) {
								System.err.println("\tSkipping comparison with "+models[compIndex].getName()+", calculations not found");
								continue;
							}
							GriddedGeoDataSet pDiff = pDiff(maps[i], maps[compIndex]);
							DiffStats stats = new DiffStats(pDiff);
							
							GriddedGeoDataSet zoomPDiff = null;
							double mean = stats.mean;
							double meanAbs = stats.meanAbs;
							double medianAbs = stats.medianAbs;
							double maxSigned = stats.maxSigned;
							if (zoomMaps[i] != null && zoomMaps[compIndex] != null) {
								zoomPDiff = pDiff(zoomMaps[i], zoomMaps[compIndex]);
								DiffStats zoomStats = new DiffStats(zoomPDiff);
								if (Math.abs(zoomStats.maxSigned) > Math.abs(maxSigned))
									maxSigned = zoomStats.maxSigned;
							}
							
							CPT cpt = Math.abs(maxSigned) > maxForTight ? pDiffCPT : pDiffTightCPT;
							
//							System.out.println("\t\t"+stats);
							System.out.println("\t\tmean="+twoDF.format(mean)
								+"%;\tmeanAbs="+twoDF.format(meanAbs)
								+"%;\tmedianAbs="+twoDF.format(medianAbs)
								+"%;\tmaxSigned="+twoDF.format(maxSigned)+"%");
							
							if (plotMaps) {
								String prefix = model.name()+"_vs_"+comp.name();
//								String label = model.getName()+" vs "+comp.getName()+", "+mapDiffLabel;
								String label = mapDiffLabel;
								if (compLabels.contains(model, comp))
									label = compLabels.get(model, comp)+", "+label;
								
								if (replot || !new File(mapCompDir, prefix+".pdf").exists()) {
									mapMaker.plotXYZData(pDiff, cpt, label);
									mapMaker.plot(mapCompDir, prefix, " ");
								}
								if (zoomPDiff != null && (replot || !new File(zoomMapCompDir, prefix+".pdf").exists())) {
									zoomMapMaker.plotXYZData(zoomPDiff, cpt, label);
									zoomMapMaker.plot(zoomMapCompDir, prefix, " ");
								}
							}
							
							compCSV.addLine(model.getName(), comp.getName(), twoDF.format(mean)+"%",
									twoDF.format(meanAbs)+"%", twoDF.format(medianAbs)+"%", twoDF.format(maxSigned)+"%");
							
							if (texFW != null) {
								String texPrefix = model.texName+"Vs"+comp.texName+perTexPrefix;
								texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"Mean", oneDF.format(mean)+"%")+"\n");
								texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MeanAbs", oneDF.format(meanAbs)+"%")+"\n");
								texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MedianAbs", oneDF.format(medianAbs)+"%")+"\n");
								texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"Max", oneDF.format(maxSigned)+"%")+"\n");
								
								if (model == Models.AS_PUBLISHED && (comp == REF_FINITE_MODEL || comp == PROPOSED_FULL_MODEL)) {
									// add extra rounded values
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MeanRounded", roundedDF.format(mean)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MeanAbsRounded", roundedDF.format(meanAbs)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MedianAbsRounded", roundedDF.format(medianAbs)+"%")+"\n");
									texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MaxRounded", roundedDF.format(maxSigned)+"%")+"\n");
								}
							}
						}
					}
				}
				
				compCSV.writeToFile(new File(subDir, "hazard_change_stats.csv"));
				
				if (texFW != null)
					texFW.close();
			}
		}
	}
	
	static GriddedGeoDataSet loadXYZ(ZipFile zip, GriddedRegion gridReg, String entryName, GriddedGeoDataSet mask) throws IOException {
		System.out.println("Loading "+entryName+" from "+zip.getName());
		ZipEntry mapEntry = zip.getEntry(entryName);
		InputStream is = zip.getInputStream(mapEntry);
		Preconditions.checkNotNull(is, "IS is null for %s", entryName);
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(mapEntry)));
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
		
		if (mask != null) {
			Preconditions.checkState(mask.size() == xyz.size());
			for (int i=0; i<mask.size(); i++)
				xyz.set(i, xyz.get(i)*mask.get(i));
		}
		
		return xyz;
	}
	
	private static GriddedGeoDataSet pDiff(GriddedGeoDataSet map, GriddedGeoDataSet refMap) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map.getRegion());
		for (int i=0; i<ret.size(); i++) {
			double val1 = map.get(i);
			double val2 = refMap.get(i);
			if (!Double.isFinite(val1) || !Double.isFinite(val2)) {
				ret.set(i, Double.NaN);
			} else if (val1 == 0d && val2 == 0d) {
				ret.set(i, 0d);	
			} else {
				ret.set(i, 100d*(val1 - val2)/val2);
			}
		}
		return ret;
	}

	private static DecimalFormat roundedDF = new DecimalFormat("0");
	private static DecimalFormat oneDF = new DecimalFormat("0.0");
	private static DecimalFormat twoDF = new DecimalFormat("0.00");
	
	private static class DiffStats {
		
		public final double mean;
		public final double meanAbs;
		public final double medianAbs;
		public final double maxSigned;
		public final int numFinite;
		
		public DiffStats(GriddedGeoDataSet pDiff) {
			double sum = 0d;
			int numFinite = 0;
			double sumAbs = 0d;
			double largest = 0;
			double largestAbs = 0d;
			
			List<Double> absVals = new ArrayList<>(pDiff.size());
			
			for (int i=0; i<pDiff.size(); i++) {
				double val = pDiff.get(i);
				if (Double.isFinite(val)) {
					numFinite++;
					sum += val;
					double absVal = Math.abs(val);
					absVals.add(absVal);
					sumAbs += absVal;
					if (absVal > largestAbs) {
						largestAbs = absVal;
						largest = val;
					}
					
				}
			}
			
			mean = sum/(double)numFinite;
			meanAbs = sumAbs/(double)numFinite;
			medianAbs = DataUtils.median(Doubles.toArray(absVals));
			maxSigned = largest;
			this.numFinite = numFinite;
		}
		
		public String toString() {
			return "mean="+twoDF.format(mean)+"%;\tmeanAbs="+twoDF.format(meanAbs)+"%;\tmaxSigned="+twoDF.format(maxSigned)+"%";
		}
	}
	
	private static void compAdd(Map<Models, List<Models>> comparisons, Models from, Models to) {
		List<Models> comps = comparisons.get(from);
		if (comps == null) {
			comps = new ArrayList<>(1);
			comparisons.put(from, comps);
		}
		if (!comps.contains(to))
			comps.add(to);
	}

}
