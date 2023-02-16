package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;

public class ExternalHazardCompPageGen {
	
	private static enum CompType {
			TOTAL,
			FAULTS,
			GRIDDED,
			SUBDUCTION
	};
	
	public static void main(String[] args) throws IOException {
//		String imtName = "PGA";
//		String imtDir = "PGA";
		
		String imtName = "1s SA";
		String imtDir = "SA1P0";
		
		ReturnPeriods[] rps = ReturnPeriods.values();
		
		String name23 = "NSHM23";
		String name18 = "NSHM18";
		
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "nshmp-haz-comparisons-"+imtDir);
		File sourcesDir23 = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
				+ "conus-2023-erf-6a6-all-vs760-0p1-20230203-fb449bf7b93c9c/vs30-760/"+imtDir+"/source");
		File sourcesDir18 = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/ext_hazard_calcs/"
				+ "conus-2018-530-all-vs760-0p1-20230203-3e0458e688c967/vs30-760/"+imtDir+"/source");

		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		GriddedRegion mapReg = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousUS(),
				0.1d, GriddedRegion.ANCHOR_0_0);
		
		List<Location> curveLocs = new ArrayList<>();
		
		DiscretizedFunc[] fault23 = add(
				loadCurves(new File(sourcesDir23, "FAULT/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir23, "FAULT_CLUSTER/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir23, "FAULT_SYSTEM/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir23, "ZONE/curves.csv"), curveLocs)
				);
		DiscretizedFunc[] grid23 = loadCurves(new File(sourcesDir23, "GRID/curves.csv"), curveLocs);
		DiscretizedFunc[] sub23 = add(
				loadCurves(new File(sourcesDir23, "INTERFACE/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir23, "SLAB/curves.csv"), curveLocs)
				);
		DiscretizedFunc[] full23 = add(fault23, grid23, sub23);
		
		DiscretizedFunc[] fault18 = add(
				loadCurves(new File(sourcesDir18, "FAULT/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir18, "FAULT_CLUSTER/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir18, "FAULT_SYSTEM/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir18, "ZONE/curves.csv"), curveLocs)
				);
		DiscretizedFunc[] grid18 = loadCurves(new File(sourcesDir18, "GRID/curves.csv"), curveLocs);
		DiscretizedFunc[] sub18 = add(
				loadCurves(new File(sourcesDir18, "INTERFACE/curves.csv"), curveLocs),
				loadCurves(new File(sourcesDir18, "SLAB/curves.csv"), curveLocs)
				);
		DiscretizedFunc[] full18 = add(fault18, grid18, sub18);
		
		// these maps substitute a single nshm18 ingredient at a time, keeping everything else at 23
		DiscretizedFunc[] withFaults18 = add(fault18, grid23, sub23);
		DiscretizedFunc[] withGrid18 = add(fault23, grid18, sub23);
		DiscretizedFunc[] withSub18 = add(fault23, grid23, sub18);
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(List.of(), mapReg);
		mapMaker.setDefaultPlotWidth(1000);
		
		Color transparent = new Color(255, 255, 255, 0);
		
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 1);
		hazCPT.setNanColor(transparent);
		
		CPT pDiffCPT = CA_HazardChangeFigures.getCenterMaskedCPT(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance(), 10d, 50d);
		pDiffCPT.setNanColor(transparent);
		
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2d, 0.2d);
		diffCPT.setNanColor(transparent);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Hazard Comparisons, "+name23+" vs "+name18);
		lines.add("");
		
		lines.add("This page compares "+name23+" and "+name18+" hazard maps, computed externally with nshmp-haz. "
				+ "Each individual ERF contributor to hazard changes are plotted separately, but GMMs are held constant.");
		lines.add("");
		
		String csvStr = "Download map CSVs:";
		for (ReturnPeriods rp : rps) {
			System.out.println("Building CSV for "+rp);
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Location Index", "Latitude", "Longitude", name23, name18,
					name18+" Faults", name18+" Gridded", name18+" Subduction");
			
			GriddedGeoDataSet[] maps = {
					curvestoMap(full23, curveLocs, mapReg, rp),
					curvestoMap(full18, curveLocs, mapReg, rp),
					curvestoMap(withFaults18, curveLocs, mapReg, rp),
					curvestoMap(withGrid18, curveLocs, mapReg, rp),
					curvestoMap(withSub18, curveLocs, mapReg, rp)
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
		
		for (CompType type : CompType.values()) {
			DiscretizedFunc[] curves23, curves18;
			String label, description;
			String mapLabelAdd;
			
			switch (type) {
			case TOTAL:
				curves23 = full23;
				curves18 = full18;
				label = "Full Model Comparison";
				description = "Full model comparison, including fault, gridded, and subduction model changes.";
				mapLabelAdd = "";
				break;
			case FAULTS:
				curves23 = full23;
				curves18 = withFaults18;
				label = "Fault-Only Comparison";
				description = "Fault-only comparison, holding gridded and subduction sources constant (using those from "+name23+").";
				mapLabelAdd = "Faults, ";
				break;
			case GRIDDED:
				curves23 = full23;
				curves18 = withGrid18;
				label = "Gridded-Only Comparison";
				description = "Gridded-only comparison, holding fault and subduction sources constant (using those from "+name23+").";
				mapLabelAdd = "Gridded, ";
				break;
			case SUBDUCTION:
				curves23 = full23;
				curves18 = withSub18;
				label = "Subduction-Only Comparison";
				description = "Subduction-only comparison, holding fault and gridded sources constant (using those from "+name23+").";
				mapLabelAdd = "Subduction, ";
				break;

			default:
				throw new IllegalStateException();
			}
			
			lines.add("## "+label);
			lines.add(topLink); lines.add("");
			
			lines.add(description);
			lines.add("");
			
			for (ReturnPeriods rp : rps) {
				System.out.println("Plotting "+label+", "+rp.label);
				
				lines.add("### "+label+", "+rp.label);
				lines.add(topLink); lines.add("");
				
				String hazLabel = imtName+", "+rp.label;
				String prefix = type.name()+"_"+rp.name();
				
				GriddedGeoDataSet map23 = curvestoMap(curves23, curveLocs, mapReg, rp);
				GriddedGeoDataSet map18 = curvestoMap(curves18, curveLocs, mapReg, rp);
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				
				table.addLine(name23, name18);
				
				table.initNewLine();
				
				mapMaker.plotXYZData(asLog10(map23), hazCPT, mapLabelAdd+name23+", "+hazLabel+" (g)");
				mapMaker.plot(resourcesDir, prefix+"_nshm23", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_nshm23.png)");
				mapMaker.plotXYZData(asLog10(map18), hazCPT, mapLabelAdd+name18+", "+hazLabel+" (g)");
				mapMaker.plot(resourcesDir, prefix+"_nshm18", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_nshm18.png)");
				
				table.finalizeLine();
				
				table.addLine(MarkdownUtils.boldCentered("Ratio"), MarkdownUtils.boldCentered("Difference"));
				
				GriddedGeoDataSet pDiff = WUS_HazardChangePageGen.mapPDiff(map23, map18);
				GriddedGeoDataSet diff = WUS_HazardChangePageGen.mapDiff(map23, map18);
				
				table.initNewLine();
				
				mapMaker.plotXYZData(pDiff, pDiffCPT, mapLabelAdd+name23+" vs "+name18+", % Change, "+hazLabel);
				mapMaker.plot(resourcesDir, prefix+"_pDiff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_pDiff.png)");
				mapMaker.plotXYZData(diff, diffCPT, mapLabelAdd+name23+" - "+name18+", "+hazLabel+" (g)");
				mapMaker.plot(resourcesDir, prefix+"_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_diff.png)");
				
				table.finalizeLine();
				
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static DiscretizedFunc[] loadCurves(File csvFile, List<Location> locs) throws IOException {
		boolean locTest = !locs.isEmpty();
		
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		double[] xVals = new double[csv.getNumCols()-2];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = csv.getDouble(0, i+2);
		
		DiscretizedFunc[] curves = new DiscretizedFunc[csv.getNumRows()-1];
		
		for (int row=1; row<csv.getNumRows(); row++) {
			int index = row-1;
			double lon = csv.getDouble(row, 0);
			double lat = csv.getDouble(row, 1);
			Location loc = new Location(lat, lon);
			if (locTest) {
				Preconditions.checkState(locs.size() > index, "Input has %s locs, but encountered index %s", locs.size(), index);
				Location oLoc = locs.get(index);
				Preconditions.checkState(LocationUtils.areSimilar(oLoc, loc),
						"Locatio mismatch at index %s:\n\tInput: %s\n\tCSV:%s", index, oLoc, loc);
			} else {
				locs.add(loc);
			}
			double[] yVals = new double[xVals.length];
			for (int i=0; i<yVals.length; i++)
				yVals[i] = csv.getDouble(row, i+2);
			curves[index] = new LightFixedXFunc(xVals, yVals);
		}
		
		return curves;
	}
	
	private static DiscretizedFunc[] add(DiscretizedFunc[]... allCurves) {
		DiscretizedFunc[] ret = new DiscretizedFunc[allCurves[0].length];
		
		double[] xVals = new double[allCurves[0][0].size()];
		for (int i=0; i<xVals.length; i++)
			xVals[i] = allCurves[0][0].getX(i);
		
		for (int i=0; i<ret.length; i++) {
			double[] yVals = new double[xVals.length];
			for (DiscretizedFunc[] curves : allCurves) {
				DiscretizedFunc curve = curves[i];
				Preconditions.checkState(curve.size() == xVals.length);
				for (int j=0; j<xVals.length; j++) {
					Preconditions.checkState((float)xVals[j] == (float)curve.getX(j));
					yVals[j] += curve.getY(j);
				}
			}
			ret[i] = new LightFixedXFunc(xVals, yVals);
		}
		
		return ret;
	}
	
	private static GriddedGeoDataSet curvestoMap(DiscretizedFunc[] curves, List<Location> curveLocs,
			GriddedRegion gridRegion, ReturnPeriods rp) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(gridRegion);
		for (int i=0; i<ret.size(); i++)
			ret.set(i, Double.NaN);
		
		double rate = -Math.log(1d - rp.oneYearProb);
		
		for (int i=0; i<curveLocs.size(); i++) {
			Location loc = curveLocs.get(i);
			int index = gridRegion.indexForLocation(loc);
			if (index >= 0) {
				DiscretizedFunc curve = curves[i];
				double val;
				
				if (rate > curve.getMaxY())
					val = 0d;
				else if (rate < curve.getMinY())
					// saturated
					val = curve.getMaxX();
				else
					val = curve.getFirstInterpolatedX_inLogXLogYDomain(rate);
				ret.set(index, val);
			}
		}
		
		return ret;
	}
	
	private static GriddedGeoDataSet asLog10(GriddedGeoDataSet xyz) {
		xyz = zerosToNaNs(xyz);
		xyz.log10();
		return xyz;
	}
	
	private static GriddedGeoDataSet zerosToNaNs(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		for (int i=0; i<xyz.size(); i++)
			if (xyz.get(i) == 0d)
				xyz.set(i, Double.NaN);
		return xyz;
	}

}
