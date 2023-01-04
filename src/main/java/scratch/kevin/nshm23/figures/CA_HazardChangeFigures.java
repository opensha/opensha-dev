package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;

import com.google.common.base.Preconditions;

public class CA_HazardChangeFigures {
	
	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/u3_haz_change_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemRupSet rupSetU3 = FaultSystemRupSet.load(new File(
				"/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip"));
		FaultSystemRupSet rupSet23 = FaultSystemRupSet.load(new File(invsDir,
				"2022_12_07-nshm23_branches-no_paleo_slip-mod_dm_weights-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		File u3HazFile = new File(invsDir,
				"2022_12_08-u3-FM3_1-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File u3ConvergedHazFile = new File(invsDir,
				"2022_03_24-u3_branches-FM3_1-2000ip-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File u3_23GridHazFile = new File(invsDir,
				"2022_12_13-u3-FM3_1-ba_only-nshm23_gridded/"
				+ "results_hazard_include_0.1deg.zip");
		
		File methodsU3GridHazFile = new File(invsDir,
				"2022_12_06-nshm23_u3_hybrid_branches-no_paleo_slip-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		File methods23GridHazFile = new File(invsDir,
				"2022_12_06-nshm23_u3_hybrid_branches-no_paleo_slip-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only-nshm23_gridded/"
				+ "results_hazard_include_0.1deg.zip");
		
		File modelHazFile = new File(invsDir,
				"2022_12_07-nshm23_branches-no_paleo_slip-mod_dm_weights-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only/"
				+ "results_hazard_include_0.1deg.zip");
		
		String entryName = "mean_map_pga_TWO_IN_50.txt";
		String hazLabel = "PGA, 2% in 50 yrs";
		
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		Region[] caRegions = new Region[caOutlines.length];
		for (int i=0; i<caOutlines.length; i++) {
			LocationList outline = new LocationList();
			for (Point2D pt : caOutlines[i])
				outline.add(new Location(pt.getY(), pt.getX()));
			caRegions[i] = new Region(outline, BorderType.MERCATOR_LINEAR);
		}
		
		GriddedGeoDataSet u3Map = mask(caRegions, loadXYZ(u3HazFile, entryName));
		GriddedGeoDataSet u3ConvergedMap = mask(caRegions, loadXYZ(u3ConvergedHazFile, entryName));
		GriddedGeoDataSet u3_23GridMap = mask(caRegions, loadXYZ(u3_23GridHazFile, entryName));
		GriddedGeoDataSet methodsU3GridMap = mask(caRegions, loadXYZ(methodsU3GridHazFile, entryName));
		GriddedGeoDataSet methods23GridMap = mask(caRegions, loadXYZ(methods23GridHazFile, entryName));
		GriddedGeoDataSet modelMap = mask(caRegions, loadXYZ(modelHazFile, entryName));
		
		GriddedRegion refReg = u3Map.getRegion();
		RupSetMapMaker mapMakerU3 = new RupSetMapMaker(rupSetU3, refReg);
		RupSetMapMaker mapMaker23 = new RupSetMapMaker(rupSet23, refReg);
		
		for (RupSetMapMaker mapMaker : new RupSetMapMaker[] {mapMakerU3, mapMaker23}) {
			mapMaker.setSectOutlineChar(null);
			mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(100, 100, 100, 127)));
		}
		
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		pDiffCPT.setNanColor(new Color(255, 255, 255, 0));
		plotHazardChange(outputDir, "u3_converged_vs_u3", mapMakerU3, pDiffCPT, refReg, u3ConvergedMap, u3Map,
				"Converged UCERF3 vs UCERF3, % Change, "+hazLabel);
		plotHazardChange(outputDir, "methodology_vs_u3", mapMakerU3, pDiffCPT, refReg, methodsU3GridMap, u3Map,
				"NSHM23 Methodology vs UCERF3, % Change, "+hazLabel);
		plotHazardChange(outputDir, "ingredients", mapMaker23, pDiffCPT, refReg, modelMap, methods23GridMap,
				"NSHM23 vs UCERF3 Ingredients, % Change, "+hazLabel);
		plotHazardChange(outputDir, "full_change", mapMaker23, pDiffCPT, refReg, modelMap, u3_23GridMap,
				"NSHM23 vs UCERF3, % Change, "+hazLabel);
		
		// attribution
		CPT methodsCPT = getHalfCPT(GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance(), 10d, 50d);
		CPT modelCPT = getHalfCPT(GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse(), 10d, 50d);
		
		CPT attributionCPT = new CPT();
		for (int i=methodsCPT.size(); --i>=0;) {
			CPTVal val = methodsCPT.get(i);
			attributionCPT.add(new CPTVal(-val.end, val.maxColor, -val.start, val.minColor));
		}
		for (CPTVal val : modelCPT)
			attributionCPT.add(val);
		attributionCPT.setBelowMinColor(methodsCPT.getAboveMaxColor());
		attributionCPT.setAboveMaxColor(modelCPT.getAboveMaxColor());
		attributionCPT.setNanColor(new Color(255, 255, 255, 0));
//		System.out.println(combCPT);
		
		GriddedGeoDataSet attXYZ = new GriddedGeoDataSet(u3Map.getRegion(), false);
		
		int numMethod = 0;
		int numIngredient = 0;
		int numWithin = 0;
		for (int i=0; i<attXYZ.size(); i++) {
			double u3Val = u3Map.get(i);
			Location loc = attXYZ.getLocation(i);
			
			double u3_23GridVal = u3_23GridMap.get(u3_23GridMap.indexOf(loc));
			double methodsU3GridVal = methodsU3GridMap.get(methodsU3GridMap.indexOf(loc));
			double methods23GridVal = methods23GridMap.get(methods23GridMap.indexOf(loc));
			double modelVal = modelMap.get(modelMap.indexOf(loc));
			
			double totalDiff = Math.abs(modelVal - u3_23GridVal);
			double modelDiff = Math.abs(modelVal - methods23GridVal);
			double methodsDiff = Math.abs(methodsU3GridVal - u3Val);
			
			double pDiff = 100d*totalDiff/u3Val;
			if (modelDiff < methodsDiff)
				// when methods dominate, go negative
				pDiff = -pDiff;
			attXYZ.set(i, pDiff);
			if (Double.isFinite(pDiff)) {
				if (Math.abs(pDiff) > 10d) {
					if (modelDiff > methodsDiff)
						numIngredient++;
					else
						numMethod++;
				} else {
					numWithin++;
				}
			}
		}
		
		int totNum = numMethod + numIngredient + numWithin;
		System.out.println("Attribution stats:");
		System.out.println("\tWithin 10%: "+pDF.format((double)numWithin/(double)totNum));
		System.out.println("\tMethodology: "+pDF.format((double)numMethod/(double)totNum));
		System.out.println("\tIngredients: "+pDF.format((double)numIngredient/(double)totNum));
		
		mapMaker23.plotXYZData(attXYZ, attributionCPT, "Methodology     ←     |Hazard % Change|     →     Ingredients");
		mapMaker23.plot(outputDir, "comb_hazard_attribution", " ");
	}
	
	private static GriddedGeoDataSet loadXYZ(File zipFile, String entryName) throws IOException {
		ZipFile zip = new ZipFile(zipFile);
		
		ZipEntry regEntry = zip.getEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(regEntry)));
		GriddedRegion gridReg = GriddedRegion.fromFeature(Feature.read(bRead));
		
		ZipEntry mapEntry = zip.getEntry(entryName);
		InputStream is = zip.getInputStream(mapEntry);
		Preconditions.checkNotNull(is, "IS is null for %s", entryName);
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(mapEntry)));
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
		
		zip.close();
		
		return xyz;
	}
	
	private static CPT getHalfCPT(CPT cpt, double maskRange, double fullRange) {
		cpt = cpt.rescale(-1d, 1d);
		CPT ret = new CPT();
		Color zeroColor = cpt.getColor(0f);
		ret.add(new CPTVal(0f, zeroColor, (float)maskRange, zeroColor));
		// now scale up to full range
		for (double x=0.01; x<1d; x+=0.01) {
			double mappedX = maskRange + (fullRange-maskRange)*x;
			
			float prevX = ret.getMaxValue();
			Color prevColor = ret.getMaxColor();
			
			ret.add(new CPTVal(prevX, prevColor, (float)mappedX, cpt.getColor((float)x)));
		}
		ret.add(new CPTVal(ret.getMaxValue(), ret.getMaxColor(), (float)fullRange, cpt.getMaxColor()));
		ret.setBelowMinColor(ret.getMinColor());
		ret.setAboveMaxColor(ret.getMaxColor());
		ret.setNanColor(Color.WHITE);
		return ret;
	}
	
	static GriddedGeoDataSet mask(Region[] includeRegions, GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		GriddedRegion gridReg = xyz.getRegion();
		double halfLatSpacing = gridReg.getLatSpacing()*0.5;
		double halfLonSpacing = gridReg.getLonSpacing()*0.5;
		for (int i=0; i<xyz.size(); i++) {
			boolean inside = false;
			Location center = xyz.getLocation(i);
			
			// first just check center
			for (Region reg : includeRegions) {
				if (reg.contains(center)) {
					inside = true;
					break;
				}
			}
			if (!inside) {
				// test full cell intersection
				Location upLeft = new Location(center.getLatitude()+halfLatSpacing, center.getLongitude()-halfLonSpacing);
				Location botRight = new Location(center.getLatitude()-halfLatSpacing, center.getLongitude()+halfLonSpacing);
				Region cell = new Region(upLeft, botRight);
				
				for (Region reg : includeRegions) {
					double dist = reg.distanceToLocation(center);
					if (dist < 30d && Region.intersect(cell, reg) != null) {
						inside = true;
						break;
					}
				}
			}
			if (!inside)
				xyz.set(i, Double.NaN);
		}
		return xyz;
	}
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static void plotHazardChange(File outputDir, String prefix, RupSetMapMaker mapMaker, CPT pDiffCPT,
			GriddedRegion refReg, GriddedGeoDataSet numerator, GriddedGeoDataSet denominator, String label) throws IOException {
		GriddedGeoDataSet pDiff = new GriddedGeoDataSet(refReg, false);
		
		MinMaxAveTracker meanAbsTrack = new MinMaxAveTracker();
		MinMaxAveTracker meanTrack = new MinMaxAveTracker();
		int numWithin10 = 0;
		int numWithin5 = 0;
		int numWithin1 = 0;
		int numValid = 0;
		for (int i=0; i<pDiff.size(); i++) {
			Location loc = pDiff.getLocation(i);
			
			double v1 = numerator.get(numerator.indexOf(loc));
			double v2 = denominator.get(denominator.indexOf(loc));
			
			if (Double.isFinite(v1) && Double.isFinite(v2)) {
				double val = 100d*(v1-v2)/v2;
				pDiff.set(i, val);
				
				numValid++;
				if (Math.abs(val) < 10d)
					numWithin10++;
				if (Math.abs(val) < 5d)
					numWithin5++;
				if (Math.abs(val) < 1d)
					numWithin1++;
				
				meanAbsTrack.addValue(Math.abs(val));
				meanTrack.addValue(val);
			} else {
				pDiff.set(i, Double.NaN);
			}
		}
		
		System.out.println("Plotting "+prefix+", "+label);
		System.out.println("\tRange: ["+twoDigits.format(meanTrack.getMin())+"%, "+twoDigits.format(meanTrack.getMax())+"%]");
		System.out.println("\tAverage: "+twoDigits.format(meanTrack.getAverage())+"%");
		System.out.println("\tAverage Absolute: "+twoDigits.format(meanAbsTrack.getAverage())+"%");
		System.out.println("\tWithin 1%: "+pDF.format((double)numWithin1/(double)numValid));
		System.out.println("\tWithin 5%: "+pDF.format((double)numWithin5/(double)numValid));
		System.out.println("\tWithin 10%: "+pDF.format((double)numWithin10/(double)numValid));
		
		
		mapMaker.plotXYZData(pDiff, pDiffCPT, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}

}