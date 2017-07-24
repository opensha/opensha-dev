package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.kevin.ucerf3.GriddedSeisImportanceHazardMapCalc.CalcType;

public class GriddedSeisImportanceHazardMapPlotter {
	
	public static void main(String[] args) throws Exception {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_07_14-ucerf3-gridded-tests");
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		// the point on the hazard curve we are plotting
		boolean isProbAtIML = false;
		double level = 0.0004;
		String durationLabel = "2% in 50 yrs";
		String fileLabel = "2pin50";
		boolean mapParallel = true;
		
		List<CalcType[]> ratioSets = new ArrayList<>();
		ratioSets.add(new CalcType[] {CalcType.FULL, CalcType.SUPRA_ONLY});
		ratioSets.add(new CalcType[] {CalcType.FULL, CalcType.SUPRA_PLUS_SUB_ONLY});
		ratioSets.add(new CalcType[] {CalcType.SUPRA_ONLY, CalcType.GRIDDED_ONLY});
		ratioSets.add(new CalcType[] {CalcType.SUPRA_PLUS_SUB_ONLY, CalcType.OFF_FAULT_ONLY});
		
		Map<CalcType, Map<Location, ? extends DiscretizedFunc>> curvesMap = new HashMap<>();
		Map<CalcType, GriddedGeoDataSet> xyzMap = new HashMap<>();
		
		for (CalcType type : CalcType.values()) {
			File curveFile = new File(mainDir, type.name().toLowerCase()+"/curves/imrs1.bin");
			
			BinaryHazardCurveReader curveReader = new BinaryHazardCurveReader(curveFile.getAbsolutePath());
			Map<Location, ArbitrarilyDiscretizedFunc> curves = curveReader.getCurveMap();
			
			curvesMap.put(type, curves);
			
			GriddedGeoDataSet data = new GriddedGeoDataSet(gridReg, false);
			
			for (Location loc : curves.keySet()) {
				DiscretizedFunc curve = curves.get(loc);
				double val = HazardDataSetLoader.getCurveVal(curve, isProbAtIML, level);
				data.set(loc, val);
			}
			
			xyzMap.put(type, data);
		}

		List<Thread> threads = new ArrayList<>();
		for (CalcType[] ratioSet : ratioSets) {
			Runnable run = new Runnable() {
				
				@Override
				public void run() {
					try {
						plotRatio(xyzMap, ratioSet[0], ratioSet[1], mainDir, durationLabel, fileLabel);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			};
			if (mapParallel) {
				Thread thread = new Thread(run);
				thread.start();
				threads.add(thread);
			} else {
				run.run();
			}
		}
		
		for (Thread thread : threads)
			thread.join();
	}
	
	private static void plotRatio(Map<CalcType, GriddedGeoDataSet> xyzMap,
			CalcType numerator, CalcType denominator, File outputDir, String durationLabel, String fileLabel)
					throws IOException, GMT_MapException {
		GriddedGeoDataSet xyz1 = xyzMap.get(numerator);
		GriddedGeoDataSet xyz2 = xyzMap.get(denominator);
		GriddedGeoDataSet ratioData = new GriddedGeoDataSet(xyz1.getRegion(), false);
		for (int i=0; i<xyz1.size(); i++)
			ratioData.set(i, xyz1.get(i) / xyz2.get(i));
		ratioData.log10();
		
		double min = -0.5;
		double max = 0.5;
		
		CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance();
		ratioCPT = ratioCPT.rescale(min, max);
		
		System.out.println("Creating map instance...");
		GMT_Map map = new GMT_Map(xyz1.getRegion(), ratioData, xyz1.getRegion().getLatSpacing(), ratioCPT);
		
		String label = "Log10("+numerator.getLabel()+"/"+denominator.getLabel()+")";
		String prefix = "ratio_"+numerator.name()+"_vs_"+denominator.name()+"_"+fileLabel;
		map.setCustomLabel(label);
//		map.setTopoResolution(TopographicSlopeFile.US_SIX);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setCoast(new CoastAttributes(Color.GRAY, 1));
		map.setLogPlot(false);
		map.setDpi(300);
		map.setCustomScaleMin(min);
		map.setCustomScaleMax(max);
		map.setRescaleCPT(false);
		map.setBlackBackground(false);
		
		System.out.println("Making map...");
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
	}

}
