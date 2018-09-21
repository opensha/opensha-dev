package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.io.Files;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.kevin.simulators.hazard.HazardMapComparePlotter;
import scratch.kevin.ucerf3.GriddedSeisImportanceHazardMapCalc.CalcType;

public class GriddedSeisImportanceHazardMapPlotter {
	
	public static void main(String[] args) throws Exception {
//		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_29-ucerf3-geol-gridded-tests-pga");
//		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_29-ucerf3-geol-gridded-tests-sa-1.0s");
//		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_29-ucerf3-geol-gridded-tests-sa-5.0s");
//		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_29-ucerf3-geol-gridded-tests-sa-10.0s");
//		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_30-ucerf3-geol-gridded-tests-sa-0.1s");
//		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_30-ucerf3-geol-gridded-tests-sa-0.2s");
		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/maps/2018_08_30-ucerf3-geol-gridded-tests-sa-0.5s");
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		// the point on the hazard curve we are plotting
		boolean isProbAtIML = false;
		
		int[] rps = { 1000, 2500, 10000 };
		double[] levels = new double[rps.length];
		String[] durationLabels = new String[rps.length];
		String[] fileLabels = new String[rps.length];
		for (int i=0; i<rps.length; i++) {
			int rp = rps[i];
			levels[i] = 1d/(double)rp;
			durationLabels[i] = rp+"yr";
			fileLabels[i] = rp+"yr";
		}
//		double level = 1d/1000d;
//		String durationLabel = "1000yr";
//		String fileLabel = "1000yr";
//		double level = 1d/10000d;
//		String durationLabel = "10000yr";
//		String fileLabel = "10000yr";
//		double[] levels = { 0.0004 };
//		String[] durationLabels = { "2% in 50 yrs" };
//		String[] fileLabels = { "2pin50" };
		boolean mapParallel = true;
		
		List<CalcType[]> ratioSets = new ArrayList<>();
		ratioSets.add(new CalcType[] {CalcType.FULL, CalcType.SUPRA_ONLY});
		ratioSets.add(new CalcType[] {CalcType.FULL, CalcType.SUPRA_PLUS_SUB_ONLY});
//		ratioSets.add(new CalcType[] {CalcType.SUPRA_ONLY, CalcType.GRIDDED_ONLY});
//		ratioSets.add(new CalcType[] {CalcType.SUPRA_PLUS_SUB_ONLY, CalcType.OFF_FAULT_ONLY});
		
		HashSet<CalcType> types = new HashSet<>();
		for (CalcType[] ratioTypes : ratioSets)
			for (CalcType type : ratioTypes)
				types.add(type);
		
		CPT hazardCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		hazardCPT = hazardCPT.rescale(0d, 1.2d);
		String imtLabel;
		if (mainDir.getName().toLowerCase().contains(SA_Param.NAME.toLowerCase())) {
			String periodStr = mainDir.getName().toLowerCase();
			periodStr = periodStr.substring(0, periodStr.lastIndexOf("s"));
			periodStr = periodStr.substring(periodStr.lastIndexOf("-")+1);
			double period = Double.parseDouble(periodStr);
			if (period < 1d)
				hazardCPT = hazardCPT.rescale(0d, 2d);
			else if (period < 2d)
				hazardCPT = hazardCPT.rescale(0d, 1.2);
			else if (period < 5d)
				hazardCPT = hazardCPT.rescale(0d, 1d);
			else
				hazardCPT = hazardCPT.rescale(0d, 0.6);
			imtLabel = (float)period+"s SA (g)";
		} else {
			imtLabel = "PGA (g)";
		}
		hazardCPT.setNanColor(Color.WHITE);
		
		Map<CalcType, Map<Location, ? extends DiscretizedFunc>> curvesMap = new HashMap<>();
		Table<CalcType, Double, GriddedGeoDataSet> xyzMap = HashBasedTable.create();
		
		for (CalcType type : types) {
			File curveFile = new File(mainDir, type.name().toLowerCase()+"/curves/imrs1.bin");
			
			BinaryHazardCurveReader curveReader = new BinaryHazardCurveReader(curveFile.getAbsolutePath());
			Map<Location, ArbitrarilyDiscretizedFunc> curves = curveReader.getCurveMap();
			
			curvesMap.put(type, curves);
			
			for (int i=0; i<levels.length; i++) {
				GriddedGeoDataSet data = new GriddedGeoDataSet(gridReg, false);
				
				for (Location loc : curves.keySet()) {
					DiscretizedFunc curve = curves.get(loc);
					double val = HazardDataSetLoader.getCurveVal(curve, isProbAtIML, levels[i]);
					data.set(loc, val);
				}
				
				plotMap(data, hazardCPT, mainDir, type.getLabel()+", "+durationLabels[i], type.name()+"_"+fileLabels[i]);
				
				xyzMap.put(type, levels[i], data);
			}			
		}

		List<Thread> threads = new ArrayList<>();
		for (CalcType[] ratioSet : ratioSets) {
			for (int i=0; i<levels.length; i++) {
				final double level = levels[i];
				final String durationLabel = durationLabels[i];
				final String fileLabel = fileLabels[i];
				Runnable run = new Runnable() {
					
					@Override
					public void run() {
						try {
							plotRatio(xyzMap.column(level), ratioSet[0], ratioSet[1], mainDir, durationLabel, fileLabel);
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
			if (rps != null) {
				HazardMapComparePlotter.plotHists(curvesMap.get(ratioSet[0]), curvesMap.get(ratioSet[1]), ratioSet[1].getLabel(), gridReg,
					rps, -1, mainDir, false, imtLabel);
				String prefix = "hist_"+ratioSet[0].name()+"_vs_"+ratioSet[1].name();
				for (File file : mainDir.listFiles()) {
					String name = file.getName();
					if (name.startsWith("hist_") && !name.contains("_vs_")) {
						File dest = new File(mainDir, name.replaceAll("hist", prefix));
						System.out.println("Moving "+name+" to "+dest.getName());
						Files.move(file, dest);
					}
				}
			}
		}
		
		for (Thread thread : threads)
			thread.join();
	}
	
	private static void plotMap(GriddedGeoDataSet xyz, CPT hazardCPT, File outputDir, String label, String prefix)
					throws IOException, GMT_MapException {
		System.out.println("Creating map instance...");
		GMT_Map map = new GMT_Map(xyz.getRegion(), xyz, xyz.getRegion().getLatSpacing(), hazardCPT);
		
		map.setCustomLabel(label);
//		map.setTopoResolution(TopographicSlopeFile.US_SIX);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setCoast(new CoastAttributes(Color.GRAY, 1));
		map.setLogPlot(false);
		map.setDpi(300);
		map.setCustomScaleMin((double)hazardCPT.getMinValue());
		map.setCustomScaleMax((double)hazardCPT.getMaxValue());
		map.setRescaleCPT(false);
		map.setBlackBackground(false);
		
		System.out.println("Making map...");
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
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
