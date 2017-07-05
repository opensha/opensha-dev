package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.hazardMap.MakeXYZFromHazardMapDir;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class MapPlotGen {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws GMT_MapException 
	 */
	public static void main(String[] args) throws IOException, GMT_MapException {
//		File timeIndepDir = new File("/home/kevin/Simulators/maps/2013_04_10-rsqsim-long-la-box-0.05" +
//					"-pga-cb2008-30yrs/");
//		File timeDepDir = new File("/home/kevin/Simulators/maps/2013_04_10-rsqsim-long-la-box-0.05" +
//				"-pga-cb2008-30yrs-quietmojave156/");
//		
//		File writeDir = new File("/home/kevin/Simulators/maps/2013_04_10-plots");
		
//		File timeIndepDir = new File("/home/kevin/Simulators/maps/2013_04_11-rsqsim-long-la-box-0.05" +
//					"-1sec-cb2008-30yrs/");
//		File timeDepDir = new File("/home/kevin/Simulators/maps/2013_04_11-rsqsim-long-la-box-0.05" +
//					"-1sec-cb2008-30yrs-quietmojave156/");
//
//		File writeDir = new File("/home/kevin/Simulators/maps/2013_04_11-1sec-plots");
		
//		File timeIndepDir = new File("/home/kevin/Simulators/maps/2013_04_10-ucerf2-compare-time-indep");
//		File timeDepDir = new File("/home/kevin/Simulators/maps/2013_04_10-ucerf2-compare-time-dep-empirical");
//	
//		File writeDir = new File("/home/kevin/Simulators/maps/2013_04_10-ucerf2-plots-empirical");
		
		File timeIndepDir = new File("/home/kevin/Simulators/maps/2013_04_11-ucerf2-compare-time-1sec-indep");
		File timeDepDir = new File("/home/kevin/Simulators/maps/2013_04_11-ucerf2-compare-time-1sec-dep");
	
		File writeDir = new File("/home/kevin/Simulators/maps/2013_04_11-ucerf2-plots-1sec");
	
		if (!writeDir.exists())
			writeDir.mkdir();
		
		double[] probs = { 0.02, 0.1, 0.5 };
		double[] scaleMaxes = { 1d, 1d, 1d };
		
		double spacing = 0.01;
		
//		String imt = "PGA";
		String imt = "1sec SA";
		
		double duration = 30;
		
		Region region = new CaliforniaRegions.LA_BOX();
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		CPT ratioCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		ratioCPT = ratioCPT.rescale(0, 2);
		
		File timeIndepZip = new File(timeIndepDir, timeIndepDir.getName()+"_curves.zip");
		File timeDepZip = new File(timeDepDir, timeDepDir.getName()+"_curves.zip");
		
		for (int i=0; i<probs.length; i++) {
			double prob = probs[i];
			double scaleMax = scaleMaxes[i];
			
			cpt = cpt.rescale(0d, 1d);
			
			ArbDiscrGeoDataSet indepData = new MakeXYZFromHazardMapDir(
					timeIndepZip.getAbsolutePath(), true, true).getXYZDataset(false, prob);
			
			ArbDiscrGeoDataSet depData = new MakeXYZFromHazardMapDir(
					timeDepZip.getAbsolutePath(), true, true).getXYZDataset(false, prob);
			
			String probStr = imt+" (G), "+(int)(prob*100)+"% in "+(int)duration+" yrs";
			String fileName = imt.toLowerCase().replaceAll(" ", "_")+"_"+(int)(prob*100)+"in"+(int)duration;
			
			GMT_Map indepMap = new GMT_Map(region, indepData, spacing, cpt);
			
			indepMap.setCustomLabel("Time Independent "+probStr);
			indepMap.setTopoResolution(TopographicSlopeFile.CA_THREE);
			indepMap.setLogPlot(false);
			indepMap.setDpi(300);
			indepMap.setCustomScaleMin(0d);
			indepMap.setCustomScaleMax(scaleMax);
			
			FaultBasedMapGen.plotMap(writeDir, fileName+"_indep", false, indepMap);
			
			GMT_Map depMap = new GMT_Map(region, depData, spacing, cpt);
			
			depMap.setCustomLabel("Time Dependent "+probStr);
			depMap.setTopoResolution(TopographicSlopeFile.CA_THREE);
			depMap.setLogPlot(false);
			depMap.setDpi(300);
			depMap.setCustomScaleMin(0d);
			depMap.setCustomScaleMax(scaleMax);
			
			FaultBasedMapGen.plotMap(writeDir, fileName+"_dep", false, depMap);
			
			ArbDiscrGeoDataSet ratioData = new ArbDiscrGeoDataSet(true);

			for (int refInd=0; refInd<indepData.size(); refInd++) {
				Location refLoc = indepData.getLocation(refInd);
				double refZ = indepData.get(refInd);
				// if they're both sorted the same this should work
				if (depData.size() > refInd && refLoc.equals(depData.getLocation(refInd))) {
					double modZ = depData.get(refInd);
					double gain = modZ / refZ;

					ratioData.set(refLoc, gain);
				} else {
					for (int modInd=0; modInd<depData.size(); modInd++) {
						Location modLoc = depData.getLocation(modInd);
						double modZ = depData.get(modInd);

						if (!refLoc.equals(modLoc))
							continue;

						double gain = modZ / refZ;

						ratioData.set(modLoc, gain);
					}
				}
			}
			
			GMT_Map ratioMap = new GMT_Map(region, ratioData, spacing, ratioCPT);
			
			ratioMap.setCustomLabel("Dependent/Independent Ratio "+probStr);
			ratioMap.setTopoResolution(TopographicSlopeFile.CA_THREE);
			ratioMap.setLogPlot(false);
			ratioMap.setDpi(300);
			ratioMap.setCustomScaleMin((double)ratioCPT.getMinValue());
			ratioMap.setCustomScaleMax((double)ratioCPT.getMaxValue());
			
			FaultBasedMapGen.plotMap(writeDir, fileName+"_ratio", false, ratioMap);
		}
	}

}
