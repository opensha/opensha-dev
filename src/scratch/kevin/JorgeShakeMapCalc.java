package scratch.kevin;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.xyz.AbstractGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.ScenarioShakeMapCalculator;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.mod.ModAttenRelRef;
import org.opensha.sha.imr.mod.ModAttenuationRelationship;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.Component;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;
import org.opensha.sha.mapping.GMT_MapGeneratorForShakeMaps;
import org.opensha.sha.util.component.ComponentConverter;
import org.opensha.sha.util.component.ComponentTranslation;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class JorgeShakeMapCalc {

	public static void main(String[] args) throws IOException, DocumentException, GMT_MapException {
		int fssIndex = 227088;
		File mainDir = new File("/home/kevin/OpenSHA/UCERF3/jorge_shakemap/");
		File vs30File = new File(mainDir, "thompson_hyb_tj.xyz");
		boolean directivity = false;
		boolean rotd100 = true;
		String dirName = "results";
		if (directivity)
			dirName += "_directivity";
		if (rotd100)
			dirName += "_rotd100";
		File resultsDir = new File(mainDir, dirName);
		Preconditions.checkState(resultsDir.exists() || resultsDir.mkdir());
		// -1 = PGV
		// 0 = PGA
		double[] periods = { -1d, 0d, 0.3, 1d, 3d };
		double[] cptMax = { 100d, 0.5d, 1d, 0.6d, 0.4d };
		Location hypo = new Location(33.01, -117.32, 7.7);
		double inc = 0.008333333;
		
		ModAttenRelRef directivityModel = ModAttenRelRef.BAYLESS_SOMERVILLE_2013_DIRECTIVITY;
		
		FaultSystemSolution fss = FaultSystemIO.loadSol(new File(
				"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fss);
		erf.updateForecast();
		ProbEqkRupture rup = erf.getSource(erf.getSrcIndexForFltSysRup(fssIndex)).getRupture(0);
		System.out.println("Rupture: M="+rup.getMag());
		if (directivity) {
			double minDist = Double.POSITIVE_INFINITY;
			Location closest = null;
			for (Location loc : rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface()) {
				double dist = LocationUtils.linearDistanceFast(hypo, loc);
				if (dist < minDist) {
					minDist = dist;
					closest = loc;
				}
			}
			System.out.println("Closest hypo on surface is "+(float)minDist+" km away: "+closest);
			rup.setHypocenterLocation(closest);
		}
		
		FileWriter metaFW = new FileWriter(new File(mainDir, "metadata.txt"));
		
		metaFW.write("Calculated with OpenSHA class "+JorgeShakeMapCalc.class.getName()+" at "+new Date().toString()+"\n\n");
		metaFW.write("UCERF3 Scenario Rupture: FM3.1 Rupture "+fssIndex+"\n");
		metaFW.write("Mag: "+rup.getMag()+"\n");
		metaFW.write("Rake: "+rup.getAveRake()+"\n");
		if (directivity) {
			metaFW.write("Directivity model: "+directivityModel.getName()+"\n");
			metaFW.write("Hypocenter for Directivity Calcs: "+rup.getHypocenterLocation()+"\n");
		}
		metaFW.write("UCERF3 Subsection For Rupture:\n");
		for (FaultSectionPrefData sect : fss.getRupSet().getFaultSectionDataForRupture(fssIndex))
			metaFW.write("\t"+sect.getName()+"\n");
		metaFW.write("Fault Trace:\n");
		for (Location loc : rup.getRuptureSurface().getUpperEdge())
			metaFW.write("\t"+loc+"\n");
		metaFW.close();
		
		List<AttenRelRef> refs = Lists.newArrayList();
		
		refs.add(AttenRelRef.ASK_2014);
		refs.add(AttenRelRef.BSSA_2014);
		refs.add(AttenRelRef.CB_2014);
		refs.add(AttenRelRef.CY_2014);
		refs.add(AttenRelRef.IDRISS_2014);
		refs.add(AttenRelRef.NGAWest_2014_AVG);
		refs.add(AttenRelRef.NGAWest_2014_AVG_NOIDRISS);
		
		List<Site> sites = loadSites(vs30File);
		
//		MinMaxAveTracker latTrack = new MinMaxAveTracker();
//		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
//		for (Site site : sites) {
//			Location loc = site.getLocation();
//			latTrack.addValue(loc.getLatitude());
//			lonTrack.addValue(loc.getLongitude());
//		}
//		Region reg = new Region(new Location(latTrack.getMax(), lonTrack.getMax()),
//				new Location(latTrack.getMin(), lonTrack.getMin()));
		Region reg = new Region(new Location(33.35, -117.8), new Location(32.2, -116.6));
		
		System.out.println("Loaded "+sites.size()+" sites");
		
		CPT cpt = GMT_CPT_Files.SHAKEMAP.instance();
		
		ComponentTranslation trans = null;
		if (rotd100)
			trans = ComponentConverter.getConverter(Component.RotD50, Component.RotD100);
		
		for (AttenRelRef ref : refs) {
			AttenuationRelationship imr;
			
			if (directivity) {
				imr = new ModAttenuationRelationship(ref, directivityModel);
			} else {
				imr = ref.instance(null);
			}
			imr.setParamDefaults();
			
			for (int i=0; i<periods.length; i++) {
				double period = periods[i];
				if (directivity && period <= 0) {
					System.out.println("Skipping non S(a) for directivity");
					continue;
				}
				// use ref because directivity case changes to ModAttenRel
				String prefix = ref.getShortName();
				String label = ref.getShortName().replaceAll("_", " ");
				CPT myCPT = cpt.rescale(0d, cptMax[i]);
				if (period == 0) {
					imr.setIntensityMeasure(PGA_Param.NAME);
					prefix += "_pga";
					label += " PGA";
				} else if (period == -1) {
					if (!imr.isIntensityMeasureSupported(PGV_Param.NAME)) {
						System.out.println(imr.getShortName()+" doesn't support PGV, skipping...");
						continue;
					}
					imr.setIntensityMeasure(PGV_Param.NAME);
					prefix += "_pgv";
					label += " PGV";
				} else {
					Preconditions.checkState(period > 0d);
					imr.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), period);
					prefix += "_sa"+(float)period;
					label += " "+(float)period+"s SA";
				}
				
				if (directivity) {
					prefix += "_directivity";
					label += " w/ Directivity";
				}
				if (rotd100) {
					prefix += "_rotd100";
					label += ", RotD100";
				}
				
				System.out.println("Calculating: "+prefix);
				
				ScenarioShakeMapCalculator calc = new ScenarioShakeMapCalculator();
				GeoDataSet xyz = calc.getScenarioShakeMapData(
						Lists.newArrayList(imr), Lists.newArrayList(1d), sites, rup, false, 0.5);
				xyz.exp();
				if (trans != null && period > 0) {
					double factor = trans.getScalingFactor(period);
					System.out.println("Scaling by factor of "+factor+" to RotD100");
					xyz.scale(factor);
				}
				File xyzFile = new File(resultsDir, prefix+".txt");
				System.out.println("Writing: "+xyzFile.getName());
				AbstractGeoDataSet.writeXYZFile(xyz, xyzFile);
				
				GMT_Map map = new GMT_Map(reg, xyz, inc, myCPT);
				map.setCustomScaleMin((double)myCPT.getMinValue());
				map.setCustomScaleMax((double)myCPT.getMaxValue());
				map.setCustomLabel(label);
				CompoundSurface surf = (CompoundSurface)rup.getRuptureSurface();
				for (RuptureSurface gridSurf : surf.getSurfaceList())
					GMT_MapGeneratorForShakeMaps.addRupture(map, (EvenlyGriddedSurface)gridSurf, rup.getHypocenterLocation(),
							GMT_MapGeneratorForShakeMaps.RUP_PLOT_PARAM_PERIMETER);
				
				System.out.println("Generating maps");
				FaultBasedMapGen.plotMap(resultsDir, prefix, false, map);
			}
		}
	}
	
	private static List<Site> loadSites(File vs30File) throws IOException {
		List<Site> sites = Lists.newArrayList();
		
		BufferedReader read = new BufferedReader(new FileReader(vs30File));
		
		String line;
		
		while ((line = read.readLine()) != null) {
			line = line.trim();
			if (line.isEmpty())
				continue;
			String[] split = line.split("\t");
			Preconditions.checkState(split.length == 3);
			
			double lon = Double.parseDouble(split[0]);
			double lat = Double.parseDouble(split[1]);
			double vs30 = Double.parseDouble(split[2]);
			
			Site site = new Site(new Location(lat, lon));
			
			Vs30_Param vs30Param = new Vs30_Param();
			vs30Param.setValue(vs30);
			site.addParameter(vs30Param);
			
			Vs30_TypeParam vs30TypeParam = new Vs30_TypeParam();
			vs30TypeParam.setValue(Vs30_TypeParam.VS30_TYPE_INFERRED);
			site.addParameter(vs30TypeParam);
			
			DepthTo1pt0kmPerSecParam z10Param = new DepthTo1pt0kmPerSecParam(null, true);
			z10Param.setValue(null);
			site.addParameter(z10Param);
			
			DepthTo2pt5kmPerSecParam z25Param = new DepthTo2pt5kmPerSecParam(null, 0d, 1000000, true);
			z25Param.setValue(null);
			site.addParameter(z25Param);
			
			sites.add(site);
		}
		
		read.close();
		
		return sites;
	}

}
