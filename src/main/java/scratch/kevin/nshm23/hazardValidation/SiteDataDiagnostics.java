package scratch.kevin.nshm23.hazardValidation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.data.siteData.CONUS_Downloader;
import org.opensha.commons.data.siteData.CONUS_Versions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CONUS_SiteDataProvider;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.LocalRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.NSHM23_BaseRegion;

import com.google.common.base.Preconditions;

import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.nshmp.shaded.model.NshmpSiteData.Values;

public class SiteDataDiagnostics {
	
	public static void main(String[] args) throws IOException {
		String[] types = {
				SiteData.TYPE_DEPTH_TO_1_0,
				SiteData.TYPE_DEPTH_TO_2_5,
				SiteData.TYPE_SEDIMENT_THICKNESS
		};
		String[] typePrefixes = {
				"z10",
				"z25",
				"zSed"
		};
		double[] cptMaxs = {
				1d,
				5d,
				10d
		};
		
		NSHM23_BaseRegion[] regions = {
				AnalysisRegions.CONUS_EAST,
				AnalysisRegions.CONUS_IMW,
				AnalysisRegions.CONUS_U3_RELM,
				LocalRegions.CONUS_LA_BASIN,
				LocalRegions.CONUS_SF_BAY
				
		};
		
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/site_data_debug");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NshmpHazardModel model = NshmpHazardModel.load(Path.of("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0"));
		org.opensha.nshmp.shaded.model.NshmpSiteData modelSiteData = model.siteData();
		
		for (int t=0; t<types.length; t++) {
			String type = types[t];
			String typePrefix = typePrefixes[t];
			CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, cptMaxs[t]);
			cpt.setNanColor(new Color(0, 0, 0, 0));
			CONUS_SiteDataProvider data = new CONUS_SiteDataProvider(type, CONUS_Versions.NSHM23);
			
			for (NSHM23_BaseRegion analReg : regions) {
				System.out.println("Plotting for "+type+" ("+analReg+")");
				Region reg = analReg.load();
				GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
				
				List<Region> plotRegList = new ArrayList<>();
				List<PlotCurveCharacterstics> outlineChars = new ArrayList<>();
				for (Region plotReg : data.getRegions()) {
					plotRegList.add(plotReg);
					outlineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				}
				for (Region plotReg : data.getRegions()) {
					plotRegList.add(plotReg);
					outlineChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, Color.WHITE));
				}
				mapMaker.plotInsetRegions(plotRegList, outlineChars, null, 0d);
				
				GriddedRegion dataReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
				GriddedGeoDataSet xyz = new GriddedGeoDataSet(dataReg);
				ArrayList<Double> values = data.getValues(dataReg.getNodeList());
				boolean anyGood = false;
				for (int i=0; i<xyz.size(); i++) {
					Double value = values.get(i);
					if (value == null) {
						xyz.set(i, Double.NaN);
					} else {
						anyGood |= value > 0d;
						xyz.set(i, value);
					}
				}
				if (!anyGood)
					continue;
				
				mapMaker.plotXYZData(xyz, cpt, type);
				
				mapMaker.plot(outputDir, typePrefix+"_"+analReg.name()+"_opensha", "OpenSHA");
				
				// load NSHMP-haz
				GriddedGeoDataSet modelXYZ = new GriddedGeoDataSet(dataReg);
				for (int i=0; i<xyz.size(); i++) {
					Location loc = modelXYZ.getLocation(i);
					Values modelVals = modelSiteData.get(org.opensha.nshmp.shaded.geo.NshmpLocation.create(loc.lon, loc.lat));
					OptionalDouble val = getFromModel(modelVals, type);
					if (val.isPresent())
						modelXYZ.set(i, val.getAsDouble());
					else
						modelXYZ.set(i, Double.NaN);
					
					if (!Precision.equals(xyz.get(i), modelXYZ.get(i), 0.001)
							&& (Double.isFinite(xyz.get(i)) || Double.isFinite(modelXYZ.get(i)))) {
						System.out.println("MISMATCH at "+loc);
						System.out.println("\tOpenSHA:\t"+(float)xyz.get(i));
						System.out.println("\tNSHMP-Haz:\t"+(float)modelXYZ.get(i));
					}
				}
				
				mapMaker.plotXYZData(modelXYZ, cpt, type);
				
				mapMaker.plot(outputDir, typePrefix+"_"+analReg.name()+"_nshmp_lib", "NSHMP-Lib");
			}
		}
	}
	
	private static OptionalDouble getFromModel(Values modelVals, String type) {
		return switch (type) {
		case SiteData.TYPE_DEPTH_TO_1_0:
			yield modelVals.z1p0;
		case SiteData.TYPE_DEPTH_TO_2_5:
			yield modelVals.z2p5;
		case SiteData.TYPE_SEDIMENT_THICKNESS:
			yield modelVals.zSed;
		default:
			throw new IllegalArgumentException("Unexpected value: " + type);
		};
	}

}
