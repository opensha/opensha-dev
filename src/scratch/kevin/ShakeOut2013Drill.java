package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Shaw_2009_ModifiedMagAreaRel;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class ShakeOut2013Drill {
	
	private static Location calcCenter(Location loc1, Location loc2) {
		// vector from loc1 to loc2
		LocationVector vector = LocationUtils.vector(loc1, loc2);
		// scale vector to half
		vector.set(vector.getAzimuth(), vector.getHorzDistance()*0.5, vector.getVertDistance()*0.5);
		return LocationUtils.location(loc1, vector);
	}

	public static void main(String[] args) throws FileNotFoundException, IOException, GMT_MapException {
		Location faultP1 = new Location(34.3163, -117.549);
		Location faultP2 = new Location(33.957096, -116.81772);
		double mag = 6.85;
		double ddw = 13; // km
		double dip = 90;
		double rake = 180;
		double strike = LocationUtils.azimuth(faultP1, faultP2);
		double actualDist = LocationUtils.linearDistance(faultP1, faultP2);
		
		Shaw_2009_ModifiedMagAreaRel shaw09 = new Shaw_2009_ModifiedMagAreaRel();
		double calcArea = shaw09.getMedianArea(mag, rake);
		double calcLen = calcArea / ddw;
		
		Location topCenter = calcCenter(faultP1, faultP2);
		System.out.println("Top Center: "+topCenter);
		System.out.println("Actual length: "+actualDist+" km");
		System.out.println("Calculated length (from Mag, ddw): "+calcLen+" km");
		System.out.println("Strike: "+strike);
		
		File dir = new File("/home/kevin/Documents/2013 ShakeOut");
		File xyzFile = new File(dir, "/shakeout-site-peaks.txt");
		
		//
		ArbDiscrGeoDataSet pgvData = new ArbDiscrGeoDataSet(false);
		ArbDiscrGeoDataSet pgaData = new ArbDiscrGeoDataSet(false);
		ArbDiscrGeoDataSet mmiData = new ArbDiscrGeoDataSet(false);
		for (String line : FileUtils.loadFile(xyzFile.getAbsolutePath())) {
			while (line.contains("  "))
				line = line.replaceAll("  ", " ");
			String[] elems = line.split(" ");
			double lon = Double.parseDouble(elems[0]);
			double lat = Double.parseDouble(elems[1]);
			double pgv = Double.parseDouble(elems[3]);
			double pga = Double.parseDouble(elems[4]) / 980.665; // convert to g
			
			Location loc = new Location(lat, lon);
			pgvData.set(loc, pgv);
			pgaData.set(loc, pga);
			
			double mmi = Wald_MMI_Calc.getMMI(pga, pgv);
			mmiData.set(loc, mmi);
		}
		
		CPT cpt = GMT_CPT_Files.SHAKEMAP.instance();
		makeMap(pgaData, cpt.rescale(0, 0.5), "PGA (g)", dir, "bb_pga");
		makeMap(pgvData, cpt.rescale(0, 50), "PGV (cm/s)", dir, "bb_pgv");
		makeMap(mmiData, cpt.rescale(0, 10), "MMI", dir, "bb_mmi");
	}
	
	private static void makeMap(GeoDataSet griddedData, CPT cpt, String label, File dir, String prefix)
			throws GMT_MapException, IOException {
		Region region = new Region(new Location(35.0, -118.5), new Location(33.2, -116));
//		GriddedRegion reg = new GriddedRegion(region, 0.02, null);
//		double horzSpacing = LocationUtils.horzDistance(reg.getNodeList().get(0), reg.getNodeList().get(1));
//		int indexChangeLat = -1;
//		for (int i=1; i<reg.getNodeList().size(); i++) {
//			Location loc1 = reg.getNodeList().get(i-1);
//			Location loc2 = reg.getNodeList().get(i);
//			if (loc2.getLatitude() != loc1.getLatitude()) {
//				indexChangeLat = i;
//				break;
//			}
//		}
//		double vertSpacing = LocationUtils.horzDistance(reg.getNodeList().get(0), reg.getNodeList().get(indexChangeLat));
//		System.out.println("horz: "+horzSpacing);
//		System.out.println("vert: "+vertSpacing);
//		System.out.println(reg.getNumLocations());
//		System.exit(0);
		double spacing = 0.1;
		
		GMT_Map map = new GMT_Map(region, griddedData, spacing, cpt);
		
		map.setBlackBackground(false);
		map.setRescaleCPT(false);
		map.setCustomScaleMin((double)cpt.getMinValue());
		map.setCustomScaleMax((double)cpt.getMaxValue());
//		map.setCoast(new CoastAttributes(Color.BLACK, 2));
//		map.s
		map.setCustomLabel(label);
//		map.setUseGMTSmoothing(true);
//		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
		map.setUseGMTSmoothing(false);
		map.setTopoResolution(null);
		
		if (!map.isUseGMTSmoothing())
			prefix += "_nosmooth";
		
		FaultBasedMapGen.plotMap(dir, prefix, true, map);
	}

}
