package scratch.kevin.nshm27;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM27_GridSourceBuilder;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import net.mahdilamb.colormap.Colors;

public class BathymetryConvert {

	public static void main(String[] args) throws IOException {
		// data from https://www.ncei.noaa.gov/maps/grid-extract/
		// ETOPO 2022 layer (ice surface)
		// converted from geotiff to xyz via gdal_translate -of XYZ exportImage gnmi_bathy.xyz
		
		NSHM27_SeismicityRegions seisReg = NSHM27_SeismicityRegions.GNMI;
		File inFile = new File("/home/kevin/Downloads/gnmi_bathy.xyz");
		
//		NSHM27_SeismicityRegions seisReg = NSHM27_SeismicityRegions.AMSAM;
//		File inFile = new File("/home/kevin/Downloads/amsam_bathy.xyz");
		
		GriddedRegion gridReg = NSHM27_GridSourceBuilder.initGridReg(seisReg);
		GriddedGeoDataSet maxDists = new GriddedGeoDataSet(gridReg);
		GriddedGeoDataSet depths = new GriddedGeoDataSet(gridReg);
		for (int i=0; i<maxDists.size(); i++) {
			maxDists.set(i, Double.POSITIVE_INFINITY);
			depths.set(i, Double.NaN);
		}
		
		try (BufferedReader br = new BufferedReader(new FileReader(inFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				String[] split = line.split(" ");
				Preconditions.checkState(split.length == 3);
				double lon = Double.parseDouble(split[0]);
				if (lon < 0d)
					lon += 360d;
				double lat = Double.parseDouble(split[1]);
				double elevation = Double.parseDouble(split[2]);
				double depth = -elevation / 1000d;
				Location loc = new Location(lat, lon);
				int index = gridReg.indexForLocation(loc);
				if (index >= 0) {
					Location gridLoc = gridReg.getLocation(index);
					double dist = LocationUtils.cartesianDistanceSq(loc, gridLoc);
					if (dist < maxDists.get(index)) {
						maxDists.set(index, dist);
						depths.set(index, depth);
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		GriddedGeoDataSet.writeXYZFile(depths, new File("/tmp/"+seisReg.name()+"_depths.xyz"));
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		CPT cpt = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().reverse().rescale(0d, 10d);
		cpt.setNanColor(Colors.tab_orange);
		mapMaker.plotXYZData(depths, cpt, "NOAA ETOPO (2022) bathymetric depth (km)");
		mapMaker.plot(new File("/tmp"), seisReg.name()+"_depths", " ");
	}

}
