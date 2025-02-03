package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

public class BranchAveragedHazardFigure {

	public static void main(String[] args) throws IOException {
		File baseDir = new File(INV_DIR, COMBINED_DIR.getName()+"-ba_only-wider_region-vs760");
		File outputDir = new File(FIGURES_DIR, "hazard_map_ba");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File mapResourcesFile = new File(baseDir, "results/hazard_0.02deg_grid_seis_INCLUDE/map_1.0s_TEN_IN_50.txt");
		Region reg = PRVI25_RegionLoader.loadPRVI_IntermediateModelMapExtents();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.02, GriddedRegion.ANCHOR_0_0);
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-1.5, 0);
		cpt.setPreferredTickInterval(0.25);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
		
		String label = "Log10 1s SA (g), 2% in 50 years";
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		BufferedReader bRead = new BufferedReader(new FileReader(mapResourcesFile));
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
		bRead.close();
		
		xyz.log10();
		
		mapMaker.plotXYZData(xyz, cpt, label);
		
		mapMaker.plot(outputDir, "hazard_map", " ");
		
		// now annotate cities mentioned
		
		Font font = new Font(Font.SANS_SERIF, Font.PLAIN, 20);
		List<LocationList> arrows = new ArrayList<>();
		
		Location sjLoc = new Location(18.47, -66.12);
		Location sjAnnLoc = new Location(sjLoc.lat+0.25, sjLoc.lon+0.25);
		XYTextAnnotation ann = new XYTextAnnotation("San Juan", sjAnnLoc.getLongitude()+0.1, sjAnnLoc.getLatitude()+0.1);
		ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
		ann.setFont(font);
		mapMaker.addAnnotation(ann);
		arrows.add(LocationList.of(sjAnnLoc, sjLoc));
		
		Location stcLoc = new Location(17.67, -64.92);
		Location stcAnnLoc = new Location(stcLoc.lat-0.25, stcLoc.lon-0.15);
		ann = new XYTextAnnotation("St. Croix", stcAnnLoc.getLongitude()-0.1, stcAnnLoc.getLatitude());
		ann.setTextAnchor(TextAnchor.TOP_CENTER);
		ann.setFont(font);
		mapMaker.addAnnotation(ann);
		arrows.add(LocationList.of(stcAnnLoc, stcLoc));
		
		Location mgzLoc = new Location(18.2, -67.14);
		Location mgzAnnLoc = new Location(mgzLoc.lat-0.25, mgzLoc.lon-0.35);
		ann = new XYTextAnnotation("Mayagüez", mgzAnnLoc.getLongitude()-0.1, mgzAnnLoc.getLatitude());
		ann.setTextAnchor(TextAnchor.TOP_CENTER);
		ann.setFont(font);
		mapMaker.addAnnotation(ann);
		arrows.add(LocationList.of(mgzAnnLoc, mgzLoc));
		
		Location salLoc = new Location(17.96, -66.30);
		Location salAnnLoc = new Location(salLoc.lat-0.25, salLoc.lon-0.08);
		ann = new XYTextAnnotation("Salinas", salAnnLoc.getLongitude()-0.1, salAnnLoc.getLatitude());
		ann.setTextAnchor(TextAnchor.TOP_CENTER);
		ann.setFont(font);
		mapMaker.addAnnotation(ann);
		arrows.add(LocationList.of(salAnnLoc, salLoc));
		
		for (LocationList arrow : arrows) {
			double az = LocationUtils.azimuth(arrow.first(), arrow.last());
			double dist = LocationUtils.horzDistanceFast(arrow.first(), arrow.last());
			arrow.set(1, LocationUtils.location(arrow.first(), Math.toRadians(az), dist*0.90));
		}
		
		mapMaker.plotArrows(arrows, 7d, Color.BLACK, 1f);
		mapMaker.setFillArrowheads(true);
		
		mapMaker.plot(outputDir, "hazard_map_annotated", " ");
	}

}
