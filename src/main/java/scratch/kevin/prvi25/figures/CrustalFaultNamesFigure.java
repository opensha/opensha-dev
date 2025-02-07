package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

class CrustalFaultNamesFigure {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "crustal_fm");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		PRVI25_CrustalFaultModels fm = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1;
		
//		List<? extends FaultSection> sects = fm.getFaultSections();
		List<? extends FaultSection> sects = fm.getDefaultDeformationModel().build(fm);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(SlipRateFigures.CRUSTAL_FAULT_MAP_REG, sects);
		mapMaker.setWriteGeoJSON(false);
		
		addStandardFaultLabels(mapMaker, sects);
		
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
		
		mapMaker.plot(outputDir, "crustal_sect_names", " ");
	}
	
	static void addStandardFaultLabels(GeographicMapMaker mapMaker, List<? extends FaultSection> sects) {
		List<XYTextAnnotation> anns = new ArrayList<>();
		List<LocationList> arrows = new ArrayList<>();
		Font labelFont = new Font(Font.SANS_SERIF, Font.BOLD, 14);
//		buildLabelsAndArrows(sects, "Bunce", false, new Location(20, -66),
//				labelFont, TextAnchor.BASELINE_CENTER, 0d, anns, arrows);
		buildLabelsAndArrows(sects, List.of("Bunce 1", "Bunce 2", "Bunce 3", "Bunce 4", "Bunce 5"), "Bunce 1-5", false, new Location(19.75, -67.7),
				labelFont, TextAnchor.BASELINE_CENTER, -8d, anns, arrows);
		buildLabelsAndArrows(sects, List.of("Bunce 6", "Bunce 7", "Bunce 8"), "Bunce 6-8", false, new Location(19.45, -64),
				labelFont, TextAnchor.TOP_CENTER, 8d, anns, arrows);
		buildLabelsAndArrows(sects, "Anegada Passage", false, new Location(18.7, -64.7),
				labelFont, TextAnchor.BASELINE_CENTER, -27, anns, arrows);
		buildLabelsAndArrows(sects, "Main Ridge", false, new Location(19.2, -66.3),
				labelFont, TextAnchor.TOP_CENTER, 0d, anns, arrows);
		buildLabelsAndArrows(sects, "Septentrional", false, new Location(19.0, -69.5),
				labelFont, TextAnchor.TOP_CENTER, 9d, anns, arrows);
		buildLabelsAndArrows(sects, "South Lajas", false, new Location(17.75, -67.35),
				labelFont, TextAnchor.CENTER, 0d, anns, arrows);
		buildLabelsAndArrows(sects, "Mona Passage", false, new Location(17.5, -68.25),
				labelFont, TextAnchor.CENTER, 0d, anns, arrows);
		buildLabelsAndArrows(sects, "Great Southern Puerto Rico", false, new Location(17.25, -66.25),
				labelFont, TextAnchor.CENTER, 0d, anns, arrows);
		buildLabelsAndArrows(sects, "Salinas", false, new Location(18.15, -66.35),
				labelFont, TextAnchor.BASELINE_CENTER, 0d, anns, arrows);
		
		mapMaker.addAnnotations(anns);
		mapMaker.plotArrows(arrows, 6d, new Color(0, 0, 0, 180), 1f);
		mapMaker.setFillArrowheads(true);
	}
	
	static void buildLabelsAndArrows(List<? extends FaultSection> sects, String name, boolean stitch,
			Location textLoc, Font font, TextAnchor textAnchor, double rotationAngle,
			List<? super XYTextAnnotation> anns, List<LocationList> arrows) {
		buildLabelsAndArrows(sects, List.of(name), name, stitch, textLoc, font, textAnchor, rotationAngle, anns, arrows);
	}
	
	static void buildLabelsAndArrows(List<? extends FaultSection> sects, List<String> searchNames, String plotName, boolean stitch,
			Location textLoc, Font font, TextAnchor textAnchor, double rotationAngle,
			List<? super XYTextAnnotation> anns, List<LocationList> arrows) {
		Map<Integer, List<FaultSection>> byParent = stitch ? null : new HashMap<>();
		List<FaultSection> allCombined = stitch ? new ArrayList<>() : null;
		for (FaultSection sect : sects) {
			for (String name : searchNames) {
				if (sect.getName().contains(name)) {
					if (stitch) {
						allCombined.add(sect);
					} else {
						int parent = sect.getParentSectionId();
						if (parent == -1)
							// already parents
							parent = sect.getSectionId();
						if (!byParent.containsKey(parent))
							byParent.put(parent, new ArrayList<>());
						byParent.get(parent).add(sect);
					}
					break;
				}
			}
		}
		XYTextAnnotation ann = new XYTextAnnotation(plotName, textLoc.getLongitude(), textLoc.getLatitude());
		ann.setTextAnchor(textAnchor);
		ann.setFont(font);
		if (rotationAngle != 0d)
			ann.setRotationAngle(GeographicMapMaker.getRotationAngleCorrectedForAspectRatio(Math.toRadians(rotationAngle), textLoc));
		anns.add(ann);
		
		Collection<List<FaultSection>> bundles = stitch ? List.of(allCombined) : byParent.values();
		
		for (List<FaultSection> bundle : bundles) {
			Location middle = findMiddleLoc(bundle);
			
			// don't draw the arrow all the way, leave a buffer on each end
			double dist = LocationUtils.horzDistanceFast(middle, textLoc);
			double az = LocationUtils.azimuthRad(textLoc, middle);
			
			double buffer;
			if (dist > 40)
				buffer = 10d;
			else if (dist > 20)
				buffer = 5;
			else if (dist > 10)
				buffer = 2;
			else
				buffer = 0.15*dist;
			
			Location start = LocationUtils.location(textLoc, az, buffer);
			Location end = LocationUtils.location(textLoc, az, dist-buffer);
			
			arrows.add(LocationList.of(start, end));
		}
	}
	
	private static Location findMiddleLoc(List<? extends FaultSection> sects) {
		FaultTrace combTrace = new FaultTrace(null);
		
		for (FaultSection sect : sects)
			combTrace.addAll(sect.getFaultTrace());
		
		combTrace = FaultUtils.resampleTrace(combTrace, 2);
		return combTrace.get(1);
	}

}
