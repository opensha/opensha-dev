package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class TectonicSettingPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = FIGURES_DIR;
		Preconditions.checkState(outputDir.exists() || outputDir.mkdirs());
		
		Region modelReg = PRVI25_RegionLoader.loadPRVI_ModelBroad();
		
		double bufferKM = 30d;
		Region buffered = new Region(
				LocationUtils.location(new Location(modelReg.getMinLat(), modelReg.getMinLon()), 5d*Math.PI/4d, bufferKM),
				LocationUtils.location(new Location(modelReg.getMaxLat(), modelReg.getMaxLon()), Math.PI/4d, bufferKM));
		GeographicMapMaker mapMaker = new GeographicMapMaker(buffered);
		mapMaker.setWriteGeoJSON(false);
		
		List<String> legendItems = new ArrayList<>();
		List<PlotCurveCharacterstics> legendChars = new ArrayList<>();
		
		List<FaultSection> allSects = new ArrayList<>();
		List<PlotCurveCharacterstics> allSectChars = new ArrayList<>();
		PlotCurveCharacterstics crustalChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY);
		PlotCurveCharacterstics subChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, Color.DARK_GRAY);
		allSects.addAll(PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1.getFaultSections());
		while (allSectChars.size() < allSects.size())
			allSectChars.add(crustalChar);
		allSects.addAll(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE.getFaultSections());
		while (allSectChars.size() < allSects.size())
			allSectChars.add(subChar);
		mapMaker.setFaultSections(allSects);
		mapMaker.setSectOutlineChar(null);
		mapMaker.setPlotProxySectPolys(false);
		mapMaker.plotSectChars(allSectChars);
		legendItems.add("Crustal traces");
		legendChars.add(crustalChar);
		legendItems.add("Inteface traces");
		legendChars.add(subChar);
		
		List<Region> regions = new ArrayList<>();
		List<PlotCurveCharacterstics> outlines = new ArrayList<>();
		List<Color> fills = new ArrayList<>();
		double fillTrans = 0.1;
		
		
		regions.add(modelReg);
		outlines.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 8f, Color.BLACK));
		fills.add(null);
//		legendItems.add("Model Bounds");
//		legendChars.add(outlines.get(outlines.size()-1));
		Font modelRegFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		XYTextAnnotation modelRegAnn = new XYTextAnnotation(" Model bounds", modelReg.getMinLon(), modelReg.getMaxLat());
		modelRegAnn.setTextAnchor(TextAnchor.TOP_LEFT);
		modelRegAnn.setFont(modelRegFont);
		mapMaker.addAnnotation(modelRegAnn);
		
		Region mapRegion = PRVI25_RegionLoader.loadPRVI_MapExtents();
		regions.add(mapRegion);
		outlines.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		fills.add(null);
//		legendItems.add("Map Region");
//		legendChars.add(outlines.get(outlines.size()-1));
		Font mapRegFont = new Font(Font.SANS_SERIF, Font.PLAIN, 20);
		XYTextAnnotation mapRegAnn = new XYTextAnnotation("Map region", mapRegion.getMaxLon(), mapRegion.getMaxLat());
		mapRegAnn.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
		mapRegAnn.setFont(mapRegFont);
		mapMaker.addAnnotation(mapRegAnn);
		
		PlotCurveCharacterstics interfaceChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue);
		PlotCurveCharacterstics slabChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green);
		
		regions.add(PRVI25_SeismicityRegions.CAR_INTRASLAB.load());
		outlines.add(slabChar);
		fills.add(null);
		
		regions.add(PRVI25_SeismicityRegions.MUE_INTRASLAB.load());
		outlines.add(slabChar);
		fills.add(null);
		legendItems.add("Intraslab regions");
		legendChars.add(slabChar);
		
		regions.add(PRVI25_SeismicityRegions.CAR_INTERFACE.load());
		outlines.add(interfaceChar);
		fills.add(null);
		legendItems.add("Interface regions");
		legendChars.add(interfaceChar);
		
		regions.add(PRVI25_SeismicityRegions.MUE_INTERFACE.load());
		outlines.add(interfaceChar);
		fills.add(null);
		
		Font interfaceFont = new Font(Font.SANS_SERIF, Font.BOLD, 22);
		XYTextAnnotation carAnn = new XYTextAnnotation("Caribbean Trench", -65.5, 20.05);
		carAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
		carAnn.setFont(interfaceFont);
		carAnn.setPaint(subChar.getColor());
		mapMaker.addAnnotation(carAnn);
		XYTextAnnotation mueAnn = new XYTextAnnotation("Muertos Trough", -68, 17.2);
		mueAnn.setTextAnchor(TextAnchor.TOP_CENTER);
		mueAnn.setFont(interfaceFont);
		mueAnn.setPaint(subChar.getColor());
		mapMaker.addAnnotation(mueAnn);
		
		Location arrowStart = LocationUtils.location(new Location(modelReg.getMaxLat(), modelReg.getMaxLon()), 5d*Math.PI/4d, bufferKM);
		Location refLoc = new Location(20.5, -62.5);
		List<? extends FaultSection> subSects = PRVI25_SubductionDeformationModels.FULL.build(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE);
		Vector3D slipVect = null;
		double closestDist = Double.POSITIVE_INFINITY;
		for (FaultSection sect : subSects) {
			if (sect.getName().contains("Muertos"))
				continue;
			Vector3D sectVector = FaultSystemLineIntegralCalculator.calcHangingWallSlipVector(sect);
//			if (slipVect == null)
//				slipVect = sectVector;
//			else
//				slipVect = slipVect.add(sectVector);
			Location start = sect.getFaultTrace().first();
			Location end = sect.getFaultTrace().last();
			Location middleLoc = new Location(0.5*(start.lat + end.lat), 0.5*(start.lon + end.lon));
			double dist = LocationUtils.horzDistanceFast(middleLoc, refLoc);
			if (dist < closestDist) {
				slipVect = sectVector;
				closestDist = dist;
			}
		}
		double slipAzimuth = FaultSystemLineIntegralCalculator.vectorAzimuth(slipVect, arrowStart);
		System.out.println("CAR slip azimuth: "+slipAzimuth);
		Location arrowEnd = LocationUtils.location(arrowStart, Math.toRadians(slipAzimuth+180d), 200d);
		LocationList convergenceArrow = LocationList.of(arrowStart, arrowEnd);
		mapMaker.plotArrows(List.of(convergenceArrow), 30d, Color.BLACK, 6f);
		mapMaker.setFillArrowheads(true);
		
		Font convergenceFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		Location middleLoc = new Location(0.5*(arrowEnd.lat + arrowStart.lat), 0.5*(arrowEnd.lon + arrowStart.lon));
		XYTextAnnotation convergenceAnn = new XYTextAnnotation("15-20 mm/yr", middleLoc.lon, middleLoc.lat);
		convergenceAnn.setTextAnchor(TextAnchor.TOP_CENTER);
		convergenceAnn.setRotationAnchor(TextAnchor.TOP_CENTER);
		convergenceAnn.setFont(convergenceFont);
		double convergenceRotationAngle = 0.5*Math.PI + LocationUtils.azimuthRad(arrowEnd, arrowStart);
		convergenceAnn.setRotationAngle(mapMaker.getRotationAngleCorrectedForAspectRatio(convergenceRotationAngle)-0.004*Math.PI);
		mapMaker.addAnnotation(convergenceAnn);
		
		mapMaker.plotInsetRegions(regions, outlines, fills, fillTrans);
		mapMaker.setPlotRegionsAboveFaults(true);
		
		mapMaker.setCustomLegendItems(legendItems, legendChars);
		
		mapMaker.plot(outputDir, "tectonic_setting", " ");
		PlotSpec plot = mapMaker.buildPlot(" ");
		plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT, 0.05, 0.05, 0.55, false);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, false, mapMaker.getXRange(), mapMaker.getYRange());
		
		PlotUtils.setXTick(gp, 1d);
		PlotUtils.setYTick(gp, 1d);
		
		PlotUtils.writePlots(outputDir, "tectonic_setting", gp, 800, true, true, true, false);;
	}

}
