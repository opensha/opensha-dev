package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCouplingModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;
import com.google.common.collect.MapMaker;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class PRVI_SubductionSubSectPlots {
	
	public static Region plotReg = new Region(new Location(16.5, -71.5), new Location(21, -60));
	private static Location regCenter = new Location(0.5*(plotReg.getMinLat() + plotReg.getMaxLat()),
			0.5*(plotReg.getMinLon() + plotReg.getMaxLon()));

	public static void main(String[] args) throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = PRVI25_LogicTreeBranch.DEFAULT_SUBDUCTION_INTERFACE;
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.AVERAGE;
		String scaleLabel = "average scaling";
		branch = branch.copy();
		branch.setValue(scale);
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		
//		File outputDir = new File("/home/kevin/Documents/papers/2024_PRVI_Subduction/figures/fault_model");
		File outputDir = new File(FIGURES_DIR, "sub_fm");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		DecimalFormat magDF = new DecimalFormat("0.0");
		
		CPT magCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(7.5d, 9d);
		CPT minMagCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(7.5d, 8.400001d);
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 8d);
		CPT slipUncertCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 2d);
		CPT rakeCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, 90d);
		
		PlotUtils.writeScaleLegendOnly(outputDir, "slip_cpt",
				GeographicMapMaker.buildCPTLegend(slipCPT, "Slip Deficit Rate (mm/yr)"),
				GeographicMapMaker.PLOT_WIDTH_DEFAULT, true, true);
		for (PRVI25_SubductionFaultModels fm : PRVI25_SubductionFaultModels.values()) {
			if (fm.getNodeWeight(branch) == 0d)
				continue;
			
			String fmName = fm.getShortName()+" Fault Model";
			
			branch.setValue(fm);
			branch.requireValue(PRVI25_SubductionDeformationModels.class).build(fm);
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
			
			List<? extends FaultSection> sects = rupSet.getFaultSectionDataList();
			GeographicMapMaker mapMaker = new GeographicMapMaker(plotReg);
			mapMaker.setWriteGeoJSON(false);
			mapMaker.setFaultSections(sects);
			mapMaker.setFillSurfaces(true);
			
			List<Double> minMags = new ArrayList<>();
			List<Double> maxMags = new ArrayList<>();
			for (int s=0; s<sects.size(); s++) {
				minMags.add(rupSet.getMinMagForSection(s));
				maxMags.add(rupSet.getMaxMagForSection(s));
			}
			
			double minMin = minMags.stream().mapToDouble(D->D).min().getAsDouble();
			double maxMin = minMags.stream().mapToDouble(D->D).max().getAsDouble();
			double minMax = maxMags.stream().mapToDouble(D->D).min().getAsDouble();
			double maxMax = maxMags.stream().mapToDouble(D->D).max().getAsDouble();
			
			Range xRange = mapMaker.getXRange();
			Range yRange = mapMaker.getYRange();
			
			double annX = xRange.getLowerBound() + 0.97*xRange.getLength();
			double annY = yRange.getLowerBound() + 0.95*yRange.getLength();

			mapMaker.plotSectScalars(minMags, magCPT, "Minimum Magnitude ("+scaleLabel+")");
			XYTextAnnotation rangeAnn = new XYTextAnnotation("["+magDF.format(minMin)+", "+magDF.format(maxMin)+"]", annX, annY);
			rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			rangeAnn.setFont(font);
			mapMaker.addAnnotation(rangeAnn);
			mapMaker.plot(outputDir, "subduction_min_mag_"+fm.getFilePrefix(), fmName);
			
			List<XYTextAnnotation> labelAnns = getLabelAnns(mapMaker);

			mapMaker.plotSectScalars(minMags, minMagCPT, "Minimum Magnitude ("+scaleLabel+")");
			mapMaker.setAnnotations(labelAnns);
			mapMaker.plot(outputDir, "subduction_min_mag_"+fm.getFilePrefix()+"_names", " ");
			
			// now highlight the largest and smallest  supra-seis ruptures
			int smallestOverallRupIndex = -1;
			double smallestOverallMag = Double.POSITIVE_INFINITY;
			int largestSmallestSectRupIndex = -1;
			double largestSmallestMag = 0d;
			for (int s=0; s<sects.size(); s++) {
				double mySmallest = Double.POSITIVE_INFINITY;
				int mySmallestIndex = -1;
				for (int rupIndex : rupSet.getRupturesForSection(s)) {
					double mag = rupSet.getMagForRup(rupIndex);
					if (mag < mySmallest) {
						mySmallest = mag;
						mySmallestIndex = rupIndex;
					}
				}
				if (mySmallest < smallestOverallMag) {
					smallestOverallMag = mySmallest;
					smallestOverallRupIndex = mySmallestIndex;
				}
				if (mySmallest > largestSmallestMag) {
					largestSmallestMag = mySmallest;
					largestSmallestSectRupIndex = mySmallestIndex;
				}
			}
			
			Region smallRupRegion = getRupRegion(rupSet, smallestOverallRupIndex);
			Region largeRupRegion = getRupRegion(rupSet, largestSmallestSectRupIndex);
//			List<Region> rupRegions = List.of(smallRupRegion, largeRupRegion);
			List<Region> rupRegions = List.of(smallRupRegion, largeRupRegion, smallRupRegion, largeRupRegion);
			List<PlotCurveCharacterstics> rupChars = List.of(
//					new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_blue),
//					new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_orange));
					new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.WHITE),
					new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.WHITE),
					new PlotCurveCharacterstics(PlotLineType.DASHED, 5f, Color.BLACK),
					new PlotCurveCharacterstics(PlotLineType.DASHED, 5f, Color.BLACK));
			mapMaker.plotInsetRegions(rupRegions, rupChars, null, 0f);
			mapMaker.plot(outputDir, "subduction_min_mag_"+fm.getFilePrefix()+"_names_outlines", " ");
			mapMaker.clearInsetRegions();
			
			mapMaker.plotSectScalars(maxMags, magCPT, "Maximum Magnitude ("+scaleLabel+")");
			mapMaker.clearAnnotations();
			rangeAnn = new XYTextAnnotation("["+magDF.format(minMax)+", "+magDF.format(maxMax)+"]", annX, annY);
			rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			rangeAnn.setFont(font);
			mapMaker.addAnnotation(rangeAnn);
			mapMaker.plot(outputDir, "subduction_max_mag_"+fm.getFilePrefix(), fmName);
			
			mapMaker.setAnnotations(labelAnns);
			mapMaker.plot(outputDir, "subduction_max_mag_"+fm.getFilePrefix()+"_names", fmName);
			
			mapMaker.clearAnnotations();
			
			for (PRVI25_SubductionDeformationModels dm : PRVI25_SubductionDeformationModels.values()) {
				branch.setValue(dm);
				if (dm.getNodeWeight(branch) == 0d)
					continue;
				
				for (PRVI25_SubductionCouplingModels coupling : PRVI25_SubductionCouplingModels.values()) {
					branch.setValue(coupling);
					if (coupling.getNodeWeight(branch) == 0d)
						continue;
					
					String dmName = dm.getName()+", "+coupling.getShortName()+" Coupling";
					String dmPrefix = fm.getFilePrefix()+"_"+dm.getFilePrefix()+"_"+coupling.getFilePrefix();
					sects = dm.build(fm, branch);
					
					mapMaker.setFaultSections(sects);
					List<Double> slips = new ArrayList<>(sects.size());
					List<Double> slipUncerts = new ArrayList<>(sects.size());
					List<Double> rakes = new ArrayList<>(sects.size());
					for (FaultSection sect : sects) {
						slips.add(sect.getOrigAveSlipRate());
						slipUncerts.add(sect.getOrigSlipRateStdDev());
						rakes.add(sect.getAveRake());
					}
					
//					mapMaker.plotSectScalars(slipUncerts, slipUncertCPT, dm.getShortName()+" Slip Rate Uncertainty (mm/yr)");
//					mapMaker.plot(outputDir, "subduction_slip_uncert_"+fm.getFilePrefix()+"_"+dm.getFilePrefix(), fmName+", "+dmName);
					
					// draw rake lines
					int rakeMod = 3;
					List<LocationList> rakeArrows = new ArrayList<>();
					for (FaultSection sect : sects)
						if (sect.getSectionId() % rakeMod == 0)
							rakeArrows.addAll(buildRakeArrows(sect, slipCPT.getMaxValue()));
					mapMaker.plotArrows(rakeArrows, 20d, Colors.tab_red.darker(), 2f);
					mapMaker.setFillArrowheads(true, Color.WHITE, 0.5f);
//					mapMaker.plotLines(rakeArrows, Color.BLACK, 2f);
					
					mapMaker.plotSectScalars(slips, slipCPT, dm.getShortName()+" Slip Rate (mm/yr)");
					mapMaker.plot(outputDir, "subduction_slip_"+dmPrefix, fmName+", "+dmName);
					
					mapMaker.plotSectScalars(slips, slipCPT, null);
					mapMaker.plot(outputDir, "subduction_slip_"+dmPrefix+"_no_cpt", fmName+", "+dmName);
					
					mapMaker.plotSectScalars(rakes, rakeCPT, dm.getShortName()+" Rake");
					String prefix = "subduction_rake_"+fm.getFilePrefix()+"_"+dm.getFilePrefix();;
					System.out.println("Plotting "+prefix);
					mapMaker.plot(outputDir, prefix, fmName+", "+dmName);
					
					mapMaker.clearLines();
					mapMaker.clearArrows();
				}
			}
		}
	}
	
	private static Region getRupRegion(FaultSystemRupSet rupSet, int rupIndex) {
//		RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 5d);
//		LocationList border = surf.getEvenlyDiscritizedPerimeter();
//		if (!border.first().equals(border.last()))
//			border.add(border.first());
//		System.out.println(border);
		LocationList border = new LocationList();
		List<FaultSection> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
		for (FaultSection sect : sects) {
			for (Location loc : sect.getFaultTrace())
				if (border.isEmpty() || !LocationUtils.areSimilar(loc, border.last()))
					border.add(loc);
		}
		for (int i=sects.size(); --i>=0;) {
			FaultTrace trace = sects.get(i).getLowerFaultTrace();
			for (int j=trace.size(); --j>=0;) {
				Location loc = trace.get(j);
				if (border.isEmpty() || !LocationUtils.areSimilar(loc, border.last()))
					border.add(loc);
			}
		}
		border.add(border.first());
		return new Region(border, BorderType.MERCATOR_LINEAR);
	}

	public static List<XYTextAnnotation> getLabelAnns(GeographicMapMaker mapMaker) {
		List<XYTextAnnotation> labelAnns = new ArrayList<>();
		Font subInterfaceFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		XYTextAnnotation hspAnn = new XYTextAnnotation("Northern Hispaniola", -69.8, 20.05);
		hspAnn.setRotationAngle(mapMaker.getRotationAngleCorrectedForAspectRatio(Math.toRadians(18)));
		hspAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
		hspAnn.setFont(subInterfaceFont);
		labelAnns.add(hspAnn);
		
		XYTextAnnotation prviAnn = new XYTextAnnotation("Puerto Rico-Virgin Islands", -65.5, 20.05);
//			prviAnn.setRotationAngle(mapMaker.getRotationAngleCorrectedForAspectRatio(Math.toRadians(20)));
		prviAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
		prviAnn.setFont(subInterfaceFont);
		labelAnns.add(prviAnn);
		
		XYTextAnnotation laAnn = new XYTextAnnotation("Lesser Antilles", -62, 19.55);
		laAnn.setRotationAngle(mapMaker.getRotationAngleCorrectedForAspectRatio(Math.toRadians(23)));
		laAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
		laAnn.setFont(subInterfaceFont);
		labelAnns.add(laAnn);
		
		Font interfaceFont = new Font(Font.SANS_SERIF, Font.BOLD, 26);
		XYTextAnnotation carAnn = new XYTextAnnotation("Caribbean Trench", -65.5, 20.4);
		carAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
		carAnn.setFont(interfaceFont);
		labelAnns.add(carAnn);
		XYTextAnnotation mueAnn = new XYTextAnnotation("Muertos Trough", -68, 17.2);
		mueAnn.setTextAnchor(TextAnchor.TOP_CENTER);
		mueAnn.setFont(interfaceFont);
		labelAnns.add(mueAnn);
		return labelAnns;
	}
	
	public static List<LocationList> buildRakeArrows(FaultSection sect, double slipForMaxLen) {
		return buildRakeArrows(sect, regCenter, 200d, slipForMaxLen);
	}
	
	public static List<LocationList> buildRakeArrows(FaultSection sect, Location refLoc, double maxLen, double slipForMaxLen) {
		List<LocationList> rakeArrows = new ArrayList<>(2);
		Vector3D vect = FaultSystemLineIntegralCalculator.calcHangingWallSlipVector(sect);
		double len = maxLen * sect.getOrigAveSlipRate()/slipForMaxLen;
		System.out.println("Rake arrow for "+sect.getSectionName()+" "+vect+", plotLen="+(float)len+", slipRate="
				+(float)sect.getOrigAveSlipRate()+", rake="+(float)sect.getAveRake()+", dip="+(float)sect.getAveDip());
		
//		FaultTrace upperTrace = sect.getFaultTrace();
//		FaultTrace lowerTrace = sect.getLowerFaultTrace();
//		Location center = new Location(0.25*(upperTrace.first().lat + lowerTrace.first().lat + upperTrace.last().lat + lowerTrace.last().lat),
//				0.25*(upperTrace.first().lon + lowerTrace.first().lon + upperTrace.last().lon + lowerTrace.last().lon));
		
		RuptureSurface surf = sect.getFaultSurface(1d);
		LocationList locs = surf.getEvenlyDiscritizedListOfLocsOnSurface();
		double avgLat = locs.stream().mapToDouble(L->L.lat).average().getAsDouble();
		double avgLon = locs.stream().mapToDouble(L->L.lon).average().getAsDouble();
		Location center = new Location(avgLat, avgLon);
		
		double az = FaultSystemLineIntegralCalculator.vectorAzimuth(vect, center);
		Location lineEnd = FaultSystemLineIntegralCalculator.locConstPlotDist(center, refLoc,
				az, len);
		rakeArrows.add(LocationList.of(center, lineEnd));
//		Location arrowStart = FaultSystemLineIntegralCalculator.locConstPlotDist(lineEnd, center, az - 135, 0.2*len);
//		Location arrowEnd = FaultSystemLineIntegralCalculator.locConstPlotDist(lineEnd, center, az + 135, 0.2*len);
//		rakeArrows.add(LocationList.of(arrowStart, lineEnd, arrowEnd));
		return rakeArrows;
	}

}
