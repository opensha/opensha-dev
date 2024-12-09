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
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator;

public class PRVI_SubductionSubSectPlots {
	
	public static Region plotReg = new Region(new Location(16.5, -71.5), new Location(21, -60));
	private static Location regCenter = new Location(0.5*(plotReg.getMinLat() + plotReg.getMaxLat()),
			0.5*(plotReg.getMinLon() + plotReg.getMaxLon()));

	public static void main(String[] args) throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = PRVI25_LogicTreeBranch.DEFAULT_SUBDUCTION_INTERFACE;
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.LOGA_C4p0;
		branch = branch.copy();
		branch.setValue(scale);
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		
//		File outputDir = new File("/home/kevin/Documents/papers/2024_PRVI_Subduction/figures/fault_model");
		File outputDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures/sub_fm");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		DecimalFormat magDF = new DecimalFormat("0.00");
		
		CPT magCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(7.5d, 9d);
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 5d);
		CPT slipUncertCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 2d);
		CPT rakeCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, 90d);
		
		PlotUtils.writeScaleLegendOnly(outputDir, "slip_cpt",
				GeographicMapMaker.buildCPTLegend(slipCPT, "Slip Rate (mm/yr)"),
				GeographicMapMaker.PLOT_WIDTH_DEFAULT, true, true);
		for (PRVI25_SubductionFaultModels fm : PRVI25_SubductionFaultModels.values()) {
			if (fm.getNodeWeight(branch) == 0d)
				continue;
			
			String fmName = fm.getShortName()+" Fault Model";
			
			branch = branch.copy();
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

			mapMaker.plotSectScalars(minMags, magCPT, "Minimum Magnitude ("+scale.getShortName()+")");
			XYTextAnnotation rangeAnn = new XYTextAnnotation("["+magDF.format(minMin)+", "+magDF.format(maxMin)+"]", annX, annY);
			rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			rangeAnn.setFont(font);
			mapMaker.addAnnotation(rangeAnn);
			mapMaker.plot(outputDir, "subduction_min_mag_"+fm.getFilePrefix(), fmName);
			
			mapMaker.plotSectScalars(maxMags, magCPT, "Maximum Magnitude ("+scale.getShortName()+")");
			mapMaker.clearAnnotations();
			rangeAnn = new XYTextAnnotation("["+magDF.format(minMax)+", "+magDF.format(maxMax)+"]", annX, annY);
			rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			rangeAnn.setFont(font);
			mapMaker.addAnnotation(rangeAnn);
			mapMaker.plot(outputDir, "subduction_max_mag_"+fm.getFilePrefix(), fmName);
			
			mapMaker.clearAnnotations();
			
			for (PRVI25_SubductionDeformationModels dm : PRVI25_SubductionDeformationModels.values()) {
				if (dm.getNodeWeight(branch) == 0d)
					continue;
				String dmName = dm.getName();
				sects = dm.build(fm);
				
				mapMaker.setFaultSections(sects);
				List<Double> slips = new ArrayList<>(sects.size());
				List<Double> slipUncerts = new ArrayList<>(sects.size());
				List<Double> rakes = new ArrayList<>(sects.size());
				for (FaultSection sect : sects) {
					slips.add(sect.getOrigAveSlipRate());
					slipUncerts.add(sect.getOrigSlipRateStdDev());
					rakes.add(sect.getAveRake());
				}
				
//				mapMaker.plotSectScalars(slipUncerts, slipUncertCPT, dm.getShortName()+" Slip Rate Uncertainty (mm/yr)");
//				mapMaker.plot(outputDir, "subduction_slip_uncert_"+fm.getFilePrefix()+"_"+dm.getFilePrefix(), fmName+", "+dmName);
				
				// draw rake lines
				int rakeMod = 3;
				List<LocationList> rakeArrows = new ArrayList<>();
				for (FaultSection sect : sects)
					if (sect.getSectionId() % rakeMod == 0)
						rakeArrows.addAll(buildRakeArrows(sect, slipCPT.getMaxValue()));
				mapMaker.plotArrows(rakeArrows, 20d, Color.BLACK, 2f);
				mapMaker.setFillArrowheads(true);
				mapMaker.plotLines(rakeArrows, Color.BLACK, 2f);
				
				mapMaker.plotSectScalars(slips, slipCPT, dm.getShortName()+" Slip Rate (mm/yr)");
				mapMaker.plot(outputDir, "subduction_slip_"+fm.getFilePrefix()+"_"+dm.getFilePrefix(), fmName+", "+dmName);
				
				mapMaker.plotSectScalars(slips, slipCPT, null);
				mapMaker.plot(outputDir, "subduction_slip_"+fm.getFilePrefix()+"_"+dm.getFilePrefix()+"_no_cpt", fmName+", "+dmName);
				
				mapMaker.plotSectScalars(rakes, rakeCPT, dm.getShortName()+" Rake");
				String prefix = "subduction_rake_"+fm.getFilePrefix()+"_"+dm.getFilePrefix();;
				System.out.println("Plotting "+prefix);
				mapMaker.plot(outputDir, prefix, fmName+", "+dmName);
				
				mapMaker.clearLines();
				mapMaker.clearArrows();
			}
		}
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
