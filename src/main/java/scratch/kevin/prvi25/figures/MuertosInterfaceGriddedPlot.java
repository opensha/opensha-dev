package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class MuertosInterfaceGriddedPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "sub_grid");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution subSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		levels.addAll(PRVI25_LogicTree.levelsSubduction);
		levels.addAll(PRVI25_LogicTree.levelsSubductionGridded);
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
		
		branch.setValue(PRVI25_SubductionScalingRelationships.AVERAGE);
		branch.setValue(PRVI25_SubductionMuertosSeismicityRate.PREFFERRED);
		branch.setValue(PRVI25_SubductionCaribbeanSeismicityRate.PREFFERRED);
		branch.setValue(PRVI25_SeisSmoothingAlgorithms.AVERAGE);
		branch.setValue(PRVI25_DeclusteringAlgorithms.AVERAGE);
		
		Location targetLoc = new Location(17.96, -66.30);
//		Location targetLoc = new Location(17.8, -66.30);
		double targetMag = 7.45;
		
		Range latRange = new Range(17.2, 18.2);
		Range depthRange = new Range(-2d, 30d);
		
		double latLenKm = LocationUtils.horzDistanceFast(new Location(latRange.getLowerBound(), 0d),
				new Location(latRange.getUpperBound(), 0d));
		
		PRVI25_GridSourceBuilder.INTERFACE_USE_SECT_PROPERTIES = true;
		GridSourceList sectSurfaceGridList = PRVI25_GridSourceBuilder.buildInterfaceGridSourceList(
				subSol, branch, PRVI25_SeismicityRegions.MUE_INTERFACE);
		File sectFile = new File("/tmp/sect_surf_mue_interface.zip");
		if (!sectFile.exists()) {
			subSol.setGridSourceProvider(sectSurfaceGridList);
			subSol.write(sectFile);
		}
		PRVI25_GridSourceBuilder.INTERFACE_USE_SECT_PROPERTIES = false;
		GridSourceList slabSurfaceGridList = PRVI25_GridSourceBuilder.buildInterfaceGridSourceList(
				subSol, branch, PRVI25_SeismicityRegions.MUE_INTERFACE);
		File slabFile = new File("/tmp/slab_surf_mue_interface.zip");
		if (!slabFile.exists()) {
			subSol.setGridSourceProvider(slabSurfaceGridList);
			subSol.write(slabFile);
		}
		
		int closestLocIndex = -1;
		Location mappedLoc = null;
		double closestDist = Double.POSITIVE_INFINITY;
		for (int l=0; l<sectSurfaceGridList.getNumLocations(); l++) {
			Location loc = sectSurfaceGridList.getLocation(l);
			double dist = LocationUtils.cartesianDistanceSq(targetLoc, loc);
//			double dist = LocationUtils.horzDistanceFast(targetLoc, loc);
			if (dist < closestDist) {
				closestDist = dist;
				closestLocIndex = l;
				mappedLoc = loc;
			}
		}
		
		System.out.println("Mapped "+targetLoc+" to closest grid loc: "+mappedLoc);
		
		GriddedRupture slabSurfRup = null;
		for (GriddedRupture rup : slabSurfaceGridList.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, closestLocIndex)) {
			if ((float)rup.properties.magnitude == (float)targetMag) {
				slabSurfRup = rup;
				break;
			}
		}
		Preconditions.checkNotNull(slabSurfRup, "Couldn't find slab surface rupture for M%s", (float)targetMag);
		
		GriddedRupture sectSurfRup = null;
		for (GriddedRupture rup : sectSurfaceGridList.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, closestLocIndex)) {
			if ((float)rup.properties.magnitude == (float)targetMag) {
				sectSurfRup = rup;
				break;
			}
		}
		Preconditions.checkNotNull(sectSurfRup, "Couldn't find sect surface rupture for M%s", (float)targetMag);
		
		System.out.println("Slab rup has dip="+(float)slabSurfRup.properties.dip+", depth range=["
				+(float)slabSurfRup.properties.upperDepth+", "+(float)slabSurfRup.properties.lowerDepth+"]");
		System.out.println("Sect rup has dip="+(float)sectSurfRup.properties.dip+", depth range=["
				+(float)sectSurfRup.properties.upperDepth+", "+(float)sectSurfRup.properties.lowerDepth+"]");
		
		DiscretizedFunc slabDepthProfile = new ArbitrarilyDiscretizedFunc();
		for (int l=0; l<slabSurfaceGridList.getNumLocations(); l++) {
			Location loc = slabSurfaceGridList.getLocation(l);
			if ((float)loc.lon == (float)mappedLoc.lon) {
				GriddedRupture firstRup = slabSurfaceGridList.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, l).get(0);
				double depth = firstRup.properties.getHypocentralDepth();
				slabDepthProfile.set(loc.lat, depth);
			}
		}
		System.out.println("Slab depth profile has "+slabDepthProfile.size()+" points");
		
		FaultSection mappedSect = null;
		EvenlyGriddedSurface mappedSectSurf = null;
		int mappedSectSurfRow = -1;
		int mappedSectSurfCol = -1;
		closestDist = Double.POSITIVE_INFINITY;
		for (FaultSection sect : subSol.getRupSet().getFaultSectionDataList()) {
			if (!sect.getName().contains("Muertos"))
				continue;
			EvenlyGriddedSurface surf = (EvenlyGriddedSurface)sect.getFaultSurface(5d);
			for (int row=0; row<surf.getNumRows(); row++) {
				for (int col=0; col<surf.getNumCols(); col++) {
					Location loc = surf.getLocation(row, col);
					double dist = LocationUtils.cartesianDistanceSq(targetLoc, loc);
					if (dist < closestDist) {
						closestDist = dist;
						mappedSect = sect;
						mappedSectSurf = surf;
						mappedSectSurfRow = row;
						mappedSectSurfCol = col;
					}
				}
			}
		}
		
		DiscretizedFunc sectDepthProvile = new ArbitrarilyDiscretizedFunc();
		
		for (int row = 0; row < mappedSectSurf.getNumRows(); row++) {
			Location loc = mappedSectSurf.getLocation(row, mappedSectSurfCol);
			double depth = mappedSectSurf.get(row, mappedSectSurfCol).getDepth();
			sectDepthProvile.set(loc.lat, depth);
		}
		System.out.println("Sect depth profile has "+sectDepthProvile.size()+" points");
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		slabDepthProfile.setName("Slab2 surface");
		funcs.add(slabDepthProfile);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
		
		DiscretizedFunc slabRupDepthProfile = rupSizeProfile(slabSurfRup, mappedLoc);
		slabRupDepthProfile.setName("Slab2 surface rupture");
		funcs.add(slabRupDepthProfile);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
		
		sectDepthProvile.setName("Section surface");
		funcs.add(sectDepthProvile);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		
		DiscretizedFunc sectRupDepthProfile = rupSizeProfile(sectSurfRup, mappedLoc);
		sectRupDepthProfile.setName("Section surface rupture");
		funcs.add(sectRupDepthProfile);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Latitude", "Depth (km)");
		plot.setLegendVisible(true);
		plot.setYAxisInverted(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, false, latRange, depthRange);
		
		int width = 800;
		
		double aspectRatio = latLenKm/depthRange.getLength();
		System.out.println("Aspect ratio: "+(float)aspectRatio);
		int height = PlotUtils.calcHeight(gp.getChartPanel(), width, aspectRatio);
		
		PlotUtils.writePlots(outputDir, "muertos_side_profile", gp, width, height, true, true, false);
	}
	
	private static DiscretizedFunc rupSizeProfile(GriddedRupture rup, Location refLoc) {
		EvenlyGriddedSurface surf = GridSourceList.surfBuilderForRup(rup).buildGriddedSurface();
		
		int closestLocCol = -1;
		double minDist = Double.POSITIVE_INFINITY;
		for (int row=0; row<surf.getNumRows(); row++) {
			for (int col = 0; col < surf.getNumCols(); col++) {
				Location loc = surf.getLocation(row, col);
				double dist = LocationUtils.cartesianDistanceSq(refLoc, loc);
				if (dist < minDist) {
					minDist = dist;
					closestLocCol = col;
				}
			}
		}
		
		DiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		for (int r=0; r<surf.getNumRows(); r++) {
			Location loc = surf.get(r, closestLocCol);
			func.set(loc.lat, loc.depth);
		}
		
		return func;
	}

}
