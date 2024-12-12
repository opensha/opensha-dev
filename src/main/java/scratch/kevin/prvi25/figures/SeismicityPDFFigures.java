package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.netlib.util.intW;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

class SeismicityPDFFigures {

	public static void main(String[] args) throws IOException {
		File crustalOutputDir = new File(FIGURES_DIR, "crustal_grid");
		Preconditions.checkState(crustalOutputDir.exists() || crustalOutputDir.mkdir());
		File subOutputDir = new File(FIGURES_DIR, "sub_grid");
		Preconditions.checkState(subOutputDir.exists() || subOutputDir.mkdir());
		
		PRVI25_DeclusteringAlgorithms avgDecluster = PRVI25_DeclusteringAlgorithms.AVERAGE;
		PRVI25_SeisSmoothingAlgorithms avgSmooth = PRVI25_SeisSmoothingAlgorithms.AVERAGE;
		GriddedRegion fullGrid = new GriddedRegion(PRVI25_SeismicityRegions.CRUSTAL.load(), 0.1, GriddedRegion.ANCHOR_0_0);
		
		GridSourceList crustalGridProv = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED).requireModule(GridSourceList.class);
		
		plotPDFs(crustalOutputDir, "crustal_pdf", 0, avgDecluster, avgSmooth, fullGrid,
				"Crustal", PRVI25_SeismicityRegions.CRUSTAL);
//		plotPDFs(crustalOutputDir, "crustal_m5", 5d, avgDecluster, avgSmooth, fullGrid,
//				"Crustal", PRVI25_SeismicityRegions.CRUSTAL);
		plotNucleationRates(crustalOutputDir, "crustal_m5", crustalGridProv, 5d, TectonicRegionType.ACTIVE_SHALLOW, fullGrid, "Crustal Gridded", null);
		plotNucleationRates(crustalOutputDir, "crustal_m7", crustalGridProv, 7d, TectonicRegionType.ACTIVE_SHALLOW, fullGrid, "Crustal Gridded", null);
		
		plotPDFs(subOutputDir, "sub_interface_pdf", 0, avgDecluster, avgSmooth, fullGrid,
				"Interface Gridded", PRVI25_SeismicityRegions.CAR_INTERFACE, PRVI25_SeismicityRegions.MUE_INTERFACE);
		plotPDFs(subOutputDir, "sub_interface_m5", 5, avgDecluster, avgSmooth, fullGrid,
				"Interface Gridded", PRVI25_SeismicityRegions.CAR_INTERFACE, PRVI25_SeismicityRegions.MUE_INTERFACE);
		plotPDFs(subOutputDir, "sub_slab_pdf", 0, avgDecluster, avgSmooth, fullGrid,
				"Intraslab", PRVI25_SeismicityRegions.CAR_INTRASLAB, PRVI25_SeismicityRegions.MUE_INTRASLAB);
		plotPDFs(subOutputDir, "sub_slab_m5", 5, avgDecluster, avgSmooth, fullGrid,
				"Intraslab", PRVI25_SeismicityRegions.CAR_INTRASLAB, PRVI25_SeismicityRegions.MUE_INTRASLAB);
	}
	
	private static EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(PRVI25_GridSourceBuilder.OVERALL_MMIN, 7.95);
	private static DecimalFormat oDF = new DecimalFormat("0.#");
	
	private static void plotPDFs(File outputDir, String prefix, double magForRate, PRVI25_DeclusteringAlgorithms decluster,
			PRVI25_SeisSmoothingAlgorithms smooth, GriddedRegion gridReg, String name, PRVI25_SeismicityRegions... seisRegs) throws IOException {
		if (gridReg == null) {
			Preconditions.checkState(seisRegs.length == 1);
			gridReg = new GriddedRegion(seisRegs[0].load(), 0.1, GriddedRegion.ANCHOR_0_0);
		}
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		List<Region> plotRegions = new ArrayList<>();
		BitSet everSets = new BitSet(xyz.size());
		for (PRVI25_SeismicityRegions seisReg : seisRegs) {
			GriddedGeoDataSet pdf = smooth.loadXYZ(seisReg, decluster);
			
			boolean remap = !pdf.getRegion().equalsRegion(gridReg);
			
			double rateScalar = 1;
			if (magForRate > 0d)
				rateScalar = PRVI25_RegionalSeismicity.PREFFERRED.build(seisReg, refMFD,
						refMFD.getMaxX()).getCumRate(refMFD.getClosestXIndex(magForRate+0.01));
			
			Region reg = seisReg.load();
			if (!reg.isRectangular())
				plotRegions.add(reg);
			
			for (int i=0; i<pdf.size(); i++) {
				int index = remap ? gridReg.indexForLocation(pdf.getLocation(i)) : i;
				Preconditions.checkState(index >= 0);
				everSets.set(index);
				xyz.set(index, xyz.get(index) + rateScalar*pdf.get(i));
			}
		}
		for (int i=0; i<xyz.size(); i++)
			if (!everSets.get(i))
				xyz.set(i, Double.NaN);
		
		CPT cpt;
		String label;
		if (magForRate > 0d) {
			cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-7, -2);
			xyz.log10();
			label = "Log10 "+name+" M>"+oDF.format(magForRate)+" Nucleation Rate";
		} else {
			cpt = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(-7d, -2d);
			xyz.log10();
			label = "Log10 "+name+" Fractional Rate";
		}
		cpt.setNanColor(new Color(255, 255, 255));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setWriteGeoJSON(false);
		
		if (!plotRegions.isEmpty())
			mapMaker.plotInsetRegions(plotRegions, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY), null, 1);
		mapMaker.plotXYZData(xyz, cpt, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}
	
	private static void plotNucleationRates(File outputDir, String prefix, GridSourceList gridSources,
			double minMag, TectonicRegionType trt, GriddedRegion gridReg, String name,
			List<Region> plotRegions) throws IOException {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		BitSet everSets = new BitSet(xyz.size());
		boolean remap = gridSources.getGriddedRegion() == null || !gridReg.equalsRegion(gridSources.getGriddedRegion());
		
		for (int l=0; l<gridSources.getNumLocations(); l++) {
			int index = remap ? gridReg.indexForLocation(gridSources.getLocation(l)) : l;
			if (index < 0)
				continue;
			everSets.set(index);
			xyz.set(index, gridSources.getCumulativeNucleationRate(trt, l, minMag));
		}
		for (int i=0; i<xyz.size(); i++)
			if (!everSets.get(i))
				xyz.set(i, Double.NaN);
		xyz.log10();
		CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-7, -2);
		String label = "Log10 "+name+" M>"+oDF.format(minMag)+" Nucleation Rate";
		cpt.setNanColor(new Color(255, 255, 255));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setWriteGeoJSON(false);
		
		if (plotRegions != null && !plotRegions.isEmpty())
			mapMaker.plotInsetRegions(plotRegions, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY), null, 1);
		mapMaker.plotXYZData(xyz, cpt, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}

}
