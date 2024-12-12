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
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;

import com.google.common.base.Preconditions;

class SeismicityPDFFigures {

	public static void main(String[] args) throws IOException {
		File figsDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures");
		File crustalOutputDir = new File(figsDir, "crustal_grid");
		Preconditions.checkState(crustalOutputDir.exists() || crustalOutputDir.mkdir());
		File subOutputDir = new File(figsDir, "sub_grid");
		Preconditions.checkState(subOutputDir.exists() || subOutputDir.mkdir());
		
		PRVI25_DeclusteringAlgorithms avgDecluster = PRVI25_DeclusteringAlgorithms.AVERAGE;
		PRVI25_SeisSmoothingAlgorithms avgSmooth = PRVI25_SeisSmoothingAlgorithms.AVERAGE;
		GriddedRegion fullGrid = new GriddedRegion(PRVI25_SeismicityRegions.CRUSTAL.load(), 0.1, GriddedRegion.ANCHOR_0_0);
		
		plotPDFs(crustalOutputDir, "crustal_pdf", 0, avgDecluster, avgSmooth, fullGrid,
				"Crustal", PRVI25_SeismicityRegions.CRUSTAL);
		plotPDFs(crustalOutputDir, "crustal_m5", 5d, avgDecluster, avgSmooth, fullGrid,
				"Crustal", PRVI25_SeismicityRegions.CRUSTAL);
		
		plotPDFs(subOutputDir, "sub_interface_pdf", 0, avgDecluster, avgSmooth, fullGrid,
				"Interface", PRVI25_SeismicityRegions.CAR_INTERFACE, PRVI25_SeismicityRegions.MUE_INTERFACE);
		plotPDFs(subOutputDir, "sub_interface_m5", 5, avgDecluster, avgSmooth, fullGrid,
				"Interface", PRVI25_SeismicityRegions.CAR_INTERFACE, PRVI25_SeismicityRegions.MUE_INTERFACE);
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

}
