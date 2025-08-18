package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.netlib.util.intW;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionSlabMMax;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

class SeismicityPDFFigures {

	public static void main(String[] args) throws IOException {
		File figsDir = FIGURES_DIR;
//		File figsDir = new File("/tmp/test_pdfs");
		Preconditions.checkState(figsDir.exists() || figsDir.mkdir());
		File crustalOutputDir = new File(figsDir, "crustal_grid");
		Preconditions.checkState(crustalOutputDir.exists() || crustalOutputDir.mkdir());
		File subOutputDir = new File(figsDir, "sub_grid");
		Preconditions.checkState(subOutputDir.exists() || subOutputDir.mkdir());
		
		PRVI25_DeclusteringAlgorithms avgDecluster = PRVI25_DeclusteringAlgorithms.AVERAGE;
		PRVI25_SeisSmoothingAlgorithms avgSmooth = PRVI25_SeisSmoothingAlgorithms.AVERAGE;
		GriddedRegion fullGrid = new GriddedRegion(PRVI25_SeismicityRegions.CRUSTAL.load(), 0.1, GriddedRegion.ANCHOR_0_0);
		
		GridSourceList crustalGridProv = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED).requireModule(GridSourceList.class);
		
		plotPDFs(crustalOutputDir, "crustal_pdf", 0, avgDecluster, avgSmooth, fullGrid,
				"Crustal", PRVI25_SeismicityRegions.CRUSTAL);
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
		
		HazardModel prevModel = HazardModel.load(Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-prvi-2003-main"));
		// needs to be include for slab and interface
		NshmErf subERF = new NshmErf(prevModel, Set.of(TectonicRegionType.SUBDUCTION_SLAB, TectonicRegionType.SUBDUCTION_INTERFACE),
				IncludeBackgroundOption.INCLUDE);
		subERF.getTimeSpan().setDuration(1d);
		subERF.updateForecast();
		
		plotPrevNucleationRates(subOutputDir, "sub_slab_2003_m5", subERF, 5d,
				TectonicRegionType.SUBDUCTION_SLAB, fullGrid, "2003 Intraslab", null);
		plotPrevNucleationRates(subOutputDir, "sub_slab_2003_m5_50km", subERF, 5d,
				TectonicRegionType.SUBDUCTION_SLAB, fullGrid, "2003 Intraslab (50km)", Range.closedOpen(0d, 79d));
		plotPrevNucleationRates(subOutputDir, "sub_slab_2003_m5_80km", subERF, 5d,
				TectonicRegionType.SUBDUCTION_SLAB, fullGrid, "2003 Intraslab (80km)", Range.closedOpen(80d, 119d));
		plotPrevNucleationRates(subOutputDir, "sub_slab_2003_m5_120km", subERF, 5d,
				TectonicRegionType.SUBDUCTION_SLAB, fullGrid, "2003 Intraslab (120km)", Range.closedOpen(120d, 1000d));
		plotPrevNucleationRates(subOutputDir, "sub_interface_2003_m5", subERF, 5d,
				TectonicRegionType.SUBDUCTION_INTERFACE, fullGrid, "2003 Interface", null);
		plotPrevNucleationRates(subOutputDir, "sub_interface_2003_m7", subERF, 7d,
				TectonicRegionType.SUBDUCTION_INTERFACE, fullGrid, "2003 Interface", null);
		
		GridSourceList subGridProv = FaultSystemSolution.load(SUBDUCTION_SOLS_COMBINED).requireModule(GridSourceList.class);
		plotSlabDepths(subOutputDir, "sub_slab_depths", subGridProv, fullGrid,
				List.of(PRVI25_SeismicityRegions.CAR_INTRASLAB.load(), PRVI25_SeismicityRegions.MUE_INTRASLAB.load()));
		
		NshmErf crustalERF = new NshmErf(prevModel, Set.of(TectonicRegionType.ACTIVE_SHALLOW),
				IncludeBackgroundOption.ONLY);
		crustalERF.getTimeSpan().setDuration(1d);
		crustalERF.updateForecast();
		
		plotPrevNucleationRates(crustalOutputDir, "crustal_2003_m5", crustalERF, 5d, TectonicRegionType.ACTIVE_SHALLOW, fullGrid, "2003 Crustal Gridded", null);
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
		double totalRate = 0d;
		for (PRVI25_SeismicityRegions seisReg : seisRegs) {
			GriddedGeoDataSet pdf = smooth.loadXYZ(seisReg, decluster);
			
			boolean remap = !pdf.getRegion().equalsRegion(gridReg);
			
			double rateScalar = 1;
			if (magForRate > 0d) {
				double sumWeight = 0d;
				rateScalar = 0d;
				for (PRVI25_SeismicityRateEpoch epoch : PRVI25_SeismicityRateEpoch.values()) {
					double weight = epoch.getNodeWeight(null);
					if (weight == 0d)
						continue;
					IncrementalMagFreqDist mfd;
					switch (seisReg) {
					case CRUSTAL:
						mfd = PRVI25_CrustalSeismicityRate.PREFFERRED.build(epoch, refMFD, refMFD.getMaxX());
						break;
					case CAR_INTERFACE:
						mfd = PRVI25_SubductionCaribbeanSeismicityRate.PREFFERRED.build(epoch, refMFD, refMFD.getMaxX(), Double.NaN, false);
						break;
					case CAR_INTRASLAB:
						mfd = PRVI25_SubductionCaribbeanSeismicityRate.PREFFERRED.build(epoch, refMFD,
								PRVI25_SubductionSlabMMax.MAG_8p0.getIncrementalMmax(), PRVI25_GridSourceBuilder.SLAB_M_CORNER, true);
						break;
					case MUE_INTERFACE:
						mfd = PRVI25_SubductionMuertosSeismicityRate.PREFFERRED.build(epoch, refMFD, refMFD.getMaxX(), Double.NaN, false);
						break;
					case MUE_INTRASLAB:
						mfd = PRVI25_SubductionMuertosSeismicityRate.PREFFERRED.build(epoch, refMFD,
								PRVI25_SubductionSlabMMax.MAG_8p0.getIncrementalMmax(), PRVI25_GridSourceBuilder.SLAB_M_CORNER, true);
						break;

					default:
						throw new IllegalStateException();
					}
					rateScalar += weight * mfd.getCumRate(refMFD.getClosestXIndex(magForRate+0.01));
					sumWeight += weight;
				}
				rateScalar /= sumWeight;
				totalRate += rateScalar;
			}
			
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
			if (!everSets.get(i) || xyz.get(i) == 0d)
				xyz.set(i, Double.NaN);
		
		if (magForRate > 0d) {
			System.out.println(name);
			System.out.println("Rate M>"+(float)magForRate+": "+(float)totalRate);
		}
		
		CPT cpt;
		String label;
		if (magForRate > 0d) {
			cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-7, -2);
			cpt.setLog10(true);
			label = name+" M>"+oDF.format(magForRate)+" Nucleation Rate";
		} else {
			cpt = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(-7d, -2d);
			cpt.setLog10(true);
			label = name+" Fractional Rate";
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
		CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-7, -2);
		cpt.setLog10(true);
		String label = name+" M>"+oDF.format(minMag)+" Nucleation Rate";
		cpt.setNanColor(new Color(255, 255, 255));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setWriteGeoJSON(false);
		
		if (plotRegions != null && !plotRegions.isEmpty())
			mapMaker.plotInsetRegions(plotRegions, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY), null, 1);
		mapMaker.plotXYZData(xyz, cpt, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}
	
	private static void plotPrevNucleationRates(File outputDir, String prefix, NshmErf erf, double minMag,
			TectonicRegionType trt, GriddedRegion gridReg, String name, Range<Double> depthRange) throws IOException {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		HashSet<Double> slabDepths = trt == TectonicRegionType.SUBDUCTION_SLAB ? new HashSet<>() : null;
		double totalRate = 0d;
		for (ProbEqkSource source : erf) {
			if (source.getTectonicRegionType() == trt) {
//				if (trt == TectonicRegionType.SUBDUCTION_INTERFACE)
//					System.out.println("Including source "+source.getName()+" ("+source+") for trt="+trt.name()
//							+" with "+source.getNumRuptures()+" ruptures");
//				RuptureSurface sourceSurf = source.getRupture(0).getRuptureSurface();
//				LocationList sourceLocs = sourceSurf.getEvenlyDiscritizedListOfLocsOnSurface();
//				Preconditions.checkState(sourceLocs.size() == 1, "Not a point surface (%s locs): %s", sourceLocs.size(), sourceSurf);
//				Location loc = sourceLocs.get(0);
//				Preconditions.checkState(sourceSurf instanceof PointSurface, "Not a point surface: %s", sourceSurf);
//				Location loc = ((PointSurface)sourceSurf).getLocation();
//				int index = gridReg.indexForLocation(loc);
//				if (index >= 0) {
//					double rate = xyz.get(index);
//					for (ProbEqkRupture rup : source)
//						if (rup.getMag() >= minMag)
//							rate += rup.getMeanAnnualRate(1d);
//					xyz.set(index, rate);
//				}
				double sourceMinMag = Double.POSITIVE_INFINITY;
				double sourceMaxMag = 0d;
				for (ProbEqkRupture rup : source) {
					sourceMinMag = Math.min(sourceMinMag, rup.getMag());
					sourceMaxMag = Math.max(sourceMaxMag, rup.getMag());
					RuptureSurface sourceSurf = rup.getRuptureSurface();
					LocationList sourceLocs = sourceSurf.getEvenlyDiscritizedListOfLocsOnSurface();
					double rupRate = rup.getMeanAnnualRate(1d);
					totalRate += rupRate;
					double rupRateEach = rupRate / sourceLocs.size();
					for (Location loc : sourceLocs) {
						if (depthRange != null && !depthRange.contains(loc.depth))
							continue;
						if (slabDepths != null)
							slabDepths.add(loc.depth);
						int index = gridReg.indexForLocation(loc);
						if (index >= 0)
							xyz.add(index, rupRateEach);
					}
				}
//				if (trt == TectonicRegionType.SUBDUCTION_INTERFACE)
//					System.out.println("\tMag range: ["+(float)sourceMinMag+", "+(float)sourceMaxMag+"]");
			}
		}
		System.out.println(name);
		System.out.println("Rate M>"+(float)minMag+": "+(float)totalRate);
		if (slabDepths != null)
			System.out.println("Slab depths: "+slabDepths);
		for (int i=0; i<xyz.size(); i++)
			if (xyz.get(i) == 0d)
				xyz.set(i, Double.NaN);
		CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-7, -2);
		String label = name+" M>"+oDF.format(minMag)+" Nucleation Rate";
		cpt.setLog10(true);
		cpt.setNanColor(new Color(255, 255, 255));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setWriteGeoJSON(false);
		
		mapMaker.plotXYZData(xyz, cpt, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}
	
	private static void plotSlabDepths(File outputDir, String prefix, GridSourceList gridSources,
			GriddedRegion gridReg, List<Region> plotRegions) throws IOException {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<xyz.size(); i++)
			xyz.set(i, Double.POSITIVE_INFINITY);
		BitSet everSets = new BitSet(xyz.size());
		boolean remap = gridSources.getGriddedRegion() == null || !gridReg.equalsRegion(gridSources.getGriddedRegion());
		
		for (int l=0; l<gridSources.getNumLocations(); l++) {
			int index = remap ? gridReg.indexForLocation(gridSources.getLocation(l)) : l;
			if (index < 0)
				continue;
			everSets.set(index);
			
			double depth = xyz.get(index);
			for (GriddedRupture rup : gridSources.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, l))
				depth = Math.min(depth, rup.properties.upperDepth);
			xyz.set(index, depth);
		}
		for (int i=0; i<xyz.size(); i++)
			if (!everSets.get(i))
				xyz.set(i, Double.NaN);
//		CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance();
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
		cpt = cpt.rescale(0d, 80);
		String label = "Slab Depth (km)";
		cpt.setNanColor(new Color(255, 255, 255));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setWriteGeoJSON(false);
		
		if (plotRegions != null && !plotRegions.isEmpty())
			mapMaker.plotInsetRegions(plotRegions, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY), null, 1);
		mapMaker.plotXYZData(xyz, cpt, label);
		
		mapMaker.plot(outputDir, prefix, " ");
	}

}
