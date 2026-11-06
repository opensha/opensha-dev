package scratch.kevin.nshm23.griddedUncert;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.ColorUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class GriddedUncertExplorationFigures {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/grid_uncert");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<Region> regions = new ArrayList<>();
		List<String> regNames = new ArrayList<>();
		List<String> regPrefixes = new ArrayList<>();
		
		NSHM23_RegionLoader.SeismicityRegions seisReg = NSHM23_RegionLoader.SeismicityRegions.CONUS_WEST;
		
		double mMax = 7.65;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, mMax-0.1);
		
		regions.add(NSHM23_RegionLoader.loadFullConterminousWUS());
		regNames.add("WUS");
		regPrefixes.add("wus");
		
		regions.add(NSHM23_RegionLoader.AnalysisRegions.CONUS_U3_RELM.load());
		regNames.add("CA");
		regPrefixes.add("ca");
		
		regions.add(NSHM23_RegionLoader.LocalRegions.CONUS_SF_BAY.load());
		regNames.add("SF");
		regPrefixes.add("sf");
		
		regions.add(new Region(new Location(34.05, -118.35), new Location(34.15, -118.25)));
		regNames.add("SF 0.1deg cell");
		regPrefixes.add("sf_cell");
		
		regions.add(NSHM23_RegionLoader.LocalRegions.CONUS_LA_BASIN.load());
		regNames.add("LA");
		regPrefixes.add("la");
		
		regions.add(new Region(new Location(37.75, -122.55), new Location(37.85, -122.45)));
		regNames.add("LA 0.1deg cell");
		regPrefixes.add("la_cell");
		
//		NSHM23_DeclusteringAlgorithms.
		GriddedGeoDataSet avgPDF = NSHM23_SeisSmoothingAlgorithms.AVERAGE.loadXYZ(seisReg, NSHM23_DeclusteringAlgorithms.AVERAGE);
		WeightedList<GriddedGeoDataSet> branchPDFs = new WeightedList<>();
		for (NSHM23_DeclusteringAlgorithms decluster : NSHM23_DeclusteringAlgorithms.values()) {
			if (decluster.getNodeWeight() == 0d)
				continue;
			for (NSHM23_SeisSmoothingAlgorithms smooth : NSHM23_SeisSmoothingAlgorithms.values()) {
				if (smooth.getNodeWeight() == 0d)
					continue;
				double weight = decluster.getNodeWeight() * smooth.getNodeWeight();
				branchPDFs.add(smooth.loadXYZ(seisReg, decluster), weight);
			}
		}
		if (!branchPDFs.isNormalized())
			branchPDFs.normalize();
		
		IncrementalMagFreqDist prefTotal = NSHM23_RegionalSeismicity.PREFFERRED.build(seisReg, refMFD, mMax);
		IncrementalMagFreqDist highTotal = NSHM23_RegionalSeismicity.HIGH.build(seisReg, refMFD, mMax);
		IncrementalMagFreqDist lowTotal = NSHM23_RegionalSeismicity.LOW.build(seisReg, refMFD, mMax);
		
		UncertainBoundedIncrMagFreqDist totalBounds = new UncertainBoundedIncrMagFreqDist(
				prefTotal, lowTotal, highTotal, UncertaintyBoundType.CONF_95);
		
		for (int r=0; r<regions.size(); r++) {
			Region reg = regions.get(r);
			String name = regNames.get(r);
			String prefix = regPrefixes.get(r);
			
			System.out.println("Doing "+name);
			
			List<IncrementalMagFreqDist> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			totalBounds.setName(null);
			funcs.add(totalBounds);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 0, 0, 80)));
			
			prefTotal.setName("Full regional mean & 95%");
			funcs.add(prefTotal);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.BLACK));
			
			UncertainBoundedIncrMagFreqDist remapped = NSHM23_RegionalSeismicity.getRemapped(reg,
					NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, mMax);
			
			double avgPDFinReg = 0d;
			for (int i=0; i<avgPDF.size(); i++)
				if (reg.contains(avgPDF.getLocation(i)))
					avgPDFinReg += avgPDF.get(i);
			System.out.println("PDF in region:\t"+(float)avgPDFinReg);
			
			IncrementalMagFreqDist prefScaled = prefTotal.deepClone();
			prefScaled.scale(avgPDFinReg);
			IncrementalMagFreqDist highScaled = highTotal.deepClone();
			highScaled.scale(avgPDFinReg);
			IncrementalMagFreqDist lowScaled = lowTotal.deepClone();
			lowScaled.scale(avgPDFinReg);
			
			UncertainBoundedIncrMagFreqDist scaledBounds = new UncertainBoundedIncrMagFreqDist(
					prefScaled, lowScaled, highScaled, UncertaintyBoundType.CONF_95);
			scaledBounds.setName(null);
			funcs.add(scaledBounds);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, ColorUtils.transparent(Colors.tab_red, 80)));
			
			prefScaled.setName("Scaled regional mean & 95%");
			funcs.add(prefScaled);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_red));
			
			IncrementalMagFreqDist remappedUpper = remapped.getUpper().deepClone();
			IncrementalMagFreqDist remappedLower = remapped.getLower().deepClone();
			remappedUpper.setName("Sqrt(N) expanded bounds");
			remappedLower.setName(null);
			
			funcs.add(remappedUpper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_blue));
			funcs.add(remappedLower);
			chars.add(chars.getLast());
			
			for (int i=0; i<branchPDFs.size(); i++) {
				GriddedGeoDataSet pdf = branchPDFs.getValue(i);
				double pdfInReg = 0d;
				for (int j=0; j<pdf.size(); j++)
					if (reg.contains(pdf.getLocation(j)))
						pdfInReg += pdf.get(j);
				
				IncrementalMagFreqDist prefPDF = prefTotal.deepClone();
				prefPDF.scale(pdfInReg);
				IncrementalMagFreqDist highPDF = highTotal.deepClone();
				highPDF.scale(pdfInReg);
				IncrementalMagFreqDist lowPDF = lowTotal.deepClone();
				lowPDF.scale(pdfInReg);
				
				if (i == 0)
					prefPDF.setName("Indv. PDF mean & 95% bounds");
				else
					prefPDF.setName(null);
				highPDF.setName(null);
				lowPDF.setName(null);
				
				funcs.add(prefPDF);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.DARK_GRAY));
				
				funcs.add(highPDF);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.GRAY));
				funcs.add(lowPDF);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.GRAY));
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, name, "Magnitude", "Incremental rate (1/yr)");
			plot.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			gp.drawGraphPanel(plot, false, true, new Range(5d, 8d), new Range(1e-8, 1e1));
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 800, true, true, false);
		}
	}

}
