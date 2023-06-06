package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_MaxMagOffFault;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class SeisGridPlotter {

	public static void main(String[] args) throws IOException {
		
		File outputDir = new File("/home/kevin/markdown/nshm23-misc/seis_pdfs");
		
		List<String> lines = new ArrayList<>();

		double maxMagOff = NSHM23_MaxMagOffFault.MAG_7p6.getMaxMagOffFault();
		List<? extends FaultSection> subSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v2);
		
//		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		Region region = Region.union(SeismicityRegions.CONUS_WEST.load(), SeismicityRegions.CONUS_EAST.load());
		List<SeismicityRegions> seisRegions = NSHM23_InvConfigFactory.getSeismicityRegions(region);
		GriddedRegion gridReg = NSHM23_InvConfigFactory.getGriddedSeisRegion(seisRegions); 
		RupSetMapMaker mapMaker = new RupSetMapMaker(subSects, region);
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, new Color(0, 0, 0, 80)));
		mapMaker.setSectOutlineChar(null);
		mapMaker.setWritePDFs(true);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8d);
		
		CPT linearCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 0.01);
		CPT logCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-7, -1);
		
		CPT rawRatioCPT = GMT_CPT_Files.DIVERGING_BLUE_RED_UNIFORM.instance().rescale(0d, 1d);
		CPT ratioCPT = new CPT();
		double prevX = Double.NaN;
		Color prevColor = null;
		for (double x=0.5d; (float)x<1f; x+=0.01) {
			double remappedX = x-0.5d;
			Color color = rawRatioCPT.getColor((float)remappedX);
			
			if (prevColor != null)
				ratioCPT.add(new CPTVal((float)prevX, prevColor, (float)x, color));
			
			prevColor = color;
			prevX = x;
		}
		for (double x=1d; (float)x<=2f; x+=0.01) {
			double remappedX = 0.5*x;
			Color color = rawRatioCPT.getColor((float)remappedX);
			
			ratioCPT.add(new CPTVal((float)prevX, prevColor, (float)x, color));
			
			prevColor = color;
			prevX = x;
		}
		ratioCPT.setBelowMinColor(rawRatioCPT.getBelowMinColor());
		ratioCPT.setAboveMaxColor(rawRatioCPT.getAboveMaxColor());
		ratioCPT.setNanColor(rawRatioCPT.getNanColor());
//		CPT ratioCPT = GMT_CPT_Files.DIVERGING_BLUE_RED_UNIFORM.instance().rescale(0.5, 2d);
		CPT logRatioCPT = GMT_CPT_Files.DIVERGING_BLUE_RED_UNIFORM.instance().rescale(-1d, 1d);
		
		lines.add("# NSHM23 Spatial Seismicity PDFs");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		for (NSHM23_RegionalSeismicity seisRates : NSHM23_RegionalSeismicity.values()) {
			lines.add("## "+seisRates.getName());
			lines.add(topLink); lines.add("");
			
			GriddedGeoDataSet averagePDF = buildPDF(seisRates, NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE,
					maxMagOff, gridReg, refMFD, seisRegions);
			GriddedGeoDataSet logAveragePDF = averagePDF.copy();
			logAveragePDF.log10();
			
			lines.add("### Average PDF, "+seisRates.getName());
			lines.add(topLink); lines.add("");
			
			TableBuilder linearRatioTable = null;
			int linearTableIndex = -1;
			TableBuilder logRatioTable = null;
			int logTableIndex = -1;
			
			for (boolean log : new boolean[] {false,true}) {
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine(log ? "Log10, Average PDF" : "Linear, Average PDF");
				
				String prefix = seisRates.getFilePrefix()+"_average";
				if (log)
					prefix += "_log";
				if (log) {
					mapMaker.plotXYZData(logAveragePDF, logCPT, "Log10 (Rate M≥5)");
					mapMaker.plot(outputDir, prefix, "Average PDF, "+seisRates.getName());
				} else {
					mapMaker.plotXYZData(averagePDF, linearCPT, "Rate M≥5");
					mapMaker.plot(outputDir, prefix, "Average PDF, "+seisRates.getName());
				}
				table.addLine("![Map]("+prefix+".png)");
				
				lines.addAll(table.build());
				lines.add("");
				if (log) {
					logTableIndex = lines.size();
					logRatioTable = MarkdownUtils.tableBuilder();
				} else {
					linearTableIndex = lines.size();
					linearRatioTable = MarkdownUtils.tableBuilder();
				}
				lines.add("");
			}
			
			// plot moment rate
			GriddedGeoDataSet moRates = calcMoRate(seisRates, NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE,
					maxMagOff, gridReg, refMFD, seisRegions);
			
			double minMoRate = Double.POSITIVE_INFINITY;
			double maxMoRate = Double.NEGATIVE_INFINITY;
			for (int i=0; i<moRates.size(); i++) {
				double val = moRates.get(i);
				if (val > 0d) {
					minMoRate = Math.min(minMoRate, val);
					maxMoRate = Math.max(maxMoRate, val);
				}
			}
			
			double logMaxMo = Math.ceil(Math.log10(maxMoRate));
			double logMinMo = Math.max(logMaxMo-10d, Math.floor(Math.log10(minMoRate)));
			
			CPT moCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMinMo, logMaxMo);
			moCPT.setNanColor(Color.WHITE);
			
			String prefix = seisRates.getFilePrefix()+"_average";
			
			File moRatesXYZ = new File(outputDir, prefix+"_mo_rates.xyz");
			GriddedGeoDataSet.writeXYZFile(moRates, moRatesXYZ);
			
			// log10
			moRates.log10();
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			table.addLine("Total Moment Rate");
			
			mapMaker.plotXYZData(moRates, moCPT, "Log10 Total Moment Rate (N-m/yr)");
			mapMaker.plot(outputDir, prefix+"_mo_rates", " ");
			
			table.addLine("![Map]("+prefix+"_mo_rates.png)");
			
			lines.add("");
			lines.add("Download moment rate XYZ data: "+ "["+moRatesXYZ.getName()+"]("+moRatesXYZ.getName()+")");
			lines.add("");
			lines.addAll(table.build());
			lines.add("");

			for (NSHM23_SeisSmoothingAlgorithms smooth : NSHM23_SeisSmoothingAlgorithms.values()) {
				if (smooth.getNodeWeight(null) == 0d)
					continue;
				for (NSHM23_DeclusteringAlgorithms declustering : NSHM23_DeclusteringAlgorithms.values()) {
					if (declustering.getNodeWeight(null) == 0d)
						continue;
					lines.add("### "+declustering.getName()+", "+smooth.getName()+
							(seisRates == NSHM23_RegionalSeismicity.PREFFERRED ? "" : ", "+seisRates.getName()));
					lines.add(topLink); lines.add("");
					System.out.println("Processing "+smooth.getName()+", "+seisRates.getName());
					table = MarkdownUtils.tableBuilder();
					table.addLine(MarkdownUtils.boldCentered(smooth.getName()), "Ratio to Mean");
					
					GriddedGeoDataSet combinedPDF = buildPDF(seisRates, declustering, smooth, maxMagOff, gridReg, refMFD, seisRegions);
					
					String title = declustering.getName()+", "+smooth.getName();
					if (seisRates != NSHM23_RegionalSeismicity.PREFFERRED)
						title += ", "+seisRates.getName();
					
					table.initNewLine();
					mapMaker.plotXYZData(combinedPDF, linearCPT, "Rate M≥5");
					prefix = declustering.getFilePrefix()+"_"+seisRates.getFilePrefix()+"_"+smooth.getFilePrefix();
					mapMaker.plot(outputDir, prefix, title);
					table.addColumn("![Map]("+prefix+".png)");
					
					GeoDataSet ratio = GeoDataSetMath.divide(combinedPDF, averagePDF);
					mapMaker.plotXYZData(ratio, ratioCPT, "Rate / Avg Rate");
					mapMaker.plot(outputDir, prefix+"_ratio", title);
					table.addColumn("![Map]("+prefix+"_ratio.png)");
					table.finalizeLine();
					
					linearRatioTable.addColumn("![Map]("+prefix+"_ratio.png)");
					
					combinedPDF.log10();
					table.initNewLine();
					mapMaker.plotXYZData(combinedPDF, logCPT, "Log10(Rate M≥5)");
					mapMaker.plot(outputDir, prefix+"_log", title);
					table.addColumn("![Map]("+prefix+"_log.png)");
					
					ratio.log10();
					mapMaker.plotXYZData(ratio, logRatioCPT, "Log10(Rate / Avg Rate)");
					mapMaker.plot(outputDir, prefix+"_log_ratio", title);
					table.addColumn("![Map]("+prefix+"_log_ratio.png)");
					table.finalizeLine();
					
					logRatioTable.addColumn("![Map]("+prefix+"_log_ratio.png)");
					
					lines.addAll(table.build());
					lines.add("");
				}
			}
			
			// add ratio tables
			
			// log first as it's further down
			lines.addAll(logTableIndex, logRatioTable.wrap(3, 0).build());
			lines.addAll(linearTableIndex, linearRatioTable.wrap(3, 0).build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static GriddedGeoDataSet buildPDF(NSHM23_RegionalSeismicity seisRates, NSHM23_DeclusteringAlgorithms declusteringAlg, NSHM23_SeisSmoothingAlgorithms smooth,
			double maxMagOff, GriddedRegion gridReg, EvenlyDiscretizedFunc refMFD, List<SeismicityRegions> seisRegions) throws IOException {
		GriddedGeoDataSet combinedPDF = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<combinedPDF.size(); i++)
			combinedPDF.set(i, Double.NaN);
		
		for (SeismicityRegions primaryReg : seisRegions) {
			IncrementalMagFreqDist totalGR = seisRates.build(primaryReg, refMFD, maxMagOff);
			double rateM5 = totalGR.getCumRate(totalGR.getClosestXIndex(5.01));
			
			GriddedGeoDataSet regPDF = smooth.loadXYZ(primaryReg, declusteringAlg);
			regPDF.scale(rateM5);
			
			for (int i=0; i<regPDF.size(); i++) {
				Location loc = regPDF.getLocation(i);
				int index = gridReg.indexForLocation(loc);
				if (index >= 0)
					combinedPDF.set(index, regPDF.get(i));
			}
		}
		return combinedPDF;
	}
	
	private static GriddedGeoDataSet calcMoRate(NSHM23_RegionalSeismicity seisRates, NSHM23_DeclusteringAlgorithms declusteringAlg, NSHM23_SeisSmoothingAlgorithms smooth,
			double maxMagOff, GriddedRegion gridReg, EvenlyDiscretizedFunc refMFD, List<SeismicityRegions> seisRegions) throws IOException {
		GriddedGeoDataSet moRates = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<moRates.size(); i++)
			moRates.set(i, Double.NaN);
		
		for (SeismicityRegions primaryReg : seisRegions) {
			IncrementalMagFreqDist totalGR = seisRates.build(primaryReg, refMFD, maxMagOff);
			double totMoment = totalGR.getTotalMomentRate();
			GriddedGeoDataSet regPDF = smooth.loadXYZ(primaryReg, declusteringAlg);
			
			for (int i=0; i<regPDF.size(); i++) {
				Location loc = regPDF.getLocation(i);
				int index = gridReg.indexForLocation(loc);
				if (index >= 0)
					moRates.set(index, totMoment*regPDF.get(i));
			}
		}
		return moRates;
	}

}
