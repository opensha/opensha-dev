package scratch.kevin.ucerf3;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;

public class CatalogMFDInRegionTest {
	
	private static final double minX = 5.05d;
	private static final double maxX = 9.05d;
	private static final double delta = 0.1d;

	public static void main(String[] args) throws IOException, DocumentException {
		File csvFile = new File("/home/kevin/Downloads/u3ti_100k_Catalog.csv");
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		Region[] regions = {
				new CaliforniaRegions.LA_BOX(),
				new CaliforniaRegions.SF_BOX(),
				new CaliforniaRegions.NORTHRIDGE_BOX(),
		};
		File compareFile = new File("/home/kevin/.opensha/ucerf3_erf/"
				+ "cached_dep100.0_depMean_rakeMean.zip");
		FaultSystemSolution compSol = FaultSystemIO.loadSol(compareFile);
		FaultSystemRupSet compRupSet = compSol.getRupSet();
		GridSourceProvider gridProv = compSol.getGridSourceProvider();
		
		String baseURL = "http://opensha.usc.edu/ftp/kmilner/ucerf3/2013_05_10-ucerf3p3-production-10runs/erf_mfd_plots/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_MFD_ERF_";
		
		File outputDir = new File("/tmp");
		
		double duration = 102020 - 2018;
		double rateEach = 1d/duration;
		
		for (Region region : regions) {
			double[] fractsInside = compRupSet.getFractRupsInsideRegion(region, true);
			
			IncrementalMagFreqDist compNuclMFD =
					compSol.calcNucleationMFD_forRegion(region, minX, maxX, delta, true);
			IncrementalMagFreqDist compPartMFD = new IncrementalMagFreqDist(minX, maxX, compNuclMFD.size());
			for (int fssIndex=0; fssIndex<fractsInside.length; fssIndex++)
				if (fractsInside[fssIndex] > 0)
					compPartMFD.add(compPartMFD.getClosestXIndex(compRupSet.getMagForRup(fssIndex)),
							compSol.getRateForRup(fssIndex));
			for (int index=0; index<gridProv.size(); index++) {
				Location loc = gridProv.getGriddedRegion().getLocation(index);
				if (!region.contains(loc))
					continue;
				IncrementalMagFreqDist nodeMFD = gridProv.getNodeMFD(index);
				if (nodeMFD == null)
					continue;
				for (Point2D pt : nodeMFD) {
					double mag = pt.getX();
					if (mag < minX-0.5*delta)
						continue;
					int mfdIndex = compNuclMFD.getClosestXIndex(mag);
					compNuclMFD.add(mfdIndex, pt.getY());
					compPartMFD.add(mfdIndex, pt.getY());
				}
			}
			IncrementalMagFreqDist catalogHypoMFD =
					new IncrementalMagFreqDist(minX, maxX, compNuclMFD.size());
			IncrementalMagFreqDist catalogFractMFD =
					new IncrementalMagFreqDist(minX, maxX, compNuclMFD.size());
			IncrementalMagFreqDist catalogPartMFD =
					new IncrementalMagFreqDist(minX, maxX, compNuclMFD.size());
			
			String regPrefix = region.getName().replaceAll(" ", "_");
			
			String onlineURL = baseURL+regPrefix+".txt";
			File targetFile = new File(outputDir, "orig_"+regPrefix+".txt");
			FileUtils.downloadURL(onlineURL, targetFile);
			boolean correctMFD = false;
			boolean insideData = false;
			ArbitrarilyDiscretizedFunc onlineCumulative = new ArbitrarilyDiscretizedFunc();
			for (String line : Files.readLines(targetFile, Charset.defaultCharset())) {
				line = line.trim();
				if (line.startsWith("Name:"))
					correctMFD = line.contains("UCERF3 MFDs (weighted mean)");
				if (line.isEmpty() || line.contains("DATASET"))
					insideData = false;
				if (line.startsWith("Data[x,y]:")) {
					insideData = true;
					continue;
				}
				if (correctMFD && insideData) {
					StringTokenizer tok = new StringTokenizer(line);
					double mag = Double.parseDouble(tok.nextToken());
					double rate = Double.parseDouble(tok.nextToken());
					onlineCumulative.set(mag, rate);
				}
			}
			
			for (int row=1; row<csv.getNumRows(); row++) {
				double lat = Double.parseDouble(csv.get(row, 6));
				double lon = Double.parseDouble(csv.get(row, 7));
				double depth = Double.parseDouble(csv.get(row, 8));
				double mag = Double.parseDouble(csv.get(row, 9));
				if (mag < minX-0.5*delta)
					continue;
				int mfdIndex = catalogFractMFD.getClosestXIndex(mag);
				int fssIndex = Integer.parseInt(csv.get(row, 16));
				Location loc = new Location(lat, lon, depth);
				if (region.contains(loc))
					catalogHypoMFD.add(mfdIndex, rateEach);
				
				if (fssIndex >= 0) {
					double solMag = compRupSet.getMagForRup(fssIndex);
					Preconditions.checkState((float)mag == (float)solMag);
					catalogFractMFD.add(mfdIndex, rateEach*fractsInside[fssIndex]);
					if (fractsInside[fssIndex] > 0)
						catalogPartMFD.add(mfdIndex, rateEach);
				} else if (region.contains(loc)) {
					catalogFractMFD.add(mfdIndex, rateEach);
					catalogPartMFD.add(mfdIndex, rateEach);
				}		
			}
			
			for (boolean cumulative : new boolean[] {false, true}) {
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				String yAxisLabel;
				if (cumulative) {
					EvenlyDiscretizedFunc cumNuclCompMFD = compNuclMFD.getCumRateDistWithOffset();
					EvenlyDiscretizedFunc cumPartCompMFD = compPartMFD.getCumRateDistWithOffset();
					EvenlyDiscretizedFunc cumCatalogHypoMFD = catalogHypoMFD.getCumRateDistWithOffset();
					EvenlyDiscretizedFunc cumCatalogFractMFD = catalogFractMFD.getCumRateDistWithOffset();
					EvenlyDiscretizedFunc cumCatalogPartMFD = catalogPartMFD.getCumRateDistWithOffset();
					
					onlineCumulative.setName("Online MFD");
					funcs.add(onlineCumulative);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
					
					cumNuclCompMFD.setName("FSS Zip File (Nucl)");
					funcs.add(cumNuclCompMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
					
					cumPartCompMFD.setName("FSS Zip File (Part)");
					funcs.add(cumPartCompMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.ORANGE));

					cumCatalogHypoMFD.setName("Catalog Hypocenter In Region");
					funcs.add(cumCatalogHypoMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));

					cumCatalogFractMFD.setName("Catalog Fract Surface In Region");
					funcs.add(cumCatalogFractMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));

					cumCatalogPartMFD.setName("Catalog Participation");
					funcs.add(cumCatalogPartMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.MAGENTA));
					
					yAxisLabel = "Cumulative Annual Frequency";
				} else {
					compNuclMFD.setName("FSS Zip File (Nucl)");
					funcs.add(compNuclMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
					
					compPartMFD.setName("FSS Zip File (Part)");
					funcs.add(compPartMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.ORANGE));

					catalogHypoMFD.setName("Catalog Hypocenter In Region");
					funcs.add(catalogHypoMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));

					catalogFractMFD.setName("Catalog Fract Surface In Region");
					funcs.add(catalogFractMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));

					catalogPartMFD.setName("Catalog Participation");
					funcs.add(catalogPartMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.MAGENTA));
					
					yAxisLabel = "Incremental Annual Frequency";
				}
				
				PlotSpec plot = new PlotSpec(funcs, chars, region.getName(), "Magnitude", yAxisLabel);
				plot.setLegendVisible(true);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.setLegendFontSize(18);
				gp.setBackgroundColor(Color.WHITE);
				
				gp.drawGraphPanel(plot, false, true, null, new Range(1e-8, 1e0));
				gp.getChartPanel().setSize(1000, 800);
				String prefix = region.getName().replaceAll(" ", "_")+"_comparison";
				if (cumulative)
					prefix += "_cumulative";
				gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
				gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
				gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
			}
		}
	}

}
