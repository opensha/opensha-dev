package scratch.kevin.ucerf3;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

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
		
		File outputDir = new File("/tmp");
		
		double duration = 102020 - 2018;
		double rateEach = 1d/duration;
		
		for (Region region : regions) {
			IncrementalMagFreqDist compMFD =
					compSol.calcNucleationMFD_forRegion(region, minX, maxX, delta, true);
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
					int mfdIndex = compMFD.getClosestXIndex(mag);
					compMFD.add(mfdIndex, pt.getY());
				}
			}
			IncrementalMagFreqDist catalogHypoMFD =
					new IncrementalMagFreqDist(minX, maxX, compMFD.size());
			IncrementalMagFreqDist catalogFractMFD =
					new IncrementalMagFreqDist(minX, maxX, compMFD.size());
			
			double[] fractsInside = compRupSet.getFractRupsInsideRegion(region, true);
			
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
				} else if (region.contains(loc)) {
					catalogFractMFD.add(mfdIndex, rateEach);
				}
			}
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			compMFD.setName("FSS Zip File");
			funcs.add(compMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));

			catalogHypoMFD.setName("Catalog Hypocenter In Region");
			funcs.add(catalogHypoMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));

			catalogFractMFD.setName("Catalog Fract Surface In Region");
			funcs.add(catalogFractMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			
			PlotSpec plot = new PlotSpec(funcs, chars, region.getName(), "Magnitude", "Annual Frequency");
			plot.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setLegendFontSize(18);
			gp.setBackgroundColor(Color.WHITE);
			
			gp.drawGraphPanel(plot, false, true, null, new Range(1e-6, 1e0));
			gp.getChartPanel().setSize(1000, 800);
			String prefix = region.getName().replaceAll(" ", "_")+"_comparison";
			gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
		}
	}

}
