package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.simulators.RSQSimEvent;

import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class NGAW2MagAreaCompare {

	public static void main(String[] args) throws IOException {
		File csvFile = new File("/home/kevin/Simulators/nga_w2_flatfile_rotd50_d050.csv");
		CSVFile<String> csv = CSVFile.readFile(csvFile, false);
		
		int magCol = 9;
		int areaCol = 34;
		
		DefaultXY_DataSet ngaScatter = new DefaultXY_DataSet();
		ngaScatter.setName("NGA-W2");
		for (int row=1; row<csv.getNumRows(); row++) {
			if (csv.get(row, 0).isEmpty())
				continue;
			double mag = csv.getDouble(row, magCol);
			double area = csv.getDouble(row, areaCol);
			if (area > 0 && mag > 0)
				ngaScatter.set(area, mag);
		}
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
		
		DefaultXY_DataSet rsScatter = new DefaultXY_DataSet();
		rsScatter.setName("RSQSim");
		for (RSQSimEvent event : catalog.loader().iterable()) {
			double area = event.getArea() * 1e-6;
			double mag = event.getMagnitude();
			if (area > 0 && mag > 0)
				rsScatter.set(area, mag);
		}
		
		DefaultXY_DataSet rsDownsampled = new DefaultXY_DataSet();
		Random r = new Random();
		for (int i=0; i<ngaScatter.size(); i++)
			rsDownsampled.set(rsScatter.get(r.nextInt(rsScatter.size())));
		
		System.out.println("Done loading catalog");

		double minArea = Math.min(rsScatter.getMinX(), ngaScatter.getMinX());
		double maxArea = Math.max(rsScatter.getMaxX(), ngaScatter.getMaxX());
		double minMag = Math.min(rsScatter.getMinY(), ngaScatter.getMinY());
		double maxMag = Math.max(rsScatter.getMaxY(), ngaScatter.getMaxY());
		
		System.out.println("Raw area range: ["+minArea+", "+maxArea+"]");
		System.out.println("Raw mag range: ["+minMag+", "+maxMag+"]");
		
		maxMag = Math.ceil(maxMag);
		minMag = Math.floor(minMag);
		maxArea = Math.pow(10, Math.ceil(Math.log10(maxArea)));
		minArea = Math.pow(10, Math.floor(Math.log10(minArea)));
		
		System.out.println("Plot area range: ["+minArea+", "+maxArea+"]");
		System.out.println("Plot mag range: ["+minMag+", "+maxMag+"]");
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		PlotCurveCharacterstics ngaChar =
				new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLUE.darker());
		PlotCurveCharacterstics rsChar =
				new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.GRAY);
		PlotCurveCharacterstics rsDownChar =
				new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.RED.darker());
		
		funcs.add(ngaScatter);
		chars.add(ngaChar);
		PlotSpec ngaSpec = new PlotSpec(funcs, chars,
				"Mag-Area Scatters", "Area (km^2)", "Magnitude");
		
		funcs = new ArrayList<>();
		chars = new ArrayList<>();
		funcs.add(rsScatter);
		chars.add(rsChar);
		funcs.add(rsDownsampled);
		chars.add(rsDownChar);
		PlotSpec rsSpec = new PlotSpec(funcs, chars,
				"Mag-Area Scatters", "Area (km^2)", "Magnitude");
		
		funcs = new ArrayList<>();
		chars = new ArrayList<>();
		funcs.add(rsScatter);
		chars.add(rsChar);
		funcs.add(rsDownsampled);
		chars.add(rsDownChar);
		funcs.add(ngaScatter);
		chars.add(ngaChar);
		PlotSpec combSpec = new PlotSpec(funcs, chars,
				"Mag-Area Scatters", "Area (km^2)", "Magnitude");

		List<PlotSpec> specs = new ArrayList<>();
		specs.add(ngaSpec);
		specs.add(rsSpec);
		specs.add(combSpec);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(minArea, maxArea);
		Range yRange = new Range(minMag, maxMag);
		
		List<Range> xRanges = new ArrayList<>();
		List<Range> yRanges = new ArrayList<>();
		
		xRanges.add(xRange);
		for (int i=0; i<specs.size(); i++)
			yRanges.add(yRange);
		
		System.out.println("Plotting....");
		gp.drawGraphPanel(specs, true, false, xRanges, yRanges);
		gp.getYAxis().setTickLabelsVisible(false);
		
		File file = new File("/tmp/nga_w2_ma_comparision");
		gp.getChartPanel().setSize(800, 1200);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
//		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		System.out.println("DONE");
	}

}
