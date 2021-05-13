package scratch.kevin.ucerf3.eal.spatialCorr;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

public class SpatiallyCorrelatedLossPlot {

	public static void main(String[] args) throws IOException {
		File resultsDir = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "2020_06_30-ucerf3-ngaw2-cea-100pct-spatial/results");
		
		int numZeros = 0;
		
		for (File file : resultsDir.listFiles()) {
			if (!file.getName().endsWith(".csv"))
				continue;
			CSVFile<String> csv = CSVFile.readFile(file, true);
			String name = file.getName();
			name = name.substring(0, name.indexOf(".csv"));
			name = name.replaceAll("_", " ");
			
			// between, mean, indv
			Table<Double, Double, List<Double>> resultsTable = HashBasedTable.create();
			
			for (int row=1; row<csv.getNumRows(); row++) {
				double meanLoss = csv.getDouble(row, 6);
				double between = csv.getDouble(row, 7);
				List<Double> vals = resultsTable.get(between, meanLoss);
				if (vals == null) {
					vals = new ArrayList<>();
					resultsTable.put(between, meanLoss, vals);
				}
				double val = csv.getDouble(row, 9);
				vals.add(val);
				if (val == 0d)
					numZeros++;
			}
			System.out.println("have "+numZeros+" zeroes");
			
			for (Double between : resultsTable.rowKeySet()) {
				Map<Double, List<Double>> map = resultsTable.row(between);
				DefaultXY_DataSet scatter = new DefaultXY_DataSet();
				DefaultXY_DataSet means = new DefaultXY_DataSet();
				
				for (Double meanLoss : map.keySet()) {
					List<Double> withinVals = map.get(meanLoss);
					double meanWithin = StatUtils.mean(Doubles.toArray(withinVals));
					
					for (double within : withinVals)
						scatter.set(meanLoss, within);
					
					means.set(meanLoss, meanWithin);
				}
				
				double maxVal = Math.max(scatter.getMaxX(), scatter.getMaxY());
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				Range range = new Range(0d, maxVal);
				
				DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
				oneToOne.set(0d,  0d);
				oneToOne.set(maxVal, maxVal);
				funcs.add(oneToOne);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY));
				
				scatter.setName("Indiv. within-events");
				funcs.add(scatter);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, new Color(0, 0, 0, 80)));
				
				means.setName("Mean within-event");
				funcs.add(means);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.RED));
				
				String title = name+", Between="+between.floatValue();
				System.out.println("Doing "+title);
				PlotSpec spec = new PlotSpec(funcs, chars, title,
						"Mean Loss", "Within-Event Variability Loss");
				spec.setLegendVisible(true);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(24);
				gp.setLegendFontSize(28);
				gp.setBackgroundColor(Color.WHITE);
				
				gp.drawGraphPanel(spec, false, false, range, range);
				
				String prefix = file.getName();
				prefix = prefix.substring(0, prefix.indexOf(".csv"));
				prefix += "_between_"+between.floatValue();
				File outFile = new File(resultsDir, prefix);
				gp.getChartPanel().setSize(800, 800);
				gp.saveAsPNG(outFile.getAbsolutePath()+".png");
			}
		}
	}

}
