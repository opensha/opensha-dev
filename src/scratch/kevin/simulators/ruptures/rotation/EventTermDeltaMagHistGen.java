package scratch.kevin.simulators.ruptures.rotation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;

public class EventTermDeltaMagHistGen {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/git/rsqsim-analysis/catalogs/"
				+ "rundir2585_1myr/bbp_LA_BASIN_500");
		String csvName = "between_events_m7.2_50km_std_dev_3s_medians_table.csv";
		File outputFile = new File("/tmp/event_term_hists.png");
		
		Map<Double, File> dmDirs = new HashMap<>();
		
		dmDirs.put(0d, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface"));
		dmDirs.put(0.05, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface_rnd_mag_0p05"));
		dmDirs.put(0.1, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface_rnd_mag"));
		dmDirs.put(0.15, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface_rnd_mag_0p15"));
		dmDirs.put(0.2, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface_rnd_mag_0p2"));
		dmDirs.put(0.25, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface_rnd_mag_0p25"));
		dmDirs.put(0.3, new File(baseDir, "rotated_ruptures_m7p2_vert_ss_surface_rnd_mag_0p3"));
		
		Map<Double, List<Double>> termLists = new HashMap<>();
		
		MinMaxAveTracker eventTermTrack = new MinMaxAveTracker();
		
		for (Double dm : dmDirs.keySet()) {
			File dmDir = dmDirs.get(dm);
			File resourcesDir = new File(dmDir, "resources");
			File csvFile = new File(resourcesDir, csvName);
			CSVFile<String> csv = CSVFile.readFile(csvFile, true);
			List<Double> terms = new ArrayList<>();
			for (int row=1; row<csv.getNumRows(); row++) {
				double eventTerm = csv.getDouble(row, 1);
				eventTermTrack.addValue(eventTerm);
				terms.add(eventTerm);
			}
			termLists.put(dm, terms);
		}
		
		HistogramFunction bins = HistogramFunction.getEncompassingHistogram(
				eventTermTrack.getMin(), eventTermTrack.getMax(), 0.05);
		
		List<Double> dms = new ArrayList<>(dmDirs.keySet());
		Collections.sort(dms);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		CPT cpt = new CPT(0d, dms.size()-1, Color.BLACK, Color.LIGHT_GRAY);
		
		List<HistogramFunction> hists = new ArrayList<>();
		
		for (int j=0; j<dms.size(); j++) {
			Double dm = dms.get(j);
			HistogramFunction hist = new HistogramFunction(
					bins.getMinX(), bins.getMaxX(), bins.size());
			for (Double eventTerm : termLists.get(dm))
				hist.add(hist.getClosestXIndex(eventTerm), 1d);
			
			XY_DataSet xy = new DefaultXY_DataSet();
			double delta = hist.getDelta();
			double halfDelta = 0.5*delta;
			for (int i=0; i<hist.size(); i++) {
				double center = hist.getX(i);
				double y = hist.getY(i);
				if (i == 0)
					xy.set(center-halfDelta, 0d);
				xy.set(center-halfDelta, y);
				xy.set(center+halfDelta, y);
				if (i == hist.size()-1)
					xy.set(center+halfDelta, 0d);
			}
			
			if (dm > 0)
				xy.setName("dM="+dm.floatValue());
			else
				xy.setName("Closest");
			
			hists.add(hist);
			
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f,
					cpt.getColor((float)j)));
		}
		
//		for (int i=0; i<hists.size(); i++) {
//			funcs.add(hists.get(i));
//			Color c = cpt.getColor((float)i);
////			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f,
////					new Color(c.getRed(), c.getGreen(), c.getBlue(), 127)));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, c));
//		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "", "Event Term", "Count");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		gp.drawGraphPanel(spec, false, false, null, null);
		gp.getYAxis().setTickLabelsVisible(false);
		
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(outputFile.getAbsolutePath());
	}

}
