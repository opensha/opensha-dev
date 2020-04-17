package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimState;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class TransTotalPatchSlipCompare {

	public static void main(String[] args) throws IOException {
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_4950.instance();
		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/singleSS"),
				"Single SS", null, null);
//		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/singleSS_tri"),
//				"Single SS Triangles", null, null);
		
		System.out.println("Loading list slips...");
		double[] patchSlips = new double[catalog.getElements().size()];
		for (RSQSimEvent event : catalog.loader().iterable()) {
			int[] ids = event.getAllElementIDs();
			double[] slips = event.getAllElementSlips();
			
			for (int i=0; i<ids.length; i++)
				patchSlips[ids[i]-1] += slips[i];
		}
		
		RSQSimStateTransitionFileReader transReader = catalog.getTransitions();
		System.out.println("Trans version: "+transReader.getVersion());
		transReader.setQuiet(true);
		
		double[] patchTransSlips = new double[catalog.getElements().size()];
		
		int numSlips = 0;
		double totalTimeSlipping = 0d;
		
		System.out.println("Loading trans slips...");
		for (RSQSimStateTime trans : transReader.getTransitionsIterable(0d, Double.POSITIVE_INFINITY)) {
			if (trans.state == RSQSimState.EARTHQUAKE_SLIP) {
				double slip = trans.velocity*trans.getDuration();
				patchTransSlips[trans.patchID-1] += slip;
			}
		}
		
		System.out.println("Total number of slips: "+numSlips);
		System.out.println("Total time slipping: "+totalTimeSlipping);
		
		System.out.println("Plotting...");
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		double maxSlip = Math.max(StatUtils.max(patchTransSlips), StatUtils.max(patchSlips));
		double minSlip = Double.POSITIVE_INFINITY;
		for (int i=0; i<patchSlips.length; i++) {
			xy.set(patchSlips[i], patchTransSlips[i]);
			if (patchSlips[i] > 0)
				minSlip = Math.min(minSlip, patchSlips[i]);
			if (patchTransSlips[i] > 0)
				minSlip = Math.min(minSlip, patchTransSlips[i]);
		}
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Total Slip Comparison", "From List Files", "From Transitions");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(0, Math.ceil(maxSlip));
		Range yRange = xRange;
		
		XY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(xRange.getLowerBound(), xRange.getLowerBound());
		oneToOne.set(xRange.getUpperBound(), xRange.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File file = new File(catalog.getCatalogDir(), "trans_patch_total_slip_compare");
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		
		// now log
		
		xRange = new Range(Math.pow(10, Math.floor(Math.log10(minSlip))),
				Math.pow(10, Math.ceil(Math.log10(maxSlip))));
		yRange = xRange;
		
		oneToOne = new DefaultXY_DataSet();
		oneToOne.set(xRange.getLowerBound(), xRange.getLowerBound());
		oneToOne.set(xRange.getUpperBound(), xRange.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		file = new File(catalog.getCatalogDir(), "trans_patch_total_slip_compare_log");
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}

}
