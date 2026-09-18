package scratch.kevin.nshm23.hazardValidation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.ui.RectangleEdge;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm.DistanceCorrection2013;
import org.opensha.sha.util.TectonicRegionType;

public class CorrRJBPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		
		CPT magCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(6d, 8d);
		
		EvenlyDiscretizedFunc magRange = new EvenlyDiscretizedFunc(6.05, 19, 0.1);
		EvenlyDiscretizedFunc distRange = new EvenlyDiscretizedFunc(0d, 100, 1d);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DistanceCorrection2013 corr = new DistanceCorrection2013();
		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
		
		for (int m=0; m<magRange.size(); m++) {
			double mag = magRange.getX(m);
			EvenlyDiscretizedFunc rJBfunc = distRange.deepClone();
			for (int d=0; d<rJBfunc.size(); d++)
				rJBfunc.set(d, corr.getCorrectedDistanceJB(mag, rJBfunc.getX(d), trt));
			
			Color color = magCPT.getColor(mag);
			
			funcs.add(rJBfunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
		}
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		PaintScaleLegend cptLegend = GraphPanel.getLegendForCPT(magCPT, "Magnitude",
				gp.getPlotPrefs(), 0.1, RectangleEdge.BOTTOM);
		PlotSpec plot = new PlotSpec(funcs, chars, "NSHM Dist-Corr", "Epicentral Distance (km)", "Distance J-B (km)");
		plot.addSubtitle(cptLegend);
		
		gp.drawGraphPanel(plot);
		
		PlotUtils.writePlots(outputDir, "nshm_rjb_corr", gp, 1000, 1000, true, true, false);
	}

}
