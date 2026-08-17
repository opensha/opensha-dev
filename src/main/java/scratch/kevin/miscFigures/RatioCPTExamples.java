package scratch.kevin.miscFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jfree.chart.ui.RectangleEdge;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

public class RatioCPTExamples {

	public static void main(String[] args) throws IOException {
		DefaultXY_DataSet fakeXY = new DefaultXY_DataSet(0d, 0d);
		
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		PlotPreferences prefs = gp.getPlotPrefs();
		
		PlotSpec plot = new PlotSpec(List.of(fakeXY), List.of(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.WHITE)), " ", null, null);
		
		List<CPT> cpts = new ArrayList<>();
		List<String> labels = new ArrayList<>();
		
		cpts.add(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging VIK");
		
		cpts.add(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1).expandAt(1d, 0.9, 1.1));
		labels.add("Diverging VIK (expanded)");
		
		cpts.add(GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1)
				.asDiscrete(0.2, true).mask(cpts.getLast().getColor(1d), 0.8, 1.2));
		labels.add("Diverging VIK (discrete)");
		
		cpts.add(GMT_CPT_Files.DIVERGENT_RYB.instance().reverse().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging RYB");
		
		cpts.add(GMT_CPT_Files.DIVERGENT_RYB.instance().reverse().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1).expandAt(Color.LIGHT_GRAY, 1, 0.9, 1.1));
		labels.add("Diverging RYB (expanded)");
		
		cpts.add(GMT_CPT_Files.DIVERGING_DARK_BLUE_RED_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging Dark Blue-Red");
		
		cpts.add(GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging BAM");
		
		cpts.add(GMT_CPT_Files.DIVERGING_BROC_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging BROC");
		
		cpts.add(GMT_CPT_Files.DIVERGING_CORK_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging CORK");
		
		cpts.add(GMT_CPT_Files.DIVERGING_BLUE_RED_UNIFORM.instance().rescale(0, 1).anchoredRescale(0, 3, 0.5, 1));
		labels.add("Diverging Blue-Red");
		
		Collections.reverse(cpts);
		Collections.reverse(labels);
		
		for (int i=0; i<cpts.size(); i++)
			plot.addSubtitle(GraphPanel.getLegendForCPT(cpts.get(i), labels.get(i), prefs, 0.2, RectangleEdge.BOTTOM));
		
		gp.drawGraphPanel(plot);
		
		PlotUtils.writePlots(new File("/tmp"), "ratio_cpts", gp, 800, 100+100*cpts.size(), true, false, false);
	}

}
