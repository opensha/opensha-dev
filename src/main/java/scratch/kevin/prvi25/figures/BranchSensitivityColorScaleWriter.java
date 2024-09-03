package scratch.kevin.prvi25.figures;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

public class BranchSensitivityColorScaleWriter {

	public static void main(String[] args) throws IOException {
		CPT cpt = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-50d, 50d);
		
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(101, 101, 0d, 0d, 0.01);
		
		XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, null, null, null, "Branch Choice / Overall Mean, % Change");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec);
		
		PlotUtils.writePlots(new File("/tmp"), "branch_cpt", gp, 700, 650, true, true, false);
	}

}
