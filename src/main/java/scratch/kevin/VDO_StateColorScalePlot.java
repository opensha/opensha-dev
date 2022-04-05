package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.util.XMLUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;

import com.google.common.base.Preconditions;

public class VDO_StateColorScalePlot {

	public static void main(String[] args) throws DocumentException, IOException {
		File dir = new File("/home/kevin/Documents/papers/2021_UCERF4_Plausibility/figures/cff_vdo_figure");
		File file = new File(dir, "vdo_figure.xml");
		String cptLabel = "ΔCFF (MPa)";
		Document doc = XMLUtils.loadDocument(file);
		
		Element root = doc.getRootElement();
		
		Element pluginEl = root.element("UCERF3-Fault-System-Ruptures-Plugin");
		
		Element colorerEl = pluginEl.element("Colorer");
		
		Element cptEl = colorerEl.element("CPT");
		
		CPT cpt = CPT.fromXMLMetadata(cptEl);
		
		cpt.writeCPTFile(new File(dir, "cpt.cpt"));
		
		// make a fake plot
		
		int width = 400;
		int height = 1000;
		
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(10, 10, 0d, 0d, 0.1d);
		
		XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, " ", " ", " ", cptLabel);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(xyzSpec);
		
		PlotUtils.writePlots(dir, "cpt", gp, width, height, true, true, false);
		
		CPT logUpperCPT = new CPT();
		CPT logLowerCPT = new CPT();
		
		for (CPTVal val : cpt) {
			if (val.start == 0f || val.end == 0f)
				continue;
			if (val.start < 0f) {
				Preconditions.checkState(val.end < 0f);
				logLowerCPT.add(0, new CPTVal((float)Math.log10(-val.end), val.maxColor, (float)Math.log10(-val.start), val.minColor));
			} else {
				Preconditions.checkState(val.end > 0f);
				logUpperCPT.add(new CPTVal((float)Math.log10(val.start), val.minColor, (float)Math.log10(val.end), val.maxColor));
			}
		}
		
		logLowerCPT.add(0, new CPTVal(-5f, Color.WHITE, logLowerCPT.getMinValue(), logLowerCPT.getMinColor()));
		logUpperCPT.add(0, new CPTVal(-5f, Color.WHITE, logUpperCPT.getMinValue(), logUpperCPT.getMinColor()));
		
		System.out.println(logUpperCPT);
		System.out.println(logLowerCPT);
		
		xyzSpec = new XYZPlotSpec(xyz, logLowerCPT, " ", " ", " ", cptLabel);
		gp.drawGraphPanel(xyzSpec);
		PlotUtils.writePlots(dir, "cpt_log_lower", gp, width, height, true, true, false);
		
		xyzSpec = new XYZPlotSpec(xyz, logUpperCPT, " ", " ", " ", cptLabel);
		gp.drawGraphPanel(xyzSpec);
		PlotUtils.writePlots(dir, "cpt_log_upper", gp, width, height, true, true, false);
	}

}
