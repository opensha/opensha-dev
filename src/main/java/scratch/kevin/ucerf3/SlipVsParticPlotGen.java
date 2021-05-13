package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class SlipVsParticPlotGen {

	public static void main(String[] args) throws IOException, DocumentException {
		// TODO Auto-generated method stub
		File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
		
		
		File outputDir = new File("/tmp");
		double[] mags = {6d, 6.5d, 7d, 7.5d, 8d, 8.5d };
		
		for (int i=1; i<mags.length; i++) {
			double minMag = mags[i-1];
			double maxMag = mags[i];
			
			double[] particRates = sol.calcParticRateForAllSects(minMag, maxMag);
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			
			for (int j=0; j<particRates.length; j++)
				if (particRates[j] > 0)
					scatter.set(sol.getRupSet().getSlipRateForSection(j)*1000d, particRates[j]);
			
			List<XY_DataSet> elems = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			elems.add(scatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
			
			String xAxisLabel = "Slip Rate (mm/yr)";
			String yAxisLabel = "M"+(float)minMag+" => M"+(float)maxMag+" Partic Rate";
			
			PlotSpec spec = new PlotSpec(elems, chars, "Partic Rate vs Slip Rate", xAxisLabel, yAxisLabel);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setBackgroundColor(Color.WHITE);
			
			gp.setXLog(true);
			gp.setYLog(true);
//			gp.setUserBounds(1e-2, 1e0, 1e-7, 1e-2);
			gp.drawGraphPanel(spec);
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(new File(outputDir, "partic_vs_slip_"+(float)minMag+"_"+(float)maxMag+".png").getAbsolutePath());
		}
	}

}
