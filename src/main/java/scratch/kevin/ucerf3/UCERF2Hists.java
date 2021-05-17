package scratch.kevin.ucerf3;

import java.awt.Color;
import java.util.ArrayList;

import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.collect.Lists;

public class UCERF2Hists {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MeanUCERF2 ucerf2 = new MeanUCERF2();
		ucerf2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
		ucerf2.updateForecast();
		
		
		
		HistogramFunction magHist = new HistogramFunction(5.05, 36, 0.1);
		HistogramFunction lengthHist = new HistogramFunction(0d, 63, 20d);
		int num = 0;
		for (ProbEqkSource s : ucerf2) {
			for (ProbEqkRupture r : s) {
				magHist.add(r.getMag(), 1.0);
				lengthHist.add(r.getRuptureSurface().getAveLength(), 1.0);
				num++;
			}
		}
		magHist.normalizeBySumOfY_Vals();
		magHist.setName("Mag Histogram for UCERF2");
		magHist.setInfo("(based on "+num+" ruptures)");
		lengthHist.normalizeBySumOfY_Vals();
		lengthHist.setName("Length Histogram for UCERF2");
		lengthHist.setInfo("(based on "+num+" ruptures)");
		
		
		ArrayList<HistogramFunction> histList = new ArrayList<HistogramFunction>();
		histList.add(magHist);
		PlotCurveCharacterstics histChar = new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 0.1f, Color.BLUE);
		
		GraphWindow graph = new GraphWindow(histList, "Rupture Magnitude Histogram",Lists.newArrayList(histChar)); 
		graph.setX_AxisLabel("Magnitude");
		graph.setY_AxisLabel("Normalized Number");
		
		histList.clear();
		histList.add(lengthHist);
		
		graph = new GraphWindow(histList, "Ruptre Length Histogram",Lists.newArrayList(histChar)); 
		graph.setX_AxisLabel("Length");
		graph.setY_AxisLabel("Normalized Number");
	}

}
