package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class BranchScratch {
	
	private static void test1() throws IOException {
		MeanUCERF2 u2 = new MeanUCERF2();
		u2.getParameter(UCERF2.PROB_MODEL_PARAM_NAME).setValue(UCERF2.PROB_MODEL_POISSON);
		u2.getTimeSpan().setDuration(1d);
		u2.updateForecast();
		
		int sourceID = 128;
		ProbEqkSource source = u2.getSource(sourceID);
		
		System.out.println("Source has "+source.getNumRuptures()+" ruptures");
		
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(5.05, 39, 0.1);
		int[] counts = new int[mfd.size()];
		for (ProbEqkRupture rup : source) {
			int index = mfd.getClosestXIndex(rup.getMag());
			mfd.add(index, rup.getMeanAnnualRate(1d));
			counts[index]++;
		}
		for (int i=0; i<mfd.size(); i++)
			System.out.println("M"+(float)mfd.getX(i)+"\t"+counts[i]+" rups, rate="+(float)mfd.getY(i));
		
		PlotSpec plot = new PlotSpec(List.of(mfd),
				List.of(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK)),
				"UCERF2 Source "+sourceID+": "+source.getName(), "Magnitude", "Annual Rate");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, true, new Range(6d, 9d), new Range(1e-6, 1e0));
		
		PlotUtils.writePlots(new File("/tmp"), "u2_source_"+sourceID, gp, 800, 800, true, false, false);
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		try {
			test1();
		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(1);
		}
	}

}
