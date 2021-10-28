package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.AnnealingProgress;

import com.google.common.base.Preconditions;

public class CARC_InversionThreadTestPlotter {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/data/kevin/ucerf4/batch_inversions/2021_10_18-reproduce-ucerf3-ref_branch-uniform-thread_test");
		
		String prefix = "u3_reproduce_";
		int maxThreads = 20;
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(1d, maxThreads);
		
		List<XY_DataSet> timeFuncs = new ArrayList<>();
		List<XY_DataSet> iterFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int t=2; t<=maxThreads; t++) {
			File subDir = new File(mainDir, prefix+t+"_threads");
			Preconditions.checkState(subDir.exists());
			
			File solFile = new File(subDir, "solution.zip");
			if (!solFile.exists())
				continue;
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			
			AnnealingProgress progress = sol.requireModule(AnnealingProgress.class);
			
			Color color = cpt.getColor((float)t);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
			
			XY_DataSet timeFunc = new DefaultXY_DataSet();
			XY_DataSet iterFunc = new DefaultXY_DataSet();
			
			timeFunc.setName(t+" Threads");
			iterFunc.setName(t+" Threads");
			
			for (int i=0; i<progress.size(); i++) {
				long time = progress.getTime(i);
				double secs = time/1000d;
				double mins = secs/60d;
				long iters = progress.getIterations(i);
				double e = progress.getEnergies(i)[0];
				
				timeFunc.set(mins, e);
				iterFunc.set((double)iters, e);
			}
			
			timeFuncs.add(timeFunc);
			iterFuncs.add(iterFunc);
		}
		
		PlotSpec iterSpec = new PlotSpec(iterFuncs, chars, "Energy vs Iterations", "Iteration Count", "Energy");
		iterSpec.setLegendVisible(true);
		PlotSpec timeSpec = new PlotSpec(timeFuncs, chars, "Energy vs Time", "Inversion Time (m)", "Energy");
		timeSpec.setLegendVisible(true);
		
		Range xRange = null;
		Range yRange = new Range(65d, 100d);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		gp.drawGraphPanel(timeSpec, false, false, xRange, yRange);
		PlotUtils.writePlots(mainDir, "energy_vs_time", gp, 1000, 800, true, false, false);
		gp.drawGraphPanel(iterSpec, false, false, xRange, yRange);
		PlotUtils.writePlots(mainDir, "energy_vs_iters", gp, 1000, 800, true, false, false);
	}

}
