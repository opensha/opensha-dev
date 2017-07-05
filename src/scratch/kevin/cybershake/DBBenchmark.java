package scratch.kevin.cybershake;

import java.awt.Color;
import java.util.List;
import java.util.Random;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeIM.CyberShakeComponent;
import org.opensha.sha.cybershake.db.CybershakeIM.IMType;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.PeakAmplitudesFromDB;

import com.google.common.collect.Lists;

public class DBBenchmark {

	public static void main(String[] args) {
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		PeakAmplitudesFromDB amps2db = new PeakAmplitudesFromDB(db);
		
		CybershakeIM im = new CybershakeIM(146, IMType.SA, 3d, null, CyberShakeComponent.RotD100);
		
		List<Integer> runIDs = Lists.newArrayList(2657,3037,2722,3022,3030,3027,2636,2638,2660,2703,3504,2988,2965,3007);
		Random r = new Random();
		
		HistogramFunction hist = new HistogramFunction(5d, 20, 10d);
		
		for (int i=0; i<200; i++) {
			int runID = runIDs.get(r.nextInt(runIDs.size()));
			long curTime = System.currentTimeMillis();
			amps2db.countAmps(runID, im);
			long delta = System.currentTimeMillis()-curTime;
			// convert to seconds
			double secs = (double)delta/1000d;
			System.out.println("Took "+(float)secs+" for "+runID);
			hist.add(hist.getClosestXIndex(secs), 1d);
		}
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		PlotSpec spec = new PlotSpec(funcs, chars, "Count Amps time Distribution", "Time (s)", "Number");
		
		new GraphWindow(spec);
	}

}
