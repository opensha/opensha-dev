package scratch.kevin;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.rupForecastImpl.Frankel96.Frankel96_AdjustableEqkRupForecast;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

public class SampedeKNLThreadTest {
	
	private static final String imt = PGA_Param.NAME;
	private static final DiscretizedFunc xVals = IMT_Info.getUSGS_PGA_Function();
	private static DiscretizedFunc logXVals;
	static {
		logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			logXVals.set(Math.log(xVals.getX(i)), 1d);
	}
	
	public static void main(String args[]) {
		int numEach = 5000;
		int[] threads = { 68/4, 68/2, 68, 68*2, 68*3, 68*4 };
		
		AbstractERF erf = new Frankel96_AdjustableEqkRupForecast();
		erf.updateForecast();
		AttenRelRef gmpeRef = AttenRelRef.BSSA_2014;
		
		List<Task> tasks = Lists.newArrayList();
		System.out.println("Building tasks");
		for (int i=0; i<numEach; i++)
			tasks.add(new Task(erf, gmpeRef));
		
		for (int num : threads) {
			System.gc();
			System.out.println("Calculating with "+num+" threads");
			ExecutorService exec = Executors.newFixedThreadPool(num);
			
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			Stopwatch watch = Stopwatch.createStarted();
			List<Future<?>> futures = Lists.newArrayList();
			for (Task task : tasks)
				futures.add(exec.submit(task));
			exec.shutdown();
			for (Future<?> future : futures) {
				try {
					future.get();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (ExecutionException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			watch.stop();
			
			double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
			System.out.println("Took "+secs+" secs");
			double rate = (double)numEach/(double)secs;
			System.out.println("Rate: "+(float)rate+" tasks/sec");
		}
	}
	
	private static class Task implements Runnable {
		
		private AbstractERF erf;
		private HazardCurveCalculator calc;
		private AttenuationRelationship gmpe;
		private Site site;

		private Task(AbstractERF erf, AttenRelRef gmpeRef) {
			this.erf = erf;
			calc = new HazardCurveCalculator();
			gmpe = gmpeRef.instance(null);
			gmpe.setParamDefaults();
			gmpe.setIntensityMeasure(imt);
			site = new Site();
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter(param);
		}

		@Override
		public void run() {
			double lat = 35 + 0.05*Math.random();
			double lon = -120 + 0.05*Math.random();
			Location loc = new Location(lat, lon);
			site.setLocation(loc);
			site.getParameter(Double.class, Vs30_Param.NAME).setValue(new Double(300d + 300d*Math.random()));
			calc.getHazardCurve(logXVals.deepClone(), site, gmpe, erf);
		}
		
	}

}
