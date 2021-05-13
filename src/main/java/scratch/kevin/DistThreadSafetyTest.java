package scratch.kevin;

import java.util.ArrayList;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.RuptureSurface;

public class DistThreadSafetyTest {
	
	private static class TestRunnable implements Runnable {
		private double expectedJB;
		private RuptureSurface surf;
		private Location loc;
		public TestRunnable(double expectedJB, RuptureSurface surf, Location loc) {
			this.expectedJB = expectedJB;
			this.surf = surf;
			this.loc = loc;
		}
		
		public void run() {
			long cnt = 0;
			while (true) {
				cnt++;
				double val = surf.getDistanceJB(loc);
				if (val != expectedJB) {
					System.out.println("AHA it happened after "+cnt+" gets! expected: "+expectedJB+", got: "+val);
					System.exit(1);
				}
			}
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ERF erf = new MeanUCERF2();
		erf.updateForecast();
		
		ProbEqkRupture rup = erf.getRupture(128, 1000);
		Location loc1 = new Location(34, -118);
		Location loc2 = new Location(36, -120);
		
		RuptureSurface surf = rup.getRuptureSurface();
		
		System.out.println("Surf type: "+surf.getClass().getName());
		
		double jbDist1 = surf.getDistanceJB(loc1);
		System.out.println("dist for loc1: "+jbDist1);
		double jbDist2 = surf.getDistanceJB(loc2);
		System.out.println("dist for loc2: "+jbDist2);
		
		int num = 12;
		
		for (int i=0; i<num; i++) {
			TestRunnable run;
			if (i % 2 == 0)
				run = new TestRunnable(jbDist1, surf, loc1);
			else
				run = new TestRunnable(jbDist2, surf, loc2);
			Thread thread = new Thread(run);
			thread.start();
		}
	}

}
