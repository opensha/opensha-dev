package scratch.kevin;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils;

public class DistSpeedTest {
	
	private Location startLoc;
	private LocationVector vector;
	private int iterations;
	
	private LocationList locs;
	private double[] fastDists;
	private double[] slowDists;
	private long[] fastTimes;
	private long[] slowTimes;
	
	private long fastMilis;
	private long slowMilis;
	
	public DistSpeedTest(Location startLoc, LocationVector vector, int iterations) {
		this.startLoc = startLoc;
		this.vector = vector;
		this.iterations = iterations;
	}
	
	public void doCalc() {
		if (locs == null) {
			Location testLoc = startLoc;
			locs = new LocationList();
			
			for (int i=0; i<iterations; i++) {
				testLoc = LocationUtils.location(testLoc, vector);
				locs.add(testLoc);
			}
		}
		
		slowDists = new double[locs.size()];
		slowTimes = new long[locs.size()];
		slowMilis = doTest(false);
		fastDists = new double[locs.size()];
		fastTimes = new long[locs.size()];
		fastMilis = doTest(true);
	}
	
	public void printReport() {
		double fastSecs = (double)fastMilis / 1000d;
		double slowSecs = (double)slowMilis / 1000d;
		double avgFast = (double)fastMilis / (double)locs.size();
		double avgSlow = (double)slowMilis / (double)locs.size();
		double speedup = slowSecs / fastSecs;
		System.out.println("TOTAL TIMES");
		System.out.println("Fast: " + fastMilis + " milis = " + fastSecs + " secs (" + avgFast + " milis/calc)");
		System.out.println("Slow: " + slowMilis + " milis = " + slowSecs + " secs (" + avgSlow + " milis/calc)");
		System.out.println("Speedup: " + speedup + "x");
		
		int targetNum = 100;
		
		int modulus = 1;
		if (iterations > targetNum)
			modulus = iterations / targetNum;
		
		System.out.println("Printing with modulus " + modulus);
		
		System.out.println();
		System.out.println("slow dist,slow time,fast dist,fast time,fast-slow,percent diff,start lat,start lon,end lat, end lon");
		for (int i=0; i<iterations; i++) {
			if ((i+1) % modulus != 0)
				continue;
			double fastDist = fastDists[i];
			double slowDist = slowDists[i];
			Location loc = locs.get(i);
			
			double diff = fastDist - slowDist;
			double pDiff = DataUtils.getPercentDiff(fastDist, slowDist);
			
			System.out.println(slowDist+","+slowTimes[i]+","+fastDist+","+fastTimes[i]+","+diff+","+pDiff+","
					+startLoc.getLatitude()+","+startLoc.getLongitude()+","+loc.getLatitude()+","+loc.getLongitude());
		}
	}
	
	private long doTest(boolean fast) {
		long start = System.currentTimeMillis();
		for (int i=0; i<locs.size(); i++) {
			if (fast) {
				fastDists[i] = LocationUtils.linearDistanceFast(startLoc, locs.get(i));
				fastTimes[i] = System.currentTimeMillis()-start;
			} else {
				slowDists[i] = LocationUtils.linearDistance(startLoc, locs.get(i));
				slowTimes[i] = System.currentTimeMillis()-start;
			}
		}
		long end = System.currentTimeMillis();
		return end - start;
	}
	
	public static void main(String[] args) {

		Location testLoc = new Location(34, -118);
		double azimuth = Math.toRadians(45);
		double horizDist = 0.001d;
		double vertDist = 0d;
		LocationVector vector = new LocationVector(azimuth, horizDist, vertDist);
		
		int iterations = 5000000;
		
		DistSpeedTest test = new DistSpeedTest(testLoc, vector, iterations);
		
		test.doCalc();
		test.doCalc();
		test.doCalc();
		test.printReport();
	}

}
