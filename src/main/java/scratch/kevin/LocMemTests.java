package scratch.kevin;

import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Stopwatch;

class LocMemTests {
	
	private static class LocWithExtra {
		final double latDeg;
		final double latRad;
		final double lonDeg;
		final double lonRad;
		final double depth;
		
		LocWithExtra() {
			latDeg = Math.random();
			latRad = Math.random();
			lonDeg = Math.random();
			lonRad = Math.random();
			depth = Math.random();
		}
	}
	
	private static class LocExtraWithArray {
		final double[] values;
		
		static final int latDegIndex = 0;
		static final int latRadIndex = 1;
		static final int lonDegIndex = 2;
		static final int lonRadIndex = 3;
		static final int depthIndex = 4;
		
		LocExtraWithArray() {
			values = new double[5];
			for (int i=0; i<values.length; i++)
				values[i] = Math.random();
		}
	}
	
	private static class LocBasic {
		final double latRad;
		final double lonRad;
		final double depth;
		
		LocBasic() {
			latRad = Math.random();
			lonRad = Math.random();
			depth = Math.random();
		}
	}
	
	private static class LocBasicWithArray {
		final double[] values;
		
		static final int latRadIndex = 0;
		static final int lonRadIndex = 1;
		static final int depthIndex = 2;
		
		LocBasicWithArray() {
			values = new double[3];
			for (int i=0; i<values.length; i++)
				values[i] = Math.random();
		}
	}
	
	private static class WrappedDegOnlyLoc extends Location {

		private double lat;
		private double lon;

		public WrappedDegOnlyLoc(double lat, double lon) {
			super(lat, lon);
			this.lat = lat;
			this.lon = lon;
		}

		@Override
		public double getLatRad() {
			return Math.toRadians(lat);
		}

		@Override
		public double getLonRad() {
			return Math.toRadians(lon);
		}
		
	}

	public static void main(String[] args) throws InterruptedException {
		int num = 100000000;
		Object[] locs = new Object[num];
		
		long initialMem = memUsage();
		System.out.println("Initial: "+initialMem);
		for (int i=0; i<num; i++)
//			locs[i] = new LocWithExtra();
			locs[i] = new Location(Math.random(), Math.random());
//			locs[i] = new LocExtraWithArray();
//			locs[i] = new LocBasic();
//			locs[i] = new LocBasicWithArray();
		long newMem = memUsage();
		System.out.println("New: "+newMem);
		long memUsed = newMem - initialMem;
		double memUsedMB = (double)memUsed/(double)(1024*1024);
		System.out.println("Used "+memUsed+" bytes = "+memUsedMB+" MB for "+num+" instances");
		double memEach = (double)memUsed/(double)num;
		System.out.println(ClassUtils.getClassNameWithoutPackage(locs[0].getClass())+": "+memEach+" each");
		
		double[] latDegs = new double[num];
		double[] lonDegs = new double[num];
		double[] latRads = new double[num];
		double[] lonRads = new double[num];
		for (int i=0; i<num; i++) {
			latDegs[i] = Math.random()*180d-90d;
			lonDegs[i] = Math.random()*360d;
		}
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<num; i++) {
			latRads[i] = Math.toRadians(latDegs[i]);
			lonRads[i] = Math.toRadians(lonDegs[i]);
		}
		watch.stop();
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Took "+secs+" s to convert "+num+" locs from degrees to radians");
		for (boolean degOnly : new boolean[] {false, true}) {
			System.out.println();
			if (degOnly)
				System.out.println("Using DEGREES-ONLY location class");
			else
				System.out.println("Using original location class");
			Location[] origLocs = new Location[num];
			for (int i=0; i<num; i++) {
				if (degOnly)
					origLocs[i] = new WrappedDegOnlyLoc(latDegs[i], lonDegs[i]);
				else
					origLocs[i] = new Location(latDegs[i], lonDegs[i]);
			}
			double sumDist = 0d;
			watch = Stopwatch.createStarted();
			for (int i=0; i<origLocs.length; i++)
				sumDist += LocationUtils.horzDistanceFast(origLocs[0], origLocs[1]);
			watch.stop();
			secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
			System.out.println("Took "+secs+" s to calculate "+num+" fast distances");
			sumDist = 0d;
			watch = Stopwatch.createStarted();
			for (int i=0; i<origLocs.length; i++)
				sumDist += LocationUtils.horzDistance(origLocs[0], origLocs[1]);
			watch.stop();
			secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
			System.out.println("Took "+secs+" s to calculate "+num+" slow distances");
			System.out.println(sumDist);
		}
	}
	
	private static long memUsage() throws InterruptedException {
		System.gc();
		Thread.sleep(1000);
		System.gc();
		Thread.sleep(1000);
		System.gc();
		Thread.sleep(1000);
		Runtime runtime = Runtime.getRuntime();
		long total = runtime.totalMemory();
		long free = runtime.freeMemory();
		long used = total - free;
		System.out.println(used+" = "+total+" - "+free);
		return used;
	}

}
