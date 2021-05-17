package scratch.kevin;

import java.util.ArrayList;

import org.opensha.commons.geo.Location;
import org.opensha.sha.faultSurface.FaultTrace;

public class StrikeDirTest {
	
	public static FaultTrace getNorthLineTrace(int segs) {
		FaultTrace trace = new FaultTrace("north line");
		
		for (int i=0; i<=segs; i++)
			trace.add(new Location(i, 0.0));
		
		return trace;
	}
	
	public static FaultTrace getEastLineTrace(int segs) {
		FaultTrace trace = new FaultTrace("east line");
		
		for (int i=0; i<=segs; i++)
			trace.add(new Location(0.0, i));
		
		return trace;
	}
	
	public static FaultTrace getRightAngle() {
		FaultTrace trace = new FaultTrace("right angle");
		
		trace.add(new Location(0.0, 0.0));
		trace.add(new Location(1.0, 0.0));
		trace.add(new Location(1.0, 1.0));
		
		return trace;
	}
	
	public static FaultTrace getZigZagTrace(double width, double skew) {
		FaultTrace trace = new FaultTrace("zig zag, width="+width+", skew="+skew);
		
		trace.add(new Location(0.0, 0.0));
		trace.add(new Location(1.0, width+skew));
		trace.add(new Location(2.0, -width));
		trace.add(new Location(3.0, 0.0));
		
		return trace;
	}
	
	private static void addCircle(FaultTrace trace, double radius, double xMax) {
		double radSquared = radius*radius;
		for (double x=0.0; x<=xMax; x+=0.1) {
			double y = Math.sqrt(Math.abs(radSquared - Math.pow(x, 2)));
//			System.out.println(x+", "+y);
			trace.add(new Location(x, y));
		}
	}
	
	public static FaultTrace getQuarterCircleTrace() {
		FaultTrace trace = new FaultTrace("quarter circle");
		
		double radius = 10.0;
		
		addCircle(trace, radius, radius);
		
		return trace;
	}
	
	public static FaultTrace getHalfCircleTrace() {
		FaultTrace trace = new FaultTrace("half circle");
		
		double radius = 10.0;
		
		addCircle(trace, radius, radius*2+0.05);
		
		return trace;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int calcTries = 100000;
		
		ArrayList<FaultTrace> traces = new ArrayList<FaultTrace>();
		
		traces.add(getNorthLineTrace(1));
		traces.add(getNorthLineTrace(20));
		traces.add(getEastLineTrace(1));
		traces.add(getEastLineTrace(20));
		traces.add(getRightAngle());
		traces.add(getZigZagTrace(10.0, 0.0));
		traces.add(getZigZagTrace(10.0, 10.0));
		traces.add(getQuarterCircleTrace());
		traces.add(getHalfCircleTrace());
		
		for (FaultTrace trace : traces) {
			double strikeDir = trace.getStrikeDirection();
			long start = System.currentTimeMillis();
			for (int i=0; i<calcTries; i++) {
				trace.getStrikeDirection();
			}
			long strikeDirMilis = System.currentTimeMillis()-start;
			
			double aveStrike = trace.getAveStrike();
			start = System.currentTimeMillis();
			for (int i=0; i<calcTries; i++) {
				trace.getAveStrike();
			}
			long aveStrikeMilis = System.currentTimeMillis()-start;
			
			System.out.println("---------------------------------------------------");
			System.out.println("Trace:\t\t\t\t" + trace.getName());
			System.out.println("Trace Points:\t\t\t" + trace.getNumLocations());
			System.out.println("getStrikeDirection() time:\t" + strikeDirMilis);
			System.out.println("getAveStrike() time:\t\t" + aveStrikeMilis);
			System.out.println("getStrikeDirection() val:\t" + strikeDir);
			System.out.println("getAveStrike() val:\t\t" + aveStrike);
			System.out.println("---------------------------------------------------");
		}
	}

}
