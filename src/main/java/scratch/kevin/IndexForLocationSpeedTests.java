package scratch.kevin;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;

import com.google.common.base.Stopwatch;

public class IndexForLocationSpeedTests {

	public static void main(String[] args) {
		int num = 100000000;
		
		double spacing = 0.1/5d;
		GriddedRegion reg = new GriddedRegion(new CaliforniaRegions.RELM_COLLECTION(), spacing, null);
		
		double halfSpacing = reg.getLatSpacing()*0.5;
		
		Location[] testLocs = new Location[num];
		int[] nodes = new int[num];
		int[] results = new int[num];
		
		Random r = new Random(0);
		System.out.println("Preparing locations...");
		for (int i=0; i<num; i++) {
			nodes[i] = r.nextInt(reg.getNodeCount());
			testLocs[i] = reg.locationForIndex(nodes[i]);
			testLocs[i] = new Location(testLocs[i].getLatitude() + r.nextDouble()*spacing - halfSpacing,
					testLocs[i].getLongitude() + r.nextDouble()*spacing - halfSpacing);
		}
		
		System.out.println("Testing...");
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<num; i++)
			results[i] = reg.indexForLocation(testLocs[i]);
		watch.stop();
		System.out.println("Took "+(float)(watch.elapsed(TimeUnit.MILLISECONDS)/1000d)+" seconds");
		int numEqual = 0;
		for (int i=0; i<num; i++)
			if (nodes[i] == results[i])
				numEqual++;
		System.out.println(numEqual+"/"+num+" matched exactly");
	}

}
