package scratch.kevin.pointSources;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Named;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class SuperSamplingDistExperiments {

	public static void main(String[] args) {
//		Location centerLoc = new Location(0d, 0d);
		Location centerLoc = new Location(37d, -110d);
		
		double discr = 0.1;
		Region cell = new Region(new Location(centerLoc.getLatitude()-0.5*discr, centerLoc.getLongitude()-0.5*discr),
			new Location(centerLoc.getLatitude()+0.5*discr, centerLoc.getLongitude()+0.5*discr));
//		Region cell = new Region(centerLoc, 11d);
		
		List<TestCase> testCases = new ArrayList<>();
		
		GriddedRegion superDuperSampled = new GriddedRegion(cell, 0.0001, new Location(0.00005, 0.00005));
		System.out.println("SuperDuperSampled has "+superDuperSampled.getNodeCount()+" nodes");
		TestCase refCase = new LocListTestCase(superDuperSampled.getNodeList(), "Ideal sample");
		
		testCases.add(new GriddedRegionTestCase(cell, 21));
		testCases.add(new GriddedRegionTestCase(cell, 20));
		testCases.add(new GriddedRegionTestCase(cell, 12));
		testCases.add(new GriddedRegionTestCase(cell, 11));
		testCases.add(new GriddedRegionTestCase(cell, 10));
		testCases.add(new GriddedRegionTestCase(cell, 6));
		testCases.add(new GriddedRegionTestCase(cell, 5));
		
		testCases.add(new ExteriorNodesTestCase(cell, 12));
		testCases.add(new ExteriorNodesTestCase(cell, 11));
		testCases.add(new ExteriorNodesTestCase(cell, 10));
		testCases.add(new ExteriorNodesTestCase(cell, 6));
		testCases.add(new ExteriorNodesTestCase(cell, 5));
		
		LocationList border = cell.getBorder();
		if (LocationUtils.areSimilar(border.first(), border.last())) {
			// remove the closing point (don't want any duplicates)
			border = new LocationList(border);
			border.remove(border.size()-1);
		}
		
		testCases.add(new LocListTestCase(border, "Corners"));
		
		testCases.add(new ResampledBorderTestCase(border, 1d));
		testCases.add(new ResampledBorderTestCase(border, 2d));
		testCases.add(new ResampledBorderTestCase(border, 5d));
		
		LocationList sideCenterBorder = new LocationList();
		sideCenterBorder.add(new Location(centerLoc.getLatitude(), cell.getMaxLon()));
		sideCenterBorder.add(new Location(centerLoc.getLatitude(), cell.getMinLon()));
		sideCenterBorder.add(new Location(cell.getMaxLat(), centerLoc.getLongitude()));
		sideCenterBorder.add(new Location(cell.getMinLat(), centerLoc.getLongitude()));
		testCases.add(new LocListTestCase(sideCenterBorder, "Side centers"));
		LocationList centersAndCorners = new LocationList();
		centersAndCorners.addAll(sideCenterBorder);
		centersAndCorners.addAll(border);
		testCases.add(new LocListTestCase(centersAndCorners, "Corners and centers"));
		
		LocationList centerAndCorners = new LocationList();
		centerAndCorners.add(centerLoc);
		centerAndCorners.addAll(border);
		testCases.add(new LocListTestCase(centerAndCorners, "Corners and center"));
		
		LocationList centerOnly = new LocationList();
		centerOnly.add(centerLoc);
		testCases.add(new LocListTestCase(centerOnly, "Center only"));
		
		double[] dists = {
			0d, 10d, 20d, 30d, 40d, 50d, 75d, 100d, 150d, 200d, 300d
		};
		
		int longestName = 0;
		for (TestCase test : testCases)
			longestName = Integer.max(longestName, test.getName().length());
		
		for (double dist : dists) {
			System.out.println("Center Distance: "+(float)dist);
			int numTests = dist == 0d ? 1 : 50;
			Location[] testLocs = new Location[numTests];
			if (numTests == 1) {
				testLocs[0] = centerLoc;
			} else {
				for (int i=0; i<numTests; i++) {
					double az = i*Math.PI/(double)numTests;
					testLocs[i] = LocationUtils.location(centerLoc, az, dist);
				}
			}

			MinMaxAveTracker refTrack = refCase.run(centerLoc, testLocs);
			System.out.println("\t"+paddedName(refCase.getName()+":", longestName+1)+"\t"+trackString(refTrack, null));
			
			for (TestCase testCase : testCases) {
				MinMaxAveTracker track = testCase.run(centerLoc, testLocs);
				System.out.println("\t"+paddedName(testCase.getName()+":", longestName+1)+"\t"+trackString(track, refTrack));
			}
			System.out.println();
		}
	}
	
	private static String paddedName(String name, int longest) {
		while (name.length() < longest)
			name += " ";
		return name;
	}
	
	private static final DecimalFormat distDF = new DecimalFormat("0.00");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static String trackString(MinMaxAveTracker track, MinMaxAveTracker ref) {
		String ret = "["+distDF.format(track.getMin())+", "+distDF.format(track.getMax())+"]; avg="+distDF.format(track.getAverage());
		if (ref != null) {
			double diff = track.getAverage() - ref.getAverage();
			ret += " (";
			if (diff >= 0d)
				ret += "+";
			ret += pDF.format(diff/ref.getAverage())+")";
		}
		return ret;
	}
	
	private static class LocListTestCase implements TestCase {
		
		protected LocationList locs;
		private String name;

		public LocListTestCase(LocationList locs, String name) {
			this.locs = locs;
			this.name = name;
		}

		@Override
		public String getName() {
			return name+" ["+locs.size()+"]";
		}

		@Override
		public MinMaxAveTracker run(Location centerLoc, Location[] testLocs) {
			MinMaxAveTracker track = new MinMaxAveTracker();
			for (Location testLoc : testLocs)
				for (Location loc : locs)
					track.addValue(LocationUtils.linearDistanceFast(testLoc, loc));
			return track;
		}
		
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static class ResampledBorderTestCase extends LocListTestCase {

		public ResampledBorderTestCase(LocationList border, double discrKM) {
			super(null, oDF.format(discrKM)+"km border");
			if (LocationUtils.areSimilar(border.first(), border.last())) {
				// remove the closing point (don't want any duplicates)
				border = new LocationList(border);
				border.remove(border.size()-1);
			}
			FaultTrace borderAsTrace = new FaultTrace(null);
			borderAsTrace.addAll(border);
			borderAsTrace.add(border.first()); // need to close it, will remove at the end
			int samples = (int)Math.round(borderAsTrace.getTraceLength()/discrKM);
			locs = FaultUtils.resampleTrace(borderAsTrace, samples);
			locs.remove(locs.size()-1); // last is a duplicate
		}
		
	}
	
	private static class GriddedRegionTestCase extends LocListTestCase {

		public GriddedRegionTestCase(Region cell, int numSamples) {
			super(getSampled(cell, numSamples).getNodeList(), numSamples+"x sampled");
		}
		
	}
	
	private static class ExteriorNodesTestCase extends LocListTestCase {

		public ExteriorNodesTestCase(Region cell, int numSamples) {
			super(null, numSamples+"x exterior");
			GriddedRegion sampled = getSampled(cell, numSamples);
			int numLats = sampled.getNumLatNodes();
			int[] minLonIndexes = new int[numLats];
			int[] maxLonIndexes = new int[numLats];
			for (int i=0; i<numLats; i++)
				minLonIndexes[i] = Integer.MAX_VALUE;
			int minLatIndex = Integer.MAX_VALUE;
			int maxLatIndex = -1;
			LocationList nodes = sampled.getNodeList();
			for (Location loc : nodes) {
				int latIndex = sampled.getLatIndex(loc);
				int lonIndex = sampled.getLonIndex(loc);
				Preconditions.checkState(latIndex >= 0 && lonIndex >= 0);
				// overall min
				minLatIndex = Integer.min(latIndex, minLatIndex);
				maxLatIndex = Integer.max(latIndex, maxLatIndex);
				// min/max lon for this lat
				minLonIndexes[latIndex] = Integer.min(minLonIndexes[latIndex], lonIndex);
				maxLonIndexes[latIndex] = Integer.max(maxLonIndexes[latIndex], lonIndex);
			}
			System.out.println("LatIndex range: ["+minLatIndex+", "+maxLatIndex+"] for numLat="+numLats);
			this.locs = new LocationList();
			for (Location loc : nodes) {
				int latIndex = sampled.getLatIndex(loc);
				int lonIndex = sampled.getLonIndex(loc);
				if (latIndex == minLatIndex || latIndex == maxLatIndex
						|| lonIndex == minLonIndexes[latIndex] || lonIndex == maxLonIndexes[latIndex])
					locs.add(loc);
			}
			System.out.println("Kept "+locs.size()+" exterior nodes (of "+sampled.getNodeCount()+" nodes) for "+numSamples+"x sampled");
		}
		
	}
	
	private static GriddedRegion getSampled(Region cell, int numSamples) {
		double discrX = cell.getMaxLon() - cell.getMinLon();
		double discrY = cell.getMaxLat() - cell.getMinLat();
		return new GriddedRegion(cell, discrY/(double)numSamples, discrX/(double)numSamples,
				new Location(discrY/(2d*numSamples), discrX/(2d*numSamples)));
	}
	
	private static interface TestCase extends Named {
		
		public MinMaxAveTracker run(Location centerLoc, Location[] testLocs);
	}

}
