package scratch.peter.nshmp;

import java.util.Arrays;
import java.util.List;

import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.nshmp2.erf.NSHMP2008;
import org.opensha.nshmp2.erf.source.FixedStrikeSource;
import org.opensha.nshmp2.erf.source.GridERF;
import org.opensha.nshmp2.erf.source.NSHMP_ERF;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkSource;

import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;

/**
 * 2008 NSHMP grid source generator that restricts sources to those that fall
 * within the CA state border.
 * 
 * @author Peter Powers
 */
public class NSHMP08_ClippedGridSourceGen {
	
	private static final GriddedRegion ca_region;
	
	private NSHMP2008 gridListERF;
	
	// array of summed erf source counts; facilitates indexing of nested erfs
	private int[] erfStartIndices;
	
	// lists of actual indices inside CA for each ERF
	private List<List<Integer>> erfIndices = Lists.newArrayList();
	
	// total number of sources
	private int srcCount = 0;
	
	public NSHMP08_ClippedGridSourceGen() {
		gridListERF = NSHMP2008.createCaliforniaGridded();

		
		List<Integer> indexList = Lists.newArrayList();
		for (ERF erf : gridListERF) {
			
			// need to include an initial 0 value and can skip last value
			indexList.add(srcCount);
			
			int count = 0; // valid index counter
			List<Integer> indices = Lists.newArrayList();
			erfIndices.add(indices);

			// this creates a lot of sources unnecessarily but should only 
			// be called once per ERF init
			for (ProbEqkSource src : erf) {
				if (src instanceof PointSource13b) {
					Location loc = ((PointSource13b) src).getLocation();
					if (ca_region.indexForLocation(loc) != -1) indices.add(count);
				} else if (src instanceof FixedStrikeSource) {
					Location loc = ((FixedStrikeSource) src).getLocation();
					if (ca_region.indexForLocation(loc) != -1) indices.add(count);
				} else {
					throw new IllegalStateException(
						"Illegal grid source: " + src.getClass());
				}
				count++;
			}
			srcCount += indices.size();
		}
		erfStartIndices = Ints.toArray(indexList);
		
		// NOTE:
		// adjust mfd rates of sources; this is possible becuase the mfds can
		// be accessed publicly without cloning (and therefore adjusted),
		// this isn't good
		
		for (NSHMP_ERF erf : gridListERF) {
			GridERF gerf = (GridERF) erf;
			gerf.scaleRatesToWeight();
		}
	}
	
	public ProbEqkSource getSource(int srcIdx) {
		int erfIdx = Arrays.binarySearch(erfStartIndices, srcIdx);
		// for index matches, select the next highest erf
		erfIdx = (erfIdx < 0) ? -(erfIdx + 2) : erfIdx;
		srcIdx = srcIdx - erfStartIndices[erfIdx];
		return gridListERF.getERF(erfIdx).getSource(srcIdx);
	}
	
	public int getNumSources() {
		return srcCount;
	}
	
	public void setForecastDuration(double duration) {
		gridListERF.getTimeSpan().setDuration(duration);
	}
	
	static {
		LocationList locs = new LocationList();
		locs.add(new Location(39.000, -119.999));
		locs.add(new Location(35.000, -114.635));
		locs.add(new Location(34.848, -114.616));
		locs.add(new Location(34.719, -114.482));
		locs.add(new Location(34.464, -114.371));
		locs.add(new Location(34.285, -114.122));
		locs.add(new Location(34.097, -114.413));
		locs.add(new Location(33.934, -114.519));
		locs.add(new Location(33.616, -114.511));
		locs.add(new Location(33.426, -114.636));
		locs.add(new Location(33.401, -114.710));
		locs.add(new Location(33.055, -114.676));
		locs.add(new Location(33.020, -114.501));
		locs.add(new Location(32.861, -114.455));
		locs.add(new Location(32.741, -114.575));
		locs.add(new Location(32.718, -114.719));
		locs.add(new Location(32.151, -120.861));
		locs.add(new Location(39.000, -126.000));
		locs.add(new Location(42.001, -126.000));
		locs.add(new Location(42.001, -119.999));
		locs.add(new Location(39.000, -119.999));
		ca_region = new GriddedRegion(locs, BorderType.MERCATOR_LINEAR, 0.1, GriddedRegion.ANCHOR_0_0);
	}

//	public static void main(String[] args) {
//		NSHMP08_ClippedGridSourceGen grdGen = new NSHMP08_ClippedGridSourceGen();
//		System.out.println(Arrays.toString(grdGen.erfStartIndices));
//		for (List<Integer> erfIndices : grdGen.erfIndices) {
//			System.out.println(erfIndices);
//		}
//		System.out.println(grdGen.srcCount);
//	}

}
