package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import org.junit.Test;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import scratch.aftershockStatisticsETAS.griddedInterpGMPE.NDimArrayCalc;

public class NDimArrayCalcTest {
	
	private static Random r = new Random();
	private static Joiner join = Joiner.on(",");
	
	private static int[] randomDims(int nDims, int min, int max) {
		int delta = max - min;
		Preconditions.checkState(delta > 0);
		int[] dimensions = new int[nDims];
		for (int i=0; i<nDims; i++)
			dimensions[i] = min+r.nextInt(delta);
		return dimensions;
	}
	
	private static int[] randomIndex(int[] dimensions) {
		int[] index = new int[dimensions.length];
		for (int i=0; i<dimensions.length; i++)
			index[i] = r.nextInt(dimensions[i]);
		return index;
	}

	@Test
	public void testFirstIndexAgainstRegular() {
		int maxDim = 5;
		int numPerDim = 1000;
		for (int nDims=1; nDims<maxDim; nDims++) {
			int[] dimensions = randomDims(nDims, 1, 5);
			NDimArrayCalc calc = new NDimArrayCalc(dimensions);
			
			for (int i=0; i<numPerDim; i++) {
				int[] indexes = randomIndex(dimensions);
				int index1 = calc.getIndex(indexes);
				int index2 = calc.getIndex(indexes[0], Arrays.copyOfRange(indexes, 1, indexes.length));
				assertEquals("", index1, index2);
			}
		}
	}
	
	private static int[] nextIndex(int[] current, int[] dimensions) {
		int[] next = Arrays.copyOf(current, current.length);
		for (int i=next.length; --i>=0;) {
			if (next[i] < dimensions[i]-1) {
				next[i]++;
				for (int j=i+1; j<next.length; j++)
					next[j] = 0;
				return next;
			}
		}
		return null;
	}
	
	@Test
	public void testAllIndexesUsed() {
		int maxDim = 5;
		int numPerDim = 100;
		for (int nDims=1; nDims<maxDim; nDims++) {
			for (int i=0; i<numPerDim; i++) {
				int[] dimensions = randomDims(nDims, 1, 5);
				NDimArrayCalc calc = new NDimArrayCalc(dimensions);
				
//				System.out.println("Test for dims: ["+join.join(Ints.asList(dimensions))+"]");
				
				int numIndexes = calc.rawArraySize();
				
				HashSet<Integer> indexTracker = new HashSet<>();
				int[] indexes = new int[nDims];
				while (indexes != null) {
//					System.out.println(join.join(Ints.asList(indexes)));
					int index = calc.getIndex(indexes);
					assertTrue("Duplicate index found: "+index, indexTracker.add(index));
					indexes = nextIndex(indexes, dimensions);
				}
				assertEquals("Didn't process all indexes for dims: ["+join.join(Ints.asList(dimensions))+"]",
						numIndexes, indexTracker.size());
			}
		}
	}
	
	@Test
	public void testRawArraySize() {
		int maxDim = 5;
		int numPerDim = 100;
		for (int nDims=1; nDims<maxDim; nDims++) {
			for (int i=0; i<numPerDim; i++) {
				int[] dimensions = randomDims(nDims, 1, 5);
				NDimArrayCalc calc = new NDimArrayCalc(dimensions);
				int numIndexes = calc.rawArraySize();
				
				int calculated = 1;
				for (int size : dimensions)
					calculated *= size;
				assertEquals("Unexpected array size", calculated, numIndexes);
			}
		}
	}
	
	@Test
	public void testCollapsed() {
		int maxDim = 5;
		int numPerDim = 100;
		for (int nDims=2; nDims<maxDim; nDims++) {
			for (int i=0; i<numPerDim; i++) {
				int[] dimensions = randomDims(nDims, 2, 5);
				NDimArrayCalc calc = new NDimArrayCalc(dimensions);
				
				int index = r.nextInt(nDims);
				int value = r.nextInt(dimensions[index]);
				
				NDimArrayCalc collapsed = calc.getCollapsedView(index, value);
				assertEquals("Raw rray size shouldn't change when collapsing", calc.rawArraySize(), collapsed.rawArraySize());
				int[] indexes = randomIndex(dimensions);
				indexes[index] = value;
				
				int rawIndex = calc.getIndex(indexes);
				int[] collapsedIndexes = new int[indexes.length-1];
				for (int j=0; j<collapsedIndexes.length; j++) {
					if (j < index)
						collapsedIndexes[j] = indexes[j];
					else
						collapsedIndexes[j] = indexes[j+1];
				}
				int collapsedRawIndex = collapsed.getIndex(collapsedIndexes);
				assertEquals("Bad collapses index calc", rawIndex, collapsedRawIndex);
			}
		}
	}

	@Test
	public void testCollapsedFirstIndexAgainstRegular() {
		int maxDim = 5;
		int numPerDim = 1000;
		for (int nDims=2; nDims<maxDim; nDims++) {
			int[] dimensions = randomDims(nDims, 1, 5);
			NDimArrayCalc calc = new NDimArrayCalc(dimensions);
			
			int collapseIndex = r.nextInt(nDims);
			calc = calc.getCollapsedView(collapseIndex, r.nextInt(dimensions[collapseIndex]));
			dimensions = calc.getDimensions();
			
			for (int i=0; i<numPerDim; i++) {
				int[] indexes = randomIndex(dimensions);
				int index1 = calc.getIndex(indexes);
				int index2 = calc.getIndex(indexes[0], Arrays.copyOfRange(indexes, 1, indexes.length));
				assertEquals("", index1, index2);
			}
		}
	}

}
