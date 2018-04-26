package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import com.google.common.base.Preconditions;

public class NDimArrayCalc {
	
	private int[] dimensions;
	private int[] strides;
	private int size;
	
	public NDimArrayCalc(int[] dimensions) {
		Preconditions.checkArgument(dimensions.length > 0);
		strides = new int[dimensions.length];
		this.dimensions = dimensions;
		
		strides[0] = 1;
		size = dimensions[0];
		
		for (int i=1; i<dimensions.length; i++) {
			strides[i] = strides[i-1]*dimensions[i-1];
			size *= dimensions[i];
		}
	}
	
	private NDimArrayCalc(int[] dimensions, int[] strides, int size) {
		this.dimensions = dimensions;
		this.strides = strides;
		this.size = size;
	}
	
	public int getNumDimensions() {
		return dimensions.length;
	}
	
	public int[] getDimensions() {
		return dimensions;
	}
	
	public int getIndex(int firstIndex, int[] nextIndexes) {
		int index = firstIndex;
		for (int i=0; i<nextIndexes.length; i++)
			index += nextIndexes[i]*strides[i+1];
		return index;
	}
	
	public int getIndex(int... indexes) {
		int index = indexes[0];
		for (int i=1; i<indexes.length; i++)
			index += indexes[i]*strides[i];
		return index;
	}
	
	public int rawArraySize() {
		return size;
	}
	
	public NDimArrayCalc getCollapsedView(int index, int value) {
		Preconditions.checkState(getNumDimensions() > 1, "Can't collapse with only 1 dimension!");
		return new CollapsedNDimArrayCalc(this, index, value);
	}
	
	private static int[] collapseRemove(int[] orig, int index) {
		Preconditions.checkState(index < orig.length && index >= 0);
		int[] collapsed = new int[orig.length-1];
		for (int i=0; i<collapsed.length; i++) {
			if (i < index)
				collapsed[i] = orig[i];
			else
				collapsed[i] = orig[i+1];
		}
		return collapsed;
	}
	
	private static int[] expandAdd(int[] orig, int index, int value) {
		Preconditions.checkState(index <= orig.length && index >= 0, "Bad expand index with orig len=%s, index=%s, value=%s",
				orig.length, index, value);
		int[] expanded = new int[orig.length+1];
		for (int i=0; i<expanded.length; i++) {
			if (i == index)
				expanded[i] = value;
			else if (i < index)
				expanded[i] = orig[i];
			else
				expanded[i] = orig[i-1];
		}
		return expanded;
	}
	
	private class CollapsedNDimArrayCalc extends NDimArrayCalc {

		private NDimArrayCalc orig;
		private int index;
		private int value;

		public CollapsedNDimArrayCalc(NDimArrayCalc orig, int index, int value) {
			super(collapseRemove(orig.dimensions, index), null, orig.size);
			this.orig = orig;
			this.index = index;
			this.value = value;
		}

		@Override
		public int getIndex(int firstIndex, int[] nextIndexes) {
			if (this.index == 0)
				return orig.getIndex(value, expandAdd(nextIndexes, 0, firstIndex));
			return orig.getIndex(firstIndex, expandAdd(nextIndexes, index-1, value));
		}

		@Override
		public int getIndex(int... indexes) {
			return orig.getIndex(expandAdd(indexes, index, value));
		}
		
	}

}
