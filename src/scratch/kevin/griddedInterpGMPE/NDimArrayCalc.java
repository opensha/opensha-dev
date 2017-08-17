package scratch.kevin.griddedInterpGMPE;

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

}
