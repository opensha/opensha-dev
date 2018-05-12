package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

public class NDimensionalLinearInterpolation {
	
	private int numDimensions;
	
	private NDimArrayCalc[] arrayCalcs;
	private List<List<int[]>> subInterpIndexes;
	
	public NDimensionalLinearInterpolation(int numDimensions) {
		this.numDimensions = numDimensions;
		
		arrayCalcs = new NDimArrayCalc[numDimensions];
		subInterpIndexes = new ArrayList<>();
		int[] subDimensions = new int[numDimensions];
		for (int i=0; i<subDimensions.length; i++)
			subDimensions[i] = 2;
		for (int i=0; i<numDimensions; i++) {
			arrayCalcs[i] = new NDimArrayCalc(subDimensions);
			subDimensions = Arrays.copyOfRange(subDimensions, 1, subDimensions.length);
			if (i == 0) {
				subInterpIndexes.add(null);
			} else {
				subInterpIndexes.add(getSubInterpIndexes(numDimensions-(i), 2));
//				System.out.println("Sub interp index count for "+i+": "+subInterpIndexes.get(i).size());
			}
		}
	}
	
	private NDimensionalLinearInterpolation(int numDimensions, NDimArrayCalc[] arrayCalcs, List<List<int[]>> subInterpIndexes) {
		this.numDimensions = numDimensions;
		this.arrayCalcs = arrayCalcs;
		this.subInterpIndexes = subInterpIndexes;
	}
	
	private NDimensionalLinearInterpolation collapsed;
	
	private synchronized NDimensionalLinearInterpolation getCollapsedInterpolator() {
		Preconditions.checkState(numDimensions > 1, "Can't collapse with only 1 dimension!");
		if (collapsed == null) {
			collapsed = new NDimensionalLinearInterpolation(numDimensions-1,
					Arrays.copyOfRange(arrayCalcs, 1, arrayCalcs.length), subInterpIndexes.subList(1, subInterpIndexes.size()));
		}
		return collapsed;
	}
	
	public double interpolate(double[] flatDataArray, NDimArrayCalc arrayCalc, double[] indexes) {
		Preconditions.checkState(arrayCalc.getNumDimensions() == numDimensions);
		
		// check for any fixed values and collapse to avoid unnecessary interpolation
		if (indexes.length > 1) {
			for (int i=0; i<indexes.length; i++) {
				if ((int)indexes[i] == (float)indexes[i]) {
					// remove that index
					double[] collapsed = new double[indexes.length-1];
					for (int j=0; j<collapsed.length; j++) {
						if (j < i)
							collapsed[j] = indexes[j];
						else
							collapsed[j] = indexes[j+1];
					}
					NDimensionalLinearInterpolation cInterp = getCollapsedInterpolator();
					return cInterp.interpolate(flatDataArray, arrayCalc.getCollapsedView(i, (int)indexes[i]), collapsed);
				}
			}
		} else {
			// only one dimension
			if ((int)indexes[0] == (float)indexes[0]) {
				return flatDataArray[arrayCalc.getIndex((int)indexes[0])];
			}
		}
		
		int[] dimensions = arrayCalc.getDimensions();
		for (int j=0; j<dimensions.length; j++) {
			int dimension = dimensions[j];
			Preconditions.checkArgument(dimension > 1,
					"Only dimensions with size >1 are supported, encountered %s at dim %s/%s", dimension, j+1, dimensions.length);
		}
		
		// convert to array where each dimension has 2 values, and delta is fraction between those 2 values
		// handles edge cases where index is on the greater edge
		
		int[] baseIndexes = new int[indexes.length];
		double[] deltas = new double[indexes.length];
		for (int i=0; i<baseIndexes.length; i++) {
			baseIndexes[i] = (int)indexes[i];
			if (baseIndexes[i] == dimensions[i] - 1) {
				// we're at an edge, move index down and set delta to 1
				Preconditions.checkState(indexes[i] == dimensions[i] - 1,
						"Index is outside of range for dimension %s, size=%s, index=%s", i, dimensions[i], indexes[i]);
				baseIndexes[i]--;
				deltas[i] = 1d;
			} else {
				deltas[i] = indexes[i] - baseIndexes[i];
			}
			Preconditions.checkState(baseIndexes[i] < dimensions[i] - 1 && baseIndexes[i]  >= 0,
					"Bad index for dimension %s, index=%s = %s+%s with size %s", i, indexes[i], baseIndexes[i], deltas[i], dimensions[i]);
		}
		
		double[] subFlatArray = new double[arrayCalcs[0].rawArraySize()];
		
		if (numDimensions == 1) {
			subFlatArray[0] = flatDataArray[arrayCalc.getIndex(baseIndexes[0])];
			subFlatArray[1] = flatDataArray[arrayCalc.getIndex(baseIndexes[0]+1)];
		} else {
			for (int[] subIndexes : subInterpIndexes.get(1)) {
				int[] origNextIndexes = new int[subIndexes.length];
				for (int i=0; i<origNextIndexes.length; i++)
					origNextIndexes[i] = baseIndexes[i+1] + subIndexes[i];
				int origIndex0 = arrayCalc.getIndex(baseIndexes[0], origNextIndexes);
				int origIndex1 = arrayCalc.getIndex(baseIndexes[0]+1, origNextIndexes);
				int newIndex0 = arrayCalcs[0].getIndex(0, subIndexes);
				int newIndex1 = arrayCalcs[0].getIndex(1, subIndexes);
				subFlatArray[newIndex0] = flatDataArray[origIndex0];
				subFlatArray[newIndex1] = flatDataArray[origIndex1];
			}
		}
		
		return interpolateTrimmed(subFlatArray, deltas);
	}
	
	public double interpolateTrimmed(double[] flatDataArray, double[] deltas) {
		Preconditions.checkArgument(flatDataArray.length == arrayCalcs[0].rawArraySize(),
				"Array size inconsistent, should be already trimmed such that each dimension is of size=2");
		for (int i=0; i<deltas.length; i++)
			Preconditions.checkState(deltas[i] >= 0d && deltas[i] <= 1d, "deltas values not in range [0 1]: %s", deltas[i]);
		
		return interpRecursive(flatDataArray, deltas, 0);
	}
	
	private double interpRecursive(double[] data, double[] deltas, int interpIndex) {
		if (deltas.length == 1) {
			Preconditions.checkState(data.length == 2);
			return (1d - deltas[0])*data[0] + deltas[0]*data[1];
		}
		
		NDimArrayCalc curArrayCalc = arrayCalcs[interpIndex];
		NDimArrayCalc nextArrayCalc = arrayCalcs[interpIndex+1];
		
		double myDelta = deltas[0];
		double[] subData = new double[nextArrayCalc.rawArraySize()];
		for (int[] indexes : subInterpIndexes.get(interpIndex+1)) {
			double v0 = data[curArrayCalc.getIndex(0, indexes)];
			double v1 = data[curArrayCalc.getIndex(1, indexes)];
			subData[nextArrayCalc.getIndex(indexes)] = (1d-myDelta)*v0 + myDelta*v1;
		}
		
		double[] subDeltas = Arrays.copyOfRange(deltas, 1, deltas.length);
		
		return interpRecursive(subData, subDeltas, interpIndex+1);
	}
	
	static List<int[]> getSubInterpIndexes(int numDimensions, int numEach) {
		int[] dimensions = new int[numDimensions];
		for (int i=0; i<dimensions.length; i++)
			dimensions[i] = numEach;
		return getSubInterpIndexes(new ArrayList<>(), new int[0], 0, dimensions.length, dimensions);
	}
	
	static List<int[]> getSubInterpIndexes(int[] dimensions) {
		return getSubInterpIndexes(new ArrayList<>(), new int[0], 0, dimensions.length, dimensions);
	}
	
	private static List<int[]> getSubInterpIndexes(List<int[]> indexes, int[] parent, int startDimension, int numDimensions, int[] dimensions) {
		int[] nextDimensions = Arrays.copyOfRange(dimensions, 1, dimensions.length);
		for (int i=0; i<dimensions[0]; i++) {
			int[] child = Arrays.copyOf(parent, parent.length+1);
			child[child.length-1] = i;
			
			if (child.length == numDimensions)
				indexes.add(child);
			else
				indexes.addAll(getSubInterpIndexes(indexes, child, startDimension+1, numDimensions, nextDimensions));
		}
		return indexes;
	}
	
//	public static void main(String[] args) {
//		// 1d case
//		int[] dimensions = { 2 };
//		double[] data = { 0d, 10d };
//		double[] deltas = { 0.5 };
//		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation( dimensions.length);
//		System.out.println(interp.interpolateTrimmed(data, deltas));
//		
//		// 2d case
//		dimensions = new int[] { 2, 2 };
//		NDimArrayCalc arrayCalc = new NDimArrayCalc(dimensions);
//		data = new double[arrayCalc.rawArraySize()];
//		for (int x=0; x<2; x++)
//			for (int y=0; y<2; y++)
//				data[arrayCalc.getIndex(x,y)] = x*y;
////		deltas = new double[] { 0.0, 0.0 };
//		deltas = new double[] { 1.0, 1.0 };
//		interp = new NDimensionalLinearInterpolation(dimensions.length);
//		System.out.println(interp.interpolateTrimmed(data, deltas));
//		
//		// speed tests
//		int numPerDim = 20;
//		// 2D
//		arrayCalc = new NDimArrayCalc(new int[] {numPerDim, numPerDim});
//		data = new double[arrayCalc.rawArraySize()];
//		for (int x=0; x<numPerDim; x++)
//			for (int y=0; y<numPerDim; y++)
//				data[arrayCalc.getIndex(x, y)] = Math.random();
//		speedTest(data, arrayCalc, (int)4e7, numPerDim);
//		// 3D
//		arrayCalc = new NDimArrayCalc(new int[] {numPerDim, numPerDim, numPerDim});
//		data = new double[arrayCalc.rawArraySize()];
//		for (int x=0; x<numPerDim; x++)
//			for (int y=0; y<numPerDim; y++)
//				for (int z=0; z<numPerDim; z++)
//					data[arrayCalc.getIndex(x, y, z)] = Math.random();
//		speedTest(data, arrayCalc, (int)1e7, numPerDim);
//		// 4D
//		arrayCalc = new NDimArrayCalc(new int[] {numPerDim, numPerDim, numPerDim, numPerDim});
//		data = new double[arrayCalc.rawArraySize()];
//		for (int x=0; x<numPerDim; x++)
//			for (int y=0; y<numPerDim; y++)
//				for (int z=0; z<numPerDim; z++)
//					for (int t=0; t<numPerDim; t++)
//						data[arrayCalc.getIndex(x, y, z, t)] = Math.random();
//		speedTest(data, arrayCalc, (int)1e6, numPerDim);
//		// 5D
//		arrayCalc = new NDimArrayCalc(new int[] {numPerDim, numPerDim, numPerDim, numPerDim, numPerDim});
//		data = new double[arrayCalc.rawArraySize()];
//		for (int x=0; x<numPerDim; x++)
//			for (int y=0; y<numPerDim; y++)
//				for (int z=0; z<numPerDim; z++)
//					for (int t=0; t<numPerDim; t++)
//						for (int p=0; p<numPerDim; p++)
//							data[arrayCalc.getIndex(x, y, z, t, p)] = Math.random();
//		speedTest(data, arrayCalc, (int)1e5, numPerDim);
//	}
//	
//	private static void speedTest(double[] flatData, NDimArrayCalc arrayCalc, int numTests, int numPerDim) {
//		// prepare search indexes first
//		System.out.println("Preparing for "+arrayCalc.getNumDimensions()+"-D test");
//		List<double[]> searchIndexes = new ArrayList<>();
//		for (int i=0; i<numTests; i++) {
//			double[] indexes = new double[arrayCalc.getNumDimensions()];
//			for (int j=0; j<indexes.length; j++)
//				indexes[j] = Math.random()*(numPerDim - 1);
//			searchIndexes.add(indexes);
//		}
//		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation(arrayCalc.getNumDimensions());
//		
//		System.out.println("Doing "+arrayCalc.getNumDimensions()+"-D test");
//		Stopwatch watch = Stopwatch.createStarted();
//		for (double[] indexes : searchIndexes)
//			interp.interpolate(flatData, arrayCalc, indexes);
//		watch.stop();
//		long millis = watch.elapsed(TimeUnit.MILLISECONDS);
//		double secs = millis/1000d;
//		System.out.println("\tTime to "+numTests+": "+(float)secs+"s");
//		double calcsPerSec = numTests/secs;
//		System.out.println("\tInterp/s: "+(float)calcsPerSec+"s");
//	}

}
