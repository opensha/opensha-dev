package scratch.kevin.markov;

import java.util.HashMap;
import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * This is an n-dimensional dataset without any set limits on dataset or dimensional size. Only non zero values
 * are stored. Solution space is assumed to be all zero before any set operations.
 * @author kevin
 *
 */
public class SparseNDimensionalHashDataset<E> {
	
	private static double[] asArray(double val, int num) {
		double[] ret = new double[num];
		for (int i=0; i<num; i++)
			ret[i] = val;
		return ret;
	}
	
	private static final double PRECISION_SCALE = 1 + 1e-14;
	
	private int nDims;
	private double[] zeroIndLocPerDim;
	private double[] spacingsPerDim;
	
	private HashMap<IndicesKey, E> data;
	
	public SparseNDimensionalHashDataset(int nDims, double minValsPerDim, double spacing) {
		this(nDims, asArray(minValsPerDim, nDims), asArray(spacing, nDims));
	}
	
	public SparseNDimensionalHashDataset(int nDims, double[] zeroIndLocPerDim, double[] spacingsPerDim) {
		this.nDims = nDims;
		this.zeroIndLocPerDim = zeroIndLocPerDim;
		this.spacingsPerDim = spacingsPerDim;
		
		this.data = Maps.newHashMap();
	}
	
	public int getNDims() {
		return nDims;
	}
	
	public E get(int[] indices) {
		Preconditions.checkArgument(indices.length == this.nDims);
		E val = data.get(new IndicesKey(indices));
		return val;
	}
	
	public void set(int[] indices, E val) {
		Preconditions.checkArgument(indices.length == this.nDims);
		data.put(new IndicesKey(indices), val);
	}
	
	public int indexForDimVal(int nDim, double val) {
		double iVal = PRECISION_SCALE * (val - zeroIndLocPerDim[nDim]) / spacingsPerDim[nDim];
		int i = (spacingsPerDim[nDim] == 0) ? 0 : (int) Math.round(iVal);
		return i;
	}
	
	public int[] indexesForDimVals(double[] vals) {
		Preconditions.checkArgument(vals.length == this.nDims);
		int[] indexes = new int[nDims];
		for (int i=0; i<nDims; i++)
			indexes[i] = indexForDimVal(i, vals[i]);
		return indexes;
	}
	
	public List<int[]> getPopulatedIndices() {
		List<int[]> indicesList = Lists.newArrayList();
		for (IndicesKey key : data.keySet())
			indicesList.add(key.indices);
		return indicesList;
	}

}
