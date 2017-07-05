package scratch.kevin.matBench;

import org.opensha.commons.util.ClassUtils;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class ColtDoubleWrapper extends AbstractBenchableMatrix2D {
	
	private DoubleMatrix2D mat;
	
	public ColtDoubleWrapper(DoubleMatrix2D mat) {
		this.mat = mat;
	}

	@Override
	public int rows() {
		return mat.rows();
	}

	@Override
	public int cols() {
		return mat.columns();
	}

	@Override
	public double get(int row, int col) {
		return mat.get(row, col);
	}

	@Override
	public void set(int row, int col, double value) {
		mat.set(row, col, value);
	}

	@Override
	protected BenchableMatrix2D doMult(BenchableMatrix2D mat) {
		return new ColtDoubleWrapper(this.mat.zMult(((ColtDoubleWrapper)mat).mat, null));
	}
	
	@Override
	public String toString() {
		return ClassUtils.getClassNameWithoutPackage(mat.getClass());
	}

}
