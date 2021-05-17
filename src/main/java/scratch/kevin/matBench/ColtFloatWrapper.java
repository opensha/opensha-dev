package scratch.kevin.matBench;

import org.opensha.commons.util.ClassUtils;

import cern.colt.matrix.tfloat.FloatMatrix2D;

public class ColtFloatWrapper extends AbstractBenchableMatrix2D {
	
	private FloatMatrix2D mat;
	
	public ColtFloatWrapper(FloatMatrix2D mat) {
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
		return (double)mat.get(row, col);
	}

	@Override
	public void set(int row, int col, double value) {
		mat.set(row, col, (float)value);
	}

	@Override
	protected BenchableMatrix2D doMult(BenchableMatrix2D mat) {
		return new ColtFloatWrapper(this.mat.zMult(((ColtFloatWrapper)mat).mat, null));
	}
	
	@Override
	public String toString() {
		return ClassUtils.getClassNameWithoutPackage(mat.getClass());
	}

}
