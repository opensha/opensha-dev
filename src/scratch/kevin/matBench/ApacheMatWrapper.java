package scratch.kevin.matBench;

import org.apache.commons.math3.linear.RealMatrix;
import org.opensha.commons.util.ClassUtils;

public class ApacheMatWrapper extends AbstractBenchableMatrix2D {
	
	private RealMatrix mat;
	
	public ApacheMatWrapper(RealMatrix mat) {
		this.mat = mat;
	}

	@Override
	public int rows() {
		return mat.getRowDimension();
	}

	@Override
	public int cols() {
		return mat.getColumnDimension();
	}

	@Override
	public double get(int row, int col) {
		return mat.getEntry(row, col);
	}

	@Override
	public void set(int row, int col, double value) {
		mat.setEntry(row, col, value);
	}

	@Override
	protected BenchableMatrix2D doMult(BenchableMatrix2D mat) {
		return new ApacheMatWrapper(this.mat.multiply(((ApacheMatWrapper)mat).mat));
	}
	
	@Override
	public String toString() {
		return ClassUtils.getClassNameWithoutPackage(mat.getClass());
	}

}
