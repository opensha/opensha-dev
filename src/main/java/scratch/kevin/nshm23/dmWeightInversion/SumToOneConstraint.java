package scratch.kevin.nshm23.dmWeightInversion;

import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class SumToOneConstraint extends InversionConstraint {

	protected SumToOneConstraint(double weight) {
		super("Sum to 1", "SumTo1", weight, false, ConstraintWeightingType.UNNORMALIZED);
	}

	@Override
	public int getNumRows() {
		return 1;
	}

	@Override
	public long encode(DoubleMatrix2D A, double[] d, int startRow) {
		double aScalar = getWeight()*getWeightingType().getA_Scalar(1d, Double.NaN);
		for (int i=0; i<A.columns(); i++)
			setA(A, startRow, i, aScalar);
		d[startRow] = getWeight()*getWeightingType().getD(1d, Double.NaN);
		return A.columns();
	}

}
