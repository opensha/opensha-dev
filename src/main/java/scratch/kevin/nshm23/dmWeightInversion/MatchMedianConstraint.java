package scratch.kevin.nshm23.dmWeightInversion;

import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;

import com.google.common.base.Preconditions;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class MatchMedianConstraint extends InversionConstraint {

	private List<double[]> dms;
	private double[] medians;

	protected MatchMedianConstraint(List<double[]> dms, double[] medians, double weight, ConstraintWeightingType weightType) {
		super("Match Median Slip Rates", "MatchMed", weight, false, weightType);
		this.dms = dms;
		this.medians = medians;
	}

	@Override
	public int getNumRows() {
		return medians.length;
	}

	@Override
	public long encode(DoubleMatrix2D A, double[] d, int startRow) {
		ConstraintWeightingType type = getWeightingType();
		Preconditions.checkState(A.columns() == dms.size());
		
		double weight = getWeight();
		
		for (int i=0; i<medians.length; i++) {
			type.getA_Scalar(medians[i], Double.NaN);
			for (int m=0; m<dms.size(); m++) {
				double modelSlip = dms.get(m)[i];
				setA(A, startRow+i, m, modelSlip*weight);
			}
			d[startRow+i] = weight*type.getD(medians[i], Double.NaN);
		}
		return medians.length*dms.size();
	}

}
