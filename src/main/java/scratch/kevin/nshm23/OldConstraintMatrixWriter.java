package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;

import cern.colt.function.tdouble.IntIntDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class OldConstraintMatrixWriter {

	public static void main(String[] args) throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
				new File("/home/kevin/markdown/inversions/fm3_1_u3ref_uniform_reproduce_ucerf3.zip"));
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/experiments/2021_10-paleo-uncertainty-refactor-verification");
		
//		String prefix = "paleo_rate";
//		String suffix = "new"
//		PaleoRateInversionConstraint constraint = new PaleoRateInversionConstraint(rupSet, 1d,
//				UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList()),
//				UCERF3_PaleoProbabilityModel.load());
		
		String prefix = "paleo_slip";
		String suffix = "mod2";
		PaleoSlipInversionConstraint constraint = new PaleoSlipInversionConstraint(rupSet, 1d,
				U3AveSlipConstraint.load(rupSet.getFaultSectionDataList()), U3AveSlipConstraint.slip_prob_model, false);
		
		int numRups = rupSet.getNumRuptures();
		int numRows = constraint.getNumRows();
		
		System.out.println("A size: "+numRows+" x "+numRups+" = "+(numRows*numRups));
		
		DoubleMatrix2D A = new SparseDoubleMatrix2D(numRows, numRups);
		double[] d = new double[numRows];
		
		constraint.encode(A, d, 0);
		
		CSVFile<Number> aCSV = new CSVFile<>(true);
		A.forEachNonZero(new IntIntDoubleFunction() {
			
			@Override
			public double apply(int row, int col, double val) {
				aCSV.addLine(row, col, val);
				return val;
			}
		});
		System.out.println("A has "+aCSV.getNumRows()+" values");
		
		CSVFile<Double> dCSV = new CSVFile<>(true);
		for (double v : d)
			dCSV.addLine(v);
		
		aCSV.writeToFile(new File(outputDir, prefix+"_A_"+suffix+".csv"));
		dCSV.writeToFile(new File(outputDir, prefix+"_d_"+suffix+".csv"));
	}

}
