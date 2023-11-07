package scratch.kevin.nshm23.dmCovarianceTests;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import Jama.Matrix;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;

class DecompBenchmark {

	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		
		// disable the upper bound
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = 0d;
		
		List<? extends FaultSection> subSects = dm.build(fm);
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(subSects);
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		Preconditions.checkState(subSects.size() == sol.getRupSet().getNumSections());
		SectionCovarianceSampler sampler = new PrecomputedConnectivityCorrelationSampler(
				subSects, sol, distAzCalc, 100d, 0.95, 30d);
		
		subSects = sampler.getSubsampled(subSects, 1);
		System.out.println("Downsampled to "+subSects.size()+" subsects");
		
		double[][] corrs = sampler.calcCorrs(subSects);
		int numEach = 5;
		
//		benchmarkJama(corrs, numEach);
		benchmarkApache(corrs, numEach);
//		benchmarkColt(corrs, numEach);
	}
	
	private static void benchmarkJama(double[][] corrs, int num) {
		Matrix C = new Matrix(corrs);
		
		Stopwatch totWatch = Stopwatch.createStarted();
		
		System.out.println("Jama eigen test");
		for (int i=0; i<num; i++) {
			System.out.println("Jama eig "+i);
			Stopwatch singleWatch = Stopwatch.createStarted();
			C.eig();
			singleWatch.stop();
			System.out.println("\tTook "+elapsed(singleWatch));
		}
		totWatch.stop();
		System.out.println("\tTook "+elapsed(totWatch)+" overall");
		System.out.println("\tTook "+elapsed(totWatch, num)+" each");
	}
	
	private static void benchmarkApache(double[][] corrs, int num) {
		RealMatrix C = new Array2DRowRealMatrix(corrs, false);
		
		Stopwatch totWatch = Stopwatch.createStarted();
		
		System.out.println("Apache eigen test");
		for (int i=0; i<num; i++) {
			System.out.println("Apache eig "+i);
			Stopwatch singleWatch = Stopwatch.createStarted();
			EigenDecomposition eig = new EigenDecomposition(C);
			eig.getD();
			singleWatch.stop();
			System.out.println("\tTook "+elapsed(singleWatch));
		}
		totWatch.stop();
		System.out.println("\tTook "+elapsed(totWatch)+" overall");
		System.out.println("\tTook "+elapsed(totWatch, num)+" each");
	}
	
	private static void benchmarkColt(double[][] corrs, int num) {
		DoubleMatrix2D C = new DenseDoubleMatrix2D(corrs);
		
		Stopwatch totWatch = Stopwatch.createStarted();
		
		System.out.println("Colt eigen test");
		for (int i=0; i<num; i++) {
			System.out.println("Colt eig "+i);
			Stopwatch singleWatch = Stopwatch.createStarted();
			DenseDoubleEigenvalueDecomposition eig = new DenseDoubleEigenvalueDecomposition(C);
			eig.getD();
			singleWatch.stop();
			System.out.println("\tTook "+elapsed(singleWatch));
		}
		totWatch.stop();
		System.out.println("\tTook "+elapsed(totWatch)+" overall");
		System.out.println("\tTook "+elapsed(totWatch, num)+" each");
	}
	
	private static DecimalFormat df = new DecimalFormat("0.00");
	
	private static String elapsed(Stopwatch watch) {
		return elapsed(watch, 1);
	}
	
	private static String elapsed(Stopwatch watch, int instances) {
		double millis = watch.elapsed(TimeUnit.MILLISECONDS);
		if (instances > 1)
			millis /= (double)instances;
		double secs = millis/1000d;
		if (secs < 90d)
			return df.format(secs)+" s";
		double mins = secs/60d;
		return df.format(mins)+" m";
	}

}
