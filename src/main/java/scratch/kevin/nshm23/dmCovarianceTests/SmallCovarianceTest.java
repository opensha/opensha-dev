package scratch.kevin.nshm23.dmCovarianceTests;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
//import org.opensha.commons.calc.cholesky.CholeskyDecomposition;
//import org.opensha.commons.calc.cholesky.NearPD;
import org.opensha.commons.calc.cholesky.NearPD;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;

import com.google.common.base.Preconditions;

import Jama.Matrix;
import net.mahdilamb.colormap.Colors;

public class SmallCovarianceTest {

	private static double getCorr(int i, int j) {
		if (i == j)
			return 1d;
		return -0.5;
	}
	
	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException {
		RandomGenerator rng = new Well19937c(System.currentTimeMillis());
		int N = 4;
		
		double[] means = new double[N];
		double[] sigmas = new double[N];
		for (int i=0; i<N; i++) {
			means[i] = 10;
			sigmas[i] = 10;
		}
		
//		Double conditionalFirstSlip = 2d;
		Double conditionalFirstSlip = 20d;
//		Double conditionalFirstSlip = 30d;
//		Double conditionalFirstSlip = null;
		
//		int numSamples = 100;
//		int samplePrintMod = 1;
		int numSamples = 10000;
		int samplePrintMod = 1000;
		
		double[][] corrs = new double[N][N];
		for (int i=0; i<N; i++)
			for (int j=0; j<N; j++)
				corrs[i][j] = getCorr(i, j);
		
		Preconditions.checkState(corrs.length == corrs[0].length);
		RealMatrix C = MatrixUtils.createRealMatrix(corrs);
		
		System.out.println("COV matrix:");
		printMatrix(C);
		
		System.out.println("Doing CholeskyDecomposition");
		CholeskyDecomposition chol;
		try {
			chol = new CholeskyDecomposition(C);
		} catch (MathIllegalArgumentException ex) {
			System.out.println("Matrix is not positive definite, will try to iteratively find nearPD solution");
			NearPD nearPD = new NearPD();
			nearPD.setKeepDiag(true);
			nearPD.setUseApache(true);
//			nearPd.setEigTol(1.e-6);
			boolean success = nearPD.calcNearPD(new Matrix(C.getData()), false);
			double normFrob = nearPD.getFrobNorm();
			if (!success) {
				System.out.println("B is size "+C.getRowDimension()+"x"+C.getColumnDimension()+". normFrob="+normFrob);
//				printMatrix(B);
//				throw new RuntimeException("Error: nearPD failed to converge, the correlation matrix maybe" +
//						" significantly different from a PD matrix, check that the correlation equations" +
//						" used are reasonable");
				System.err.println("WARNING: nearPD failed to converge, the correlation matrix maybe" +
						" significantly different from a PD matrix, check that the correlation equations" +
						" used are reasonable. Convergence: "+(float)nearPD.getConvergedTolerence());
			}
			
			
			C = MatrixUtils.createRealMatrix(nearPD.getX().getArray());
			System.out.println("Nearest PD:");
			printMatrix(C);
			// Now get the CholDecomp of this nearest matrix
			chol = new CholeskyDecomposition(C);
//			Preconditions.checkState(cholPD.isSPD(), "Error: Even after NearPD the matrix is not PD");
//			chol = cholPD;
		}
		
		List<double[]> samples = new ArrayList<>(numSamples);
		
		if (conditionalFirstSlip == null) {
			RealMatrix L = chol.getL();
			
			System.out.println("L");
			printMatrix(L);
			
			for (int s=0; s<numSamples; s++) {
				boolean print = s % samplePrintMod == 0;
				if (print) System.out.println("Sample "+s);
				// random gaussian
				double[] rands = new double[N];
				for (int i=0; i<N; i++)
					rands[i] = rng.nextGaussian();
				RealMatrix X = MatrixUtils.createColumnRealMatrix(rands);
				
				// sample z-scores
				RealMatrix zScores = L.multiply(X);
//				Matrix zScores = L.times(X);
				if (print) {
					System.out.println("Sampled z-scores:");
					printMatrix(zScores);
				}
				
				// now convert back to real values
				double[] sample = new double[N];
				for (int i=0; i<N; i++)
					sample[i] = slip(means[i], sigmas[i], zScores.getEntry(i, 0));
				
				samples.add(sample);
				
				if (print) {
					System.out.println("Sampled slip rates:");
					printArray(sample);
					System.out.println("Sum:\t"+(float)StatUtils.sum(sample));
					System.out.println();
				}
	
			}
		} else {
			ConditionalNormalSampler cond = new ConditionalNormalSampler(C, rng);
			
			// convert that physical value to the latent z-score
			double a = (conditionalFirstSlip - means[0]) / sigmas[0];   // z = (x-μ)/σ
			
			System.out.println("Mean z-scores:");
			double[] meanZ = cond.getMeanZ(a);
			printArray(meanZ);
			System.out.println("Mean slips:");
			double[] meanSlips = new double[N];
			for (int i=0; i<N; i++)
				meanSlips[i] = slip(means[i], sigmas[i], meanZ[i]);
			printArray(meanSlips);

			// ------------------------------------------------------------------
			// Monte-Carlo loop
			// ------------------------------------------------------------------
			for (int s = 0; s < numSamples; s++) {
				boolean print = s % samplePrintMod == 0;
				if (print) System.out.println("Sample "+s);

			    // 1) latent draw with z[0] = a
			    double[] z = cond.sample(a);      // z[0] == a by construction
			    
				if (print) {
					System.out.println("Sampled z-scores:");
					printArray(z);
				}

			    // 2) map latent normals -> physical normals
			    double[] sample = new double[N];
			    for (int i = 0; i < N; i++)
			        sample[i] = slip(means[i], sigmas[i], z[i]);   // x = μ + σ z

			    // (optional) make absolutely sure:
			    sample[0] = conditionalFirstSlip;
			    
			    samples.add(sample);
				
				if (print) {
					System.out.println("Sampled slip rates:");
					printArray(sample);
					System.out.println("Sum:\t"+(float)StatUtils.sum(sample));
					System.out.println();
				}
			}

		}
		
		double maxSum = 0d;
		for (double[] sample : samples)
			maxSum = Math.max(maxSum, StatUtils.sum(sample));
		// now plot
		HistogramFunction sumHist = HistogramFunction.getEncompassingHistogram(0d, maxSum+1, 1d);
		HistogramFunction[] indvHists = new HistogramFunction[N];
		for (int i=0; i<N; i++)
			indvHists[i] = new HistogramFunction(sumHist.getMinX(), sumHist.size(), sumHist.getDelta());
		for (double[] sample : samples) {
			double sum = StatUtils.sum(sample);
			sumHist.add(sumHist.getClosestXIndex(sum), 1d);
			for (int i=0; i<N; i++)
				indvHists[i].add(indvHists[i].getClosestXIndex(sample[i]), 1d);
		}
		
		List<PlotSpec> plots = new ArrayList<>();
		List<Range> xRanges = List.of(new Range(0d, maxSum+2));
		List<Range> yRanges = new ArrayList<>();
		
		for (HistogramFunction hist : indvHists) {
			plots.add(buildHistPlot(hist, Colors.tab_blue));
			yRanges.add(new Range(0d, hist.getMaxY()+2));
		}
		
		plots.add(buildHistPlot(sumHist, Colors.tab_orange));
		yRanges.add(new Range(0d, sumHist.getMaxY()+2));
		
		GraphWindow gw = new GraphWindow(new GraphWidget(plots, PlotUtils.getDefaultFigurePrefs(), false, false, xRanges, yRanges));
		
		gw.setVisible(true);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}
	
	private static double slip(double mean, double sigma, double z) {
		return Math.max(0, mean + sigma * z);
	}
	
	private static PlotSpec buildHistPlot(HistogramFunction hist, Color color) {
		return new PlotSpec(List.of(hist), List.of(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color)),
				" ", "Slip Rate (mm/yr)", "Count");
	}
	
	static void printMatrix(RealMatrix matrix) {
		for (int row=0; row<matrix.getRowDimension(); row++) {
			for (int col=0; col<matrix.getColumnDimension(); col++)
				System.out.print((float)matrix.getEntry(row, col)+"\t");
			System.out.println();
		}
	}

	
	static void printArray(double[] array) {
		for (double val : array)
			System.out.print((float)val+"\t");
		System.out.println();
	}
	
	/**
	 * Samples an (n-1)-vector z₂ ~ N( μ_cond , Σ_cond )
	 * given that the first latent-normal component z₁ is fixed to 'a'.
	 *
	 *  z = [ z₁ ; z₂ ]  ~  N(0, R)
	 *  z₂ | z₁=a  ~  N( R21 * a ,  R22 − R21 R21ᵀ )
	 *
	 *  All conditioning is done in latent space; you map every component
	 *  through its inverse-CDF afterwards.
	 */
	public static class ConditionalNormalSampler {

		private final RealVector r21;     // (n-1) × 1   ==  R21
		private final RealMatrix L;       // lower-triangular chol( Σ_cond )
		private final GaussianRandomGenerator gauss;

		/**
		 * @param R   full n×n correlation (or covariance) matrix, assumed SPD
		 * @param rng your preferred Commons-Math RandomGenerator
		 */
		public ConditionalNormalSampler(RealMatrix R, RandomGenerator rng) {
			int n = R.getRowDimension();
			// blocks: R = [ 1  r12ᵀ ; r12  R22 ]
			r21 = R.getSubMatrix(1, n - 1, 0, 0).getColumnVector(0);   // (n-1)×1
			RealMatrix R22 = R.getSubMatrix(1, n - 1, 1, n - 1);

			// Σ_cond = R22 − r21 r21ᵀ   (Schur complement)
			RealMatrix schur = R22.subtract(r21.outerProduct(r21));
			RealMatrix Lmat  = new CholeskyDecomposition(schur, 1e-12, 1e-9).getL();

			this.L = Lmat;
			this.gauss = new GaussianRandomGenerator(rng);
		}

		/** Draw the full latent vector z given z₁ = a. */
		public double[] sample(double a) {
			int m = r21.getDimension();          // n-1
			double[] zRand = new double[m];
			for (int i = 0; i < m; i++) zRand[i] = gauss.nextNormalizedDouble();

			RealVector z2 = L.operate(new ArrayRealVector(zRand))
					.add(r21.mapMultiply(a));

			double[] out = new double[m + 1];
			out[0] = a;
			for (int i = 0; i < m; i++) out[i + 1] = z2.getEntry(i);
			return out;
		}
		
		public double[] getMeanZ(double a) {
			int m = r21.getDimension();
			// r21 is exactly the vector of ρ_{j1} in ConditionalNormalSampler
			RealVector condMeanZ = r21.mapMultiply(a);          // ⟨ρ_{j1}·a⟩
			double[] ret = new double[m+1];
			ret[0] = a;
			for (int i=0; i<m; i++)
				ret[i+1] = condMeanZ.getEntry(i);
			return ret;
		}
	}

}
