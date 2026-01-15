package scratch.kevin.tdProbModelPlayground;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.TimeDepFaultSystemSolutionERF;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class TD_CompareProbsToOrig {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = TD_ERF_Example.fetchU3_BA();
		FaultSystemSolutionERF origERF = new FaultSystemSolutionERF(sol);
		TimeDepFaultSystemSolutionERF newERF = new TimeDepFaultSystemSolutionERF(sol);
		
		TD_Benchmark.parameterizeERF(origERF, true);
		TD_Benchmark.parameterizeERF(newERF, true);
		
		double[] origPoisson = getProbs(origERF);
		double[] newPoisson = getProbs(newERF);
		
		TD_Benchmark.parameterizeERF(origERF, false);
		TD_Benchmark.parameterizeERF(newERF, false);
		
		double[] origTD = getProbs(origERF);
		double[] newTD = getProbs(newERF);
		
		report("Poisson: original vs new", origPoisson, newPoisson);
		
		report("Time-dependent: original vs new", origTD, newTD);
		
		// Optional: compare TD gain relative to each impl's own Poisson (often the most diagnostic)
		double[] origGain = gains(origPoisson, origTD);
		double[] newGain = gains(newPoisson, newTD);
		report("TD gain (TD/Poisson): original vs new", origGain, newGain);
	}
	
	private static double[] getProbs(BaseFaultSystemSolutionERF erf) {
		erf.updateForecast();
		double[] ret = new double[erf.getSolution().getRupSet().getNumRuptures()];
		for (int i=0; i<erf.getNumFaultSystemSources(); i++)
			ret[erf.getFltSysRupIndexForSource(i)] = erf.getSource(i).computeTotalProb();
		return ret;
	}
	
	/**
	 * Returns elementwise TD/Poisson gains with safe handling for Poisson==0.
	 * - If poisson==0 and td==0 => gain=1 (no information, treat as neutral)
	 * - If poisson==0 and td>0 => gain=+Inf (flag)
	 * - If poisson>0 => td/poisson
	 */
	private static double[] gains(double[] poisson, double[] td) {
		int n = poisson.length;
		double[] ret = new double[n];
		for (int i = 0; i < n; i++) {
			double p = poisson[i];
			double t = td[i];
			if (p > 0d) {
				ret[i] = t / p;
			} else if (t == 0d) {
				ret[i] = 1d;
			} else {
				ret[i] = Double.POSITIVE_INFINITY;
			}
		}
		return ret;
	}

	private static void report(String label, double[] a, double[] b) {
		if (a.length != b.length)
			throw new IllegalStateException("Length mismatch: " + a.length + " vs " + b.length);

		int n = a.length;

		// Basic counts
		int bothZero = 0;
		int aZeroOnly = 0;
		int bZeroOnly = 0;
		int bothFinitePos = 0;
		int eitherNaN = 0;
		int eitherInf = 0;

		// Error metrics (for pairs where both are finite)
		ArrayList<Double> absErr = new ArrayList<>();
		ArrayList<Double> relErr = new ArrayList<>();
		ArrayList<Double> log10RatioAbs = new ArrayList<>();

		// Track worst cases
		int maxAbsIdx = -1;
		double maxAbs = -1d;

		int maxRelIdx = -1;
		double maxRel = -1d;

		int maxLogIdx = -1;
		double maxLog = -1d;

		double sumA = 0d;
		double sumB = 0d;

		for (int i = 0; i < n; i++) {
			double x = a[i];
			double y = b[i];

			if (Double.isNaN(x) || Double.isNaN(y)) {
				eitherNaN++;
				continue;
			}
			if (Double.isInfinite(x) || Double.isInfinite(y)) {
				eitherInf++;
				continue;
			}

			sumA += x;
			sumB += y;

			boolean xZero = x == 0d;
			boolean yZero = y == 0d;

			if (xZero && yZero) {
				bothZero++;
				continue;
			} else if (xZero) {
				aZeroOnly++;
				// abs error still meaningful
				double ae = Math.abs(y);
				absErr.add(ae);
				if (ae > maxAbs) {
					maxAbs = ae;
					maxAbsIdx = i;
				}
				continue;
			} else if (yZero) {
				bZeroOnly++;
				double ae = Math.abs(x);
				absErr.add(ae);
				if (ae > maxAbs) {
					maxAbs = ae;
					maxAbsIdx = i;
				}
				continue;
			}

			// Both positive/nonzero
			bothFinitePos++;

			double ae = Math.abs(x - y);
			absErr.add(ae);

			// Symmetric relative error: |x-y| / ((|x|+|y|)/2)
			double denom = (Math.abs(x) + Math.abs(y)) * 0.5;
			double re = denom > 0d ? (ae / denom) : Double.NaN;
			relErr.add(re);

			// Log-ratio distance: |log10(y/x)|
			double lr = Math.abs(Math.log10(y / x));
			log10RatioAbs.add(lr);

			if (ae > maxAbs) {
				maxAbs = ae;
				maxAbsIdx = i;
			}
			if (!Double.isNaN(re) && re > maxRel) {
				maxRel = re;
				maxRelIdx = i;
			}
			if (lr > maxLog) {
				maxLog = lr;
				maxLogIdx = i;
			}
		}

		double[] absArr = toArray(absErr);
		double[] relArr = toArray(relErr);
		double[] logArr = toArray(log10RatioAbs);

		System.out.println();
		System.out.println("============================================================");
		System.out.println(label);
		System.out.println("N ruptures: " + n);
		System.out.println("NaN pairs skipped: " + eitherNaN + ", Inf pairs skipped: " + eitherInf);
		System.out.println("Both zero: " + bothZero
				+ ", a==0 only: " + aZeroOnly
				+ ", b==0 only: " + bZeroOnly
				+ ", both nonzero finite: " + bothFinitePos);
		System.out.println();

		System.out.println("Sum(a): " + fmt(sumA) + ", Sum(b): " + fmt(sumB)
				+ ", Sum ratio b/a: " + (sumA > 0d ? fmt(sumB / sumA) : "NaN"));

		printDist("Absolute error |a-b|", absArr, "prob");
		printDist("Symmetric relative error", relArr, "");
		printDist("|log10(b/a)|", logArr, "dex");

		if (maxAbsIdx >= 0) {
			System.out.println();
			System.out.println("Worst absolute error:");
			System.out.println("\tidx=" + maxAbsIdx + " a=" + sci(a[maxAbsIdx]) + " b=" + sci(b[maxAbsIdx])
					+ " |a-b|=" + sci(Math.abs(a[maxAbsIdx] - b[maxAbsIdx])));
		}
		if (maxRelIdx >= 0) {
			System.out.println("Worst symmetric relative error:");
			System.out.println("\tidx=" + maxRelIdx + " a=" + sci(a[maxRelIdx]) + " b=" + sci(b[maxRelIdx])
					+ " rel=" + fmt(symRel(a[maxRelIdx], b[maxRelIdx])));
		}
		if (maxLogIdx >= 0) {
			System.out.println("Worst |log10(b/a)|:");
			System.out.println("\tidx=" + maxLogIdx + " a=" + sci(a[maxLogIdx]) + " b=" + sci(b[maxLogIdx])
					+ " |log10(b/a)|=" + fmt(Math.abs(Math.log10(b[maxLogIdx] / a[maxLogIdx]))));
		}
	}

	private static void printDist(String name, double[] vals, String units) {
		if (vals.length == 0) {
			System.out.println(name + ": (no comparable values)");
			return;
		}
		Arrays.sort(vals);

		double min = vals[0];
		double p50 = DataUtils.median(vals);
		double p90 = percentileSorted(vals, 0.90);
		double p95 = percentileSorted(vals, 0.95);
		double p99 = percentileSorted(vals, 0.99);
		double max = vals[vals.length - 1];

		String suffix = units == null || units.isEmpty() ? "" : " " + units;

		System.out.println(name + " (n=" + vals.length + "):");
		System.out.println("\tmin=" + sci(min) + suffix
				+ "  p50=" + sci(p50) + suffix
				+ "  p90=" + sci(p90) + suffix
				+ "  p95=" + sci(p95) + suffix
				+ "  p99=" + sci(p99) + suffix
				+ "  max=" + sci(max) + suffix);
	}

	private static double percentileSorted(double[] sorted, double p) {
		// p in [0,1]. Linear interpolation between nearest ranks.
		if (p <= 0d)
			return sorted[0];
		if (p >= 1d)
			return sorted[sorted.length - 1];

		double idx = p * (sorted.length - 1);
		int i0 = (int)Math.floor(idx);
		int i1 = (int)Math.ceil(idx);
		if (i0 == i1)
			return sorted[i0];
		double f = idx - i0;
		return sorted[i0] * (1d - f) + sorted[i1] * f;
	}

	private static double[] toArray(ArrayList<Double> list) {
		double[] ret = new double[list.size()];
		for (int i = 0; i < ret.length; i++)
			ret[i] = list.get(i);
		return ret;
	}

	private static double symRel(double a, double b) {
		double denom = (Math.abs(a) + Math.abs(b)) * 0.5;
		if (denom == 0d)
			return Double.NaN;
		return Math.abs(a - b) / denom;
	}

	private static String fmt(double val) {
		if (Double.isNaN(val))
			return "NaN";
		if (Double.isInfinite(val))
			return val > 0 ? "Inf" : "-Inf";
		return String.format("%.4f", val);
	}

	private static String sci(double val) {
		if (Double.isNaN(val))
			return "NaN";
		if (Double.isInfinite(val))
			return val > 0 ? "Inf" : "-Inf";
		return String.format("%.6e", val);
	}

}
