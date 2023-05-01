package scratch.kevin.spatialVar;

import com.google.common.base.Preconditions;

public class LothBaker2013_SpatialVarCalc {

	/*
	 *  b1 b2 b3 values from Loth&Baker 2013
	 */
	static final double[][] b1, b2, b3;
	
	static final double[] periods = { 0.01, 0.1, 0.2, 0.5, 1, 2, 5, 7.5, 10};
	
	static {
		b1 = new double[periods.length][];
		b2 = new double[periods.length][];
		b3 = new double[periods.length][];
		
		int row = 0;
		b1[row++] = new double[] { 0.30, 0.24, 0.23, 0.22, 0.16, 0.07, 0.03, 0, 0 };
		b1[row++] = new double[] { 0.24, 0.27, 0.19, 0.13, 0.08, 0, 0, 0, 0 };
		b1[row++] = new double[] { 0.23, 0.19, 0.26, 0.19, 0.12, 0.04, 0, 0, 0 };
		b1[row++] = new double[] { 0.22, 0.13, 0.19, 0.32, 0.23, 0.14, 0.09, 0.06, 0.04 };
		b1[row++] = new double[] { 0.16, 0.08, 0.12, 0.23, 0.32, 0.22, 0.13, 0.09, 0.07 };
		b1[row++] = new double[] { 0.07, 0, 0.04, 0.14, 0.22, 0.33, 0.23, 0.19, 0.16 };
		b1[row++] = new double[] { 0.03, 0, 0, 0.09, 0.13, 0.23, 0.34, 0.29, 0.24 };
		b1[row++] = new double[] { 0, 0, 0, 0.06, 0.09, 0.19, 0.29, 0.30, 0.25 };
		b1[row++] = new double[] { 0, 0, 0, 0.04, 0.07, 0.16, 0.24, 0.25, 0.24 };
		
		row = 0;
		b2[row++] = new double[] { 0.31, 0.26, 0.27, 0.24, 0.17, 0.11, 0.08, 0.06, 0.05 };
		b2[row++] = new double[] { 0.26, 0.29, 0.22, 0.15, 0.07, 0, 0, 0, -0.03 };
		b2[row++] = new double[] { 0.27, 0.22, 0.29, 0.24, 0.15, 0.09, 0.03, 0.02, 0 };
		b2[row++] = new double[] { 0.24, 0.15, 0.24, 0.33, 0.27, 0.23, 0.17, 0.14, 0.14 };
		b2[row++] = new double[] { 0.17, 0.07, 0.15, 0.27, 0.38, 0.34, 0.23, 0.19, 0.21 };
		b2[row++] = new double[] { 0.11, 0, 0.09, 0.23, 0.34, 0.44, 0.33, 0.29, 0.32 };
		b2[row++] = new double[] { 0.08, 0, 0.03, 0.17, 0.23, 0.33, 0.45, 0.42, 0.42 };
		b2[row++] = new double[] { 0.06, 0, 0.02, 0.14, 0.19, 0.29, 0.42, 0.47, 0.47 };
		b2[row++] = new double[] { 0.05, -0.03, 0, 0.14, 0.21, 0.32, 0.42, 0.47, 0.54 };
		
		row = 0;
		b3[row++] = new double[] { 0.38, 0.36, 0.35, 0.17, 0.04, 0.04, 0, 0.03, 0.08 };
		b3[row++] = new double[] { 0.36, 0.43, 0.35, 0.13, 0, 0.02, 0, 0.02, 0.08 };
		b3[row++] = new double[] { 0.35, 0.35, 0.45, 0.11, -0.04, -0.02, -0.04, -0.02, 0.03 };
		b3[row++] = new double[] { 0.17, 0.13, 0.11, 0.35, 0.20, 0.06, 0.02, 0.04, 0.02 };
		b3[row++] = new double[] { 0.04, 0, -0.04, 0.20, 0.30, 0.14, 0.09, 0.12, 0.04 };
		b3[row++] = new double[] { 0.04, 0.02, -0.02, 0.06, 0.14, 0.22, 0.12, 0.13, 0.09 };
		b3[row++] = new double[] { 0, 0, -0.04, 0.02, 0.09, 0.12, 0.21, 0.17, 0.13 };
		b3[row++] = new double[] { 0.03, 0.02, -0.02, 0.04, 0.12, 0.13, 0.17, 0.23, 0.10 };
		b3[row++] = new double[] { 0.08, 0.08, 0.03, 0.02, 0.04, 0.09, 0.13, 0.10, 0.22 };
	}
	
	/**
	 * Calculates covariance between a log intensity measure at two sides that are the given distance apart, and between
	 * intensity measures of the given periods. The correlation coefficients will be bilinearly interpolated if periods
	 * are supplied between those defined in Loth & Baker (2013)
	 * 
	 * @param distance distance between 2 sites, in KM
	 * @param period1 period of site 1
	 * @param period2 period of site 2
	 * @return
	 */
	public static double calcCovariance(double distance, double period1, double period2) {
		validatePeriod(period1);
		validatePeriod(period2);
		int ind1 = getPeriodIndexBefore(period1);
		int ind2 = getPeriodIndexBefore(period2);
		
		double B1 = getInterpB(b1, period1, period2, ind1, ind2);
		double B2 = getInterpB(b2, period1, period2, ind1, ind2);
		
		double neg3H = -3d*distance;
		
		double ret = B1 * Math.exp(neg3H/20d) + B2*Math.exp(neg3H/70d);
		if (distance == 0d) {
			double B3 = getInterpB(b3, period1, period2, ind1, ind2);
			ret += B3;
		}
		return ret;
	}
	
	static void validatePeriod(double period) {
		Preconditions.checkState(period >= 0d, "Period must be >=%s: %s", periods[0], period);
		Preconditions.checkState(period <= periods[periods.length-1], "Period must be <=%s: %s",
				periods[periods.length-1], period);
	}
	
	static double getInterpB(double[][] b, double period1, double period2, int ind1, int ind2) {
		// X here corresponds to period/ind1, Y to period/ind2
		int x0, x1, y0, y1;
		if (ind1 < 0) {
			// must be exactly on the bottom edge (checked externally)
			Preconditions.checkState(ind1 == -1);
			Preconditions.checkState((float)period1 <= (float)periods[0]);
			x0 = 0;
			x1 = 0;
			if (period1 < periods[0])
				System.err.println("Warning: can't interpolate for period="
						+(float)period1+", using lower bound of "+(float)periods[0]);
		} else {
			Preconditions.checkState(ind1 <= periods.length);
			x0 = ind1;
			x1 = ind1+1;
			if (period1 > periods[periods.length-1]) {
				System.err.println("Warning: can't interpolate for period="
						+(float)period1+", using upper bound of "+(float)periods[periods.length-1]);
			} else {
				Preconditions.checkState((float)period1 >= (float)periods[x0],
						"period1=%s, periods[%s]=%s, periods[%s]=%s", period2, x0, periods[x0], x1, periods[x1]);
				Preconditions.checkState((float)period1 <= (float)periods[x1],
						"period1=%s, periods[%s]=%s, periods[%s]=%s", period2, x0, periods[x0], x1, periods[x1]);
			}
			if (x1 == periods.length)
				x1 = x0;
		}
		if (ind2 < 0) {
			// must be exactly on the bottom edge (checked externally)
			Preconditions.checkState(ind2 == -1);
			Preconditions.checkState((float)period2 <= (float)periods[0]);
			y0 = 0;
			y1 = 0;
			if (period2 < periods[0])
				System.err.println("Warning: can't interpolate for period="
						+(float)period2+", using lower bound of "+(float)periods[0]);
		} else {
			Preconditions.checkState(ind2 <= periods.length);
			y0 = ind2;
			y1 = ind2+1;
			if (period2 > periods[periods.length-1]) {
				System.err.println("Warning: can't interpolate for period="
						+(float)period2+", using upper bound of "+(float)periods[periods.length-1]);
			} else {
				Preconditions.checkState((float)period2 >= (float)periods[y0],
						"period2=%s, periods[%s]=%s, periods[%s]=%s", period2, y0, periods[y0], y1, periods[y1]);
				Preconditions.checkState((float)period2 <= (float)periods[y1],
						"period2=%s, periods[%s]=%s, periods[%s]=%s", period2, y0, periods[y0], y1, periods[y1]);
			}
			if (y1 == periods.length)
				y1 = y0;
		}
		
		// "central"
		double s00 = b[x0][y0];
		// to the right
		double s01 = b[x1][y0];
		// below
		double s10 = b[x0][y1];
		// below and to the right
		double s11 = b[x1][y1];
		
		double xfrac = x0 == x1 ? 0d : (period1 - periods[x0])/(periods[x1] - periods[x0]);
		double yfrac = y0 == y1 ? 0d : (period2 - periods[y0])/(periods[y1] - periods[y0]);
		
		return (1 - yfrac) * ((1 - xfrac)*s00 + xfrac*s01) + 
			    yfrac * ((1 - xfrac)*s10 + xfrac*s11);
	}
	
	static int getPeriodIndexBefore(double period) {
		if (period <= periods[0])
			return -1;
		int ind = -1;
//		System.out.println("period="+period);
		for (double p : periods) {
//			System.out.println("checking p="+p+" with ind="+ind+" and periods[ind]="+periods[ind]);
			if (p >= period)
				break;
			ind++;
		}
		Preconditions.checkState(ind < periods.length);
		Preconditions.checkState((float)period >= (float)periods[ind],
				"getPeriodIndexBefore failed with period=%s, ind=%s, periods[%s]=%s", period, ind, ind, periods[ind]);
		if (ind < periods.length-1)
			Preconditions.checkState((float)period <= (float)periods[ind+1],
				"getPeriodIndexBefore failed with period=%s, ind=%s, periods[%s]=%s", period, ind, ind, periods[ind]);
		return ind;
	}

	public static void main(String[] args) {
		System.out.println(calcCovariance(10d, 1d, 0d));
		System.out.println(calcCovariance(10d, 1d, 0.5d));
		System.out.println(calcCovariance(10d, 1d, 0.6d));
		System.out.println(calcCovariance(10d, 1d, 0.7d));
		System.out.println(calcCovariance(10d, 1d, 0.8d));
		System.out.println(calcCovariance(10d, 1d, 0.9d));
		System.out.println(calcCovariance(10d, 1d, 1d));
		System.out.println(calcCovariance(10d, 1d, 2d));
		System.out.println(calcCovariance(10d, 1d, 3d));
//		System.out.println(calcCovariance(10d, 1d, 11d));
	}

}
