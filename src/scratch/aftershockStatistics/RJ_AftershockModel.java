package scratch.aftershockStatistics;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.distribution.PoissonDistribution;
//import org.mongodb.morphia.annotations.Transient;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;


/**
 * This represents a Reasenberg-Jones (1989, 1994) aftershock model.
 * 
 * Subclasses should normalize Likelihood values so they sum to 1.0 over the range of parameter-values specified;
 * this can be accomplished using the convertArrayToLikelihood(double maxLogLikeVal) method.
 *
 * @author field
 *
 * Modified by Michael Barall.
 *
 * This model includes both aleatory uncertainty (uncertainty due to randomness in nature)
 * and epistemic uncertainty (uncertainty due to our lack of knowledge of the correct
 * parameter values).
 *
 * The aleatory uncertainty is represented by a Reasenberg-Jones probability distribution.
 * According to R&J, the rate of aftershocks of magnitude >= magMin is
 *  lambda(t) = k * (t + c)^(-p)
 * where
 *  k = 10^(a + b*(magMain - magMin))
 * Here time t is measured in days, and lambda(t) is events per day.  The parameters are:
 *  a = Reasenberg-Jones productivity
 *  b = Gutenberg-Richter exponent
 *  c = Omori offset
 *  p = Omori exponent
 *
 * In a time interval [tMin, tMax], the expected number of aftershocks is the integral
 * of lambda(t) from tMin to tMax:
 *  lambda(tMin, tMax) = INTEGRAL(lambda(t) * dt, tMin <= t <= tMax)
 * The actual number of aftershocks in the time interval is assumed to be a Poisson random
 * variable with expected value equal to lambda(tMin, tMax).  So the probability of there
 * being n aftershocks in the time interval [tMin, tMax] is
 *  PA(n) = exp(-lambda(tMin, tMax)) * (lambda(tMin, tMax))^n / n!
 * The notation PA(n) indicates this is probability due to aleatory uncertainty.
 *
 * Epistemic uncertainty is introduced by treating the parameter triple (a,p,c) as a
 * random variable, with probability density function rho(a,p,c).  The parameter b is
 * considered to be known, and so there is no probability distribution associated with b.
 * Then, the probability of there being n aftershocks in the time interval [tMin, tMax],
 * taking into account both aleatory and epistemic uncertainty, is
 *  PAE(n) = INTEGRAL(PA(n) * rho(a,p,c) * da * dp * dc)
 * where the integration runs over all allowed values of a, p, and c.
 * Notice that this is a non-Poissonian distribution.
 *
 * In this model, the probability density rho(a,p,c) is taken to be a discrete function,
 * which is defined on a uniformly-spaced grid of points in (a,p,c) space.  Each of the
 * three parameters a, p, and c is allowed to assume a set of evenly-spaced discrete values,
 * ranging from min_a to max_a, min_p to max_p, and min_c to max_c respectively.  The
 * probabilities are stored in a three-dimensional array, and the array elements are
 * normalized so they sum to 1.  Then the integral is replaced with a sum:
 *  PAE(n) = SUM(PA(n) * rho(a,p,c), min_a <= a <= max_a, min_p <= p <= max_p, min_c <= c <= max_c)
 *
 * This is an abstract class.  The main job of a subclass is to supply the probability
 * distribution rho(a,p,c).
 */
public abstract class RJ_AftershockModel {

	//----- Fields that must be set by the subclass -----

	// Debug flag.
	// When set to true, some values are written to System.out, which are logged in the OAF server.

	//@Transient
	protected boolean D = true;	// debug flag

	// The Gutenberg-Richter b-value
	// This is set to a fixed value by the subclass.
	
	protected double b = 0.0; 

	// The magnitude of the mainshock.
	// This is set to a fixed value by the subclass.
	
	protected double magMain = 0.0; 

	// The parameter space for the Reasenberg-Jones productivity parameter a and the Omori parameters p and c.
	// This class considers a discrete set of a-values, ranging from min_a to max_a.
	// There are num_a different values, evenly spaced with spacing delta_a, hence delta_a = (max_a - min_a)/(num_a - 1).
	// (If num_a == 1 then max_a == min_a and delta_a == 0.)
	// The p-values and c-values are discretized in the same way.

	protected double min_a = 0.0;
	protected double max_a = 0.0;
	protected double delta_a = 0.0;
	protected int num_a = 1;

	protected double min_p = 0.0;
	protected double max_p = 0.0;
	protected double delta_p = 0.0;
	protected int num_p = 1;

	protected double min_c = 0.0;
	protected double max_c = 0.0;
	protected double delta_c = 0.0;
	protected int num_c = 1;

	// Likelihood values for each parameter triple (a, p, c).
	// The dimensions are apc_likelihood[num_a][num_p][num_c].
	// This array may contain either likelihood or log-likelihood depending on the context.
	// Note: A subclass or user of this class must supply the likelihood values.
	// They are not computed in this class.

	//@Transient
	protected double[][][] apc_likelihood = null;

	// The fraction of the (a,p,c) probability distribution that can be ignored as negligably small.

	protected double apc_tail_fraction = 0.0001;

	//----- Fields that are set and used by this class -----

	// Index values for the maximum likelihood parameter triple (a, p, c).
	// In other words, apc_likelihood[max_a_index][max_p_index][max_c_index] is the largest element in apc_likelihood.
	// Note: This identifies the maximum likelihood values, considering the parameters to form a triple (a, p, c).
	// It does NOT identify the maximum likelihood values of each parameter considered separately.

	protected int max_a_index = -1;
	protected int max_p_index = -1;
	protected int max_c_index = -1;

	// Total size of the apc_likelihood matrix.
	// This is num_a * num_p * num_c.

	protected int apc_total_size = -1;

	// Total number of apc_likelihood elements in the support.
	// The support consists of elements not in the tail, and so not negligable.

	protected int apc_support_size = -1;

	// Total of all the apc_likelihood elements in the support.
	// This will be very close to 1.0.

	protected double apc_support_total = 1.0;

	// The maximum element in the tail of the apc_likelihood distribution.
	// Elements are in the tail, and considered negligable, if <= apc_max_tail_element.

	protected double apc_max_tail_element = 0.0;

	// The range of index values that contain the support of the (a,p,c) probability distribution.
	// (The support is the portion that is outside the tail, that is, non-negligable.)

	protected int a_support_lo = -1;
	protected int a_support_hi = -1;

	protected int p_support_lo = -1;
	protected int p_support_hi = -1;

	protected int c_support_lo = -1;
	protected int c_support_hi = -1;

	// The calculated mean, standard deviation, and maximum likelihood value of each parameter.

	protected double stat_a_mean = 0.0;
	protected double stat_a_sdev = 0.0;
	protected double stat_a_like = 0.0;

	protected double stat_p_mean = 0.0;
	protected double stat_p_sdev = 0.0;
	protected double stat_p_like = 0.0;

	protected double stat_c_mean = 0.0;
	protected double stat_c_sdev = 0.0;
	protected double stat_c_like = 0.0;

	// This is a discrete function that represents likelihood as a function of the number of aftershocks
	// of magnitude >= 5 during the time interval tMinDaysCurrent <= t <= tMaxDaysCurrent.
	// Specifically, this function is a collection of points (x,y) where:
	//  x = Expected number of aftershocks, as computed by the R&J formula.
	//  y = Likelihood, obtained from apc_likelihood[aIndex][pIndex][cIndex].
	// The set of points is obtained by iterating aIndex, pIndex, and cIndex over their ranges.
	// If two points have exactly the same x-value, they are combined by adding their y-values.
	// (The combination is done in class EmpiricalPoint2DToleranceSortedList.)
	//
	// The use of magnitude 5 is arbitrary.  According to the R&J formula, results for any other
	// magnitude threshold M can by obtained by multiplying all x values by 10^(b*(5 - M)).
	//
	// Conceptually, this function captures the epistemic uncertainty associated with the
	// uncertainty in the parameters (a,p,c).  The triple (a,p,c) is considered to be a random
	// variable, whose probability density is given by apc_likelihood.  The R&J expected number
	// of aftershocks is then a function of (a,p,c) and so has an induced probability
	// distribution.  This function is the induced probability density.
	//
	// Note that the x-values of this function are themselves merely expected numbers of
	// aftershocks.  Even if "correct" values of (a,p,c) were known, there would still be
	// aleatory uncertainty due to randomness in nature, which causes the observed number
	// of aftershocks to differ from the expected number.
	//
	// This function must be used with care, because the x-values are neither equally-spaced
	// nor binned.  You can use numMag5_DistributionFunc.getMean(),
	// numMag5_DistributionFunc.getStdDev(), and numMag5_DistributionFunc.getInterpolatedFractile(fractile)
	// to get the mean, standard deviation, and fractile of the R&J expected number of
	// aftershocks.  But, for example, a naive computation of the mode would be wrong.

	//@Transient
	protected ArbDiscrEmpiricalDistFunc numMag5_DistributionFunc = null;

	// The time interval used to calculate numMag5_DistributionFunc,
	// measured in days after the mainshock.
	// Note: The value of numMag5_DistributionFunc is cached so it need not be recomputed
	// every time it is needed.  The values of tMinDaysCurrent and tMaxDaysCurrent are
	// used to check whether recomputation is needed.

	//@Transient
	protected double tMinDaysCurrent = -1.0;
	//@Transient
	protected double tMaxDaysCurrent = -1.0;




	// Return the maximum-likelihood value of a, p, or c.
	// Note: These are the values of a, p, and c in the maximum-likelihood triple (a, p, c).
	// They are NOT the maximum likelihood values of a, p, and c considered separately.
	// Note: These are the a, p, and c values for the largest element in apc_likelihood.
	
	public double getMaxLikelihood_a() {return get_a(max_a_index);}
	
	public double getMaxLikelihood_p() {return get_p(max_p_index);}
	
	public double getMaxLikelihood_c() {return get_c(max_c_index);}

	// Return the floating-point value of a, p, or c, given its index number.
		
	protected double get_a(int aIndex) {return min_a+aIndex*delta_a;}
	
	protected double get_p(int pIndex) {return min_p+pIndex*delta_p;}
	
	protected double get_c(int cIndex) {return min_c+cIndex*delta_c;}

	// Return the magnitude of the mainshock.
	
	public double getMainShockMag() {return magMain;}

	// Return the Gutenberg-Richter b value.
	
	public double get_b() {return b;}

	// Return the mean value of a, p, or c.

	public double getMean_a() {return stat_a_mean;}
	public double getMean_p() {return stat_p_mean;}
	public double getMean_c() {return stat_c_mean;}

	// Return the standard deviation of a, p, or c.

	public double getStdDev_a() {return stat_a_sdev;}
	public double getStdDev_p() {return stat_p_sdev;}
	public double getStdDev_c() {return stat_c_sdev;}

	// Return the maximum likelihood value of a, p, or c (alternate implementation).

	public double getMaxLike_a() {return stat_a_like;}
	public double getMaxLike_p() {return stat_p_like;}
	public double getMaxLike_c() {return stat_c_like;}

	// Return parameter space information.

	public double getMin_a() {return min_a;}
	public double getMin_p() {return min_p;}
	public double getMin_c() {return min_c;}

	public double getMax_a() {return max_a;}
	public double getMax_p() {return max_p;}
	public double getMax_c() {return max_c;}

	public double getDelta_a() {return delta_a;}
	public double getDelta_p() {return delta_p;}
	public double getDelta_c() {return delta_c;}

	public int getNum_a() {return num_a;}
	public int getNum_p() {return num_p;}
	public int getNum_c() {return num_c;}




	/**
	 * Turn verbose mode on or off.
	 */
	public void set_verbose(boolean f_verbose) {
		D = f_verbose;
		return;
	}




	/**
	 * Set the tail fraction that is used for clipping the (a,p,c) probability distribution.
	 * It could be set to zero (or a very small value) if the (a,p,c) distribution
	 * has already been clipped external to this class.
	 * For example, the default value of 0.0001 means to throw out the smallest elements
	 * in apc_likelihood (after normalizing it) until their total reaches apc_tail_fraction.
	 */
	public void set_tail_fraction(double apc_tail_fraction) {
		this.apc_tail_fraction = apc_tail_fraction;
		return;
	}




	/**
	 * Set a fixed value of parameter a.
	 */
	protected void set_fixed_a(double a) {
		num_a = 1;
		min_a = a;
		max_a = a;
		delta_a = 0.0;
		return;
	}




	/**
	 * Set a fixed value of parameter p.
	 */
	protected void set_fixed_p(double p) {
		num_p = 1;
		min_p = p;
		max_p = p;
		delta_p = 0.0;
		return;
	}




	/**
	 * Set a fixed value of parameter c.
	 */
	protected void set_fixed_c(double c) {
		num_c = 1;
		min_c = c;
		max_c = c;
		delta_c = 0.0;
		return;
	}




	/**
	 * Return the name of this model.
	 */
	public abstract String getModelName();

	


	/**
	 * Finish setting up apc_likelihood, which contains the (a,p,c) probability distribution.
	 * @param f_log = True if apc_likelihood contains log-likelihood, false if it contains likelihood.
	 * @return
	 * Convert log-likelihood to likelihood if necessary.
	 * Identify the tail of the distribution.
	 * Normalize the likelihood values so they sum to 1.0.
	 */
	protected void apcFinish(boolean f_log) {

		// Error if matrix is so large it cannot be stored in a one-dimensional matrix.
		// (Maybe should use a lower limit than Integer.MAX_VALUE)

		if (((Integer.MAX_VALUE / num_a) / num_p) / num_c == 0) {
			throw new RuntimeException("RJ_AftershockModel: Parameter likelihood matrix is too large");
		}

		// Invalidate the event count likelihood function

		numMag5_DistributionFunc = null;
		tMinDaysCurrent = -1.0;
		tMaxDaysCurrent = -1.0;

		// Find the biggest element in the matrix

		double max_element = apc_likelihood[0][0][0];
		max_a_index = 0;
		max_p_index = 0;
		max_c_index = 0;

		for (int aIndex = 0; aIndex < num_a; aIndex++) {
			for (int pIndex = 0; pIndex < num_p; pIndex++) {
				for (int cIndex = 0; cIndex < num_c; cIndex++) {
					if (apc_likelihood[aIndex][pIndex][cIndex] > max_element) {
						max_element = apc_likelihood[aIndex][pIndex][cIndex];
						max_a_index = aIndex;
						max_p_index = pIndex;
						max_c_index = cIndex;
					}
				}
			}
		}

		// If matrix contains log likelihood, convert it to likelihood

		if (f_log) {
			for (int aIndex = 0; aIndex < num_a; aIndex++) {
				for (int pIndex = 0; pIndex < num_p; pIndex++) {
					for (int cIndex = 0; cIndex < num_c; cIndex++) {
						apc_likelihood[aIndex][pIndex][cIndex] = Math.exp(apc_likelihood[aIndex][pIndex][cIndex] - max_element);
									// subtract max_element to avoid overflows
					}
				}
			}
			max_element = 1.0;
		}

		// Check for nonzero probabilities

		if (max_element < Double.MIN_NORMAL * 1.0e16) {
			throw new RuntimeException("RJ_AftershockModel: Parameter likelihood matrix effective zero");
		}

		// Initialize for statistics computation

		stat_a_mean = 0.0;
		stat_a_sdev = 0.0;
		stat_a_like = getMaxLikelihood_a();

		stat_p_mean = 0.0;
		stat_p_sdev = 0.0;
		stat_p_like = getMaxLikelihood_p();

		stat_c_mean = 0.0;
		stat_c_sdev = 0.0;
		stat_c_like = getMaxLikelihood_c();

		// Calculate total weight of the matrix, and calculate the means

		double total_weight = 0.0;
		for (int aIndex = 0; aIndex < num_a; aIndex++) {
			double a = get_a(aIndex);
			for (int pIndex = 0; pIndex < num_p; pIndex++) {
				double p = get_p(pIndex);
				for (int cIndex = 0; cIndex < num_c; cIndex++) {
					double c = get_c(cIndex);
					double w = apc_likelihood[aIndex][pIndex][cIndex];
					total_weight += w;
					stat_a_mean += a * w;
					stat_p_mean += p * w;
					stat_c_mean += c * w;
				}
			}
		}

		stat_a_mean /= total_weight;
		stat_p_mean /= total_weight;
		stat_c_mean /= total_weight;

		// Normalize the matrix so it sums to 1.0, and dump the normalized values to a one-dimensional array,
		// and compute the standard deviations

		apc_total_size = num_a * num_p * num_c;
		double[] apc_sorted = new double[apc_total_size];
		int sIndex = 0;

		for (int aIndex = 0; aIndex < num_a; aIndex++) {
			double a = get_a(aIndex);
			for (int pIndex = 0; pIndex < num_p; pIndex++) {
				double p = get_p(pIndex);
				for (int cIndex = 0; cIndex < num_c; cIndex++) {
					double c = get_c(cIndex);
					double w = apc_likelihood[aIndex][pIndex][cIndex] / total_weight;
					apc_likelihood[aIndex][pIndex][cIndex] = w;
					apc_sorted[sIndex++] = w;
					stat_a_sdev += (a - stat_a_mean) * (a - stat_a_mean) * w;
					stat_p_sdev += (p - stat_p_mean) * (p - stat_p_mean) * w;
					stat_c_sdev += (c - stat_c_mean) * (c - stat_c_mean) * w;
				}
			}
		}

		stat_a_sdev = Math.sqrt(stat_a_sdev);
		stat_p_sdev = Math.sqrt(stat_p_sdev);
		stat_c_sdev = Math.sqrt(stat_c_sdev);

		if (num_a == 1) {
			stat_a_mean = min_a;
			stat_a_sdev = 0.0;
		}

		if (num_p == 1) {
			stat_p_mean = min_p;
			stat_p_sdev = 0.0;
		}

		if (num_c == 1) {
			stat_c_mean = min_c;
			stat_c_sdev = 0.0;
		}

		// Sort the array from low to high

		Arrays.sort (apc_sorted, 0, apc_total_size);

		// Scan the sorted array to find the largest tail element

		apc_max_tail_element = 0.0;
		double tail_weight = apc_sorted[0];		// sum of all elements prior to sIndex
		for (sIndex = 1; sIndex < apc_total_size && tail_weight <= apc_tail_fraction; ++sIndex) {

			// If greater than the prior element, then the prior element could be the last element of the tail
			// (Don't let the tail end in the middle of a run of equal elements)

			if (apc_sorted[sIndex] > apc_sorted[sIndex - 1]) {
				apc_max_tail_element = apc_sorted[sIndex - 1];
			}

			// Add current element to tail weight

			tail_weight += apc_sorted[sIndex];
		}

		// Get the support bounds and total

		apc_support_size = 0;
		apc_support_total = 0.0;

		a_support_lo = num_a;
		a_support_hi = 0;

		p_support_lo = num_p;
		p_support_hi = 0;

		c_support_lo = num_c;
		c_support_hi = 0;

		for (int aIndex = 0; aIndex < num_a; aIndex++) {
			for (int pIndex = 0; pIndex < num_p; pIndex++) {
				for (int cIndex = 0; cIndex < num_c; cIndex++) {
					if (apc_likelihood[aIndex][pIndex][cIndex] > apc_max_tail_element) {
						++apc_support_size;
						apc_support_total += apc_likelihood[aIndex][pIndex][cIndex];
						a_support_lo = Math.min(a_support_lo, aIndex);
						a_support_hi = Math.max(a_support_hi, aIndex + 1);
						p_support_lo = Math.min(p_support_lo, pIndex);
						p_support_hi = Math.max(p_support_hi, pIndex + 1);
						c_support_lo = Math.min(c_support_lo, cIndex);
						c_support_hi = Math.max(c_support_hi, cIndex + 1);
					}
				}
			}
		}

		// Verbose output if desired
		
		if(D) {
			System.out.println(String.format("b=%.4g  magMain=%.4g  apcTot=%d  apcSup=%d",
				b, magMain, apc_total_size, apc_support_size));
			System.out.println(String.format("a: like=%.4g  mean=%.4g  sdev=%.4g  min=%.4g  max=%.4g  delta=%.4g  num=%d  lo=%d  hi=%d",
				stat_a_like, stat_a_mean, stat_a_sdev, min_a, max_a, delta_a, num_a, a_support_lo, a_support_hi));
			System.out.println(String.format("p: like=%.4g  mean=%.4g  sdev=%.4g  min=%.4g  max=%.4g  delta=%.4g  num=%d  lo=%d  hi=%d",
				stat_p_like, stat_p_mean, stat_p_sdev, min_p, max_p, delta_p, num_p, p_support_lo, p_support_hi));
			System.out.println(String.format("c: like=%.4g  mean=%.4g  sdev=%.4g  min=%.4g  max=%.4g  delta=%.4g  num=%d  lo=%d  hi=%d",
				stat_c_like, stat_c_mean, stat_c_sdev, min_c, max_c, delta_c, num_c, c_support_lo, c_support_hi));
		}

		return;
	}

	

	
	/**
	 * This computes the distribution of the number of M >= 5.0 events given all a, p, and c values, as well as the associated
	 * weight for each set of values.  This is used as a reference function that can be scaled to other magnitudes for greater
	 * efficiency.
	 * @param tMinDays = Beginning of the time interval, in days since the mainshock.
	 * @param tMaxDays = End of the time interval, in days since the mainshock.
	 * @return
	 * See comments for numMag5_DistributionFunc above.
	 * This uses only the elements in the support of apc_likelihood.
	 */
	public ArbDiscrEmpiricalDistFunc computeNumMag5_DistributionFunc(double tMinDays, double tMaxDays) {
		
		// If we already have a function computed for this time interval, then just return it

		if(tMinDaysCurrent == tMinDays && tMaxDaysCurrent == tMaxDays && numMag5_DistributionFunc != null) { // already computed
			return numMag5_DistributionFunc;
		}

		// Allocate a new function for this time interval
		
		tMinDaysCurrent = tMinDays;
		tMaxDaysCurrent = tMaxDays;
		numMag5_DistributionFunc = new ArbDiscrEmpiricalDistFunc();

		// Add points to the function, x = expected number of M5 aftershocks, y = probability of (a,p,c)

		for (int aIndex = a_support_lo; aIndex < a_support_hi; aIndex++) {
			for (int pIndex = p_support_lo; pIndex < p_support_hi; pIndex++) {
				for (int cIndex = c_support_lo; cIndex < c_support_hi; cIndex++) {
					if (apc_likelihood[aIndex][pIndex][cIndex] > apc_max_tail_element) {
						double numM5 = AftershockStatsCalc.getExpectedNumEvents(get_a(aIndex), b, magMain, 5.0, get_p(pIndex), get_c(cIndex), tMinDays, tMaxDays);
						numMag5_DistributionFunc.set(numM5, apc_likelihood[aIndex][pIndex][cIndex] / apc_support_total);
					}
				}
			}
		}

		// Debug or verbose output

		if(D) {
			System.out.println("M>=5 mean = "+numMag5_DistributionFunc.getMean());
			//System.out.println("M>=5 mode (caution) = "+numMag5_DistributionFunc.getApparentMode());	// don't do this because it can throw
			System.out.println("M>=5 median = "+numMag5_DistributionFunc.getMedian());
			System.out.println("M>=5 2.5 Percentile = "+numMag5_DistributionFunc.getInterpolatedFractile(0.025));
			System.out.println("M>=5 97.5 Percentile = "+numMag5_DistributionFunc.getInterpolatedFractile(0.975));
		}

		return numMag5_DistributionFunc;
	}
	


	
	/**
	 * This gives the expected number of aftershocks associated with the maximum likelihood a/p/c parameters (which represents the mode
	 * in the number of events space) above the given minimum magnitude and over the specified time span.  A GR distribution with
	 * no upper bound is assumed.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return
	 * According to R&J, the rate of aftershocks of magnitude >= magMin is
	 *  lambda(t) = k * (t + c)^(-p)
	 * where
	 *  k = 10^(a + b*(magMain - magMin))
	 *  Parmeters a, p, and c are set to their maximum-likelihood values.
	 *  Parameter b is set to the fixed value in this model.
	 * The value returned by this function is the integral of lambda(t) from t=tMin to t=tMax.
	 *
	 * Note (MJB): References here and elsewhere to "the mode in the number of events space"
	 * don't make sense to me.  There is a mapping from parameter triples (a,p,c) to expected
	 * number of events, and so the probability density on parameter triples induces a
	 * probability density on number of events.  But the mapping is non-linear and non-injective.
	 * So, the mode of the former probability density does not necessarily map into the mode
	 * of the latter probability density.  In other words: the most likely parameter triple
	 * (a,p,c) does not necessarily yield the most likely expected number of events.  Functions
	 * referring to "the mode in the number of events space" are actually computing values that
	 * correspond to the maximum likelihood (a,p,c) triple.
	 */
	public double getModalNumEvents(double magMin, double tMinDays, double tMaxDays) {
		return AftershockStatsCalc.getExpectedNumEvents(getMaxLikelihood_a(), b, magMain, magMin, getMaxLikelihood_p(), getMaxLikelihood_c(), tMinDays, tMaxDays);
	}
	


	
	/**
	 * This gives the probability of one or more aftershocks above the given minimum magnitude
	 * and over the specified time span.  A GR distribution with no upper bound is assumed.
	 * Epistemic uncertainty is included.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return
	 * According to R&J, the rate of aftershocks of magnitude >= magMin is
	 *  lambda(t) = k * (t + c)^(-p)
	 * where
	 *  k = 10^(a + b*(magMain - magMin))
	 *  Parmeters a, p, and c vary over their support.
	 *  Parameter b is set to the fixed value in this model.
	 * For given a/p/c, the number of aftershocks is assumed to be a Poisson distribution with
	 * expected value N(a,p,c) equal to the integral of lambda(t) from t=tMin to t=tMax.
	 * Then the probability P(a,p,c) of one or more aftershocks is 1 - exp(-N(a,p,c)).
	 * This function returns the weighted average of P(a,p,c) over all a/p/c triples,
	 * which represents the epistemic uncertainty.
	 */
	public double getProbOneOrMoreEvents(double magMin, double tMinDays, double tMaxDays) {
		double result = 0.0;

		for (int aIndex = a_support_lo; aIndex < a_support_hi; aIndex++) {
			for (int pIndex = p_support_lo; pIndex < p_support_hi; pIndex++) {
				for (int cIndex = c_support_lo; cIndex < c_support_hi; cIndex++) {
					if (apc_likelihood[aIndex][pIndex][cIndex] > apc_max_tail_element) {
						double expectedVal = AftershockStatsCalc.getExpectedNumEvents(get_a(aIndex), b, magMain, magMin, get_p(pIndex), get_c(cIndex), tMinDays, tMaxDays);
						double poissonProb = 1.0 - Math.exp(-expectedVal);
						result += (poissonProb * apc_likelihood[aIndex][pIndex][cIndex] / apc_support_total);
					}
				}
			}
		}

		if (result > 1.0) {
			result = 1.0;		// in case rounding produces a result a little larger than 1.0
		}

		return result;
	}
	


	
	/**
	 * This gives the expected cumulative number of aftershocks as a function of time for the maximum likelihood a/p/c parameters 
	 * (which represents the mode in the number of events space) above the given minimum magnitude and over 
	 * the specified time span.  A GR distribution with no upper bound is assumed.
	 * @param magMin = Minimum magnitude of aftershocks to consider.
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @param tDelta = Spacing between time values in the returned function.
	 * @return
	 * Returns a discrete function with:
	 *  x = time (in days after mainshock)
	 *  y = cumulative expected number of aftershocks for the corresponding time interval.
	 * The range [tMinDays, tMaxDays] is partitioned into equal-sized intervals, with the
	 * width of each interval equal to approximately tDelta.  The interval width is
	 * adjusted so that a whole number of intervals fit within the range [tMinDays, tMaxDays].
	 * For each such interval, the x value is the midpoint of the interval (and so
	 * the first and last x values are tMin+x_delta/2 and tMax-x_delta/2,
	 * where x_delta is the adjusted interval width).  The y value is the expected
	 * number of aftershocks within that interval plus all preceding intervals according
	 * to the R&J formula.  In evaluating the R&J formula, parameters a, p, and c are set
	 * to their maximum-likelihood values, and parameter b is the fixed value in this model.
	 */
	public EvenlyDiscretizedFunc getModalCumNumEventsWithTime(double magMin, double tMinDays, double tMaxDays, double tDelta) {
		return AftershockStatsCalc.getExpectedCumulativeNumWithTimeFunc(getMaxLikelihood_a(), b, magMain, magMin, 
				getMaxLikelihood_p(), getMaxLikelihood_c(), tMinDays, tMaxDays, tDelta);
	}



	
	/**
	 * This returns the cumulative MFD associated with the maximum likelihood a/p/c parameters (which represents the mode
	 * in the number of events space) and over the specified time span.   A GR distribution with no upper bound is assumed.
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return
	 * Returns a discrete function containing the cumulative magnitude-frequency distribution:
	 *  x = Magnitude.
	 *  y = Expected number of aftershocks with magnitude >= x, in the given time span.
	 * The discrete function contains evenly-spaced x values, with the smallest x-value equal
	 * to minMag and the largest x-value equal to maxMag.  The number of x-values is numMag.
	 * If numMag==1 then minMag and maxMag must be exactly equal.
	 * The expected number of aftershocks is computed using the R&J formula, with a, p, and c
	 * set to their maximum-likelihood values, and b set to the fixed value in this model.
	 *
	 * Note: The result will be a Gutenberg-Richter distribution, with y proportional to 10^(-b*x).
	 */
	public EvenlyDiscretizedFunc getModalCumNumMFD(double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc mfd = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
		for(int i=0;i<mfd.size();i++) {
			mfd.set(i, getModalNumEvents(mfd.getX(i), tMinDays, tMaxDays));
		}
		mfd.setName("Modal Num Events");
		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
//		double totExpNum = getExpectedNumEvents(minMag, tMinDays, tMaxDays);
//		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(b, totExpNum, minMag, maxMag, numMag);
//		mfd.setName("Expected Num Incr. MFD");
//		mfd.setInfo("Total Expected Num = "+totExpNum);
		return mfd;
	}

	
	
	
	/**
	 * This returns the mean cumulative MFD given each a/p/c parameter set, and the associated likelihood, 
	 * for the specified time span.   A GR distribution with no upper bound is assumed.
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return
	 * Returns a discrete function containing the cumulative magnitude-frequency distribution:
	 *  x = Magnitude.
	 *  y = Scaled expected number of aftershocks with magnitude >= x, in the given time span.
	 * The discrete function contains evenly-spaced x values, with the smallest x-value equal
	 * to minMag and the largest x-value equal to maxMag.  The number of x-values is numMag.
	 * If numMag==1 then minMag and maxMag must be exactly equal.
	 * The expected number of aftershocks is computed using the R&J formula, with a, p, and c
	 * set to their maximum-likelihood values, and b set to the fixed value in this model.
	 * Then, all the numbers of aftershocks are scaled, with a scale factor chosen so that the
	 * expected number of magnitude >= 5 aftershocks equals the mean of numMag5_DistributionFunc
	 * (which is the mean expected number of magnitude >= 5 aftershocks implied by the probability
	 * distribution of the parameters (a,p,c)).
	 * The choice of magnitude 5 is arbitrary; any other choice of magnitude would yield the
	 * same result; which implies that every value in this MFD equals the mean expected number
	 * of aftershocks implied by the probability distribution of (a,p,c).
	 *
	 * Note: The result will be a Gutenberg-Richter distribution, with y proportional to 10^(-b*x).
	 * The constant of proportionality is chosen so that each y is the mean expected number of
	 * aftershocks with magnitude >= x, according to the probability distribution on (a,p,c).
	 */
	public EvenlyDiscretizedFunc getMeanCumNumMFD(double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		// get the reference MFD, which we will scale to the mean
		EvenlyDiscretizedFunc mfd = getModalCumNumMFD(minMag, maxMag, numMag, tMinDays, tMaxDays);

//		double m5val = mfd.getInterpolatedY(5.0);	// fails if minMag > 5 || maxMag < 5
		double m5val = getModalNumEvents(5.0, tMinDays, tMaxDays);

		computeNumMag5_DistributionFunc(tMinDays, tMaxDays);
		mfd.scale(numMag5_DistributionFunc.getMean()/m5val);	// scale MFD to the mean at M5
		mfd.setName("Mean Num Events");
		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		return mfd;
	}



	
	/**
	 * This returns the fractile MFD implied by each a/p/c parameter set and their associated likelihoods, 
	 * for the specified time span.   A GR distribution with no upper bound is assumed.  Only epistemic uncertainty
	 * is considered (no aleatory variability).
	 * @param fractile - the fractile (percentile/100) for the distribution
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return
	 * Returns a discrete function containing the fractile magnitude-frequency distribution:
	 *  x = Magnitude.
	 *  y = Scaled expected number of aftershocks with magnitude >= x, in the given time span.
	 * The discrete function contains evenly-spaced x values, with the smallest x-value equal
	 * to minMag and the largest x-value equal to maxMag.  The number of x-values is numMag.
	 * If numMag==1 then minMag and maxMag must be exactly equal.
	 * The expected number of aftershocks is computed using the R&J formula, with a, p, and c
	 * set to their maximum-likelihood values, and b set to the fixed value in this model.
	 * Then, all the numbers of aftershocks are scaled, with a scale factor chosen so that the
	 * expected number of magnitude >= 5 aftershocks equals the fractile of numMag5_DistributionFunc
	 * (which is the fractile expected number of magnitude >= 5 aftershocks implied by the probability
	 * distribution of the parameters (a,p,c)).
	 * In other words:  Given the probability distribution on (a,p,c), the parameter "fractile"
	 * is the probability that the expected number of aftershocks of magnitude >= x is <= y.
	 * The choice of magnitude 5 is arbitrary; any other choice of magnitude would yield the
	 * same result; which implies that every value in this MFD equals the fractile expected
	 * number of aftershocks implied by the probability distribution of (a,p,c).
	 *
	 * Note: The result will be a Gutenberg-Richter distribution, with y proportional to 10^(-b*x).
	 * The constant of proportionality is chosen so that each y is the fractile expected number of
	 * aftershocks with magnitude >= x, according to the probability distribution on (a,p,c).
	 */
	public EvenlyDiscretizedFunc getCumNumMFD_Fractile(double fractile, double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		// get the reference MFD, which we will scale
		EvenlyDiscretizedFunc mfd = getModalCumNumMFD(minMag, maxMag, numMag, tMinDays, tMaxDays);

		// get the modal value at M5
//		double m5val;
//		if (minMag > 5 || maxMag < 5)	// in case requested range does not include M5
//			m5val = getModalCumNumMFD(5d, 5d, 1, tMinDays, tMaxDays).getY(0);
//		else
//			m5val = mfd.getInterpolatedY(5.0);
		double m5val = getModalNumEvents(5.0, tMinDays, tMaxDays);

		computeNumMag5_DistributionFunc(tMinDays, tMaxDays);
		mfd.scale(numMag5_DistributionFunc.getInterpolatedFractile(fractile)/m5val);
		mfd.setName(fractile+" Fractile for Num Events");
		mfd.setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		return mfd;
	}



	
	/**
	 * This returns the fractile MFDs implied by each a/p/c parameter set and their associated likelihoods, 
	 * for the specified time span.   A GR distribution with no upper bound is assumed.  Both epistemic uncertainty
	 * and aleatory variability are considered.
	 * @param fractileArray - Desired fractiles (percentile/100) of the probability distribution.
	 * @param minMag - the minimum magnitude considered
	 * @param maxMag - the maximum magnitude considered
	 * @param numMag - number of mags in the MFD
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return EvenlyDiscretizedFunc[]
	 * The return value is an array whose length equals fractileArray.length.
	 * Each element is a discrete function containing the fractile magnitude-frequency distribution:
	 *  x = Magnitude.
	 *  y = Fractile expected number of aftershocks with magnitude >= x, in the given time span.
	 * Specifically: y is an integer value, which is the smallest number such that the probability
	 * of there being y or fewer aftershocks of magnitude >= x is >= fractileArray[i].
	 * More succinctly: y is the fractileArray[i] fractile of the number of aftershocks of magnitude >= x.
	 * This computation takes into account both epistemic uncertainty (represented by the probability
	 * distribution on the parameter triple (a,p,c)) and aleatory uncertainty (represented by a
	 * Poisson distribution whose expected value equals the R&J expected number of aftershocks).
	 * See getCumNumFractileWithAleatory for further info on the probability distribution.
	 * Each discrete function contains evenly-spaced x values, with the smallest x-value equal
	 * to minMag and the largest x-value equal to maxMag.  The number of x-values is numMag.
	 * If numMag==1 then minMag and maxMag must be exactly equal.
	 */
	public EvenlyDiscretizedFunc[] getCumNumMFD_FractileWithAleatoryVariability(double[] fractileArray, double minMag, double maxMag, int numMag, double tMinDays, double tMaxDays) {
		EvenlyDiscretizedFunc[] mfdArray = new EvenlyDiscretizedFunc[fractileArray.length];
		for(int i=0;i<fractileArray.length;i++) {
			mfdArray[i] = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
			mfdArray[i].setName(fractileArray[i]+" Fractile for Num Events, including aleatory variability");
			mfdArray[i].setInfo("Cumulative distribution (greater than or equal to each magnitude)");
		}
		for(int i=0;i<numMag;i++) {
			double mag = mfdArray[0].getX(i);	// any MFD will do, as they all have the same x-axis values
			double[] valsArray = getCumNumFractileWithAleatory(fractileArray, mag, tMinDays, tMaxDays);
			for(int j=0;j<fractileArray.length;j++) {
				mfdArray[j].set(i,valsArray[j]);
//				System.out.println("\tworking on "+mag);
			}

		}
		return mfdArray;
	}



	
	/**
	 * This provides the cumulative number for the given fractiles, where aleatory variability
	 * is included in the result based on a Poisson distribution.
	 * @param fractileArray = Desired fractiles (percentile/100) of the probability distribution.
	 * @param mag = Minimum magnitude of aftershocks considered.
	 * @param tMinDays = Start of time range, in days after the mainshock.
	 * @param tMaxDays = End of time range, in days after the mainshock.
	 * @return
	 * This function constructs a probability distribution for number of aftershocks
	 * of magnitude mag or greater, in the time interval from tMinDays to tMagDays,
	 * which combines epistemic and aleatory uncertainty.
	 * Epistemic uncertainty is represented by numMag5_DistributionFunc, which is the
	 * probability distribution of the R&J expected number of aftershocks induced by the
	 * probability distribution of the parameter triple (a,p,c).
	 * Aleatory uncertainty is represented by a Poisson distribution whose expected value
	 * equals the R&J expected number of aftershocks.
	 * This function constructs a Poisson distribution for each set of parameter values (a,p,c),
	 * and then sums all the Poisson distributions with weight equal to the probability
	 * (i.e., likelihood) of the correponding (a,p,c) triple.
	 * The return value is an array whose length equals fractileArray.length.
	 * The i-th element of the return value is the fractileArray[i] fractile of the
	 * probability distribution, which is defined to be the minimum number n of aftershocks
	 * such that the probability of n or fewer aftershocks is >= fractileArray[i].
	 * Note that, although the return type is double[], the return values are integers.
	 *
	 * Implementation notes:
	 * Documentation for PoissonDistribution is here:
	 * http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/distribution/PoissonDistribution.html
	 * Since we do not generate samples from PoissonDistribution, it is better to use the constructor:
	 * PoissonDistribution(RandomGenerator rng, double p, double epsilon, int maxIterations)
	 * Set rng = null so that it does not create a random number generator.
	 * Set epsilon and maxIterations to their default values.
	 * However, this constructor is apparently planned for removal from Apache Math 4.0,
	 * so leave the original constructor commented-out in case it needs to be restored.
	 */
	public double[] getCumNumFractileWithAleatory(double[] fractileArray, double mag, double tMinDays, double tMaxDays) {
		// compute the distribution for the expected num aftershocks with M >= 5 (which we will scale to other magnitudes)
		computeNumMag5_DistributionFunc(tMinDays, tMaxDays);

		// get the maximum expected num, which we will use to set the maximum num in the distribution function
//System.out.print("\tworking on M "+mag+"\nunm="+numMag5_DistributionFunc.size()+"\n");
		double maxExpNum = numMag5_DistributionFunc.getMaxX()*Math.pow(10d, b*(5-mag));

//		PoissonDistribution poissDist = new PoissonDistribution(maxExpNum);
		PoissonDistribution poissDist = new PoissonDistribution(null, maxExpNum, PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		int maxAleatoryNum = poissDist.inverseCumulativeProbability(0.999);
		
		HistogramFunction cumDistFunc = new HistogramFunction(0d, (double)maxAleatoryNum,maxAleatoryNum+1);
		double[] distFunc = new double[cumDistFunc.size()];
		double totWt=0;
		
		for(int i=0;i<numMag5_DistributionFunc.size();i++) {
//System.out.print(", "+i);
			double expNum = numMag5_DistributionFunc.getX(i)*Math.pow(10d, b*(5-mag));
			double wt = numMag5_DistributionFunc.getY(i);
//			poissDist = new PoissonDistribution(expNum);
			poissDist = new PoissonDistribution(null, expNum, PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
			totWt+=wt;
			
			int minLoopVal = poissDist.inverseCumulativeProbability(0.0001);
			int maxLoopVal = poissDist.inverseCumulativeProbability(0.9999);
			if(maxLoopVal>cumDistFunc.size()-1)
				maxLoopVal=cumDistFunc.size()-1;
			if(minLoopVal < 0)
				minLoopVal = 0;
			for(int j=minLoopVal;j<=maxLoopVal;j++) {
				distFunc[j] += poissDist.probability(j)*wt;
			}
			
		}
		double sum=0;
		for(int j=0;j<distFunc.length;j++) {
			sum+=distFunc[j];
			cumDistFunc.set(j,sum);
		}
//System.out.print("\n");
		double[] fractValArray = new double[fractileArray.length];
		for(int i=0;i<fractileArray.length;i++) {
			double fractVal = (int)Math.round(cumDistFunc.getClosestXtoY(fractileArray[i]));
			if(cumDistFunc.getY(fractVal)<fractVal)
				fractVal += 1;	// this is how PoissonDistribution class does it	
			fractValArray[i]=fractVal;
		}

		
//		System.out.println("totWt="+totWt);
//		System.out.println("cumDistFunc.getMaxY()="+cumDistFunc.getMaxY());
//		System.out.println("fractVal="+fractVal+"\tfractile="+fractile);
//		GraphWindow graph = new GraphWindow(cumDistFunc, "cumDistFunc"); 

		return fractValArray;
	}



	
	/**
	 * This returns the PDF of a, which is a marginal distribution if either c or p 
	 * are unconstrained (either num_p or num_c not equal to 1). Null is returned if
	 * a is constrained (num_a=1).
	 * @return
	 */
	public HistogramFunction getPDF_a() {
		if(num_a == 1) {
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_a, num_a, delta_a);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						hist.add(get_a(aIndex), apc_likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
			String name = "PDF of a-value";
			if(num_p !=1 || num_c != 1)
				name += " (marginal)";
			hist.setName(name);
//			if(D) {
//				System.out.println("PDF of a-value:  "+hist);
//				System.out.println("PDF of a-value: totalTest = "+hist.calcSumOfY_Vals());
//			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}



	
	/**
	 * This returns the PDF of p, which is a marginal distribution if either a or c 
	 * are unconstrained (either num_a or num_c not equal to 1). Null is returned if
	 * p is constrained (num_p=1).
	 * @return
	 */
	public HistogramFunction getPDF_p() {
		if(num_p == 1) {
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_p, num_p, delta_p);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						hist.add(get_p(pIndex), apc_likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
			String name = "PDF of p-value";
			if(num_a !=1 || num_c != 1)
				name += " (marginal)";
			hist.setName(name);
//			if(D) {
//				System.out.println("PDF of p-value: totalTest = "+hist.calcSumOfY_Vals());
//			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}



	
	/**
	 * This returns the PDF of c, which is a marginal distribution if either a or p 
	 * are unconstrained (either num_a or num_p not equal to 1). Null is returned if
	 * c is constrained (num_c=1).
	 * @return
	 */
	public HistogramFunction getPDF_c() {
		if(num_c == 1) {
			return null;
		}
		else {
			HistogramFunction hist = new HistogramFunction(min_c, num_c, delta_c);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						hist.add(get_c(cIndex), apc_likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
			String name = "PDF of c-value";
			if(num_a !=1 || num_p != 1)
				name += " (marginal)";
			hist.setName(name);
//			if(D) {
//				System.out.println("PDF of c-value: totalTest = "+hist.calcSumOfY_Vals());
//			}
			hist.scale(1d/hist.getDelta());
			return hist;
		}
	}



	
	/**
	 * This returns a 2D PDF for a and p, which is a marginal distribution if c 
	 * is unconstrained (num_c not equal to 1). Null is returned if either
	 * a or p are constrained (num_a=1 or num_p=1).
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_a_and_p() {
		if(num_a == 1 || num_p == 1) {
			return null;
		}
		else {
			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_a, num_p, min_a, min_p, delta_a, delta_p);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						double prevVal = hist2D.get(aIndex,pIndex);
						hist2D.set(aIndex,pIndex, prevVal+apc_likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
//			String name = "2D PDF of a vs p";
//			if(num_c != 1)
//				name += " (marginal)";
//			if(D) {
//				System.out.println("2D PDF of a vs p: totalTest = "+hist2D.getSumZ());
//			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}

	

	
	/**
	 * This returns a 2D PDF for a and c, which is a marginal distribution if p 
	 * is unconstrained (num_p not equal to 1). Null is returned if either
	 * a or c are constrained (num_a=1 or num_c=1).
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_a_and_c() {
		if(num_a == 1 || num_c == 1) {
			return null;
		}
		else {
			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_a, num_c, min_a, min_c, delta_a, delta_c);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						double prevVal = hist2D.get(aIndex,cIndex);
						hist2D.set(aIndex,cIndex, prevVal+apc_likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
//			String name = "2D PDF of a vs c";
//			if(num_p != 1)
//				name += " (marginal)";
//			if(D) {
//				System.out.println("2D PDF of a vs c: totalTest = "+hist2D.getSumZ());
//			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}




	/**
	 * This returns a 2D PDF for c and p, which is a marginal distribution if a 
	 * is unconstrained (num_a not equal to 1). Null is returned if either
	 * c or p are constrained (num_c=1 or num_p=1).
	 * @return
	 */
	public EvenlyDiscrXYZ_DataSet get2D_PDF_for_c_and_p() {
		if(num_c == 1 || num_p == 1) {
			return null;
		}
		else {
			EvenlyDiscrXYZ_DataSet hist2D = new EvenlyDiscrXYZ_DataSet(num_c, num_p, min_c, min_p, delta_c, delta_p);
			for(int aIndex=0;aIndex<num_a;aIndex++) {
				for(int pIndex=0;pIndex<num_p;pIndex++) {
					for(int cIndex=0;cIndex<num_c;cIndex++) {
						double prevVal = hist2D.get(cIndex,pIndex);
						hist2D.set(cIndex,pIndex, prevVal+apc_likelihood[aIndex][pIndex][cIndex]);
					}
				}
			}
//			String name = "2D PDF of c vs p";
//			if(num_a != 1)
//				name += " (marginal)";
//			if(D) {
//				System.out.println("2D PDF of c vs p: totalTest = "+hist2D.getSumZ());
//			}
			hist2D.scale(1d/(hist2D.getGridSpacingX()*hist2D.getGridSpacingY()));
			return hist2D;
		}
	}



	
	/**
	 * This sets the apc_likelihood and maximum likelihood values from the a-values in the given discretized function,
	 *  and holding the other parameters fixed at the values given.  The apc_likelihood is normalized so values sum
	 *  to 1.0.
	 * @param aValueFunc = Function giving likelihood of Reasenberg-Jones productivity parameter a at equally-spaced values.
	 * @param p = Omori p-parameter (exponent).
	 * @param c = Omori c-parameter (time offset), in days.
	 * Note: The caller must set up this.b and this.magMain.
	 * This function sets up all the fields related to a, p, and c, including apc_likelihood.
	 * The likelihood values in aValueFunc need not be normalized.
	 * The values of p and c are considered to be known.
	 */
	protected void setArrayAndMaxLikelyValuesFrom_aValueFunc(EvenlyDiscretizedFunc aValueFunc, double p, double c) {
		
		delta_a = aValueFunc.getDelta();
		min_a = aValueFunc.getMinX();
		max_a = aValueFunc.getMaxX();
		num_a = aValueFunc.size();

		set_fixed_p(p);
		set_fixed_c(c);
		
		apc_likelihood = new double[num_a][num_p][num_c];
		for (int aIndex = 0; aIndex < num_a; aIndex++) {
			double wt = aValueFunc.getY(aIndex);
			apc_likelihood[aIndex][0][0] = wt;
		}

		// Complete the likelihood setup

		apcFinish (false);
		
		return;
	}




	public static void main(String[] args) {
		// TODO Auto-generated method stub

		return;
	}

}
