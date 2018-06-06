package scratch.aftershockStatistics;

import java.util.ArrayList;
import java.util.List;

//import org.mongodb.morphia.annotations.Transient;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;
import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;


/**
 * This computes a Reasenberg-Jones (1989, 1994) aftershock model from aftershock data and with an assumed time-dependent 
 * magnitude completeness model described in the constructor.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests in order to cover all cases.
 *
 * @author field
 *
 * Modified by Michael Barall.
 *
 * According to R&J, the rate of aftershocks of magnitude >= magMin is
 *  lambda(t) = k * (t + c)^(-p)
 * where
 *  k = 10^(a + b*(magMain - magMin))
 * According to Page et al. the time-dependent magnitude of completeness is
 *  magMin(t) = Max(magMain/2 - G - H*log10(t), magCat)
 * In these formulas, t is measured in days.
 *
 * Suppose that in a time interval tMinDays <= t <= tMaxDays, there is a sequence of aftershocks
 * at times t_i with magnitudes M_i, for 1 <= i <= N.  Assume that M_i >= magMin(t_i).
 * Then, for a parameter triple (a,p,c), the log-likelihood is defined as:
 *  log L(a,p,c) = SUM(log(lambda(t_i)), 1 <= i <= N) - INTEGRAL(lambda(t)*dt, tMinDays <= t <= tMaxDays)
 *
 * This model uses the likelihood L(a,p,c) as a probability density for the parameter triple (a,p,c).
 */
public class RJ_AftershockModel_SequenceSpecific extends RJ_AftershockModel {

//	//@Transient
//	boolean D=true;	// debug flag (inherited)

	// Parameters for the time-dependent magnitude of completeness.

	protected double capG = 10.0;
	protected double capH = 0.0;
	protected double magCat = 0.0;

	// The list of aftershocks used to construct the model.

	//@Transient
	protected List<ObsEqkRupture> aftershockList = null;

	// The mainshock.

	//@Transient
	protected ObsEqkRupture mainShock = null;

	// The time interval covered by the list of aftershocks, in days since the mainshock.

	protected double dataStartTimeDays = 0.0;
	protected double dataEndTimeDays = 0.0;

	// The number of aftershocks that were used to determine parameters.

	protected int numAftershocks = 0;




	/**
	 * Return the name of this model.
	 */
	@Override
	public String getModelName() {
		return "Reasenberg-Jones (1989, 1994) aftershock model (Sequence Specific)";
	}



	
	/**
	 * Use this constructor to initialize from parameter holders.
	 * @param mainShock - the mainshock
	 * @param aftershockList - list of aftershocks; events with mag below magCat will be filtered out
	 * @param dataStartTimeDays - start time for data, in days since the mainshock
	 * @param dataEndTimeDays - end time for data, in days since the mainshock
	 * @param mcParam - magnitude of completeness parameters
	 * @param sqParam - sequence-specific range parameters
	 */
	public RJ_AftershockModel_SequenceSpecific(ObsEqkRupture mainShock, List<ObsEqkRupture> aftershockList,
				double dataStartTimeDays, double dataEndTimeDays,
				MagCompPage_Parameters mcParam, SeqSpecRJ_Parameters sqParam) {
		
		this(mainShock, aftershockList,
				mcParam.get_magCat(), mcParam.get_capG(), mcParam.get_capH(),
				sqParam.get_b(), dataStartTimeDays, dataEndTimeDays,
				sqParam.get_min_a(), sqParam.get_max_a(), sqParam.get_num_a(),
				sqParam.get_min_p(), sqParam.get_max_p(), sqParam.get_num_p(),
				sqParam.get_min_c(), sqParam.get_max_c(), sqParam.get_num_c());

	}



	
	/**
	 * Use this constructor to apply a time-independent magnitude of completeness.
	 * @param mainShock - the mainshock
	 * @param aftershockList - list of aftershocks; events with mag below magCat will be filtered out
	 * @param magCat - the magnitude of completeness (time independent)
	 * @param b - assumed b value
	 * @param dataStartTimeDays - start time for data, in days since the mainshock
	 * @param dataEndTimeDays - end time for data, in days since the mainshock
	 * @param min_a \
	 * @param max_a  | - range of a-values for grid search (set min=max and num=1 to constrain to single value)
	 * @param num_a /
	 * @param min_p \
	 * @param max_p  | - range of p-values for grid search (set min=max and num=1 to constrain to single value)
	 * @param num_p /
	 * @param min_c \
	 * @param max_c  | - range of c-values for grid search (set min=max and num=1 to constrain to single value)
	 * @param num_c /
	 * Note: In previous versions, this constructor was faster than the constructor below because it
	 * used an analytic solution whereas the constructor below used numerical integration.  There is no
	 * longer any performance difference between the two constructors, because the code auto-detects
	 * when an analytic solution can be used.
	 */
	public RJ_AftershockModel_SequenceSpecific(ObsEqkRupture mainShock, List<ObsEqkRupture> aftershockList,
				double magCat, double b, double dataStartTimeDays, double dataEndTimeDays,
				double min_a, double max_a, int num_a, 
				double min_p, double max_p, int num_p, 
				double min_c, double max_c, int num_c) {
		
		this(mainShock, aftershockList, magCat, Double.NaN, Double.NaN, b, dataStartTimeDays, dataEndTimeDays,
				min_a, max_a, num_a, min_p, max_p, num_p, min_c, max_c, num_c);

	}

		

	
	/**
	 * This solves for the Reasenberg-Jones parameters from the given mainShock, aftershockList,
	 * and other specifications as described below, and for a time-dependent magnitude of completeness
	 * model defined as Mc(t,Mm) = Max(Mm/2-G-H*log10(t); Mcat), where Mm is the main-shock magnitude, G and H 
	 * are model parameters, and Mcat is the magnitude of completeness for the network during normal times. 
	 * Likelihood values are normalized so they sum to 1.0 over the range of parameter-values specified.
	 * @param mainShock - the mainshock
	 * @param aftershockList - list of aftershocks; events with mag below magCat will be filtered out
	 * @param magCat - "Mcat" in the in the time-dependent magnitude of completeness model defined above
	 * @param capG - the "G" parameter in the time-dependent magnitude of completeness model defined above; 
	 *               As a special case, if capG == 10.0 then the magnitude of completeness is always magCat.
	 * @param capH - the "H" parameter in the time-dependent magnitude of completeness model defined above
	 * @param b - assumed b value
	 * @param dataStartTimeDays - start time for data, in days since the mainshock
	 * @param dataEndTimeDays - end time for data, in days since the mainshock
	 * @param min_a \
	 * @param max_a  | - range of a-values for grid search (set min=max and num=1 to constrain to single value)
	 * @param num_a /
	 * @param min_p \
	 * @param max_p  | - range of p-values for grid search (set min=max and num=1 to constrain to single value)
	 * @param num_p /
	 * @param min_c \
	 * @param max_c  | - range of c-values for grid search (set min=max and num=1 to constrain to single value)
	 * @param num_c /
	 * Note: For compatibility, if either capG or capH is equal to Double.NaN, then it is treated
	 * as if capG==10.0 and capH==0.0, that is, the magnitude of completeness is always magCat.
	 * New code should not rely on this behavior.
	 */
	public RJ_AftershockModel_SequenceSpecific(ObsEqkRupture mainShock, List<ObsEqkRupture> aftershockList,
			 								double magCat, double capG, double capH,
											double b, double dataStartTimeDays, double dataEndTimeDays,
											double min_a, double max_a, int num_a, 
											double min_p, double max_p, int num_p, 
											double min_c, double max_c, int num_c) {
		
		// check range values
		if(num_a == 1 && min_a != max_a) {
			throw new RuntimeException("RJ_AftershockModel_SequenceSpecific: num_a == 1 && min_a != max_a");
		}
		if(num_p == 1 && min_p != max_p) {
			throw new RuntimeException("RJ_AftershockModel_SequenceSpecific: num_p == 1 && min_p != max_p");
		}
		if(num_c == 1 && min_c != max_c) {
			throw new RuntimeException("RJ_AftershockModel_SequenceSpecific: num_c == 1 && min_c != max_c");
		}
		if(min_a > max_a) {
			throw new RuntimeException("RJ_AftershockModel_SequenceSpecific: min_a > max_a");
		}
		if(min_p > max_p) {
			throw new RuntimeException("RJ_AftershockModel_SequenceSpecific: min_p > max_p");
		}
		if(min_c > max_c) {
			throw new RuntimeException("RJ_AftershockModel_SequenceSpecific: min_c > max_c");
		}

		this.min_a = min_a;
		this.max_a = max_a;
		this.num_a = num_a;
		this.min_p = min_p;
		this.max_p = max_p;
		this.num_p = num_p;
		this.min_c = min_c;
		this.max_c = max_c;
		this.num_c = num_c;
		this.b = b;
		this.magCat = magCat;

		this.aftershockList=aftershockList;
//		this.aftershockList = new ObsEqkRupList();
		this.mainShock=mainShock;
		this.dataStartTimeDays=dataStartTimeDays;
		this.dataEndTimeDays=dataEndTimeDays;

		if(Double.isNaN(capG) || Double.isNaN(capH)) {
			this.capG = 10.0;
			this.capH = 0.0;
		} else {
			this.capG = capG;
			this.capH = capH;
		}
		
		this.magMain = mainShock.getMag();

		if(num_a>1) {
			this.delta_a = (max_a-min_a)/((double)num_a - 1.0);
		} else {
			this.delta_a = 0.0;
		}
		if(num_p>1) {
			this.delta_p = (max_p-min_p)/((double)num_p - 1.0);
		} else {
			this.delta_p = 0.0;
		}
		if(num_c>1) {
			this.delta_c = (max_c-min_c)/((double)num_c - 1.0);
		} else {
			this.delta_c = 0.0;
		}

		apc_build(mainShock, aftershockList, dataStartTimeDays, dataEndTimeDays);
		
	}




	/**
	 * This default constructor creates an empty model.
	 * This is intended for use in database retrieval.
	 */
    public RJ_AftershockModel_SequenceSpecific() {
		// When retrieving from database, remain quiet by default
		D = false;
    }




	/**
	 * Build the apc_likelihood matrix, that gives the probability distribution of (a,p,c).
	 * @param mainShock - the mainshock
	 * @param aftershockList - list of aftershocks; events with mag below magCat will be filtered out
	 * @param dataStartTimeDays - start time for data, in days since the mainshock
	 * @param dataEndTimeDays - end time for data, in days since the mainshock
	 */
    public void apc_build(ObsEqkRupture mainShock, List<ObsEqkRupture> aftershockList, double dataStartTimeDays, double dataEndTimeDays) {

		// Save the parameters

		this.aftershockList = aftershockList;
		this.mainShock = mainShock;
		this.dataStartTimeDays = dataStartTimeDays;
		this.dataEndTimeDays = dataEndTimeDays;

		this.numAftershocks = 0;

		// Allocate the array

		apc_likelihood = new double[num_a][num_p][num_c];
		double ln10 = Math.log(10);

		// Loop over c first, so we can accumulate log(t+c)

		for(int cIndex = 0; cIndex < num_c; cIndex++) {
			double c = get_c(cIndex);

			// Sum of magMain - magMin(t_i)

			double sum1 = 0.0;

			// Sum of log(t_i + c)

			double sum2 = 0.0;

			// Number of aftershocks

			int numEvents = 0;

			// Scan list of aftershocks

			for(ObsEqkRupture rup:aftershockList) {

				// Get time since the mainshock in days, skip it if it is outside our time interval

				double timeSinceMainDays = (double)(rup.getOriginTime()-mainShock.getOriginTime()) / (double)AftershockStatsCalc.MILLISEC_PER_DAY;
				if(timeSinceMainDays < dataStartTimeDays || timeSinceMainDays > dataEndTimeDays) { // not necessary if list already filtered
					continue;
				}

				// Get the magnitude of completeness at this time

				double magMin = AftershockStatsCalc.getPageMagCompleteness(
									magMain, magCat, capG, capH, timeSinceMainDays);

				// If the aftershock magnitude is at least the magnitude of completeness, accumulate it

				if(rup.getMag() >= magMin) {
					numEvents += 1;
					sum1 += (magMain - magMin);
					sum2 += Math.log(timeSinceMainDays + c);
					++numAftershocks;
				}
			}

			// Now loop over p and a

			for(int pIndex=0;pIndex<num_p;pIndex++) {
				double p = get_p(pIndex);
				for(int aIndex=0;aIndex<num_a;aIndex++) {
					double a = get_a(aIndex);

					// Compute the integral of the aftershock rate over the time interval

					double integral = AftershockStatsCalc.getPageExpectedNumEvents(
						a, b, magMain, magCat, capG, capH, p, c, dataStartTimeDays, dataEndTimeDays);

					// Form the log likelihood

					double logLike = numEvents*a*ln10 + b*ln10*sum1 - p*sum2 - integral;

					// Save it as the array element

					apc_likelihood[aIndex][pIndex][cIndex] = logLike;
				}
			}
		}

		// Complete the likelihood setup

		apcFinish (true);	// true means array contains log-likelihood

		if(D) {
			System.out.println(String.format("G=%.4g  H=%.4g  magCat=%.4g  tStart=%.8g  tEnd=%.8g  nEvents=%d",
				capG, capH, magCat, dataStartTimeDays, dataEndTimeDays, numAftershocks));
		}
		
		return;
	}




	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("RJ_AftershockModel_SequenceSpecific : Missing subcommand");
			return;
		}


		// Subcommand : Test #1
		// Command format:
		//  test1
		// Generate a simulated aftershock sequence.
		// Then, construct the model and see if it can recover the parameters a and p.

		if (args[0].equalsIgnoreCase ("test1")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("RJ_AftershockModel_SequenceSpecific : Invalid 'test1' subcommand");
				return;
			}

			// Parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
			double capG = 1.25;
			double capH = 0.75;
			double dataStartTimeDays = 0.0;
			double dataEndTimeDays = 30.0;
		
			double min_a = -2.0;
			double max_a = -1.0;
			int num_a = 101;

			double min_p = 0.9; 
			double max_p = 1.2; 
			int num_p = 31;
		
			double min_c=0.05;
			double max_c=0.05;
			int num_c=1;

			// Run the simulation

			ObsEqkRupList aftershockList = AftershockStatsCalc.simAftershockSequence(a, b, magMain, magCat, capG, capH, p, c, dataStartTimeDays, dataEndTimeDays);

			// Make the mainshock

			ObsEqkRupture mainShock = new ObsEqkRupture("0", 0L, null, magMain);

			// Make the model, it will output some information

			RJ_AftershockModel_SequenceSpecific gen =
				new RJ_AftershockModel_SequenceSpecific(mainShock, aftershockList,
			 								magCat, capG, capH,
											b, dataStartTimeDays, dataEndTimeDays,
											min_a, max_a, num_a, 
											min_p, max_p, num_p, 
											min_c, max_c, num_c);

			// A few calculations
		
			EvenlyDiscretizedFunc lowFract = gen.getCumNumMFD_Fractile(0.025, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("2.5%: "+lowFract.getX(0)+"\t"+lowFract.getY(0));
			EvenlyDiscretizedFunc hiFract = gen.getCumNumMFD_Fractile(0.975, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("97.5%: "+hiFract.getX(0)+"\t"+hiFract.getY(0));
		
			double[] fractArray = {0.025, 0.975};
			EvenlyDiscretizedFunc[] fractalWithAleatoryMFDArray = gen.getCumNumMFD_FractileWithAleatoryVariability(fractArray, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("2.5% With Aleatory: "+fractalWithAleatoryMFDArray[0].getX(0)+"\t"+fractalWithAleatoryMFDArray[0].getY(0));
			System.out.println("97.5% With Aleatory: "+fractalWithAleatoryMFDArray[1].getX(0)+"\t"+fractalWithAleatoryMFDArray[1].getY(0));

			return;
		}


		// Subcommand : Test #2
		// Command format:
		//  test2
		// Generate a simulated aftershock sequence.
		// Then, construct the model and see if it can recover the parameters a and p.
		// It's done twice, once with ObsEqkRupList and again with CompactEqkRupList.
		// The second time we use the new constructor that takes parameter holders.
		// Note: Most of the time the two sets of results are the same.  Occasionally the
		// second set has one less aftershock (because magnitude rounding causes one event just
		// above the magnitude of completeness to fall below it), causing the results to differ.

		if (args[0].equalsIgnoreCase ("test2")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("RJ_AftershockModel_SequenceSpecific : Invalid 'test2' subcommand");
				return;
			}

			// Parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
			double capG = 1.25;
			double capH = 0.75;
			double dataStartTimeDays = 0.0;
			double dataEndTimeDays = 30.0;
		
			double min_a = -2.0;
			double max_a = -1.0;
			int num_a = 101;

			double min_p = 0.9; 
			double max_p = 1.2; 
			int num_p = 31;
		
			double min_c=0.05;
			double max_c=0.05;
			int num_c=1;

			// Run the simulation

			ObsEqkRupList aftershockList = AftershockStatsCalc.simAftershockSequence(a, b, magMain, magCat, capG, capH, p, c, dataStartTimeDays, dataEndTimeDays);

			// Make the mainshock

			ObsEqkRupture mainShock = new ObsEqkRupture("0", 0L, null, magMain);

			// Make the model, it will output some information

			RJ_AftershockModel_SequenceSpecific gen =
				new RJ_AftershockModel_SequenceSpecific(mainShock, aftershockList,
			 								magCat, capG, capH,
											b, dataStartTimeDays, dataEndTimeDays,
											min_a, max_a, num_a, 
											min_p, max_p, num_p, 
											min_c, max_c, num_c);

			// A few calculations
		
			EvenlyDiscretizedFunc lowFract = gen.getCumNumMFD_Fractile(0.025, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("2.5%: "+lowFract.getX(0)+"\t"+lowFract.getY(0));
			EvenlyDiscretizedFunc hiFract = gen.getCumNumMFD_Fractile(0.975, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("97.5%: "+hiFract.getX(0)+"\t"+hiFract.getY(0));
		
			double[] fractArray = {0.025, 0.975};
			EvenlyDiscretizedFunc[] fractalWithAleatoryMFDArray = gen.getCumNumMFD_FractileWithAleatoryVariability(fractArray, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("2.5% With Aleatory: "+fractalWithAleatoryMFDArray[0].getX(0)+"\t"+fractalWithAleatoryMFDArray[0].getY(0));
			System.out.println("97.5% With Aleatory: "+fractalWithAleatoryMFDArray[1].getX(0)+"\t"+fractalWithAleatoryMFDArray[1].getY(0));

			// Compact form of list

			CompactEqkRupList compactList = new CompactEqkRupList (aftershockList);

			// Make the model, it will output some information

			//RJ_AftershockModel_SequenceSpecific compactGen =
			//	new RJ_AftershockModel_SequenceSpecific(mainShock, compactList,
			// 								magCat, capG, capH,
			//								b, dataStartTimeDays, dataEndTimeDays,
			//								min_a, max_a, num_a, 
			//								min_p, max_p, num_p, 
			//								min_c, max_c, num_c);

			MagCompPage_Parameters mcParam = new MagCompPage_Parameters (magCat, capG, capH);
			SeqSpecRJ_Parameters sqParam = new SeqSpecRJ_Parameters (b,
				min_a, max_a, num_a, min_p, max_p, num_p, min_c, max_c, num_c);

			RJ_AftershockModel_SequenceSpecific compactGen =
				new RJ_AftershockModel_SequenceSpecific(mainShock, compactList,
			 								dataStartTimeDays, dataEndTimeDays,
											mcParam, sqParam);

			// A few calculations
		
			EvenlyDiscretizedFunc lowFract2 = compactGen.getCumNumMFD_Fractile(0.025, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("2.5%: "+lowFract2.getX(0)+"\t"+lowFract2.getY(0));
			EvenlyDiscretizedFunc hiFract2 = compactGen.getCumNumMFD_Fractile(0.975, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("97.5%: "+hiFract2.getX(0)+"\t"+hiFract2.getY(0));
		
			double[] fractArray2 = {0.025, 0.975};
			EvenlyDiscretizedFunc[] fractalWithAleatoryMFDArray2 = compactGen.getCumNumMFD_FractileWithAleatoryVariability(fractArray2, 5.0, 8.0, 31, 0d, 7d);
			System.out.println("2.5% With Aleatory: "+fractalWithAleatoryMFDArray2[0].getX(0)+"\t"+fractalWithAleatoryMFDArray2[0].getY(0));
			System.out.println("97.5% With Aleatory: "+fractalWithAleatoryMFDArray2[1].getX(0)+"\t"+fractalWithAleatoryMFDArray2[1].getY(0));

			return;
		}


		// Subcommand : Test #3
		// Command format:
		//  test3
		// Generate a simulated aftershock sequence.
		// Convert to compact form.
		// Construct the parameter holder objects.
		// Construct generic, sequence-specific, and bayesian models.
		// Construct summary objects for each of the three models.
		// Marshal the summary objects.
		// Display summaries twice: once direct, once unmarshaled.

		if (args[0].equalsIgnoreCase ("test3")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("RJ_AftershockModel_SequenceSpecific : Invalid 'test3' subcommand");
				return;
			}

			// Parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
			double capG = 1.25;
			double capH = 0.75;
			double dataStartTimeDays = 0.0;
			double dataEndTimeDays = 30.0;
		
			double min_a = -2.0;
			double max_a = -1.0;
			int num_a = 101;

			double min_p = 0.9; 
			double max_p = 1.2; 
			int num_p = 31;
		
			double min_c=0.05;
			double max_c=0.05;
			int num_c=1;

			// Run the simulation

			ObsEqkRupList aftershockList = AftershockStatsCalc.simAftershockSequence(a, b, magMain, magCat, capG, capH, p, c, dataStartTimeDays, dataEndTimeDays);

			// Make the mainshock

			ObsEqkRupture mainShock = new ObsEqkRupture("0", 0L, null, magMain);

			// Compact form of list

			CompactEqkRupList compactList = new CompactEqkRupList (aftershockList);

			// Parameter holders

			double a_sigma = 1.76;
			double a_sigma1 = 750;
			double a_sigma0 = 0.49;
			double a_delta = (max_a - min_a)/((double)(num_a - 1));

			GenericRJ_Parameters rjParam = new GenericRJ_Parameters (a, a_sigma, a_sigma0, a_sigma1, 
					b, p, c, min_a, max_a, a_delta);

			MagCompPage_Parameters mcParam = new MagCompPage_Parameters (magCat, capG, capH);

			SeqSpecRJ_Parameters sqParam = new SeqSpecRJ_Parameters (b,
				min_a, max_a, num_a, min_p, max_p, num_p, min_c, max_c, num_c);

			// Make the models

			RJ_AftershockModel_SequenceSpecific seqModel =
				new RJ_AftershockModel_SequenceSpecific(mainShock, compactList,
			 								dataStartTimeDays, dataEndTimeDays,
											mcParam, sqParam);

			rjParam = new GenericRJ_Parameters (a, a_sigma, a_sigma0, a_sigma1, 
					b, seqModel.getMaxLikelihood_p(), c, min_a, max_a, a_delta);	// same p so bayesian model can be formed

			RJ_AftershockModel_Generic genModel = new RJ_AftershockModel_Generic (magMain, rjParam);

			RJ_AftershockModel_Bayesian bayModel = new RJ_AftershockModel_Bayesian(genModel, seqModel);

			// Make the summaries

			RJ_Summary_Generic genSummary = new RJ_Summary_Generic (genModel);

			RJ_Summary_SequenceSpecific seqSummary = new RJ_Summary_SequenceSpecific (seqModel);

			RJ_Summary_Bayesian baySummary = new RJ_Summary_Bayesian (bayModel);

			// Marshal the summaries

			MarshalImpArray store = new MarshalImpArray();
			store.marshalMapBegin (null);

			//RJ_Summary_Generic.marshal_poly (store, "Generic", genSummary);
			genSummary.marshal (store, "Generic");

			//RJ_Summary_SequenceSpecific.marshal_poly (store, "SeqSpec", seqSummary);
			seqSummary.marshal (store, "SeqSpec");

			//RJ_Summary_Bayesian.marshal_poly (store, "Bayesian", baySummary);
			baySummary.marshal (store, "Bayesian");

			// Prepare to unmarshal

			store.marshalMapEnd ();
			store.check_write_complete ();
			store.unmarshalMapBegin (null);

			// Write the summaries

			System.out.println ("Generic, Direct:\n" + genSummary.toString());

			//System.out.println ("Generic, Marshaled:\n" + RJ_Summary_Generic.unmarshal_poly(store, "Generic").toString());
			System.out.println ("Generic, Marshaled:\n" + (new RJ_Summary_Generic()).unmarshal(store, "Generic").toString());

			System.out.println ("Sequence-Specific, Direct:\n" + seqSummary.toString());

			//System.out.println ("Sequence-Specific, Marshaled:\n" + RJ_Summary_SequenceSpecific.unmarshal_poly(store, "SeqSpec").toString());
			System.out.println ("Sequence-Specific, Marshaled:\n" + (new RJ_Summary_SequenceSpecific()).unmarshal(store, "SeqSpec").toString());

			System.out.println ("Bayesian, Direct:\n" + baySummary.toString());

			//System.out.println ("Bayesian, Marshaled:\n" + RJ_Summary_Bayesian.unmarshal_poly(store, "Bayesian").toString());
			System.out.println ("Bayesian, Marshaled:\n" + (new RJ_Summary_Bayesian()).unmarshal(store, "Bayesian").toString());

			store.unmarshalMapEnd ();
			store.check_read_complete ();

			return;
		}


		// Subcommand : Test #4
		// Command format:
		//  test3
		// Generate a simulated aftershock sequence.
		// Convert to compact form.
		// Construct the parameter holder objects.
		// Construct generic, sequence-specific, and bayesian models.
		// Construct summary objects for each of the three models.
		// Marshal the summary objects.
		// Display summaries twice: once direct, once unmarshaled.
		// This version uses marshaling to JSON, and displays the JSON string.

		if (args[0].equalsIgnoreCase ("test4")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("RJ_AftershockModel_SequenceSpecific : Invalid 'test4' subcommand");
				return;
			}

			// Parameter values
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
			double capG = 1.25;
			double capH = 0.75;
			double dataStartTimeDays = 0.0;
			double dataEndTimeDays = 30.0;
		
			double min_a = -2.0;
			double max_a = -1.0;
			int num_a = 101;

			double min_p = 0.9; 
			double max_p = 1.2; 
			int num_p = 31;
		
			double min_c=0.05;
			double max_c=0.05;
			int num_c=1;

			// Run the simulation

			ObsEqkRupList aftershockList = AftershockStatsCalc.simAftershockSequence(a, b, magMain, magCat, capG, capH, p, c, dataStartTimeDays, dataEndTimeDays);

			// Make the mainshock

			ObsEqkRupture mainShock = new ObsEqkRupture("0", 0L, null, magMain);

			// Compact form of list

			CompactEqkRupList compactList = new CompactEqkRupList (aftershockList);

			// Parameter holders

			double a_sigma = 1.76;
			double a_sigma1 = 750;
			double a_sigma0 = 0.49;
			double a_delta = (max_a - min_a)/((double)(num_a - 1));

			GenericRJ_Parameters rjParam = new GenericRJ_Parameters (a, a_sigma, a_sigma0, a_sigma1, 
					b, p, c, min_a, max_a, a_delta);

			MagCompPage_Parameters mcParam = new MagCompPage_Parameters (magCat, capG, capH);

			SeqSpecRJ_Parameters sqParam = new SeqSpecRJ_Parameters (b,
				min_a, max_a, num_a, min_p, max_p, num_p, min_c, max_c, num_c);

			// Make the models

			RJ_AftershockModel_SequenceSpecific seqModel =
				new RJ_AftershockModel_SequenceSpecific(mainShock, compactList,
			 								dataStartTimeDays, dataEndTimeDays,
											mcParam, sqParam);

			rjParam = new GenericRJ_Parameters (a, a_sigma, a_sigma0, a_sigma1, 
					b, seqModel.getMaxLikelihood_p(), c, min_a, max_a, a_delta);	// same p so bayesian model can be formed

			RJ_AftershockModel_Generic genModel = new RJ_AftershockModel_Generic (magMain, rjParam);

			RJ_AftershockModel_Bayesian bayModel = new RJ_AftershockModel_Bayesian(genModel, seqModel);

			// Make the summaries

			RJ_Summary_Generic genSummary = new RJ_Summary_Generic (genModel);

			RJ_Summary_SequenceSpecific seqSummary = new RJ_Summary_SequenceSpecific (seqModel);

			RJ_Summary_Bayesian baySummary = new RJ_Summary_Bayesian (bayModel);

			// Marshal the summaries

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			store.marshalMapBegin (null);

			RJ_Summary_Generic.marshal_poly (store, "Generic", genSummary);

			RJ_Summary_SequenceSpecific.marshal_poly (store, "SeqSpec", seqSummary);

			RJ_Summary_Bayesian.marshal_poly (store, "Bayesian", baySummary);

			// Prepare to unmarshal

			store.marshalMapEnd ();
			store.check_write_complete ();
			String json_string = store.get_json_string();

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			retrieve.unmarshalMapBegin (null);

			// Write the summaries

			System.out.println ("Generic, Direct:\n" + genSummary.toString());

			System.out.println ("Generic, Marshaled:\n" + RJ_Summary_Generic.unmarshal_poly(retrieve, "Generic").toString());

			System.out.println ("Sequence-Specific, Direct:\n" + seqSummary.toString());

			System.out.println ("Sequence-Specific, Marshaled:\n" + RJ_Summary_SequenceSpecific.unmarshal_poly(retrieve, "SeqSpec").toString());

			System.out.println ("Bayesian, Direct:\n" + baySummary.toString());

			System.out.println ("Bayesian, Marshaled:\n" + RJ_Summary_Bayesian.unmarshal_poly(retrieve, "Bayesian").toString());

			retrieve.unmarshalMapEnd ();
			retrieve.check_read_complete ();

			// Display the JSON string

			System.out.println (json_string);

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("RJ_AftershockModel_SequenceSpecific : Unrecognized subcommand : " + args[0]);
		return;

	}

}
