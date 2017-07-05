package scratch.peter.nshmp;

import static org.opensha.nshmp2.util.GaussTruncation.ONE_SIDED;

import java.util.concurrent.TimeUnit;

import org.opensha.commons.calc.GaussianDistCalc;
import org.opensha.commons.exceptions.IMRException;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.param.Parameter;
import org.opensha.nshmp2.util.Utils;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;

import com.google.common.base.Stopwatch;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class PE_Tester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//   
		double mean = -0.6133377010039655;
		double std = 0.564;
		double iml = 0.10754923586540488;
		System.out.println(getExceedProbability(iml, mean, std, false, 0.0));
		System.out.println(getExceedProbability(mean, std, iml, 1));
		
		Stopwatch sw = Stopwatch.createStarted();
		for (int i=0; i<1000000; i++) {
			double r = (Math.random() - 0.5) * 0.1;
			getExceedProbability(iml, mean+r, std, false, 0.0);
		}
		sw.stop();
		System.out.println("nshmp: " + sw.elapsed(TimeUnit.MILLISECONDS));
		sw.reset().start();
		for (int i=0; i<1000000; i++) {
			double r = (Math.random() - 0.5) * 0.1;
			getExceedProbability(mean+r, std, iml, 1);
		}
		sw.stop();
		System.out.println("  sha: " + sw.elapsed(TimeUnit.MILLISECONDS));

	}
	
	public static double getExceedProbability(double iml,
			double mean, double sigma, boolean clamp, double clampVal) {

		double clip = mean + 3 * sigma;
		if (clamp) {
			double clip3s = Math.exp(clip);
			double clipPer = clampVal;
			if (clipPer < clip3s && clipPer > 0) clip = Math.log(clipPer);
		}
		double Pclip = Utils.gaussProbExceed(mean, sigma, clip);
		
		return Utils.gaussProbExceed(mean, sigma, iml, Pclip, ONE_SIDED);
	}

	public static double getExceedProbability(double mean, double stdDev, double iml, int truncType) throws
	ParameterException, IMRException {

		double truncLevel = 3.0;
		if (stdDev != 0) {
			double stRndVar = (iml - mean) / stdDev;
			// compute exceedance probability based on truncation type
			if (truncType == 0) {
				return GaussianDistCalc.getExceedProb(stRndVar);
			}
			else {
				double numSig = truncLevel;
				if (truncType == 1) {
					return GaussianDistCalc.getExceedProb(stRndVar, 1, numSig);
				}
				else {
					return GaussianDistCalc.getExceedProb(stRndVar, 2, numSig);
				}
			}
		}
		else {
			if (iml > mean) {
				return 0;
			}
			else {
				return 1;
			}
		}
	}


}
