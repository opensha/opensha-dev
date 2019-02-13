package scratch.kevin.util;

import com.google.common.base.Preconditions;

public class ReturnPeriodUtils {
	
	/**
	 * Uses equivalence p1Star/t1 = p2Star/t2 where pStar = -Ln(1-p)
	 * @param prob1
	 * @param duration1
	 * @param duration2
	 * @return
	 */
	public static double calcExceedanceProb(double prob1, double duration1, double duration2) {
		Preconditions.checkArgument(prob1 >= 0 && prob1 <= 1);
		Preconditions.checkArgument(duration1 >= 0 && Double.isFinite(duration1));
		Preconditions.checkArgument(duration2 >= 0 && Double.isFinite(duration2));
		double p1star = calcProbStar(prob1);
//		System.out.println("p1star="+(float)p1star);
		double p2star = p1star*duration2/duration1;
//		System.out.println("p2star="+(float)p2star);
		return calcProbFromPorbStar(p2star);
	}
	
	private static double calcProbStar(double prob) {
		// no input validation, assumed done externally
		return -Math.log(1d - prob);
	}
	
	private static double calcProbFromPorbStar(double probStar) {
		return 1d - Math.exp(-probStar);
	}
	
	public static double calcDurationWithExceedanceProb(double exceedProb, double referenceProb, double referenceDuration) {
		double targetProbStar = calcProbStar(exceedProb);
//		targetProbStar = 1;
		double referenceProbStar = calcProbStar(referenceProb);
		
		return referenceDuration * targetProbStar / referenceProbStar;
	}
	
	public static double calcDurationWithExceedanceProb(double exceedProb, double referenceReturnPeriod) {
		double targetProbStar = calcProbStar(exceedProb);
		
		return referenceReturnPeriod * targetProbStar;
	}
	
	public static double calcReturnPeriod(double exceedProb, double duration) {
		return duration / calcProbStar(exceedProb);
	}
	
	public static void main(String[] args) {
		System.out.println(calcReturnPeriod(0.5, 30));
		System.out.println(calcExceedanceProb(0.5, 30, 1d));
//		System.out.println(calcExceedanceProb(0.02, 50, 1d));
		System.exit(0);
//		double r1 = 0.02;
//		double t1 = 50;
//		
//		double[] rps = { 1d, 50d, 500d, 1547.0297, 1733, 2474, 2500 };
//		
//		for (double t2 : rps) {
//			double r2 = calcExceedanceProb(r1, t1, t2);
//			double r2Star = calcProbFromPorbStar(r1)*t2/t1;
//			
//			System.out.println("T2="+(float)t2+"\tR2*="+(float)r2Star+"\tR2="+(float)r2);
//		}
		
		System.out.println("Return periods for probability levels");
		for (double p : new double[] {0.2, 0.1, 0.05, 0.02, 0.01}) {
			System.out.println("\t"+(float)(p*100d)+"% in 50:\t"+(float)calcReturnPeriod(p, 50d));
			System.out.println("\t\tVerificiation: "
					+(float)calcExceedanceProb(0.5, calcDurationWithExceedanceProb(0.5, calcReturnPeriod(p, 50d)), 50d));
			System.out.println("\t\tProb for 1yr curves: "+(float)calcExceedanceProb(p, 50d, 1d));
		}
		System.out.println();
		System.out.println("duration with p=0.5 for 2% in 50");
		System.out.println(calcDurationWithExceedanceProb(0.5, 0.02, 50));
		System.out.println("duration with p=0.5 for 2474.9yr RP");
		System.out.println(calcDurationWithExceedanceProb(0.5, calcReturnPeriod(0.02, 50d)));
		System.out.println("duration with p=0.5 for 2500yr RP");
		System.out.println(calcDurationWithExceedanceProb(0.5, 2500));
		System.out.println();
		System.out.println("duration with p=0.5 for 10% in 50");
		System.out.println(calcDurationWithExceedanceProb(0.5, 0.10, 50));
		System.out.println("duration with p=0.5 for 474.6yr RP");
		System.out.println(calcDurationWithExceedanceProb(0.5, calcReturnPeriod(0.1, 50d)));
		System.out.println("duration with p=0.5 for 500yr RP");
		System.out.println(calcDurationWithExceedanceProb(0.5, 500));
	}

}
