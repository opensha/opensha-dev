package scratch.ned.nshm23;

import java.awt.Color;
import java.util.ArrayList;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

public class PoissonRateFromNinT_Calc {
	
	final static double[] millerRice83_5pt_cumVals = {0.034893, 0.211702, 0.5, 0.788298, 0.965107};
	final static double[] millerRice83_5pt_wts = {0.101, 0.244, 0.31, 0.244, 0.101};
	
	
	public static GammaDistribution getGammaDistribution(int n, double t) {
		double alpha = n+1;  // also called shape paramter
		double beta = 1/t; 
		return new GammaDistribution(alpha,beta);
	}
	
	public static double getRateForCumulativeProb(int n, double t, double prob) {
		GammaDistribution gd = getGammaDistribution(n,t);
//		System.out.println(gd.inverseCumulativeProbability(prob)+"\t"+(float)gd.cumulativeProbability(gd.inverseCumulativeProbability(prob)));
		return gd.inverseCumulativeProbability(prob);
	}
	
	public static void getCumRateDistribution(int n, double t, XY_DataSet distFunc) {
		GammaDistribution gd = getGammaDistribution(n,t);
		for(int i=0;i<distFunc.size();i++) {
			distFunc.set(i,gd.cumulativeProbability(distFunc.getX(i)));
		}
	}
	
	public static void getRateDistribution(int n, double t, XY_DataSet distFunc) {
		GammaDistribution gd = getGammaDistribution(n,t);
		for(int i=0;i<distFunc.size();i++) {
			distFunc.set(i,gd.density(distFunc.getX(i)));
		}
	}
	
	/**
	 * This gives the five-point rates for a range of T values (uniform distribution assumed between T1 and T2)
	 * @param n		N - number of events
	 * @param t1	T1 - lower bound on time period
	 * @param t2	T2 - upper bound on time period
	 */
	public static double[] getFivePointRates(int n, double t1, double t2) {
		if(t1>t2) throw new RuntimeException("t1 must be less that t2");
		
		double[] rateArray = new double[millerRice83_5pt_cumVals.length];
		
		// find and adequate range of rates for the function
		double minRate = getRateForCumulativeProb(n, t2, 0.01); // the 1% rate for the lowest possible rate branch
		double maxRate = getRateForCumulativeProb(n, t1, 0.99); // the 99% rate for the lowest possible rate branch
		int num=1000;
		
		// to match Allison Excel spreadsheet
//		minRate = 1e-5;
//		maxRate = 6.29E-03;
//		num = (int)Math.round(maxRate/minRate);
		
		// mkae the cumulative distribution function
		EvenlyDiscretizedFunc cumDist = new EvenlyDiscretizedFunc(minRate, maxRate, num);
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			EvenlyDiscretizedFunc tempCumDist = new EvenlyDiscretizedFunc(minRate, maxRate, num);
			double t = t1 + (t2-t1)*millerRice83_5pt_cumVals[i];
			getCumRateDistribution(n, t, tempCumDist);
//			if(i==0)System.out.println(tempCumDist);
//			System.out.println("T:\t"+t+"\t"+millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]);
			for(int j=0;j<tempCumDist.size();j++) {
				cumDist.add(j, millerRice83_5pt_wts[i]*tempCumDist.getY(j));
			}
		}
		
//		System.out.println(cumDist);
		
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			rateArray[i] = cumDist.getFirstInterpolatedX(millerRice83_5pt_cumVals[i]);
//			double rate = cumDist.getClosestXtoY(millerRice83_5pt_cumVals[i]);
		}
		return rateArray;

	}
	
	/**
	 * This one samples the uniform T distribution densely, and 
	 * gives a result very similar to other implementation
	 * @param n
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static double[] getFivePointRatesAlt(int n, double t1, double t2) {
		if(t1>t2) throw new RuntimeException("t1 must be less that t2");
		
		double[] rateArray = new double[millerRice83_5pt_cumVals.length];
		
		// find and adequate range of rates for the function
		double minRate = getRateForCumulativeProb(n, t2, 0.01); // the 1% rate for the lowest possible rate branch
		double maxRate = getRateForCumulativeProb(n, t1, 0.99); // the 99% rate for the lowest possible rate branch
		int num=1000;
		
		// to match Allison Excel spreadsheet
//		minRate = 1e-5;
//		maxRate = 6.29E-03;
//		num = (int)Math.round(maxRate/minRate);
		
		// make the cumulative distribution function
		EvenlyDiscretizedFunc cumDist = new EvenlyDiscretizedFunc(minRate, maxRate, num);
		
		int numT=1000;
		double deltaT = (t2-t1)/numT;
		EvenlyDiscretizedFunc tFunc= new EvenlyDiscretizedFunc(t1+deltaT/2d, t2-deltaT/2d, numT);
		
		for(int i=0;i<tFunc.size();i++) {
			EvenlyDiscretizedFunc tempCumDist = new EvenlyDiscretizedFunc(minRate, maxRate, num);
			double t = tFunc.getX(i);
			getCumRateDistribution(n, t, tempCumDist);
//			if(i==0)System.out.println(tempCumDist);
//			System.out.println("T:\t"+t+"\t"+millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]);
			for(int j=0;j<tempCumDist.size();j++) {
				cumDist.add(j, tempCumDist.getY(j)/numT);
			}
		}
		
//		System.out.println(cumDist);
		
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			rateArray[i] = cumDist.getFirstInterpolatedX(millerRice83_5pt_cumVals[i]);
//			double rate = cumDist.getClosestXtoY(millerRice83_5pt_cumVals[i]);
		}
		return rateArray;

	}

	
	
	public static double[] getFivePointRates(int n, double t) {
		double[] rateArray = new double[millerRice83_5pt_cumVals.length];
		GammaDistribution gd = getGammaDistribution(n,t);
		for(int i=0;i<rateArray.length;i++)
			rateArray[i] = gd.inverseCumulativeProbability(millerRice83_5pt_cumVals[i]);
		return rateArray;
	}

	
		
	/**
	 * This demonstrates that the Poisson rate likelihood is equivalent to the Gamma distribution
	 * with alpha=N+1 and beta=1/T
	 * @param n
	 * @param t
	 * @param distFunc
	 */
	public static void testDistribution(int n, double t, XY_DataSet distFunc) {
		GammaDistribution gd = getGammaDistribution(n,t);
		double nFactorial = 1;
		for(int i =n; i>0; i--) nFactorial *= i;
		for(int i=0;i<distFunc.size();i++) {
			double rate = distFunc.getX(i);
			double val1 = gd.density(rate);
			double val2 = t*Math.pow(rate*t,n)*Math.exp(-rate*t)/nFactorial;
			distFunc.set(i,val1/val2);
		}
	}

	
	public static double solveForT(int n, double rate, double cumProb) {
		double maxT = 10*n/rate;
		double minT = 0.1*n/rate;
		int num = 10000;
		EvenlyDiscretizedFunc rateDist = new EvenlyDiscretizedFunc(minT, maxT, num);
		for(int i=0;i <rateDist.size();i++) {
			double t = rateDist.getX(i);
			GammaDistribution gd = getGammaDistribution(n,t);
			rateDist.set(i,gd.inverseCumulativeProbability(cumProb));
		}
		return rateDist.getFirstInterpolatedX(rate);
	}
	
	
	public static double[] writeFivePointRates(String name, int n, double t1, double t2, double[] testValsArray) {
		double[] rateArray;
		if(!Double.isNaN(t2))
			rateArray = getFivePointRates(n, t1, t2);
		else
			rateArray = getFivePointRates(n, t1);
		String result = "";
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			result += String.format("%.3f",millerRice83_5pt_cumVals[i])+"\t"+
					  String.format("%.3f",millerRice83_5pt_wts[i])+"\t"+String.format("%6.3e",rateArray[i]);
			if(testValsArray != null)
				result += "\t"+String.format("%.2f",(rateArray[i]/testValsArray[i]));
			result += "\n";
		}
		System.out.println("\n"+name+":\n"+"prob\tweight\trate\t\ttestRatio\n"+result);
		return rateArray;
	}
	
	

	public static void main(String[] args) {
		
		// Notes:
		// 1) CEUS SSCn (& some of Allison's tables) list  values in reverse order of what's here
		// 2) getFivePointRatesAlt() produces almost identical results (diff approximation of the uniform distribution for T)
		
		// Following test values come from Allison's Electronic Supplement; first one off because her rate range was not wide enough
		writeFivePointRates("Saline River 1 in 322-5100", 1, 322, 5100, new double[] {1.00E-04,3.10E-04,6.80E-04,1.52E-03,4.84E-03});
		
		// Following test values come from Allison's Electronic Supplement; Values differ a bit from Allison's because her's are accidently offset (first two values are the same)
		writeFivePointRates("Central Virginia 1 in 1800-2800", 1, 1800, 2800, new double[] {1.4E-04,3.9E-04,7.5E-04,1.3E-03,2.3E-03});

		// Following test values come from Table H-5.4-2 from CEUS SSCn (2012); another branch (20% wt) has 3 RIs
		writeFivePointRates("Meers 1 in 2153-2968", 1, 2153, 2968, new double[] {1.2E-04,3.4E-04,6.7E-04,1.2E-03,2.1E-03});

		// Following test values come from Table H-5.9-2 from CEUS SSCn (2012)
		writeFivePointRates("Wabash Valley 1 in 11k-13k", 1, 11e3, 13e3, new double[] {2.4E-05,7.2E-05,1.4E-04,2.5E-04,4.4E-04});

		// Following test values come from Table H-5.7-2 from CEUS SSCn (2012) 
		writeFivePointRates("Marianna: 2 in 9.6k-10.2k",2, 9.6e3, 10.2e3, new double[] {7.2E-05,1.6E-04,2.7E-04,4.2E-04,6.9E-04});
		
		// Following test values come from Table H-5.7-3 from CEUS SSCn (2012) 
		writeFivePointRates("Marianna: 3 in 9.6k-10.2k",3, 9.6e3, 10.2e3, new double[] {1.2E-04,2.4E-04,3.7E-04,5.5E-04,8.4E-04});
		
		// Following test values come from Table H-5.1-2 from CEUS SSCn (2012); Table 6.1-1 says N=1 T=348 but other table says "Recurrence Intervals"
		writeFivePointRates("Charlevoix 1 in 348", 1, 348, Double.NaN, new double[] {7.7E-04,2.2E-03,4.2E-03,6.7E-03,9.3E-03});
		//this one not matching
		double t = solveForT(1, 4.2e-3, 0.5); System.out.println("T solved for: "+t+"??");
		writeFivePointRates("Charlevoix attemped fix of 1 in "+Float.toString((float)t), 1, t, Double.NaN, new double[] {7.7E-04,2.2E-03,4.2E-03,6.7E-03,9.3E-03});
	
		// Following test values come from Table H-5.1-3 from CEUS SSCn (2012)
		writeFivePointRates("Charlevoix 3 in 6k-7k", 3, 6e3, 7e3, new double[] {1.9E-04,3.7E-04,5.7E-04,8.4E-04,1.3E-03});

		// Following test values come from Table H-5.1-4 from CEUS SSCn (2012)
		writeFivePointRates("Charlevoix 4 in 9.5k-10.2k", 4, 9.5e3, 10.2e3, new double[] {1.8E-04,3.2E-04,4.7E-04,6.7E-04,9.8E-04});

		// Following test values come from first (& 2nd) column of Table 6.1.2-4 from CEUS SSCn (2012); Range of T is the CON bounds 95% bounds for for Earthquake C in Table 6.1-1
		writeFivePointRates("Charleston 4 in (1569+1854)/2", 4, (1569+1854)/2, Double.NaN, new double[] {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03});
//		writeFivePointRates("Charleston 4 in 1569", 4, 1569, Double.NaN, new double[] {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03});
//		writeFivePointRates("Charleston 4 in 1854", 4, 1854, Double.NaN, new double[] {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03});
		t = solveForT(4, 2.1E-03, 0.5); System.out.println("T solved for: "+t+"??");
		writeFivePointRates("Charleston attemped fix of 4 in "+Float.toString((float)t), 4, t, Double.NaN, new double[] {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03});
		writeFivePointRates("Charleston attempted fix of 3 in (1569+1854)/2", 3, (1569+1854)/2, Double.NaN, new double[] {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03});

		
		
//		writeFivePointRates("", , new double[] {});

		
		// Cannot reproduce what's below:

		
//		// Charleston
//		double t = solveForT(6, 1.1e-3, 0.5);
//		System.out.println("\nCharleston T: "+t);
//		double[] rateArray = getFivePointRates(6, t);	// 6 events since event E
//		for(int i=0;i<millerRice83_5pt_cumVals.length;i++)
//			System.out.println(millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]+"\t"+(float)rateArray[i]);

//		double t = solveForT(4, 2.1e-3, 0.5);
//		System.out.println("\nCharleston T: "+t);
//		double[] rateArray = getFivePointRates(4, t);	// 4 events since event C
//		for(int i=0;i<millerRice83_5pt_cumVals.length;i++)
//			System.out.println(millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]+"\t"+(float)rateArray[i]);
		
//		System.out.println("\nCharleston: 1569");
//		double[] rateArray = getFivePointRates(4, 1569);	// 4 events since event C
//		for(int i=0;i<millerRice83_5pt_cumVals.length;i++)
//			System.out.println(millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]+"\t"+(float)rateArray[i]);
//		System.out.println("\nCharleston: 1854");
//		rateArray = getFivePointRates(4, 1854);	// 4 events since event C
//		for(int i=0;i<millerRice83_5pt_cumVals.length;i++)
//			System.out.println(millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]+"\t"+(float)rateArray[i]);
//		System.out.println("\nCharleston: 2000");
//		rateArray = getFivePointRates(4, 2000);	// 4 events since event C
//		for(int i=0;i<millerRice83_5pt_cumVals.length;i++)
//			System.out.println(millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]+"\t"+(float)rateArray[i]);

		
		
		
//		int n=2;
//		double t=2000;
//		
//		n=1;
//		t=488.718754;
//		
//		EvenlyDiscretizedFunc func1 = new EvenlyDiscretizedFunc(1e-4, 500, 1e-5);
//		getCumRateDistribution(n, t, func1);
//		for(int i=0;i<func1.size();i++)
//			System.out.println((float)func1.getX(i)+"\t"+(float)func1.getY(i));
//		
//		EvenlyDiscretizedFunc func2 = new EvenlyDiscretizedFunc(1e-4, 500, 1e-5);
//		getRateDistribution(n, t, func2);
//		func2.scale(1.0/func2.getMaxY());
//
//		
//    	ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
//    	funcs.add(func1);
//    	funcs.add(func2);
//    	
//    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
//    	
//		GraphWindow graph = new GraphWindow(funcs, "",plotChars); 
//		graph.setX_AxisLabel("rate");
//		graph.setY_AxisLabel("likelihood");
//		graph.setXLog(true);
////		graph.setX_AxisRange(50, 2e4);
////		graph.setY_AxisRange(5, 9);
//		graph.setPlotLabelFontSize(18);
//		graph.setAxisLabelFontSize(18);
//		graph.setTickLabelFontSize(16);


		
		
//		// test equivalence
//		EvenlyDiscretizedFunc func1 = new EvenlyDiscretizedFunc(1e-4, 500, 1e-5);
//		testDistribution(n, t, func1);
//		for(int i=0;i<func1.size();i++)
//			System.out.println((float)func1.getX(i)+"\t"+(float)func1.getY(i));

		

		
	}

}
