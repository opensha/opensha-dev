package scratch.ned.nshm23;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

public class PoissonRateFromNinT_Calc {
	
	final static double[] millerRice83_5pt_cumVals = {0.034893, 0.211702, 0.5, 0.788298, 0.965107};
	final static double[] millerRice83_5pt_wts = {0.101, 0.244, 0.31, 0.244, 0.101};
	
	// these are used to define the incremental x-axis rate values for PDFs and CDFs.
	final static double minProb = 0.001;
	final static double maxProb = 0.999;
	final static int numX = 10000;
	
	

	
	/**
	 * This provides and instance of a GammaDistribution for the given N and T
	 * @param n
	 * @param t
	 * @return
	 */
	private static GammaDistribution getGammaDistribution(int n, double t) {
		double alpha = n+1;  // also called shape paramter
		double beta = 1/t; 
		return new GammaDistribution(alpha,beta);
	}
	
	/**
	 * This fills in the incremental probability (or density, on y-axis ) for each rate on the
	 * x-axis for the provided distFunc.
	 * @param n
	 * @param t
	 * @param minRate - min x-axis value for distribution function
	 * @param maxRate - max x-axis value for distribution function
	 * @param numRate - num points on x-axis
	 */
	public static EvenlyDiscretizedFunc getRateDistribution(int n, double t, double minRate, double maxRate, int numRate) {
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(minRate,maxRate,numRate);
		GammaDistribution gd = getGammaDistribution(n,t);
		for(int i=0;i<distFunc.size();i++) {
			distFunc.set(i,gd.density(distFunc.getX(i)));
		}
		return distFunc;
	}
	
	/**
	 * this provides the rate for the given cumulative probability
	 * @param n
	 * @param t
	 * @param cumProb
	 * @return
	 */
	public static double getRateForCumulativeProb(int n, double t, double cumProb) {
		GammaDistribution gd = getGammaDistribution(n,t);
		return gd.inverseCumulativeProbability(cumProb);
	}
	
	
	/**
	 * This fills in the cumulative probability (y-axis value) for each rate on the
	 * x-axis for the provided distFunc.
	 * @param n
	 * @param t
	 * @param minRate - min x-axis value for distribution function
	 * @param maxRate - max x-axis value for distribution function
	 * @param numRate - num points on x-axis
	 */
	public static EvenlyDiscretizedFunc getCumRateDistribution(int n, double t, double minRate, double maxRate, int numRate) {
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(minRate,maxRate,numRate);
		GammaDistribution gd = getGammaDistribution(n,t);
		for(int i=0;i<distFunc.size();i++) {
			distFunc.set(i,gd.cumulativeProbability(distFunc.getX(i)));
		}
		return distFunc;
	}
	
	/**
	 * This fills in the cumulative probability (y-axis value) for each rate on the
	 * x-axis for the provided distFunc, assuming a uniform probability of T being
	 * between T1 and T2.  This uses a 5-point approximation for this uniform T distribution.
	 * @param n
	 * @param t1
	 * @param t2
	 * @param minRate - min x-axis value for distribution function
	 * @param maxRate - max x-axis value for distribution function
	 * @param numRate - num points on x-axis
	 */
	public static EvenlyDiscretizedFunc getCumRateDistribution(int n, double t1, double t2, double minRate, double maxRate, int numRate) {
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(minRate,maxRate,numRate);
		
		// make the cumulative distribution function sampling the range of T with 5 points
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			double t = t1 + (t2-t1)*millerRice83_5pt_cumVals[i];
			EvenlyDiscretizedFunc tempCumDist = getCumRateDistribution(n, t, minRate, maxRate, numRate);
//			if(i==0)System.out.println(tempCumDist);
//			System.out.println("T:\t"+t+"\t"+millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]);
			for(int j=0;j<tempCumDist.size();j++) {
				distFunc.add(j, millerRice83_5pt_wts[i]*tempCumDist.getY(j));
			}
		}
		return distFunc;
	}
	
	/**
	 * This gets the cumulative distribution for the given NinT_Data object.  This does not apply the 
	 * weight therein.
	 * @param nInT_Data
	 * @param minRate - min x-axis value for distribution function
	 * @param maxRate - max x-axis value for distribution function
	 * @param numRate - num points on x-axis
	 */
	public static EvenlyDiscretizedFunc getCumRateDistribution(NinT_Data nInT_Data, double minRate, double maxRate, int numRate) {
		if(Double.isNaN(nInT_Data.t2))
			return getCumRateDistribution(nInT_Data.n, nInT_Data.t, minRate, maxRate, numRate);
		else
			return getCumRateDistribution(nInT_Data.n, nInT_Data.t, nInT_Data.t2, minRate, maxRate, numRate);
	}
	
	
	
	/**
	 * This gets the cumulative distribution for the given list of NinT_Data objects by stacking
	 * the associated CDFs with the weight value in each NinT_Data object.  This method solves for
	 * the minRate, maxRate, and numRate required in the other like-named method.
	 * @param List<nInT_Data>
	 */
	public static EvenlyDiscretizedFunc getCumRateDistribution(List<NinT_Data> nInT_DataList) {

		// find min and max rate for the distribution
		double minRate = Double.POSITIVE_INFINITY;
		double maxRate = Double.NEGATIVE_INFINITY;
		for(NinT_Data data:nInT_DataList) {
			double testMinRate = getRateForCumulativeProb(data.n, data.t, minProb);
			double testMaxRate = getRateForCumulativeProb(data.n, data.t, maxProb);
			if(!Double.isNaN(data.t2)) {
				testMinRate = getRateForCumulativeProb(data.n, data.t2, minProb); // the 1% rate for the lowest possible rate branch
			}
			if(minRate>testMinRate)
				minRate = testMinRate;
			if(maxRate<testMaxRate)
				maxRate = testMaxRate;
		}
		return getCumRateDistribution(nInT_DataList, minRate, maxRate, numX);	
	}

	
	/**
	 * This gets the cumulative distribution for the given list of NinT_Data objects by stacking
	 * the associated CDFs with the weight value in each NinT_Data object.
	 * @param List<nInT_Data>
	 * @param minRate - min x-axis value for distribution function
	 * @param maxRate - max x-axis value for distribution function
	 * @param numRate - num points on x-axis
	 */
	public static EvenlyDiscretizedFunc getCumRateDistribution(List<NinT_Data> nInT_DataList, double minRate, double maxRate, int numRate) {
	
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(minRate,maxRate,numRate);
		double weightTest = 0;
		
//ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		for(NinT_Data nInT_Data : nInT_DataList) {
			EvenlyDiscretizedFunc tempFunc = null;
			if(Double.isNaN(nInT_Data.t2))
				tempFunc = getCumRateDistribution(nInT_Data.n, nInT_Data.t, minRate, maxRate, numRate);
			else
				tempFunc = getCumRateDistribution(nInT_Data.n, nInT_Data.t, nInT_Data.t2, minRate, maxRate, numRate);
			
			
//// test means
//			double mean1 = getMeanFromFivePointRates(getFivePointRatesFromDistribution(tempFunc));
//			HistogramFunction pdf = new HistogramFunction(tempFunc.getMinX(),tempFunc.getMaxX(),tempFunc.size());
//			for(int i=tempFunc.size()-1; i>0; i--)
//				pdf.set(i,tempFunc.getY(i)-tempFunc.getY(i-1));
//			double mean2 = pdf.computeMean();
//			funcs.add(pdf);
//			System.out.println("\tMeanTest:\t5pt = "+(float)mean1+"\tpdf = "+(float)mean2+"\tratio = "+(float)(mean1/mean2));			
//// end test	
			
//			tempFunc.scale(nInT_Data.weight); // scale by the weight
			weightTest+=nInT_Data.weight;
			for(int i=0;i<distFunc.size();i++) {
				distFunc.add(i, nInT_Data.weight*tempFunc.getY(i));
			}
		}
		if(Math.abs(weightTest-1.0) > 0.001)
			throw new RuntimeException("weights do not add to 1.0");
		
//// test means
//		double mean1 = getMeanFromFivePointRates(getFivePointRatesFromDistribution(distFunc));
//		HistogramFunction pdf = new HistogramFunction(distFunc.getMinX(),distFunc.getMaxX(),distFunc.size());
//		for(int i=distFunc.size()-1; i>0; i--)
//			pdf.set(i,distFunc.getY(i)-distFunc.getY(i-1));
//		double mean2 = pdf.computeMean();
//		funcs.add(pdf);
//		System.out.println("\tMeanTestList:\t5pt = "+(float)mean1+"\tpdf = "+(float)mean2+"\tratio = "+(float)(mean1/mean2));			
//
//		plotRateDistributions(funcs);
//// end test
		
		
		return distFunc;
		
	}



	
	/**
	 * This fills in the cumulative probability (y-axis value) for each rate on the
	 * x-axis for the provided distFunc, assuming a uniform probability of T being
	 * between T1 and T2.  This samples the uniform T distribution more densely.
	 * @param n
	 * @param t1
	 * @param t2
	 * @param minRate - min x-axis value for distribution function
	 * @param maxRate - max x-axis value for distribution function
	 * @param numRate - num points on x-axis
	 */
	public static EvenlyDiscretizedFunc getCumRateDistributionAlt(int n, double t1, double t2, double minRate, double maxRate, int numRate) {
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(minRate,maxRate,numRate);
		
		// discrtize the unfirm distribution between T1 and T2
		int numT=1000;
		double deltaT = (t2-t1)/numT;
		EvenlyDiscretizedFunc tFunc= new EvenlyDiscretizedFunc(t1+deltaT/2d, t2-deltaT/2d, numT);
		
		for(int i=0;i<tFunc.size();i++) {
			double t = tFunc.getX(i);
			EvenlyDiscretizedFunc tempCumDist = getCumRateDistribution(n, t, minRate,maxRate,numRate);
//			if(i==0)System.out.println(tempCumDist);
//			System.out.println("T:\t"+t+"\t"+millerRice83_5pt_cumVals[i]+"\t"+millerRice83_5pt_wts[i]);
			for(int j=0;j<tempCumDist.size();j++) {
				distFunc.add(j, tempCumDist.getY(j)/numT);
			}
		}
		return distFunc;
	}

	/**
	 * This gives the millerRice83 five-point values from the given cumulative distribution
	 * @param cumDist
	 * @return
	 */
	public static double[] getFivePointRatesFromDistribution(EvenlyDiscretizedFunc cumDist) {
		double[] rateArray = new double[millerRice83_5pt_cumVals.length];
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			rateArray[i] = cumDist.getFirstInterpolatedX(millerRice83_5pt_cumVals[i]);
		}
		
		return rateArray;
	}
	
	
	/**
	 * This gives the five-point rates for a range of T values (uniform distribution assumed between T1 and T2)
	 * @param n		N - number of events
	 * @param t1	T1 - lower bound on time period
	 * @param t2	T2 - upper bound on time period
	 */
	public static double[] getFivePointRates(int n, double t1, double t2) {
		if(t1>t2) throw new RuntimeException("t1 must be less that t2");
		
		
		// find and adequate range of rates for the function
		double minRate = getRateForCumulativeProb(n, t2, minProb); // the 1% rate for the lowest possible rate branch
		double maxRate = getRateForCumulativeProb(n, t1, maxProb); // the 99% rate for the lowest possible rate branch
		
		// to match Allison Excel spreadsheet
//		minRate = 1e-5;
//		maxRate = 6.29E-03;
//		num = (int)Math.round(maxRate/minRate);

		// make the cumulative distribution function sampling the range of T with 5 points
		EvenlyDiscretizedFunc cumDist = getCumRateDistribution(n, t1, t2, minRate, maxRate, numX);
		
//		System.out.println(cumDist);
		return getFivePointRatesFromDistribution(cumDist);
	}
	
	
	/**
	 * This one samples the uniform T distribution densely, and 
	 * gives a result very similar to other implementation
	 * @param n
	 * @param t1
	 * @param t2
	 * @return double[]
	 */
	public static double[] getFivePointRatesAlt(int n, double t1, double t2) {
		if(t1>t2) throw new RuntimeException("t1 must be less that t2");
		
		// find and adequate range of rates for the function
		double minRate = getRateForCumulativeProb(n, t2, minProb); // the 1% rate for the lowest possible rate branch
		double maxRate = getRateForCumulativeProb(n, t1, maxProb); // the 99% rate for the lowest possible rate branch
		
		// to match Allison Excel spreadsheet
//		minRate = 1e-5;
//		maxRate = 6.29E-03;
//		num = (int)Math.round(maxRate/minRate);
		
		// make the cumulative distribution function
		EvenlyDiscretizedFunc cumDist = getCumRateDistributionAlt(n, t1, t2, minRate, maxRate, numX);
		
//		System.out.println(cumDist);
		return getFivePointRatesFromDistribution(cumDist);
	}

	/**
	 * This gives an array of rates for the Miller and Rice (1985) five-point distribution
	 * @param n
	 * @param t
	 * @return double[]
	 */
	public static double[] getFivePointRates(int n, double t) {
		double[] rateArray = new double[millerRice83_5pt_cumVals.length];
		GammaDistribution gd = getGammaDistribution(n,t);
		for(int i=0;i<rateArray.length;i++)
			rateArray[i] = gd.inverseCumulativeProbability(millerRice83_5pt_cumVals[i]);
		return rateArray;
	}

	
	/**
	 * This gives an array of rates for the Miller and Rice (1985) five-point distribution
	 * @param NinT_Data
	 * @return double[]
	 */
	public static double[] getFivePointRates(NinT_Data nInT_Data) {
		if(Double.isNaN(nInT_Data.t2))
			return getFivePointRates(nInT_Data.n, nInT_Data.t);
		else
			return getFivePointRates(nInT_Data.n, nInT_Data.t, nInT_Data.t2);
	}

	
	/**
	 * This gives an array of rates for the Miller and Rice (1985) five-point distribution
	 * @param List<NinT_Data>
	 * @return double[]
	 */
	public static double[] getFivePointRates(List<NinT_Data> nInT_DataList) {
		
		EvenlyDiscretizedFunc cumDist = getCumRateDistribution(nInT_DataList);
		
		return getFivePointRatesFromDistribution(cumDist);
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

	/**
	 * This solves for T given n, rate, and cumProb for that rate.
	 * @param n
	 * @param rate
	 * @param cumProb
	 * @return
	 */
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
	
	
	/**
	 * this gives the mean rate from a five-point distribution
	 * @param rateArray
	 * @return
	 */
	public static double getMeanFromFivePointRates(double[] rateArray) {
		double mean = 0;
		if(rateArray.length != 5)
			throw new RuntimeException("array must have five elements");
		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			mean += millerRice83_5pt_wts[i]*rateArray[i];
		}
		return mean;
	}
	
	
	public static double getMeanFromCDF(EvenlyDiscretizedFunc cdf) {
		// convert cdf to pdf
		HistogramFunction pdf = new HistogramFunction(cdf.getMinX(),cdf.getMaxX(),cdf.size());
		for(int i=cdf.size()-1; i>0; i--)
			pdf.set(i,cdf.getY(i)-cdf.getY(i-1));
		return pdf.computeMean();
	}

	
	
	/**
	 * This writes out the five-point rates and compares them to the given target testValsArray
	 * @param name
	 * @param n
	 * @param t1
	 * @param t2
	 * @param testValsArray
	 * @return
	 */
	public static double[] writeFivePointRates(String name, int n, double t1, double t2, double[] testValsArray) {
		double[] rateArray;
		if(!Double.isNaN(t2))
			rateArray = getFivePointRates(n, t1, t2);
		else
			rateArray = getFivePointRates(n, t1);
		String result = "";
		double mean = 0;
		double testValsMean = 0;

		for(int i=0;i<millerRice83_5pt_cumVals.length;i++) {
			result += String.format("%.3f",millerRice83_5pt_cumVals[i])+"\t"+
					  String.format("%.3f",millerRice83_5pt_wts[i])+"\t"+String.format("%6.3e",rateArray[i]);
			if(testValsArray != null) {
				result += "\t"+String.format("%.2f",(rateArray[i]/testValsArray[i]));
				testValsMean+=testValsArray[i]*millerRice83_5pt_wts[i];
			}
			result += "\n";
			mean += millerRice83_5pt_wts[i]*rateArray[i];
		}
		String appendString="";
		if(testValsArray != null)
			appendString = "\ttestValMean="+(float)testValsMean+"\tratio="+(float)(mean/testValsMean);
		System.out.println("\n"+name+":\n"+"prob\tweight\trate\t\ttestRatio\n"+result+"mean="+(float)mean+appendString+"\n");
		return rateArray;
	}
	
	
	
	public static void plotRateDistributions(int n, double t) {
		
		double minRate = getRateForCumulativeProb(n, t, minProb);
		double maxRate = getRateForCumulativeProb(n, t, maxProb);
		
		EvenlyDiscretizedFunc func1 = getCumRateDistribution(n, t, minRate,maxRate,1000);
//		for(int i=0;i<func1.size();i++)
//			System.out.println((float)func1.getX(i)+"\t"+(float)func1.getY(i));
		
		EvenlyDiscretizedFunc func2 = getRateDistribution(n, t, minRate,maxRate,1000);
		func2.scale(1.0/func2.getMaxY());

		
    	ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
    	funcs.add(func1);
    	funcs.add(func2);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
    	
		GraphWindow graph = new GraphWindow(funcs, "",plotChars); 
		graph.setX_AxisLabel("rate");
		graph.setY_AxisLabel("likelihood");
//		graph.setXLog(true);
//		graph.setX_AxisRange(50, 2e4);
//		graph.setY_AxisRange(5, 9);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);

	}
	
	
	public static void plotRateDistributions(XY_DataSet func) {
				
    	ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
    	funcs.add(func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
    	
		GraphWindow graph = new GraphWindow(funcs, "",plotChars); 
		graph.setX_AxisLabel("rate");
		graph.setY_AxisLabel("likelihood");
//		graph.setXLog(true);
//		graph.setX_AxisRange(50, 2e4);
//		graph.setY_AxisRange(5, 9);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);

	}

	
	
	public static void plotRateDistributions(ArrayList<XY_DataSet> funcs) {
		
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
    	
		GraphWindow graph = new GraphWindow(funcs, "",plotChars); 
		graph.setX_AxisLabel("rate");
		graph.setY_AxisLabel("likelihood");
//		graph.setXLog(true);
//		graph.setX_AxisRange(50, 2e4);
//		graph.setY_AxisRange(5, 9);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);

	}

	
	public static void test2023_CEUS_results() {
		// Notes:
		// 1) CEUS SSCn (& some of Allison's tables) list  values in reverse order of what's here
		// 2) getFivePointRatesAlt() produces almost identical results (diff approximation of the uniform distribution for T)
		
		// Following test values come from Table H-5.7-2 from CEUS SSCn (2012) 
		writeFivePointRates("Marianna: 2 in 9.6k-10.2k",2, 9.6e3, 10.2e3, new double[] {7.2E-05,1.6E-04,2.7E-04,4.2E-04,6.9E-04});
		// Following test values come from Table H-5.7-3 from CEUS SSCn (2012) 
		writeFivePointRates("Marianna: 3 in 9.6k-10.2k",3, 9.6e3, 10.2e3, new double[] {1.2E-04,2.4E-04,3.7E-04,5.5E-04,8.4E-04});
		// Note that Figure 28 of Petersen et al. (2014) has the N values erroneously increased by 1
		// Above check against mean value in Peter's files


		// Following test values come from Table H-5.9-2 from CEUS SSCn (2012)
		writeFivePointRates("Wabash Valley 1 in 11k-13k", 1, 11e3, 13e3, new double[] {2.4E-05,7.2E-05,1.4E-04,2.5E-04,4.4E-04});
		// the mean here is 1.4% below that in Peter's file

		// Following test values come from Allison's Electronic Supplement; first one off because her rate range was not wide enough
		writeFivePointRates("Saline River 1 in 322-5100", 1, 322, 5100, new double[] {1.00E-04,3.10E-04,6.80E-04,1.52E-03,4.84E-03});
				
		// Following test values come from Allison's Electronic Supplement; Values differ a bit from Allison's because her's are accidently offset (first two values are the same)
//		writeFivePointRates("Central Virginia 1 in 1800-2800", 1, 1800, 2800, new double[] {1.4E-04,3.9E-04,7.5E-04,1.3E-03,2.3E-03});
		// updated, apparently to fix Allison's previous error
		writeFivePointRates("Central Virginia 1 in 1800-2800", 1, 1800, 2800, new double[] {1.3E-04,3.8E-04,7.4E-04,1.3E-03,2.3E-03});
		
		// Following test values come from Table H-5.4-2 from CEUS SSCn (2012); another branch (20% wt) has 3 RIs
		writeFivePointRates("Meers 1 in 2153-2968", 1, 2153, 2968, new double[] {1.2E-04,3.4E-04,6.7E-04,1.2E-03,2.1E-03});


		// Charlevoix:
		// Summing the following means with their branch wts (0.2, 0.6, 0.2), From Figure ?? of Petersen et al (2014), matches the
		// rate in Peter's file.
		
		// Following test values come from Table H-5.1-2 from CEUS SSCn (2012); Table 6.1-1 says N=1 T=348 but other table 
		// says "Recurrence Intervals".  Here, N is the number of intervals (not quakes) and so it is 1.0
		writeFivePointRates("Charlevoix 1 interval in 348", 1, 348, Double.NaN, new double[] {7.7E-04,2.2E-03,4.2E-03,6.7E-03,9.3E-03});
		//this one not matching, the following matches the median target
//		double t = solveForT(1, 4.2e-3, 0.5); System.out.println("T solved for: "+t+"??");
//		writeFivePointRates("Charlevoix attemped fix of 1 in "+Float.toString((float)t), 1, t, Double.NaN, new double[] {7.7E-04,2.2E-03,4.2E-03,6.7E-03,9.3E-03});
		// this one gets the right mean target, which is what we want to match NSHM
		int n=1; // intervals, not quakes
		double t=(n+1)/0.00449067; // the latter is the target
		writeFivePointRates("Charlevoix attemped fix of 1 in "+Float.toString((float)t), n, t, Double.NaN, new double[] {7.7E-04,2.2E-03,4.2E-03,6.7E-03,9.3E-03});

		// Following test values come from Table H-5.1-3 from CEUS SSCn (2012)
		writeFivePointRates("Charlevoix 3 in 6k-7k", 3, 6e3, 7e3, new double[] {1.9E-04,3.7E-04,5.7E-04,8.4E-04,1.3E-03});
		
		// Following test values come from Table H-5.1-4 from CEUS SSCn (2012)
		writeFivePointRates("Charlevoix 4 in 9.5k-10.2k", 4, 9.5e3, 10.2e3, new double[] {1.8E-04,3.2E-04,4.7E-04,6.7E-04,9.8E-04});

		
//		System.exit(0);
	


		
		// For Charleston, it appears they randomly sampled event dates from distributions in Figure 6.1.2-19 and averaged results.  
		// Here I compute effective T from their branch weighted mean for each column in Table 6.1.2-4
		// Note that to match their results, I had to set N as the number of intervals (1 minus num quakes); doing so also leads to 
		// implied T that matches the PDF for oldest event in each case (their figure 6.1.2-19).
		// Wt averaging the 5 cases using branch wts in Figure 6.1.2-1a (0.8,0.04,0.06,0.04,0.06, respectively) and then multiplying by 0.9
		// (effective prob of existance weight in Figure 30 of Petersen et al. (2014)) equals the ave rate implied by Peter's file

		// Post-2000 years Earthquakes 1886, A, B, and C 
		double[] targetBranchRates1 = {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03};
		double meanRate =0;
		for(int i=0; i<5; i++) meanRate += millerRice83_5pt_wts[i]*targetBranchRates1[i];
		n=3; // num intervals, not num events
		t = (n+1)/meanRate;
		writeFivePointRates("Charleston Post-2000 years Earthquakes 1886, A, B, and C; T = "+(float)t, n, t, Double.NaN, targetBranchRates1);

		// Post-5,500 years Earthquakes 1886, A, B, and C; 
		double[] targetBranchRates2 = {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03};
		meanRate =0;
		for(int i=0; i<5; i++) meanRate += millerRice83_5pt_wts[i]*targetBranchRates2[i];
		n=3; // num intervals, not num events
		t = (n+1)/meanRate;
		writeFivePointRates("Charleston Post-5,500 years Earthquakes 1886, A, B, and C; T = "+(float)t, n, t, Double.NaN, targetBranchRates2);

		// Post-5,500 years Earthquakes 1886, A, B, C and D; 
		double[] targetBranchRates3 = {5.0E-04,8.8E-04,1.3E-03,1.9E-03,2.7E-03};
		meanRate =0;
		for(int i=0; i<5; i++) meanRate += millerRice83_5pt_wts[i]*targetBranchRates3[i];
		n=4; // num intervals, not num events
		t = (n+1)/meanRate;
		writeFivePointRates("Charleston Post-5,500 years Earthquakes 1886, A, B, C and D; T = "+(float)t, n, t, Double.NaN, targetBranchRates3);

		
		// Post-5,500 years Earthquakes 1886, A, B, C and E; 
		double[] targetBranchRates4 = {3.4E-04,6.4E-04,9.2E-04,1.3E-03,1.9E-03};
		meanRate =0;
		for(int i=0; i<5; i++) meanRate += millerRice83_5pt_wts[i]*targetBranchRates4[i];
		n=4; // num intervals, not num events
		t = (n+1)/meanRate;
		writeFivePointRates("Charleston Post-5,500 years Earthquakes 1886, A, B, C and E; T = "+(float)t, n, t, Double.NaN, targetBranchRates4);
		
		// Post-5,500 years Earthquakes 1886, A, B, C, D and E; 
		double[] targetBranchRates5 = {4.6E-04,7.8E-04,1.1E-03,1.5E-03,2.2E-03};
		meanRate =0;
		for(int i=0; i<5; i++) meanRate += millerRice83_5pt_wts[i]*targetBranchRates5[i];
		n=5; // num intervals, not num events
		t = (n+1)/meanRate;
		writeFivePointRates("Charleston Post-5,500 years Earthquakes 1886, A, B, C, D and E; T = "+(float)t, n, t, Double.NaN, targetBranchRates5);


	}
	
	public static void testFivePointRatesAlt(String name, int n, double t1, double t2, double[] testVals) {
		double[]  rates1 = getFivePointRates(n, t1, t2);
		double[]  rates2 = getFivePointRatesAlt(n, t1, t2);
		System.out.println(name);
		for(int i=0;i<rates1.length;i++)
			System.out.println((float)rates1[i]+"\t"+(float)rates2[i]+"\t"+(float)(rates1[i]/rates2[i]));
	}
	
	

	public static void main(String[] args) {
		
		test2023_CEUS_results();
		
		
//		plotRateDistributions(2,10000d);

		
		// This tests that getFivePointRates() and getFivePointRatesAlt() produce same results for 2023 NSHM models
//		testFivePointRatesAlt("Marianna: 2 in 9.6k-10.2k",2, 9.6e3, 10.2e3, new double[] {7.2E-05,1.6E-04,2.7E-04,4.2E-04,6.9E-04});
//		testFivePointRatesAlt("Marianna: 3 in 9.6k-10.2k",3, 9.6e3, 10.2e3, new double[] {1.2E-04,2.4E-04,3.7E-04,5.5E-04,8.4E-04});
//		testFivePointRatesAlt("Wabash Valley 1 in 11k-13k", 1, 11e3, 13e3, new double[] {2.4E-05,7.2E-05,1.4E-04,2.5E-04,4.4E-04});
//		testFivePointRatesAlt("Saline River 1 in 322-5100", 1, 322, 5100, new double[] {1.00E-04,3.10E-04,6.80E-04,1.52E-03,4.84E-03});
//		testFivePointRatesAlt("Central Virginia 1 in 1800-2800", 1, 1800, 2800, new double[] {1.4E-04,3.9E-04,7.5E-04,1.3E-03,2.3E-03});
//		testFivePointRatesAlt("Meers 1 in 2153-2968", 1, 2153, 2968, new double[] {1.2E-04,3.4E-04,6.7E-04,1.2E-03,2.1E-03});
//		testFivePointRatesAlt("Charlevoix 3 in 6k-7k", 3, 6e3, 7e3, new double[] {1.9E-04,3.7E-04,5.7E-04,8.4E-04,1.3E-03});
//		testFivePointRatesAlt("Charlevoix 4 in 9.5k-10.2k", 4, 9.5e3, 10.2e3, new double[] {1.8E-04,3.2E-04,4.7E-04,6.7E-04,9.8E-04});
		
		
	}

}
