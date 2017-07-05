package scratch.aftershockStatisticsETAS;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.LegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;


/**
 * This represents an ETAS aftershock model where a-values are assumed to 
 * be Gaussian distributed and p and c are held fixed at the values given.
 * 
 * a-value discretization is hard-coded as 0.01 (the resolution of value given for various regions by
 * Page et al. (2016).
 * 
 * Note also that the Gaussian distribution is renormalized so that values sum to 1.0 over the range of
 * a-values represented.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests that cover all cases.
 *
 * @author field
 * @author van der Elst
 *
 */
public class ETAS_AftershockModel_Generic extends ETAS_AftershockModel {
	
	Boolean D=false;	// debug flag
	
	/**
	 * This instantiates a generic ETAS model from a GenericETAS_Parameters object, where aValueMin and aValueMax
	 * are set as -4.5 and -0.5, respectively.
	 * @param mainShock
	 * @param aftershockList
	 * @param genericETAS_Parameters
	 */
//	public ETAS_AftershockModel_Generic(ETAScatalog simulatedCatalog, ObsEqkRupture mainShock, ObsEqkRupList aftershockList,
//			GenericETAS_Parameters genericETAS_Parameters) {
	public ETAS_AftershockModel_Generic(
			ObsEqkRupture mainShock, ObsEqkRupList aftershockList, GenericETAS_Parameters genericETAS_Parameters,
			double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, double Mc, 
			double maxSimMag, int maxNumGenerations, int nSims) {
		
		this(mainShock, aftershockList,
				genericETAS_Parameters.get_aValueMean(), genericETAS_Parameters.get_aValueSigma(), 
				genericETAS_Parameters.get_bValue(), genericETAS_Parameters.get_pValue(), genericETAS_Parameters.get_pValueSigma(),
				genericETAS_Parameters.get_cValue(), genericETAS_Parameters.get_logcValueSigma(), genericETAS_Parameters.get_alpha(),
				genericETAS_Parameters.get_refMag(),
				dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, Mc,
				maxSimMag, maxNumGenerations, nSims); 
				
				
	}
	
	
	/**
	 * This instantiates a generic ETAS model for the values given
	 * @param magMain - main shock magnitude
	 * @param mean_a - mean a-value for the Gaussian distribution
	 * @param sigma_a - a-value standard deviation for the Gaussian distribution
	 * @param b - b-value
	 * @param p - p-value
	 * @param c - c-value
	 */
	public ETAS_AftershockModel_Generic(ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double mean_a, double sigma_a, double b, double mean_p, double sigma_p, double mean_c, double sigma_logc, double alpha, double refMag,
			double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, double Mc,
			double maxMag, int maxGenerations, int nSims) {
		
//		this.simulatedCatalog = simulatedCatalog;
		
		//this.magMain= magMain;
		this.magMain = mainshock.getMag();
		this.magAftershocks = this.getAftershockMags(aftershocks);
		this.b=b;

		//correct the a-value if Mc is not the same as refMag
		mean_a += Math.log10((maxMag - refMag)/(maxMag - Mc));
		
		this.mean_a=mean_a;
		this.sigma_a=sigma_a;
		this.mean_p = mean_p;
		this.sigma_p = sigma_p;
		this.mean_c = mean_c;
		this.sigma_logc = sigma_logc;
		this.alpha = alpha;
		this.refMag = refMag;
		this.magComplete = Mc;
		this.nSims = nSims;
		this.maxMag = maxMag;
		this.maxGenerations = maxGenerations;
		this.dataStartTimeDays = dataMinDays;
		this.dataEndTimeDays = dataMaxDays;
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		
		
		// compute a_value likelihoods from Gaussian distribution (need to convert a_mean, a_sigma into a vector to be parallel with seq specific model
		this.min_a = mean_a - 3*sigma_a;
		this.max_a = mean_a + 3*sigma_a;
		
		System.out.println(min_a +" "+ max_a +" "+ num_a);
		double[] a_vec = ETAS_StatsCalc.linspace(min_a, max_a, num_a);
		
		this.min_p = mean_p - 3*sigma_p;
		this.max_p = mean_p + 3*sigma_p;
		
		double[] p_vec;
		if(sigma_p == 0)
			p_vec = new double[]{mean_p};
		else
			p_vec = ETAS_StatsCalc.linspace(min_p, max_p, num_p);
		
		this.min_c = Math.pow(10, Math.log10(mean_c) - 3*sigma_logc);
		this.max_c = Math.pow(10, Math.log10(mean_c) + 3*sigma_logc);
		
		double[] c_vec;
		if(sigma_logc == 0)
			c_vec = new double[]{mean_c};
		else
			c_vec = ETAS_StatsCalc.logspace(min_c, max_c, num_c);
		
		double[][][] likelihood = get_likelihoodMatrix(a_vec, p_vec, c_vec);		
			
		this.a_vec = a_vec;
		this.p_vec = p_vec;
		this.c_vec = c_vec;
		this.likelihood = likelihood;
		
		this.max_a_index = Math.round((a_vec.length-1)/2);
		this.max_p_index = Math.round((p_vec.length-1)/2);
		this.max_c_index = Math.round((c_vec.length-1)/2);
		
		
		
//		this.min_a = a_vec[0];
//		this.max_a = a_vec[a_vec.length-1];
//		this.num_a = a_vec.length;
		this.min_p = p_vec[0];
		this.max_p = p_vec[p_vec.length-1];
		this.num_p = p_vec.length;
		this.min_c = c_vec[0];
		this.max_c = c_vec[c_vec.length-1];
		this.num_c = c_vec.length;
		this.delta_a = (max_a - min_a)/(num_a-1);
		this.delta_c = (max_c - min_c)/(num_c-1);
		this.delta_p = (max_p - min_p)/(num_p-1);
		
		this.aftershockList=aftershocks;
		this.mainShock=mainshock;
		this.magMain = mainshock.getMag();
		
		if(min_a>max_a) {
			throw new RuntimeException("Problem: aValueMin > aValueMax");
		}
		if(sigma_a<=0){
			throw new RuntimeException("Problem: sigma_a must be greater than 0");
		}
		
		// get aftershock times and mags and store as simple doubles[]
		double[] relativeEventTimes = ETAS_StatsCalc.getDaysSinceMainShockArray(mainShock, aftershockList.getRupsAboveMag(magComplete));
		double[] magAftershocks = getAftershockMags(aftershockList.getRupsAboveMag(magComplete));

		List<double[]> sortedEQlist = new ArrayList<double[]>();

		System.out.println(relativeEventTimes.length + " "+ dataEndTimeDays +" "+ magComplete);

		for(int i = 0; i < relativeEventTimes.length; i++){
			double[] temp = new double[]{relativeEventTimes[i], magAftershocks[i]};
			if(temp[0] < dataEndTimeDays)
				sortedEQlist.add(temp);
		}

		//sort double[] of times and magnitudes
		Collections.sort(sortedEQlist, new java.util.Comparator<double[]>() {
			public int compare(double[] a, double[] b) {
				return Double.compare(a[0], b[0]);
			}
		});

		for(int i = 0; i < relativeEventTimes.length; i++){
			double[] temp = sortedEQlist.get(i);
			relativeEventTimes[i] = temp[0];
			magAftershocks[i] = temp[1];
			//					sortedEQlist.add(temp);
		}

		this.magAftershocks = magAftershocks;
		this.relativeTimeAftershocks = relativeEventTimes;
		
		// generate simulated ETAS catalogs
//		ETAScatalog simulatedCatalog = new ETAScatalog(mean_a, a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
//				mainshock, aftershocks, dataMaxDays+forecastMinDays, dataMaxDays+forecastMaxDays, Mc, 9.5, 100, nSims); //maxMag = 9.5, maxGeneratons = 100;
//		
//		this.simulatedCatalog = simulatedCatalog;
		computeNewForecast(dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, nSims);
		
	}		

	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
		System.out.println("Data/Forecast duration: " + dataMinDays +" "+ dataMaxDays +" "+ forecastMinDays +" "+ forecastMaxDays +" "+ nSims);
		System.out.println("Params: "+ mean_a +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete);
		
		ETAScatalog simulatedCatalog = new ETAScatalog(mean_a, a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
				mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
		
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.simulatedCatalog = simulatedCatalog;
		this.nSims = nSims;
	}

	/**
	 * Returns likelihood for given [a,p,c] assuming a Gaussian distribution on each parameter (specified by mean_a, sigma_a, mean_p, sigma_p, etc.)
	 */
	public double get_likelihood(double a, double p, double c){
		double aLike, pLike, cLike;
		double like;
		
		if(sigma_a == 0)
			if(mean_a == a)
				aLike = 1;
			else
				aLike = 0;
		else
			aLike = 1/Math.sqrt(2*Math.PI)/sigma_a * Math.exp( -(a - mean_a) * (a - mean_a ) / (2*sigma_a*sigma_a) );
		
		if(sigma_p == 0)
			if(mean_p == p)
				pLike = 1;
			else
				pLike = 0;
		else
			pLike = 1/Math.sqrt(2*Math.PI)/sigma_p * Math.exp( -(p - mean_p) * (p - mean_p ) / (2*sigma_p*sigma_p) );

			
		double logc = Math.log10(c);
		double mean_logc = Math.log10(mean_c);
		
		if(sigma_logc == 0)
			if(mean_logc == logc)
				cLike = 1;
			else
				cLike = 0;
		else
			cLike = 1/Math.sqrt(2*Math.PI)/sigma_logc * Math.exp( -(logc - mean_logc) * (logc - mean_logc ) / (2*sigma_logc*sigma_logc) );

		
		like = aLike*pLike*cLike;
		return like;
		
	}
	
	/**
	 * Returns likelihood matrix for vectors [a,p,c] assuming a Gaussian distribution on each parameter (specified by mean_a, sigma_a, mean_p, sigma_p, etc.)
	 */
	
	private double[][][] get_likelihoodMatrix(double[] a_vec, double[] p_vec, double[] c_vec){
		double[][][] likelihood = new double[a_vec.length][p_vec.length][c_vec.length];
		double aLike, pLike, cLike;	
		double cumSum = 0;
		double mean_logc = Math.log10(mean_c);
		double logc;

		
		
		
		for(int i = 0; i < a_vec.length ; i++ ){
			if(sigma_a == 0)
				if(mean_a == a_vec[i])
					aLike = 1;
				else
					aLike = 0;
			else
				aLike = 1/Math.sqrt(2*Math.PI)/sigma_a * Math.exp( -(a_vec[i] - mean_a) * (a_vec[i] - mean_a ) / (2*sigma_a*sigma_a) );
			
			for(int j = 0; j < p_vec.length ; j++ ){
				if(sigma_p == 0)
					if(mean_p == p_vec[j])
						pLike = 1;
					else
						pLike = 0;
				else
					pLike = 1/Math.sqrt(2*Math.PI)/sigma_p * Math.exp( -(p_vec[j] - mean_p) * (p_vec[j] - mean_p ) / (2*sigma_p*sigma_p) );
				
				for(int k = 0; k < c_vec.length ; k++ ){
					logc = Math.log10(c_vec[k]);
					if(sigma_logc == 0)
						if(mean_logc == logc)
							cLike = 1;
						else
							cLike = 0;
					else
						cLike = 1/Math.sqrt(2*Math.PI)/sigma_logc * Math.exp( -(logc - mean_logc) * (logc - mean_logc ) / (2*sigma_logc*sigma_logc) );
					
					likelihood[i][j][k] = aLike*pLike*cLike;
					cumSum += likelihood[i][j][k];
				}
			}
		}

		for(int i = 0; i < a_vec.length ; i++ ){
			for(int j = 0; j < p_vec.length ; j++ ){
				for(int k = 0; k < c_vec.length ; k++ ){
					likelihood[i][j][k] /= cumSum;
				}
			}
		}
		
		
		return likelihood;
	}
		
		

}
