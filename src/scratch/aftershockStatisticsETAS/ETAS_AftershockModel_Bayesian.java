package scratch.aftershockStatisticsETAS;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.LegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
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

import com.google.common.base.Preconditions;


/**
 * This represents an ETAS aftershock model where the solution is the Bayesian combination of the generic and sequence-specific models.
 * The generic model is treated as if the covariance between a, p, and c were zero. This means that the generic model uncertainties given 
 * in the generic model data file should be based on the marginal distributions of the parameters. 
 * 
 * a-value discretization is the same as for the sequence-specific model
 * 
 * Note also that the Gaussian distribution is renormalized so that values sum to 1.0 over the final range of
 * a-values represented.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests that cover all cases.
 *
 * @author field, vdElst
 *
 */
public class ETAS_AftershockModel_Bayesian extends ETAS_AftershockModel{
	
	Boolean D=true;	// debug flag
	
	
//	double maxMag;
//	double magComplete;
//	int maxGenerations;
//	int nSims;
	
	/**
	 * This instantiates a Bayesian combination of the two given ETAS models
	 * @param g_model
	 * @param ss_model
	 */
	public ETAS_AftershockModel_Bayesian(ETAS_AftershockModel_SequenceSpecific ss_model, ETAS_AftershockModel_Generic g_model) {
		
		this.magMain = g_model.getMainShockMag();
		this.b = g_model.get_b();
		
		// check similarity of these with the second model.
		Preconditions.checkArgument(areModelsEquivalent(g_model, ss_model),
				"Models are not equivalent so Bayesian combination impossible.");

		
		double like1, like2;
		double[][][] likelihood = new double[ss_model.a_vec.length][ss_model.p_vec.length][ss_model.c_vec.length];
		double cumSum = 0;
		
		for(int i = 0; i < ss_model.a_vec.length; i++ ){
			for(int j = 0; j < ss_model.p_vec.length; j++ ){
				for(int k = 0; k < ss_model.c_vec.length; k++){
					like1 = g_model.get_likelihood(ss_model.a_vec[i], ss_model.p_vec[j], ss_model.c_vec[k]);
					like2 = ss_model.likelihood[i][j][k];
					likelihood[i][j][k] = like1*like2;
					
					cumSum += likelihood[i][j][k];
				}
			}
		}
		
		//normalize likelihood and find maximum
		double maxVal = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < ss_model.a_vec.length; i++ ){
			for(int j = 0; j < ss_model.p_vec.length; j++ ){
				for(int k = 0; k < ss_model.c_vec.length; k++){
					likelihood[i][j][k] /= cumSum;

					if(maxVal<likelihood[i][j][k]){
						maxVal=likelihood[i][j][k];
						this.max_a_index=i;
						this.max_p_index=j;
						this.max_c_index=k;
					}

				}
			}
		}
		
		
		
		this.likelihood = likelihood;
		this.mainShock = ss_model.mainShock;
		this.magMain = ss_model.magMain;
		this.aftershockList = ss_model.aftershockList;
		this.dataEndTimeDays = ss_model.dataEndTimeDays;
		this.dataStartTimeDays = ss_model.dataStartTimeDays;
		this.forecastMinDays = ss_model.forecastMinDays;
		this.forecastMaxDays = ss_model.forecastMaxDays;
		this.maxMag = ss_model.maxMag;
		this.magComplete = ss_model.magComplete;
		this.maxGenerations = ss_model.maxGenerations;
		this.nSims = ss_model.nSims;
		this.a_vec = ss_model.a_vec;
		this.p_vec = ss_model.p_vec;
		this.c_vec = ss_model.c_vec;
		this.alpha = ss_model.alpha;
		this.refMag = ss_model.refMag;
		this.min_a = a_vec[0];
		this.max_a = a_vec[a_vec.length-1];
		this.num_a = a_vec.length;
		this.min_p = p_vec[0];
		this.max_p = p_vec[p_vec.length-1];
		this.num_p = p_vec.length;
		this.min_c = c_vec[0];
		this.max_c = c_vec[c_vec.length-1];
		this.num_c = c_vec.length;
		this.magAftershocks = ss_model.magAftershocks;
		this.relativeTimeAftershocks = ss_model.relativeTimeAftershocks;
		
//		// generate simulated ETAS catalogs
//		ETAScatalog simulatedCatalog = new ETAScatalog(a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
//				mainShock, aftershockList, dataEndTimeDays+forecastMinDays, dataEndTimeDays+forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
//		this.simulatedCatalog = simulatedCatalog;
				
			
		computeNewForecast(dataStartTimeDays,dataEndTimeDays, forecastMinDays,  forecastMaxDays, nSims);
		
	}		

	public void computeNewForecast(double dataMinDays, double dataMaxDays, double forecastMinDays, double forecastMaxDays, int nSims){
		System.out.println("Data/Forecast duration: " + dataMinDays +" "+ dataMaxDays +" "+ forecastMinDays +" "+ forecastMaxDays +" "+ nSims);
		System.out.println("Params: "+ mean_a +" "+ getMaxLikelihood_a() +" "+ getMaxLikelihood_p() +" "+ getMaxLikelihood_c() +" "+ alpha +" "+ b +" "+ magComplete);
		
		
		ETAScatalog simulatedCatalog = new ETAScatalog(a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
				mainShock, aftershockList, dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, magComplete, maxMag, maxGenerations, nSims); //maxMag = 9.5, maxGeneratons = 100;
		
		this.forecastMinDays = forecastMinDays;
		this.forecastMaxDays = forecastMaxDays;
		this.simulatedCatalog = simulatedCatalog;
		this.nSims = nSims;

		
	}

	
	private boolean areModelsEquivalent(ETAS_AftershockModel g_model, ETAS_AftershockModel ss_model) {
		// TODO Auto-generated method stub
		
//		if(g_model.get_b() != ss_model.get_b())
//			return false;
//		else
			return true;
	}
	
		

//	public static void main(String[] args) {
//		RJ_AftershockModel_Generic gen1 = new RJ_AftershockModel_Generic(7.0, -2, 0.3, -4.5, -0.5, 1.0, 1.12, 0.018);
//		RJ_AftershockModel_Generic gen2 = new RJ_AftershockModel_Generic(7.0, -3, 0.3, -4.5, -0.5, 1.0, 1.12, 0.018);
//		RJ_AftershockModel_Bayesian bayes = new RJ_AftershockModel_Bayesian(gen1,gen2);
//		
//		ArrayList<HistogramFunction> funcList = new ArrayList<HistogramFunction>();
//		funcList.add(gen1.getPDF_a());
//		funcList.add(gen2.getPDF_a());
//		funcList.add(bayes.getPDF_a());
//		ArrayList<PlotCurveCharacterstics> curveCharList = new ArrayList<PlotCurveCharacterstics>();
//		curveCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//		curveCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//		curveCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//		GraphWindow numVsTimeGraph = new GraphWindow(funcList, "PDF"); 
//		numVsTimeGraph.setX_AxisLabel("a-value");
//		numVsTimeGraph.setY_AxisLabel("Density");
//	}

}
