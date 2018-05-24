package scratch.aftershockStatistics;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.data.Range;
//import org.mongodb.morphia.annotations.Transient;
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
 * This represents a Reasenberg-Jones (1989, 1994) aftershock model where the a-value distribution is the Bayesian
 * combination of two given models, and where b, p, and c are held fixed (and they must be the same in each model).
 * 
 * a-value discretization is the smallest between the two given models, and the a-value range is the intersection of the two.
 * 
 * Note also that the probability distribution is renormalized so that values sum to 1.0 over the final range of
 * a-values represented.
 *
 * The probability density function of a in the Bayesian model, P(a), is proportional to:
 *  P(a) = P1(a)*P2(a)
 * where P1(a) and P2(a) are the probability density functions of a in the two models.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests that cover all cases.
 *
 * @author field
 *
 * Modified by Michael Barall.
 *
 */
public class RJ_AftershockModel_Bayesian extends RJ_AftershockModel {

//	//@Transient
//	protected boolean D=true;	// debug flag (inherited)




	/**
	 * Return the name of this model.
	 */
	@Override
	public String getModelName() {
		return "Reasenberg-Jones (1989, 1994) aftershock model (Bayesian Combination)";
	}




	/**
	 * Return true if x1 and x2 are approximately equal (to 7 digits).
	 */
    public static boolean areParamsEquivalent(double x1, double x2) {
		if (Math.abs(x1 - x2) <= 1.0e-7 * Math.max(Math.abs(x1), Math.abs(x2)) + Double.MIN_NORMAL) {
			return true;
		}
		return false;
    }




    /**
	 * Check similarity of both models, return true if the parameters are equivalent
	 * 
	 * @param model1 = First model to check.
	 * @param model2 = Second model to check.
	 * @return
	 * Returns true if the models are similar, so that the Bayesian combination can be formed.
	 */
	public static boolean areModelsEquivalent(RJ_AftershockModel model1, RJ_AftershockModel model2) {
		if (!( areParamsEquivalent (model1.getMainShockMag(), model2.getMainShockMag())
			&& areParamsEquivalent (model1.get_b(), model2.get_b())
			&& areParamsEquivalent (model1.getMaxLikelihood_p(), model2.getMaxLikelihood_p())
			&& areParamsEquivalent (model1.getMaxLikelihood_c(), model2.getMaxLikelihood_c()) )) {
			return false;
		}
		if (!( Math.max(model1.min_a, model2.min_a) < Math.min(model1.max_a, model2.max_a)
			&& model1.num_a > 1 && model2.num_a > 1 )) {
			return false;
		}
		return true;
	}



	
	/**
	 * This instantiates a Bayesian combination of the two given RJ models
	 * @param model1 = First model to combine.
	 * @param model2 = Second model to combine.
	 */
	public RJ_AftershockModel_Bayesian(RJ_AftershockModel model1, RJ_AftershockModel model2) {
		setup_model(model1, model2);
	}




	/**
	 * This default constructor creates an empty model.
	 * This is intended for use in database retrieval.
	 */
    public RJ_AftershockModel_Bayesian() {
		// When retrieving from database, remain quiet by default
		D = false;
    }




	/**
	 * Set up the model.
	 * @param model1 = First model to combine.
	 * @param model2 = Second model to combine.
	 */
    public void setup_model(RJ_AftershockModel model1, RJ_AftershockModel model2) {

		// check similarity of the two models
		Preconditions.checkArgument(areModelsEquivalent(model1, model2),
				"Models are not equivalent so Bayesian combination impossible."
						+"\n\tMain shock mags: %s, %s"
						+"\n\tb-values: %s, %s"
						+"\n\tp-values: %s, %s"
						+"\n\tc-values: %s, %s"
						+"\n\tmin a-values: %s, %s"
						+"\n\tmax a-values: %s, %s"
						+"\n\tnum a-values: %s, %s",
						model1.getMainShockMag(), model2.getMainShockMag(),
						model1.get_b(), model2.get_b(),
						model1.getMaxLikelihood_p(), model2.getMaxLikelihood_p(),
						model1.getMaxLikelihood_c(), model2.getMaxLikelihood_c(),
						model1.min_a, model2.min_a,
						model1.max_a, model2.max_a,
						model1.num_a, model2.num_a);
		
		// Use averages here so that the construction is symmetric
		this.magMain = 0.5*(model1.getMainShockMag() + model2.getMainShockMag());
		this.b = 0.5*(model1.get_b() + model2.get_b());
		double p = 0.5*(model1.getMaxLikelihood_p() + model2.getMaxLikelihood_p());
		double c = 0.5*(model1.getMaxLikelihood_c() + model2.getMaxLikelihood_c());

		set_fixed_p(p);
		set_fixed_c(c);
		
		
		HistogramFunction aValFunc1 =  model1.getPDF_a();
		if (aValFunc1 == null) {
			throw new RuntimeException("RJ_AftershockModel_Bayesian: aValFunc1 == null");
		}

		HistogramFunction aValFunc2 =  model2.getPDF_a();
		if (aValFunc2 == null) {
			throw new RuntimeException("RJ_AftershockModel_Bayesian: aValFunc2 == null");
		}

		this.min_a = Math.max(aValFunc1.getMinX(), aValFunc2.getMinX());
		this.max_a = Math.min(aValFunc1.getMaxX(), aValFunc2.getMaxX());
		if(min_a >= max_a) {
			throw new RuntimeException("RJ_AftershockModel_Bayesian: aValueMin >= aValueMax");
		}

		double minDelta = Math.min(aValFunc1.getDelta(), aValFunc2.getDelta());
		this.num_a = (int)Math.ceil((max_a-min_a)/minDelta) + 1;
		this.delta_a = (max_a - min_a)/((double)(num_a - 1));
		
		apc_likelihood = new double[num_a][num_p][num_c];
		for (int aIndex = 0; aIndex < num_a; aIndex++) {
			double a = get_a(aIndex);
			double wt = aValFunc1.getInterpolatedY(a)*aValFunc2.getInterpolatedY(a);
			apc_likelihood[aIndex][0][0] = wt;
		}

		// Complete the likelihood setup

		apcFinish (false);

		return;
	}
	
		


	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("RJ_AftershockModel_Bayesian : Missing subcommand");
			return;
		}


		// Subcommand : Test #1
		// Command format:
		//  test1
		// This is the pre-existing test for the class.

		if (args[0].equalsIgnoreCase ("test1")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("RJ_AftershockModel_Bayesian : Invalid 'test1' subcommand");
				return;
			}

			// Run the test
		
			RJ_AftershockModel_Generic gen1 = new RJ_AftershockModel_Generic(7.0, -2, 0.3, -4.5, -0.5, 0.01, 1.0, 1.12, 0.018);
			RJ_AftershockModel_Generic gen2 = new RJ_AftershockModel_Generic(7.0, -3, 0.3, -4.5, -0.5, 0.01, 1.0, 1.12, 0.018);
			RJ_AftershockModel_Bayesian bayes = new RJ_AftershockModel_Bayesian(gen1,gen2);
		
			ArrayList<HistogramFunction> funcList = new ArrayList<HistogramFunction>();
			funcList.add(gen1.getPDF_a());
			funcList.add(gen2.getPDF_a());
			funcList.add(bayes.getPDF_a());
			ArrayList<PlotCurveCharacterstics> curveCharList = new ArrayList<PlotCurveCharacterstics>();
			curveCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			curveCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			curveCharList.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			GraphWindow numVsTimeGraph = new GraphWindow(funcList, "PDF"); 
			numVsTimeGraph.setX_AxisLabel("a-value");
			numVsTimeGraph.setY_AxisLabel("Density");

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("RJ_AftershockModel_Bayesian : Unrecognized subcommand : " + args[0]);
		return;

	}

}
