package scratch.aftershockStatistics;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.LegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.jfree.data.Range;
import org.mongodb.morphia.annotations.Transient;
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
 * a-value discretization is the smallest between the two given models, and the a-value range if union of the two.
 * 
 * Note also that the Gaussian distribution is renormalized so that values sum to 1.0 over the final range of
 * a-values represented.
 * 
 * TODO:
 * 
 *  1) Carefully define jUnit tests that cover all cases.
 *
 * @author field
 *
 */
public class RJ_AftershockModel_Bayesian extends RJ_AftershockModel {

	@Transient
	Boolean D=true;	// debug flag

    public RJ_AftershockModel_Bayesian() {

    }

    /**
	 * Check similarity of both models, currently to floating point precision
	 * 
	 * @param model1
	 * @param model2
	 * @return
	 */
	public static boolean areModelsEquivalent(RJ_AftershockModel model1, RJ_AftershockModel model2) {
		if ((float)model1.getMainShockMag() != (float)model2.getMainShockMag())
			return false;
		if ((float)model1.get_b() != (float)model2.get_b())
			return false;
		if ((float)model1.getMaxLikelihood_p() != (float)model2.getMaxLikelihood_p())
			return false;
		if ((float)model1.getMaxLikelihood_c() != (float)model2.getMaxLikelihood_c())
			return false;
		return true;
	}
	
	/**
	 * This instantiates a Bayesian combination of the two given RJ models
	 * @param model1
	 * @param model2
	 */
	public RJ_AftershockModel_Bayesian(RJ_AftershockModel model1, RJ_AftershockModel model2) {
		
		this.magMain = model1.getMainShockMag();
		this.b = model1.get_b();
		double p = model1.getMaxLikelihood_p();
		double c = model1.getMaxLikelihood_c();
		
		// check similarity of these with the second model
		Preconditions.checkArgument(areModelsEquivalent(model1, model2),
				"Models are not equivalent so Bayesian combination impossible."
						+"\n\tMain shock mags: %s, %s"
						+"\n\tb-values: %s, %s"
						+"\n\tp-values: %s, %s"
						+"\n\tc-values: %s, %s",
						model1.getMainShockMag(), model2.getMainShockMag(), model1.get_b(), model2.get_b(),
						model1.getMaxLikelihood_p(), model2.getMaxLikelihood_p(),
						model1.getMaxLikelihood_c(), model2.getMaxLikelihood_c());
		
		this.min_p=p;
		this.max_p=p;
		this.num_p=1;
		this.min_c=c;
		this.max_c=c;
		this.num_c=1;
		
		
		HistogramFunction aValFunc1 =  model1.getPDF_a();
		HistogramFunction aValFunc2 =  model2.getPDF_a();
		this.min_a = Math.max(aValFunc1.getMinX(), aValFunc2.getMinX());
		this.max_a = Math.min(aValFunc1.getMaxX(), aValFunc2.getMaxX());
		double minDelta = Math.min(aValFunc1.getDelta(), aValFunc2.getDelta());
		this.num_a = (int)Math.ceil((max_a-min_a)/minDelta) + 1;
		EvenlyDiscretizedFunc aValueFuncBayes = new EvenlyDiscretizedFunc(min_a, max_a, num_a);
		this.delta_a = aValueFuncBayes.getDelta();

		if(min_a>max_a) {
			throw new RuntimeException("Problem: aValueMin > aValueMax");
		}
		
		for(int i=0;i<aValueFuncBayes.size();i++) {
			double a = aValueFuncBayes.getX(i);
			double wt = aValFunc1.getInterpolatedY(a)*aValFunc2.getInterpolatedY(a);
			aValueFuncBayes.set(i,wt);
		}
		
		setArrayAndMaxLikelyValuesFrom_aValueFunc(aValueFuncBayes, b, p, c);
		
	}
	
		

	public static void main(String[] args) {
		RJ_AftershockModel_Generic gen1 = new RJ_AftershockModel_Generic(7.0, -2, 0.3, -4.5, -0.5, 1.0, 1.12, 0.018);
		RJ_AftershockModel_Generic gen2 = new RJ_AftershockModel_Generic(7.0, -3, 0.3, -4.5, -0.5, 1.0, 1.12, 0.018);
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
	}

}
